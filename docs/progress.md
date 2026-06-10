# H_psi

H_psi 作为一个Functor（即实现了operator()的类）
构造是读取 GKK，ATOM，pseudo.ncpp文件们，Functor 对 WG 的数据进行变换。

Functor 对 WG 变换过程中，累积到自己内部的与 WG 同样尺寸的vector `H_psi`上，包含以下几项


```
H_psi[:] = 0
```

## 动能
### 已实现
- Callable 构造时加载 GKK k 点数据 (`loadKPoint`) + 验证 `kinetic = 0.5*(Kx²+Ky²+Kz²)` (`validateKineticConsistency`)
- 直接使用 GKK 的 kinetic 数据，逐 G-矢量累积到 hpsi：
```
hpsi[ig] += kinetic[ig] * psi[ig]
```
- `dim()` 返回 `ng_`（该 k 点的 G-矢量数量）

### 设计决策
- `gkk_` 标记为 `mutable`，`gkk() const` 返回 `GKK&`，允许在 `Callable::operator()`（const）中调用 `loadKPoint`
- `Callable` 保持 `const Hamiltonian*`，`at_k` 保持 const

## 局域势
### 已实现
通过傅里叶变换快速实现：
```
psi[G] --[FFT:G2R]--> psi(r) --> psi(r) * VR(r) --[FFT:R2G]--> (V_loc|ψ⟩)[G]
hpsi[G] += (V_loc|ψ⟩)[G]
```

- Callable 构造时启用 Integer view，捕获 g_idx 和 FFT 网格尺寸 n1,n2,n3
- `operator()` 中：放置系数 → G2R → 乘 VR → R2G → 提取回 hpsi
- FFT 变换使用 `math.FFT3D`（基于 pocketfft）
- G2R 后不含 1/N 因子（forward R2G 已包含缩放）


## Pseudo Potential 非局域项

### 结构因子
- [x] 实现 StructureFactor 类（`H_psi.structure_factor`）
  - 两种模式：`CacheMode::None`（直接 exp 计算，不缓存）和 `CacheMode::Separable`（分离式 1D 相位表缓存，默认模式）
  - Separable 模式利用恒等式：
    `S(g,k) = exp(-i·2π·(k+g)·τ) = exp(-i·2π·k·τ) × p_x[g.x] × p_y[g.y] × p_z[g.z]`
    其中 `p_x[i] = exp(-i·2π·(i - n1/2)·τ_x)`（p_y, p_z 同理）
  - 缓存：三个 1D 复数向量，大小 `n_i + 2*CACHE_BUFFER`，`CACHE_BUFFER = 4`，构造时预计算
  - `cacheIndex(g, n) = g + n/2 + CACHE_BUFFER` 映射 G 矢量分量到缓存索引
  - `reset_tau()` / `reset_frac_atomic_position()` 重建缓存


### 傅里叶-贝塞尔变换
在目前 fourier_bessel.cppm 中，实现的是用 $e^{i\bm{q}\cdot\bm{r}}$ 变换 $\beta(r)$。

然而平面波 $\ket{\bm{q}}=e^{i\bm{q}\cdot\bm{r}}$ 的归一化系数依赖于晶格体积。
$$
\int_\Omega \mathrm{d}\bm{r} \braket{\bm{q}|\bm{q}} = 1 \cdot \Omega
$$
所以归一化系数应该是 $1/\sqrt{\Omega}$。

我们在构造函数里多引入一个归一化系数参数。


### 非局域势计算

实现了非局域项在 `Callable::operator()` 中的计算，位置在动能和局域势项之后。

利用 `RealSphericalHarmonics`、`StructureFactor`、`BetaqInterpolator` 的 reset 接口避免循环中重复创建：


#### 数据流


```cpp
Y_lm = SphericalHarmonics.reset(GKK, l_max)
structure_factor = StructureFactor.init(n1, n2, n3, mode=Separable)
for(ityp: ATOM.eachtype())
    auto& V_psp = ncpp[ityp];
    V_psp.diagnoalizeNonlocal()
    beta_q_interpolaror = BetaqInterpolator.init(V_psp.mesh.r, V_psp.mesh.rab, l=0, q_max)
    for(iatm: ATOM.eachatom(ityp))
        structure_factor.reset_tau(tau)
        for(l=0; l<=l_max; ++l)
            V_NL_l = V_psp.projectorBlock(l)
            for(ib=0; ib<=nb; ++ib)
                beta_q_interpolaror.reset_beta(V_NL_l.beta[ib], l)  // handle q_max 
                beta_q = beta_q_interpolaror(GKK::KVecs.q)
                for(m=-l; m<=l; ++m)
                    projector = (beta_q .* Y_lm .* structure_factor)
                    H_psi[:] +=
                        innerProduct(conj(projector), WG) * 
                        V_NL_l.B[ib, ib] * 
                        projector[:]
```



# 未来
我们会增加 $H^{\alpha} \psi$ 和 $H^{\alpha\beta} \psi$，alpha和beta指的是对于  以 k 为参数的哈密顿量 H(k)，对于 k_alpha 和 k_beta 分量的一阶偏导数和二阶偏导数

