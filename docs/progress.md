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
  - 单原子对象，`set_tau`/`set_frac_atomic_position` 设分数坐标，默认 τ=0 → S=1.0
  - 三种访问入口：不缓存 `operator()(g, k)`（直接 exp）、缓存 `prepare(gs)` + `operator()(ig, k)`（三数组查表）、批处理 `operator()(gs, k, out)`
- [ ] 命名（已完成，使用 StructureFactor 类）

### 非局域势计算

- 在 SphericalHarmonics 中加入 reset(l_max), 更改常驻 Y_lm(GKK) 的数量
- [x] RealSphericalHarmonics::get(l,m) → operator()(l,m) API 统一

- simpson积分 meshtype 从 ncpp 到 bessel_fourier 转换
```cpp
  auto mt = (ncpp.mesh.type == MeshType::Uniform)
            ? RadialMeshType::Uniform : RadialMeshType::General;
  simpson(f, rab, mt);
```


```cpp
for(ityp: ATOM.eachtype())
    for(iatm: ATOM.eachatom(ityp))
        V_psp = ncpp(ityp)
        Y_lm = SphericalHarmonics.reset(GKK, l_max)
        structure_factor = structure_factor(GKK)
        for(l=0; l<=l_max; ++l)
            V_NL_l = V_psp.projectorBlock(l)
            for(ib=0; ib<=nb; ++ib)
                beta_q_interpolaror = BetaQInterpolator(V_psp.mesh, V_NL_l.beta[ib])  // handle q_max 
                beta_q = (GKK => beta_q_interpolaror)
                for(m=-l; m<=l; ++m)
                    projector = (beta_q .* Y_lm .* structure_factor)
                    H_psi[:] +=
                        innerProduct(conj(projector), WG) * 
                        V_NL_l.B[ib, ib] * 
                        projector[:]
```




# 未来
我们会增加 $H^{\alpha} \psi$ 和 $H^{\alpha\beta} \psi$，alpha和beta指的是对于  以 k 为参数的哈密顿量 H(k)，对于 k_alpha 和 k_beta 分量的一阶偏导数和二阶偏导数

