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

---

# HPSI_DEBUG 诊断系统

HPSI_DEBUG 是一个编译期开关，在 `Callable::operator()` 中插入诊断输出，用于验证 H|ψ⟩ 各分项的正确性。

## 编译选项

CMake 变量，支持三级（0/1/2）：

```bash
# 关闭（默认）
cmake -B build -G Ninja -DHPSI_DEBUG=0

# 级别 1：逐测试项期望值汇总
cmake -B build -G Ninja -DHPSI_DEBUG=1

# 级别 2：逐原子 逐(l,m,ib) 分解
cmake -B build -G Ninja -DHPSI_DEBUG=2
```

配置后重新编译：

```bash
cmake --build build
```

## 输出说明

### 级别 1（`HPSI_DEBUG=1`）

在每个 H|ψ⟩ 作用后打印逐项期望值汇总表：

```
[HPSI_DEBUG] ikpt=1  iband(WG)=24  ng=1714
[HPSI_DEBUG]   ⟨ψ|ψ⟩       =       0.9999999535  (should = 1)
[HPSI_DEBUG]   ─────────────── expectation values ───────────────
[HPSI_DEBUG]   ⟨ψ|T|ψ⟩     =       1.5888611304
[HPSI_DEBUG]   ⟨ψ|V_loc|ψ⟩ =      -0.5431430021
[HPSI_DEBUG]   ⟨ψ|V_NL|ψ⟩  =       0.0586388123
[HPSI_DEBUG]   ─────────────── eigenvalues (÷⟨ψ|ψ⟩) ─────────────
[HPSI_DEBUG]   E_kin       =       1.5888612043
[HPSI_DEBUG]   E_loc       =      -0.5431430274
[HPSI_DEBUG]   E_NL        =       0.0586388150
[HPSI_DEBUG]   ──────────────────────────────────────────────────
[HPSI_DEBUG]   E_sum       =       1.1043569919
[HPSI_DEBUG]   E_total     =       1.1043569919  (Rayleigh)
[HPSI_DEBUG]   Im(E_total) =  -8.2650116077e-18
```

列字段：
- `⟨ψ|ψ⟩` × Ω — 波函数模方×体积，应=1
- `⟨ψ|T|ψ⟩` — 动能期望值
- `⟨ψ|V_loc|ψ⟩` — 局域势期望值
- `⟨ψ|V_NL|ψ⟩` — 非局域势期望值
- `E_sum` = `(E_kin + E_loc + E_NL)`，拼接检验
- `E_total` = `⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩` 完整 Rayleigh 商
- `Im(E_total)` — 接近 0 说明厄米性

配套脚本 `tools/parse_hpsi_debug.sh`（依赖 `tools/parse_hpsi_debug.awk`）可从 ctest 输出中提取对齐表格：

```bash
bash tools/parse_hpsi_debug.sh
```

输出类似：
```
ikpt iband             <T>        <V_loc>         <V_NL>         E_file         E_calc         DE
   0     0      2.55493117    -5.68917223    -1.82281428    -3.47498702    -4.95705549   1.482068
   0     1      2.55513613    -5.68914516    -1.82310777    -3.47487044    -4.95711749   1.482247
```

`DE = E_file − E_calc` 应在容差内接近 0。较大偏差表示对应能带的 H|ψ⟩ 结果有问题。

### 级别 2（`HPSI_DEBUG=2`）

在级别 1 的基础上，额外输出 V_NL 循环内部的逐原子逐投影仪分解。

#### 格式

```
[HPSI_DEBUG:2]   atom ityp=0 tau=(0.666667,0.333333,0.880713)
[HPSI_DEBUG:2]     l=1 m=-1 ib=0  lambda=+0.9621964223  |<p|psi>|^2=8.089e-06  <p|p>=0.7618
[HPSI_DEBUG:2]     l=1 m=+0 ib=0  lambda=+0.9621964223  |<p|psi>|^2=1.559e-11  <p|p>=0.7595
[HPSI_DEBUG:2]     l=1 m=+1 ib=0  lambda=+0.9621964223  |<p|psi>|^2=5.039e-05  <p|p>=0.7618
[HPSI_DEBUG:2]   atom ityp=0 tau=(0.666667,0.333333,0.880713)  total_V_NL=5.091e-05  ×Ω=1.461e-02
```

每个字段：
- `ityp` — 原子类型索引
- `tau=(x,y,z)` — 该原子在分数坐标下的位置
- `l,m,ib` — 角量子数/磁量子数/投影仪通道索引
- `lambda` — `diagonalizeNonlocal()` 对角化后 B 矩阵的本征值 λ_i
- `|<p|psi>|^2` — 已乘体积 Ω 的投影仪与波函数内积模方 `|⟨p|ψ⟩|² × Ω`
- `<p|p>` — 投影仪自身模方 `Σ_g |p(g)|²`，用于归一化校验
- `total_V_NL(atom)` — 原始值 Σ(λ · |⟨p|ψ⟩|²)；`×Ω` — 乘体积后的期望值，与 L1 `⟨ψ|V_NL|ψ⟩` 直接对比

#### 用途

级别 2 适用于定位 V_NL 实现中的具体故障：

| 现象 | 可能的根因 |
|------|-----------|
| ⟨p|p⟩ 偏离 1 且分散 | β(q) 傅里叶-贝塞尔变换归一化系数错误 |
| 同 l 不同 m 的 ⟨p|p⟩ 不等 | 球谐函数归一化与 β(q) 的交互问题 |
| λ 值异常（如比预期大 100×） | `diagonalizeNonlocal()` 中 B 矩阵对角化问题 |
| |<p|psi>|^2 为 0 的通道 | k 点/能带对称性导致的正确行为，或投影仪方向错误 |
| 不同 tau 的同 (l,m,ib) 贡献相同 | 结构因子相位未正确区分原子位置 |

## 性能影响

- `HPSI_DEBUG=1`：多一次 `dbg_vloc_contrib` 存储 + 期望值后处理，开销约 5-10%
- `HPSI_DEBUG=2`：多一次额外循环计算 ⟨p|p⟩ + 大量 `std::println`，在非局域项较多的体系上开销可达 30-50%
- `HPSI_DEBUG=0`：编译期完全移除，零开销

