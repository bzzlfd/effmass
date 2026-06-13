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
# 重构
## 1. advance 内嵌数据结构
ncpp.cppm 当中，advance 成员的设计                                                                                                           
  1. 我的本意是，因为ncpp.cppm 仅仅是对 ncpp 数据的表示，所以涉及对它数据处理的代码应该是 advance 的部分。顺着当时设计 alias                              
  的思路，我采用了类内内嵌结构的写法                                                                                                                      
  2. 然而在我开始写代码的时候，我发觉这样很怪。ncpp.advance.XXX 的写法1. 不自然 2.甚至会让人误解成这是某个迭代容器中的 advance 向前迭代                   


当前设计的问题根源

  Alias 作为内嵌 struct 很自然，因为：
  - alias 是名词 → 读作 ncpp 的别名视图                                                   
  - 它是只读访问 → 没有副作用
  - 含义清晰，不会与其他概念混淆

  而 Advance 作为内嵌 struct 别扭，因为：
  - advance 是动词 → 读作 ncpp 去推进，但 . 运算符使它在语法上像名词属性
  - 和 std::advance / 迭代器概念重名 → 确实会误解
  - 内嵌一个执行突变操作的结构体 → 数据对象内部携带"操作空间"的模式本身就很罕见

  关键矛盾在于：你用同一模式（内嵌 struct）承载了两种截然不同的语义——视图（名词）和操作（动词/祈使句）。

### 重构方案
自由函数 + 独立模块
```cpp
  // src/pseudo/ncpp_advance.cppm
  export module pseudo.ncpp.advance;
  import pseudo.ncpp;

  export void diagonalizeNonlocal(NCPP& ncpp);
  export void upf2upfSO(NCPP& ncpp);
```
  调用：
```cpp
  diagonalizeNonlocal(ncpp);
```


## 2. Betaq table

重构 hamiltonian_callable.cpp 中计算nonlocal pseudopotentail 的部分，涉及的代码包括 
  hamiltonian.cppm, hamiltonian_callable.cpp, fourier_bessel.cppm, ncpp_advance.cppm, ncpp.cppm
  1. 现在的计算瓶颈很大程度上在于每次要把ncpp中实空间的径向beta变换到q空间，而这里计算大量的simpson积分，这些计算被重复计算。然而，BetaqIntepolator不是
  每次at_k 得到 Callable对象之后需要重新计算的，它可以放在 Hamiltonian类内部，甚至直接在ncpp对象里面（具体我打算放到 ncpp_advance.cppm 下）
  2. 目前的fourier_bessel.cppm 的框架设计，对于计算一个type 的原子，即计算一个ncpp数据，是一个不错的框架。
    2.1. 一个 原子类型/ncpp 只有一种 mesh(r, rab)
    2.2. 球贝塞尔函数采用递推方式，仅仅在 l=0 初始化时因三角函数涉及较多计算量
    2.3. 我们可以让 beta 从 r 到 q 的计算（积分）从 l 小的beta开始，一步步升高，这样可以搭球贝塞尔函数递推的便车
  3. 目前的fourier_bessel.cppm 的框架设计，对于 多l,多beta table_缓存 来说，不是好设计，让我们先讨论这部分的设计
    3.1. BetaqInterpolator 降级成函数，无需保留内部状态
    3.2. BetaqInterpolator 升级成 多l, 多table 记录
    3.3. Betaq 和 interpolator 函数

### 重构方案
1. 添加 sort beta