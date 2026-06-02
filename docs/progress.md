# H_psi

H_psi 作为一个Functor（即实现了operator()的类）
构造是读取 GKK，ATOM，pseudo.ncpp文件们，Functor 对 WG 的数据进行变换。

Functor 对 WG 变换过程中，累积到自己内部的与 WG 同样尺寸的vector `H_psi`上，包含以下几项
- [ ] `H_psi` 的命名

```
H_psi[:] = 0
```

## 动能
验证 GKK中 kinetic 是否可靠，然后直接使用这个kinetic数据。进行数据累加
```
H_psi[g] += kinetic[g] * WG[g]
```

## 局域势
通过傅里叶变换快速实现：
```
WG --[FFT:G2R]--> WR --> WR * VR --[FFT:R2G]--> WVG
H_psi[g] += WVG[g]
```


## Pseudo Potential 非局域项

### 结构因子
- [ ] 命名
Pseudo Potential 和原子位置有关。对于平面波 exp(iqr), 和原子位置 tau，
相对于原子位于原点，pseudopotential对平面波的作用会额外多一个 exp(iq tau)的系数

具体实现上，
1. 我们 GKK IntegerView 时把其中的数据向量看成 -K=-(k+G)
2. exp函数内部的加法可以拆分成外面的乘法
所以我们不采用 n1*n2*n3 大小的数组记录，而是使用三个数组，大小分别是 n1 n2 n3 记录 exp(-i b1 tau)^iG exp(-i b2 tau)^jG  exp(-i b3 tau)^kG，然后再乘 exp(-i k tau)

### 非局域势计算

- 在 SphericalHarmonics 中加入 reset(l_max), 更改常驻 Y_lm(GKK) 的数量

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
- [ ] 命名



# 球谐函数
## buildYlm
buildYlm 做了 批量化 计算实球谐函数：对所有 ng_ 个 G 点同时算 Y_lm，结果存成 y_lm_[block][ig]，其中 block = l*l + l + m。

  外层循环是 固定 m，遍历 l（column-wise），比固定 l 遍历 m 更高效。全程分 6 步：

  ---
  1. Seed P_m^m（line 477-485）
  
  对每个 G 点算 l=m 的 Legendre，封闭公式：
  P_m^m(x) = (2m-1)!! * sin(θ)^m
  
  m=0 时 P_0^0 = 1，m>0 时直接乘上去。

  ---
  2. P_{m+1}^m 递推一步（line 487-491）
  
  从 P_m^m 到 P_{m+1}^m 用初值公式：
  P_{m+1}^m(x) = (2m+1) * cos(θ) * P_m^m

  ---
  3. 三项递推到 l_max（line 493-501）
  
  更高 l 用标准 Legendre 三项递推：
  (l-m) * P_l^m = (2l-1) * cos(θ) * P_{l-1}^m - (l+m-1) * P_{l-2}^m
  
  对每个 G 点各自独立算。p_buf 拍平成 [l_idx * ng_ + ig]。

  ---
  4. 保存 top_p_ 状态（line 503-516）
  
  为 reset() 增量展开做准备。把当前 m 的最后一层 P_{l_max} 和倒数第二层 P_{l_max-1} 复制到 top_p_[m]，以便下次 expand 时不用重头算 Legendre。

  ---
  5. 组装 Y_lm（line 518-534）
  
  实球谐函数的完整公式：

  m = 0:    Y_{l0} = N_l^0 * P_l^0
  m > 0:    Y_{l,+m} = √2 * N_l^m * P_l^m * cos(mφ)   → block[l*l + l + m]
            Y_{l,-m} = √2 * N_l^m * P_l^m * sin(mφ)   → block[l*l + l - m]

  cos(mφ)/sin(mφ) 来自外层的 cm/sm 数组。

  ---
  6. 递进 cos(mφ), sin(mφ)（line 536-543）
  
  用三角恒等式从 m 递推到 m+1：
  cos((m+1)φ) = cos(mφ)cos(φ) - sin(mφ)sin(φ)
  sin((m+1)φ) = sin(mφ)cos(φ) + cos(mφ)sin(φ)

  ---
  整体数据流：

  theta_, phi_ ──→ cos(θ), sin(θ), cos(φ), sin(φ)   (构造函数，只跑一次)
                        │
                        ↓
                 m 外层循环:
                    1. P_m^m seed ──────────────┐
                    2. P_{m+1}^m                │
                    3. 三项递推 P_l^m           │
                    4. 存 top_p_  ←─────────────┘ (供 reset expand 用)
                    5. 组装 Y_lm → y_lm_[]
                    6. 更新 cos(mφ)/sin(mφ)



## Legendre 递推里实际值的范围：

  ┌─────┬──────────────────┬──────────────┬────────────┐
  │  m  │ P_m^m = (2m-1)!! │ 归一化 N_m^m │ 乘积 Y_m^m │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 1   │ 1                │ 0.345        │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 5   │ 945              │ 4.6e-4       │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 10  │ 6.5e8            │ 6.4e-10      │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 20  │ 3.2e21           │ 1.3e-22      │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 50  │ 5.1e62           │ 2.2e-63      │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 100 │ 7.3e144          │ 1.5e-145     │ O(1)       │
  ├─────┼──────────────────┼──────────────┼────────────┤
  │ 170 │ 1.1e258          │ ~1e-259      │ O(1)       │
  └─────┴──────────────────┴──────────────┴────────────┘

  P_l^m 和 N_l^m 都是各自往极端走，乘积始终 O(1)。对于当前代码的目标范围（MAX_L_SAFE = 10，实际用到的 l ≤ 3），P_l^m ~ 1e9 完全不溢出，不需要合并。

  但如果要支撑 l > 50 的场景，确实该做归一化 Legendre 递推——把 N_l^m 吸收到递推关系里，中间量始终 O(1)。代码里像 top_p_ 存的就是未归一化的
  P，如果要扩展，这里也得跟着改。

## top_p_ 
  - top_p_[m].curr = $P_{l_{\text{max}}}^{,m}$ （最高的）
  - top_p_[m].prev = $P_{l_{\text{max}}-1}^{,m}$ （次高的，$l_{\text{max}} = m$ 时等于 curr）

  存这两个的目的，是给 reset() 增量展开时做三项递推往下走。比如旧 $l_{\text{max}}=4$、新 $l_{\text{max}}=6$，就从 top_p_[m] 存的 $P_4^m$、$P_3^m$
  出发，递推算出 $P_5^m$、$P_6^m$，不用重头开始。

