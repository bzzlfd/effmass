# RealSphericalHarmonics 接口设计

## 概述

`RealSphericalHarmonics` 是向量化的实球谐函数求值器，对一组 K 点的 G 矢量（theta, phi）同时计算 $Y_{lm}$。类位于模块 `math.sph_harmonics`，是 `math` 模块的一部分。

## 两种缓存策略，Ylm数据获得方式

### 1. None Cache：compute 模式

```
构造 → 只存 trig 数组 → compute(l, m) 是唯一访问路径
```

None 模式下，构造器只做一件事：**预计算三角函数**（cosθ, sinθ, cosφ, sinφ）。`setLMax()` 不会被调用，`y_lm_` 和 `top_column_Q_` 保持空。

访问完全依赖 `compute(l, m)`：

```
compute(l, m)
  ├── seedPmm(m, q_curr)       ← Q_m^m = N_m^m · sin^m(θ) · 1/(2√π)
  ├── stepPmm1 / advanceColumn  ← 递推到 Q_l^m
  └── 组装: Q_l^m → Y_lm       ← 乘以 √2 · {cos(mφ), sin(mφ)}
```

每次都是完整计算，O(l · ng)，结果不缓存。适合访问稀疏、l 上限不可预测、或内存受限的场景。

---

### 2. Full Cache：setLMax → get

这是设计的核心，分为冷热两条路径。

#### get() — 热路径

```
setLMax(l_max)
  └── 逐列递推: m=0..l_max
        ├── seedPmm  →  Q_m^m
        ├── stepPmm1 →  Q_{m+1}^m
        ├── advanceColumn → Q_{m+2..l_max}^m
        └── assembleYlm → 写入 y_lm_[l*l + (m+l)]

get(l, m)
  └── y_lm_[l*l + (m+l)]   ← O(1)，零拷贝
```

`setLMax()` 完成后，所有 `l ≤ l_max_resident` 的 Y_lm 在 `y_lm_` 数组中立即可用。`get()` 就是一次数组索引。

`setLMax()` 的增量设计：**保存 `top_column_Q_`**（每个 m 在 `l_max_resident` 处的 Legendre 递推顶端状态）。后续扩展不需要从头再来，只需从当前顶端继续递推，复杂度 O(Δl · l_max · ng)。缩小同样通过 `retreatColumn` 反向递推恢复状态。

> [!NOTE] **`top_column_Q_` 是 Full Cache 的专属状态**：None Cache 模式下保持为空，只有 Full Cache 的 `setLMax()` 才会写入和更新。
> 它是一个按 m 索引的结构（m = 0..l_max_resident）——对每个 m 保存 `l = l_max_resident` 处的 Q 矢量（大小 ng），记录递推中断处的顶端行。后续 `setLMax(l)` 扩展或收缩时，算法从这些顶端值继续递推，而非从头重算。

#### compute() — 冷路径 (bonus of NoneCache)

当 `l > l_max_resident_`，`get()` 不工作。调用者有两个选择：

| 路径 | 开销 | 适用场景 |
|------|------|---------|
| `setLMax(l); get(l, m)` | 一次性 O(Δl·l_max·ng)，后续 O(1) | 后续会反复访问高 l |
| `compute(l, m)` | 每次 O(l·ng)，不占用缓存 | 偶发一次查询，不确定是否会复用 |

`compute()` 在 Full 模式下的角色是：**缓存的补充机制**，覆盖缓存范围之外的即时查询。它不是主路径，而是处理边界情况的逃生口。


#### 必要之物：setLMax()
初始化构建缓存过程中，会设置 l_max 以确定缓存范围。这段代码被抽象出来，并添加了 reset l_max 逻辑，作为公共接口暴露出来，给使用者带来更多灵活性。
新创建的缓存向量，会以 cos_theta 的 capacity 作为参考，设置自己的 capacity。

---

### 3. operator() — 统一接口

```
operator()(l, m)
  ├── 校验量子数 (l, m)
  ├── Full 模式
  │     ├── l > l_max_resident_ → 报错，提示 setLMax 或 compute
  │     └── 正常 → return get(l, m)   ← 拷贝一次
  └── None 模式 → return compute(l, m)
```

定位：**便利入口**，调用者可以不关心 CacheMode 直接用 `sh(l, m)` 拿到 Y_lm。代价是返回值而非引用（Full 模式下多做一次拷贝），性能敏感场合应直接调用 `get()` 或 `compute()`。

`operator()` 独自做校验而非委托给 `get()`，保证了报错信息中的函数名与实际调用一致。

### 总结

| 方法 | 返回类型 | CacheMode | l 范围 | 开销特征 |
|------|---------|-----------|--------|---------|
| `get(l, m)` | `const vector<double>&` | Full 必须 | ≤ l_max_resident | O(1)，零拷贝 |
| `operator()(l, m)` | `vector<double>` | 任意 | 任意 | Full: 拷贝一次；None: 完整计算 |
| `compute(l, m)` | `vector<double>` | 任意 | 任意 | 每次 O(l·ng) |

---

## 构造器的工作

构造器做层次清晰的初始化，每层有不同的生命周期：

```
构造器
  ├── 1. 存 theta_/phi_ 的 span  ← 生命周期由调用者保证
  ├── 2. 预计算 trig 数组        ← 全生命周期不变
  │     ├── cos_theta_, sin_theta_
  │     └── cos_phi_, sin_phi_
  ├── 3. Full 模式 → setLMax(l)    ← 可后续通过 setLMax() 调整
```



---

## 存储复用

`reinit(theta, phi, l_max_resident)` 和 `reserveNg(max_ng)` 用于支持 `Hamiltonian::Callable` 中 Ylm 缓存在 k 点切换时的复用。

### reserveNg — 预分配每个 ylm 向量最大容量

切换 k 点时，GKK 数据长度可能变长超过 capacity，resize 触发新建空间和数据转移，所以应该有一个函数可以设置缓存向量的 capacity 。
与 setLMax 相比，它在 ng 维度扩容，并不承担数据构造职责。

```cpp
auto reserveNg(int max_ng) -> void;
```

将内部所有 vector 的 capacity 提升到 `max_ng`。不改变 `ng_` 和当前缓存值。

分模式执行：
- **Full 模式**：trig 数组 (`cos_theta_`, `sin_theta_`, `cos_phi_`, `sin_phi_`)、`y_lm_` 各条目、`top_column_Q_` 各条目的 `prev`/`curr` 全部分配到 `max_ng` 容量
- **None 模式**：仅 trig 数组分配（`y_lm_` 和 `top_column_Q_` 不存在）

调用 `std::vector::reserve`（`reserveNg` 内部），缩小请求无效果。

### reinit

实际代码等价于构造器。差别在于检查 cos/sin_theta/phi 的 capacity，如果不能容纳则 throw error；临时销毁 l_max_resident_（=-1） 以强制出发 setLMax 的构造程序。
因为使用过程中 theta，phi 与之前没有什么重合，所以权当是构造了一个新球谐函数缓存，只不过复用了之前的内存，避免销毁+重新构造的开销。


```cpp
auto reinit(theta, phi, l_max_resident) -> void;
```

**容量复用**：若之前调用了 `reserveNg(ng_max_)` 且 `ng_ <= ng_max_`，`resize` 不触发堆重分配。`setLMax()` 中对 `y_lm_` 和 `top_column_Q_` 的 `resize` 同理。

### 关键设计点

- **三角函数是全类共享的**：不论 Mode，不论访问路径，cos/sin 只算一次。`compute()` 直接复用这些数组，不需要传入 theta/phi 参数。
- **`theta_`/`phi_` 存的是 span**：类不拥有数据，调用者需保证生命周期覆盖 `RealSphericalHarmonics` 对象的使用期。
- **缓存构建**：构造器只在 Full 模式下调用 `setLMax()`；None 模式跳过，缓存数组保持空。

三层初始化对应三个不同的可变层级：trig 数组不可变、θ/φ 数据不可变但引用可失效、缓存可伸缩。

---

## 命名说明：`get_ang_grad_phi` / `get_ang_grad_theta`

`ang_grad` = angular gradient on the unit sphere（单位球面角梯度），即表面梯度 `∇_Ω`。
在球坐标中，完整 3D 梯度算子为 `∇ = (∂/∂r, (1/r)∂/∂θ, (1/(r sinθ))∂/∂φ)`。
球谐函数 Y_lm 只含角度部分，不含径向 `1/r` 因子，因此：
- `get_ang_grad_theta(l,m) = ∂Y_lm/∂θ`（球面梯度 θ 分量，r=1 时度量因子为 1）
- `get_ang_grad_phi(l,m)   = (1/sinθ)·∂Y_lm/∂φ`（球面梯度 φ 分量，含 sinθ 度量因子）

`1/q` 径向因子由调用方处理（记为 `β(q)/q`，q→0 时取 `β'(0)` 极限）。

## Q_lm_ang_grad_phi_ 缓存设计

`Q_lm_ang_grad_phi_` 是 `get_ang_grad_phi()` 的惰性缓存，存储 `R_l^m = Q_l^m / sinθ`，用于计算 `(1/sinθ)·∂Y_lm/∂φ`。

### 存储格式

```cpp
mutable std::optional<std::vector<std::vector<double>>> Q_lm_ang_grad_phi_;
```

- **`optional` 作为有效性标志**：`std::nullopt` 表示未初始化或已失效；engaged 状态表示缓存数据就绪。
- **索引方案**：与 `Q_lm_` 相同，`l*(l+1)/2 + m_abs`。
- **m=0 条目**：空 `vector`（默认构造，无堆分配），`fillAngGradPhiCache()` 跳过不填充。
- 没有单独的 `valid_` 布尔标志——`optional` 的 engaged/disengaged 状态就是有效性指示器。

### 生命周期

| 函数 | Q_lm_ang_grad_phi_ 操作 | 原因 |
|------|---------------------|------|
| `reinit()` | `.reset()` |  θ/φ 数据改变，缓存内容全部失效 |
| `reserveNg()` | 如果缓存已填充，对 inner vectors 做 `reserve()`；否则 no-op | 只影响 capacity，不影响数据 |
| `setLMax()` 扩展 | **no-op**（惰性） | 下次 `get_ang_grad_phi()` 发现 `optional` size 与 Q_lm_ 不匹配 → 自动填充 |
| `setLMax()` 收缩 | **outer vector resize 至新大小**（inner vectors 数据保留） | 与 Q_lm_ 设计一致——低 l/m 的 R_l^m 数据保留，避免惰性重建误判 size 不匹配而触发全量重建 |

`get_ang_grad_phi()` 首次调用或检测到 size 不匹配时，调用 `fillAngGradPhiCache()`：

```
fillAngGradPhiCache()
  ├── Q_lm_ang_grad_phi_.emplace(total_size)   ← 分配 outer vector
  ├── m_abs=0 的 inner vector 保持空（无堆分配）
  ├── for l = 1..l_max_resident_:
  │     for m_abs = 1..l:
  │       reserve(ng_capacity) + resize(ng_)  ← 只为 m≥1 预分配
  └── for m = 1..l_max_resident_:
        └── AngGradPhiRecurrence::seed → step1 → advance 填充整列
```

### AngGradPhiRecurrence 结构体

与 `QRecurrence` 接口相同，`step1`/`advance` 代码相同，仅有 `seed` 不同：

| | QRecurrence::seed | AngGradPhiRecurrence::seed |
|---|---|---|
| m=0 | `1/(2√π)` | `0`（结果始终为 0） |
| m=1 | `-1/(2√π) · √(3/2) · sinθ` | `-1/(2√π) · √(3/2)`（常数，θ=0 处有限） |
| m>=2 | `(-1)^m/(2√π) · ∏_{k=1}^{m} √((2k+1)/(2k)) · sin^m(θ)` | `(-1)^m/(2√π) · ∏_{k=1}^{m} √((2k+1)/(2k)) · sin^(m-1)(θ)`（θ=0 → 0） |

### 与纯 ∂Y/∂φ 的比较

| | ∂Y/∂φ (`get_dphi`，已移除) | (1/sinθ)·∂Y/∂φ (`get_ang_grad_phi`) |
|---|---|---|
| 需要的数据 | Q_l^m | R_l^m = Q_l^m / sinθ |
| 缓存 | 无（可直接用 Q_lm_） | 有（Q_lm_ang_grad_phi_，惰性分配） |
| 递推 | 无 | AngGradPhiRecurrence（seed 少一个 sinθ） |
| θ=0，m=0 | 0 | 0 |
| θ=0，|m|=1 | 0 | **有限非零** |
| θ=0，|m|>=2 | 0 | 0 |

---








---

## 内部实现

### Legendre 递推公式

所有访问路径共享同一组递推核心。定义归一化的连带 Legendre 函数（含自带的归一化系数）：

$$Q_l^m = N_l^m \cdot P_l^m(\cos\theta)$$

#### `seedPmm(m, curr)` — 种子值 $Q_m^m$

包含 Condon-Shortley 相位因子 $(-1)^m$：

$$Q_0^0 = \frac{1}{2\sqrt{\pi}}$$

$$Q_m^m = (-1)^m \cdot Q_0^0 \cdot \prod_{k=1}^{m} \sqrt{\frac{2k+1}{2k}} \cdot \sin^m(\theta)$$

#### `stepPmm1(m, prev, curr)` — $Q_{m+1}^m$ 第一步

$$Q_{m+1}^m = \sqrt{2m+3} \cdot \cos\theta \cdot Q_m^m$$

#### `advanceColumn(m, l, prev, curr)` — 向前三项递推

$$
r_1 = \sqrt{\frac{(2l-1)(2l+1)(l-m)}{l+m}}, \quad
r_2 = (l+m-1) \sqrt{\frac{(2l+1)(l-m)(l-m-1)}{(2l-3)(l+m)(l+m-1)}}
$$

$$
Q_l^m = \frac{\cos\theta \cdot r_1 \cdot Q_{l-1}^m - r_2 \cdot Q_{l-2}^m}{l-m}
$$

#### `retreatColumn(m, l, prev, curr)` — 反向递推（`setLMax` 缩小用）

系数 $c_1, c_2$ 与 `advanceColumn` 中的 $r_1, r_2$ 相同：

$$
Q_{l-2}^m = \frac{\cos\theta \cdot c_1 \cdot Q_{l-1}^m - (l-m) \cdot Q_l^m}{c_2}
$$

归一化系数 $N_l^m$ 已全部折叠进递推系数中，计算 $Q_l^m$ 时无需额外乘法。

### 缓存索引

`y_lm_` 是一维 `vector<vector<double>>`，按 `l*l + (m+l)` 平铺索引。`m = 0` 在每块中间，`±m` 对称分布两侧。

`assembleYlm` 将 `Q_l^m`（已含归一化系数和 Condon-Shortley 因子）乘以 $\sqrt{2} \cdot \{\cos(m\phi), \sin(m\phi)\}$ 写入缓存。

#### `top_column_Q_` — 增量递推的顶端缓存

`top_column_Q_` 按 `m` 索引（`m = 0..l_max_resident`），每个条目保存两个 vector：

| 成员 | 含义 |
|------|------|
| `curr` | $Q_{l_{\text{max}}}^{m}$ — 当前 `l_max_resident` 处的 Q 值 |
| `prev` | $Q_{l_{\text{max}}-1}^{m}$ — 上一行的 Q 值，用于继续递推 |

`top_column_Q_` 是 **Full Cache 模式专属**，None 模式下保持为空。它的作用是将 `setLMax()` 的扩展/收缩从 O(l_max · ng) 降为 O(Δl · l_max · ng)：

- **扩展**：从当前 `curr`/`prev` 继续向前递推 $\Delta l$ 步
- **收缩**：通过 `retreatColumn` 反向递推 $\Delta l$ 步恢复低 l_max 的顶端状态

当 `l_max_resident = m`（m 列只有种子值）时，`prev = curr`。
