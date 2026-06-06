# RealSphericalHarmonics 接口设计

## 概述

`RealSphericalHarmonics` 是向量化的实球谐函数求值器，对一组 K 点的 G 矢量（theta, phi）同时计算 $Y_{lm}$。类位于模块 `math.sph_harmonics`，是 `math` 模块的一部分。

核心设计思想：**两种缓存策略 × 三种访问方法**，调用者根据自己的访问模式选择合适路径。

---

## 1. None Cache：compute 模式

```
构造 → 只存 trig 数组 → compute(l, m) 是唯一访问路径
```

None 模式下，构造器只做一件事：**预计算三角函数**（cosθ, sinθ, cosφ, sinφ）。`reset()` 不会被调用，`y_lm_` 和 `top_column_Q_` 保持空。

访问完全依赖 `compute(l, m)`：

```
compute(l, m)
  ├── seedPmm(m, q_curr)       ← Q_m^m = N_m^m · sin^m(θ) · 1/(2√π)
  ├── stepPmm1 / advanceColumn  ← 递推到 Q_l^m
  └── 组装: Q_l^m → Y_lm       ← 乘以 √2 · {cos(mφ), sin(mφ)}
```

每次都是完整计算，O(l · ng)，结果不缓存。适合访问稀疏、l 上限不可预测、或内存受限的场景。

---

## 2. Full Cache：reset → get

这是设计的核心，分为冷热两条路径。

### reset() + get() — 热路径

```
reset(l_max)
  └── 逐列递推: m=0..l_max
        ├── seedPmm  →  Q_m^m
        ├── stepPmm1 →  Q_{m+1}^m
        ├── advanceColumn → Q_{m+2..l_max}^m
        └── assembleYlm → 写入 y_lm_[l*l + (m+l)]

get(l, m)
  └── y_lm_[l*l + (m+l)]   ← O(1)，零拷贝
```

`reset()` 完成后，所有 `l ≤ l_max_resident` 的 Y_lm 在 `y_lm_` 数组中立即可用。`get()` 就是一次数组索引。

`reset()` 的增量设计：**保存 `top_column_Q_`**（每个 m 在 `l_max_resident` 处的 Legendre 递推顶端状态）。后续扩展不需要从头再来，只需从当前顶端继续递推，复杂度 O(Δl · l_max · ng)。缩小同样通过 `retreatColumn` 反向递推恢复状态。

> [!NOTE] **`top_column_Q_` 是 Full Cache 的专属状态**：None Cache 模式下保持为空，只有 Full Cache 的 `reset()` 才会写入和更新。
> 它是一个按 m 索引的结构（m = 0..l_max_resident）——对每个 m 保存 `l = l_max_resident` 处的 Q 矢量（大小 ng），记录递推中断处的顶端行。后续 `reset(l)` 扩展或收缩时，算法从这些顶端值继续递推，而非从头重算。

### compute() — 冷路径 (bonus of NoneCache)

当 `l > l_max_resident_`，`get()` 不工作。调用者有两个选择：

| 路径 | 开销 | 适用场景 |
|------|------|---------|
| `reset(l); get(l, m)` | 一次性 O(Δl·l_max·ng)，后续 O(1) | 后续会反复访问高 l |
| `compute(l, m)` | 每次 O(l·ng)，不占用缓存 | 偶发一次查询，不确定是否会复用 |

`compute()` 在 Full 模式下的角色是：**缓存的补充机制**，覆盖缓存范围之外的即时查询。它不是主路径，而是处理边界情况的逃生口。

---

## 3. operator() — 统一接口

```
operator()(l, m)
  ├── 校验量子数 (l, m)
  ├── Full 模式
  │     ├── l > l_max_resident_ → 报错，提示 reset 或 compute
  │     └── 正常 → return get(l, m)   ← 拷贝一次
  └── None 模式 → return compute(l, m)
```

定位：**便利入口**，调用者可以不关心 CacheMode 直接用 `sh(l, m)` 拿到 Y_lm。代价是返回值而非引用（Full 模式下多做一次拷贝），性能敏感场合应直接调用 `get()` 或 `compute()`。

`operator()` 独自做校验而非委托给 `get()`，保证了报错信息中的函数名与实际调用一致。

### 三种方法对比

| 方法 | 返回类型 | CacheMode | l 范围 | 开销特征 |
|------|---------|-----------|--------|---------|
| `get(l, m)` | `const vector<double>&` | Full 必须 | ≤ l_max_resident | O(1)，零拷贝 |
| `operator()(l, m)` | `vector<double>` | 任意 | 任意 | Full: 拷贝一次；None: 完整计算 |
| `compute(l, m)` | `vector<double>` | 任意 | 任意 | 每次 O(l·ng) |

---

## 4. 构造器的工作

构造器做层次清晰的初始化，每层有不同的生命周期：

```
构造器
  ├── 1. 存 theta_/phi_ 的 span  ← 生命周期由调用者保证
  ├── 2. 预计算 trig 数组        ← 全生命周期不变
  │     ├── cos_theta_, sin_theta_
  │     └── cos_phi_, sin_phi_
  └── 3. Full 模式 → reset(l)    ← 可后续通过 reset() 调整
```

### 关键设计点

- **三角函数是全类共享的**：不论 Mode，不论访问路径，cos/sin 只算一次。`compute()` 直接复用这些数组，不需要传入 theta/phi 参数。
- **`theta_`/`phi_` 存的是 span**：类不拥有数据，调用者需保证生命周期覆盖 `RealSphericalHarmonics` 对象的使用期。
- **缓存构建是延迟的**：构造器只在 Full 模式下调用 `reset()`；None 模式跳过，缓存数组保持空。

三层初始化对应三个不同的可变层级：trig 数组不可变、θ/φ 数据不可变但引用可失效、缓存可伸缩。

---

## 5. 内部实现机制

### Legendre 递推

所有访问路径共享同一组递推核心：

- `seedPmm(m, curr)` — 种子值 Q_m^m
- `stepPmm1(m, prev, curr)` — Q_{m+1}^m 第一步
- `advanceColumn(m, l, prev, curr)` — 向前三项递推
- `retreatColumn(m, l, prev, curr)` — 反向递推（reset 缩小用）

归一化系数已折叠进递推系数中，计算 Q_l^m 时无需额外乘法。

### 缓存索引

`y_lm_` 是一维 `vector<vector<double>>`，按 `l*l + (m+l)` 平铺索引。`m = 0` 在每块中间，`±m` 对称分布两侧。

`assembleYlm` 将 Q_l^m（含归一化）乘以 √2 · {cos(mφ), sin(mφ)} 写入缓存。
