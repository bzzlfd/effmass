# Hamiltonian 模块设计与一致性检查架构

## 接口设计

Hamiltonian 作为纯粹的数据综合，通过 `.at_k(ikpt)` 返回内嵌 functor 类 Hamiltonian::Callable，可作为参数传入其他函数。
callable functor 访问/使用 Hamiltonian 中的数据。

此外，还有两个 `Gradient`、`Hessian` 内嵌类，有别于 $H \ket{\psi}$，它们代表操作 $\frac{\partial H}{\partial k_\alpha} \ket{\psi}$、$\frac{\partial^2 H}{\partial k_{\alpha}\partial k_{\beta} } \ket{\psi}$，二者内部也有类似 `at_k` 成员函数返回内嵌 `Callable` 类作为 functor。

### H|ψ⟩ 计算接口

```cpp
auto op = h.at_k(ikpt);      // → Callable
op(psi, hpsi);                // H|ψ⟩

auto grad = h.gradient();
grad.at_k_a(ikpt, KDir::X);   // ∂H/∂k_x |ψ⟩

auto hess = h.hessian();
hess.at_k_ab(ikpt, KDir::X, KDir::Y); // ∂²H/∂k_x∂k_y |ψ⟩
```


### 数据加载
#### `.loadXxx()` — 逐步加载模式

`Hamiltonian` 采用**逐步加载**（incremental loading）模式，每个数据文件通过独立的 `loadXxx()` 方法加载：

```cpp
auto h = Hamiltonian("path/to/calc");
h.loadGKK("OUT.GKK");   // 加载波函数 G 矢量数据
h.loadWG("OUT.WG");     // 加载波函数系数
h.loadATOM("atom.config");
h.loadNCPPs("UPF");     // 加载赝势
// ...
```

**设计要点：**

- **路径解析**：相对路径基于构造时指定的 `base_dir` 解析，绝对路径原样使用。调用者无需要 `chdir` 即可加载不同目录的数据。
- **自动一致性检查**：每个 `loadXxx()` 调用末尾自动触发 `checkConsistency()`，验证新加载的数据与已有数据是否来自同一次计算。错误尽早暴露。
- **可选依赖**：部分文件是可选的（如 `VR`、`RHO`），未加载时对应的访问器会抛出异常，但不会影响其他功能。

#### `loadFromDirectory()` — 批量加载

```cpp
auto h = Hamiltonian("path/to/calc");
h.loadFromDirectory();
```

PWmat 的计算结果存放在一个文件夹下，所以也应该可以 load Directory。
按固定顺序加载标准文件，每步之后都触发一致性检查，让检查逐步覆盖更多交叉验证对：

1. `atom.config` → 设定晶格、原子数等基础信息
2. `UPF/`（或 `.`）→ 加载赝势，验证元素覆盖
3. `OUT.GKK` → 设定 k 网格、FFT 网格、Ecut 等核心参数
4. `OUT.WG` → 波函数系数，验证 nband 等
5. `OUT.VR` → 局域势
6. `OUT.RHO` → 电荷密度
7. `OUT.EIGEN` → 本征值，验证 k 点坐标
8. `OUT.OCC` → 占据数

#### 数据访问

```cpp
auto& gkk = h.gkk();    // 返回 const GKK&
auto& wg  = h.wg();
auto& occ = h.occ();
```

未加载时抛出 `std::runtime_error`。支持通过 `hasGKK()`、`hasWG()` 等查询是否存在。


#### checkConsistency 三部分架构

一致性检查分三部分，按计算成本和语义区分：

```
checkConsistency()  ← 每个 loadXxx() 后自动调用
├── Part 1  规范物理量（量级检查，O(1)）
├── Part 2  文件间完整性检查（pair 风格，O(1)）
└── checkConsistencyExtended()  ← 用户按需调用
    └── Part 3  高级计算性检查（FFT/积分，重量级）
```

##### Part 1: 规范物理量（canonical quantities）

**定义**：在多个文件中出现、且应该在 Hamiltonian 层面提供统一接口的物理量。

**判定原则**：如果一个变量出现在多个文件中，且 H|ψ⟩ 计算（Callable／gradient／hessian）需要从某处读取它，那么这个变量就应该进入 Part 1 作为规范值存储在 Hamiltonian 中。内部计算通过 `canonical_*` 访问，而非从 `gkk_`、`atom_` 等文件对象中临时获取。

**机制**：首个加载文件设定规范值，后续文件必须与之匹配。

```cpp
// 错误：直接从文件对象取用（如果它出现在多个文件中）
use(gkk_->meta.lattice);

// 正确：从 Hamiltonian 规范值取用（统一来源）
use(*canonical_lattice_);
```

**当前 Part 1 量集（9 个）：**

| 规范量 | 出现在 | 用途 |
|--------|--------|------|
| `lattice` | GKK, WG, VR, RHO, ATOM | 坐标变换、倒格子构造 |
| `nkpt` | GKK, WG, EIGEN, OCC | k 点循环维度 |
| `nband` | WG, EIGEN, OCC | 能带循环维度 |
| `n1,n2,n3` | GKK, WG, VR, RHO | FFT 网格 |
| `is_SO` | GKK, WG, EIGEN | 自旋轨道耦合分支 |
| `islda` | GKK, WG, EIGEN, OCC | 自旋极化分支 |
| `Ecut` | GKK, WG | 截断能 |
| `natom` | EIGEN, ATOM | 原子数 |
| `kpt_vec` | GKK, EIGEN, OCC | k 点坐标（分数坐标） |

**新量放置指南**：当一个物理量在 ≥ 2 个文件中出现，且 H|ψ⟩ 计算需要访问它时，应：
1. 在 `hamiltonian.cppm` 中添加 `std::optional<T> canonical_` 成员
2. 在 `checkPart1()` 中添加设定/比对逻辑
3. 内部代码通过 `canonical_*` 而非文件对象的元数据访问该量

##### Part 2: 文件间完整性检查

**定义**：纯验证文件间数据一致性的检查，不涉及计算语义。

**判定原则**：该量只用于"文件是否来自同一次计算"的验证，H|ψ⟩ 计算不需要直接访问它。这类检查保留 pair 风格，每次 `checkConsistency()` 全量重查（均为 O(1) 整数/浮点比较）。

**当前 Part 2 检查（4 组）：**

| 组 | 检查内容 | 涉及文件 | 说明 |
|----|---------|---------|------|
| 1 | `mg_nx`, `ng_tot_per_kpt` | GKK ↔ WG | 最大 G 矢量数、每 k 点 G 矢量总数 |
| 2 | `nnode` | GKK ↔ WG ↔ EIGEN | 并行节点数三方一致（VR/RHO 豁免） |
| 3 | `nstate` | VR ↔ RHO | 状态数 |
| 4 | 元素覆盖 | NCPP ↔ ATOM | 每类原子有对应 UPF |

**新量放置指南**：若一个量：
- 只出现在 1–2 个文件中（大概率），或
- 其唯一用途是验证文件一致性

则放入 Part 2。如果后续发现 H|ψ⟩ 计算需要它，再提升到 Part 1。

##### Part 3: 高级计算性检查

`checkConsistencyExtended()` — 用户按需调用，不自动触发。

**选择性执行**：通过 `ExtendedCheck` 枚举指定要跑的检查，内部用 `std::uint64_t` bitmask 去重：

```cpp
// 全部重检查（数据未加载的 check 自动跳过）
h.checkConsistencyExtended();

// 只跑子集
h.checkConsistencyExtended({
    ExtendedCheck::RHOReconstruct,
    ExtendedCheck::NCPPZVal,
});
```

**定义**：需要通过 FFT、积分、或其他复杂计算才能验证物理自洽性的检查。

**判定原则**：
- 计算量远超 O(1)，不适合每次 `loadXxx()` 后自动执行
- 验证的是物理自洽性而非文件一致性（例如"从 OCC+WG+GKK 重建的 RHO 是否与文件 RHO 一致"）
- 通常是端到端的数据链路验证

**Part 3 检查清单：**

| 检查 | 依赖数据 | 算法 | 容差 |
|------|---------|------|------|
| RHO 重建 | GKK+WG+OCC+RHO | OCC × |WG|² → FFT → 实空间 → vs 文件 RHO | RMSE < 1e-2 |
| 价电子数 | ATOM+RHO+NCPPs | Σ(NCPP.z_val × count) ≈ ∫RHO d³r | diff ≤ 1.0 |
| NCPP 价电子合理性 | NCPPs（及 ATOM 元素名） | 0 < z_val ≤ atomic_number | — |

**实现模块**：重量级计算逻辑（`reconstructRHO`、`integrateRHO`、`compareRHO`）封装在 `H_psi.density_reconstruction` 模块中，可独立导入复用：

##### 三部分判定速查表

| 新检查的特征 | 放置 |
|-------------|------|
| 同一个量出现在不同文件 | **Part 2**（ pair 检查） |
| 量出现在 ≥ 2 个文件；H|ψ⟩ 等计算需要访问 | **Part 1**|
| 验证涉及复杂计算 | **Part 3**（单独方法） |

## 设计要点总结

- **逐步加载 + 自动检查**：尽早暴露数据不一致问题，且不强制所有文件同时加载
- **规范量集中存储**：避免 H|ψ⟩ 计算在多个文件对象间查找所需数据
- **三部分分离**：轻量级自动检查与重量级计算验证解耦
- **可扩展**：新数据类型只需在对应 Part 的检查块中增加条目，无需新增 N 个 pair
