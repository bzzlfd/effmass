# IO 模块设计细节与接口

本文档介绍 `io` 模块及其子模块 `io.GKK`、`io.WG` 的设计思想与对外接口。关于文件二进制格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md)。

## 模块结构

```
io (aggregate)
├── io.GKK     →  src/io/GKK.cppm
├── io.WG      →  src/io/WG.cppm
└── io.lattice →  src/io/lattice.cppm
```

`io.cppm` 是聚合模块（aggregate module），通过 `export import io.GKK; export import io.WG; export import io.lattice;` 重新导出三个子模块的所有公开接口。外部代码只需 `import io;` 即可使用全部 IO 功能。

## 数据结构

### `GKKMetadata`

存储 `OUT.GKK` 文件的元数据：

```cpp
export struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;
    double Ecut;              // cutoff energy (Hartree)
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};
```

### `KVecsView`

控制 `KVecs` 中哪些表示被计算和暴露的位掩码枚举：

```cpp
export enum class KVecsView : unsigned int {
    Cartesian = 1 << 0,  // kinetic, Kx, Ky, Kz
    Spherical = 1 << 1,  // r, theta, phi
    Integer   = 1 << 2,  // iG, jG, kG, kPoint, reciprocalLattice
};
```

支持位运算 `|`, `&`, `~` 以及 `hasView(flags, view)` 辅助函数，可自由组合。例如：

```cpp
auto view = KVecsView::Cartesian | KVecsView::Spherical;
```

### `KVecs`

某个 k 点的 G 向量数据视图，通过指针指向类内部预分配的连续缓冲区。各字段是否非空取决于 `GKK::currentView()` 和已就绪的表示：

```cpp
export struct KVecs {
    // Cartesian representation: K = G - k
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, G_x-k_x, G_y-k_y, G_z-k_z
    
    // Spherical representation of K = G - k
    std::span<const double> r, theta, phi;        // |K|, polar angle [0,π], azimuthal angle [-π,π]

    // Integer indices of G vector in reciprocal lattice basis
    std::span<const int>    iG, jG, kG;           // G = iG*b1 + jG*b2 + kG*b3
    // Per-k-point metadata (valid whenever Integer view is enabled)
    std::array<double, 3>               kPoint{};            // fractional coordinate of current k-point
    std::array<std::array<double, 3>, 3> reciprocalLattice{}; // reciprocal lattice vectors: row n = b_n
};
```

- **Cartesian**（始终可加载）：`kinetic` 为 `|G+k|²/2`，`Kx/Ky/Kz` 为 `G-k` 的 Cartesian 分量。
- **Spherical**：由 `Kx/Ky/Kz` 转换而来，`r = |K|`，`theta = acos(Kz/r)`，`phi = atan2(Ky, Kx)`。
- **Integer**：利用 `inferCurrent_k()` 推断 k 点分数坐标，通过公式 `(iG,jG,kG) = round(A^T * K_cart / (2π) + k_frac)` 计算 G 向量的整数米勒指数。同时填充该 k 点相关的全局元数据：
  - `kPoint`：当前 k 点的分数坐标（由 `inferCurrent_k()` 推断，范围 `(-0.5, 0.5]`）
  - `reciprocalLattice`：倒格子矩阵（行向量 `b1, b2, b3`，单位 Bohr⁻¹），由 `Lattice::B()` 在每次调用时动态计算。

### `Lattice`

晶格（实空间 + 倒空间）的物理量封装，由 `io.lattice` 子模块导出。

```cpp
export class Lattice {
public:
    Lattice() = default;  // 默认构造，A_ 为全零
    explicit Lattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit);
    explicit Lattice(std::span<const double, 9> flat, LengthUnit unit);
    Lattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit);

    auto setLattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit) -> void;
    auto setLattice(std::span<const double, 9> flat, LengthUnit unit) -> void;
    auto setLattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit) -> void;

    auto A(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;
    auto B(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;
};
```

#### 单位

`Lattice` 内部统一采用 **Hartree 原子单位制** 存储：

通过构造或 `setLattice()` 传入 `LengthUnit`（`Bohr` 或 `Angstrom`），转换成 Bohr 单位。
读取时（`A()`, `B()`）再由 Bohr 单位转换。

#### 构造函数与 `setLattice`

-  `Lattice() = default;` 
-  `Lattice(span<const double, 9>, LengthUnit)` 
-  `Lattice(const array<array>& A, LengthUnit)` 
-  `Lattice(std::initializer_list<initializer_list> A, LengthUnit)` 
-  `setLattice(span<const double, 9>, LengthUnit)` 
-  `setLattice(const array<array>& A, LengthUnit)` 
-  `setLattice(std::initializer_list<initializer_list> A, LengthUnit)` 

转发路径

```
Lattice(span) - setLattice(span) - setLatticeFromFlat
Lattice(array) - setLattice(array) - setLatticeFromFlat
Lattice(ilist) - setLattice(ilist) - setLatticeFromFlat
```

`setLatticeFromFlat` 完成：将 flat 数据写入 `A_` → 单位转换（Angstrom→Bohr）→ `computeReciprocalLattice()` 提前验证体积非零。

#### 访问接口

- `A(LengthUnit unit = LengthUnit::Bohr)`：返回实空间晶格矩阵，可选指定输出单位。
- `B(LengthUnit unit = LengthUnit::Bohr)`：返回倒空间晶格矩阵，计算后可选指定输出单位。

`A()` 和 `B()` 均采用按值返回，因为当请求的单位与内部存储单位不同时需要进行数值转换。

内部数据存储, `A()`/`B()` 均由计算得到，其数据不由 lattice 存储。

## `GKK` 类

`GKK` 类封装了对 `OUT.GKK` 文件的读取操作。

### 构造函数与生命周期

```cpp
explicit GKK(const std::string& filename);
~GKK();
```

- 构造函数打开文件并读取元数据
- 析构函数关闭文件句柄

### 移动语义

```cpp
GKK(GKK&& other) noexcept;
GKK& operator=(GKK&& other) noexcept;
```

支持移动构造和移动赋值。移动后会自动更新 `KVecs` 内部的指针，使其指向新的缓冲区。

> 拷贝被显式删除，因为类持有 `FILE*` 句柄和大量缓冲区。

### 公共接口

```cpp
const GKKMetadata& metadata() const;                          // 获取元数据
void setDataView(KVecsView view);                             // 设置需要计算的数据表示
KVecsView currentView() const;                                // 查询当前设置的数据表示
const KVecs& loadKPoint(int ikpt);                            // 加载指定 k 点数据（带缓存）
int current_ikpt() const;                                     // 当前缓存的 k 点索引
const KVecs& currentData() const;                             // 当前缓存的数据视图（按 currentView 过滤）
std::array<double, 3> inferCurrent_k() const;                 // 从 K=G-k 数据推断 k 点分数坐标 (fractional coordinate)
```

- `setDataView`
    1. 设置 `desired_views_`，即 `GKK` 对象返回的 `KVecs` 中有哪些数据表示类型（笛卡尔坐标，球坐标，整数索引）。
    2. 分配/销毁对应数据表示的存储空间
- `loadKPoint`
    1. 从文件中加载第 `ikpt` 个 k 点的 K 向量。
    2. 按照 `desired_views_`，按需计算其他数据表示。


## 设计要点

### 延迟加载

构造函数**不**加载全部 k 点数据，只读取元数据并预计算每个 k 点在文件中的偏移量。实际数据按需通过 `loadKPoint` 加载。

#### 偏移预计算

`computeOffsets()` 在构造时遍历所有 record，记录每个 k 点的起始文件偏移。这样 `loadKPoint` 可以直接 `fseek` 到目标位置，而不需要顺序跳过前面的数据。

#### 缓冲区复用/惰性分配

类在构造时根据 `mg_nx * nnodes` 预分配 **Cartesian** 必需的缓冲区：
- `kinetic_buf_`, `Kx_buf_`, `Ky_buf_`, `Kz_buf_`

**Spherical** 和 **Integer** 的缓冲区（`r_buf_`、`theta_buf_`、`phi_buf_`、`iG_buf_`、`jG_buf_`、`kG_buf_`）采用**惰性分配**：仅在 `setDataView` 启用对应视图时分配，禁用时通过 `clear()` + `shrink_to_fit()` 释放堆内存。


## 错误处理

- **文件打开失败**：抛出 `std::runtime_error("Cannot open file: ...")`
- **Record 长度不匹配**：`checkRecordLength()` 抛出 `std::runtime_error("Record length mismatch")`
- **Record 大小与预期不符**：`readRecord()` 抛出 `std::runtime_error("context: record size mismatch")`
- **读取失败**：`readRecord()` 抛出 `std::runtime_error("context: read failed")`
- **k 点索引越界**：`loadKPoint()` / `seekToKPoint()` 抛出 `std::out_of_range`
- **元数据不一致**：`nnodes mismatch` 等异常

所有异常信息都包含上下文描述（如 `header`、`Ecut`、`gkk` 等），便于定位出错的 record。

## `WG` 类

`WG` 类封装了对 `OUT.WG` 波函数系数文件的读取操作。其设计思想与 `GKK` 类基本一致，重复部分不再赘述，仅记录差异点。

### 数据结构

#### `WGMetadata`

与 `GKKMetadata` 结构相同，但多一个 `mx` 字段表示能带数：

```cpp
export struct WGMetadata {
    int n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda;
    double Ecut;
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
    std::vector<int> ng_tot_per_kpt;
};
```

#### `WGCoeffs`

存储单个 k 点、单个能带的波函数系数视图：

```cpp
export struct WGCoeffs {
    std::span<const std::complex<double>> up;     // 自旋向上分量
    std::span<const std::complex<double>> down;   // 仅在 is_SO == 1 时有效
};
```

### 公共接口

```cpp
const WGMetadata& metadata() const;
const WGCoeffs& loadBand(int ikpt, int iband);  // 加载指定 k 点、指定能带（带缓存）
int current_ikpt() const;
int current_iband() const;
const WGCoeffs& currentData() const;
```

### 与 `GKK` 的差异

1. **偏移预计算维度不同**：
   `OUT.WG` 的数据按 **k 点 → 节点 → 能带** 的三层嵌套顺序存储。因此 `computeOffsets()` 预计算的是每个 `(ikpt, iband)` 组合的文件偏移，而非仅每个 k 点。

2. **自旋轨道耦合的数据布局**：
   当 `is_SO == 1` 时，一个 Fortran record 中连续存放了向上和向下两个 `complex*16` 数组。`loadBand` 读取后将其拆分：用成员 `tmp_buf_` 作为临时缓冲区读入完整 record，再分别拷贝到缓存条目的 `up` 和 `down` 向量中。

3. **多 band 缓存**：
   `WG` 支持缓存**多个** `(k-point, band)` 组合，数量由构造函数参数 `cache_capacity` 控制（默认 `4`）。

   缓存实现要点：
   - **独立存储**：每个缓存条目（`CacheEntry`）拥有独立的 `std::vector<std::complex<double>>`，存储该 band 完整的 up/down 系数数据，以及一个 `WGCoeffs` 视图（`std::span`）指向自身数据。
   - **地址稳定**：使用 `std::vector<std::unique_ptr<CacheEntry>>` 管理条目。`unique_ptr` 保证 `CacheEntry` 对象地址在容器重新分配时不发生变化，因此 `WGCoeffs` 内部的 `std::span` 始终有效。
   - **替换策略**：采用 **FIFO 环形替换**。当缓存已满且未命中时，复用 `cache_next_slot_` 指向的旧槽位，然后该索引循环前进。
   - **缓存命中**：`loadBand` 先在缓存中线性查找；命中后直接返回对应条目的视图，避免文件 IO。

   这意味着遍历同一 k 点的多个能带时，只要在缓存容量范围内，就不会重复读取磁盘。返回的 `const WGCoeffs&` 只要该 band 未被逐出缓存就保持有效。

4. **移动语义简化**：
   由于每个 `CacheEntry` 通过 `unique_ptr` 拥有独立堆内存，移动 `WG` 时各条目地址不变，内部 `std::span` 无需像 `GKK` 那样手动修复指针。移动构造/赋值只需直接转移 `unique_ptr` 所有权即可。

### 元数据一致性

`OUT.WG` 的头部元数据（`n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda, Ecut, lattice, ngtotnod`）必须与 `OUT.GKK` 保持一致（仅多出 `mx`）。当前类内部**暂不做交叉校验**，由调用方负责验证。
