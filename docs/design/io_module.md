# IO 模块设计细节与接口

本文档介绍 `io` 模块及其子模块 `io.GKK`、`io.WG` 的设计思想与对外接口。关于文件二进制格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md)。

## 模块结构

```
io (aggregate)
├── io.ATOM    →  src/io/ATOM.cppm
├── io.EIGEN   →  src/io/EIGEN.cppm
├── io.GKK     →  src/io/GKK.cppm
├── io.lattice →  src/io/lattice.cppm
├── io.OCC     →  src/io/OCC.cppm
├── io.RHO     →  src/io/RHO.cppm
├── io.VR      →  src/io/VR.cppm
└── io.WG      →  src/io/WG.cppm
```

`io.cppm` 是聚合模块（aggregate module），通过 `export import` 重新导出所有子模块的公开接口。外部代码只需 `import io;` 即可使用全部 IO 功能。

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
    Spherical = 1 << 1,  // q, theta, phi
    Integer   = 1 << 2,  // g_idx, kPoint, reciprocalLattice
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
    // Cartesian: semantically interpreted as -K = -(G+k)
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, -(G_x+k_x), -(G_y+k_y), -(G_z+k_z)
    
    // Spherical representation of -K
    std::span<const double> q, theta, phi;        // |K|, polar angle [0,π], azimuthal angle [-π,π]

    // Integer indices of G vector in reciprocal lattice basis: G = g_idx.x*b1 + g_idx.y*b2 + g_idx.z*b3
    std::span<const vector3d<int>>  g_idx;
    // Per-k-point metadata (valid whenever Integer view is enabled)
    kVec                                kPoint{};            // fractional coordinate of current k-point
    array2d<double, 3, 3>               reciprocalLattice{}; // reciprocal lattice vectors: row n = b_n
};
```

- **Cartesian**（始终可加载）：`kinetic` 为 $|G+k|^2/2$，`Kx/Ky/Kz` 语义解释为 $-(G+k)$ 的 Cartesian 分量。
- **Spherical**：由 `Kx/Ky/Kz` 转换而来，`q = |K|`，`theta = acos(Kz/q)`，`phi = atan2(Ky, Kx)`。
- **Integer**：使用点取反恢复 $G+k$ 后，通过公式 `(g_idx.x, g_idx.y, g_idx.z) = round(A^T·(G+k) / (2π) - k_frac)` 计算 G 向量的整数米勒指数。同时填充该 k 点相关的全局元数据：
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
GKKMetadata meta;                                                // 元数据（公有成员）
void setDataView(KVecsView view);                                // 设置需要计算的数据表示
KVecsView currentView() const;                                   // 查询当前设置的数据表示
const KVecs& loadKPoint(int ikpt);                               // 加载指定 k 点数据（带缓存）
int current_ikpt() const;                                        // 当前缓存的 k 点索引
const KVecs& currentData() const;                                // 当前缓存的数据视图（按 currentView 过滤）
std::array<double, 3> inferCurrent_k() const;                    // 从 -K = -(G+k) 数据推断 k 点分数坐标
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
- `kinetic_`, `Kx_`, `Ky_`, `Kz_`

**Spherical** 和 **Integer** 的缓冲区（`q_`、`theta_`、`phi_`、`g_idx_`）采用**惰性分配**：仅在 `setDataView` 启用对应视图时分配，禁用时通过 `clear()` + `shrink_to_fit()` 释放堆内存。


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

#### `WGCoeffsView`

存储单个 k 点、单个能带的波函数系数视图（非 owning，`std::span`）：

```cpp
export struct WGCoeffsView {
    std::span<const std::complex<double>> up;     // 自旋向上分量
    std::span<const std::complex<double>> down;   // 仅在 is_SO == 1 时有效
};
```

### 公共接口

```cpp
WGMetadata meta;                                          // 元数据（公有成员）
const WGCoeffsView& loadBand(int ikpt, int iband);            // 加载指定 k 点、指定能带（带缓存）
int current_ikpt() const;
int current_iband() const;
const WGCoeffsView& currentData() const;
```

### 与 `GKK` 的差异

1. **偏移预计算维度不同**：
   `OUT.WG` 的数据按 **k 点 → 节点 → 能带** 的三层嵌套顺序存储。因此 `computeOffsets()` 预计算的是每个 `(ikpt, iband)` 组合的文件偏移，而非仅每个 k 点。

2. **自旋轨道耦合的数据布局**：
   当 `is_SO == 1` 时，一个 Fortran record 中连续存放了向上和向下两个 `complex*16` 数组。`loadBand` 读取后将其拆分：用成员 `tmp_buf_` 作为临时缓冲区读入完整 record，再分别拷贝到缓存条目的 `up` 和 `down` 向量中。

3. **多 band 缓存**：
   `WG` 支持缓存**多个** `(k-point, band)` 组合，数量由构造函数参数 `cache_capacity` 控制（默认 `4`）。

   缓存实现要点：
   - **独立存储**：每个缓存条目（`CachedState`）拥有独立的 `std::vector<std::complex<double>>`，存储该 band 完整的 up/down 系数数据，以及一个 `WGCoeffsView` 视图（`std::span`）指向自身数据。
   - **地址稳定**：使用 `std::vector<std::unique_ptr<CachedState>>` 管理条目。`unique_ptr` 保证 `CachedState` 对象地址在容器重新分配时不发生变化，因此 `WGCoeffsView` 内部的 `std::span` 始终有效。
   - **替换策略**：采用 **FIFO 环形替换**。当缓存已满且未命中时，复用 `cache_next_slot_` 指向的旧槽位，然后该索引循环前进。
   - **缓存命中**：`loadBand` 先在缓存中线性查找；命中后直接返回对应条目的视图，避免文件 IO。

   这意味着遍历同一 k 点的多个能带时，只要在缓存容量范围内，就不会重复读取磁盘。返回的 `const WGCoeffsView&` 只要该 band 未被逐出缓存就保持有效。

4. **移动语义简化**：
   由于每个 `CachedState` 通过 `unique_ptr` 拥有独立堆内存，移动 `WG` 时各条目地址不变，内部 `std::span` 无需像 `GKK` 那样手动修复指针。移动构造/赋值只需直接转移 `unique_ptr` 所有权即可。

### 元数据一致性

`OUT.WG` 的头部元数据（`n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda, Ecut, lattice, ngtotnod`）必须与 `OUT.GKK` 保持一致（仅多出 `mx`）。当前类内部**暂不做交叉校验**，由调用方负责验证。

## `EIGEN` 类

`EIGEN` 类封装了对 `OUT.EIGEN` 文件的读取操作。

### 数据结构

#### `EIGENMetadata`

```cpp
export struct EIGENMetadata {
    int islda;        // 自旋极化标志
    int nkpt;         // k 点数
    int nband;        // 能带数（Fortran 原变量名 mx）
    int nref_tot_8;   // （含义尚不明确，原样保留）
    int natom;        // 原子数
    int nnode;        // 节点数（Fortran 原变量名 nnodes）
    int is_SO;        // 自旋轨道耦合标志（旧格式文件可能缺失，此时设为 0）
};
```

#### `kVec`

k 点分数坐标向量类型，io 各子模块内部统一使用 `kVec`（`using kVec = vector3d<double>;`），但**不 export**：

```cpp
struct vector3d<double> {
    double x{}, y{}, z{};
    auto operator[](int i) -> double&;           // 0→x, 1→y, 2→z
    auto data() -> double*;                      // returns &x
    auto norm_squared() const -> double;
    auto norm() const -> double;
    // operator+=, -=, *=, /=, +, -, *, /, dot, cross, ==, !=
};
```

外部代码通过 `vector3d<double>` 或成员名（`kpt_vec[ik].x`）访问，`kVec` 仅为模块内部别名。

#### `Array3d`

用于存储能量本征值的三维数组封装，一维连续存储，顺序 `[ispin][ikpt][iband]`：

```cpp
export class Array3d {
    // 构造
    Array3d(int nband, int nkpt, int islda);

    // 查询
    auto dims() const -> std::array<int, 3>;  // {nband, nkpt, islda}
    auto empty() const -> bool;
    auto size() const -> std::size_t;  // 总元素数
    auto data() -> double*;
    auto data() const -> const double*;

    // 多参数直接访问
    auto operator[](int iband, int ikpt, int ispin) -> double&;  // 3 参数
    auto operator[](int iband, int ikpt, int ispin) const -> double;

    auto operator[](int iband, int ikpt) -> double&;  // 2 参数，仅 islda==1 时可用
    auto operator[](int iband, int ikpt) const -> double;

    // 链式访问 [ispin][ikpt][iband]
    auto operator[](int ispin)       -> Slice2d;
    auto operator[](int ispin) const -> ConstSlice2d;
    // Slice2d[ikpt] 返回 std::span<double> (len = nband)
    // ConstSlice2d[ikpt] 返回 std::span<const double> (len = nband)
};
```

### `EIGEN` 类

```cpp
export class EIGEN {
public:
    explicit EIGEN(const std::string& filename);
    ~EIGEN();

    // 禁用拷贝，启用移动
    EIGEN(const EIGEN&) = delete;
    auto operator=(const EIGEN&) -> EIGEN& = delete;
    EIGEN(EIGEN&&) noexcept;
    auto operator=(EIGEN&&) noexcept -> EIGEN&;

    // 打印元数据 + sum(kpt_weight)
    auto print_info() const -> void;

    // 公有数据成员
    EIGENMetadata meta;
    std::vector<double> kpt_weight;    // k 点权重，size = nkpt
    std::vector<vector3d<double>> kpt_vec;    // k 点坐标，size = nkpt
    Array3d eigenvalue;                // 本征值，dims = [nband, nkpt, islda]
};
```

### 设计要点

1. **一次性读取**：EIGEN 文件通常较小（数 KB），构造函数中读取所有元数据和数据，不采用延迟加载策略。文件读取完毕立即关闭，不持有成员 `FILE*`。Rule of Zero。

2. **文件格式兼容性**：header record 通过前导长度标记区分新/旧格式：
   - 28 字节 → 7 int（含 `is_SO`）
   - 24 字节 → 6 int（无 `is_SO`，设为 0）
   - 其他长度 → 抛异常

3. **验证**：读取 k 点循环时验证 `islda_tmp` 和 `ikpt_tmp` 是否与循环计数器一致（Fortran 1-based），不一致时抛异常。

4. **数据冗余**：文件中每个 `(ispin, ikpt)` 对都有一份 weight/ak 数据（4 doubles）。这些量在语义上按 kpt 索引（而非按自旋），因此最终保留最后一次读入的值。

5. **移动语义**：Rule of Zero。无 `FILE*` 成员，默认移动/拷贝即可满足需求。

## `OCC` 类

`OCC` 类封装了对 `OUT.OCC` 文本格式文件的读取操作。关于文件格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md#outocc)。

### `OCCMetadata`

```cpp
export struct OCCMetadata {
    int islda{};      // 自旋数（1 或 2），自动检测
    int nkpt{};       // k 点数
    int nband{};      // 每个 k 点的能带数
};
```

### `OCC` 类

```cpp
export class OCC {
public:
    explicit OCC(const std::string& filename);

    // Rule of zero
    OCC() = default;
    OCC(const OCC&) = default;
    auto operator=(const OCC&) -> OCC& = default;
    OCC(OCC&&) noexcept = default;
    auto operator=(OCC&&) noexcept -> OCC& = default;

    auto energy(int iband, int ikpt) const -> double;              // islda=1
    auto energy(int iband, int ikpt, int ispin) const -> double;   // islda=2
    auto occupation(int iband, int ikpt) const -> double;          // islda=1
    auto occupation(int iband, int ikpt, int ispin) const -> double; // islda=2
    auto print_info() const -> void;

    OCCMetadata meta;
    std::vector<kVec> kpt_vec;          // k 点分数坐标，size = nkpt

private:
    std::vector<double> energies_;       // eV，平铺 [ispin][ikpt][iband]
    std::vector<double> occupations_;    // 平铺 [ispin][ikpt][iband]
};
```

参数顺序 `(iband, ikpt)` 与 `EIGEN` 一致。

### 设计要点

1. **文本解析**：与 `ATOM` 类一致，使用 `fopen("r")` + `fgets` 逐行读取，`strtol`/`strtod` 链式解析数值。与其它二进制 IO 类（GKK/WG/EIGEN）的最显著区别。

2. **Spin=2 自动检测**：读取首行，若以 `=` 开头则进入 spin 模式，按 `==========  SPIN N  ==========` 分隔两段数据；否则视为单自旋（islda=1）。第二自旋的 k 点坐标自动与第一自旋校验一致。

3. **一次性读取**：OCC 文件通常较小（数百行），构造函数中读取所有元数据和数据。不采用延迟加载策略。

4. **局部 FILE***：构造函数内使用局部 FILE*，解析完成后立即关闭。try/catch 保护解析过程。Rule of Zero。

5. **数据验证**：
   - 每个 k 点的能带数必须一致（跨 k 点验证）
   - 总数据条目必须等于 islda × nkpt × nband
   - k 点索引必须为 1-based 递增正整数

6. **占据数语义**：文件中存储的占据数已包含 k 点权重，即 $\sum_{\text{band}} f_{n\mathbf{k}} = N_{\text{valence}} \cdot w_{\mathbf{k}}$。

7. **与 OUT.EIGEN 的交叉验证**：`OUT.OCC` 与 `OUT.EIGEN` 共享相同的 k 点网格和能带结构。测试代码通过对比两个文件的能量本征值（OCC 四舍五入到 4 位小数，与 EIGEN 原始值差距 ≤ 0.001 eV），以及 k 点坐标（≤ 1e-4 误差），来验证文件读取的正确性。

8. **内部类型别名**：OCC 内部使用 `using kVec = vector3d<double>;` 表示 k 点向量（不 export），与 `EIGEN`、`GKK` 保持一致。

## `RHO` / `VR` 类

`RHO` 类封装了对 `OUT.RHO`（电荷密度）和 `OUT.VR`（价电子势能）实空间 FFT 网格数据的读取操作。两者文件格式完全相同，`VR` 是 `RHO` 的类型别名（`using VR = RHO`）。

### 数据结构

#### `RHOMetadata`

```cpp
export struct RHOMetadata {
    int n1, n2, n3;       // FFT 网格维度
    int nnodes;           // 节点数
    int nstate;           // （含义尚不明确，原样保留）
};
```

### `RHO` 类

```cpp
export class RHO {
public:
    explicit RHO(const std::string& filename);
    ~RHO();

    // 禁用拷贝，启用移动
    RHO(const RHO&) = delete;
    auto operator=(const RHO&) -> RHO& = delete;
    RHO(RHO&&) noexcept;
    auto operator=(RHO&&) noexcept -> RHO&;

    // 公有数据成员
    RHOMetadata meta;
    Lattice lattice;

    // 三维网格访问
    auto operator[](int i, int j, int k) -> double&;              // state 0，无运行时检查
    auto operator[](int i, int j, int k) const -> double;
    auto operator[](int i, int j, int k, int state) -> double&;   // 显式 state
    auto operator[](int i, int j, int k, int state) const -> double;

    auto print_info() const -> void;
};
```

### 索引公式

0-based 坐标 `(i, j, k)` 对应的扁平索引（文件存储顺序与访问顺序一致，无需重排）：

```
idx = state × n1 × n2 × n3 + i × n2 × n3 + j × n3 + k
```

其中 `k` 变化最快，`i` 变化最慢，与 Fortran 参考代码中的映射一致。

### 设计要点

1. **索引设计——与 EIGEN 的差异**：
   - `RHO::operator[](int i, int j, int k)` **不做任何运行时检查**，直接计算索引并返回引用。这是有意为之——`[i, j, k]` 是最常用的访问模式，每次判断 `nstate` 会拖慢高频访问的热路径。
   - 与 `EIGEN` 的 `Array3d::operator[](iband, ikpt)` 不同：后者在 `islda != 1` 时抛异常，因为语义上是"省略最后一个自旋维度"；而 `RHO` 的 `[i, j, k]` 语义明确为"访问第一个 state"，不需要防御性检查。
   - 显式 state 的 `operator[](i, j, k, state)` 不做 nstate 越界检查，由调用方保证。
   
2. **nstate 说明**：
   - `nstate` 是早期文件格式字段。实际中，当存在第二个自旋时，数据写在独立的 `OUT.RHO_2` 文件中（而非在同一文件增加 state），因此 `nstate` 一般为 1。
   
3. **文件格式兼容性**：header record 通过前导长度标记区分新/旧格式：
   - 20 字节 → 5 int（含 `nstate`）
   - 16 字节 → 4 int（旧格式，`nstate = 1`）

4. **一次性读取**：文件约 100KB（实测），构造函数中读取所有元数据和数据。文件读取完毕立即关闭。Rule of Zero。

5. **VR = RHO 别名**：`VR` 与 `RHO` 是同一类型，两者在代码中可互换使用。`src/io/VR.cppm` 仅包含：

```cpp
export using VR = RHO;
```

## `ATOM` 类

`ATOM` 类封装了对 `atom.config` 文本格式文件的读取操作。与其它 IO 类不同，`atom.config` 是**纯文本**（非 Fortran binary）格式，因此采用 `fopen("r")` + `fgets` 的行式解析方式，而非 `fread` + Fortran record 格式。

关于文件格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md#atomconfig)。

### `ATOM` 类

```cpp
export class ATOM {
public:
    explicit ATOM(const std::string& filename);
    ~ATOM();

    ATOM(const ATOM&) = delete;
    auto operator=(const ATOM&) -> ATOM& = delete;
    ATOM(ATOM&&) noexcept;
    auto operator=(ATOM&&) noexcept -> ATOM&;

    auto print_info() const -> void;

    // parsed data (as-read order, unsorted)
    int natom{};
    Lattice lattice;
    std::vector<int> species;
    std::vector<double> x, y, z;

    // species analysis (computed during construction)
    int ntyp{};
    std::vector<int> zvals;       // ntyp
    std::vector<int> type_counts; // ntyp
    std::vector<int> atom_types;  // ntyp
    std::vector<int> sorted_idx; // natom
};
```

### 设计要点

1. **文本解析**：与其它类的最大区别。使用 `fgets` 逐行读取，`strtol`/`strtod` 链式解析数值。第一行读 `natom`，然后进入卡片循环匹配关键词。

2. **卡片循环（card-loop）**：文件以关键词为驱动，卡片（LATTICE、POSITION 等）可自由排列。解析器逐行去除前后空格，匹配已知关键词：
   - `LATTICE` → 3 行晶格矢量，构造 `Lattice(A, LengthUnit::Angstrom)`
   - `POSITION` → `natom` 行原子数据（忽略 `imove_*` 标志）
   - 其它行 → 跳过

3. **解析辅助函数**：共同的解析逻辑放在模块内匿名命名空间（anonymous namespace）的 `parseAtomConfigFile` 函数中。匿名命名空间保证符号不跨模块泄露，避免 ODR 冲突。主实现和 `archived` 实现均调用此函数。

4. **species 分析**：构造中解析完成后调用 `analyzeSpecies()`，计算：
   - `ntyp` / `zvals[itype]` / `type_counts[itype]`：去重排序后的 species 种类、原子序数、每类数量
   - `atom_types[iatom]`：第 `ia` 个原子属于哪一类型（按 `zvals` 的排序）
   - `sorted_idx[new] = old`：将原子按类型分组的排列（`stable_sort`，同类型内保持原序）
   
   数据存储选型上，采用**分离的 vector**（`species`、`x`、`y`、`z`）而非 `vector<Atom>` 结构体数组。利弊：
   - 分离 vector 允许调用方直接传递单个数组（如 `species.data()`），无需先拆包
   - 只关心 species 时无需附带位置数据
   - 排序需要额外维护一个排列数组（`sorted_idx`），不能直接对容器 `sort`
   - 若改为结构体，排序更直接但失去扁平数组的灵活性。当前项目倾向于扁平数组风格。

5. **文件句柄生命周期**：与 OCC 一致——使用局部 `FILE*`，解析完成后立即关闭。try/catch 保护解析过程。Rule of Zero。

6. **单位**：晶格矢量从 Å 读入，立即通过 `Lattice` 转换为 Bohr。分数坐标不涉及单位转换。

### 元素符号 ↔ 原子序数转换

```cpp
static auto elementName(int atomic_number) -> std::string_view;
static auto atomicNumber(std::string_view name) -> int;
```

用法：

```cpp
ATOM::elementName(6)   // → "C"
ATOM::atomicNumber("Au")  // → 79
```

实现要点：

1. **数据源**：内部维护 `constexpr element_symbols[113]` 数组（下标 0 空置，1–112 对应 IUPAC 元素符号），数据与参考 Fortran 代码 `gen_element_name_number.f90` 的 112 个元素表一致（仅 112 号由 "Ch" 改为现代符号 "Cn"）。
2. **`elementName(z)`**：下标直接查表，Z 超出 1–112 范围抛 `std::out_of_range`。
3. **`atomicNumber(name)`**：线性扫描 1–112 匹配，未匹配抛 `std::invalid_argument`。ntyp 通常很小（实际算例 ≤ 5 种），线性扫描的性能开销可忽略。

### 迭代视图

ATOM 提供三种 range-for 可遍历的视图，用于访问物种类型和原子数据，统一用 `auto&&` 接收元素（详见下方「为什么用 `auto&&`」）：

```cpp
// TypeView：遍历物种类型
for (auto&& t : atom.eachType())          // ATOM::TypeEntry
    std::println("Z = {}, count = {}", t.z, t.count);

// AtomView：遍历指定类型的所有原子（按 Z 分组后的顺序）
for (auto&& a : atom.eachAtom(ityp))      // ATOM::AtomEntry
    std::println("species = {}, frac = ({}, {}, {})", a.species, a.x, a.y, a.z);

// SpecieView：按原始读取顺序遍历所有原子
for (auto&& a : atom.eachSpecie())        // ATOM::AtomEntry
    std::println("{}", a.species);
```

#### 为什么用 `auto&&`

`operator*()` 按值返回 `TypeEntry` 或 `AtomEntry`，产生临时对象（右值）：

```cpp
auto operator*() const -> TypeEntry { return {a_->zvals[it_], a_->type_counts[it_]}; }
```

| 写法               | 能否编译 | 说明 |
|-------------------|---------|------|
| `auto& a`         | ❌      | 非 const 左值引用不能绑右值 |
| `const auto& a`   | ✅      | const 引用可延寿，但只读 |
| `auto&& a`        | ✅      | 转发引用，万能绑定 |

range-for 展开后 `auto&& a = *__it`，`auto&&`（转发引用）可以绑定到任何值类别，因此是通用写法。

#### 实现方式

三个视图类均定义为 `ATOM` 的**嵌套类**（随 `export class ATOM` 自动导出），每个视图内部嵌套自己的 `Iterator` 类：

```
ATOM
├── TypeView          — eachType() 的返回类型
│   └── Iterator      — 内部迭代器，*it → TypeEntry
├── AtomView          — eachAtom(ityp) 的返回类型
│   └── Iterator      — 内部迭代器，*it → AtomEntry
└── SpecieView        — eachSpecie() 的返回类型
    └── Iterator      — 内部迭代器，*it → AtomEntry
```

每个迭代器实现 range-for 所需的三个操作：

```cpp
auto operator*()  const -> EntryType;   // 解引用
auto operator++() -> Iterator&;         // 前进
bool operator!=(const Iterator& o);     // 判等
```

View 类提供 `begin()` / `end()` 返回迭代器对。

#### `friend` 的作用

```cpp
class TypeView {
    friend class ATOM;                 // 只允许 ATOM 调用私有构造
    TypeView(const ATOM* a);
};
```

**不是**让 `TypeView` 能访问 `ATOM` 的私有成员——`ATOM` 的 `zvals`、`ntyp` 等本来就是 `public` 的。而是**让 `ATOM` 能构造 `TypeView`**：把构造设为 `private` 防止用户绕过 `eachType()` 直接实例化视图。

嵌套类读外围类的 `public` 成员不需要 `friend`；外围类读嵌套类的 `private` 成员才需要。详见 [`cpp_conventions.md`](cpp_conventions.md) 的迭代器设计风格约定。

#### 与 eachSpecie 的关系

```
每个 atom.eachType() + eachAtom 组合  →  所有原子恰好一次，按 Z 分组
atom.eachSpecie()                     →  所有原子恰好一次，原始顺序
```

覆盖的原子集合相同，**顺序不同**。以下例 `species = [13, 13, 7, 7]`（原始顺序），`sorted_idx = [2, 3, 0, 1]`：

```cpp
eachSpecie:         13, 13, 7, 7     ← 原始顺序
Type+Atom:           7,  7, 13, 13   ← Z 分组排序
```

如果原始文件已按 Z 排序，两种遍历碰到的顺序恰巧相同——但这是巧合，不是保证。
