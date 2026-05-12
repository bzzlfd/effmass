# IO 模块设计细节与接口

本文档介绍 `io` 模块及其子模块 `io.GKK`、`io.WG` 的设计思想与对外接口。关于文件二进制格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md)。

## 模块结构

```
io (aggregate)
├── io.GKK  →  src/io/GKK.cppm
└── io.WG   →  src/io/WG.cppm
```

`io.cppm` 是聚合模块（aggregate module），通过 `export import io.GKK; export import io.WG;` 重新导出两个子模块的所有公开接口。外部代码只需 `import io;` 即可使用全部 IO 功能。

## 数据结构

### `GKKMetadata`

存储 `OUT.GKK` 文件的元数据：

```cpp
export struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;
    double Ecut;              // cutoff energy (Hartree)
    double AL[3][3];          // lattice vectors (Bohr)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};
```

### `KVecs`

某个 k 点的 G 向量数据视图，通过指针指向类内部预分配的连续缓冲区：

```cpp
export struct KVecs {
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, G_x-k_x, G_y-k_y, G_z-k_z
};
```

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
const KVecs& loadKPoint(int ikpt);                      // 加载指定 k 点数据（带缓存）
int current_ikpt() const;                                     // 当前缓存的 k 点索引
const KVecs& currentData() const;                       // 当前缓存的数据视图
std::array<double, 3> inferCurrent_k() const;                 // 从 K=G-k 数据推断 k 点分数坐标 (fractional coordinate)
```

- `loadKPoint` 会自动管理一个单 k 点缓存：如果请求的是当前已缓存的 k 点，直接返回已有视图，避免重复文件 IO。
- 读取时会将各 node 的数据合并到内部预分配的连续缓冲区中。

## 设计要点

### 延迟加载

构造函数**不**加载全部 k 点数据，只读取元数据并预计算每个 k 点在文件中的偏移量。实际数据按需通过 `loadKPoint` 加载。

#### 偏移预计算

`computeOffsets()` 在构造时遍历所有 record，记录每个 k 点的起始文件偏移。这样 `loadKPoint` 可以直接 `fseek` 到目标位置，而不需要顺序跳过前面的数据。

#### 缓冲区复用

类在构造时根据 `mg_nx * nnodes` 预分配四个 `std::vector<double>` 缓冲区。每次 `loadKPoint` 复用这些缓冲区，避免频繁的内存分配。

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
    double AL[3][3];
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

`OUT.WG` 的头部元数据（`n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda, Ecut, AL, ngtotnod`）必须与 `OUT.GKK` 保持一致（仅多出 `mx`）。当前类内部**暂不做交叉校验**，由调用方负责验证。
