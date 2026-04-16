# `io.cppm` 设计细节与接口

本文档介绍 `src/io.cppm` 中 `GKK` 类的设计思想与对外接口。关于 `OUT.GKK` 等文件的二进制格式定义，请参阅 [`note/file_formats.md`](../note/file_formats.md)。

## 模块导出

```cpp
export module io;
```

`io.cppm` 是一个 C++23 模块接口文件，导出与文件 IO 相关的结构和类。

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

### `KPointGVecs`

某个 k 点的 G 向量数据视图，通过指针指向类内部预分配的连续缓冲区：

```cpp
export struct KPointGVecs {
    std::size_t ng;           // number of G-vectors
    const double *g, *gx, *gy, *gz;  // |G+k|²/2, Gx, Gy, Gz
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

支持移动构造和移动赋值。移动后会自动更新 `KPointGVecs` 内部的指针，使其指向新的缓冲区。

> 拷贝被显式删除，因为类持有 `FILE*` 句柄和大量缓冲区。

### 公共接口

```cpp
const GKKMetadata& metadata() const;           // 获取元数据
const KPointGVecs& loadKPoint(int ikpt);       // 加载指定 k 点数据（带缓存）
int currentKPoint() const;                     // 当前缓存的 k 点索引
const KPointGVecs& currentData() const;        // 当前缓存的数据视图
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
