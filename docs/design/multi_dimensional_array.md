# 多维数组索引约定

## 背景

项目代码中同时存在两种二维数组的访问语法：

1. **`A[i][j]`** —— 嵌套 `std::vector`（或原生数组）的链式下标。
2. **`A[i,j]`** —— C++23 多参数 `operator[]` 配合 flat storage 的一次下标。

## 现状盘点

### `A[i][j]` 风格（嵌套容器 / 原生数组）

| 位置 | 类型 | 说明 |
|------|------|------|
| `ncpp-upf.cppm` | `std::vector<std::vector<double>> beta` | 非局域势的 beta 投射函数，每行被截断为不同长度 |
| `ncpp-upf.cppm` | `std::vector<std::vector<double>> chi` | 赝原子波函数，每行被截断为不同长度 |
| `io.cppm` | `std::vector<std::vector<int>> ngtotnod_` | 每 k 点、每 node 的 G-vector 计数，语义上是二维表但由嵌套 vector 承载 |
| `lattice.cppm` | `std::array<std::array<double, 3>, 3>` (in `Lattice`) | 晶格矢量，固定 3×3 尺寸，封装在 `Lattice` 类中 |
| `test_ncpp_upf.cpp` | `nl.beta[0][0]`、`wfc.chi[i][ir]` | 测试代码中对上述结构的访问 |

### `A[i,j]` 风格（flat storage + 多参数下标）

| 位置 | 类型 | 说明 |
|------|------|------|
| `ncpp-upf.cppm` | `Matrix dion` | NCPP 的非局域势 D_ij 矩阵，固定 `nb × nb` 方阵 |

其中 `Matrix` 的实现：

```cpp
struct Matrix {
    std::vector<double> data;
    int rows = 0;
    int cols = 0;

    auto operator[](int i, int j) const -> const double& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    auto operator[](int i, int j) -> double& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    // ...
};
```

## 讨论：是否需要统一？

**结论：不需要强制统一，但需按数据拓扑明确选择。**

两种风格对应两种本质上不同的二维数据结构：

| 维度 | `A[i][j]`（嵌套 vector） | `A[i,j]`（flat Matrix） |
|------|--------------------------|------------------------|
| **数据拓扑** | Jagged array（每行长度可变） | Dense matrix（行列固定） |
| **内存布局** | 行间不连续，每行独立堆分配 | 整块连续内存，缓存友好 |
| **典型场景** | `beta`、`chi` 等需要逐行截断 trailing zeros 的数据 | `dion` 等固定大小的矩阵 |
| **性能特征** | 访问多一次间接寻址，不适合大规模数值计算 | 适合密集线性代数运算 |

若强行统一，会带来反效果：

- 若把 jagged array 强行塞进 flat `Matrix`，需要额外维护每行长度/偏移，过度设计。
- 若把 dense matrix 改用嵌套 `std::vector<std::vector<T>>`，会丧失内存连续性，影响后续可能引入的密集线性代数操作性能。

## 约定

根据数据拓扑选择访问语法：

1. **Jagged array（每行长度不固定）**
   - 使用 `std::vector<std::vector<T>>`。
   - 访问语法为 `A[i][j]`。
   - 适用场景：需要逐行截断、每行有效长度不同的数据（如赝势的 `beta`、`chi`）。

2. **Dense matrix（行列尺寸固定）**
   - 使用 flat storage 的自定义 `Matrix`（或后续引入的类似结构）。
   - 访问语法为 `A[i,j]`（利用 C++23 多参数 `operator[]`）。
   - 适用场景：固定大小的数值矩阵（如 `dion`、平面波系数矩阵等）。

3. **极小固定尺寸（如 3×3）**
   - 允许直接使用原生数组 `T[N][M]` 或嵌套 `std::array`。
   - 访问语法为 `A[i][j]`。
   - 适用场景：晶格矢量 `lattice.A()[n][c]` 等。

## 示例

```cpp
// 1. Jagged array —— 每行长度不同
std::vector<std::vector<double>> beta;
beta[i][ir] = ...;          // OK: A[i][j]

// 2. Dense matrix —— 固定行列
Matrix dion;
dion[i, j] = ...;           // OK: A[i,j]

// 3. 极小固定尺寸 —— 封装在 Lattice 类中
Lattice lattice;
lattice.A()[n][c] = ...;    // OK: A[i][j]
```
