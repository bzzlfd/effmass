# 通用工具类型：`array2d` 与 `vector3d`

## `array2d<T, Rows, Cols>` — 2D 数组

### 模板签名

```cpp
export template<typename T, int Rows = -1, int Cols = -1>
struct array2d;
```

### 工作原理

`array2d` 根据 `Rows`/`Cols` 模板参数选择不同的内部存储实现：

- **动态版本** `array2d<T>`（`Rows = Cols = -1`，默认）：使用 `std::vector<T>` 存储，`rows`/`cols` 为运行时变量。**堆分配**，适用于编译期不确定大小的场景。
- **静态版本** `array2d<T, R, C>`（`Rows > 0 && Cols > 0`）：使用 `std::array<T, R*C>` 存储，`rows`/`cols` 为编译期常量。**栈分配**，适用于固定大小的 2D 数组。

两个版本共享相同的接口，切换只需改变模板参数：

```cpp
array2d<double>       mat;         // 动态，堆分配
array2d<double, 3, 3> mat;         // 静态，栈分配，3×3 固定大小
```

### 实现方式

通过**主模板 + 两个偏特化**实现：

```cpp
// 主模板（仅声明）
template<typename T, int Rows = -1, int Cols = -1>
struct array2d;

// 偏特化 1：动态版本（默认）
template<typename T>
struct array2d<T, -1, -1> { ... };

// 偏特化 2：静态版本（requires 约束确保 Rows/Cols > 0）
template<typename T, int Rows, int Cols>
    requires (Rows > 0 && Cols > 0)
struct array2d<T, Rows, Cols> { ... };
```

`requires` 子句确保编译器能正确区分两个偏特化：当 `Rows` 和 `Cols` 都是正整数时匹配静态版本，否则（默认 `-1`）匹配动态版本。

### 使用场景

| 场景 | 适用版本 | 原因 |
|------|----------|------|
| 编译器已知固定维度，如 3×3 矩阵 | 静态 `array2d<double, 3, 3>` | 栈分配，零开销 |
| 运行时确定维度，如 `n × n` 特征向量矩阵 | 动态 `array2d<double>` | 维度依赖输入数据 |

### 示例

```cpp
// 静态：倒格子矩阵，固定 3×3
array2d<double, 3, 3> B;
B[0, 1] = 1.0;            // 双参数下标
auto row0 = B[0];          // 行视图：std::span<double>

// 动态：投影算符 D 矩阵，大小由 UPF 文件决定
array2d<double> D(nb, nb);
D[i, j] = value;
```

---

## `vector3d<T>` — 3 维向量

### 模板签名

```cpp
export template<typename T>
struct vector3d;
```

### 成员与接口

```cpp
T x{}, y{}, z{};                      // 公开成员，支持结构化绑定

auto operator[](int i) -> T&;         // 索引访问（0→x, 1→y, 2→z）
auto data() -> T*;                   // 连续内存指针

auto norm() -> T;                     // 模长
auto norm_squared() -> T;             // 模长的平方

auto operator+=(vector3d&) -> vector3d&;
auto operator-=(vector3d&) -> vector3d&;
auto operator*=(T) -> vector3d&;
auto operator/=(T) -> vector3d&;
```

以及自由函数：`+`, `-`, `*`（标量乘法）, `/`（标量除法）, `-`（取反）, `dot()`, `cross()`, `==`, `!=`, `<<`。

### 使用场景

替换散落的 3 元素数组和结构体。目前已替换位置：
- `KVecs::kPoint`（原 `std::array<double, 3>`）
- `EIGEN::kVec`（原自定义 `struct kVec`）

### 示例

```cpp
vector3d<double> v{1.0, 2.0, 3.0};
v.x = 4.0;                       // 成员访问
v[0] = 5.0;                      // 索引访问
auto [x, y, z] = v;              // 结构化绑定

auto n  = v.norm();
auto d  = dot(v, w);
auto c  = cross(v, w);
auto v2 = v * 2.0;
```
