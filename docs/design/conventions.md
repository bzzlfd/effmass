# 项目约定（Conventions）

## [Convention-001] 晶格矢量索引语义 - 向量紧密安排

### 问题

从 Fortran 读取存储晶格矢量的 `AL(3,3)` 到 C++ 时有两种约定可选，新来者无法从代码中直接判断采用了哪种。

> 注：在 Fortran 代码中，这三个晶格矢量以 `AL(3,3)` 的变量名、按 column-major 的方式存储。

### 本项目的约定（约定 B：向量紧密安排）

- `lattice.A()[n][c]` 中，**n 对应晶格矢量编号 (a1/a2/a3)，c 对应分量 (x/y/z)**
- 即 `lattice.A()[n][0..2]` = 第 n 个晶格矢量的三个分量（内存连续）
- 这种约定符合 C 语言习惯，`A()[n]` 可直接作为 `double[3]` 使用

### 关键代码

```cpp
// 约定: A[n][c] 其中 n=vector, c=component
std::span<const double, 9> flat;
for (int n = 0; n < 3; ++n) {
    for (int c = 0; c < 3; ++c) {
        A_[n][c] = flat[n * 3 + c];
    }
}
```

### 使用示例

```cpp
// 获取第 n 个晶格矢量 (a1/a2/a3) 的 x, y, z 分量
double an[3] = {lattice.A()[n][0], lattice.A()[n][1], lattice.A()[n][2]};

// 获取所有晶格矢量的 x 分量
double ax[3] = {lattice.A()[0][0], lattice.A()[1][0], lattice.A()[2][0]};
```

---

## [Convention-002] 倒空间物理量命名

| 符号 | 一般含义 | 代码命名 |
|------|----------|----------|
| `k` | k 点分数坐标 | `k_frac` |
| `G` | 倒格子矢量 | `g_vec_cart` |
| `K` | 语义解释为 $-(G+k)$，由 `io.GKK` 管理 | `Kx`, `Ky`, `Kz` |
| 动能 | 平面波动能项 $\|G+k\|^2 / 2$ | `kinetic` |

### 整数索引前缀

| 含义 | 前缀 |
|------|------|
| k 点索引 | `ikpt` |
| 能带索引 | `iband` |
| 节点索引 | `inode` |
| G 向量索引 | `ig` |
