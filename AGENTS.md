# effmass 项目文档

## 项目概述

`effmass` 是一个第一性原理电子结构计算工具，用于读取平面波基组下的波函数文件（`OUT.GKK`、`OUT.WG`、`OUT.EIGEN` 等），计算得到有效质量。目前项目处于早期开发阶段，已完成 IO 模块和赝势读取模块的基础实现。

## 技术栈

- **主要语言**: C++23
- **编译器**: `clang++`（必须，依赖 libc++ 的实验性 `import std` 支持）
- **构建工具**: CMake ≥ 3.31、Ninja
- **文档**: Markdown

### C++23 模块说明

本项目使用 C++23 模块（`.cppm`）和实验性 `import std`，通过 CMake 的 `CMAKE_EXPERIMENTAL_CXX_IMPORT_STD` 机制支持。所有标准库内容通过 `import std;` 引入，不使用传统 `#include <...>` 方式引入标准库头文件（模块文件内部允许 `#include` 第三方库头文件，如 `pugixml.hpp`）。

## 构建与测试命令

```bash
# 配置构建（必须使用 Ninja generator）
cmake -B build -G Ninja

# 编译
cmake --build build

# 运行测试（CTest）
ctest --test-dir build
```

> **注意**: C++23 modules + `import std` 为实验性功能，依赖 libc++。若构建失败，删除 `build/` 目录后重新配置通常可解决问题。

## 项目结构

```
.
├── AGENTS.md         # 本文件 - AI 编码代理指南
├── CMakeLists.txt    # 根 CMake 配置
├── docs/             # 开发文档
├── src/              # 源代码
├── vendor/           # 第三方依赖
└── test/             # 测试代码与测试数据
```

### `src/` 目录

| 文件/目录 | 说明 |
|-----------|------|
| `io.cppm` | IO 模块：导出 `GKK`（`OUT.GKK` 读取器）、`WG`（`OUT.WG` 读取器） |
| `pseudo.cppm` | 赝势总模块：重新导出 `pseudo.ncpp_upf`、`pseudo.uspp_upf`、`pseudo.paw_upf` |
| `pseudo/ncpp-upf.cppm` | `NCPPUPF` 类：完整实现 UPF v.2 模守恒赝势解析 |
| `pseudo/uspp-upf.cppm` | `USPP` 类：骨架，构造函数抛出未实现异常 |
| `pseudo/paw-upf.cppm` | `PAW` 类：骨架，构造函数抛出未实现异常 |

### `docs/` 目录

| 文件/目录 | 说明 |
|-----------|------|
| `index.md` | `docs/` 目录索引 |
| `roadmap.md` | 项目整体开发路线图 |
| `progress.md` | 当前进度与修改说明 |
| `ideas.md` | 潜在改进方向（AI 编辑） |
| `proposal/` | 技术提案与改进方案 |
| `design/` | 代码设计决策与约定 |
| `note/` | 先验知识笔记（公式、文件格式等） |


## C++ 编码约定

### 命名规范

| 类型 | 命名规范 | 示例 |
|------|----------|------|
| 常量 | 全大写，下划线分隔 | `BOHR_RADIUS_ANGSTROM` |
| 类 | 大驼峰（PascalCase） | `GKK`, `NCPPUPF` |
| 结构体 | 大驼峰（PascalCase） | `GKKMetadata`, `KPointGVecs` |
| 成员变量 | 后缀下划线 | `meta_`, `fp_`, `current_kpt_` |
| 函数 | 小驼峰（camelCase） | `loadKPoint()`, `readMetadata()` |
| 局部变量 | 小驼峰 | `record_len`, `ikpt` |

### 函数返回类型写法

除构造函数和析构函数外，**所有函数**统一采用尾随返回类型写法：

```cpp
auto function_name(args...) -> ReturnType;
```

运算符重载函数也适用此写法。

### 物理常量与单位

- 禁止使用 Magic Number。所有物理量和单位转换系数必须使用 `constexpr` 定义，全大写下划线分隔命名。
- 内部计算统一使用 Hartree 原子单位制。

### 晶格矢量索引约定（约定-001）

从 Fortran 文件读取的晶格矢量 `AL[3][3]` 采用**向量紧密安排**：
- `AL[n][c]` 中 `n` 为晶格矢量编号（0=a1, 1=a2, 2=a3），`c` 为分量（0=x, 1=y, 2=z）。
- `AL[n]` 在内存中连续，可直接作为 `double[3]` 使用。

## 模块级设计约定

各模块的详细设计决策（错误处理策略、内存管理、性能优化等）已移至对应的设计文档：

- **IO 模块**（`GKK`、`WG`）：详见 `docs/design/io_module.md`
  - 错误处理、RAII 文件句柄、延迟加载、偏移预计算、缓存策略、缓冲区复用等
- **赝势模块**（`NCPPUPF` 等）：详见 `docs/design/pseudo_module.md`
  - UPF 解析错误前缀、Rule of Zero 等

## 测试

- 测试框架：不使用第三方测试框架，采用 `main()` 函数 + `try/catch` + `std::runtime_error` 的方式编写断言式测试。
- 测试运行：通过 CTest 调用（`ctest --test-dir build`）。
- 测试工作目录：`${CMAKE_SOURCE_DIR}`（即项目根目录），因此测试代码中使用相对路径 `"test/test_ncpp_upf/..."` 访问数据文件。

## 文档维护约定

- `docs/index.md` 是 `docs/` 目录的索引，每添加一个新文档，应在 `index.md` 中插入简洁概括。
- `docs/roadmap.md` 使用任务列表标记进度：`- [x]` 已完成，`- [/]` 进行中，`- [ ]` 未开始。
- `docs/progress.md` 记录当前正在做的工作和已完成的修改说明。
- `docs/ideas.md` 记录潜在的改进方向，由 AI 和开发者共同编辑。
- 每次完成一个进度后，评估是否需要更新 `docs/` 目录。

## Git commit messages

采用 Conventional Commits 格式：

```
<type>(<scope>): <subject>
```

允许的类型：`feat`, `fix`, `test`, `refactor`, `chore`, `style`, `docs`, `perf`, `build`, `ci`, `revert`。
