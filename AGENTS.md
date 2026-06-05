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

## 项目结构

```
.
├── AGENTS.md         # 本文件 - AI 编码代理指南
├── CMakeLists.txt    # 根 CMake 配置
├── docs/             # 开发文档
├── src/              # 源代码
├── test/             # 测试代码与测试数据
├── tools/            # 辅助工具与脚本（尽量脱离本项目维护）
└── vendor/           # 第三方依赖
```

### `src/` 目录

| 文件/目录 | 说明 |
|-----------|------|
| `io.cppm` + `io/` | 读取和表示 GKK、WG、EIGEN、OCC、RHO、ATOM、VR 等 |
| `pseudo.cppm` + `pseudo/` | 赝势模块：文件解析层：`pseudo/io/*.cppm`；赝势数据接口层：`pseudo/*.cppm` |
| `H_psi.cppm` + `H_psi/` | 与 `H_psi` 计算有关的代码。Hamiltonian 模块：Hamiltonian（对整个计算目录/文件的抽象），内含内嵌对象和 Callable Functor 作用算符； H_psi 作用、Gradient、Hessian |
| `math.cppm` + `math/` | 数学工具模块：FFT、Fourier-Bessel、球谐/球贝塞尔函数、线性代数 |
| `support/` | 支持模块（密度计算等） |
| `utils/` | 实用工具：vector3d、array2d、logger |
| `physical_constants.hpp` | 物理常量定义（头文件） |

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

详见：`docs/design/cpp_conventions.md` | `docs/design/conventions.md`

关键规则速览：
- **命名**：常量 `UPPER_SNAKE`，类/结构体 `PascalCase`，成员变量 `suffix_`，函数/局部变量 `camelCase`
- **尾随返回类型**：所有函数（除构造/析构外）统一使用 `auto fn() -> RetType`
- **物理单位**：禁止 Magic Number，内部计算统一使用 Hartree 原子单位制
- **晶格矢量索引**：`AL[n][c]`，`n` 为晶格矢量编号，`c` 为分量
- **向量命名的通常约定**：`k` 为分数坐标，`G` 为倒格子矢量，`K`为`k+G`，`kinetic = |G+k|²/2`

## 阅读路径

接手任务时，先阅读 `docs/progress.md` 了解当前进度和最近修改。其余设计文档和编码约定在实现过程中按需查阅。

## 问题处理

- **编译失败**：先尝试 `rm -rf build && cmake -B build -G Ninja && cmake --build build`
- **测试失败**：分析失败原因。如果是修改引入的 bug 则修复；如果是预期行为变化，更新测试断言
- **代码与文档矛盾**：以代码为准，但同步更新相关文档
- **需求不明确**：向用户提问澄清，不要自行假设
- **禁止改动**：不要修改 `vendor/` 下的第三方代码；不要引入新的测试框架依赖；大范围重构前先与用户确认

## 测试

- 测试框架：不使用第三方测试框架，采用 `main()` 函数 + `try/catch` + `std::runtime_error` 的方式编写断言式测试。
- 测试运行：通过 CTest 调用（`ctest --test-dir build`）。
- 测试代码组织：位于 `test/` 下，目录结构模仿 `src/`。如果测试涉及多个 `src` 模块，放在这些模块向上最先遇到的公共父节点上。
- 测试工作目录：`${CMAKE_SOURCE_DIR}`（即项目根目录），因此测试代码中使用相对路径 `"test/test_io_ncpp/..."` 访问数据文件。

## 文档维护约定

- `docs/index.md` 是 `docs/` 目录的索引，每添加一个新文档，应在 `index.md` 中插入简洁概括。
- `docs/roadmap.md` 使用任务列表标记进度：`- [x]` 已完成，`- [/]` 进行中，`- [ ]` 未开始。
- `docs/progress.md` 记录当前正在做的工作和已完成的修改说明。
- `docs/ideas.md` 记录潜在的改进方向，由 AI 和开发者共同编辑。
- 每次完成一个进度后，评估是否需要更新 `docs/` 目录。

## Git commit messages

> Do NOT run `git add`, `git commit`, `git push`, `git reset`, `git rebase` or any other git mutations without explicit permission. Always ask for confirmation before performing git mutations, even if the user has confirmed in earlier conversations.

### Commit Message Format (header)

Use the Conventional Commits format:
```
<type>(<scope>): <subject>
```

Allowed types: `feat`, `fix`, `test`, `refactor`, `chore`, `style`, `docs`, `perf`, `build`, `ci`, `revert`.

### body
If multiple changes are involved, list them with items.

### Co-author Declaration (footer)

Append a `Co-Authored-By` trailer at the end of the commit message.

Determine your identity from the system prompt, distinguishing two concepts:
- **Agent name** (`<agent>`): Which agent you are. The system prompt typically declares your identity at the beginning, e.g. "You are Claude Code", "You are GitHub Copilot", etc. Fill in the agent name — **do not** fill in the model name.
- **Model name** (`<model>`): Which LLM drives you. Usually accompanied by words like "model" or "powered by", e.g. "claude-opus-4-7", "deepseek-v4-flash", "gpt-5", etc.
- **Email** `[<email>]`: Optional.

Format:
```
Co-Authored-By: <agent> (<model>) [<email>]
```

If identity cannot be determined, use `AI Agent (unknown)`.
