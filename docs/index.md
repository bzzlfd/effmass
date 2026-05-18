# docs/ 目录索引

## 开发规划

- `roadmap.md`：项目整体开发路线图（IO → 线性代数 → 非局域势）
- `progress.md`：当前进度与修改说明
- `ideas.md`：和 AI 工作中临时笔记
- `proposal/`：潜在的改进方向与技术债务记录

## 知识笔记（note/）

- `note/atomic_units.md`：Hartree 原子单位制的约定与单位转换
- `note/fortran_record_format.md`：Fortran unformatted binary 的 record 格式
- `note/file_formats.md`：OUT.GKK、OUT.WG、OUT.EIGEN、OUT.KPT 的文件格式解析
- `note/upf_format.md`：UPF v.2 赝势格式的字段含义、变量说明、束缚态波函数归一化条件
- `note/clangd_cpp23_modules.md`：clangd 如何通过 `compile_commands.json` 自动解析 `import std`，无需硬编码 `std.pcm` 路径

## 设计决策（design/）

- `design/conventions.md`：**[Convention-001/002]** 晶格矢量索引约定与倒空间物理量命名规范
- `design/physical_constants.md`：物理常量命名规范与 `constexpr` 使用约定
- `design/io_module.md`：`io` 总模块及子模块 `io.GKK`、`io.WG`、`io.EIGEN`、`io.lattice` 的接口设计与实现细节
- `design/pseudo_module.md`：赝势读取模块的接口设计、`NCPPUPF` 类结构与 UPF 文件 tag 的映射
- `design/upf_vs_binary_io.md`：UPF（pugixml）与二进制 IO（GKK/WG/RHO/EIGEN）文件生命周期管理的对比分析
- `design/cpp_conventions.md`：C++ 编码约定（命名规范、函数返回类型写法等）

## 实现的模块

| 模块 | 路径 | 功能 |
|------|------|------|
| `io` | `src/io.cppm` + `src/io/` | 读取 OUT.GKK、OUT.WG 等文件，提供 K 矢量、晶格、波函数系数 |
| `pseudo` | `src/pseudo.cppm` + `src/pseudo/` | UPF 赝势文件解析（NCPPUPF） |
| `math` | `src/math.cppm` + `src/math/` | 球贝塞尔函数、实球谐函数、1D/3D FFT |
