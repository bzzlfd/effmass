# docs/ 目录索引

## 开发规划

- `roadmap.md`：项目整体开发路线图（IO → 线性代数 → 非局域势）
- `progress.md`：当前进度与修改说明
- `ideas.md`：潜在的改进方向与技术债务记录

## 知识笔记（note/）

- `note/atomic_units.md`：Hartree 原子单位制的约定与单位转换
- `note/fortran_record_format.md`：Fortran unformatted binary 的 record 格式
- `note/file_formats.md`：OUT.GKK、OUT.WG、OUT.EIGEN、OUT.KPT 的文件格式解析

## 设计决策（design/）

- `design/lattice_vectors_indexConvention.md`：**[约定-001]** 晶格矢量采用向量紧密安排（`AL[n][c]`）的索引约定
- `design/physical_constants.md`：物理常量命名规范与 `constexpr` 使用约定
- `design/io_cppm_design.md`：`io.cppm` 中 `GKK` 与 `WG` 类的接口设计与实现细节
- `design/cpp_conventions.md`：C++ 编码约定（命名规范、函数返回类型写法等）
