# io.EIGEN 已完成

`io.EIGEN` 模块已实现并测试通过（`test/test_eigen.cpp`），代码位于 `src/io/EIGEN.cppm`。

实现内容：
- `EIGENMetadata` 结构体（含 islda, nkpt, nband, nref_tot_8, natom, nnode, is_SO）
- `KPointVec` 结构体（k 点坐标，`.x/.y/.z` 命名成员 + `[0]/[1]/[2]` 下标访问）
- `Array3d` 类（三维本征值数组，一维连续存储，支持多参数 `[iband, ikpt, ispin]`、`[iband, ikpt]` 及链式 `[ispin][ikpt][iband]` 访问）
- `EIGEN` 类（一次性读取，文件格式兼容新/旧 header，验证 `islda_tmp`/`ikpt_tmp`）
- `print_info()` 成员函数

文档：
- `docs/note/file_formats.md` 补充 OUT.EIGEN 格式说明（包括 nref_tot_8 含义不明确的说明）
- `docs/design/io_module.md` 添加 EIGEN 类设计

# io.RHO / io.VR 已完成

`io.RHO`（及别名 `io.VR`）模块已实现并测试通过（`test/test_rho.cpp`），代码位于 `src/io/RHO.cppm`、`src/io/VR.cppm`。

实现内容：
- `RHOMetadata` 结构体（含 n1, n2, n3, nnodes, nstate）
- `RHO` 类：实空间网格数据读取，flat vector 存储 `[state][i][j][k]`，通过 `operator[]` 提供 `[i, j, k]` 和 `[i, j, k, state]` 两种索引方式
- `VR` 是 `RHO` 的类型别名（`using VR = RHO;`）
- 文件格式兼容新/旧 header（4 int 旧格式 vs 5 int 新格式）
- `print_info()` 成员函数

设计要点：
- `operator[](i, j, k)` 不做运行时 nstate 检查（与 EIGEN 的 Array3d 不同），保持热路径性能
- nstate 一般情况为 1；多自旋时由独立文件 `OUT.RHO_2` 处理

文档：
- `docs/note/file_formats.md` 补充 OUT.VR/OUT.RHO 格式说明
- `docs/design/io_module.md` 添加 RHO/VR 类设计

# 文档

文件格式分析写入 note/file_formats.md，EIGEN/VR/RHO 类设计写入 design/io_module
