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

# OUT.VR

最终目标是当我使用 [i,j,k] 的时候，我访问的是 xyz 的坐标

讨论1. 用 `vector` 连续存储，还是 `vector<vector<vector>>` 存储 

# 文档

文件格式分析写入 note/file_formats.md，EIGEN 类设计写入 design/io_module
