# 当前进度

## 已完成

### Lattice 内部扁平化存储 + 全部路径汇聚到 `setLatticeFromFlat`

- `src/io/lattice.cppm`：
  - 新增私有嵌套类型 `array2d`（包装 `std::array<double, 9>`，重写 C++23 双参 `operator[](int i, int j)` → `A_[n, c]` 索引语法），替代 `std::array<std::array<double, 3>, 3>` 作为 `A_` 的内部存储，内存扁平连续
  - `A()` 改为从 `A_` 按元素构建 `std::array<std::array<double, 3>, 3>` 返回
  - `setLattice(array)` 改为展平到 `std::array<double, 9>` 后委托 `setLatticeFromFlat`
  - `Lattice(ilist)` 改为调用 `setLattice(A, unit)`（不再直接调 `setLatticeFromFlat`），保持三位构造一致委托三位 set
  - 六条入口路径全部汇聚到 `setLatticeFromFlat`（唯一写入 `A_` 的位置）
  - `computeReciprocalLattice()`：`A_[i].data()` → `&A_[i, 0]`
- `docs/design/io_module.md`：更新转发路径图为一元汇聚结构，补充 `array2d` 内部存储说明

### IO 模块重构 — 拆分 `io.cppm` 为子模块

- `src/io/GKK.cppm`：新建 `io.GKK` 子模块，包含 `GKKMetadata`、`KVecs`、`GKK` 类及全部实现
- `src/io/WG.cppm`：新建 `io.WG` 子模块，包含 `WGMetadata`、`WGCoeffs`、`WG` 类及全部实现
- `src/io.cppm`：重写为聚合模块，通过 `export import io.GKK; export import io.WG;` 重新导出
- `src/physical_constants.hpp`：新增共享物理常量头文件（`BOHR_RADIUS_ANGSTROM`），供 GKK/WG 共同包含
- `CMakeLists.txt`：更新 `io` 目标的 `FILE_SET cxx_modules`，加入两个新模块文件
- 文档同步更新：`AGENTS.md`、`docs/design/io_module.md`

### NCPP UPF 解析器 — beta / chi 数据截断

- `src/pseudo/ncpp-upf.cppm`：`beta` 按 `kbeta` 截断，`chi` 按尾部零扫描截断，并新增 `kchi[]` 字段保存有效长度。
- `test/test_ncpp_upf.cpp`：更新尺寸校验，验证截断后长度与 `kbeta`/`kchi` 一致，波函数归一化校验通过。
- `docs/design/pseudo_module.md`：补充"数据截断策略"设计说明。

### NCPP UPF 解析器 — mesh 类型推断与按 l 筛选非局域项

- `src/pseudo/ncpp-upf.cppm`：
  - 新增 `MeshType` 枚举（`Uniform`, `Exponential`, `Unknown`）。
  - 新增 `NCPPUPFNonlocalByL` 结构体。
  - `NCPPUPF` 新增 `meshType()` 函数：通过分析 `r` 和 `rab` 的数值关系推断网格类型（均匀或指数/对数）。
  - `NCPPUPF` 新增 `nonlocalByL(int l)` 函数：返回指定角动量 `l` 的 beta 子列表和对应的 `D_ij` 子矩阵。
- `test/test_ncpp_upf.cpp`：新增 `meshType()` 和 `nonlocalByL()` 的测试用例，包括 beta 数量、dion 尺寸、数值一致性校验。
- `docs/design/pseudo_module.md`：更新公共接口与数据结构说明，补充 mesh 类型推断和 `nonlocalByL` 的设计文档。

## 已完成

### 赝势模块重构 — IO 层与算符层分离

- `src/pseudo/io/`：新建 IO 层目录，存放 UPF 文件读取器
  - `ncpp.upf.cppm`：原 `ncpp-upf.cppm` 迁移，模块名 `pseudo.ncpp_upf` → `pseudo.io.ncpp_upf`
  - `uspp.upf.cppm`：原 `uspp-upf.cppm` 迁移，类名 `USPP` → `USPPUPF`，模块名 `pseudo.uspp_upf` → `pseudo.io.uspp_upf`
  - `paw.upf.cppm`：原 `paw-upf.cppm` 迁移，类名 `PAW` → `PAWUPF`，模块名 `pseudo.paw_upf` → `pseudo.io.paw_upf`
- `src/pseudo/ncpp.cppm`：**新增** `NCPP` 算符骨架，模块 `pseudo.ncpp`，持有 `NCPPUPF` 数据副本
- `src/pseudo.cppm`：更新重新导出的子模块列表
- `CMakeLists.txt`：更新 `pseudo` 目标的模块文件路径
- 测试重构：
  - `test/test_io_ncpp.cpp`：原 `test_ncpp_upf.cpp` 迁移，数据路径更新
  - `test/test_ncpp.cpp`：**新增** `NCPP` 算符骨架测试
  - `test/CMakeLists.txt`：更新测试目标
- 文档同步更新：`AGENTS.md`、`docs/design/pseudo_module.md`

## 已完成

### GKK — 多视角数据表示与 `KVecsView`

- `src/io/GKK.cppm`：
  - 新增 `KVecsView` 位掩码枚举（`Cartesian`, `Spherical`, `Integer`），支持 `|`, `&`, `~` 组合。
  - 扩展 `KVecs` 结构体：
    - 新增球坐标 `r, theta, phi`
    - 新增整数索引 `iG, jG, kG`（原 `Gx, Gy, Gz` 已更名）
    - 新增每 k 点元数据：`kPoint`（k 分数坐标，由 `inferCurrent_k()` 推断）、`reciprocalLattice`（倒格子矩阵，构造时计算并缓存）
  - `GKK` 新增 `setDataView(view)` / `currentView()`：内部状态控制当前启用的数据表示。
  - `loadKPoint` 支持增量计算：同一 k 点切换 view 时不重复文件 IO，仅补充计算缺失的表示。
  - `currentData()` 返回的 `KVecs` 按 `currentView` 过滤：未启用的表示为空 `span`。
  - 新增 `computeSpherical()` 和 `computeIntegerIndices()` 私有方法：
    - 球坐标：`r = |K|`, `theta = acos(Kz/r)`, `phi = atan2(Ky, Kx)`
    - 整数索引：`(iG,jG,kG) = round(A^T * K_cart / (2π) + k_frac)`，复用已有的 `inferCurrent_k()`
  - 新增 `computeReciprocalLattice()`：从 `meta_.AL` 计算倒格子矩阵 `b_n = 2π * (a_{n+1} × a_{n+2}) / V`
  - 移动语义更新：修复新增缓冲区的 `span` 指向。
- 文档同步更新：`docs/design/io_module.md`

## 进行中

1. 特殊函数
2. ~~给 ncpp-upf.cppm 中的 NCPPUPF 添加更多函数~~
    1. ~~推断 mesh 中的网格是均匀的，还是指数的，还是未知的~~
    2. ~~获得角动量是 l 的 beta 子列表，dion子矩阵~~


1. 代码规范，函数之间两行。和python规范一样
2. ~~抽象lattice~~
  1. ~~写测试。用高斯消元法~~

## 已完成

### Lattice 抽象与 IO 模块整合

- `src/io/lattice.cppm`：新建 `io.lattice` 子模块，包含 `LengthUnit` 枚举和 `Lattice` 类。
  - `Lattice` 封装实空间晶格 `A`（Bohr）与倒空间晶格 `B`（Bohr⁻¹），构造时自动计算倒格子。
  - 支持 `LengthUnit::Bohr` 和 `LengthUnit::Angstrom` 单位传入/转换。
  - `A(LengthUnit)` / `B(LengthUnit)` 支持按指定单位返回数值矩阵（默认 Bohr）。
- `src/io/GKK.cppm`：`GKKMetadata` 中 `double AL[3][3]` 替换为 `Lattice lattice`；移除内嵌的 `computeReciprocalLattice()`，所有晶格操作统一委托给 `Lattice`。
- `src/io/WG.cppm`：`WGMetadata` 中 `double AL[3][3]` 替换为 `Lattice lattice`。
- `test/test_lattice.cpp`：新增测试，用增广矩阵高斯消元法手动求 `A` 的逆并转置，验证 `lattice.B()` 与 `2π * (A⁻¹)ᵀ` 一致（覆盖立方、FCC、Angstrom 转换）。
- `CMakeLists.txt` / `test/CMakeLists.txt`：加入 `io.lattice` 模块和 `test_lattice` 测试目标。
- 文档同步更新：`docs/design/io_module.md`、`docs/design/conventions.md`、`docs/design/multi_dimensional_array.md`、`docs/note/atomic_units.md`、`docs/index.md`。


