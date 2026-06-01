# ncpp 接口重构

以 docs/reference/unified_pseudopotential_format.md 为 UPF 文件参考。
它提供了完整ncpp数据列表

## ncpp.cppm 和 ncpp.upf.cppm 角色
1. ncpp.upf.cppm 是纯粹的 upf 文件的数据接口抽象
2. ncpp.cppm 是纯粹的 ncpp 的数据接口抽象

  1. 把 *.upf.cppm 精简成 upf.cppm
  2. upf 本身可以表示三种pseudopotential，upf.cppm 也应该表示三种PsP。
  3. ncpp.cppm 可以从 upf 文件构造出 ncpp，也可以从 （比如psp8 文件）构造出ncpp，然后作为数据接口提供给别人


## ncpp.cppm 提供接口原则
1. 不使用 ncpp.upf.cppm 中变量命名方式，而是使用更贴近其含义的命名
  1. D 矩阵（Dion）在ncpp中改名为B矩阵，这是Vanderbilt论文中的命名方式
2. NCPP 类除了有对应的数据接口，还应提供若干便利的成员函数

**每个 accessor 返回标准 C++ 类型**：`string_view`, `span<const double>`, `int`, `double`, `bool`。不返回格式相关的复合类型。

未来 ncpp.cppm 中的对象不一定来自 UPF，
未来多格式时变为 `std::variant<UPF, PSP8, ...>`，所有 accessor 内部 dispatch。不设 factory——以后再说。

## ncpp.cppm 数据存储
成员变量
- Info
- Header
- Mesh
- NLCC (optional)
- Local
- Nonlocal
- Semilocal (optional, only for norm-conserving)
- PsWFc (optional)
- RhoAtom

其中，Semilocal 数据是可选的。我们后面会有程序判断，如果ncpp是semilocal的，那么报错，说我们现在只支持 Nonlocal ，为这个程序判断做好准备。
mesh 中的数据与 beta 中的数据按位置对应。但是beta中在某处往后全都是0, beta存储数据长度小于 mesh
其内部数据组织方式参考 ncpp.upf.cppm，但是内部变量名要选择更突出物理含义的



## ncpp.cppm 成员函数（不完全）
- public
  - NonlocalByL
- private
  - inferMeshType (在初始化时调用)


## 当前命名列表

### 类型

| 名称 | 用途 |
|------|------|
| `NCPP` | 模守恒赝势数据类 |
| `MeshType` | 网格类型枚举：`Uniform`, `Exponential`, `Unknown` |
| `Relativistic` | 相对论处理枚举：`None`, `Scalar`, `Full` |
| `PseudoType` | 赝势类型枚举：`NC`, `SL`, `US`, `PAW`, `Coulomb` |
| `NonlocalByL` | 按 l 分组的非局域数据 |

### NCPP 数据成员

| 成员 | 类型 | 说明 |
|------|------|------|
| `info` | `Info` | 元信息 |
| `info.element` | `std::string` | 元素符号 |
| `info.z_valence` | `double` | 价电子数 |
| `info.relativistic` | `Relativistic` | 相对论类型 |
| `info.functional` | `std::string` | 交换关联泛函 |
| `info.pseudo_type` | `PseudoType` | 赝势类型 |
| `header` | `Header` | 数值参数 |
| `header.l_max` | `int` | 最高角动量 |
| `header.l_local` | `int` | 局域势角动量 |
| `header.mesh_size` | `int` | 网格点数 |
| `header.number_of_wfc` | `int` | 波函数数量 |
| `header.number_of_proj` | `int` | 投影子数量 |
| `header.has_so` | `bool` | 是否有自旋轨道耦合 |
| `header.core_correction` | `bool` | 是否有 NLCC |
| `header.total_psenergy` | `double` | 赝势总能量 |
| `mesh` | `Mesh` | 径向网格 |
| `mesh.r` | `std::vector<double>` | 径向坐标 |
| `mesh.rab` | `std::vector<double>` | dr/di 导数 |
| `nonlocal` | `Nonlocal` | KB 非局域投影子 |
| `nonlocal.beta` | `std::vector<std::vector<double>>` | 投影子函数 βᵢ(r) |
| `nonlocal.angular_momentum` | `std::vector<int>` | 投影子角动量 |
| `nonlocal.cutoff_index` | `std::vector<int>` | 截断半径索引 |
| `nonlocal.cutoff_radius` | `std::vector<double>` | 截断半径值 |
| `nonlocal.b_matrix` | `DenseMatrix<double>` | B 矩阵 Bᵢⱼ |
| `wfc` | `PsWFc` | 赝波函数 |
| `wfc.chi` | `std::vector<std::vector<double>>` | 波函数 χᵢ(r) |
| `wfc.kchi` | `std::vector<int>` | 波函数有效长度 |
| `wfc.lchi` | `std::vector<int>` | 波函数角动量 |
| `wfc.occupation` | `std::vector<double>` | 占据数 |
| `wfc.label` | `std::vector<std::string>` | 标记 (e.g. "3S") |
| `vloc` | `std::vector<double>` | 局域势 V_loc(r) |
| `rho_atom` | `std::vector<double>` | 原子电荷密度 ρ(r) |

### NonlocalByL 成员

| 成员 | 类型 | 说明 |
|------|------|------|
| `beta` | `std::vector<std::vector<double>>` | 当前 l 的投影子 |
| `bMatrix` | `DenseMatrix<double>` | 当前 l 的 B 子矩阵 |  改称 B 

### NCPP 成员函数

| 函数 | 说明 |
|------|------|
| `meshType()` | 推断网格类型 |
| `nonlocalByL(l)` | 提取指定 l 的投影子和 B 子矩阵 |



info和header合并成 meta，
mesh改为radialMesh
b_matrix，bMatrix 改为 B
nonlocal 改为 V_nonlocal
local 改为年 V_local
wfc 改为 pseudoWaveFunction
