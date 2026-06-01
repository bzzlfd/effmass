# Pseudo-Potential Modules

本文档描述 `src/pseudo/` 下赝势模块的设计与接口，包含 **IO 层**（UPF 文件读取）与 **算符层**（计算中使用的赝势算符抽象）。

## 模块结构

```
src/pseudo.cppm          # 总集成模块: export module pseudo;
src/pseudo/
  ├── io/
  │   └── upf.cppm       # 统一 UPF v.2 读取 (pseudo.io.upf), 支持 NC/USPP/PAW
  └── ncpp.cppm          # NCPP 算符抽象 (pseudo.ncpp)
```

`pseudo.cppm` 仅负责重新导出子模块，本身不含实现。

## IO 层 — UPF 文件读取

### UPF — Unified Pseudopotential Format (v.2)

单一 `UPF` 类读取所有三种赝势类型（Norm-Conserving、Ultrasoft、PAW）的 UPF 文件。

### 参考资料

- [Quantum ESPRESSO: Unified Pseudopotential Format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
- Vanderbilt, D. (1990). *Soft self-consistent pseudopotentials in a generalized eigenvalue formalism*. Phys. Rev. B, 41, 7892. DOI: [10.1103/PhysRevB.41.7892](https://doi.org/10.1103/PhysRevB.41.7892)

### 公共接口

```cpp
export class UPF {
public:
    explicit UPF(const std::string& filename);

    auto header() const -> const UPFHeader&;
    auto mesh() const -> const UPFMesh&;
    auto localPotential() const -> std::span<const double>;
    auto nonlocal() const -> const UPFNonlocal&;
    auto wavefunctions() const -> const UPFWavefunction&;
    auto rhoAtom() const -> std::span<const double>;
};
```

### 数据结构与 UPF 文件 tag 的对应关系

```
UPF
├── UPFHeader             ← PP_HEADER (全部属性)
│   ├── pseudo_type       ← "NC" | "US" | "PAW"
│   ├── is_ultrasoft      ← USPP 标识
│   └── is_paw            ← PAW 标识
├── UPFMesh               ← PP_MESH
│   ├── r[]               ← PP_R
│   └── rab[]             ← PP_RAB
├── vloc_                 ← PP_LOCAL
├── UPFNonlocal           ← PP_NONLOCAL
│   ├── beta[][]          ← PP_BETA.1 ... PP_BETA.nbeta（截断至 kbeta）
│   ├── lll[]             ← PP_BETA.* / angular_momentum
│   ├── kbeta[]           ← PP_BETA.* / cutoff_radius_index
│   ├── rcut[]            ← PP_BETA.* / cutoff_radius
│   └── dion              ← PP_DIJ (DenseMatrix<double>)
├── UPFWavefunction       ← PP_PSWFC
│   ├── chi[][]           ← PP_CHI.1 ... PP_CHI.nwfc（截断尾部零）
│   ├── kchi[]            ← 有效长度（最后一个非零元素位置）
│   ├── lchi[]            ← PP_CHI.* / l
│   ├── oc[]              ← PP_CHI.* / occupation
│   └── labels[]          ← PP_CHI.* / label
└── rho_at_               ← PP_RHOATOM
```

**NC / US / PAW 支持现状**：
- **NC**（Norm-Conserving）：完全实现，读取全部字段
- **US**（Ultrasoft）：读取 header 后抛出 `"UPF: USPP reader not yet implemented"`，结构体已预留，读取逻辑待实现
- **PAW**（Projector-Augmented Wave）：读取 header 后抛出 `"UPF: PAW reader not yet implemented"`，结构体已预留，读取逻辑待实现

构造函数先读取 PP_HEADER，根据 `is_ultrasoft` / `is_paw` 标识判断类型。对于 NC 类型，继续读取剩余全部字段；对于 US/PAW 类型，因 PP_NONLOCAL 中额外结构（Q_ij 增广电荷等）尚未解析，暂不支持完整读取。

#### 扁平存储 vs. 嵌套容器

- **`beta`** 与 **`chi`**：使用 `std::vector<std::vector<double>>`。读取后按实际有效长度截断尾部零，因此每行长度 ≤ `mesh_size`，该尺度（mesh ~ 10³，nbeta ~ 10）下的间接寻址开销可忽略。
- **`dion`**：使用 `DenseMatrix<double>` 值类型（来自 `utils.matrix`）。$D_{ij}$ 是 $n_{\beta} \times n_{\beta}$ 对称矩阵，扁平连续存储更适合矩阵运算，且 C++23 多维 `operator[](i, j)` 提供了自然访问语法。
- **`vloc_`** 与 **`rho_at_`**：使用 `std::vector<double>` 并通过 `std::span<const double>` 暴露，方便调用方直接遍历而无需深拷贝。

### 数据截断策略

UPF 文件中的 `beta` 与 `chi` 数组在 `cutoff_radius_index`（`kbeta`）之后通常全部为零。若保留完整的 `mesh_size`（~10³）长度，会造成大量无效内存占用，并增加后续数值积分的循环开销。

**`beta` 的截断**：
- 依据 XML 属性 `cutoff_radius_index`（UPF 中为 1-based 索引）。
- 读取完整数组后执行 `resize(kbeta)`，vector 长度变为 `kbeta`（0 ~ kbeta-1）。
- `kbeta` 由文件作者设定，通常略大于实际最后一个非零值的位置，属于安全的保守截断。

**`chi` 的截断**：
- `chi` 没有现成的 cutoff 属性，因此从末尾向前扫描，去掉精确为 `0.0` 的尾部元素。
- 截断后的有效长度存入 `kchi[]`，与 `beta` 的 `kbeta[]` 形成对称设计。
- 扫描使用 `== 0.0` 比较；UPF v.2 中的尾部零是格式化输出产生的精确零值（如 `0.0000000000E+00`），不存在浮点噪声问题。

**截断后的索引一致性**：
- `mesh.r` 与 `mesh.rab` 仍保持完整 `mesh_size` 长度。
- 截断后的 `beta[i][ir]` 与 `chi[i][ir]` 的索引 `ir` 仍与 `mesh.r[ir]`、 `mesh.rab[ir]` 一一对应，只是循环上限从 `mesh_size` 变为 `beta[i].size()` 或 `chi[i].size()`。

### 解析流程

构造函数按以下顺序解析 XML 节点（UPF v.2 格式）：

```
<UPF>
  ├── <PP_HEADER/>          → 读取全部属性 → UPFHeader (含类型判断)
  ├── <PP_MESH>
  │     ├── <PP_R>          → 解析文本为 vector<double>
  │     └── <PP_RAB>        → 解析文本为 vector<double>
  ├── <PP_LOCAL>            → 解析文本 → vloc_[mesh]
  ├── <PP_NONLOCAL>
  │     ├── <PP_BETA.1>     → 属性 + 文本 → beta[0]
  │     ├── ...             → (循环 nbeta 次)
  │     └── <PP_DIJ>        → 解析文本 → dion[nbeta*nbeta]
  ├── <PP_PSWFC>
  │     ├── <PP_CHI.1>      → 属性 + 文本 → chi[0]
  │     └── ...             → (循环 nwfc 次)
  └── <PP_RHOATOM>          → 解析文本 → rho_at_[mesh]
```

**文本解析**：pugixml 加载 DOM 后，通过 `node.child_value()` 获取标签体内的文本，再用 `std::istringstream` 逐个读取 `double`。非数值词（如示例文件中的 `...`）被自动跳过。

**UPF v.2 特性**：
- 标签全大写（`PP_HEADER`, `PP_BETA.1`）
- 数组元素标签带索引后缀（`PP_BETA.1` ~ `PP_BETA.nbeta`）
- Header 数据全部存储在 XML 属性中（自闭合标签）

### 一次性全加载

UPF 文件通常只有几百 KB，pugixml 的 DOM 解析完全无压力。构造函数一次性读取全部数据到内存后即关闭文件句柄，类内部只持有 `std::vector`、`std::string` 等标准值类型成员。

因此 `UPF` 遵循 **Rule of Zero**——不显式声明析构函数、拷贝/移动构造函数或赋值运算符，编译器自动生成的语义完全正确。这也意味着该类可以自然地放入 `std::vector` 等标准容器中。

### 程序中不读取 `columns` 的原因

UPF v.2 的部分数据节点带有 `columns` 属性，例如：

```xml
<PP_R columns="4">
   1.0  2.0  3.0  4.0
   5.0  6.0  7.0  8.0
   ...
</PP_R>
```

在 `UPF` 的实现中，`parseTextToVector` 通过 `std::istringstream` 逐 token 读取，按空白字符（空格、制表符、换行）分割。因此无论每行有多少列，所有数值都会被正确提取。

**不读取 `columns` 的设计理由**：

1. **该属性是输出格式提示，不是语义约束**：`columns` 仅影响人类可读性和 Fortran `write` 语句的格式化输出，不改变数据本身的顺序和数量。

2. **pugixml 的文本节点已抹除格式信息**：`pugi::xml_node::child_value()` 返回的是节点内所有文本内容的拼接字符串，行边界信息已经丢失。即使想验证 `columns`，也无法从纯文本中可靠还原原始行结构。

3. **数值校验已足够**：每个数组读取后都会与 `mesh_size`（或 `nbeta`、`nwfc`）进行长度校验。如果格式化错误导致数据缺失或多余，`size mismatch` 异常会立即暴露问题。

### 错误处理

| 场景 | 异常类型 | 消息 |
|------|---------|------|
| XML 解析失败 | `std::runtime_error` | `UPF: failed to parse file: ...` |
| 缺少根节点 `<UPF>` | `std::runtime_error` | `UPF: missing root element <UPF>` |
| 缺少属性 | `std::runtime_error` | `UPF: missing attribute '...'` |
| 属性类型转换失败 | `std::runtime_error` | `UPF: invalid double/int/bool attribute ...` |
| 数组大小不匹配 | `std::runtime_error` | `UPF: <tag> size mismatch` |
| USPP 文件（未实现） | `std::runtime_error` | `UPF: USPP reader not yet implemented` |
| PAW 文件（未实现） | `std::runtime_error` | `UPF: PAW reader not yet implemented` |

## 算符层 — 赝势算符抽象

### NCPP — Norm-Conserving Pseudo-Potential Operator

`NCPP` 是计算过程中使用的模守恒赝势算符抽象，与 `UPF`（文件读取器）分离。可从 `UPF` 对象构造，也可从其他来源（如 psp8 文件）构造。

```cpp
export class NCPP {
public:
    explicit NCPP(const UPF& upf);

    // Access raw UPF data
    auto header() const -> const UPFHeader&;
    auto mesh() const -> const UPFMesh&;
    auto localPotential() const -> std::span<const double>;
    auto nonlocal() const -> const UPFNonlocal&;
    auto wavefunctions() const -> const UPFWavefunction&;
    auto rhoAtom() const -> std::span<const double>;
    auto upfData() const -> const UPF&;

    // Higher-level interfaces
    auto meshType() const -> MeshType;
    auto nonlocalByL(int l) const -> NonlocalByL;
};
```

未来将逐步添加物理接口，如局域势傅里叶变换、非局域 beta 投影在平面波基组下的矩阵元等。

## 测试

```
test/
├── CMakeLists.txt
test/data_io_upf/
    ├── Ge-spd-high.PD04.PBE.UPF          # NC 测试数据
    ├── Si.pbe-n-kjpaw_psl.1.0.0.UPF      # PAW 测试数据
    └── Si.pbe-nl-rrkjus_psl.1.0.0.UPF    # USPP 测试数据
test/pseudo/io/
    └── test_ncpp_upf.cpp                 # UPF 数据读取测试
test/pseudo/
    └── test_ncpp.cpp                     # NCPP 算符测试
```

`test_io_upf` 验证内容：
1. **Header 字段校验**：element, pseudo_type, z_valence, mesh_size, l_max 等
2. **数据结构尺寸校验**：mesh.r/rab, beta, chi, dion, rhoAtom 长度与 header 一致
3. **波函数归一化校验**：对每个束缚态验证 `∑ chi[i]² · rab[i] ≈ 1.0`（容差 1e-5）

测试数据 `Ge-spd-high.PD04.PBE.UPF`（ONCVPSP 生成）包含 5 个波函数（3S, 3P, 3D, 4S, 4P）和 6 个 beta 投影，全部通过校验。
