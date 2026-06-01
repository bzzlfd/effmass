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
    auto socData() const -> const UPFSOC*;
};
```

### 数据结构与 UPF 文件 tag 的对应关系

```
UPF
├── UPFHeader             ← PP_HEADER (全部属性)
│   ├── has_so            ← 是否存在 SOC 扩展数据
│   ├── pseudo_type       ← "NC" | "US" | "PAW"
│   ├── is_ultrasoft      ← USPP 标识
│   └── is_paw            ← PAW 标识
├── UPFMesh               ← PP_MESH
│   ├── r[]               ← PP_R
│   └── rab[]             ← PP_RAB
├── vloc_                 ← PP_LOCAL
├── UPFNonlocal           ← PP_NONLOCAL
│   ├── beta[][]          ← PP_BETA.1 ... PP_BETA.nbeta（完整 mesh_size 长度）
│   ├── lll[]             ← PP_BETA.* / angular_momentum
│   ├── kbeta[]           ← PP_BETA.* / cutoff_radius_index
│   ├── rcut[]            ← PP_BETA.* / cutoff_radius
│   └── dion              ← PP_DIJ (DenseMatrix<double>)
├── UPFWavefunction       ← PP_PSWFC
│   ├── chi[][]           ← PP_CHI.1 ... PP_CHI.nwfc（完整 mesh_size 长度）
│   ├── lchi[]            ← PP_CHI.* / l
│   ├── oc[]              ← PP_CHI.* / occupation
│   └── labels[]          ← PP_CHI.* / label
├── UPFSOC                ← PP_SPIN_ORB（扩展，仅 has_so 时存在）
│   ├── jjj[]             ← PP_RELBETA.* / jjj（PseudoDojo）或 PP_BETA.* / jjj（PSlibrary）
│   └── jchi[]            ← PP_RELWFC.* / jchi（PseudoDojo）或 PP_CHI.* / jchi（PSlibrary）
└── rho_at_               ← PP_RHOATOM
```

**NC / US / PAW 支持现状**：
- **NC**（Norm-Conserving）：完全实现，读取全部字段
- **US**（Ultrasoft）：读取 header 后抛出 `"UPF: USPP reader not yet implemented"`，结构体已预留，读取逻辑待实现
- **PAW**（Projector-Augmented Wave）：读取 header 后抛出 `"UPF: PAW reader not yet implemented"`，结构体已预留，读取逻辑待实现

构造函数先读取 PP_HEADER，根据 `is_ultrasoft` / `is_paw` 标识判断类型。对于 NC 类型，继续读取剩余全部字段；对于 US/PAW 类型，因 PP_NONLOCAL 中额外结构（Q_ij 增广电荷等）尚未解析，暂不支持完整读取。

#### 扁平存储 vs. 嵌套容器

- **`beta`** 与 **`chi`**：使用 `std::vector<std::vector<double>>`。每行长度为完整的 `mesh_size`，不做截断。IO 层的职责是忠实地映射文件内容，截断优化属于算符层。
- **`dion`**：使用 `DenseMatrix<double>` 值类型（来自 `utils.matrix`）。$D_{ij}$ 是 $n_{\beta} \times n_{\beta}$ 对称矩阵，扁平连续存储更适合矩阵运算，且 C++23 多维 `operator[](i, j)` 提供了自然访问语法。
- **`vloc_`** 与 **`rho_at_`**：使用 `std::vector<double>` 并通过 `std::span<const double>` 暴露，方便调用方直接遍历而无需深拷贝。

### IO 层的数据截断原则

UPF 文件中的 `beta` 与 `chi` 数组在 `cutoff_radius_index` 之后通常全部为零——但 `upf.cppm` **不做截断**。

**设计理由**：IO 层的职责是忠实映射文件内容。`upf.cppm` 存储的 `beta` 和 `chi` 为完整 `mesh_size` 长度，与文件中的数值一一对应。截断（`resize` 到有效长度）是性能优化决策，放在算符层 `ncpp.cppm` 的构造函数中执行：

- **beta**：依据 `cutoff_radius_index` 属性（文件属性，仍存储于 `UPFNonlocal::kbeta`）
- **chi**：从末尾向前扫描，去除精确为 `0.0` 的尾部元素，有效长度存入 `PsWFc::kchi`

截断后的索引仍与 `mesh.r`、`mesh.rab` 一一对应，只是循环上限变小。

**ncpp.cppm 中的实现**：

```cpp
// beta 截断：依据 cutoff_radius_index（文件属性）
for (int i = 0; i < header.number_of_proj; ++i) {
    nonlocal.beta[i].resize(nonlocal.cutoff_index[i]);
}

// chi 截断：扫描尾部精确零，有效长度记入 kchi
for (int i = 0; i < header.number_of_wfc; ++i) {
    int cutoff = wfc.chi[i].size();
    while (cutoff > 0 && wfc.chi[i][cutoff - 1] == 0.0) --cutoff;
    wfc.chi[i].resize(cutoff);
    wfc.kchi[i] = cutoff;
}
```

其中 `PsWFc::kchi` 是算符层自建的有效长度字段（UPF 文件的 `PP_CHI.*` 没有 `cutoff_radius_index` 属性），与 `Nonlocal::cutoff_index` 形成对称。

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
  ├── <PP_RHOATOM>          → 解析文本 → rho_at_[mesh]
  └── <PP_SPIN_ORB>?        → 仅 has_so 时读取（PseudoDojo 格式）
                             或从 PP_BETA.* / PP_CHI.* 读取 jjj/jchi（PSlibrary 格式）
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
| 缺少 PP_SPIN_ORB 子标签（SOC PseudoDojo） | `std::runtime_error` | `UPF: missing <PP_RELBETA.i> in PP_SPIN_ORB` |
| 缺少 PP_BETA.* 的 jjj 属性（SOC PSlibrary） | `std::runtime_error` | `UPF: missing jjj attribute on <PP_BETA.i>` |
| 缺少 PP_CHI.* 的 jchi 属性（SOC PSlibrary） | `std::runtime_error` | `UPF: missing jchi attribute on <PP_CHI.i>` |


### 自旋轨道耦合（SOC）支持

#### 设计原则

`upf.cppm` 的设计目标是忠实地表示标准 UPF v.2 XML 格式。然而 SOC 扩展数据（`jjj`、`jchi`）并非标准 UPF v.2 规范的一部分——它们来源于 **PP_SPIN_ORB** 扩展段（PseudoDojo 约定）或作为非标准属性嵌入标准标签中（PSlibrary 约定）。

因此，所有 SOC 相关字段被隔离到独立的 `UPFSOC` 结构体中，不污染核心数据结构。读者可以从核心类型（`UPFHeader`、`UPFMesh`、`UPFNonlocal`、`UPFWavefunction`）中了解标准 UPF 结构，而不被非标准的 SOC 约定干扰。

`angular_momentum`（`lll`）字段在 `UPFNonlocal` 中按原样存储，不做任何语义解释。


##### 错误做法例子

PWmat 的 upfSO 格式将 l=1/2/3 的自旋轨道投影子 β_SO 编码为 `angular_momentum=10/20/30`。
这是一种对 SOC 的支持的处理方式，说实话这很优雅，精简数据、使得非 SOC 代码使用非SOC UPF 数据不会不会察觉变化。
可是我们**不应该在 UPF 对象中注解 `angular_momentum=10/20/30` 的含义**：
这只是 UPF 文件的 SOC 扩展的一种方案，但是现在却污染了整个 UPF 对象的定义。这让新读者产生困惑。

1. 任何依赖 `angular_momentum == l` 的通用代码（如按 l 分组投影子）会被静默破坏
2. IO 层与解释层的边界模糊，换一种 SOC 方案就得改核心读取逻辑

这正是 `upf.cppm` **不做**的事情。IO 层遇到 10/20/30 只是照存，上层（如未来的 `upfso.cppm`）自行识别。这项原则与 `UPFSOC` 独立于 `UPFNonlocal` 一脉相承：**文件里有什么就存什么，怎么用留给解释层**。

#### UPFSOC 结构体

```cpp
struct UPFSOC {
    enum class Format { pp_spin_orb, ps_library };

    Format format = Format::pp_spin_orb;
    std::vector<double> jjj;   // [nbeta] 总角动量 j
    std::vector<double> jchi;  // [nwfc]  总角动量 j
};
```

- `format` 记录检测到的 SOC 数据来源格式
- `jjj` 与 beta 投影一一对应（`UPFNonlocal::lll` 平行结构）
- `jchi` 与波函数一一对应（`UPFWavefunction::lchi` 平行结构）

#### 两种 SOC 数据格式

`readSpinOrbit()` 内部支持两种格式的自动检测：

1. **PseudoDojo 格式（`pp_spin_orb`）**：`<PP_SPIN_ORB>` 段包含自闭合标签 `<PP_RELBETA.* jjj="..."/>` 和 `<PP_RELWFC.* jchi="..."/>`，分别存储对应 beta 投影和波函数的 j 值。

2. **PSlibrary 格式（`ps_library`）**：jjj 作为 `<PP_BETA.*>` 标签上的额外属性，jchi 作为 `<PP_CHI.*>` 标签上的额外属性，均嵌入标准 `PP_NONLOCAL` / `PP_PSWFC` 段中。

**检测逻辑**：优先尝试查找 `<PP_SPIN_ORB>` 段；若不存在，回退到从 `PP_BETA.*` 和 `PP_CHI.*` 标签的属性中读取。检测结果记录在 `UPFSOC::format` 中。

#### 访问方式

```cpp
auto upf = UPF(filename);
const UPFSOC* soc = upf.socData();
if (soc) {
    // 处理 SOC 数据
    auto j_values = soc->jjj;
} else {
    // 非 SOC 赝势，无需处理
}
```

`socData()` 返回指针而非引用：当 `UPFHeader::has_so == false` 时返回 `nullptr`，调用方无需检查 `has_so` 即可安全区分 SOC / 非 SOC 文件。

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
