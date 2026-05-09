# Pseudo-Potential IO Modules

本文档描述 `src/pseudo/` 下赝势文件读取模块的设计与接口。

## 模块结构

```
src/pseudo.cppm          # 总集成模块: export module pseudo;
src/pseudo/
  ├── ncpp-upf.cppm      # NCPP UPF v.2 读取 (pseudo.ncpp_upf)
  ├── uspp-upf.cppm      # USPP 骨架 (pseudo.uspp_upf) — 未实现
  └── paw-upf.cppm       # PAW 骨架 (pseudo.paw_upf) — 未实现
```

`pseudo.cppm` 仅负责重新导出子模块，本身不含实现。

## NCPPUPF — Norm-Conserving Pseudo-Potential (UPF v.2)

### 参考资料

- [Quantum ESPRESSO: Unified Pseudopotential Format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
- Vanderbilt, D. (1990). *Soft self-consistent pseudopotentials in a generalized eigenvalue formalism*. Phys. Rev. B, 41, 7892. DOI: [10.1103/PhysRevB.41.7892](https://doi.org/10.1103/PhysRevB.41.7892)

### 公共接口

```cpp
export enum class MeshType { Uniform, Exponential, Unknown };

export struct NCPPUPFNonlocalByL {
    std::vector<std::vector<double>> beta;  // [n_beta_l][mesh]
    Matrix dion;                             // D_ij submatrix for this l
};

export class NCPPUPF {
public:
    explicit NCPPUPF(const std::string& filename);

    auto header() const -> const NCPPUPFHeader&;
    auto mesh() const -> const NCPPUPFMesh&;
    auto localPotential() const -> std::span<const double>;
    auto nonlocal() const -> const NCPPUPFNonlocal&;
    auto wavefunctions() const -> const NCPPUPFWavefunction&;
    auto rhoAtom() const -> std::span<const double>;

    auto meshType() const -> MeshType;
    auto nonlocalByL(int l) const -> NCPPUPFNonlocalByL;
};
```

### 数据结构与 UPF 文件 tag 的对应关系

```
NCPPUPF
├── NCPPUPFHeader         ← PP_HEADER (全部属性)
├── NCPPUPFMesh           ← PP_MESH
│   ├── r[]               ← PP_R
│   └── rab[]             ← PP_RAB
├── MeshType              ← 由 r/rab 推断（Uniform / Exponential / Unknown）
├── vloc_                 ← PP_LOCAL
├── NCPPUPFNonlocal       ← PP_NONLOCAL
│   ├── beta[][]          ← PP_BETA.1 ... PP_BETA.nbeta（截断至 kbeta）
│   ├── lll[]             ← PP_BETA.* / angular_momentum
│   ├── kbeta[]           ← PP_BETA.* / cutoff_radius_index
│   ├── rcut[]            ← PP_BETA.* / cutoff_radius
│   └── dion              ← PP_DIJ (Matrix 类型)
├── NCPPUPFNonlocalByL    ← 按角动量 l 筛选的 beta + dion 子矩阵
├── NCPPUPFWavefunction   ← PP_PSWFC
│   ├── chi[][]           ← PP_CHI.1 ... PP_CHI.nwfc（截断尾部零）
│   ├── kchi[]            ← 有效长度（最后一个非零元素位置）
│   ├── lchi[]            ← PP_CHI.* / l
│   ├── oc[]              ← PP_CHI.* / occupation
│   └── labels[]          ← PP_CHI.* / label
└── rho_at_               ← PP_RHOATOM
```

#### 扁平存储 vs. 嵌套容器

- **`beta`** 与 **`chi`**：使用 `std::vector<std::vector<double>>`。读取后按实际有效长度截断尾部零，因此每行长度 ≤ `mesh_size`，该尺度（mesh ~ 10³，nbeta ~ 10）下的间接寻址开销可忽略。
- **`dion`**：使用自定义的 `Matrix` 值类型。$D_{ij}$ 是 $n_{\beta} \times n_{\beta}$ 对称矩阵，扁平连续存储更适合矩阵运算，且 C++23 多维 `operator[](i, j)` 提供了自然访问语法。
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
  ├── <PP_HEADER/>          → 读取全部属性 → NCPPUPFHeader
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

### 构造函数一次性全加载

UPF 文件通常只有几百 KB，pugixml 的 DOM 解析完全无压力。构造函数一次性读取全部数据到内存后即关闭文件句柄，类内部只持有 `std::vector`、`std::string` 等标准值类型成员。

因此 `NCPPUPF` 遵循 **Rule of Zero**——不显式声明析构函数、拷贝/移动构造函数或赋值运算符，编译器自动生成的语义完全正确。这也意味着该类可以自然地放入 `std::vector` 等标准容器中。

### 程序中不读取 `columns` 的原因

UPF v.2 的部分数据节点带有 `columns` 属性，例如：

```xml
<PP_R columns="4">
   1.0  2.0  3.0  4.0
   5.0  6.0  7.0  8.0
   ...
</PP_R>
```

在 `NCPPUPF` 的实现中，`parseTextToVector` 通过 `std::istringstream` 逐 token 读取，按空白字符（空格、制表符、换行）分割。因此无论每行有多少列，所有数值都会被正确提取。

**不读取 `columns` 的设计理由**：

1. **该属性是输出格式提示，不是语义约束**：`columns` 仅影响人类可读性和 Fortran `write` 语句的格式化输出，不改变数据本身的顺序和数量。

2. **pugixml 的文本节点已抹除格式信息**：`pugi::xml_node::child_value()` 返回的是节点内所有文本内容的拼接字符串，行边界信息已经丢失。即使想验证 `columns`，也无法从纯文本中可靠还原原始行结构。

3. **数值校验已足够**：每个数组读取后都会与 `mesh_size`（或 `nbeta`、`nwfc`）进行长度校验。如果格式化错误导致数据缺失或多余，`size mismatch` 异常会立即暴露问题。



### Mesh 类型推断 (`meshType()`)

`meshType()` 通过分析 `mesh.r` 和 `mesh.rab` 的数值关系推断网格类型：

- **Uniform（均匀）**：`r[i] - r[i-1]` 在所有网格点上的偏差不超过均值的 `1e-6`（相对容差）。
- **Exponential（指数/对数）**：`rab[i] / r[i]` 在所有非零 `r[i]` 上的偏差不超过均值的 `1e-6`。对于标准的对数网格 `r[i] = r_0 * exp((i-1) * dx)`，积分权重 `rab[i] = r[i] * dx`，因此该比值严格为常数。
- **Unknown**：不满足以上两种模式，或网格点少于 2 个。

推断基于数值容差，不依赖 UPF 文件中的任何额外属性。

### 按角动量筛选非局域项 (`nonlocalByL(int l)`)

`nonlocalByL(l)` 根据 `NCPPUPFNonlocal::lll` 筛选出角动量等于 `l` 的所有 beta 投影，并提取对应的 `D_ij` 子矩阵。

- 返回类型为 `NCPPUPFNonlocalByL`，仅包含 `beta`（`std::vector<std::vector<double>>`）和 `dion`（`Matrix`）。
- 若不存在该 `l` 的 beta，返回空 `beta` 列表和 0×0 的 `dion`。
- `beta` 行按原始 `nonlocal.beta` 中出现的顺序排列，`dion` 的 `[i, j]` 对应原始矩阵中 `[global_i, global_j]` 的值。

### 错误处理

| 场景 | 异常类型 | 消息 |
|------|---------|------|
| XML 解析失败 | `std::runtime_error` | `UPF: failed to parse file: ...` |
| 缺少根节点 `<UPF>` | `std::runtime_error` | `UPF: missing root element <UPF>` |
| 缺少属性 | `std::runtime_error` | `UPF: missing attribute '...'` |
| 属性类型转换失败 | `std::runtime_error` | `UPF: invalid double/int/bool attribute ...` |
| 数组大小不匹配 | `std::runtime_error` | `UPF: <tag> size mismatch` |

## USPP 与 PAW（骨架）

`USPP` 和 `PAW` 类已创建但尚未实现。构造函数直接抛出：

```cpp
throw std::runtime_error("USPP/PAW reader not yet implemented");
```

两者与 `NCPPUPF` 保持相同的接口风格（纯值类型、Rule of Zero、提供 `header()` 接口），为后续扩展预留统一的入口。

## 测试

```
test/
├── CMakeLists.txt
test/test_ncpp_upf/
    ├── test_ncpp_upf.cpp
    └── Ge-spd-high.PD04.PBE.UPF
```

`test_ncpp_upf` 验证内容：
1. **Header 字段校验**：element, pseudo_type, z_valence, mesh_size, l_max 等
2. **数据结构尺寸校验**：mesh.r/rab, beta, chi, dion, rhoAtom 长度与 header 一致
3. **波函数归一化校验**：对每个束缚态验证 `∑ chi[i]² · rab[i] ≈ 1.0`（容差 1e-5）

测试数据 `Ge-spd-high.PD04.PBE.UPF`（ONCVPSP 生成）包含 5 个波函数（3S, 3P, 3D, 4S, 4P）和 6 个 beta 投影，全部通过校验。
