# UPF (Unified Pseudopotential Format) v.2 解析笔记

## 参考资料

- [Quantum ESPRESSO: Unified Pseudopotential Format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
- Vanderbilt, D. (1990). *Soft self-consistent pseudopotentials in a generalized eigenvalue formalism*. Phys. Rev. B, 41, 7892. DOI: [10.1103/PhysRevB.41.7892](https://doi.org/10.1103/PhysRevB.41.7892)

## 单位制

UPF 文件采用 **Rydberg 原子单位制**：
- $e^2 = 2$
- $m = 1/2$
- $\hbar = 1$
- 长度单位：Bohr ($0.529177\,\text{Å}$)
- 能量单位：Ry ($13.6058\,\text{eV}$)
- 势能已乘以 $e$，故与能量同量纲

> **注意**：某些在径向网格上计算的量含有额外的 $r$ 因子，详见各字段说明。

---

## 文件结构

```xml
<UPF version="2.0.1">
  <PP_INFO>   ...  </PP_INFO>      <!-- 生成信息，人类可读，可忽略 -->
  <PP_HEADER attr1="..." ... />    <!-- 元数据头 -->
  <PP_MESH>   ...  </PP_MESH>      <!-- 径向网格 -->
  <PP_LOCAL>  ...  </PP_LOCAL>     <!-- 局域势 -->
  <PP_NONLOCAL> ... </PP_NONLOCAL> <!-- 非局域势 (beta + D_ij) -->
  <PP_PSWFC>  ...  </PP_PSWFC>     <!-- 赝波函数 -->
  <PP_RHOATOM> ... </PP_RHOATOM>   <!-- 赝原子电荷密度 -->
</UPF>
```

以下按读取顺序解释 `NCPPUPF` 目前解析的各字段。

---

## PP_HEADER

```xml
<PP_HEADER
  generated="..."
  author="..."
  date="..."
  comment="..."
  element="Ge"
  pseudo_type="NC"
  relativistic="scalar"
  is_ultrasoft="F"
  is_paw="F"
  is_coulomb="F"
  has_so="F"
  has_wfc="F"
  has_gipaw="F"
  core_correction="T"
  functional="..."
  z_valence="22.0"
  total_psenergy="..."
  wfc_cutoff="..."
  rho_cutoff="..."
  l_max="2"
  l_local="2"
  mesh_size="1560"
  number_of_wfc="5"
  number_of_proj="6"
/>
```

| 属性 | 类型 | 含义 |
|------|------|------|
| `generated` | string | 生成该赝势的代码名称 |
| `author` | string | 作者 |
| `date` | string | 生成日期 |
| `comment` | string | 描述 |
| `element` | string | 元素化学符号 |
| `pseudo_type` | string | `NC` = 模守恒 (Norm-Conserving)；`SL` = 半局域；`US` = 超软；`PAW` = PAW；`1/r` = Coulomb |
| `relativistic` | string | `scalar` / `full` / `nonrelativistic` |
| `is_ultrasoft` | bool | 是否为超软赝势 |
| `is_paw` | bool | 是否为 PAW 数据集 |
| `is_coulomb` | bool | 是否为 Coulomb 势 (测试用，无局域/非局域项) |
| `has_so` | bool | 是否包含自旋-轨道耦合项 |
| `has_wfc` | bool | 是否包含全电子波函数 (PP_FULL_WFC) |
| `has_gipaw` | bool | 是否包含 GIPAW 重构数据 |
| `core_correction` | bool | 是否有非线性芯修正 (NLCC) |
| `functional` | string | 交换关联泛函标识符 |
| `z_valence` | double | 价电子电荷数 $Z_{\text{val}}$ |
| `total_psenergy` | double | 赝价电子总能量 (Ry) |
| `wfc_cutoff` | double | 建议的波函数平面波截断 (Ry) |
| `rho_cutoff` | double | 建议的电荷密度平面波截断 (Ry) |
| `l_max` | int | 赝势中最大的角动量分量 $l_{\text{max}}$ |
| `l_local` | int | 局域势对应的角动量 $l_{\text{loc}}$ |
| `mesh_size` | int | 径向网格点数 |
| `number_of_wfc` | int | `PP_PSWFC` 中赝波函数的个数 $n_{\text{wfc}}$ |
| `number_of_proj` | int | Kleinman-Bylander 投影函数 ($\beta$) 的个数 $n_{\text{beta}}$ |

---

## PP_MESH

```xml
<PP_MESH dx="..." mesh="..." xmin="..." rmax="..." zmesh="...">
  <PP_R>
    r(1) r(2) ... r(mesh)
  </PP_R>
  <PP_RAB>
    rab(1) rab(2) ... rab(mesh)
  </PP_RAB>
</PP_MESH>
```

| 变量 | 含义 |
|------|------|
| `r(i)` | 径向网格点。常用形式：$r(i) = e^{x_{\min}+i \cdot dx} / Z_{\text{mesh}}$ 或 $r(i) = (e^{x_{\min}+i \cdot dx} - 1) / Z_{\text{mesh}}$ |
| `rab(i)` | 离散积分权重因子，满足 $\int f(r)\,dr = \sum_i f_i \cdot \text{rab}_i$ |

> 注意：部分生成器（如 ONCVPSP）使用均匀线性网格，此时 `rab` 为常数。

---

## PP_LOCAL

```xml
<PP_LOCAL>
  vloc(1) vloc(2) ... vloc(mesh)
</PP_LOCAL>
```

- `vloc(mesh)`：局域赝势，在径向网格上采样 (Ry a.u.)。
- 对于 `is_coulomb == true` 的测试势，该字段可能不存在。

---

## PP_NONLOCAL

```xml
<PP_NONLOCAL>
  <PP_BETA.1 ...> ... </PP_BETA.1>
  ...
  <PP_BETA.nbeta ...> ... </PP_BETA.nbeta>
  <PP_DIJ> ... </PP_DIJ>
</PP_NONLOCAL>
```

### PP_BETA.*

| 属性 | 含义 |
|------|------|
| `angular_momentum` | 该投影函数的角动量 $l$ |
| `cutoff_radius_index` | 截断半径对应的网格索引 $k_{\beta}$ (≤ mesh) |
| `cutoff_radius` | 截断半径 $r_{\text{cut}}$ (a.u.) |

- `beta(i)`：投影函数 $\beta_i(r)$，**已乘以 $r$**（即文件中存储的是 $r \cdot \beta(r)$）。
- 单位：Bohr$^{-1/2}$。

### PP_DIJ

- `dion(i, j)`：Kleinman-Bylander 非局域势的 $D_{ij}$ 矩阵。
- 非局域势形式：$V_{\text{NL}} = \sum_{i,j} D_{ij} \, |\beta_i\rangle \langle \beta_j|$
- 单位：Ry。
- $D_{ij}$ 是**对称矩阵**（$D_{ij} = D_{ji}$）。

> **注意**：由于 `beta` 和 `dion` 在代码中仅以乘积形式 $(\beta \cdot D \cdot \beta)$ 出现，某些转换器可能以不同单位组合输出（例如 `dion` 用 Ry$^{-1}$ 而 `beta` 用 Ry$\cdot$Bohr$^{-1/2}$），实际计算时只要整体乘积单位正确即可。

---

## PP_PSWFC

```xml
<PP_PSWFC>
  <PP_CHI.1 l="..." occupation="..." label="...">
    chi(1,1) chi(2,1) ... chi(mesh,1)
  </PP_CHI.1>
  ...
</PP_PSWFC>
```

| 属性 | 含义 |
|------|------|
| `l` | 角动量量子数 |
| `occupation` | 该轨道的占据数 |
| `label` | 轨道标签（如 "3S", "3P" 等） |

- `chi(mesh, i)`：第 $i$ 个赝原子轨道的径向部分 $\chi_i(r)$，形式为 $\chi(r) = \sqrt{4\pi} \, r \cdot R(r)$，其中 $R(r)$ 是径向波函数。
- 束缚态波函数的归一化条件（见下节）。

---

## PP_RHOATOM

```xml
<PP_RHOATOM>
  rho_at(1) rho_at(2) ... rho_at(mesh)
</PP_RHOATOM>
```

- `rho_at(mesh)`：赝原子电荷密度，**已乘以 $4\pi r^2$**（即文件中存储的是 $4\pi r^2 \rho(r)$）。
- 积分后应等于价电子数 $Z_{\text{val}}$。

---

## 束缚态波函数的归一化条件
**UPF 中的 `chi`（赝波函数）和 以 $\chi(r) = \sqrt{4\pi} \, r \cdot R(r)$ 的形式存储**，其中 $R(r)$ 是径向波函数，即

束缚态归一化条件：
$$
\int_0^\infty |\chi(r)|^2 \, dr = 1
$$

`beta` 亦受此影响。

`test_ncpp_upf` 对 $\int_0^\infty |\chi(r)|^2 \, dr = 1$ 进行了验证。

