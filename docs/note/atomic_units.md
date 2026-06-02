# Hartree 原子单位制

本项目在内部计算中采用 **Hartree 原子单位制（atomic units, a.u.）**。

## 基本约定

在 Hartree 原子单位制下，取：

- 电子质量 $m_e = 1$
- 元电荷 $e = 1$
- 约化普朗克常数 $\hbar = 1$
- 玻尔半径 $a_0 = 1$
- 真空介电常数 $4\pi\varepsilon_0 = 1$

因此，能量的单位为 **Hartree（$E_h$）**，长度的单位为 **Bohr（$a_0$）**。其它导出单位均由此定义。

## 本项目中的单位使用

| 物理量 | 文件中的单位 | 内部计算单位 |
|--------|-------------|-------------|
| 长度（晶格矢量 `A`） | Å（Angstrom） | Bohr |
| 能量 | — | Hartree |
| 波矢 $|G+k|^2/2$ | — | Hartree（因为 $\hbar = 1, m_e = 1$） |
| 赝势（UPF 文件） | Ry（Rydberg） | 读取后在 `NCPP` 构造时自动转换为 Hartree |

### 赝势的 Rydberg → Hartree 转换

UPF 格式使用 Rydberg 原子单位制（$e^2 = 2$, $m = 1/2$, $\hbar = 1$），
而项目内部使用 Hartree 原子单位制。`NCPP` 从 `UPF` 构造时对以下能量量纲的字段做 $\div 2$ 转换：

- `total_psenergy`、`wfc_cutoff`、`rho_cutoff`（元数据）
- `vloc`（局域势）
- `dion` / `B`（非局域势 D 矩阵）

径向量（`r`、`rab`、`beta` 投影函数、`chi` 波函数）和电荷密度 `rho_atom` 不受影响。

## 单位转换

### 长度：Å ↔ Bohr

```cpp
constexpr double BOHR_RADIUS_ANGSTROM = 0.52917721067;  // 1 Bohr = 0.52917721067 Å
```

- **Angstrom → Bohr**：`value_bohr = value_angstrom / BOHR_RADIUS_ANGSTROM`
- **Bohr → Angstrom**：`value_angstrom = value_bohr * BOHR_RADIUS_ANGSTROM`

### 能量：Hartree ↔ eV

如需与实验数据对比，常用的换算关系为：

$$1 \ \text{Hartree} = 27.211386245988 \ \text{eV}$$

## 代码实践

晶格矢量等带有物理单位的数据，应直接传入封装对象（如 `Lattice`），由对象内部根据传入的单位完成转换，而非在 IO 层手动转换：

```cpp
// Fortran 文件中 AL 以 Angstrom 存储
// 直接构造 Lattice，传入原始数值和单位，转换在 Lattice 内部完成
std::array<std::array<double, 3>, 3> al{};
for (int n = 0; n < 3; ++n) {
    for (int c = 0; c < 3; ++c) {
        al[n][c] = al_flat[n * 3 + c];
    }
}
Lattice lattice(al, LengthUnit::Angstrom);
```

内部所有公式（如有效质量计算、倒格子体积、动能项等）均默认使用原子单位制，不再出现 $\hbar$ 和 $m_e$。
