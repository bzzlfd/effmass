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
| 长度（晶格矢量 `AL`） | Å（Angstrom） | Bohr |
| 能量 | — | Hartree |
| 波矢 $|G+k|^2/2$ | — | Hartree（因为 $\hbar = 1, m_e = 1$） |

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

所有从输入文件读取的以 Å 为单位的长度量，应在读取时立即转换为 Bohr：

```cpp
// Fortran 文件中 AL 以 Angstrom 存储
// 读取后统一转换为 Bohr（原子单位制长度单位）
meta_.AL[n][c] = al_flat[n * 3 + c] / BOHR_RADIUS_ANGSTROM;
```

内部所有公式（如有效质量计算、倒格子体积、动能项等）均默认使用原子单位制，不再出现 $\hbar$ 和 $m_e$。
