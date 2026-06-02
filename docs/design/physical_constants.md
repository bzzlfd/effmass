# 物理常量与单位常量

## 禁止使用 Magic Number

代码中**禁止**直接出现无名的裸数值常量（magic number）。所有物理量、单位转换系数、文件格式相关的固定数值都应当定义为具有**语义化名称**的常量。

## 使用 `constexpr` 定义常量

在 C++ 中，所有编译期可确定的常量应使用 `constexpr` 定义：

```cpp
constexpr double BOHR_RADIUS_ANGSTROM = 0.52917721067;  // 1 Bohr = 0.52917721067 Å
```

### 命名规范

- 常量名称使用**全大写**、**下划线分隔**
- 名称应清晰表达其物理意义和单位

## 本项目已定义的常量

| 常量名 | 数值 | 含义 |
|--------|------|------|
| `BOHR_RADIUS_ANGSTROM` | `0.52917721067` | 玻尔半径（单位：Å） |
| `HARTREE_TO_EV` | `27.211386245988` | 1 Hartree 对应的 eV 值（CODATA 2018） |
| `HARTREE_TO_RYDBERG` | `2.0` | 1 Hartree = 2 Rydberg |

## 单位转换函数

通过 `EnergyUnit` 枚举和 `convertEnergy` 重载实现类型安全的单位转换：

```cpp
enum class EnergyUnit { Hartree, Rydberg, eV };

// 返回纯转换因子：  x [from] * factor = x [to]
constexpr auto convertEnergy(EnergyUnit from, EnergyUnit to) -> double;

// 直接对数值做转换
constexpr auto convertEnergy(double value, EnergyUnit from, EnergyUnit to) -> double;
```

用法示例：

```cpp
// UPF 文件的数据是 Rydberg 单位，内部用 Hartree
meta.total_psenergy = convertEnergy(h.total_psenergy, EnergyUnit::Rydberg, EnergyUnit::Hartree);
```

### 与现有常量的关系

`RYDBERG_TO_EV` 可自行推导：`convertEnergy(1.0, EnergyUnit::Rydberg, EnergyUnit::eV)`。
`HARTREE_TO_RYDBERG` 是基准比值，所有单位转换均基于此常数，保证一致性。

## 单位转换原则

- 所有从输入文件读取的以 Å 为单位的长度量，应在读取时立即转换为 Bohr
- 转换公式：`value_bohr = value_angstrom / BOHR_RADIUS_ANGSTROM`
- 内部计算统一使用原子单位制，不再出现 $\hbar$ 和 $m_e$
