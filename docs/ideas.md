# 潜在改进方向

## Hamiltonian 模块拆分方案对比

### 背景

`H_psi/hamiltonian.cppm` 包含 Hamiltonian 类及其三个嵌套的 Callable 类型（`Callable`、`Gradient::Callable`、`Hessian::Callable`）。待计算逻辑完整实现后，这个文件会膨胀到难以管理。有两种拆分方案。

### 方案一：平铺类 + 独立模块

类成为顶层类，各自一个 `.cppm` 文件，通过 `friend class Hamiltonian` 获得私有成员访问权限。

```
src/H_psi/callable.cppm   → export module H_psi.callable;   → class HamiltonianCallable
src/H_psi/gradient.cppm   → export module H_psi.gradient;   → class HamiltonianGradient
src/H_psi/hessian.cppm    → export module H_psi.hessian;    → class HamiltonianHessian
src/H_psi/hamiltonian.cppm → export module H_psi.hamiltonian; → class Hamiltonian（import 上述三个模块）
```

**依赖关系**：单向，`callable/gradient/hessian` 只需前置声明 `class Hamiltonian;`（指针成员），`hamiltonian` `import` 它们获取完整类型。

**关键代码模式**：
```cpp
// callable.cppm
export module H_psi.callable;
class Hamiltonian;  // 前置声明

export class HamiltonianCallable {
    void operator()(...) const;
private:
    friend class Hamiltonian;
    HamiltonianCallable(const Hamiltonian* parent, int ikpt);
};
```

### 方案二：嵌套类 + .cpp 实现分离

保留 `Hamiltonian::Callable` / `Gradient::Callable` / `Hessian::Callable` 的嵌套结构，.cppm 只放声明，实现写入 `.cpp` 模块实现文件（`module H_psi.hamiltonian;`）。

```
src/H_psi/hamiltonian.cppm               → 接口（只有声明）
src/H_psi/hamiltonian_callable.cpp       → void Hamiltonian::Callable::operator()(...) { ... }
src/H_psi/hamiltonian_gradient.cpp       → void Hamiltonian::Gradient::Callable::operator()(...) { ... }
src/H_psi/hamiltonian_hessian.cpp        → void Hamiltonian::Hessian::Callable::operator()(...) { ... }
```

**关键代码模式**：
```cpp
// hamiltonian.cppm  — 接口
export class Hamiltonian {
public:
    class Callable {
        void operator()(...) const;  // 只有声明
    private:
        friend class Hamiltonian;
        Callable(const Hamiltonian* parent, int ikpt);
    };
};
```

```cpp
// hamiltonian_callable.cpp  — 实现
module H_psi.hamiltonian;
void Hamiltonian::Callable::operator()(...) const { /* 逻辑 */ }
```

### 对比总览

| 维度 | 方案一（平铺类 + 独立模块） | 方案二（嵌套类 + .cpp 实现） |
|------|---------------------------|---------------------------|
| 文件数 | 5 `.cppm` + 可选 `.cpp` | 1 `.cppm` + 3 `.cpp` |
| 命名 | `HamiltonianCallable` | `Hamiltonian::Callable` |
| friend | 每平铺类 1 个 `friend class Hamiltonian` | 不需要 |
| 实现文件 import | 每个 `.cpp` 需手动 import 依赖 | 继承 `.cppm` 的 import，无需重复 |
| 增量编译 | 改 `.cppm` 触发模块依赖链重编 | 改 `.cpp` 只编那个文件 |
| 模块独立性 | 各模块可被单独 import | 必须通过 `Hamiltonian` 使用 |
| 语义分组 | 靠命名约定（`Hamiltonian` 前缀） | 靠 C++ 嵌套作用域 |
| | | 嵌套类语法练习 |

### 结论

选择方案二。因为 Callable 对象都依赖 Hamiltonian 内的数据，方案一给人一种可以单独使用 Callable 的错觉，但实际上没有意义。选择二用 `Hamiltonian::Callable` 的嵌套作用域比 `HamiltonianCallable` 的前缀约定更能避免歧义。