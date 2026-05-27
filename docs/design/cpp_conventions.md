# C++ 编码约定

## 命名规范

| 类型 | 命名规范 | 示例 |
|------|----------|------|
| 常量 | 全大写，下划线分隔 | `BOHR_RADIUS_ANGSTROM` |
| 类 | 大驼峰（PascalCase） | `GKK`, `GKKBufferIterator` |
| 结构体 | 大驼峰（PascalCase） | `GKKMetadata`, `KVecs` |
| 成员变量 | 后缀下划线 | `meta_`, `fp_`, `current_kpt_` |
| 函数 | 小驼峰（camelCase） | `loadKPoint()`, `readMetadata()` |
| 局部变量 | 小驼峰 | `record_len`, `ikpt` |

## 函数返回类型写法

所有非构造函数/析构函数的返回类型，统一采用 **尾随返回类型（trailing return type）** 写法：

```cpp
auto function_name(args...) -> ReturnType;
```

### 目的

- **对齐**：当函数名长度不一时，返回类型放在右侧可以使代码左侧对齐，提升可读性。
- **一致性**：项目内所有函数采用同一种写法，降低阅读成本。
- **与现代 C++ 风格一致**：`auto` + 尾随返回类型是 C++11 以来的标准语法，在模板和模块化代码中尤其自然。

### 示例

```cpp
export class GKK {
public:
    GKKMetadata meta;

    explicit GKK(const std::string& filename);
    ~GKK();

    auto operator=(GKK&& other) noexcept -> GKK&;
    auto loadKPoint(int ikpt) -> const KVecs&;
    auto current_ikpt() const -> int { return current_ikpt_; }
    auto currentData() const -> const KVecs& { return current_data_; }

private:
    auto readRecordLength() -> int;
    auto checkRecordLength(int expected) -> void;
    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void;
    auto readMetadata() -> void;
    auto seekToKPoint(int ikpt) -> void;

    std::FILE* fp_;
    GKKMetadata meta_;
    int current_ikpt_ = -1;
};
```

### 例外

- **构造函数** 和 **析构函数** 没有返回类型，不适用此写法。
- **运算符重载**（如 `operator=`）适用此写法，因为运算符本质上也是函数。

## 空行规范

借鉴 PEP 8 的空行哲学。

### 1. 导入区与正文

文件顶部的 `module;`、`#include`、`export module`、`import` 等声明属于导入区。导入区与第一个正文（类、结构体、函数、全局常量）之间使用 **两个空行**。

导入区内部：
- `module;` / `#include` 区 与 `export module` 声明之间：**一个空行**
- `export module` 声明 与 `import` 语句之间：**一个空行**
- `import` 语句按来源分组（标准库 → 第三方 → 本项目），组间用**一个空行**

```cpp
module;
#include <cstdio>
#include "pugixml.hpp"

export module io.GKK;

import io.lattice;
import std;


// 正文开始
export struct GKKMetadata {
    ...
};
```

### 2. 顶级定义之间

文件级别的类、结构体、枚举、函数定义之间使用 **两个空行**。

```cpp
export struct KVecs {
    ...
};


export class GKK {
    ...
};


auto foo() -> int {
    ...
}
```

**例外**：一组紧密相关的短内联函数（如运算符重载、简单 getter）之间可只用一个空行或不空行。

```cpp
export constexpr auto operator|(KVecsView a, KVecsView b) -> KVecsView {
    ...
}
export constexpr auto operator&(KVecsView a, KVecsView b) -> KVecsView {
    ...
}
```

### 3. 类/结构体内部

类或结构体内部的方法声明或定义之间使用 **一个空行**。数据成员通常不空行；若按逻辑分组，组间可空一行。

```cpp
export class GKK {
public:
    GKKMetadata meta;

    explicit GKK(const std::string& filename);
    ~GKK();

    auto loadKPoint(int ikpt) -> const KVecs&;

private:
    auto readMetadata() -> void;
    auto seekToKPoint(int ikpt) -> void;

    std::FILE* fp_;
```

## 导出组织规范

所有被导出的名称统一放在文件开头的一个 `export { }` 块中（位于 `import` 声明之后），以便一眼看清模块的公共 API。

```cpp
export module my.module;

import std;

export {
    struct PublicType { ... };
    enum class PublicEnum { ... };
    class PublicClass { ... };
    class PublicClass2;
    auto publicFunction() -> void;
    using Alias = ExistingType;
}

// 实现细节不加 export 关键字，放在 export 块之后
auto publicFunction() -> void { ... }
PublicClass::PublicClass() { ... }

class PublicClass2 {
public: 
    ...
private: 
    ...
}

```

- **导出的类型/函数**：放入 `export { }` 块。
- **普通实现细节**（辅助函数、常量、内部类型）：不加 `export` 关键字，放在 `export` 块之后。
- **缩进体现层次**：`export { }` 块内部通过缩进表达类型之间的从属关系——主类型顶格，其拥有的子类型缩进一级，子类型的子类型再缩进一级。
- **namespace**: 一种额外的隔离方式
    - **`namespace archived`**：仅为展示另一种可能的实现方案，保留以备参考。
    - **`namespace deprecated`**：已弃用的代码，保留以兼容旧版或过渡期使用。

### 缩进示例

```cpp
export {
    class NCPPUPF;
        struct NCPPUPFHeader;
        struct NCPPUPFMesh;
            enum class MeshType : int;
        struct NCPPUPFNonlocal;
            struct NCPPUPFNonlocalByL;
        struct NCPPUPFWavefunction;
}
```

### 示例

```cpp
module;

export module math.sph_bessel;

import std;

// 普通实现细节，不加 export
constexpr double BESSEL_EPS = 1e-15;

auto sphericalBesselJImpl(int l, double x) -> double { ... }

export {
    auto sphericalBesselJ(int l, double x) -> double;
    auto sphericalBesselJ(int l, std::span<const double> xs, std::span<double> out) -> void;
    class SphericalBesselJ { ... };
}


auto sphericalBesselJ(int l, double x) -> double { ... }

auto sphericalBesselJ(int l, std::span<const double> xs, std::span<double> out) -> void { ... }

// 仅为参考的另一种实现
namespace archived {
    auto sphericalBesselJ_Taylor(int l, double x) -> double { ... }
}

// 已弃用的接口
namespace deprecated {
    auto oldBesselInterface(int l, double x) -> double;
}
```
