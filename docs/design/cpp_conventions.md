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
    explicit GKK(const std::string& filename);
    ~GKK();

    auto operator=(GKK&& other) noexcept -> GKK&;
    auto metadata() const -> const GKKMetadata& { return meta_; }
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
    explicit GKK(const std::string& filename);
    ~GKK();

    auto metadata() const -> const GKKMetadata&;
    auto loadKPoint(int ikpt) -> const KVecs&;

private:
    auto readMetadata() -> void;
    auto seekToKPoint(int ikpt) -> void;

    std::FILE* fp_;
    GKKMetadata meta_;
};
```
