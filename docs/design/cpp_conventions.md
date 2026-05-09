# C++ 编码约定

## 命名规范

| 类型 | 命名规范 | 示例 |
|------|----------|------|
| 常量 | 全大写，下划线分隔 | `BOHR_RADIUS_ANGSTROM` |
| 类 | 大驼峰（PascalCase） | `GKK`, `GKKBufferIterator` |
| 结构体 | 大驼峰（PascalCase） | `GKKMetadata`, `KPointGVecs` |
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
    auto loadKPoint(int ikpt) -> const KPointGVecs&;
    auto current_ikpt() const -> int { return current_ikpt_; }
    auto currentData() const -> const KPointGVecs& { return current_data_; }

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
