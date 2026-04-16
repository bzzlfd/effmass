# C++ 命名规范

| 类型 | 命名规范 | 示例 |
|------|----------|------|
| 常量 | 全大写，下划线分隔 | `BOHR_RADIUS_ANGSTROM` |
| 类 | 大驼峰（PascalCase） | `GKK`, `GKKBufferIterator` |
| 结构体 | 大驼峰（PascalCase） | `GKKMetadata`, `KPointGVecs` |
| 成员变量 | 后缀下划线 | `meta_`, `fp_`, `current_kpt_` |
| 函数 | 小驼峰（camelCase） | `loadKPoint()`, `readMetadata()` |
| 局部变量 | 小驼峰 | `record_len`, `ikpt` |

## 示例

```cpp
export class GKK {
public:
    explicit GKK(const std::string& filename);
    const KPointGVecs& loadKPoint(int ikpt);

private:
    std::FILE* fp_;
    GKKMetadata meta_;
    int current_kpt_ = -1;
};
```
