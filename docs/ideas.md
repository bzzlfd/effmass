#  ideas / 潜在改进方向

##  [idea-001] 使用 `std::ifstream` 替代 `std::FILE*` 进行二进制 IO

### 背景

当前 `io.cppm` 中的 `GKK` 和 `WG` 类使用 C 风格的 `std::FILE*`（`std::fopen`、`std::fread`、`std::fseek`）读取 Fortran unformatted binary 文件。

### 可行性分析

`std::ifstream` 原生支持二进制读取所需的所有操作：

- 二进制读取：`read(char_type* s, std::streamsize n)`
- 随机定位：`seekg(pos_type pos)`、`tellg()`
- 偏移定位：`seekg(off_type off, std::ios_base::seekdir dir)`（对应 `SEEK_SET`/`SEEK_CUR`/`SEEK_END`）

Fortran record 的读取逻辑完全可以平移。

### 优势

| 方面 | `std::ifstream` |
|------|-----------------|
| **RAII** | 析构时自动关闭文件，无需手动 `fclose` |
| **异常安全** | 可设置 `exceptions()` mask，让流错误自动抛异常 |
| **类型一致性** | 更符合现代 C++ 风格，与 `std::span`/`std::vector` 配合更自然 |
| **无 C 宏依赖** | 不需要 `SEEK_SET`/`SEEK_CUR` 宏，避免 `import std` 下的兼容性问题 |
| **扩展性** | 未来若要从内存流读取，接口统一 |

### 劣势

| 方面 | `std::ifstream` |
|------|-----------------|
| **性能** | 通常比 `std::FILE*` 慢。`std::istream` 有虚函数开销和额外的 locale/sentry 开销。对于大体积 Fortran 二进制文件，连续大量 `read` 时差距可能明显 |
| **迁移成本** | 需要同时重构已有的 `GKK` 类和 `WG` 类，改动面较大 |
| **底层控制** | `std::fread` 直接调用 C 库 read，行为更可预期 |

### 结论

暂时保留 `std::FILE*` 以最小化改动并保持性能。待后续若需要更强的可移植性（如从内存流读取）或性能差距可接受时，再统一迁移到 `std::ifstream`。

---

