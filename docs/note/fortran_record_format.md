# Fortran Record 格式

本项目需要读取 Fortran unformatted binary 文件（如 `OUT.GKK`、`OUT.WG`、`OUT.EIGEN`）。这些文件采用特定的 **Fortran Record** 格式存储数据。

## Record 结构

每个 record 的结构为：

```
[4-byte length header][data][4-byte length footer]
```

- **头部（4 bytes）**：一个 32 位整数，表示后续 `data` 的字节长度。
- **数据（data）**：实际存储的数据，长度由头部指定。
- **尾部（4 bytes）**：另一个 32 位整数，通常与头部数值相同，用于一致性校验。

> **注意**：长度标记的具体字节数（4 字节或 8 字节）取决于 Fortran 编译器选项，但在本项目使用的数据中通常为 **4 字节 `int`**。

## 读取流程

```cpp
// 1. 读取 4 字节头部，获取 record 长度
int record_len;
fread(&record_len, sizeof(int), 1, fp);

// 2. 根据长度读取数据
std::vector<char> buffer(record_len);
fread(buffer.data(), 1, record_len, fp);

// 3. 读取 4 字节尾部并校验
int record_len_end;
fread(&record_len_end, sizeof(int), 1, fp);
if (record_len != record_len_end) {
    // 长度不匹配，抛出异常
    throw std::runtime_error("Fortran record length mismatch");
}
```
