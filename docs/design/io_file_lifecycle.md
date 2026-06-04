# IO 模块文件句柄生命周期策略

> 各 IO 类采用何种文件句柄策略，取决于其读取模式（一次性 vs. 延迟加载）。**Rule of Zero 是默认目标**，只在确实需要随机访问时才持有 `FILE*`。

---

## 一、一览

| 类 | 模块 | 读取模式 | 持有 FILE* | 特殊成员 | 原因 |
|----|------|---------|-----------|---------|------|
| `GKK` | `io.GKK` | 延迟加载 | ✅ 构造–析构 | 自定义析构/移动，禁止拷贝 | 需要 `fseek` 按 k 点随机读取 |
| `WG` | `io.WG` | 延迟加载 | ✅ 构造–析构 | 同上 | 同上 |
| `EIGEN` | `io.EIGEN` | 一次性 | ❌ | Rule of Zero | 数据读完即可关闭 |
| `RHO` / `VR` | `io.RHO` / `io.VR` | 一次性 | ❌ | Rule of Zero | 数据读完即可关闭 |
| `ATOM` | `io.ATOM` | 一次性 | ❌ | Rule of Zero | 数据读完即可关闭 |
| `OCC` | `io.OCC` | 一次性 | ❌ | Rule of Zero | 数据读完即可关闭 |
| `NCPPUPF` | `pseudo.io.ncpp_upf` | 一次性 | ❌ | Rule of Zero | pugixml `load_file` 内部关闭 |

**规则很简单**：如果你能在构造结束前把全部数据装进内存，就别留着文件句柄。

---

## 二、为何追求 Rule of Zero

C++ 的 Rule of Zero 主张：如果类的成员本身都是自管理类型（`std::vector`、`std::string` 等），编译器生成的默认特殊成员函数就是正确的——不需要写析构、拷贝、移动。这带来：

| 优势 | 说明 |
|------|------|
| 更少的代码 | 不需要手写析构/移动/拷贝，减少样板 |
| 更少的 bug | 不会漏掉 `fp_ = nullptr`，不会忘记 `fclose` |
| 更好的异常安全 | 构造中异常时，已构造的成员会自动析构 |
| 自然的值语义 | 拷贝/移动开箱即用，不需要额外维护 |

唯一的代价是：你必须确保**所有需要的数据在构造结束前都已读入内存**。

---

## 三、为什么 GKK 和 WG 必须持有 FILE*

GKK 和 WG 文件较大（单个 k 点可能数百 KB，nkpt 可达数百），且大多数场景只需要访问其中几个 k 点。因此它们采用**延迟加载（lazy loading）**：

```cpp
GKK gkk("OUT.GKK");           // 只读元数据 + 预计算文件偏移
// ...
auto& kpt5 = gkk.loadKPoint(5);  // fseek → fread → 返回数据
auto& kpt0 = gkk.loadKPoint(0);  // 再次 seek → 覆盖 buffer
```

如果构造完就关闭文件，`loadKPoint` 就无法工作。文件句柄是 GKK/WG 对象状态的必要组成部分，因此它们：

- 有自定义析构函数（`fclose`）
- 删除拷贝（两个对象不能共享同一个 `FILE*`）
- 自定义移动语义（转移 `FILE*` 所有权，源对象置空）

这只影响 2 个类，其余 6 个 IO 类都是 Rule of Zero。

---

## 四、局部 FILE* 的异常安全模式

所有"读完即关"的类都使用以下模式：

```cpp
auto* fp = std::fopen(filename.c_str(), "rb");
if (!fp) throw /* ... */;

try {
    readMetadata(fp);   // 可能抛异常
    readData(fp);       // 可能抛异常
} catch (...) {
    std::fclose(fp);    // 确保异常路径不泄漏
    throw;
}
std::fclose(fp);        // 正常路径关闭
```

注意这里不使用 `std::unique_ptr<FILE, decltype(&fclose)>`——手动的 try/catch 模式更直接，且不引入额外的抽象开销。

---

## 五、演进历史

1. **早期**：GKK/WG/RHO/EIGEN 全部持有 `FILE*`，认为"IO 类就应该持有文件句柄"。ATOM 也沿用此模式但实际不需要。

2. **中期**：`NCPPUPF`（UPF 解析）采用了"读完即关"的 Rule of Zero 模式。这在当时被称为"违背祖训"，因为与其它 IO 类截然不同。OCC 模块也采用了相同的局部 FILE* 模式。

3. **当前**：将 ATOM、EIGEN、RHO 也改为局部 `FILE*` + Rule of Zero。现在只有**真正需要延迟加载的 GKK 和 WG** 才持有 `FILE*`。曾经特立独行的 `NCPPUPF` 和 `OCC` 反而成了主流模式，而 ATOM、EIGEN、RHO 也加入了这一行列。

---

## 六、如果将来 GKK/WG 也想变成 Rule of Zero？

如果哪天项目决定一次性加载全部 k 点，GKK/WG 也可以改为局部 `FILE*` + Rule of Zero。但这会改变 API 语义（`loadKPoint` 不再需要文件句柄），且增加内存占用（所有 k 点常驻内存）。当前延迟加载的设计是在**内存**和**IO**之间做的合理取舍。
