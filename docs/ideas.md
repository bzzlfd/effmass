# NCPP Interface Design

## Context

`ncpp.cppm` (`pseudo.ncpp`) 是赝势的**算符层**——统一的、格式无关的计算接口。它的数据来自 `ncpp.upf.cppm`（IO 层），前者应成为外部消费者的唯一入口。

当前 `NCPP` 暴露了大量 UPF 特化类型的 pass-through（`header() -> NCPPUPFHeader&` 等），这耦合了所有调用方。上一轮已把 `meshType()` 和 `nonlocalByL()` 从 IO 层迁入算符层。现在需要设计完整的算符层接口。

## 设计原则

1. **算符层说"物理语言"**：element, angular momentum, radial mesh, projectors — 而非 "UPF header" 或 "PP_NONLOCAL"。
2. **每个 accessor 返回标准 C++ 类型**：`string_view`, `span<const double>`, `int`, `double`, `bool`。不返回格式相关的复合类型。
3. **不平行复制类型体系**：不创建与 UPF 结构体一一对应的 `NCPPHeader` / `NCPPMesh` 等包装类型——它们最终不过是改名的字节复制，没有抽象收益。

## 接口设计

### 构造

```cpp
explicit NCPP(const NCPPUPF& upf);  // 保持现状
```

未来多格式时 `upf_` 变为 `std::variant<NCPPUPF, PSP8Data, ...>`，所有 accessor 内部 dispatch。不设 factory——以后再说。

### 格式无关 accessor（新增，替换现有的 UPF 类型 pass-through）

```
auto element() const -> std::string_view;
auto zValence() const -> double;
auto lMax() const -> int;
auto lLocal() const -> int;
auto coreCorrection() const -> bool;
auto meshR() const -> std::span<const double>;
auto meshRab() const -> std::span<const double>;
auto meshSize() const -> int;
auto localPotential() const -> std::span<const double>;
auto rhoAtom() const -> std::span<const double>;
```

注意：`localPotential()` 和 `rhoAtom()` 已经返回 `span<const double>`（格式无关），不用改。

### Nonlocal 和 wavefunction accessor（新增）

```
auto numBetaProjectors() const -> int;
auto betaProjector(int i) const -> std::span<const double>;
auto betaAngularMomentum(int i) const -> int;
auto betaCutoffIndex(int i) const -> int;
auto dionMatrix() const -> std::span<const double>;   // flat, row-major, nproj x nproj

auto numWavefunctions() const -> int;
auto wavefunction(int i) const -> std::span<const double>;
auto wavefunctionAngularMomentum(int i) const -> int;
auto wavefunctionOccupation(int i) const -> double;
auto wavefunctionLabel(int i) const -> std::string_view;
```

### 移除现有的 UPF 类型 pass-through

上一轮添加的 `header()`, `mesh()`, `nonlocal()`, `wavefunctions()` **没有任何外部调用者**（test 直接通过 `upf.header()` 访问，不通过 `ncpp`），直接删除。`upfData()` 保留作为 escape hatch。

### 保留

| 成员 | 说明 |
|------|------|
| `meshType() -> MeshType` | 已迁入，不动 |
| `nonlocalByL(int) -> NonlocalByL` | 已迁入，不动 |
| `upfData() -> const NCPPUPF&` | escape hatch，保留 |

### 导出类型

`export` 块：
```
enum class MeshType { Uniform, Exponential, Unknown };
struct NonlocalByL { vector<vector<double>> beta; DenseMatrix dion; };
class NCPP;
```

`DenseMatrix` 保持模块内部。它和 `ncpp.upf.cppm` 中的 `Matrix` 是同等逻辑的重复类型——等将来 `pseudo -> math` 的依赖建立后，统一迁入 `math.linalg`。本轮不动。

### 本轮不做的

- **Fourier 计算** (`localPotentialFourier(q)`, `betaProjectorFourier(l, i, q)`)：这些需要引入 `pseudo -> math.fourier_bessel` 的依赖，属于下一轮。当前集中精力把数据访问接口做干净。
- **MeshType → RadialMeshType 桥接**：等加 Fourier 方法时一起做。

## 需修改的文件

| 文件 | 改动 |
|------|------|
| `src/pseudo/ncpp.cppm` | 新增格式无关 accessor；移除 UPF 类型 pass-through；调整 export 块 |
| `test/pseudo/io/test_ncpp_upf.cpp` | 对新的 `ncpp.*()` accessor 增加测试（与现有 `upf.*()` 测试并存） |
| `test/pseudo/test_ncpp.cpp` | 扩展为有实际断言的测试：验证所有新 accessor 返回正确值 |

## 验证

```bash
cd build && cmake --build . && ctest --output-on-failure
```
