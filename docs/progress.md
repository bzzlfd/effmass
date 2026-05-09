# 当前进度

## 已完成

### NCPP UPF 解析器 — beta / chi 数据截断

- `src/pseudo/ncpp-upf.cppm`：`beta` 按 `kbeta` 截断，`chi` 按尾部零扫描截断，并新增 `kchi[]` 字段保存有效长度。
- `test/test_ncpp_upf.cpp`：更新尺寸校验，验证截断后长度与 `kbeta`/`kchi` 一致，波函数归一化校验通过。
- `docs/design/pseudo_module.md`：补充"数据截断策略"设计说明。

### NCPP UPF 解析器 — mesh 类型推断与按 l 筛选非局域项

- `src/pseudo/ncpp-upf.cppm`：
  - 新增 `MeshType` 枚举（`Uniform`, `Exponential`, `Unknown`）。
  - 新增 `NCPPUPFNonlocalByL` 结构体。
  - `NCPPUPF` 新增 `meshType()` 函数：通过分析 `r` 和 `rab` 的数值关系推断网格类型（均匀或指数/对数）。
  - `NCPPUPF` 新增 `nonlocalByL(int l)` 函数：返回指定角动量 `l` 的 beta 子列表和对应的 `D_ij` 子矩阵。
- `test/test_ncpp_upf.cpp`：新增 `meshType()` 和 `nonlocalByL()` 的测试用例，包括 beta 数量、dion 尺寸、数值一致性校验。
- `docs/design/pseudo_module.md`：更新公共接口与数据结构说明，补充 mesh 类型推断和 `nonlocalByL` 的设计文档。

## 进行中

1. GKK
    1. current_ikpt 这个语义上像是可以得知 k 向量，可以得到吗？
    推断 k 向量
    2. 对于currentData的返回值，内部保存一个变量，用于控制currentData以那些视角返回
        1. 源数据 gg(把这个重新起一个不会误会的名字),gx,gy,gz （叫g... 语义上也不好，让人以为是 G，其实是 G+k）
        2. 球坐标视角（）
        3. 整数视角 （）

添加命名规范。
k: (0.5, 0.5]
G: (recip. latt.) * (i,j,k)
K: G - k

ik index of k

2. 特殊函数
3. ~~给 ncpp-upf.cppm 中的 NCPPUPF 添加更多函数~~
    1. ~~推断 mesh 中的网格是均匀的，还是指数的，还是未知的~~
    2. ~~获得角动量是 l 的 beta 子列表，dion子矩阵~~

- span 还是 vector，需要思考