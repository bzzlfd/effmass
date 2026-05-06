# 当前进度

## 已完成

### NCPP UPF 解析器 — beta / chi 数据截断

- `src/pseudo/ncpp-upf.cppm`：`beta` 按 `kbeta` 截断，`chi` 按尾部零扫描截断，并新增 `kchi[]` 字段保存有效长度。
- `test/test_ncpp_upf.cpp`：更新尺寸校验，验证截断后长度与 `kbeta`/`kchi` 一致，波函数归一化校验通过。
- `docs/design/pseudo_module.md`：补充“数据截断策略”设计说明。
