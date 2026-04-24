# effmass 项目文档

## 项目概述

`effmass` 是一个第一性原理电子结构计算工具，用于读取波函数文件（平面波基组下的平面波波矢`OUT.GKK`, 平面波系数`OUT.WG` 等），进行计算得到有效质量。


## 技术栈

- **主要语言**: C++
- **编译工具**：clang++, cmake, ninja
- **文档**: Markdown
- **文件格式**: Fortran unformatted binary（OUT.GKK, OUT.WG, OUT.EIGEN）、 text（OUT.KPT）

## 构建步骤

```bash
cmake -B build -G Ninja
cmake --build build
```

C++23 modules + `import std` 实验性支持，依赖 libc++。构建失败时删除 `build/` 目录重新配置。

## 项目结构

```
.
├── AGENTS.md           # 本文件 - AI 编码代理指南
├── src                 # 源代码地址
├── docs                # 开发路线、目前正在开发的进程、编码规范、设计细节和关键决定、公式笔记、参考资料
├── test                # 测试
└── .tree               # 多开发目标时，使用 git worktree 同时工作，在此目录下创建新的 worktree
```

### src 目录


### docs 目录
- `index.md`：对 `docs/` 的索引。用简洁的词语/话概括每个文件中的内容，根据这个信息，LLM足够确定去哪个目录查看具体稳定。
- `roadmap.md`：项目整体规划，*当前认可的计划（但不是绝对真理）*
  - `- [x]` 表示已完成
  - `- [/]` 表示已开始，正在进行，未完成
  - `- [ ]` 表示未开始
- `prograss.md`：当前进度，目前着手做的工作。*开发者编写*
- `ideas.md`：我或ai提出的潜在方向，改进建议。*ai编辑*
- `note/`: 存放先于代码开发存在的知识：如计算中用到的公式笔记、文件格式解析
- `design/`：存放代码开发阶段的约定、选择和权衡过程等
- `reference/`: 存放原始参考资料

每次完成一个进度之后，思考是否需要更新 `docs` 目录。
每添加一个一个文档，在 `index.md` 中插入简介的概括。




## Git commit messages

Conventional Commits format:

```
<type>(<scope>): <subject>
```

Allowed types:
`feat`, `fix`, `test`, `refactor`, `chore`, `style`, `docs`, `perf`, `build`, `ci`, `revert`.



## 开发注意事项

### 错误处理

- 文件读取失败时抛出 `std::runtime_error`
- 索引越界时抛出 `std::out_of_range`
- Record 长度不匹配时抛出异常

### 内存管理

- 使用 RAII 管理文件句柄
- 预分配缓冲区避免频繁内存分配
- 实现移动语义支持（move constructor/assignment）
- 禁止拷贝（delete copy constructor/assignment）

### 性能考虑

- 延迟加载：构造函数不加载全部数据
- 缓存策略：缓存最近访问的 k 点数据
- 偏移预计算：构造函数中预计算每个 k 点的文件偏移
- 缓冲区复用：避免重复分配内存

