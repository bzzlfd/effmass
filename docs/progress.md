```
主要参考 quantum espresso 中读取 UPF 文件的过程。

pseudo_iosys -> readpp -> read_ps_new -> read_upf_new


1. 只需要考虑 NCPP 的实现
2. 不必考虑 PAW（tpawp，gipaw），USPP （tvanp）
3. 不必考虑 metagga （with_metagga_info）代码
4. 不必考虑 SOC（has_so, spin-orbital coupling, 全相对论）情况
```

## 已完成

- [x] 总结 QE `read_upf_new.f90` 和 `xmltool.f90` 的抽象，写入 `ideas.md`
- [x] 嵌入 pugixml v1.15 作为 XML 解析依赖
- [x] 实现 `pseudo.ncpp_upf` 模块：完整读取 UPF v.2 格式（Header / Mesh / Local / Nonlocal / PSWFC / RhoAtom）
- [x] 验证 UPF 中 `chi` 与 `beta` 的存储约定：`χ(r) = r·R(r)`，归一化 `∫χ² dr = 1`
- [x] 编写 `test_ncpp_upf` 测试，使用 Ge 赝势数据验证读取正确性与波函数归一化

- NCPP 中遇到 pp_semilocal 时，报错没有实现


- UPF 具体变量含义
    - [Unified Pseudopotential Format](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
    - [Soft self-consistent pseudopotentials in a generalized eigenvalue formalism](https://doi.org/10.1103/PhysRevB.41.7892) 
- UPF 文件中，很多数据tag有columns属性，用来指示一行有多少数据。但是代码里没有明显的使用这个属性，解释原因
- UPF 束缚态波函数的归一化条件
- NCPPUPF 类内部的数据结构与UPF文件结构的对应。public接口设计原则，用于便捷取用其中特定数据


- design 中只需要添加一个 pseudo_module ，里面的内容有
    1. 模块结构 (现在很好)
        1. NCPP 
            1. 参考资料
            2. 数据结构与 UPF 文件 tag 的对应关系 （现在写的很好）
            3. 构造函数一次性全加载
            2. 程序中不读取 `columns` 的原因
            3. 束缚态波函数 $\chi(r) = \sqrt{4\pi} r \cdot R(r)$ ,$\int_0^{\infty} |\chi(r)|^2 \, dr = 1$


## 当前工作

- NCPP UPF 读取功能已完成并测试通过
- USPP / PAW 模块仅创建骨架，待后续实现

## 未来工作
存储 beta 之类的有 cutoff 的数据时， 只存储 cutoff 段的