像我们之前 写读取 GKK 做的那样，编写抽象 WG 的类
1. 读取 docs/reference/plot_wg.f90 源码，提取读取 OUT.WG 文件的代码
2. 将该文件格式分析写进 file_formats.md
3. 在 src/io.cppm 中撰写 WG 相关类。代码风格，框架设计可以参考 
    1. docs/design 的内容，
    2. io.cppm 中关于 GKK 的实现
4. 在 docs/io_cppm_design.md 中撰写设计思路，可以参考其中 GKK 的内容，重复的部分简写。
5. docs 当中的信息发生修改，同步 index.md 内容


1. 自旋轨道耦合
2. 检查两个对象中，相同数据的一致性
3. [/] 多个缓存数据。WG 类已实现多 band 缓存，容量由构造函数参数控制（默认 4），采用 FIFO 环形替换策略。