# ncpp 接口重构

以 docs/reference/unified_pseudopotential_format.md 为 UPF 文件参考。
它提供了完整ncpp数据列表

## ncpp.cppm 和 ncpp.upf.cppm 角色
1. ncpp.upf.cppm 是纯粹的 upf 文件的数据接口抽象
2. ncpp.cppm 是纯粹的 ncpp 的数据接口抽象

  1. 把 *.upf.cppm 精简成 upf.cppm
  2. upf 本身可以表示三种pseudopotential，upf.cppm 也应该表示三种PsP。
  3. ncpp.cppm 可以从 upf 文件构造出 ncpp，也可以从 （比如psp8 文件）构造出ncpp，然后作为数据接口提供给别人
  

## ncpp.cppm 提供接口原则
1. 不使用 ncpp.upf.cppm 中变量命名方式，而是使用更贴近其含义的命名
  1. D 矩阵（Dion）在ncpp中改名为B矩阵，这是Vanderbilt论文中的命名方式
2. NCPP 类除了有对应的数据接口，还应提供若干便利的成员函数

## ncpp.cppm 数据存储
成员变量
- Info
- Header
- Mesh
  - MeshType
- NLCC (optional)
- Local
- Nonlocal
- Semilocal (optional, only for norm-conserving)
- PsWFc (optional)
- RhoAtom

以上是public的。
其中，Semilocal 数据是可选的。我们后面会有程序判断，如果ncpp是semilocal的，那么报错，说我们现在只支持 Nonlocal ，为这个程序判断做好准备。
mesh 中的数据与 beta 中的数据按位置对应。但是beta中在某处往后全都是0, beta存储数据长度小于 mesh
其内部数据组织方式参考 ncpp.upf.cppm，但是内部变量名要选择更突出物理含义的



## ncpp.cppm 成员函数（不完全）
- public
  - NonlocalByL
- private
  - inferMeshType (在初始化时调用)

