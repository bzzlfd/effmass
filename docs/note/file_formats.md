# 文件格式解析

关于 `io.cppm` 中 C++ 接口的设计细节，请参阅 [`design/io_cppm_design.md`](../design/io_cppm_design.md)。

## OUT.GKK

`OUT.GKK` 文件存储每个 k 点对应的平面波波矢信息。该文件为 Fortran unformatted binary 格式，采用 [Fortran Record 格式](fortran_record_format.md)。

### 文件结构（按读取顺序）

以下是从 `docs/reference/plot_wg.f90` 中提取的 `OUT.GKK` 读取逻辑：

```fortran
open(11, file="OUT.GKK", form="unformatted")
rewind(11)

! Record 1: 元数据头
read(11) n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda

! Record 2: 截断能
read(11) Ecut

! Record 3: 晶格矢量 AL(3,3)，文件中单位为 Å
read(11) AL
AL = AL / A_AU_1   ! 转换为原子单位（Bohr）

! Record 4: 每个 k 点、每个节点上的 G 向量数量
read(11) nnodes, ((ngtotnod_9(inode, kpt), inode=1, nnodes), kpt=1, nkpt)

! 自旋轨道耦合：有效 G 向量数减半
if (is_SO .eq. 1) then
    ngtotnod_9 = ngtotnod_9 / 2
endif

! 后续：按 k 点顺序存储 G 向量数据
do kpt = 1, nkpt
    num = 0
    do inode = 1, nnodes
        read(11) gkk_n_tmp      ! |G+k|^2 / 2
        read(11) gkk_n_xtmp     ! K_x = G_x - k_x
        read(11) gkk_n_ytmp     ! K_y = G_y - k_y
        read(11) gkk_n_ztmp     ! K_z = G_z - k_z
        do ig = 1, ngtotnod_9(inode, kpt)
            gkk(num+ig)   = gkk_n_tmp(ig)
            gkk_x(num+ig) = gkk_n_xtmp(ig) + akx(kpt)
            gkk_y(num+ig) = gkk_n_ytmp(ig) + aky(kpt)
            gkk_z(num+ig) = gkk_n_ztmp(ig) + akz(kpt)
        enddo
        num = num + ngtotnod_9(inode, kpt)
    enddo
    ng_tot = num
enddo
close(11)
```

### 数据存储方式总结

| Record | 内容 | 数据类型 | 说明 |
|--------|------|----------|------|
| 1 | `n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda` | 8 × `int` | FFT 网格、单节点最大 G 数、节点数、k 点数、自旋轨道耦合标志、自旋极化标志 |
| 2 | `Ecut` | 1 × `double` | 平面波截断能 |
| 3 | `AL(3,3)` | 9 × `double` | 晶格矢量，文件中为 **Å**，读取后应转为 **Bohr** |
| 4 | `nnodes`, `ngtotnod(inode, kpt)` | `int` + `nnodes × nkpt` × `int` | 每个节点在每个 k 点上的 G 向量数量 |
| 5+ | k 点循环数据 | 每个节点 4 个 record | `gkk`, `gkk_x`, `gkk_y`, `gkk_z` |

#### k 点数据布局

对于每个 k 点 `kpt = 1, ..., nkpt`：
- 遍历每个计算节点 `inode = 1, ..., nnodes`
- 每个节点包含 **4 个独立的 Fortran record**：
  1. `gkk(ngtotnod)` — `double` 数组，存储 `|G+k|^2 / 2`
  2. `gkk_x(ngtotnod)` — `double` 数组，存储 $K_x = G_x - k_x$
  3. `gkk_y(ngtotnod)` — `double` 数组，存储 $K_y = G_y - k_y$
  4. `gkk_z(ngtotnod)` — `double` 数组，存储 $K_z = G_z - k_z$

> **注意**：尽管分量存储的是 `G - k`，但动能项 `gkk` 仍按 `|G+k|^2 / 2` 计算。

### 自旋轨道耦合处理

当 `is_SO == 1` 时，文件中每个 G 向量实际上包含自旋向上和向下两个分量，因此**有效 G 向量数量需要除以 2**。这会影响以下量：

- **`mg_nx`**：单节点最大 G 向量数需要除以 2
- **`ngtotnod(inode, kpt)`**：每个节点在每个 k 点上的 G 向量数需要除以 2

**处理原则**：
1. 读取元数据头后立即将 `mg_nx /= 2`
2. 读取 `ngtotnod` 数组后，对每个元素执行 `ngtotnod /= 2`
3. 后续按调整后的 `ngtotnod` 读取每个 record 中的 `gkk*` 数据

这样做之后，后续读取 `OUT.WG` 时也应使用同样的 `ngtotnod`（已除以 2），以保证波函数系数与 G 向量一一对应。

## OUT.WG

`OUT.WG` 文件存储平面波基组下的波函数展开系数。该文件同样为 Fortran unformatted binary 格式。

### 文件结构（按读取顺序）

以下是从 `docs/reference/plot_wg.f90` 中提取的 `OUT.WG` 读取逻辑：

```fortran
open(11, file=filename, form="unformatted")
rewind(11)

! Record 1: 元数据头（与 OUT.GKK 类似，但多了一个 mx）
read(11) n1_t, n2_t, n3_t, mx, mg_nx, nnodes, nkpt, is_SO, islda

! Record 2: 截断能
read(11) Ecut_t

! Record 3: 晶格矢量 AL(3,3)，文件中单位为 Å
read(11) AL_t
AL_t = AL_t / A_AU_1   ! 转换为原子单位（Bohr）

! Record 4: 每个 k 点、每个节点上的 G 向量数量
read(11) nnodes, ((ngtotnod_9_t(inode, kpt), inode=1, nnodes), kpt=1, nkpt)

if (is_SO .eq. 1) then
    ngtotnod_9_t = ngtotnod_9_t / 2
endif

! 后续：按 k 点 → 节点 → 能带 顺序存储波函数系数
do kpt = 1, nkpt
    do inode = 1, nnodes
        do im = 1, mx
            read(11) ug_n_tmp
            do ig = 1, ngtotnod_9_t(inode, kpt)
                ug_one(num+ig) = ug_n_tmp(ig)
            enddo
            if (is_SO .eq. 1) then
                do ig = 1, ngtotnod_9_t(inode, kpt)
                    ug_two(num+ig) = ug_n_tmp(ngtotnod_9_t(inode, kpt) + ig)
                enddo
            endif
            num = num + ngtotnod_9_t(inode, kpt)
        enddo
    enddo
enddo
close(11)
```

### 数据存储方式总结

| Record | 内容 | 数据类型 | 说明 |
|--------|------|----------|------|
| 1 | `n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda` | 8 × `int` | FFT 网格、**能带数**、单节点最大 G 数、节点数、k 点数、自旋轨道耦合标志、自旋极化标志 |
| 2 | `Ecut` | 1 × `double` | 平面波截断能 |
| 3 | `AL(3,3)` | 9 × `double` | 晶格矢量，文件中为 **Å**，读取后应转为 **Bohr** |
| 4 | `nnodes`, `ngtotnod(inode, kpt)` | `int` + `nnodes × nkpt` × `int` | 每个节点在每个 k 点上的 G 向量数量 |
| 5+ | k 点 → 节点 → 能带 循环数据 | 每个 `(inode, im)` 一个 record | 波函数系数 `ug` |

#### 波函数数据布局

对于每个 k 点 `kpt = 1, ..., nkpt`：
- 遍历每个计算节点 `inode = 1, ..., nnodes`
- 每个节点内遍历每个能带 `im = 1, ..., mx`
- 每个 `(inode, im)` 对应 **一个独立的 Fortran record**：
  - `ug(ngtotnod)` — `complex*16` 数组，存储该节点上该能带的波函数系数

### 自旋轨道耦合处理

当 `is_SO == 1` 时：

1. **`ngtotnod`** 需要除以 2（同 OUT.GKK）
2. 每个波函数 record 包含 **两个连续存储的 `complex*16` 数组**，长度各为 `ngtotnod`：
   - 前半部分：`ug_up` — 自旋向上分量
   - 后半部分：`ug_down` — 自旋向下分量

> 因此读取一个 band 数据时，若 `is_SO == 1`，record 的总长度为 `2 × ngtotnod × sizeof(complex*16)`；否则为 `ngtotnod × sizeof(complex*16)`。

### 与 OUT.GKK 的关系

`OUT.WG` 的元数据头部（`n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda, Ecut, AL, ngtotnod`）**必须与 `OUT.GKK` 保持一致**（仅额外多出 `mx` 字段）。实际使用时应验证两者元数据的一致性。

## OUT.EIGEN

`OUT.EIGEN` 文件存储能量本征值和 k 点信息。该文件为 Fortran unformatted binary 格式，采用 [Fortran Record 格式](fortran_record_format.md)。

### 文件结构（按读取顺序）

```fortran
open(23, file="OUT.EIGEN", form="unformatted")
rewind(23)

! Record 1: 元数据头（7 ints，旧格式可能为 6 ints）
read(23) islda, nkpt, mx, nref_tot_8, natom, nnodes, is_SO

allocate(weighkpt_2(nkpt))
allocate(akx_2(nkpt))
allocate(aky_2(nkpt))
allocate(akz_2(nkpt))
allocate(E_st(mx, nkpt, islda))

! 后续：按自旋 → k 点顺序存储
do iislda = 1, islda
    do kpt = 1, nkpt
        read(23) iislda_tmp, kpt_tmp, weighkpt_2(kpt), &
            &       akx_2(kpt), aky_2(kpt), akz_2(kpt)
        read(23) (E_st(i, kpt, iislda), i = 1, mx)
    end do
end do
close(23)
```

### 数据存储方式总结

| Record | 内容 | 数据类型 | 说明 |
|--------|------|----------|------|
| 1 | `islda, nkpt, mx, nref_tot_8, natom, nnodes [, is_SO]` | 6 或 7 × `int` | 元数据。旧格式无 `is_SO`，新格式含 `is_SO`。通过 record 长度区分（24 vs 28 字节） |
| 2~ | k 点循环数据 | 每个 `(iislda, kpt)` 2 个 record | 见下方 |

#### k 点数据布局

对于每个自旋 `iislda = 1, ..., islda`，每个 k 点 `kpt = 1, ..., nkpt`：

1. **k 点元数据 record**：`iislda_tmp` (int), `kpt_tmp` (int), `weighkpt_2(kpt)` (double), `akx_2(kpt)` (double), `aky_2(kpt)` (double), `akz_2(kpt)` (double)
   - 总大小：`2 × sizeof(int) + 4 × sizeof(double) = 40` 字节
   - `iislda_tmp` 应与外层循环 `iislda` 一致，`kpt_tmp` 应与 `kpt` 一致；不一致时报错
   - weighkpt_2 和 akx/aky/akz 按 kpt 索引（非自旋），后一个自旋覆盖前一个

2. **本征值 record**：`mx` 个 double，**单位为 eV**
   - 总大小：`mx × sizeof(double)` 字节
   - 按列优先（Fortran column-major）存储：能带索引变化最快

### 兼容性说明

- **is_SO 处理**：读取时先尝试读 7 int（新格式），通过 record 长度（28 字节）判断；若 record 长度为 24 字节则为旧格式，设 `is_SO = 0`
- 与 GKK/WG 不同，EIGEN 文件中 `is_SO` 仅作标记，不改变数据布局
- `weighkpt_2`、`akx_2`、`aky_2`、`akz_2` 在文件中对每个 `(iislda, kpt)` 对都有存储，但最终只保留最后读入的值（按 kpt 索引，非自旋）
- `nref_tot_8`：该字段含义尚不明确。Fortran 注释显示不同文件间（如 `eigen_all.store` 与 `bpsiifile`）该值可能不同，读取时原样保留不做解释

## OUT.OCC

`OUT.OCC` 文件以**文本格式**存储每个 k 点每个能带的能量本征值（eV）和占据数。该文件仅在**非局域势**计算中输出，局域势计算不产生此文件。

### 文件结构（按读取顺序）

```
KPOINTS      1:    0.0000    0.0000    0.0000
 NO.   ENERGY(eV) OCCUPATION
      1    -94.5592   0.04000
      2    -94.5560   0.04000
    ...
KPOINTS      2:    0.0000    0.0000    0.3314
 NO.   ENERGY(eV) OCCUPATION
      1    ...
```

### 布局说明

#### 单自旋格式（islda=1）

每个 k 点包含 2 行头信息和 `nband` 行数据：

| 行 | 内容 | 说明 |
|----|------|------|
| 1 | `KPOINTS N: x y z` | k 点序号（1-based）、分数坐标 |
| 2 | ` NO.   ENERGY(eV) OCCUPATION` | 列标题 |
| 3+ | `iband energy occupation` | 能带序号（1-based）、本征值（eV）、占据数（已包含 k 点权重） |

#### 双自旋格式（islda=2）

文件以 `==========  SPIN 1  ==========` 开头，分两段存储，每段结构与单自旋格式一致：

```
==========  SPIN 1  ==========
KPOINTS      1:    0.0000    0.0000    0.0000
 NO.   ENERGY(eV) OCCUPATION
      1    -94.5597   0.02000
      2    -94.5566   0.02000
      ...
==========  SPIN 2  ==========
KPOINTS      1:    0.0000    0.0000    0.0000
 NO.   ENERGY(eV) OCCUPATION
      1    -93.5597   0.02000
      ...
```

两自旋共享相同的 k 点网格和能带数。k 点坐标在第二段中与第一段一致，读取时会自动校验一致性。

### 与 OUT.EIGEN 的关系

`OUT.OCC` 与 `OUT.EIGEN` 共享相同的 k 点网格和能带数。两个文件中的能量本征值一致（OUT.OCC 四舍五入到 4 位小数），占据数已包含 k 点权重。

### C++ 读取说明

`OCC` 类（`io.OCC` 模块）使用 `fopen("r")` + `fgets` 逐行解析（与 `ATOM` 类的文本解析策略一致）：

1. 检测首行：若以 `=` 开头则进入双自旋模式，按 `==========  SPIN N  ==========` 分隔两段数据
2. 匹配 `KPOINTS` 关键词 → 解析 k 点索引和分数坐标
3. 遇到数字开头的行 → 解析能带序号、能量（eV）、占据数
4. 其他行（列标题、空行、SPIN 分隔线等）→ 跳过
5. 双自旋模式下，第二自旋的 k 点坐标会与第一自旋自动校验一致性
6. 所有数据一次性读取并验证一致性（islda × nkpt × nband 数据完整性）

## OUT.KPT

`OUT.KPT` 文件以**文本格式**存储 k 点坐标。

> TODO: 待补充具体的列定义和解析方法。

## OUT.VR / OUT.RHO

`OUT.VR` 和 `OUT.RHO` 文件分别存储价电子波函数势和电荷密度在实空间 FFT 网格上的值。两者文件格式**完全相同**，因此 C++ 中 `VR` 是 `RHO` 的类型别名。该文件为 Fortran unformatted binary 格式。

### 文件结构（按读取顺序）

```fortran
subroutine input_vr()
implicit double precision (a-h,o-z)
integer ierr

open(11, file=filename, form="unformatted")
rewind(11)

! Record 1: 元数据头（5 ints，旧格式可能为 4 ints）
read(11, IOSTAT=ierr) n1, n2, n3, nnodes, nstate

if (ierr .ne. 0) then
    rewind(11)
    read(11, IOSTAT=ierr) n1, n2, n3, nnodes
    nstate = 1
endif

! Record 2: 晶格矢量 AL(3,3)，文件中单位为 Å
read(11) AL

nr = n1 * n2 * n3
nr_n = nr / nnodes
allocate(vr_tmp(nr_n))
allocate(vr(n1, n2, n3))

! 后续：按 state → node 顺序存储实空间数据
do ist = 1, nstate
    do iread = 1, nnodes
        read(11) vr_tmp
        do ii = 1, nr_n
            jj = ii + (iread - 1) * nr_n
            i = (jj - 1) / (n2 * n3) + 1
            j = ((jj - 1) - (i - 1) * n2 * n3) / n3 + 1
            k = jj - (i - 1) * n2 * n3 - (j - 1) * n3
            vr(i, j, k) = vr_tmp(ii)
        end do
    end do
end do
close(11)
```

### 映射说明

Fortran 代码中的 1-based 线性索引到三维索引的映射：

```
jj = (i - 1) * n2 * n3 + (j - 1) * n3 + k     ! 1-based
```

其中 **k 变化最快**（步长 1），**i 变化最慢**（步长 n2 × n3）。C++ 中使用的 0-based 公式：

```
idx = i * n2 * n3 + j * n3 + k                  ! 0-based
```

### 数据存储方式总结

| Record | 内容 | 数据类型 | 说明 |
|--------|------|----------|------|
| 1 | `n1, n2, n3, nnodes [, nstate]` | 4 或 5 × `int` | FFT 网格维度和节点数。旧格式无 `nstate`，设为 1。通过 record 长度区分（16 vs 20 字节） |
| 2 | `AL(3,3)` | 9 × `double` | 晶格矢量，文件中为 **Å**，读取后转为 **Bohr** |
| 3~ | `vr_tmp(nr_n)` | `nr_n` × `double` | 每个 `(state, node)` 对应一个 record，共 `nstate × nnodes` 个 record |

#### 数据读取与重组

文件按 `state` 外层循环、`node` 内层循环存储。每个 record 包含 `nr_n = n1 × n2 × n3 / nnodes` 个 double，对应一个节点上的网格切片。

线性索引到网格坐标的映射（与存储顺序一致，无需重排）：

```
data[state][i][j][k] = file_data[state][node][ii]
                     → idx = state × nr + node × nr_n + ii
```

其中 `ii` 为 record 内线性索引，`idx` 为扁平 `data_` 中的位置。

### nstate 说明

- **nstate 字段**：早期版本的文件格式字段，表示数据文件中存储的"状态"数。实际使用中，当存在第二个自旋时，数据位于独立的 `OUT.RHO_2` 文件中（而非在同一文件中增加 state）。
- 旧格式文件不含此字段（record 长度为 16 字节 = 4 ints），此时 `nstate = 1`。
- 新格式文件含此字段（record 长度为 20 字节 = 5 ints）。

## atom.config

`atom.config` 文件以纯文本格式存储体系结构信息：原子数、晶格矢量、原子位置以及其它可选的卡片数据（如 FORCE、MAGNETIC 等）。

该文件**不是** Fortran unformatted binary 格式，而是行式文本（ASCII），因此 C++ 读取时使用 `fopen("r")` + `fgets`，而非 `fread` + Fortran record 格式。

### 文件结构

```
natom                   ← 第一行：原子数（整数）
LATTICE                 ← 卡片关键词（可带前导/后导空格）
  a1_x    a1_y    a1_z  ← 晶格矢量 1 (Å)
  a2_x    a2_y    a2_z  ← 晶格矢量 2 (Å)
  a3_x    a3_y    a3_z  ← 晶格矢量 3 (Å)
POSITION                ← 卡片关键词
  species  x      y      z     imove_x imove_y imove_z  ← natom 行
  species  x      y      z     imove_x imove_y imove_z
  ...
FORCE                   ← 其它可选卡片（顺序自由）
  ...
MAGNETIC
  ...
```

### 布局说明

| 行 | 内容 | 说明 |
|----|------|------|
| 1 | `natom` | 整数，原子总数 |
| 2+ | 卡片及其数据 | 按关键词分组。已知关键词：`LATTICE`、`POSITION`、`FORCE`、`MAGNETIC` 等。**卡片间顺序自由** |

### 已知卡片

| 卡片关键词 | 数据行数 | 说明 |
|-----------|---------|------|
| `LATTICE` | 3 行 | 三个晶格矢量，每行 3 个 double，单位为 **Å** |
| `POSITION` | `natom` 行 | 每行：`species x y z imove_x imove_y imove_z`。其中 `x y z` 为**分数坐标**；`imove_*` 为 `1`/`0` 表示该方向是否可移动（读取时忽略）。|

### C++ 读取说明

`ATOM` 类（`io.ATOM` 模块）使用如下策略：

1. `fgets` 读取首行 → `std::strtol` 解析 `natom`
2. 卡片循环：逐行读取、**去除前后空格**、按关键词匹配
   - 匹配 `LATTICE` → 读取接下来 3 行，解析为晶格矩阵 → `Lattice(A, LengthUnit::Angstrom)`
   - 匹配 `POSITION` → 读取接下来 `natom` 行，每行用 `strtol`/`strtod` 链式解析 `species` 和分数坐标，跳过 `imove_*`
   - 其它行（未知卡片或其数据）→ 静默跳过
3. 卡片的顺序可以任意，解析器以关键词为驱动，不依赖固定顺序

> **注意**：目前仅实现 `LATTICE` 和 `POSITION` 两个卡片。其它卡片（FORCE、MAGNETIC 等）及其数据会被静默跳过。
>
> **单位**：晶格矢量以 Å 为单位读入，立即通过 `Lattice(..., LengthUnit::Angstrom)` 转换为 Bohr（原子单位制）。分数坐标不涉及单位转换。
