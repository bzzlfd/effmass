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

`OUT.EIGEN` 文件存储能量本征值和 k 点信息。

> TODO: 待补充详细的 record 结构和读取方法。

## OUT.KPT

`OUT.KPT` 文件以**文本格式**存储 k 点坐标。

> TODO: 待补充具体的列定义和解析方法。
