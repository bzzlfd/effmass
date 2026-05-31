# TL;DR

PWmat 大概有约定，WG ($W_G$) 与 GKK ($K=G+k$) 中的数据如何还原成平面波：

$$
\psi(r) = \sum_G W_G \exp(-i K \cdot r)
$$

PWmat `convert_wg2rho` 使用 forward FFT 做 G→R 变换。

**约定错误的典型症状**：如果 GKK 解释或 FFT 方向设反了，算出的 |ψ(r)|² 会与 PWmat reference 呈中心反演关系，即

```
my_WG2RHO[i, j, k] == PWmat_WG2RHO[-i, -j, -k]
```

看到这个 pattern 就说明 GKK 符号约定与 FFT 方向不匹配。

## 本代码的处理方式

1. `io.GKK` 读到文件数据后，将其**语义解释**为 $-K = -(G+k)$
    - 数据字节不变，仅改变解释
    - `KVecs::Cartesian` / `Spherical` 维持原样（表示 $-K$）
    - `IntegerView` 求整数索引时，使用点取反恢复 $G+k$，再计算 `iG = round(A·(G+k)/(2π) - k_frac)`
2. `math.fft` 模块采用 `R2G = FORWARD` / `G2R = BACKWARD`
    - 因为 R2G（实空间→G 空间）与 Forward（时域→频域）语义对应
    - `G2R` 即 backward，使用 `exp(+iGr)`
    - 结合第 1 点中对 GKK 数据的语义解释，端到端结果与 PWmat 一致
3. 测试 `test_wg_fft.cpp` 验证了：
    - `GKK::IntegerView` 给出正确的 G 整数索引
    - `FFT3D::G2R`（backward）给出与 PWmat `convert_wg2rho` 一致的 |ψ(r)|²
    - `Σ occ·|ψ|²` 复原 `OUT.RHO`，RMSE < 1e-2
4. FFT 归一化
    - Forward 和 Backward 基组内积为 n，为保证 $R = R \xrightarrow{\text{fft}} G \xrightarrow{\text{ifft}} R$，两系数乘积须为 $1/n$
    - 希望 $\sum_r |\psi(r)|^2 \frac{\Omega}{\Pi_i n_i} = 1$ 且 $\sum_g |\psi(g)|^2 \Omega = 1$
    - 所以 `R2G`（forward）系数 $1/n$，`G2R`（backward）系数 $1$


目前来看，这个约定只影响 FFT 时与 PWmat 生成数据的比对，其余情况从 GKK 读到 $q$ 平面波就当作 $\exp(iqr)$ 计算，管你正的负的。












# FFT Forward/Backward 

pocketfft 遵循标准 DFT 约定：

```
Forward:  F[k] = Σ_j f[j] · exp(-2πi · j·k / N)
Backward: f[j] = Σ_k F[k] · exp(+2πi · j·k / N)
```

（在未指定归一化方案之前，）
Forward 用 `exp(-iGr)` 做投影提取系数，Backward 反之。
换句话说，
1. Forward 假设数据用 `exp(+iGr)` 基表示
2. Backward 用 `exp(+iGr)` 做合成重建。
3. Backward 假设数据用 `exp(-iGr)` 基表示。


# PWmat convert_wg2rho.x

`docs/reference/convert_wg2rho/` 中存放了生成 `OUT.WG2RHO_1_16` 参考数据的源码。该程序将波函数从 G 空间变换到 R 空间，计算概率密度 `|ψ(r)|²`，以 RHO 格式输出。

## 源码文件

`*fft.f90` 文件中以 mode 选择 forward (mode>0) 还是 backward (mode<0)

而在进行 G->R 的 FFT 时，传了 `mode=1` 

### cfft.f90 

`cfft.f90:24-25` 明确注释：

```fortran
!c     mode > 0   forward
!c     mode < 0   backward
```

内部调用 (`cfft.f90:47-48`)：

```fortran
if (mode .gt. 0) call dcfftf(n3,wrk,wrk(np))    ! forward → dcfftf
if (mode .lt. 0) call dcfftb(n3,wrk,wrk(np))    ! backward → dcfftb
```

`dcfftf`/`dcfftb` 来自 FFTPACK，约定如下：

| FFTPACK 函数 | 数学操作 | 归一化 |
|-------------|---------|--------|
| `dcfftf` (forward) | `X_k = Σ_j x_j · exp(-iθ)` | 无 |
| `dcfftb` (backward) | `x_j = Σ_k X_k · exp(+iθ)` | 无 |

`cfft.f90:89-98` 在 `dcfftf`/`dcfftb` 返回后**额外**施加归一化：

```fortran
if (mode .gt. 0) return          ! forward: 不归一化，直接返回
factor = 1.D0 / DBLE(n1*n2*n3)   ! backward: 除以 (n1*n2*n3)
chdr = factor * chdr
chdi = factor * chdi
```

总结 `cfft` 的完整约定：

```
mode > 0 (forward):  f(G) = Σ_r f(r) · exp(-iGr)            无归一化
mode < 0 (backward): f(r) = (1/N) Σ_G F(G) · exp(+iGr)      有归一化 1/N
```

### convert_wg2rho 的调用

`convert_wg2rho.f90:236`：

```fortran
call cfft(n1,n2,n3,wave_r,wave_i,wrk,lwrk,1)
```

`mode=1 > 0` → **forward** 变换。此时 G 空间网格上已放置了 `W_G` 系数（`wave_r + i·wave_i`），经过 `cfft(..., 1)`：

```
ψ(r) = Σ_G W_G · exp(-iGr)        (无 1/N 因子)
```

即 **PWmat 的 G→R 变换使用 `exp(-iGr)`**。

### 对 GKK 的处理, G 矢量到网格索引的映射

`convert_wg2rho.f90:201-212`：

```fortran
aa1 = AL(1,1)*gkk_x(ig) + AL(2,1)*gkk_y(ig) + AL(3,1)*gkk_z(ig)  ! a1 · K
aa2 = AL(1,2)*gkk_x(ig) + AL(2,2)*gkk_y(ig) + AL(3,2)*gkk_z(ig)  ! a2 · K
aa3 = AL(1,3)*gkk_x(ig) + AL(2,3)*gkk_y(ig) + AL(3,3)*gkk_z(ig)  ! a3 · K

aa1 = aa1/(2*pi) + n1*2 + 0.1   ! 分数坐标，+2n 保证 mod 输入为正
aa2 = aa2/(2*pi) + n2*2 + 0.1
aa3 = aa3/(2*pi) + n3*2 + 0.1
i1 = aa1                         ! Fortran 截断取整
i1 = mod(i1, n1) + 1             ! → 1-based 网格索引
```

这里 $A^T \cdot K / (2\pi)$ 将笛卡尔的 $K = G + k$ 转为分数坐标 `iG + k_frac`，取整得到 G 矢量的整数索引。`+n*2+0.1` 保证 `mod` 的输入为正以正确处理负索引。


