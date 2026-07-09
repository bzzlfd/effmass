# TL;DR

PWmat 大概有约定，WG ($W_G$) 与 GKK ($GKK:=K=-k+G$) 中的数据如何还原成平面波：

$$
\psi(r) = \sum_G W_G \exp(-i K \cdot r)
$$

PWmat `convert_wg2rho` 使用 forward FFT 做 G→R 变换。

---

正确处理此符号约定问题，涉及
1. 本程序解析数据与 PWmat `RHO`和 `VR` 文件对应。
如果 GKK 解释或 FFT 方向设反了，算出的 |ψ(r)|² 会与 PWmat reference 呈中心反演关系，即

```
my_WG2RHO[i, j, k] == PWmat_WG2RHO[-i, -j, -k]
```
2. 通过 Hellmann–Feynman theorem 计算群速度与通过对 `EIGEN` 通过有限差分计算群速度之间相差符号。

## PWmat 平面波符号约定分析
1. 当我在 `IN.KPT` 输入 k 后， gen_G_comp 会选取截断平面波基组 $\{\tilde{+}K=\tilde{+}k\tilde{+}G\}$ 。在这一步中，截断标准实际上是 $\{(i,j,k)\,|\, (-k+(i,j,k)B)^2 < \text{Ecut} \}$ 。
2. `OUT.GKK` 中 $\{K\}$ 可以认为是以 $-k$ 为中心的球内
3. 想要得到实空间 `OUT.RHO` ，需对倒空间 GKK&WG 获得的数据进行 FFT.FORWARD：直接对K约化成分数坐标，然后取模得格点坐标 $(i,j,k)$，此坐标的系数对应实空间平面波 $\exp(-i \vec{(i,j,k)}\cdot r)$。
即约定

$$
\{(i,j,k) := reduce\&mod(GKK)\}
$$ 


$$
\{G:=(i,j,k)B \}
$$ 

$$
u(r) = \sum_G W_G \exp(-i G \cdot r)
$$。

对于以上观察，一个好的理解是，PWmat 约定对于一个平面波波矢 $K:=-k+G$（$K$ 直接对应 GKK 中所存储的数据），它代表平面波：
$
\exp(-i K \cdot r)
$。



## 本代码的处理方式

1. `io.GKK` 读取文件数据后，Kx/Ky/Kz 直接表示 $K := -(GKK)=k+G$
    - 本代码约定的 $G$ 与 PWmat 约定的 $G$ 相反；且 $\psi(r) = \sum_G W_G \exp(i K \cdot r)$。
    - `IntegerView` 求整数索引时，使用 $K = (G+k)$ 计算 `iG = round(A·(G+k)/(2π) - k_frac)`。
2. `math.fft` 模块采用 `R2G = FORWARD` / `G2R = BACKWARD`
    - 与 R2G（实空间→G 空间）- Forward（时域→频域）语义对应
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


