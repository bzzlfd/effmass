# d Ylm

对于以球坐标为自变量的函数 函数 f(q, θ, φ) := β(q) Y_lm(θ, φ), 它的梯度为

(∂f/∂q, 1/q ∂f/∂θ, 1/(q sinθ) ∂f/∂φ)

这样一来，如果我对 (∂/∂x, ∂/∂y, ∂/∂z) 感兴趣，我只需要对所得的梯度进行旋转：

(∂f/∂x)   (sinθ cosφ, cosθ cosφ, -sinφ) (∂β/∂q Y_lm(θ, φ)      )
(∂f/∂y) = (sinθ sinφ, cosθ sinφ, cosφ ) (β/q   ∂Y_lm/∂θ        )
(∂f/∂z)   (cosθ,      -sinθ    , 0    ) (β/q   (1/sin0)∂Y_lm/∂φ)


1/q 可以吸收到 β(q) 处理，当 q→0 时，计算 β'(0)。
我们在这里仅处理球谐函数（变量为 θ, φ 的部分）Y_lm = Qlm e^(imφ) 。

## RealSphericalHarmonicsEngine::get

get 时对 Qlm 的缓存计算已经完成。

## RealSphericalHarmonicsEngine::get_grad_theta

计算 `∂Y_lm/∂θ`。

新缓存变量 `Q_lm_grad_theta_`（`std::optional<std::vector<std::vector<double>>>`），惰性填充。
使用 `GradThetaRecurrence` 结构体，通过**对 QRecurrence 公式逐项求导**做逐列递推（seed → step1 → advance），
依赖常驻缓存 `Q_lm_`。

### 递推公式（对 QRecurrence 求导）

记 `C_m = (-1)^m/(2√π) · ∏_{k=1}^m √((2k+1)/(2k))`：

| 步骤 | Q_l^m | dQ_l^m/dθ |
|------|-------|-----------|
| **seed** | `Q_m^m = C_m · sin^m(θ)` | `dQ_m^m/dθ = m · C_m · sin^{m-1}(θ) · cosθ`（m=0→0, m=1→C_1·cosθ 有限） |
| **step1** | `Q_{m+1}^m = √(2m+3) · cosθ · Q_m^m` | `dQ_{m+1}^m/dθ = √(2m+3) · (-sinθ·Q_m^m + cosθ·dQ_m^m/dθ)` |
| **advance** | `Q_l^m = (cosθ·r1·Q_{l-1}^m - r2·Q_{l-2}^m)/(l-m)` | `dQ_l^m/dθ = (cosθ·r1·dQ_{l-1}^m/dθ - sinθ·r1·Q_{l-1}^m - r2·dQ_{l-2}^m/dθ)/(l-m)` |

其中 r1, r2 系数与 `QRecurrence::advance` 完全相同。

`get_grad_theta` 返回 `∂Y_lm/∂θ`：

| m | ∂Y_lm/∂θ |
|---|----------|
| 0 | **dQ_l^0/dθ** |
| >0 | **√2 · dQ_l^m/dθ · cos(mφ)** |
| <0 | **√2 · dQ_l^{\|m\|}/dθ · sin(\|m\|φ)** |


## RealSphericalHarmonicsEngine::get_grad_phi

已实现（2026-07-03）。计算 `(1/sinθ)·∂Y_lm/∂φ`。

需要新缓存变量 `Q_lm_grad_phi_`（`std::optional<std::vector<std::vector<double>>>`），惰性填充。
使用 `GradPhiRecurrence` 结构体做递推（`QRecurrence` 代码复用，seed 不同）。

### 计算公式

记 `R_l^m = Q_l^m / sinθ`，`get_grad_phi` 返回 `(1/sinθ)·∂Y_lm/∂φ`：

| m | Y_lm (real) | (1/sinθ)·∂Y_lm/∂φ |
|---|-------------|---------------------|
| 0 | Q_l^0 | **0** |
| >0 | √2 · Q_l^m · cos(mφ) | **-m · √2 · R_l^m · sin(mφ)** = -m/sinθ · Y_{l,-m} |
| <0 | √2 · Q_l^{\|m\|} · sin(\|m\|φ) | **\|m\| · √2 · R_l^{\|m\|} · cos(\|m\|φ)** = \|m\|/sinθ · Y_{l,\|m\|} |

相比于旧的 `get_dphi`（返回纯 ∂Y/∂φ），区别在于 θ=0 处：|m|=1 时 `get_grad_phi` 有限非零，而纯 ∂Y/∂φ 为 0。

### 缓存策略

- 使用 `std::optional` 存储完整的 `Q_lm_grad_phi_`：无数据 = std::nullopt，有数据 = engaged。
- `reinit()` → `.reset()` 清空；`reserveNg()` → 条件式 reserve；`setLMax()` → 惰性 no-op。
- 详见 `docs/design/real_spherical_harmonics.md`。


## 讨论

- m=0 时 `get_grad_theta` 返回**非零值**（l≥1），而 `get_grad_phi` 恒 0。
- θ=0/π 边界：
    - `get_grad_theta` ：l=1,m=0 → 0；l=1,m=±1 → 有限非零（`-sqrt(3/(4π))·cosφ` 或 `-sqrt(3/(4π))·sinφ`）；|m|≥2 → 0。
    - l=1,m=±1 → 有限非零（`sqrt(3/(4π))·sinφ` 或 `-sqrt(3/(4π))·cosφ`）；|m|≥2 → 0。
    - 等等，l=1,m=±1 时岂不是 grad_theta/phi 的奇点，沿着不同经线趋向于 0/π 有不同的结果？
        1. 毛球定理。有奇点就对了。
        2. 它们组合出的笛卡尔坐标系下的梯度向量是确定的。

