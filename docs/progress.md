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

## RealSphericalHarmonicsEngine::get_dtheta

当使用 RealSphericalHarmonics 时，Qlm 是必须的。
而缓存 Q_lm_dtheta_ 时，现在采用和 SphericalBesselJ::derivValue 一样的处理方式：
  只有在首次调用时，才会生成这些数据缓存。

计算公式

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
