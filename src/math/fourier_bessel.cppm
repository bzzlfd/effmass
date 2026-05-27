export module math.fourier_bessel;

import std;
import math.sph_bessel;

export {
    enum class RadialMeshType : int;
    class BetaQInterpolator;
}


enum class RadialMeshType {
    General,
    Uniform,
};

constexpr double FOUR_PI = 4.0 * std::numbers::pi;

// Simpson integration on the index i (dr/di encoded in rab).
// For a Uniform mesh, rab(i) is constant (dr) and the formula is simplified.
// Handles even number of points without dropping the last point by combining
// Simpson 3/8 rule for the first 4 points and Simpson 1/3 for the remainder.
auto simpson(std::span<const double> f, std::span<const double> rab,
             RadialMeshType mesh_type) -> double {
    int n = static_cast<int>(f.size());
    int nrab = static_cast<int>(rab.size());
    if (n == 0) return 0.0;
    if (n == 1) return f[0] * rab[0];
    if (nrab < n) {
        throw std::invalid_argument("simpson: rab size must be >= f size");
    }

    if (mesh_type == RadialMeshType::Uniform) {
        double dr = rab[0];
        if (n % 2 == 1) {
            // Odd number of points: standard composite Simpson 1/3
            double sum = f[0] + f[n - 1];
            for (int i = 1; i < n - 1; ++i) {
                sum += (i % 2 == 0 ? 2.0 : 4.0) * f[i];
            }
            return dr * sum / 3.0;
        } else {
            if (n == 2) {
                return dr * 0.5 * (f[0] + f[1]);
            }
            // Even number of points >= 4:
            // Simpson 3/8 for indices 0..3 (3 intervals),
            // Simpson 1/3 for indices 3..n-1 (n-3 points, odd).
            double sum38 = f[0] + 3.0 * f[1] + 3.0 * f[2] + f[3];
            double sum13 = f[3] + f[n - 1];
            for (int i = 4; i < n - 1; ++i) {
                sum13 += (i % 2 == 0 ? 4.0 : 2.0) * f[i];
            }
            return dr * (3.0 / 8.0 * sum38 + sum13 / 3.0);
        }
    } else {
        // General mesh: g(i) = f(i) * rab(i), integrate over uniform index i.
        if (n % 2 == 1) {
            double sum = f[0] * rab[0] + f[n - 1] * rab[n - 1];
            for (int i = 1; i < n - 1; ++i) {
                sum += (i % 2 == 0 ? 2.0 : 4.0) * f[i] * rab[i];
            }
            return sum / 3.0;
        } else {
            if (n == 2) {
                return 0.5 * (f[0] * rab[0] + f[1] * rab[1]);
            }
            double g0 = f[0] * rab[0];
            double g1 = f[1] * rab[1];
            double g2 = f[2] * rab[2];
            double g3 = f[3] * rab[3];
            double sum38 = g0 + 3.0 * g1 + 3.0 * g2 + g3;
            double sum13 = g3 + f[n - 1] * rab[n - 1];
            for (int i = 4; i < n - 1; ++i) {
                sum13 += (i % 2 == 0 ? 4.0 : 2.0) * f[i] * rab[i];
            }
            return 3.0 / 8.0 * sum38 + sum13 / 3.0;
        }
    }
}

auto integrateBetaQ(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> beta,
    int l,
    double q,
    RadialMeshType mesh_type
) -> double {
    int n = static_cast<int>(beta.size());
    if (n == 0) return 0.0;

    int nr = static_cast<int>(r.size());
    int nrab = static_cast<int>(rab.size());
    int nmax = std::min({n, nr, nrab});
    // Defensive: beta is already truncated when the UPF object is constructed,
    // but if there are still trailing zeros within the provided span,
    // shrink nmax further to avoid unnecessary j_l evaluations.
    while (nmax > 0 && beta[static_cast<std::size_t>(nmax) - 1] == 0.0) {
        --nmax;
    }
    if (nmax == 0) return 0.0;

    // Build integrand: beta(r) * r * j_l(qr)
    // For q=0 and l=0, j_0(0)=1; for l>0, j_l(0)=0.
    std::vector<double> integrand(nmax);
    if (q < 1e-15) {
        if (l == 0) {
            for (int i = 0; i < nmax; ++i) {
                integrand[i] = beta[i] * r[i];
            }
        } else {
            return 0.0;
        }
    } else {
        SphericalBesselJ bessel{r.subspan(0, static_cast<std::size_t>(nmax)), q};
        bessel.advance(l);
        auto jl = bessel.value();
        for (int i = 0; i < nmax; ++i) {
            integrand[i] = beta[i] * r[i] * jl[i];
        }
    }

    return FOUR_PI * simpson(integrand, rab.subspan(0, nmax), mesh_type);
}

// Batch version: evaluate at multiple q points while reusing the
// SphericalBesselJ buffers (r_ / j_curr_ / j_next_) and the integrand
// vector across reset() calls, avoiding repeated heap allocations.
auto integrateBetaQ(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> beta,
    int l,
    std::span<const double> qs,
    std::span<double> out,
    RadialMeshType mesh_type
) -> void {
    int nq = static_cast<int>(qs.size());
    if (static_cast<std::size_t>(nq) != out.size()) {
        throw std::invalid_argument("integrateBetaQ: qs and out size mismatch");
    }
    if (nq == 0) return;

    int n = static_cast<int>(beta.size());
    if (n == 0) {
        std::fill(out.begin(), out.end(), 0.0);
        return;
    }

    int nr = static_cast<int>(r.size());
    int nrab = static_cast<int>(rab.size());
    int nmax = std::min({n, nr, nrab});
    while (nmax > 0 && beta[static_cast<std::size_t>(nmax) - 1] == 0.0) {
        --nmax;
    }
    if (nmax == 0) {
        std::fill(out.begin(), out.end(), 0.0);
        return;
    }

    auto r_sub = r.subspan(0, static_cast<std::size_t>(nmax));
    auto rab_sub = rab.subspan(0, static_cast<std::size_t>(nmax));

    std::vector<double> integrand(nmax);
    SphericalBesselJ bessel{r_sub, qs.front()};

    for (int iq = 0; iq < nq; ++iq) {
        bessel.reset(qs[iq]);
        bessel.advance(l);
        auto jl = bessel.value();
        for (int i = 0; i < nmax; ++i) {
            integrand[i] = beta[i] * r[i] * jl[i];
        }
        out[iq] = FOUR_PI * simpson(integrand, rab_sub, mesh_type);
    }
}

// Lagrange 4-point cubic interpolation (same as QE interp_beta)
auto lagrangeCubic(
    double f0, double f1, double f2, double f3,
    double px
) -> double {
    double ux = 1.0 - px;
    double vx = 2.0 - px;
    double wx = 3.0 - px;
    return f0 * ux * vx * wx / 6.0 +
           f1 * px * vx * wx / 2.0 -
           f2 * px * ux * wx / 2.0 +
           f3 * px * ux * vx / 6.0;
}

// Derivative of the same 4-point cubic interpolant (same as QE interp_dbeta)
auto lagrangeCubicDerivative(
    double f0, double f1, double f2, double f3,
    double px, double dq
) -> double {
    double ux = 1.0 - px;
    double vx = 2.0 - px;
    double wx = 3.0 - px;
    return (f0 * (-vx*wx - ux*wx - ux*vx) / 6.0 +
            f1 * ( vx*wx - px*wx - px*vx) / 2.0 -
            f2 * ( ux*wx - px*wx - px*ux) / 2.0 +
            f3 * ( ux*vx - px*vx - px*ux) / 6.0) / dq;
}


namespace archived {

// Single-q Fourier-Bessel transform of a radial beta function.
// Returns: 4*pi * integral_0^R  beta(r) * r * j_l(q*r) * dr
auto fourierBesselBeta(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> beta,
    int l,
    double q,
    RadialMeshType mesh_type = RadialMeshType::General
) -> double {
    double out = 0.0;
    integrateBetaQ(r, rab, beta, l,
                           std::span<const double>(&q, 1),
                           std::span<double>(&out, 1), mesh_type);
    return out;
}

// Batch version: evaluate at multiple q points.
auto fourierBesselBeta(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> beta,
    int l,
    std::span<const double> qs,
    std::span<double> out,
    RadialMeshType mesh_type = RadialMeshType::General
) -> void {
    integrateBetaQ(r, rab, beta, l, qs, out, mesh_type);
}


} // namespace archived


// Tabulated interpolator for beta(q) = 4*pi * integral beta(r) * r * j_l(qr) * dr.
// Uses 4-point Lagrange cubic interpolation, matching QE's interp_beta / interp_dbeta.
class BetaQInterpolator {
public:
    BetaQInterpolator(
        std::span<const double> r,
        std::span<const double> rab,
        std::span<const double> beta,
        int l,
        double dq = 0.01,
        double q_max = 50.0,
        RadialMeshType mesh_type = RadialMeshType::General
    ) : dq_(dq), q_max_(q_max), l_(l), mesh_type_(mesh_type) {
        if (dq <= 0.0) throw std::invalid_argument("BetaQInterpolator: dq must be positive");
        if (q_max < 0.0) throw std::invalid_argument("BetaQInterpolator: q_max must be non-negative");

        // Build table up to q_max + 3*dq so that cubic interpolation
        // has full 4-point support throughout [0, q_max].
        int n = static_cast<int>(std::ceil(q_max / dq)) + 1 + 3;
        table_.resize(n);
        for (int i = 0; i < n; ++i) {
            double q = i * dq;
            table_[i] = integrateBetaQ(r, rab, beta, l, q, mesh_type);
        }
    }

    // Function value at q using centered 4-point Lagrange cubic interpolation.
    auto operator()(double q) const -> double {
        if (q < 0.0) throw std::domain_error("BetaQInterpolator: q must be non-negative");
        if (q > q_max_) throw std::domain_error("BetaQInterpolator: q exceeds q_max");
        if (q == 0.0) return table_.front();

        double s = q / dq_;
        int i0 = static_cast<int>(std::floor(s)) - 1;
        double px = s - static_cast<double>(i0);

        // Left boundary: fall back to forward stencil when q < dq
        if (i0 < 0) {
            i0 = 0;
            px = s;
        }

        return lagrangeCubic(
            table_[i0], table_[i0 + 1], table_[i0 + 2], table_[i0 + 3], px);
    }

    // Derivative d/dq of beta(q) using centered 4-point Lagrange cubic interpolation.
    auto derivative(double q) const -> double {
        if (q < 0.0) throw std::domain_error("BetaQInterpolator: q must be non-negative");
        if (q > q_max_) throw std::domain_error("BetaQInterpolator: q exceeds q_max");
        if (q == 0.0) return 0.0;

        double s = q / dq_;
        int i0 = static_cast<int>(std::floor(s)) - 1;
        double px = s - static_cast<double>(i0);

        // Left boundary: fall back to forward stencil when q < dq
        if (i0 < 0) {
            i0 = 0;
            px = s;
        }

        return lagrangeCubicDerivative(
            table_[i0], table_[i0 + 1], table_[i0 + 2], table_[i0 + 3], px, dq_);
    }

    auto table() const -> std::span<const double> { return table_; }
    auto step() const -> double { return dq_; }
    auto maxQ() const -> double { return q_max_; }
    auto angularMomentum() const -> int { return l_; }

private:
    int l_;
    double dq_;
    double q_max_;
    RadialMeshType mesh_type_;
    std::vector<double> table_;
};


namespace archived {

// Estimate the maximum interpolation error by evaluating at half-integer
// multiples of the table step and comparing with the direct integral.
// Returns the maximum absolute deviation.
auto estimateInterpolationError(
    const BetaQInterpolator& interp,
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> beta,
    int l,
    int n_check = 200,
    RadialMeshType mesh_type = RadialMeshType::General
) -> double {
    if (n_check <= 0) return 0.0;

    double max_err = 0.0;
    double q_max = interp.maxQ();
    double dq = interp.step();

    for (int i = 0; i < n_check; ++i) {
        double q = (i + 0.5) * dq;
        if (q > q_max) break;
        double exact = integrateBetaQ(r, rab, beta, l, q, mesh_type);
        double approx = interp(q);
        double err = std::abs(exact - approx);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

} // namespace archived
