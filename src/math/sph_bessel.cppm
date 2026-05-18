export module math.sph_bessel;

import std;

namespace detail {

    constexpr double BESSEL_EPS = 1e-15;
    constexpr double SERIES_THRESH = 0.5;
    constexpr double SERIES_TOL = 1e-16;
    constexpr int SERIES_MAX_ITER = 50;

    auto sphericalBesselJImpl(int l, double x) -> double {
        if (x < 0.0) x = -x;

        if (x < BESSEL_EPS) {
            return (l == 0) ? 1.0 : 0.0;
        }

        if (x < SERIES_THRESH) {
            double x2_half = x * x / 2.0;
            double term = 1.0;
            double sum = 1.0;
            for (int k = 1; k <= SERIES_MAX_ITER; ++k) {
                term *= -x2_half / (static_cast<double>(k) * (2.0 * l + 2.0 * k + 1.0));
                sum += term;
                if (std::abs(term) < SERIES_TOL * std::abs(sum)) break;
            }
            double dfact = 1.0;
            for (int i = 3; i <= 2 * l + 1; i += 2) dfact *= i;
            return sum * std::pow(x, l) / dfact;
        }

        double j0 = std::sin(x) / x;
        if (l == 0) return j0;

        double j1 = j0 / x - std::cos(x) / x;
        if (l == 1) return j1;

        double j_prev = j0;
        double j_curr = j1;
        for (int n = 1; n < l; ++n) {
            double j_next = (2.0 * n + 1.0) / x * j_curr - j_prev;
            j_prev = j_curr;
            j_curr = j_next;
        }
        return j_curr;
    }

} // namespace detail


export auto sphericalBesselJ(int l, double x) -> double {
    return detail::sphericalBesselJImpl(l, x);
}

export auto sphericalBesselJ(int l, std::span<const double> xs, std::span<double> out) -> void {
    int n = static_cast<int>(xs.size());

    if (l == 0) {
        for (int i = 0; i < n; ++i) {
            double x = xs[i];
            if (x < 0.0) x = -x;
            out[i] = (x < detail::BESSEL_EPS) ? 1.0 : (std::sin(x) / x);
        }
        return;
    }

    if (l == 1) {
        for (int i = 0; i < n; ++i) {
            double x = xs[i];
            if (x < 0.0) x = -x;
            if (x < detail::BESSEL_EPS) { out[i] = 0.0; continue; }
            double s = std::sin(x);
            double c = std::cos(x);
            out[i] = s / (x * x) - c / x;
        }
        return;
    }

    double dfact = 1.0;
    for (int j = 3; j <= 2 * l + 1; j += 2) dfact *= j;

    for (int i = 0; i < n; ++i) {
        double x = xs[i];
        if (x < 0.0) x = -x;

        if (x < detail::BESSEL_EPS) {
            out[i] = 0.0;
            continue;
        }

        if (x < detail::SERIES_THRESH) {
            double x2_half = x * x / 2.0;
            double term = 1.0;
            double sum = 1.0;
            for (int k = 1; k <= detail::SERIES_MAX_ITER; ++k) {
                term *= -x2_half / (static_cast<double>(k) * (2.0 * l + 2.0 * k + 1.0));
                sum += term;
                if (std::abs(term) < detail::SERIES_TOL * std::abs(sum)) break;
            }
            out[i] = sum * std::pow(x, l) / dfact;
        } else {
            double j0 = std::sin(x) / x;
            double j1 = j0 / x - std::cos(x) / x;
            double j_prev = j0;
            double j_curr = j1;
            for (int nu = 1; nu < l; ++nu) {
                double j_next = (2.0 * nu + 1.0) / x * j_curr - j_prev;
                j_prev = j_curr;
                j_curr = j_next;
            }
            out[i] = j_curr;
        }
    }
}


export auto radialIntegrate(
    std::span<const double> r,
    std::span<const double> weight,
    std::span<const double> f,
    int l,
    double q
) -> double {
    double result = 0.0;
    for (int i = 0; i < static_cast<int>(r.size()); ++i) {
        double qr = q * r[i];
        double jl = sphericalBesselJ(l, qr);
        result += r[i] * r[i] * f[i] * jl * weight[i];
    }
    return result;
}

