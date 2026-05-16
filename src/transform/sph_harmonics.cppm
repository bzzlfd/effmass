export module transform.sph_harmonics;

import std;

namespace transform {

namespace detail {

    auto legendreP(int l, int m_abs, double x) -> double {
        if (std::abs(x) > 1.0) x = std::copysign(1.0, x);

        if (m_abs < 0 || m_abs > l) return 0.0;

        if (m_abs > 0 && std::abs(x) > 0.999999999999) return 0.0;

        double somx2 = (1.0 - x) * (1.0 + x);

        double p_mm = 1.0;
        if (m_abs > 0) {
            double factor = std::sqrt(somx2);
            for (int i = 1; i <= m_abs; ++i) {
                p_mm *= (2.0 * i - 1.0) * factor;
            }
        }

        if (l == m_abs) return p_mm;

        double p_mm1 = x * (2.0 * m_abs + 1.0) * p_mm;
        if (l == m_abs + 1) return p_mm1;

        double p_prev = p_mm;
        double p_curr = p_mm1;
        for (int n = m_abs + 2; n <= l; ++n) {
            double p_next = (x * (2.0 * n - 1.0) * p_curr - (n + m_abs - 1.0) * p_prev) / (n - m_abs);
            p_prev = p_curr;
            p_curr = p_next;
        }
        return p_curr;
    }

    auto normFactor(int l, int m_abs) -> double {
        double result = std::sqrt((2.0 * l + 1.0) / (4.0 * std::numbers::pi));
        for (int k = l - m_abs + 1; k <= l + m_abs; ++k) {
            result /= std::sqrt(static_cast<double>(k));
        }
        return result;
    }

} // namespace detail


export auto realSphericalHarmonic(int l, int m, double theta, double phi) -> double {
    if (l < 0 || std::abs(m) > l) return 0.0;

    double x = std::cos(theta);
    int m_abs = std::abs(m);

    double Plm = detail::legendreP(l, m_abs, x);
    double norm = detail::normFactor(l, m_abs);

    if (m == 0) {
        return norm * Plm;
    } else if (m > 0) {
        return std::sqrt(2.0) * norm * Plm * std::cos(m * phi);
    } else {
        return std::sqrt(2.0) * norm * Plm * std::sin(m_abs * phi);
    }
}

} // namespace transform
