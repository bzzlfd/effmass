module;

export module math.sph_harmonics;

import std;

export {
    class RealSphericalHarmonics;
}


// max l such that (2l)! fits in the integer type used in normFactorTiny.
//   int (32-bit):        l ≤ 6  (12! = 4.79e8)
//   long long (64-bit):  l ≤ 10 (20! = 2.43e18)
//   unsigned long long:  l ≤ 10 (20! = 2.43e18)
constexpr int MAX_L_SAFE = [] {
    unsigned long long fact = 1;
    int n = 0;
    while (fact <= std::numeric_limits<unsigned long long>::max() / (n + 1)) {
        ++n;
        fact *= n;
    }
    return n / 2;
}();

// expected l ≤ 3 (s/p/d/f orbitals)
constexpr int MAX_L_TINY = 3;

// reference table for l ≤ 3 (computed via normFactorTiny)
//
//   N_l^m = sqrt((2l+1) / (4*pi*prod)),  prod = prod_{k=l-m+1}^{l+m} k
//
//   l=0  m=0:  sqrt( 1/(4pi))
//   l=1  m=0:  sqrt( 3/(4pi))
//        m=1:  sqrt( 3/(8pi))
//   l=2  m=0:  sqrt( 5/(4pi))
//        m=1:  sqrt( 5/(24pi))
//        m=2:  sqrt( 5/(96pi))
//   l=3  m=0:  sqrt( 7/(4pi))
//        m=1:  sqrt( 7/(48pi))
//        m=2:  sqrt( 7/(480pi))
//        m=3:  sqrt(7/(2880pi))
constexpr double NORM_TABLE[4][4] = {
    {0.28209479177387814, 0.0,                  0.0,                  0.0},
    {0.48860251190291992, 0.3454941494713355,   0.0,                  0.0},
    {0.63078313050504009, 0.25751613468212636,  0.12875806734106318,  0.0},
    {0.7463526651802308,  0.21545345607610045,  0.068132365095552164, 0.027814921575518937},
};

// logarithmic fallback for l > MAX_L_SAFE
auto normFactorLog(int l, int m_abs) -> double {
    double log_res = 0.5 * (std::log(2.0 * l + 1.0) - std::log(4.0 * std::numbers::pi));
    for (int k = l - m_abs + 1; k <= l + m_abs; ++k) {
        log_res -= 0.5 * std::log(static_cast<double>(k));
    }
    return std::exp(log_res);
}

constexpr auto normFactorTiny(int l, int m_abs) -> double {
    long long prod = 1;
    for (int k = l - m_abs + 1; k <= l + m_abs; ++k) {
        prod *= k;
    }
    return std::sqrt((2.0 * l + 1.0) / (4.0 * std::numbers::pi * static_cast<double>(prod)));
}

constexpr auto& normFactorLarge = normFactorLog;

auto normFactor(int l, int m_abs) -> double {
    if (l <= MAX_L_TINY) {
        return NORM_TABLE[l][m_abs];
    } else if (l <= MAX_L_SAFE) {
        return normFactorTiny(l, m_abs);
    } else {
        return normFactorLog(l, m_abs);
    }
}

auto legendreP(int l, int m_abs, double theta) -> double {
    if (m_abs < 0 || m_abs > l) return 0.0;

    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);

    double p_mm = 1.0;
    if (m_abs > 0) {
        for (int i = 1; i <= m_abs; ++i) {
            p_mm *= (2.0 * i - 1.0) * sin_t;
        }
    }

    if (l == m_abs) return p_mm;

    double p_mm1 = cos_t * (2.0 * m_abs + 1.0) * p_mm;
    if (l == m_abs + 1) return p_mm1;

    double p_prev = p_mm;
    double p_curr = p_mm1;
    for (int n = m_abs + 2; n <= l; ++n) {
        double p_next = (cos_t * (2.0 * n - 1.0) * p_curr - (n + m_abs - 1.0) * p_prev) / (n - m_abs);
        p_prev = p_curr;
        p_curr = p_next;
    }
    return p_curr;
}


namespace archived {

auto realSphericalHarmonic(int l, int m, double theta, double phi) -> double {
    if (l < 0 || std::abs(m) > l) return 0.0;

    int m_abs = std::abs(m);

    double Plm = legendreP(l, m_abs, theta);
    double norm = normFactor(l, m_abs);

    if (m == 0) {
        return norm * Plm;
    } else if (m > 0) {
        return std::sqrt(2.0) * norm * Plm * std::cos(m * phi);
    } else {
        return std::sqrt(2.0) * norm * Plm * std::sin(m_abs * phi);
    }
}


// Compute all real spherical harmonics Y_{lm}(theta, phi) for l = 0..l_max.
// The output span must have size (l_max + 1)^2.
// Ordering within each l-block: m = -l, -l+1, ..., 0, ..., l-1, l.
// Flat index for (l, m) is: l*l + (m + l).
auto realSphericalHarmonics(double theta, double phi, int l_max, std::span<double> out) -> void {
    if (l_max < 0) {
        throw std::invalid_argument("realSphericalHarmonics: l_max must be non-negative");
    }
    std::size_t expected = static_cast<std::size_t>((l_max + 1) * (l_max + 1));
    if (out.size() != expected) {
        throw std::invalid_argument(
            std::format("realSphericalHarmonics: output span size {} != expected {}", out.size(), expected));
    }

    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);

    // Precompute cos(m*phi) and sin(m*phi)
    std::vector<double> cos_mphi(l_max + 1);
    std::vector<double> sin_mphi(l_max + 1);
    cos_mphi[0] = 1.0;
    sin_mphi[0] = 0.0;
    if (l_max > 0) {
        double c = std::cos(phi);
        double s = std::sin(phi);
        cos_mphi[1] = c;
        sin_mphi[1] = s;
        for (int m = 2; m <= l_max; ++m) {
            cos_mphi[m] = cos_mphi[m - 1] * c - sin_mphi[m - 1] * s;
            sin_mphi[m] = sin_mphi[m - 1] * c + cos_mphi[m - 1] * s;
        }
    }

    // Precompute normalization factors N_l^m
    std::vector<std::vector<double>> norms(l_max + 1);
    for (int l = 0; l <= l_max; ++l) {
        norms[l].resize(l + 1);
        for (int m = 0; m <= l; ++m) {
            norms[l][m] = normFactor(l, m);
        }
    }

    // Compute associated Legendre functions column-wise (fixed m, varying l)
    // and fill output immediately.
    for (int m = 0; m <= l_max; ++m) {
        std::vector<double> p(l_max - m + 1);

        // P_m^m(x)
        if (m == 0) {
            p[0] = 1.0;
        } else {
            double p_mm = 1.0;
            for (int i = 1; i <= m; ++i) {
                p_mm *= (2.0 * i - 1.0) * sin_t;
            }
            p[0] = p_mm;
        }

        // P_{m+1}^m(x)
        if (m < l_max) {
            p[1] = cos_t * (2.0 * m + 1.0) * p[0];
        }

        // Recurrence for higher l
        for (int l = m + 2; l <= l_max; ++l) {
            int idx = l - m;
            p[idx] = (cos_t * (2.0 * l - 1.0) * p[idx - 1] - (l + m - 1.0) * p[idx - 2]) / (l - m);
        }

        // Assemble real spherical harmonics for this m
        for (int l = m; l <= l_max; ++l) {
            int idx = l - m;
            double norm = norms[l][m];
            int base = l * l + l;

            if (m == 0) {
                out[base] = norm * p[idx];
            } else {
                double val = std::sqrt(2.0) * norm * p[idx];
                out[base - m] = val * sin_mphi[m];
                out[base + m] = val * cos_mphi[m];
            }
        }
    }
}


} // namespace archived (realSphericalHarmonic, realSphericalHarmonics)




// Vectorized spherical harmonics evaluator over K-point G-vectors.
// Precomputes all Y_lm for l <= l_max_resident on construction
// (or after a reset() call); computes larger l on demand.
class RealSphericalHarmonics {
public:
    explicit RealSphericalHarmonics(std::span<const double> theta, std::span<const double> phi, int l_max_resident)
        : theta_(theta), phi_(phi), ng_(theta_.size())
    {
        if (theta_.empty()) {
            throw std::invalid_argument(
                "RealSphericalHarmonics: theta and phi must be non-empty");
        }

        // Precompute cos(theta), sin(theta), cos(phi), sin(phi)
        // (immutable for these spans)
        cos_theta_.resize(ng_);
        sin_theta_.resize(ng_);
        cos_phi_.resize(ng_);
        sin_phi_.resize(ng_);
        for (std::size_t i = 0; i < ng_; ++i) {
            cos_theta_[i] = std::cos(theta_[i]);
            sin_theta_[i] = std::sin(theta_[i]);
        }
        for (std::size_t i = 0; i < ng_; ++i) {
            cos_phi_[i] = std::cos(phi_[i]);
            sin_phi_[i] = std::sin(phi_[i]);
        }

        reset(l_max_resident);
    }

    // Switch l_max_resident. Shrink truncates; expand continues Legendre
    // recurrence from saved top_p_ state (first call with empty top_p_
    // performs a full build from scratch).
    auto reset(int l_max_new) -> void {
        if (l_max_new < 0) {
            throw std::invalid_argument(
                "RealSphericalHarmonics::reset: l_max_resident must be non-negative");
        }
        if (l_max_new == l_max_resident_) return;

        if (l_max_new < l_max_resident_) {
            // -- Shrink: truncate buffers --
            std::size_t new_nm = static_cast<std::size_t>(l_max_new + 1);
            y_lm_.resize(new_nm * new_nm);
            top_p_.resize(new_nm);
            l_max_resident_ = l_max_new;
            return;
        }

        // -- Expand --
        int old_nm = l_max_resident_ + 1;
        int new_nm = l_max_new + 1;

        // Extend top_p_ for new m values
        top_p_.resize(static_cast<std::size_t>(new_nm));
        for (int m = old_nm; m <= l_max_new; ++m) {
            top_p_[static_cast<std::size_t>(m)].prev.resize(ng_);
            top_p_[static_cast<std::size_t>(m)].curr.resize(ng_);
        }

        // Extend y_lm_ cache
        y_lm_.resize(static_cast<std::size_t>(new_nm) * static_cast<std::size_t>(new_nm));
        for (auto& v : y_lm_) v.resize(ng_);

        std::vector<double> cm(ng_, 1.0), sm(ng_, 0.0);

        for (int m = 0; m <= l_max_new; ++m) {
            auto& tp = top_p_[static_cast<std::size_t>(m)];

            if (m > l_max_resident_) {
                // -- New m: full seed from scratch --
                seedPmm(m, tp.curr);
                assembleYlm(m, m, tp.curr, cm, sm);

                if (m == l_max_new) {
                    tp.prev = tp.curr;
                } else {
                    stepPmm1(m, tp.prev, tp.curr);
                    assembleYlm(m + 1, m, tp.curr, cm, sm);
                    if (m + 2 <= l_max_new)
                        advanceColumn(m, m + 2, l_max_new, tp.prev, tp.curr, cm, sm);
                }

            } else if (l_max_resident_ == m) {
                // -- Existing m at P_m^m only (old l_max == m) --
                if (m < l_max_new) {
                    stepPmm1(m, tp.prev, tp.curr);
                    assembleYlm(m + 1, m, tp.curr, cm, sm);
                    if (m + 2 <= l_max_new)
                        advanceColumn(m, m + 2, l_max_new, tp.prev, tp.curr, cm, sm);
                }

            } else {
                // -- Existing m with valid (P_{old-1}^m, P_{old}^m) --
                if (l_max_resident_ + 1 <= l_max_new)
                    advanceColumn(m, l_max_resident_ + 1, l_max_new, tp.prev, tp.curr, cm, sm);
            }

            // Advance cos(m*phi), sin(m*phi) to m+1
            if (m < l_max_new) {
                for (std::size_t i = 0; i < ng_; ++i) {
                    double old_cm = cm[i];
                    cm[i] = cm[i] * cos_phi_[i] - sm[i] * sin_phi_[i];
                    sm[i] = sm[i] * cos_phi_[i] + old_cm * sin_phi_[i];
                }
            }
        }

        l_max_resident_ = l_max_new;
    }

    // Return Y_lm for all G-vectors as a vector.
    // For l <= l_max_resident_: returns from the precomputed cache.
    // For l >  l_max_resident_: computes on demand from scratch
    //   (no persistent state between calls).
    auto operator()(int l, int m) -> std::vector<double> {
        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("RealSphericalHarmonics: invalid quantum numbers (l={}, m={})", l, m));
        }

        int m_abs = std::abs(m);

        // Cached path: l <= l_max_resident_
        if (l <= l_max_resident_) {
            int block = l * l + (m + l);
            return y_lm_[block];
        }

        // On-demand path: l > l_max_resident_
        // Seed Q_m^m and recurse upward to target l (no persistent state).
        // Q_l^m = N_l^m * P_l^m already includes normalization, so no extra
        // normFactor or multiplication is needed.
        std::vector<double> q_curr(ng_);
        double q0 = 0.5 / std::sqrt(std::numbers::pi);
        if (m_abs == 0) {
            for (std::size_t i = 0; i < ng_; ++i) q_curr[i] = q0;
        } else {
            double coeff = q0;
            for (int k = 1; k <= m_abs; ++k)
                coeff *= std::sqrt(static_cast<double>(2 * k + 1) / static_cast<double>(2 * k));
            for (std::size_t i = 0; i < ng_; ++i)
                q_curr[i] = coeff * std::pow(sin_theta_[i], m_abs);
        }

        if (l > m_abs) {
            std::vector<double> q_prev(ng_);
            // Q_{m+1}^m = sqrt(2m+3) * cosθ * Q_m^m
            double step_factor = std::sqrt(2.0 * m_abs + 3.0);
            for (std::size_t i = 0; i < ng_; ++i) {
                q_prev[i] = q_curr[i];
                q_curr[i] = step_factor * cos_theta_[i] * q_curr[i];
            }

            // Three-term recurrence for Q_l^m for n = m+2 .. l
            for (int n = m_abs + 2; n <= l; ++n) {
                double r1 = std::sqrt(static_cast<double>((2*n-1)*(2*n+1)*(n-m_abs))
                                     / static_cast<double>(n+m_abs));
                double r2 = (n + m_abs - 1.0)
                          * std::sqrt(static_cast<double>((2*n+1)*(n-m_abs)*(n-m_abs-1))
                                     / static_cast<double>((2*n-3)*(n+m_abs)*(n+m_abs-1)));
                double denom = 1.0 / static_cast<double>(n - m_abs);
                for (std::size_t i = 0; i < ng_; ++i) {
                    double next = (cos_theta_[i] * r1 * q_curr[i] - r2 * q_prev[i]) * denom;
                    q_prev[i] = q_curr[i];
                    q_curr[i] = next;
                }
            }
        }

        // Assemble Y_lm (Q already includes normalization)
        std::vector<double> result(ng_);
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i)
                result[i] = q_curr[i];
        } else {
            double sf = std::sqrt(2.0);
            if (m > 0) {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * q_curr[i] * std::cos(static_cast<double>(m) * phi_[i]);
            } else {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * q_curr[i] * std::sin(static_cast<double>(m_abs) * phi_[i]);
            }
        }

        return result;
    }

private:
    std::span<const double> theta_;
    std::span<const double> phi_;
    std::size_t ng_;
    int l_max_resident_ = -1;
    std::vector<double> cos_theta_;
    std::vector<double> sin_theta_;
    std::vector<double> cos_phi_;
    std::vector<double> sin_phi_;
    std::vector<std::vector<double>> y_lm_;


    // Top Q_l^m state for each m, saved so reset() expansion can
    // continue the recurrence without recomputing from scratch.
    // top_p_[m].curr = Q_{l_max_resident_}^{m},
    // top_p_[m].prev = Q_{l_max_resident_-1}^{m} (or = curr if l_max_resident_ == m).
    struct TopP {
        std::vector<double> prev;
        std::vector<double> curr;
    };
    std::vector<TopP> top_p_;

    //
    // -- Legendre recurrence helpers (shared by reset and operator()) --
    //

    // Seed Q_m^m into curr:  Q_0^0 = 1/(2√π), then
    //   Q_m^m = Q_0^0 · Π_{k=1}^{m} sqrt((2k+1)/(2k)) · sin^m(θ)
    // All intermediate values are O(1) (no overflow at high l).
    auto seedPmm(int m, std::vector<double>& curr) -> void {
        double q0 = 0.5 / std::sqrt(std::numbers::pi);
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i) curr[i] = q0;
        } else {
            double coeff = q0;
            for (int k = 1; k <= m; ++k) {
                coeff *= std::sqrt(static_cast<double>(2 * k + 1) / static_cast<double>(2 * k));
            }
            for (std::size_t i = 0; i < ng_; ++i)
                curr[i] = coeff * std::pow(sin_theta_[i], m);
        }
    }

    // One step: Q_{m+1}^m = sqrt(2m+3) · cosθ · Q_m^m
    // (the (2m+1) factor cancels when norm is folded into the recurrence)
    // On entry: curr = Q_m^m; on exit: prev = Q_m^m, curr = Q_{m+1}^m
    auto stepPmm1(int m, std::vector<double>& prev, std::vector<double>& curr) -> void {
        double factor = std::sqrt(2.0 * m + 3.0);
        for (std::size_t i = 0; i < ng_; ++i) {
            prev[i] = curr[i];
            curr[i] = factor * cos_theta_[i] * curr[i];
        }
    }

    // Assemble Y_{l,m} from Q_l^m(q), write into y_lm_ cache.
    // Q already includes the normalization N_l^m, so no extra normFactor.
    //   Y_{l0}   = Q_l^0
    //   Y_{l,±m} = √2 · Q_l^m · {cos(mφ), sin(mφ)}
    auto assembleYlm(int l, int m, const std::vector<double>& q,
                     const std::vector<double>& cm, const std::vector<double>& sm) -> void {
        int block = l * l + l;
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i)
                y_lm_[static_cast<std::size_t>(block)][i] = q[i];
        } else {
            double sf = std::sqrt(2.0);
            for (std::size_t i = 0; i < ng_; ++i) {
                double val = sf * q[i];
                y_lm_[static_cast<std::size_t>(block - m)][i] = val * sm[i];
                y_lm_[static_cast<std::size_t>(block + m)][i] = val * cm[i];
            }
        }
    }

    // Three-term recurrence for Q_l^m (normalized), fixed m, l = start_l .. end_l.
    // On entry: prev = Q_{start_l-2}^m, curr = Q_{start_l-1}^m
    // On exit:  prev = Q_{end_l-1}^m,   curr = Q_{end_l}^m
    //
    // Coefficients arise from folding the normalization ratio N_l^m / N_{l-1}^m
    // and N_l^m / N_{l-2}^m into the standard P recurrence:
    //
    //   c1 = sqrt((2l-1)(2l+1)(l-m)/(l+m))
    //   c2 = (l+m-1) · sqrt((2l+1)(l-m)(l-m-1)/((2l-3)(l+m)(l+m-1)))
    //   Q_l^m = (cosθ · c1 · Q_{l-1}^m - c2 · Q_{l-2}^m) / (l-m)
    auto advanceColumn(int m, int start_l, int end_l,
                       std::vector<double>& prev, std::vector<double>& curr,
                       const std::vector<double>& cm, const std::vector<double>& sm) -> void {
        for (int l = start_l; l <= end_l; ++l) {
            double r1 = std::sqrt(static_cast<double>((2 * l - 1) * (2 * l + 1) * (l - m))
                                 / static_cast<double>(l + m));
            double r2 = (l + m - 1.0)
                      * std::sqrt(static_cast<double>((2 * l + 1) * (l - m) * (l - m - 1))
                                 / static_cast<double>((2 * l - 3) * (l + m) * (l + m - 1)));
            double denom = 1.0 / static_cast<double>(l - m);
            for (std::size_t i = 0; i < ng_; ++i) {
                double next = (cos_theta_[i] * r1 * curr[i] - r2 * prev[i]) * denom;
                prev[i] = curr[i];
                curr[i] = next;
            }
            assembleYlm(l, m, curr, cm, sm);
        }
    }

};
