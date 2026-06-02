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


} // namespace archived


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
        // Seed P_m^m and recurse upward to target l (no persistent state)
        std::vector<double> p_curr(ng_);
        if (m_abs == 0) {
            for (std::size_t i = 0; i < ng_; ++i) p_curr[i] = 1.0;
        } else {
            double df = 1.0;
            for (int k = 1; k <= m_abs; ++k) df *= 2.0 * k - 1.0;
            for (std::size_t i = 0; i < ng_; ++i)
                p_curr[i] = df * std::pow(sin_theta_[i], m_abs);
        }

        if (l > m_abs) {
            std::vector<double> p_prev(ng_);
            // P_{m+1}^m = (2m+1) * cosθ * P_m^m
            for (std::size_t i = 0; i < ng_; ++i) {
                p_prev[i] = p_curr[i];
                p_curr[i] = (2.0 * m_abs + 1.0) * cos_theta_[i] * p_curr[i];
            }

            // General three-term recurrence for n = m+2 .. l
            for (int n = m_abs + 2; n <= l; ++n) {
                for (std::size_t i = 0; i < ng_; ++i) {
                    double next = ((2.0 * n - 1.0) * cos_theta_[i] * p_curr[i]
                                 - static_cast<double>(n + m_abs - 1) * p_prev[i])
                                 / static_cast<double>(n - m_abs);
                    p_prev[i] = p_curr[i];
                    p_curr[i] = next;
                }
            }
        }

        // Normalization factor
        double norm;
        if (l <= MAX_L_SAFE) {
            long long prod = 1;
            for (int k = l - m_abs + 1; k <= l + m_abs; ++k) prod *= k;
            norm = std::sqrt((2.0 * static_cast<double>(l) + 1.0)
                           / (4.0 * std::numbers::pi * static_cast<double>(prod)));
        } else {
            double log_norm = 0.5 * (std::log(2.0 * static_cast<double>(l) + 1.0)
                                    - std::log(4.0 * std::numbers::pi));
            for (int k = l - m_abs + 1; k <= l + m_abs; ++k)
                log_norm -= 0.5 * std::log(static_cast<double>(k));
            norm = std::exp(log_norm);
        }

        // Assemble Y_lm
        std::vector<double> result(ng_);
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i)
                result[i] = norm * p_curr[i];
        } else {
            double sf = std::sqrt(2.0) * norm;
            if (m > 0) {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * p_curr[i] * std::cos(static_cast<double>(m) * phi_[i]);
            } else {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * p_curr[i] * std::sin(static_cast<double>(m_abs) * phi_[i]);
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


    // Top P_l^m state for each m, saved so reset() expansion can
    // continue the recurrence without recomputing from scratch.
    // top_p_[m].curr = P_{l_max_resident_}^{m},
    // top_p_[m].prev = P_{l_max_resident_-1}^{m} (or = curr if l_max_resident_ == m).
    struct TopP {
        std::vector<double> prev;
        std::vector<double> curr;
    };
    std::vector<TopP> top_p_;

    //
    // -- Legendre recurrence helpers (shared by reset and operator()) --
    //

    // Seed P_m^m into curr:  P_m^m = (2m-1)!! · sin^m(θ)
    auto seedPmm(int m, std::vector<double>& curr) -> void {
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i) curr[i] = 1.0;
        } else {
            double df = 1.0;
            for (int k = 1; k <= m; ++k) df *= 2.0 * k - 1.0;
            for (std::size_t i = 0; i < ng_; ++i)
                curr[i] = df * std::pow(sin_theta_[i], m);
        }
    }

    // One step: P_{m+1}^m = (2m+1) · cosθ · P_m^m
    // On entry: curr = P_m^m; on exit: prev = P_m^m, curr = P_{m+1}^m
    auto stepPmm1(int m, std::vector<double>& prev, std::vector<double>& curr) -> void {
        for (std::size_t i = 0; i < ng_; ++i) {
            prev[i] = curr[i];
            curr[i] = (2.0 * m + 1.0) * cos_theta_[i] * curr[i];
        }
    }

    // Assemble Y_{l,m} from P_l^m(p), write into y_lm_ cache.
    auto assembleYlm(int l, int m, const std::vector<double>& p,
                     const std::vector<double>& cm, const std::vector<double>& sm) -> void {
        int block = l * l + l;
        double norm = normFactor(l, m);
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i)
                y_lm_[static_cast<std::size_t>(block)][i] = norm * p[i];
        } else {
            double sf = std::sqrt(2.0) * norm;
            for (std::size_t i = 0; i < ng_; ++i) {
                double val = sf * p[i];
                y_lm_[static_cast<std::size_t>(block - m)][i] = val * sm[i];
                y_lm_[static_cast<std::size_t>(block + m)][i] = val * cm[i];
            }
        }
    }

    // Three-term Legendre recurrence for fixed m, l = start_l .. end_l.
    // On entry: prev = P_{start_l-2}^m, curr = P_{start_l-1}^m
    // On exit:  prev = P_{end_l-1}^m,   curr = P_{end_l}^m
    // Each step assembles Y_{l,m} into the cache.
    auto advanceColumn(int m, int start_l, int end_l,
                       std::vector<double>& prev, std::vector<double>& curr,
                       const std::vector<double>& cm, const std::vector<double>& sm) -> void {
        for (int l = start_l; l <= end_l; ++l) {
            for (std::size_t i = 0; i < ng_; ++i) {
                double next = ((2.0 * l - 1.0) * cos_theta_[i] * curr[i]
                             - static_cast<double>(l + m - 1) * prev[i])
                             / static_cast<double>(l - m);
                prev[i] = curr[i];
                curr[i] = next;
            }
            assembleYlm(l, m, curr, cm, sm);
        }
    }

};
