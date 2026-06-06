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


// Vectorized spherical harmonics evaluator over K-point G-vectors.
// Precomputes all Y_lm for l <= l_max_resident on construction
// (or after a reset() call); computes larger l on demand.
class RealSphericalHarmonics {
public:
    enum class CacheMode { Full, None };

    explicit RealSphericalHarmonics(std::span<const double> theta, std::span<const double> phi,
                                     int l_max_resident, CacheMode mode = CacheMode::Full)
        : theta_(theta), phi_(phi), ng_(theta_.size()), mode_(mode)
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

        if (mode_ == CacheMode::Full) reset(l_max_resident);
    }

    // Switch l_max_resident. Shrink truncates; expand continues Legendre
    // recurrence from saved top_column_Q_ state (first call with empty top_column_Q_
    // performs a full build from scratch).
    auto reset(int l_max_new) -> void {
        if (l_max_new < 0) {
            throw std::invalid_argument(
                "RealSphericalHarmonics::reset: l_max_resident must be non-negative");
        }
        if (mode_ == CacheMode::None) return;
        if (l_max_new == l_max_resident_) return;

        // -- Shrink: truncate buffers and back-iterate top_column_Q_ --
        if (l_max_new < l_max_resident_) {
            int old_l_max = l_max_resident_;
            std::size_t new_nm = static_cast<std::size_t>(l_max_new + 1);
            y_lm_.resize(new_nm * new_nm);
            top_column_Q_.resize(new_nm);
            l_max_resident_ = l_max_new;

            // Backward Q recurrence: recover top_column_Q_ state for new l_max.
            // Forward 3-term:  Q_l^m = (cosθ·c1·Q_{l-1}^m - c2·Q_{l-2}^m) / (l-m)
            // Inverse:          Q_{l-2}^m = (cosθ·c1·Q_{l-1}^m - (l-m)·Q_l^m) / c2
            // Valid for l ≥ m+2 (c2 ≠ 0).
            for (int m = 0; m <= l_max_new; ++m) {
                auto& tp = top_column_Q_[static_cast<std::size_t>(m)];

                if (l_max_new == m) {
                    // Seed level: use direct formula to avoid cosθ division.
                    seedPmm(m, tp.curr);
                    tp.prev = tp.curr;
                } else {
                    for (int l = old_l_max; l > l_max_new; --l)
                        retreatColumn(m, l, tp.prev, tp.curr);
                }
            }
            return;
        }

        // -- Expand --
        int old_nm = l_max_resident_ + 1;
        int new_nm = l_max_new + 1;

        // Extend top_column_Q_ for new m values
        top_column_Q_.resize(static_cast<std::size_t>(new_nm));
        for (int m = old_nm; m <= l_max_new; ++m) {
            top_column_Q_[static_cast<std::size_t>(m)].prev.resize(ng_);
            top_column_Q_[static_cast<std::size_t>(m)].curr.resize(ng_);
        }

        // Extend y_lm_ cache
        y_lm_.resize(static_cast<std::size_t>(new_nm) * static_cast<std::size_t>(new_nm));
        for (auto& v : y_lm_) v.resize(ng_);

        std::vector<double> cm(ng_, 1.0), sm(ng_, 0.0);

        for (int m = 0; m <= l_max_new; ++m) {
            auto& tp = top_column_Q_[static_cast<std::size_t>(m)];

            if (m > l_max_resident_) {
                // -- New m: full seed from scratch --
                seedPmm(m, tp.curr);
                assembleYlm(m, m, tp.curr, cm, sm);

                if (m == l_max_new) {
                    tp.prev = tp.curr;
                } else {
                    stepPmm1(m, tp.prev, tp.curr);
                    assembleYlm(m + 1, m, tp.curr, cm, sm);
                    for (int l = m + 2; l <= l_max_new; ++l) {
                        advanceColumn(m, l, tp.prev, tp.curr);
                        assembleYlm(l, m, tp.curr, cm, sm);
                    }
                }

            } else if (l_max_resident_ == m) {
                // -- Existing m at P_m^m only (old l_max == m) --
                if (m < l_max_new) {
                    stepPmm1(m, tp.prev, tp.curr);
                    assembleYlm(m + 1, m, tp.curr, cm, sm);
                    for (int l = m + 2; l <= l_max_new; ++l) {
                        advanceColumn(m, l, tp.prev, tp.curr);
                        assembleYlm(l, m, tp.curr, cm, sm);
                    }
                }

            } else {
                // -- Existing m with valid (P_{old-1}^m, P_{old}^m) --
                for (int l = l_max_resident_ + 1; l <= l_max_new; ++l) {
                    advanceColumn(m, l, tp.prev, tp.curr);
                    assembleYlm(l, m, tp.curr, cm, sm);
                }
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

    // Return Y_lm for all G-vectors from the precomputed cache.
    // Requires CacheMode::Full and l <= l_max_resident.
    auto get(int l, int m) -> const std::vector<double>& {
        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("RealSphericalHarmonics::get(): invalid quantum numbers (l={}, m={})", l, m));
        }
        if (mode_ == CacheMode::None) {
            throw std::runtime_error(
                "RealSphericalHarmonics::get(): not available in CacheMode::None, use compute()");
        }
        if (l > l_max_resident_) {
            throw std::runtime_error(
                std::format("RealSphericalHarmonics::get(): l={} > l_max_resident_={}, call reset({}) first, or use compute()", l, l_max_resident_, l));
        }

        int block = l * l + (m + l);
        return y_lm_[static_cast<std::size_t>(block)];
    }

    // Convenience: dispatches based on CacheMode.
    //   Full → get() (cached, by copy)
    //   None → compute() (on-the-fly, by move)
    auto operator()(int l, int m) -> std::vector<double> {
        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("RealSphericalHarmonics::operator(): invalid quantum numbers (l={}, m={})", l, m));
        }
        if (mode_ == CacheMode::Full) {
            if (l > l_max_resident_) {
                throw std::runtime_error(
                    std::format("RealSphericalHarmonics::operator(): l={} > l_max_resident_={}, call reset({}) first, or use compute()", l, l_max_resident_, l));
            }
            return get(l, m);
        } else {
            return compute(l, m);
        }
    }

    // Compute Y_lm on the fly. Works in any mode, for any (l,m).
    // Returns by value (no persistent cache).
    // Reuses the Q recurrence helpers (seedPmm, stepPmm1, advanceColumn);
    // only the final Y_lm assembly differs (returns a new vector).
    auto compute(int l, int m) -> std::vector<double> {
        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("RealSphericalHarmonics::compute: invalid quantum numbers (l={}, m={})", l, m));
        }

        int m_abs = std::abs(m);

        // Q_l^m = N_l^m * P_l^m already includes normalization
        std::vector<double> q_curr(ng_);
        seedPmm(m_abs, q_curr);

        if (l > m_abs) {
            std::vector<double> q_prev(ng_);
            stepPmm1(m_abs, q_prev, q_curr);
            for (int n = m_abs + 2; n <= l; ++n)
                advanceColumn(m_abs, n, q_prev, q_curr);
        }

        std::vector<double> result(ng_);
        if (m == 0) {
            result = std::move(q_curr);
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
    CacheMode mode_ = CacheMode::Full;
    std::vector<double> cos_theta_;
    std::vector<double> sin_theta_;
    std::vector<double> cos_phi_;
    std::vector<double> sin_phi_;
    std::vector<std::vector<double>> y_lm_;


    // Top Q_l^m state for each m, saved so reset() expansion can
    // continue the recurrence without recomputing from scratch.
    // top_column_Q_[m].curr = Q_{l_max_resident_}^{m},
    // top_column_Q_[m].prev = Q_{l_max_resident_-1}^{m} (or = curr if l_max_resident_ == m).
    struct TopColumnQ {
        std::vector<double> prev;
        std::vector<double> curr;
    };
    std::vector<TopColumnQ> top_column_Q_;

    //
    // -- Legendre recurrence helpers --
    //

    // Seed Q_m^m into curr:  Q_0^0 = 1/(2√π), then
    //   Q_m^m = Q_0^0 · Π_{k=1}^{m} sqrt((2k+1)/(2k)) · sin^m(θ)
    // All intermediate values are O(1) (no overflow at high l).
    auto seedPmm(int m, std::vector<double>& curr) -> void {
        double Q0 = 0.5 / std::sqrt(std::numbers::pi);
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i) curr[i] = Q0;
        } else {
            double coeff = Q0;
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

    // One step of forward three-term recurrence for Q_l^m (normalized).
    // On entry: prev = Q_{l-2}^m, curr = Q_{l-1}^m
    // On exit:  prev = Q_{l-1}^m, curr = Q_l^m
    //
    //   r1 = sqrt((2l-1)(2l+1)(l-m)/(l+m))
    //   r2 = (l+m-1) · sqrt((2l+1)(l-m)(l-m-1)/((2l-3)(l+m)(l+m-1)))
    //   Q_l^m = (cosθ · r1 · Q_{l-1}^m - r2 · Q_{l-2}^m) / (l-m)
    auto advanceColumn(int m, int l,
                       std::vector<double>& prev, std::vector<double>& curr) -> void {
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
    }

    // One step of inverse three-term recurrence for Q_l^m.
    // On entry: prev = Q_{l-1}^m, curr = Q_l^m
    // On exit:  prev = Q_{l-2}^m, curr = Q_{l-1}^m
    //
    //   Q_{l-2}^m = (cosθ · c1 · Q_{l-1}^m - (l-m) · Q_l^m) / c2
    // where c1, c2 are the same coefficients as in advanceColumn.
    auto retreatColumn(int m, int l,
                       std::vector<double>& prev, std::vector<double>& curr) -> void {
        double c1 = std::sqrt(static_cast<double>((2*l-1)*(2*l+1)*(l-m))
                             / static_cast<double>(l+m));
        double c2 = (l + m - 1.0)
                  * std::sqrt(static_cast<double>((2*l+1)*(l-m)*(l-m-1))
                             / static_cast<double>((2*l-3)*(l+m)*(l+m-1)));
        for (std::size_t i = 0; i < ng_; ++i) {
            double q_lm2 = (cos_theta_[i] * c1 * prev[i] - (l - m) * curr[i]) / c2;
            curr[i] = prev[i];
            prev[i] = q_lm2;
        }
    }

    // Assemble Y_{l,m} from Q_l^m(q), write into y_lm_ cache.
    // Q already includes the normalization N_l^m, so no extra normFactor.
    //   Y_{l0}   = Q_l^0
    //   Y_{l,±m} = √2 · Q_l^m · {cos(mφ), sin(mφ)}
    auto assembleYlm(int l, int m, const std::vector<double>& Q,
                     const std::vector<double>& cm, const std::vector<double>& sm) -> void {
        int block = l * l + l;
        if (m == 0) {
            for (std::size_t i = 0; i < ng_; ++i)
                y_lm_[static_cast<std::size_t>(block)][i] = Q[i];
        } else {
            double sf = std::sqrt(2.0);
            for (std::size_t i = 0; i < ng_; ++i) {
                double val = sf * Q[i];
                y_lm_[static_cast<std::size_t>(block - m)][i] = val * sm[i];
                y_lm_[static_cast<std::size_t>(block + m)][i] = val * cm[i];
            }
        }
    }

};






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


