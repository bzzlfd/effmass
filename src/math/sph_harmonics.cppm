module;
#include <cstdio>

export module math.sph_harmonics;

import std;

export {
    class SphericalHarmonics;
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
// Precomputes all Y_lm for l <= MAX_L_TINY on construction;
// computes larger l on demand via recurrence, caching P_l^m state.
class SphericalHarmonics {
public:
    explicit SphericalHarmonics(std::span<const double> theta, std::span<const double> phi)
        : theta_(theta), phi_(phi), ng_(theta_.size())
    {
        if (theta_.empty()) {
            throw std::invalid_argument(
                "SphericalHarmonics: theta and phi must be non-empty");
        }

        // Precompute cos(theta), sin(theta)
        cos_theta_.resize(ng_);
        sin_theta_.resize(ng_);
        for (std::size_t i = 0; i < ng_; ++i) {
            cos_theta_[i] = std::cos(theta_[i]);
            sin_theta_[i] = std::sin(theta_[i]);
        }

        // Precompute cos(m*phi), sin(m*phi) for m = 0..MAX_L_TINY
        int nm = MAX_L_TINY + 1;
        cos_mphi_.resize(ng_ * static_cast<std::size_t>(nm));
        sin_mphi_.resize(ng_ * static_cast<std::size_t>(nm));

        for (std::size_t i = 0; i < ng_; ++i) {
            double c0 = std::cos(phi_[i]);
            double s0 = std::sin(phi_[i]);
            cos_mphi_[0 * ng_ + i] = 1.0;
            sin_mphi_[0 * ng_ + i] = 0.0;
            cos_mphi_[1 * ng_ + i] = c0;
            sin_mphi_[1 * ng_ + i] = s0;
            for (int m = 2; m <= MAX_L_TINY; ++m) {
                double prev_c = cos_mphi_[static_cast<std::size_t>(m - 1) * ng_ + i];
                double prev_s = sin_mphi_[static_cast<std::size_t>(m - 1) * ng_ + i];
                cos_mphi_[static_cast<std::size_t>(m) * ng_ + i] = prev_c * c0 - prev_s * s0;
                sin_mphi_[static_cast<std::size_t>(m) * ng_ + i] = prev_s * c0 + prev_c * s0;
            }
        }

        // Precompute all Y_lm for l = 0..MAX_L_TINY
        int n_blocks = nm * nm;
        cache_.resize(static_cast<std::size_t>(n_blocks) * ng_);

        for (int m = 0; m <= MAX_L_TINY; ++m) {
            // P_l^m for all l = m..MAX_L_TINY, flat: [idx * ng_ + ig]
            int n_l = MAX_L_TINY - m + 1;
            std::vector<double> p_buf(static_cast<std::size_t>(n_l) * ng_);

            // P_m^m
            if (m == 0) {
                for (std::size_t i = 0; i < ng_; ++i) p_buf[i] = 1.0;
            } else {
                double df = 1.0;
                for (int k = 1; k <= m; ++k) df *= 2.0 * k - 1.0;
                for (std::size_t i = 0; i < ng_; ++i)
                    p_buf[i] = df * std::pow(sin_theta_[i], m);
            }

            // P_{m+1}^m
            if (m < MAX_L_TINY) {
                for (std::size_t i = 0; i < ng_; ++i)
                    p_buf[ng_ + i] = (2.0 * static_cast<double>(m) + 1.0) * cos_theta_[i] * p_buf[i];
            }

            // P_l^m for l >= m+2 via three-term recurrence
            for (int l = m + 2; l <= MAX_L_TINY; ++l) {
                int idx = l - m;
                for (std::size_t i = 0; i < ng_; ++i)
                    p_buf[static_cast<std::size_t>(idx) * ng_ + i] =
                        ((2.0 * static_cast<double>(l) - 1.0) * cos_theta_[i] * p_buf[static_cast<std::size_t>(idx - 1) * ng_ + i]
                      - static_cast<double>(l + m - 1) * p_buf[static_cast<std::size_t>(idx - 2) * ng_ + i])
                      / static_cast<double>(l - m);
            }

            // Assemble Y_lm from P_l^m
            for (int l = m; l <= MAX_L_TINY; ++l) {
                int idx = l - m;
                int block = l * l + l;
                double norm = NORM_TABLE[l][m];

                if (m == 0) {
                    for (std::size_t i = 0; i < ng_; ++i)
                        cache_[static_cast<std::size_t>(block) * ng_ + i] = norm * p_buf[static_cast<std::size_t>(idx) * ng_ + i];
                } else {
                    for (std::size_t i = 0; i < ng_; ++i) {
                        double val = std::sqrt(2.0) * norm * p_buf[static_cast<std::size_t>(idx) * ng_ + i];
                        cache_[static_cast<std::size_t>(block - m) * ng_ + i] = val * sin_mphi_[static_cast<std::size_t>(m) * ng_ + i];
                        cache_[static_cast<std::size_t>(block + m) * ng_ + i] = val * cos_mphi_[static_cast<std::size_t>(m) * ng_ + i];
                    }
                }
            }
        }
    }

    // Return Y_lm for all G-vectors as a vector.
    // For l <= MAX_L_TINY: returns a copy from the precomputed cache.
    // For l >  MAX_L_TINY: computes on demand via recurrence;
    //   intermediate P vectors for the requested m_abs are kept
    //   to enable further upward recurrence.
    auto operator()(int l, int m) -> std::vector<double> {
        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("SphericalHarmonics: invalid quantum numbers (l={}, m={})", l, m));
        }

        int m_abs = std::abs(m);

        // Cached path: l <= MAX_L_TINY
        if (l <= MAX_L_TINY) {
            int block = l * l + (m + l);
            auto* src = cache_.data() + static_cast<std::ptrdiff_t>(block) * static_cast<std::ptrdiff_t>(ng_);
            return std::vector<double>(src, src + static_cast<std::ptrdiff_t>(ng_));
        }

        // On-demand path: l > MAX_L_TINY
        std::fprintf(stderr,
            "SphericalHarmonics: warning - l=%d exceeds MAX_L_TINY (%d), computing on demand\n",
            l, MAX_L_TINY);

        // Advance Legendre recurrence to target l
        auto plm = advanceP(l, m_abs);

        // Normalization factor (self-contained)
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
                result[i] = norm * plm[i];
        } else {
            double sf = std::sqrt(2.0) * norm;
            if (m > 0) {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * plm[i] * std::cos(static_cast<double>(m) * phi_[i]);
            } else {
                for (std::size_t i = 0; i < ng_; ++i)
                    result[i] = sf * plm[i] * std::sin(static_cast<double>(m_abs) * phi_[i]);
            }
        }

        return result;
    }

private:
    std::span<const double> theta_;
    std::span<const double> phi_;
    std::size_t ng_;
    std::vector<double> cos_theta_;
    std::vector<double> sin_theta_;
    std::vector<double> cos_mphi_;
    std::vector<double> sin_mphi_;
    std::vector<double> cache_;

    // Per-m_abs recurrence state for associated Legendre polynomials.
    // p_curr = P_{l_state}^{m_abs},  p_prev = P_{l_state-1}^{m_abs}
    struct RecurState {
        std::vector<double> p_prev;
        std::vector<double> p_curr;
        int l_state = -1;
    };
    std::unordered_map<int, RecurState> state_;

    // Advance Legendre recurrence for m_abs up to target_l.
    auto advanceP(int target_l, int m_abs) -> std::span<const double> {
        auto& s = state_[m_abs];

        if (s.l_state < 0) {
            // Seed from scratch
            s.p_prev.resize(ng_);
            s.p_curr.resize(ng_);

            if (m_abs == 0) {
                for (std::size_t i = 0; i < ng_; ++i) s.p_curr[i] = 1.0;
            } else {
                double df = 1.0;
                for (int k = 1; k <= m_abs; ++k) df *= 2.0 * k - 1.0;
                for (std::size_t i = 0; i < ng_; ++i)
                    s.p_curr[i] = df * std::pow(sin_theta_[i], m_abs);
            }
            s.l_state = m_abs;

            if (target_l <= m_abs) return s.p_curr;
        }

        if (target_l > s.l_state) {
            // Special first step: P_{m+1}^m from P_m^m
            if (s.l_state == m_abs) {
                for (std::size_t i = 0; i < ng_; ++i) {
                    double next = (2.0 * static_cast<double>(m_abs) + 1.0) * cos_theta_[i] * s.p_curr[i];
                    s.p_prev[i] = s.p_curr[i];
                    s.p_curr[i] = next;
                }
                ++s.l_state;
            }

            // General three-term recurrence
            while (s.l_state < target_l) {
                int n = s.l_state + 1;
                for (std::size_t i = 0; i < ng_; ++i) {
                    double next = ((2.0 * static_cast<double>(n) - 1.0) * cos_theta_[i] * s.p_curr[i]
                                 - static_cast<double>(n + m_abs - 1) * s.p_prev[i])
                                 / static_cast<double>(n - m_abs);
                    s.p_prev[i] = s.p_curr[i];
                    s.p_curr[i] = next;
                }
                ++s.l_state;
            }
        } else if (target_l < s.l_state) {
            // Requested l below current state: recompute from scratch
            state_.erase(m_abs);
            return advanceP(target_l, m_abs);
        }

        return s.p_curr;
    }
};
