module;

export module math.sph_harmonics;

import std;

export {
    class RealSphericalHarmonicsEngine;
    class RealSphericalHarmonicsData;
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


// =============================================================================
//  QRecurrence — normalized Legendre (Q_l^m) column-wise recurrence
//
//  Three-step protocol shared by Q_l^m and its derivatives:
//    seed    →  U_m^m          (top of the column)
//    step1   →  U_{m+1}^m      (first step away from the diagonal)
//    advance →  U_l^m  (l≥m+2) (upward recurrence within a column)
//
//  The struct owns no storage; it reads trig data through spans and
//  mutates externally-owned buffers (Q_lm_ or local temporaries).
// =============================================================================
struct QRecurrence {
    std::span<const double> cos_theta;
    std::span<const double> sin_theta;

    /// Seed Q_m^m = (-1)^m / (2√π) · ∏_{k=1}^m √((2k+1)/(2k)) · sin^m(θ)
    auto seed(int m, std::vector<double>& curr) const -> void {
        auto ng = cos_theta.size();
        double Q0 = 0.5 / std::sqrt(std::numbers::pi);
        if (m == 0) {
            for (std::size_t i = 0; i < ng; ++i) curr[i] = Q0;
        } else {
            double coeff = Q0;
            if (m % 2 != 0) coeff = -coeff;  // Condon-Shortley (-1)^m
            for (int k = 1; k <= m; ++k)
                coeff *= std::sqrt(static_cast<double>(2 * k + 1) / (2 * k));
            for (std::size_t i = 0; i < ng; ++i)
                curr[i] = coeff * std::pow(sin_theta[i], m);
        }
    }

    /// One step: Q_{m+1}^m = √(2m+3) · cosθ · Q_m^m
    auto step1(int m, std::vector<double>& prev,
               std::vector<double>& curr) const -> void {
        auto ng = cos_theta.size();
        double factor = std::sqrt(2.0 * m + 3.0);
        for (std::size_t i = 0; i < ng; ++i) {
            prev[i] = curr[i];
            curr[i] = factor * cos_theta[i] * curr[i];
        }
    }

    /// Forward three-term recurrence for Q_l^m (normalized).
    /// On entry: prev = Q_{l-2}^m, curr = Q_{l-1}^m
    /// On exit:  prev = Q_{l-1}^m, curr = Q_l^m
    auto advance(int m, int l, std::vector<double>& prev,
                 std::vector<double>& curr) const -> void {
        auto ng = cos_theta.size();
        double r1 = (2.0 * l - 1.0)
                  * std::sqrt(((2.0 * l + 1.0) * (l - m))
                             / ((2.0 * l - 1.0) * (l + m)));
        double r2 = (l + m - 1.0)
                  * std::sqrt(((2.0 * l + 1.0) * (l - m) * (l - m - 1))
                             / ((2.0 * l - 3.0) * (l + m) * (l + m - 1)));
        double denom = 1.0 / (l - m);
        for (std::size_t i = 0; i < ng; ++i) {
            double next = (cos_theta[i] * r1 * curr[i] - r2 * prev[i]) * denom;
            prev[i] = curr[i];
            curr[i] = next;
        }
    }
};


// Vectorized spherical harmonics evaluator over K-point G-vectors.
// Q_l^m (normalised Legendre) is the canonical permanent form;
// Y_lm is a derived quantity assembled from Q + cos/sin tables.
// Engine owns trig, Q, cos/sin tables, and Legendre recurrence.
// Cache (separate Data class) owns the lazy-assembled y_lm_ cache.
class RealSphericalHarmonicsEngine {
public:
    enum class CacheMode { Full, None };

    explicit RealSphericalHarmonicsEngine(std::span<const double> theta, std::span<const double> phi,
                                     int l_max_resident, CacheMode mode = CacheMode::Full);

    /// Pre-allocate internal buffers to at least max_ng capacity.
    /// Q_lm_, cos_mphi_ and sin_mphi_ are reserved in Full mode;
    /// trig arrays are always reserved.
    auto reserveNg(int max_ng) -> void;

    /// Switch l_max_resident.  Preserves Q_lm_ and trig tables for
    /// zero-cost re-expansion.  Expand continues the Legendre
    /// recurrence from the existing Q_lm_ state and fills Q_lm_,
    /// cos_mphi_ and sin_mphi_ for the new (l,m) range.
    auto setLMax(int l_max_new) -> void;

    /// Re-initialise with new theta/phi spans.  Preserves CacheMode.
    /// Trig arrays are recomputed from scratch; Q_lm_ and cos/sin tables
    /// are NOT rebuilt here — call setLMax() separately if needed.
    /// Existing vector allocations are reused when ng_ <= capacity.
    auto reinit(std::span<const double> theta, std::span<const double> phi) -> void;

    /// Return Y_lm for all G-vectors, assembled from the permanent
    /// Q_l^m + cos(mφ)/sin(mφ) tables.  Returns by value (no cache).
    /// Requires CacheMode::Full and l <= l_max_resident.
    auto get(int l, int m) const -> std::vector<double>;

    /// Compute Y_lm on the fly. Works in any mode, for any (l,m).
    /// Returns by value (no persistent cache).
    auto compute(int l, int m) const -> std::vector<double>;

    /// Convenience dispatcher.  In Full mode: delegates to get()
    /// (copy) when l ≤ l_max_resident_, otherwise compute().
    /// In None mode: always compute().
    auto operator()(int l, int m) const -> std::vector<double> {
        if (mode_ == CacheMode::None) return compute(l, m);
        if (l > l_max_resident_) return compute(l, m);
        return get(l, m);
    }

    // -- accessors for Data --
    auto revision() const -> std::size_t { return revision_; }
    auto lMaxResident() const -> int { return l_max_resident_; }
    auto ng() const -> std::size_t { return ng_; }
    auto ngCapacity() const -> std::size_t { return cos_theta_.capacity(); }
    auto mode() const -> CacheMode { return mode_; }

private:
    std::span<const double> theta_;
    std::span<const double> phi_;
    std::size_t ng_;
    int l_max_resident_ = -1;
    CacheMode mode_ = CacheMode::Full;
    std::size_t revision_ = 1;
    std::vector<double> cos_theta_;
    std::vector<double> sin_theta_;
    std::vector<double> cos_phi_;
    std::vector<double> sin_phi_;
    // === Q recurrence engine ===
    QRecurrence q_recur_{};

    // === Q_l^m permanent storage ===
    // Q_lm_ is the canonical Legendre storage (includes normalization
    // and Condon-Shortley phase).
    //   Q_lm_ indexing: flat = l*(l+1)/2 + m_abs  (compact lower-triangular).
    std::vector<std::vector<double>> Q_lm_;

    // cos_mphi_[m] = cos(m·φ), sin_mphi_[m] = sin(m·φ), m in [0, l_max_resident_].
    // Precomputed during setLMax() for the full m range.
    std::vector<std::vector<double>> cos_mphi_;
    std::vector<std::vector<double>> sin_mphi_;

};


// =============================================================================
//  Construction
// =============================================================================

RealSphericalHarmonicsEngine::RealSphericalHarmonicsEngine(
    std::span<const double> theta, std::span<const double> phi,
    int l_max_resident, CacheMode mode)
    : theta_(theta), phi_(phi), ng_(theta_.size()), mode_(mode)
{
    if (theta_.empty()) {
        throw std::invalid_argument(
            "RealSphericalHarmonicsEngine: theta and phi must be non-empty");
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

    // Bind recurrence engine to trig arrays.
    q_recur_ = QRecurrence{std::span<const double>(cos_theta_.data(), cos_theta_.size()),
                           std::span<const double>(sin_theta_.data(), sin_theta_.size())};

    if (mode_ == CacheMode::Full) setLMax(l_max_resident);
}


// =============================================================================
//  Public API — reserveNg / reinit
// =============================================================================

auto RealSphericalHarmonicsEngine::reinit(
    std::span<const double> theta, std::span<const double> phi) -> void
{
    // Guard against silent reallocation.
    auto new_ng = theta.size();
    if (new_ng > cos_theta_.capacity()) {
        throw std::runtime_error(
            std::format("RealSphericalHarmonicsEngine::reinit(): new ng_={} exceeds "
                        "pre-allocated capacity={}. Call reserveNg({}) before reinit().",
                        new_ng, cos_theta_.capacity(), new_ng));
    }

    theta_ = theta;
    phi_ = phi;
    ng_ = new_ng;

    // Recompute trig arrays — resize to ng_, keeping existing capacity.
    cos_theta_.resize(ng_);
    sin_theta_.resize(ng_);
    for (std::size_t i = 0; i < ng_; ++i) {
        cos_theta_[i] = std::cos(theta[i]);
        sin_theta_[i] = std::sin(theta[i]);
    }
    cos_phi_.resize(ng_);
    sin_phi_.resize(ng_);
    for (std::size_t i = 0; i < ng_; ++i) {
        cos_phi_[i] = std::cos(phi[i]);
        sin_phi_[i] = std::sin(phi[i]);
    }

    // Rebind recurrence engine to trig arrays (data ptr may have moved).
    q_recur_ = QRecurrence{std::span<const double>(cos_theta_.data(), cos_theta_.size()),
                           std::span<const double>(sin_theta_.data(), sin_theta_.size())};

    // Invalidate Q and cos/sin tables — they are no longer consistent
    // with the new trig data.
    l_max_resident_ = -1;
    ++revision_;
}

auto RealSphericalHarmonicsEngine::reserveNg(int max_ng) -> void {
    auto n = static_cast<std::size_t>(max_ng);
    cos_theta_.reserve(n);
    sin_theta_.reserve(n);
    cos_phi_.reserve(n);
    sin_phi_.reserve(n);
    if (mode_ == CacheMode::Full) {
        for (auto& v : Q_lm_) v.reserve(n);
        for (auto& v : cos_mphi_) v.reserve(n);
        for (auto& v : sin_mphi_) v.reserve(n);
    }
}


// =============================================================================
//  Public API  —  setLMax
// =============================================================================

auto RealSphericalHarmonicsEngine::setLMax(int l_max_new) -> void {
    if (l_max_new < 0) {
        throw std::invalid_argument(
            "RealSphericalHarmonicsEngine::setLMax: l_max_resident must be non-negative");
    }
    if (mode_ == CacheMode::None) return;
    if (l_max_new == l_max_resident_) return;

    // -- Shrink: keep Q_lm_, cos_mphi_ and sin_mphi_ for zero-cost re-expansion. --
    if (l_max_new < l_max_resident_) {
        l_max_resident_ = l_max_new;
        ++revision_;
        return;
    }

    // -- Expand --
    int old_nm = l_max_resident_ + 1;
    int new_nm = l_max_new + 1;
    auto ng_capacity = cos_theta_.capacity();

    // Extend Q_lm_ for new l,m pairs; existing entries keep their capacity.
    {
        std::size_t new_total = static_cast<std::size_t>(new_nm) * static_cast<std::size_t>(new_nm + 1) / 2;
        auto old_q = Q_lm_.size();
        Q_lm_.resize(new_total);
        for (std::size_t i = old_q; i < new_total; ++i) {
            Q_lm_[i].reserve(ng_capacity);
            Q_lm_[i].resize(ng_);
        }
    }

    // Extend cos_mphi_ / sin_mphi_ tables for new m values.
    {
        auto old_cm = cos_mphi_.size();
        cos_mphi_.resize(static_cast<std::size_t>(new_nm));
        sin_mphi_.resize(static_cast<std::size_t>(new_nm));
        for (std::size_t m = old_cm; m < static_cast<std::size_t>(new_nm); ++m) {
            cos_mphi_[m].reserve(ng_capacity);
            cos_mphi_[m].resize(ng_);
            sin_mphi_[m].reserve(ng_capacity);
            sin_mphi_[m].resize(ng_);
        }
    }

    // Helper: compact triangular index for Q_lm_ (m_abs in [0, l]).
    //   Q_lm_[l*(l+1)/2 + m_abs] stores Q_l^{m_abs}.
    auto qIdx = [](int l, int m_abs) -> std::size_t {
        return static_cast<std::size_t>(l * (l + 1) / 2 + m_abs);
    };

    std::vector<double> cm(ng_, 1.0), sm(ng_, 0.0);

    for (int m = 0; m <= l_max_new; ++m) {
        if (m > l_max_resident_) {
            // -- New m column: seed from scratch into local variables,
            //    save each value to Q_lm_ as we go. --
            std::vector<double> prev(ng_), curr(ng_);
            q_recur_.seed(m, curr);
            Q_lm_[qIdx(m, m)] = curr;   // save Q_m^m (copy)

            if (m < l_max_new) {
                q_recur_.step1(m, prev, curr);   // prev = Q_m^m, curr = Q_{m+1}^m
                Q_lm_[qIdx(m + 1, m)] = curr;
                for (int l = m + 2; l <= l_max_new; ++l) {
                    q_recur_.advance(m, l, prev, curr);
                    Q_lm_[qIdx(l, m)] = curr;  // save Q_l^m (copy)
                }
            }

        } else if (l_max_resident_ == m) {
            // -- Existing m at P_m^m only (old l_max == m) --
            if (m < l_max_new) {
                std::vector<double> prev(ng_);
                std::vector<double> curr = Q_lm_[qIdx(m, m)];  // copy from storage
                q_recur_.step1(m, prev, curr);
                Q_lm_[qIdx(m + 1, m)] = curr;
                for (int l = m + 2; l <= l_max_new; ++l) {
                    q_recur_.advance(m, l, prev, curr);
                    Q_lm_[qIdx(l, m)] = curr;  // save Q_l^m (copy)
                }
            }

        } else {
            // -- Existing m with full Q_{old-1}^m and Q_{old}^m in Q_lm_ storage --
            int old_l = l_max_resident_;
            std::vector<double> prev = Q_lm_[qIdx(old_l - 1, m)];  // copy from storage
            std::vector<double> curr = Q_lm_[qIdx(old_l, m)];      // copy from storage
            for (int l = old_l + 1; l <= l_max_new; ++l) {
                q_recur_.advance(m, l, prev, curr);
                Q_lm_[qIdx(l, m)] = curr;  // save Q_l^m (copy)
            }
        }

        // Save cos(mφ)/sin(mφ) for this m.
        cos_mphi_[static_cast<std::size_t>(m)] = cm;
        sin_mphi_[static_cast<std::size_t>(m)] = sm;

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
    ++revision_;
}

auto RealSphericalHarmonicsEngine::get(int l, int m) const -> std::vector<double> {
    if (l < 0 || std::abs(m) > l) {
        throw std::invalid_argument(
            std::format("RealSphericalHarmonicsEngine::get(): invalid quantum numbers (l={}, m={})", l, m));
    }
    if (mode_ == CacheMode::None) {
        throw std::runtime_error(
            "RealSphericalHarmonicsEngine::get(): not available in CacheMode::None, use compute()");
    }
    if (l > l_max_resident_) {
        throw std::runtime_error(
            std::format("RealSphericalHarmonicsEngine::get(): l={} > l_max_resident_={}, call setLMax({}) first, or use compute()", l, l_max_resident_, l));
    }

    int m_abs = std::abs(m);
    std::size_t q_idx = static_cast<std::size_t>(l * (l + 1) / 2 + m_abs);
    const auto& Q = Q_lm_[q_idx];

    if (m == 0) {
        return Q;  // Y_{l,0} = Q_l^0, no cos/sin factor
    }

    double sf = std::sqrt(2.0);
    const auto& cm = cos_mphi_[static_cast<std::size_t>(m_abs)];
    const auto& sm = sin_mphi_[static_cast<std::size_t>(m_abs)];

    std::vector<double> result(ng_);
    bool positive = (m > 0);
    for (std::size_t i = 0; i < ng_; ++i) {
        double val = sf * Q[i];
        result[i] = positive ? (val * cm[i]) : (val * sm[i]);
    }
    return result;
}

auto RealSphericalHarmonicsEngine::compute(int l, int m) const -> std::vector<double> {
    if (l < 0 || std::abs(m) > l) {
        throw std::invalid_argument(
            std::format("RealSphericalHarmonicsEngine::compute: invalid quantum numbers (l={}, m={})", l, m));
    }

    int m_abs = std::abs(m);

    // Q_l^m = N_l^m * P_l^m already includes normalization
    std::vector<double> q_curr(ng_);
    q_recur_.seed(m_abs, q_curr);

    if (l > m_abs) {
        std::vector<double> q_prev(ng_);
        q_recur_.step1(m_abs, q_prev, q_curr);
        for (int n = m_abs + 2; n <= l; ++n)
            q_recur_.advance(m_abs, n, q_prev, q_curr);
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



// =============================================================================
//  RealSphericalHarmonicsData  —  lazy y_lm cache over an Engine
// =============================================================================

/// Lazy cache of assembled Y_lm values for a given RealSphericalHarmonicsEngine.
/// Does not own the Engine — pass a const Engine& to get() or operator().
/// The revision-based syncToEngine() detects Engine changes
/// (reinit / setLMax) and resizes + invalidates the cache automatically.
class RealSphericalHarmonicsData {
public:
    RealSphericalHarmonicsData() = default;

    /// Pre-allocate capacity for y_lm_ entries.  Safe to call before
    /// any get() — the hint is stored and applied during syncToEngine.
    auto reserveNg(std::size_t max_ng) -> void {
        reserved_ng_ = max_ng;
        for (auto& v : y_lm_) v.reserve(max_ng);
    }

    /// Return Y_lm for all G-vectors from the lazy-assembled cache.
    /// On first call for each (l,m), delegates to Engine::get() (when
    /// available: Full mode + l ≤ lMaxResident) or Engine::compute()
    /// otherwise, and moves the result into the internal cache.  Subsequent
    /// calls return the cached reference at O(1) cost.
    /// Works in all CacheModes and for any valid (l,m) pair.
    auto get(const RealSphericalHarmonicsEngine& e, int l, int m) -> const std::vector<double>& {
        syncToEngine(e);

        if (l < 0 || std::abs(m) > l) {
            throw std::invalid_argument(
                std::format("RealSphericalHarmonicsData::get(): invalid quantum numbers (l={}, m={})", l, m));
        }

        std::size_t ylm_idx = static_cast<std::size_t>(l * l + (m + l));

        // Grow the outer cache array if this (l,m) ylm_idx is beyond the
        // range tracked by syncToEngine (e.g. l > lMaxResident).
        if (ylm_idx >= y_lm_.size()) {
            auto new_size = ylm_idx + 1;
            auto old_size = y_lm_.size();
            y_lm_.resize(new_size);
            y_lm_valid_.resize(new_size, false);
            for (auto i = old_size; i < new_size; ++i) {
                y_lm_[i].reserve(std::max(reserved_ng_, e.ngCapacity()));
                y_lm_[i].resize(e.ng());
            }
        }

        // Lazy assembly: delegate to Engine and move result into cache
        if (!y_lm_valid_[ylm_idx]) {
            if (e.mode() == RealSphericalHarmonicsEngine::CacheMode::Full && l <= e.lMaxResident()) {
                y_lm_[ylm_idx] = e.get(l, m);
            } else {
                y_lm_[ylm_idx] = e.compute(l, m);
            }
            y_lm_valid_[ylm_idx] = true;
        }

        return y_lm_[ylm_idx];
    }

private:
    mutable std::vector<std::vector<double>> y_lm_;
    mutable std::vector<bool> y_lm_valid_;
    mutable std::size_t engine_revision_ = 0;
    mutable std::size_t ng_cached_ = 0;
    std::size_t reserved_ng_ = 0;

    /// Invalidate the cache and resize existing entry vectors when the
    /// Engine has changed (reinit/setLMax).  O(1) hot-path when nothing
    /// changed.  The outer cache array is grown on-demand in get().
    void syncToEngine(const RealSphericalHarmonicsEngine& e) const {
        if (engine_revision_ == e.revision()) return;

        // Invalidate all existing cached values.
        std::fill(y_lm_valid_.begin(), y_lm_valid_.end(), false);

        // Resize existing entry vectors if ng changed (e.g. after reinit).
        if (ng_cached_ != e.ng()) {
            auto ng = e.ng();
            for (auto& v : y_lm_) {
                v.reserve(std::max(reserved_ng_, e.ngCapacity()));
                v.resize(ng);
            }
            ng_cached_ = e.ng();
        }

        engine_revision_ = e.revision();
    }
};


// =============================================================================
//  Archived scalar implementations
// =============================================================================

namespace archived {

auto realSphericalHarmonic(int l, int m, double theta, double phi) -> double {
    if (l < 0 || std::abs(m) > l) return 0.0;

    int m_abs = std::abs(m);

    double Plm = legendreP(l, m_abs, theta);
    double norm = normFactor(l, m_abs);

    // Condon-Shortley phase (-1)^{|m|}
    double cs = (m_abs % 2 == 0) ? 1.0 : -1.0;

    if (m == 0) {
        return norm * Plm;
    } else if (m > 0) {
        return cs * std::sqrt(2.0) * norm * Plm * std::cos(m * phi);
    } else {
        return cs * std::sqrt(2.0) * norm * Plm * std::sin(m_abs * phi);
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

        // Condon-Shortley phase (-1)^m
        double cs = (m % 2 == 0) ? 1.0 : -1.0;

        // Assemble real spherical harmonics for this m
        for (int l = m; l <= l_max; ++l) {
            int idx = l - m;
            double norm = norms[l][m];
            int base = l * l + l;

            if (m == 0) {
                out[base] = norm * p[idx];
            } else {
                double val = cs * std::sqrt(2.0) * norm * p[idx];
                out[base - m] = val * sin_mphi[m];
                out[base + m] = val * cos_mphi[m];
            }
        }
    }
}


} // namespace archived (realSphericalHarmonic, realSphericalHarmonics)
