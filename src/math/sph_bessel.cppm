export module math.sph_bessel;

import std;

export {
    class SphericalBesselJ;
}


constexpr double BESSEL_EPS = 1e-15;
constexpr double SERIES_THRESH = 0.5;
constexpr double SERIES_TOL = 1e-16;
constexpr int SERIES_MAX_ITER = 50;

auto sphericalBesselJImpl(int l, double x) -> double {
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


namespace archived {

auto sphericalBesselJ(int l, double x) -> double {
    bool neg = (x < 0.0) && (l % 2 == 1);
    double res = sphericalBesselJImpl(l, std::abs(x));
    return neg ? -res : res;
}

auto sphericalBesselJ(int l, std::span<const double> xs, std::span<double> out) -> void {
    int n = static_cast<int>(xs.size());
    bool flip_sign = (l % 2 == 1);

    if (l == 0) {
        for (int i = 0; i < n; ++i) {
            double x = std::abs(xs[i]);
            if (x < BESSEL_EPS) {
                out[i] = 1.0;
            } else if (x < SERIES_THRESH) {
                double x2_half = x * x / 2.0;
                double term = 1.0;
                double sum = 1.0;
                for (int k = 1; k <= SERIES_MAX_ITER; ++k) {
                    term *= -x2_half / (static_cast<double>(k) * (2.0 * k + 1.0));
                    sum += term;
                    if (std::abs(term) < SERIES_TOL * std::abs(sum)) break;
                }
                out[i] = sum;
            } else {
                out[i] = std::sin(x) / x;
            }
        }
        return;
    }

    if (l == 1) {
        for (int i = 0; i < n; ++i) {
            double x = xs[i];
            bool neg = (x < 0.0);
            x = std::abs(x);
            if (x < BESSEL_EPS) { out[i] = 0.0; continue; }
            if (x < SERIES_THRESH) {
                double x2_half = x * x / 2.0;
                double term = 1.0;
                double sum = 1.0;
                for (int k = 1; k <= SERIES_MAX_ITER; ++k) {
                    term *= -x2_half / (static_cast<double>(k) * (2.0 * k + 3.0));
                    sum += term;
                    if (std::abs(term) < SERIES_TOL * std::abs(sum)) break;
                }
                double res = sum * x / 3.0;
                out[i] = neg ? -res : res;
                continue;
            }
            double s = std::sin(x);
            double c = std::cos(x);
            double res = s / (x * x) - c / x;
            out[i] = neg ? -res : res;
        }
        return;
    }

    double dfact = 1.0;
    for (int j = 3; j <= 2 * l + 1; j += 2) dfact *= j;

    for (int i = 0; i < n; ++i) {
        double x = xs[i];
        bool neg = (x < 0.0) && flip_sign;
        x = std::abs(x);

        if (x < BESSEL_EPS) {
            out[i] = 0.0;
            continue;
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
            double res = sum * std::pow(x, l) / dfact;
            out[i] = neg ? -res : res;
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
            out[i] = neg ? -j_curr : j_curr;
        }
    }
}

} // namespace archived


// Forward-only spherical Bessel iterator.
//   j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x)
//
//   SphericalBesselJ(r_span,q) -- vector: j_l(q*r_i)  (q >= 0)
//
// value() returns a read-only span over the grid.
class SphericalBesselJ {
public:
    SphericalBesselJ(std::span<const double> r, double q)
        : n_(r.size()), l_(0)
    {
        r_cache_.assign(r.begin(), r.end());
        x_.resize(n_);
        j_curr_.resize(n_);
        j_next_.resize(n_);
        for (std::size_t i = 0; i < n_; ++i)
            x_[i] = q * r_cache_[i];
        initSeeds();
    }

    auto reset(double q) -> SphericalBesselJ& {
        l_ = 0;
        for (std::size_t i = 0; i < n_; ++i)
            x_[i] = q * r_cache_[i];
        initSeeds();
        return *this;
    }

    auto advance(int n = 1) -> SphericalBesselJ& {
        if (n < 0) {
            throw std::invalid_argument(
                "SphericalBesselJ::advance: negative step not supported");
        }
        if (n == 0) return *this;

        for (int step = 0; step < n; ++step) {
            for (std::size_t i = 0; i < n_; ++i) {
                if (std::abs(x_[i]) < BESSEL_EPS) {
                    // j_{l+1}(0) = 0 for all l >= 0
                    j_curr_[i] = j_next_[i];
                    j_next_[i] = 0.0;
                } else if (std::abs(x_[i]) >= SERIES_THRESH) {
                    double j_new = (2.0 * l_ + 3.0) / x_[i] * j_next_[i] - j_curr_[i];
                    j_curr_[i] = j_next_[i];
                    j_next_[i] = j_new;
                } else {
                    j_curr_[i] = j_next_[i];
                    j_next_[i] = sphericalBesselJImpl(l_ + 2, std::abs(x_[i]));
                }
            }
            ++l_;
        }
        return *this;
    }

    // Returns a "read-only" "span" of j_l(x) for all x in the grid.
    auto value() const -> std::span<const double> {
        return j_curr_;
    }

    auto currentL() const -> int { return l_; }

    auto operator++() -> double {
        advance(1);
        return value()[0];
    }

    auto operator++(int) -> double {
        double old = value()[0];
        advance(1);
        return old;
    }

private:
    std::size_t n_;
    int l_;
    std::vector<double> r_cache_;
    std::vector<double> x_;
    std::vector<double> j_curr_;
    std::vector<double> j_next_;

    auto initSeeds() -> void {
        for (std::size_t i = 0; i < n_; ++i) {
            double ax = std::abs(x_[i]);
            if (ax < BESSEL_EPS) {
                j_curr_[i] = 1.0;
                j_next_[i] = 0.0;
            } else if (ax < SERIES_THRESH) {
                double x2_half = x_[i] * x_[i] / 2.0;
                double term = 1.0;
                double sum0 = 1.0;
                for (int k = 1; k <= SERIES_MAX_ITER; ++k) {
                    term *= -x2_half / (static_cast<double>(k) * (2.0 * k + 1.0));
                    sum0 += term;
                    if (std::abs(term) < SERIES_TOL * std::abs(sum0)) break;
                }
                j_curr_[i] = sum0;
                term = 1.0;
                double sum1 = 1.0;
                for (int k = 1; k <= SERIES_MAX_ITER; ++k) {
                    term *= -x2_half / (static_cast<double>(k) * (2.0 * k + 3.0));
                    sum1 += term;
                    if (std::abs(term) < SERIES_TOL * std::abs(sum1)) break;
                }
                j_next_[i] = sum1 * x_[i] / 3.0;
            } else {
                j_curr_[i] = std::sin(x_[i]) / x_[i];
                j_next_[i] = j_curr_[i] / x_[i] - std::cos(x_[i]) / x_[i];
            }
        }
    }
};
