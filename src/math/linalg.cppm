export module math.linalg;

import std;

export {
    struct DiagonalizeResult {
        std::vector<double> eigenvalues;
        std::vector<double> eigenvectors;   // column-major, column j ↔ eigenvalues[j]
        bool converged = false;
    };

    auto diagonalize_jacobi(std::span<const double> matrix, int n) -> DiagonalizeResult;
}


namespace {

constexpr double eps = 1e-15;
constexpr int max_sweeps = 50;

auto jacobi_rotate(std::span<double> a, std::span<double> v, int n) -> bool {
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double max_off = 0.0;
        for (int i = 0; i < n - 1; ++i)
            for (int j = i + 1; j < n; ++j)
                max_off = std::max(max_off, std::abs(a[static_cast<std::size_t>(i) * n + j]));

        if (max_off < eps) return true;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                auto& apq = a[static_cast<std::size_t>(p) * n + q];
                if (std::abs(apq) < eps) continue;

                auto theta = (a[static_cast<std::size_t>(q) * n + q]
                            - a[static_cast<std::size_t>(p) * n + p])
                            / (2.0 * apq);
                auto t = std::copysign(1.0, theta)
                       / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
                auto c = 1.0 / std::sqrt(1.0 + t * t);
                auto s = t * c;

                auto app = a[static_cast<std::size_t>(p) * n + p];
                auto aqq = a[static_cast<std::size_t>(q) * n + q];
                a[static_cast<std::size_t>(p) * n + p] = app - t * apq;
                a[static_cast<std::size_t>(q) * n + q] = aqq + t * apq;
                apq = 0.0;
                a[static_cast<std::size_t>(q) * n + p] = 0.0;

                auto n_sz = static_cast<std::size_t>(n);
                auto pp = static_cast<std::size_t>(p);
                auto qq = static_cast<std::size_t>(q);
                for (std::size_t i = 0; i < n_sz; ++i) {
                    if (i == pp || i == qq) continue;
                    auto api = a[i * n_sz + pp];
                    auto aqi = a[i * n_sz + qq];
                    auto n_api = c * api - s * aqi;
                    auto n_aqi = s * api + c * aqi;
                    a[i * n_sz + pp] = n_api;
                    a[i * n_sz + qq] = n_aqi;
                    a[pp * n_sz + i] = n_api;
                    a[qq * n_sz + i] = n_aqi;
                }

                for (std::size_t i = 0; i < n_sz; ++i) {
                    auto vip = v[i * n_sz + pp];
                    auto viq = v[i * n_sz + qq];
                    v[i * n_sz + pp] = c * vip - s * viq;
                    v[i * n_sz + qq] = s * vip + c * viq;
                }
            }
        }
    }
    return false;
}

} // anonymous namespace


auto diagonalize_jacobi(std::span<const double> matrix, int n) -> DiagonalizeResult {
    auto n_sz = static_cast<std::size_t>(n);
    std::vector<double> a(matrix.begin(), matrix.end());
    std::vector<double> v(n_sz * n_sz, 0.0);
    for (std::size_t i = 0; i < n_sz; ++i) v[i * n_sz + i] = 1.0;

    bool conv = jacobi_rotate(a, v, n);

    // Sort eigenvalues ascending and permute eigenvectors
    std::vector<std::size_t> idx(n_sz);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::ranges::sort(idx, [&](auto i, auto j) {
        return a[i * n_sz + i] < a[j * n_sz + j];
    });

    std::vector<double> eigenvalues(n_sz);
    std::vector<double> sorted_v(n_sz * n_sz);
    for (std::size_t j = 0; j < n_sz; ++j) {
        eigenvalues[j] = a[idx[j] * n_sz + idx[j]];
        for (std::size_t i = 0; i < n_sz; ++i) {
            sorted_v[i * n_sz + j] = v[i * n_sz + idx[j]];
        }
    }

    return {std::move(eigenvalues), std::move(sorted_v), conv};
}
