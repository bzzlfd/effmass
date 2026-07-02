import std;
import math;


auto check(bool cond, std::string_view msg) -> void {
    if (!cond) {
        throw std::runtime_error(std::string("FAILED: ") + std::string(msg));
    }
    std::println("  PASSED: {}", msg);
}

auto near(double a, double b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}


auto test_spherical_harmonics() -> void {
    std::println("\n=== Real Spherical Harmonics ===");

    double eps = 1e-12;

    // Y_{00} = 1/(2*sqrt(pi))
    {
        double theta = 0.0, phi = 0.0;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 3};
        RealSphericalHarmonicsData data;
        const auto& y00 = data.get(eng, 0, 0);
        check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps), "Y_00 constant");
    }

    // Y_{10}(theta) = sqrt(3/(4pi)) * cos(theta)
    {
        double theta0 = 0.5, phi0 = 0.0;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        RealSphericalHarmonicsData data;
        const auto& y10 = data.get(eng, 1, 0);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta0);
        check(near(y10[0], expected, eps), "Y_10(theta)");
    }

    // Y_{1,1} = -sqrt(3/(4pi)) * sin(theta) * cos(phi)
    // (Condon-Shortley factor (-1)^1 = -1)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        RealSphericalHarmonicsData data;
        const auto& y11 = data.get(eng, 1, 1);
        double expected = -std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::cos(phi0);
        check(near(y11[0], expected, eps), "Y_11(theta,phi) (CS)");
    }

    // Y_{1,-1} = -sqrt(3/(4pi)) * sin(theta) * sin(phi)
    // (Condon-Shortley factor (-1)^1 = -1)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        RealSphericalHarmonicsData data;
        const auto& y1m1 = data.get(eng, 1, -1);
        double expected = -std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::sin(phi0);
        check(near(y1m1[0], expected, eps), "Y_1,-1(theta,phi) (CS)");
    }

    // |m| > l should throw
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        RealSphericalHarmonicsData data;
        try {
            data.get(eng, 2, 3);
            check(false, "Y_{2,3} should have thrown");
        } catch (const std::invalid_argument&) {
            check(true, "Y_{2,3} throws (|m| > l)");
        }
    }
}

auto test_spherical_harmonics_batch() -> void {
    std::println("\n=== Real Spherical Harmonics (batch) ===");

    double eps = 1e-12;
    double theta = 0.7;
    double phi = 0.3;
    int l_max = 4;

    RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};
    RealSphericalHarmonicsData data;

    for (int l = 0; l <= l_max; ++l) {
        for (int m = -l; m <= l; ++m) {
            const auto& y_lm = data.get(eng, l, m);
            check(y_lm.size() == 1, std::format("batch Y_{{{},{}}} size", l, m));
        }
    }
}

auto test_compute() -> void {
    std::println("\n=== compute() ===");

    double eps = 1e-12;
    double theta = 0.5;
    double phi = 0.3;

    // compute() matches get() for cached range
    {
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 3};
        RealSphericalHarmonicsData data;
        for (int l = 0; l <= 3; ++l) {
            for (int m = -l; m <= l; ++m) {
                const auto& cached = data.get(eng, l, m);
                auto computed = eng.compute(l, m);
                check(near(cached[0], computed[0], eps),
                      std::format("compute({},{}) matches get", l, m));
            }
        }
    }

    // compute() works for l > l_max_resident
    {
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 2};
        RealSphericalHarmonicsData data;
        auto y40 = eng.compute(4, 0);
        check(std::isfinite(y40[0]), "compute(4,0) finite");
        auto y42 = eng.compute(4, 2);
        check(std::isfinite(y42[0]), "compute(4,2) finite");
    }

    // compute() at θ=0: Y_{l,0} = sqrt((2l+1)/(4π)), Y_{l,±m} = 0
    {
        double theta_zero = 0.0;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta_zero, 1),
                                   std::span<const double>(&phi, 1), 0};
        RealSphericalHarmonicsData data;

        auto y50 = eng.compute(5, 0);
        double expected = std::sqrt(11.0 / (4.0 * std::numbers::pi));
        check(near(y50[0], expected, eps), "compute(5,0) at θ=0");

        auto y53 = eng.compute(5, 3);
        check(near(y53[0], 0.0, eps), "compute(5,3) at θ=0 (=0)");
    }
}

auto test_cache_mode_none() -> void {
    std::println("\n=== CacheMode::None ===");

    double eps = 1e-12;
    double theta = 0.5;
    double phi = 0.3;

    RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1),
                               0, RealSphericalHarmonicsEngine::CacheMode::None};
    RealSphericalHarmonicsData data;

    // compute() works
    auto y00 = eng.compute(0, 0);
    check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps),
          "None mode compute Y_00");

    auto y10 = eng.compute(1, 0);
    double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta);
    check(near(y10[0], expected, eps), "None mode compute Y_10");

    // get() falls back to compute() in None mode
    {
        const auto& y00 = data.get(eng, 0, 0);
        check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps),
              "get() works in None mode via compute fallback");
    }
    {
        const auto& y10 = data.get(eng, 1, 0);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta);
        check(near(y10[0], expected, eps), "get() Y_10 in None mode");
    }
}

auto test_get_beyond_lmax() -> void {
    std::println("\n=== get() beyond lMaxResident ===");

    double theta = 0.5;
    double phi = 0.3;
    RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 2};
    RealSphericalHarmonicsData data;

    // get() falls back to compute() for l > l_max_resident
    auto y30 = data.get(eng, 3, 0);
    check(std::isfinite(y30[0]), "get(3,0) works via compute fallback");
}

auto test_reset_shrink_expand() -> void {
    std::println("\n=== reset(shrink) → reset(expand) ===");

    double eps = 1e-12;
    double theta_arr[] = {0.5};
    double phi_arr[] = {0.3};
    int l_max = 4;

    RealSphericalHarmonicsEngine eng{theta_arr, phi_arr, l_max};
    RealSphericalHarmonicsData data;

    // Save pre-shrink references
    int n4 = 2 * l_max + 1;
    int n5 = 2 * (l_max + 1) + 1;
    auto ref_low = data.get(eng, 2, 1);
    std::vector<std::vector<double>> ref_high(static_cast<std::size_t>(n4));
    std::vector<std::vector<double>> ref_beyond(static_cast<std::size_t>(n5));
    for (int m = -l_max; m <= l_max; ++m)
        ref_high[static_cast<std::size_t>(m + l_max)] = data.get(eng, l_max, m);
    for (int m = -(l_max + 1); m <= (l_max + 1); ++m)
        ref_beyond[static_cast<std::size_t>(m + l_max + 1)] = eng.compute(l_max + 1, m);

    // Shrink and verify remaining cache intact
    eng.setLMax(2);
    check(near(data.get(eng, 2, 1)[0], ref_low[0], eps),
          "get(2,1) unchanged after shrink");
    // get() falls back to compute() even after shrink
    auto y30 = data.get(eng, 3, 0);
    check(std::isfinite(y30[0]), "get(3,0) works after shrink via compute fallback");

    // Expand back and verify get() and compute() match pre-shrink
    eng.setLMax(l_max);
    for (int m = -l_max; m <= l_max; ++m) {
        const auto& y = data.get(eng, l_max, m);
        check(near(y[0], ref_high[static_cast<std::size_t>(m + l_max)][0], eps),
              std::format("get({},{}) matches pre-shrink after expand", l_max, m));
    }
    for (int m = -(l_max + 1); m <= (l_max + 1); ++m) {
        auto y = eng.compute(l_max + 1, m);
        check(near(y[0], ref_beyond[static_cast<std::size_t>(m + l_max + 1)][0], eps),
              std::format("compute({},{}) matches pre-shrink after expand", l_max + 1, m));
    }
}

auto main() -> int {
    try {
        std::println("=== Spherical Harmonics Tests ===");
        test_spherical_harmonics();
        test_spherical_harmonics_batch();
        test_compute();
        test_cache_mode_none();
        test_get_beyond_lmax();
        test_reset_shrink_expand();
        std::println("\nAll Spherical Harmonics tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
