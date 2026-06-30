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
        RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 3};
        const auto& y00 = sh.get(0, 0);
        check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps), "Y_00 constant");
    }

    // Y_{10}(theta) = sqrt(3/(4pi)) * cos(theta)
    {
        double theta0 = 0.5, phi0 = 0.0;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        const auto& y10 = sh.get(1, 0);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta0);
        check(near(y10[0], expected, eps), "Y_10(theta)");
    }

    // Y_{1,1} = -sqrt(3/(4pi)) * sin(theta) * cos(phi)
    // (Condon-Shortley factor (-1)^1 = -1)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        const auto& y11 = sh.get(1, 1);
        double expected = -std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::cos(phi0);
        check(near(y11[0], expected, eps), "Y_11(theta,phi) (CS)");
    }

    // Y_{1,-1} = -sqrt(3/(4pi)) * sin(theta) * sin(phi)
    // (Condon-Shortley factor (-1)^1 = -1)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        const auto& y1m1 = sh.get(1, -1);
        double expected = -std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::sin(phi0);
        check(near(y1m1[0], expected, eps), "Y_1,-1(theta,phi) (CS)");
    }

    // |m| > l should throw
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        try {
            sh.get(2, 3);
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

    RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};

    for (int l = 0; l <= l_max; ++l) {
        for (int m = -l; m <= l; ++m) {
            const auto& y_lm = sh.get(l, m);
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
        RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 3};
        for (int l = 0; l <= 3; ++l) {
            for (int m = -l; m <= l; ++m) {
                const auto& cached = sh.get(l, m);
                auto computed = sh.compute(l, m);
                check(near(cached[0], computed[0], eps),
                      std::format("compute({},{}) matches get", l, m));
            }
        }
    }

    // compute() works for l > l_max_resident
    {
        RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 2};
        auto y40 = sh.compute(4, 0);
        check(std::isfinite(y40[0]), "compute(4,0) finite");
        auto y42 = sh.compute(4, 2);
        check(std::isfinite(y42[0]), "compute(4,2) finite");
    }

    // compute() at θ=0: Y_{l,0} = sqrt((2l+1)/(4π)), Y_{l,±m} = 0
    {
        double theta_zero = 0.0;
        RealSphericalHarmonics sh{std::span<const double>(&theta_zero, 1),
                                   std::span<const double>(&phi, 1), 0};

        auto y50 = sh.compute(5, 0);
        double expected = std::sqrt(11.0 / (4.0 * std::numbers::pi));
        check(near(y50[0], expected, eps), "compute(5,0) at θ=0");

        auto y53 = sh.compute(5, 3);
        check(near(y53[0], 0.0, eps), "compute(5,3) at θ=0 (=0)");
    }
}

auto test_cache_mode_none() -> void {
    std::println("\n=== CacheMode::None ===");

    double eps = 1e-12;
    double theta = 0.5;
    double phi = 0.3;

    RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1),
                               0, RealSphericalHarmonics::CacheMode::None};

    // compute() works
    auto y00 = sh.compute(0, 0);
    check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps),
          "None mode compute Y_00");

    auto y10 = sh.compute(1, 0);
    double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta);
    check(near(y10[0], expected, eps), "None mode compute Y_10");

    // get() throws
    try {
        sh.get(0, 0);
        check(false, "get() should throw in None mode");
    } catch (const std::exception&) {
        check(true, "get() throws in CacheMode::None");
    }
}

auto test_get_out_of_range() -> void {
    std::println("\n=== get() out of range ===");

    double theta = 0.5;
    double phi = 0.3;
    RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 2};

    // l > l_max_resident should throw
    try {
        sh.get(3, 0);
        check(false, "get(3,0) should have thrown (l > l_max_resident)");
    } catch (const std::runtime_error&) {
        check(true, "get(3,0) throws (l > l_max_resident)");
    }
}

auto test_reset_shrink_expand() -> void {
    std::println("\n=== reset(shrink) → reset(expand) ===");

    double eps = 1e-12;
    double theta_arr[] = {0.5};
    double phi_arr[] = {0.3};
    int l_max = 4;

    RealSphericalHarmonics sh{theta_arr, phi_arr, l_max};

    // Save pre-shrink references
    int n4 = 2 * l_max + 1;
    int n5 = 2 * (l_max + 1) + 1;
    auto ref_low = sh.get(2, 1);
    std::vector<std::vector<double>> ref_high(static_cast<std::size_t>(n4));
    std::vector<std::vector<double>> ref_beyond(static_cast<std::size_t>(n5));
    for (int m = -l_max; m <= l_max; ++m)
        ref_high[static_cast<std::size_t>(m + l_max)] = sh.get(l_max, m);
    for (int m = -(l_max + 1); m <= (l_max + 1); ++m)
        ref_beyond[static_cast<std::size_t>(m + l_max + 1)] = sh.compute(l_max + 1, m);

    // Shrink and verify remaining cache intact
    sh.setLMax(2);
    check(near(sh.get(2, 1)[0], ref_low[0], eps),
          "get(2,1) unchanged after shrink");
    try {
        sh.get(3, 0);
        check(false, "get(3,0) should throw after shrink");
    } catch (const std::runtime_error&) {
        check(true, "get(3,0) throws after shrink");
    }

    // Expand back and verify get() and compute() match pre-shrink
    sh.setLMax(l_max);
    for (int m = -l_max; m <= l_max; ++m) {
        const auto& y = sh.get(l_max, m);
        check(near(y[0], ref_high[static_cast<std::size_t>(m + l_max)][0], eps),
              std::format("get({},{}) matches pre-shrink after expand", l_max, m));
    }
    for (int m = -(l_max + 1); m <= (l_max + 1); ++m) {
        auto y = sh.compute(l_max + 1, m);
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
        test_get_out_of_range();
        test_reset_shrink_expand();
        std::println("\nAll Spherical Harmonics tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
