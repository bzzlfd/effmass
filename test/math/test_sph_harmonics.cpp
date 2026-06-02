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
        auto y00 = sh(0, 0);
        check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps), "Y_00 constant");
    }

    // Y_{10}(theta) = sqrt(3/(4pi)) * cos(theta)
    {
        double theta0 = 0.5, phi0 = 0.0;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        auto y10 = sh(1, 0);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta0);
        check(near(y10[0], expected, eps), "Y_10(theta)");
    }

    // Y_{1,1} = sqrt(3/(4pi)) * sin(theta) * cos(phi)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        auto y11 = sh(1, 1);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::cos(phi0);
        check(near(y11[0], expected, eps), "Y_11(theta,phi)");
    }

    // Y_{1,-1} = sqrt(3/(4pi)) * sin(theta) * sin(phi)
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        auto y1m1 = sh(1, -1);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::sin(phi0);
        check(near(y1m1[0], expected, eps), "Y_1,-1(theta,phi)");
    }

    // m > l should return 0
    {
        double theta0 = 0.5, phi0 = 0.3;
        RealSphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1), 3};
        try {
            auto y23 = sh(2, 3);
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
    int n = (l_max + 1) * (l_max + 1);

    RealSphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};

    for (int l = 0; l <= l_max; ++l) {
        for (int m = -l; m <= l; ++m) {
            auto y_lm = sh(l, m);
            check(y_lm.size() == 1, std::format("batch Y_{{{},{}}} size", l, m));
        }
    }
}

auto main() -> int {
    try {
        std::println("=== Spherical Harmonics Tests ===");
        test_spherical_harmonics();
        test_spherical_harmonics_batch();
        std::println("\nAll Spherical Harmonics tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
