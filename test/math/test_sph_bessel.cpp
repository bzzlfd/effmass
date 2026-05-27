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


auto test_spherical_bessel() -> void {
    std::println("\n=== Spherical Bessel Functions ===");

    double eps = 1e-12;

    double r_one = 1.0;
    auto r1 = std::span<const double>(&r_one, 1);

    // j_0(x) = sin(x)/x
    {
        SphericalBesselJ bj{r1, 0.0};
        check(near(bj.value()[0], 1.0, 1e-10), "j_0(0) = 1");
    }
    {
        SphericalBesselJ bj{r1, 1.0};
        check(near(bj.value()[0], std::sin(1.0)), "j_0(1) = sin(1)");
    }
    {
        double x = std::numbers::pi / 2;
        SphericalBesselJ bj{r1, x};
        check(near(bj.value()[0], std::sin(x) / x), "j_0(pi/2) = sin(pi/2)/(pi/2)");
    }

    // j_1(x) = sin(x)/x^2 - cos(x)/x
    {
        double x = 2.0;
        SphericalBesselJ bj{r1, x};
        bj.advance(1);
        double j1_expected = std::sin(x) / (x * x) - std::cos(x) / x;
        check(near(bj.value()[0], j1_expected, eps), "j_1(2) analytic");
    }

    // j_2(x) via recurrence
    {
        double x = 2.0;
        SphericalBesselJ bj{r1, x};
        double j0 = bj.value()[0];
        bj.advance(1);
        double j1 = bj.value()[0];
        bj.advance(1);
        double j2 = bj.value()[0];
        check(near(j2, 3.0 / x * j1 - j0, eps), "j_2(2) recurrence relation");
    }

    // Small x series
    {
        SphericalBesselJ bj{r1, 0.01};
        check(near(bj.value()[0], std::sin(0.01) / 0.01, eps), "j_0(small x) series");
    }
    {
        SphericalBesselJ bj{r1, 0.001};
        bj.advance(3);
        double j3_expected = std::pow(0.001, 3) / 105.0;
        check(near(bj.value()[0], j3_expected, 1e-18), "j_3(0.001) = x^3/105");
    }
}

auto test_spherical_bessel_negative() -> void {
    std::println("\n=== Spherical Bessel (negative x) ===");

    double eps = 1e-12;
    double x = 1.5;
    double r_one = 1.0;
    auto r1 = std::span<const double>(&r_one, 1);

    // Even l: j_l(-x) = j_l(x)
    {
        SphericalBesselJ bj_pos{r1, x};
        SphericalBesselJ bj_neg{r1, -x};
        check(near(bj_neg.value()[0], bj_pos.value()[0], eps), "j_0 even");
        bj_pos.advance(2); bj_neg.advance(2);
        check(near(bj_neg.value()[0], bj_pos.value()[0], eps), "j_2 even");
        bj_pos.advance(2); bj_neg.advance(2);
        check(near(bj_neg.value()[0], bj_pos.value()[0], eps), "j_4 even");
    }

    // Odd l: j_l(-x) = -j_l(x)
    {
        SphericalBesselJ bj_pos{r1, x};
        SphericalBesselJ bj_neg{r1, -x};
        bj_pos.advance(1); bj_neg.advance(1);
        check(near(bj_neg.value()[0], -bj_pos.value()[0], eps), "j_1 odd");
        bj_pos.advance(2); bj_neg.advance(2);
        check(near(bj_neg.value()[0], -bj_pos.value()[0], eps), "j_3 odd");
        bj_pos.advance(2); bj_neg.advance(2);
        check(near(bj_neg.value()[0], -bj_pos.value()[0], eps), "j_5 odd");
    }
}

auto main() -> int {
    try {
        std::println("=== Spherical Bessel Tests ===");
        test_spherical_bessel();
        test_spherical_bessel_negative();
        std::println("\nAll Spherical Bessel tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
