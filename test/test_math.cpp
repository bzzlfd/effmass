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

auto near_c(std::complex<double> a, std::complex<double> b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}

auto test_spherical_bessel() -> void {
    std::println("\n=== Spherical Bessel Functions ===");

    double eps = 1e-12;

    // SphericalBesselJ computes j_l(q * r_i). Use r={1} so x_i = q.
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

auto test_spherical_harmonics() -> void {
    std::println("\n=== Real Spherical Harmonics ===");

    double eps = 1e-12;

    // Y_{00} = 1/(2*sqrt(pi))
    {
        double theta = 0.0, phi = 0.0;
        SphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1)};
        auto y00 = sh(0, 0);
        check(near(y00[0], 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps), "Y_00 constant");
    }

    // Y_{10}(theta) = sqrt(3/(4pi)) * cos(theta)
    {
        double theta0 = 0.5, phi0 = 0.0;
        SphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1)};
        auto y10 = sh(1, 0);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta0);
        check(near(y10[0], expected, eps), "Y_10(theta)");
    }

    // Y_{1,1} = sqrt(3/(4pi)) * sin(theta) * cos(phi)
    {
        double theta0 = 0.5, phi0 = 0.3;
        SphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1)};
        auto y11 = sh(1, 1);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::cos(phi0);
        check(near(y11[0], expected, eps), "Y_11(theta,phi)");
    }

    // Y_{1,-1} = sqrt(3/(4pi)) * sin(theta) * sin(phi)
    {
        double theta0 = 0.5, phi0 = 0.3;
        SphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1)};
        auto y1m1 = sh(1, -1);
        double expected = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::sin(phi0);
        check(near(y1m1[0], expected, eps), "Y_1,-1(theta,phi)");
    }

    // m > l should return 0
    {
        double theta0 = 0.5, phi0 = 0.3;
        SphericalHarmonics sh{std::span<const double>(&theta0, 1), std::span<const double>(&phi0, 1)};
        try {
            auto y23 = sh(2, 3);  // throws: |m| > l
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

    SphericalHarmonics sh{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1)};

    for (int l = 0; l <= l_max; ++l) {
        for (int m = -l; m <= l; ++m) {
            auto y_lm = sh(l, m);
            // just verify it returns the right size
            check(y_lm.size() == 1, std::format("batch Y_{{{},{}}} size", l, m));
        }
    }
}

auto test_fft1d() -> void {
    std::println("\n=== 1D FFT ===");

    double eps = 1e-12;

    int n = 8;
    FFT3D fft8(1, 1, n);

    // Delta roundtrip
    std::vector<std::complex<double>> data(n);
    data[0] = 1.0;
    auto original = data;
    fft8(data, R2G);
    for (int i = 0; i < n; ++i) {
        check(near_c(data[i], {1.0 / n, 0.0}, eps), std::format("delta after R2G [{}]", i));
    }
    fft8(data, G2R);
    for (int i = 0; i < n; ++i) {
        check(near_c(data[i], original[i], 1e-14), std::format("delta FFT roundtrip [{}]", i));
    }

    // Cosine wave test
    std::vector<std::complex<double>> data2(n);
    for (int i = 0; i < n; ++i) {
        data2[i] = std::cos(2.0 * std::numbers::pi * i / n);
    }
    fft8(data2, R2G);
    double expected = 0.5;
    check(near(std::abs(data2[1]), expected, eps), "cos(2pi k/n) FFT peak at k=1");
    check(near(std::abs(data2[n - 1]), expected, eps), "cos(2pi k/n) FFT peak at k=n-1");

    // Constant array test: all ones -> peak at k=0
    std::vector<std::complex<double>> data_const(n, 1.0);
    fft8(data_const, R2G);
    check(near(std::abs(data_const[0]), 1.0, eps), "constant FFT peak at k=0");
    for (int i = 1; i < n; ++i) {
        check(near(std::abs(data_const[i]), 0.0, eps), std::format("constant FFT zero at k={}", i));
    }

    // n=2
    int n2 = 2;
    std::vector<std::complex<double>> data3 = {{1.0, 0.0}, {0.0, 0.0}};
    FFT3D fft2(1, 1, n2);
    fft2(data3, G2R);
    check(near_c(data3[0], {1.0, 0.0}, 1e-15), "n=2 inverse[0]");
    check(near_c(data3[1], {1.0, 0.0}, 1e-15), "n=2 inverse[1]");

    // Non-power-of-2: n=6
    int n6 = 6;
    std::vector<std::complex<double>> data6(n6);
    data6[0] = 1.0;
    original = data6;
    FFT3D fft6(1, 1, n6);
    fft6(data6, R2G);
    fft6(data6, G2R);
    for (int i = 0; i < n6; ++i) {
        check(near_c(data6[i], original[i], 1e-14), std::format("n=6 delta FFT roundtrip [{}]", i));
    }

    // Large non-power-of-2: n=45
    int n45 = 45;
    std::vector<std::complex<double>> data45(n45);
    for (int i = 0; i < n45; ++i) {
        data45[i] = std::complex<double>(i * 0.1, std::sin(0.1 * i));
    }
    original = data45;
    FFT3D fft45(1, 1, n45);
    fft45(data45, R2G);
    fft45(data45, G2R);
    for (int i = 0; i < n45; ++i) {
        check(near_c(data45[i], original[i], 1e-12), std::format("n=45 FFT roundtrip [{}]", i));
    }
}

auto test_fft3d() -> void {
    std::println("\n=== 3D FFT ===");

    int n1 = 6, n2 = 8, n3 = 5;
    int total = n1 * n2 * n3;
    FFT3D fft(n1, n2, n3);

    // Arbitrary grid roundtrip
    std::vector<std::complex<double>> grid(total);
    for (int i = 0; i < total; ++i) {
        grid[i] = std::complex<double>(i % 7, (i * i) % 11);
    }
    auto original = grid;

    fft(grid, G2R);
    fft(grid, R2G);

    for (int i = 0; i < total; ++i) {
        check(near_c(grid[i], original[i], 1e-12), std::format("3D FFT roundtrip [{}]", i));
    }

    // Delta function 3D
    std::vector<std::complex<double>> delta(total, 0.0);
    delta[0] = 1.0;
    auto delta_orig = delta;
    fft(delta, G2R);
    fft(delta, R2G);
    for (int i = 0; i < total; ++i) {
        check(near_c(delta[i], delta_orig[i], 1e-14), std::format("3D delta FFT roundtrip [{}]", i));
    }
}

auto test_beta_q_interpolator() -> void {
    std::println("\n=== Beta(q) Interpolator ===");

    // Gaussian-like beta(r) = exp(-r^2)
    int n = 401;
    double dr = 0.02;
    std::vector<double> r(n);
    std::vector<double> rab(n, dr);
    std::vector<double> beta(n);
    for (int i = 0; i < n; ++i) {
        r[i] = i * dr;
        beta[i] = std::exp(-r[i] * r[i]);
    }

    double dq = 0.01;
    double q_max = 10.0;

    // l=0 interpolator
    BetaQInterpolator interp0(r, rab, beta, 0, dq, q_max, RadialMeshType::Uniform);
    // Spot check: near a known value (integral of Gaussian-like function)
    double val0 = interp0(0.0);
    check(val0 > 0.0, "beta(0) > 0 for Gaussian-like l=0");

    double val1 = interp0(1.0);
    double val2 = interp0(2.0);
    check(val2 < val1, "beta(q) decreasing for l=0");

    // Derivative at q=0 should be 0 for even l
    check(near(interp0.derivative(0.0), 0.0, 1e-12), "d(beta)/dq = 0 at q=0 for l=0");

    // l=1 interpolator
    BetaQInterpolator interp1(r, rab, beta, 1, dq, q_max, RadialMeshType::Uniform);
    check(near(interp1(0.0), 0.0, 1e-14), "beta(0) = 0 for l=1 (j_1(0)=0)");

    // Table properties
    check(interp0.step() == dq, "step check");
    check(interp0.maxQ() == q_max, "q_max check");
    check(interp0.angularMomentum() == 0, "l check");
}

auto main() -> int {
    try {
        std::println("=== Math Module Tests ===");
        test_spherical_bessel();
        test_spherical_bessel_negative();
        test_spherical_harmonics();
        test_spherical_harmonics_batch();
        test_fft1d();
        test_fft3d();
        test_beta_q_interpolator();
        std::println("\nAll tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
