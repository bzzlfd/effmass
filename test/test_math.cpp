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

    // j_0(x) = sin(x)/x
    check(near(sphericalBesselJ(0, 0.0), 1.0, 1e-10), "j_0(0) = 1");
    check(near(sphericalBesselJ(0, 1.0), std::sin(1.0)), "j_0(1) = sin(1)");
    check(near(sphericalBesselJ(0, std::numbers::pi / 2), std::sin(std::numbers::pi / 2) / (std::numbers::pi / 2)),
          "j_0(pi/2) = sin(pi/2)/(pi/2)");

    // j_1(x) = sin(x)/x^2 - cos(x)/x
    double x = 2.0;
    double j1_expected = std::sin(x) / (x * x) - std::cos(x) / x;
    check(near(sphericalBesselJ(1, x), j1_expected, eps), "j_1(2) analytic");

    // j_2(x) via recurrence
    double j0 = sphericalBesselJ(0, x);
    double j1 = sphericalBesselJ(1, x);
    double j2 = sphericalBesselJ(2, x);
    check(near(j2, 3.0 / x * j1 - j0, eps), "j_2(2) recurrence relation");

    // Small x series
    check(near(sphericalBesselJ(0, 0.01), std::sin(0.01) / 0.01, eps), "j_0(small x) series");
    double j3_small = sphericalBesselJ(3, 0.001);
    double j3_expected = std::pow(0.001, 3) / 105.0;
    check(near(j3_small, j3_expected, 1e-18), "j_3(0.001) = x^3/105");
}


auto test_spherical_harmonics() -> void {
    std::println("\n=== Real Spherical Harmonics ===");

    double eps = 1e-12;

    // Y_{00} = 1/(2*sqrt(pi))
    check(near(realSphericalHarmonic(0, 0, 0.0, 0.0), 1.0 / (2.0 * std::sqrt(std::numbers::pi)), eps),
          "Y_00 constant");

    // Y_{10}(theta) = sqrt(3/(4pi)) * cos(theta)
    double theta0 = 0.5;
    double y10 = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(theta0);
    check(near(realSphericalHarmonic(1, 0, theta0, 0.0), y10, eps), "Y_10(theta)");

    // Y_{1,1} = sqrt(3/(4pi)) * sin(theta) * cos(phi)
    double phi0 = 0.3;
    double y11 = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::cos(phi0);
    check(near(realSphericalHarmonic(1, 1, theta0, phi0), y11, eps), "Y_11(theta,phi)");

    // Y_{1,-1} = sqrt(3/(4pi)) * sin(theta) * sin(phi)
    double y1m1 = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(theta0) * std::sin(phi0);
    check(near(realSphericalHarmonic(1, -1, theta0, phi0), y1m1, eps), "Y_1,-1(theta,phi)");

    // m > l should return 0
    check(near(realSphericalHarmonic(2, 3, 0.5, 0.3), 0.0, eps), "Y_{2,3} = 0 (|m| > l)");
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


auto main() -> int {
    try {
        std::println("=== Math Module Tests ===");
        test_spherical_bessel();
        test_spherical_harmonics();
        test_fft1d();
        test_fft3d();
        std::println("\nAll tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
