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

    double eps = 1e-12;
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
        std::println("=== FFT Tests ===");
        test_fft1d();
        test_fft3d();
        std::println("\nAll FFT tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
