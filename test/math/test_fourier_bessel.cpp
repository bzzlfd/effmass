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
    BetaqInterpolator interp0(r, rab, beta, 0, dq, q_max, RadialMeshType::Uniform);
    double val0 = interp0(0.0);
    check(val0 > 0.0, "beta(0) > 0 for Gaussian-like l=0");

    double val1 = interp0(1.0);
    double val2 = interp0(2.0);
    check(val2 < val1, "beta(q) decreasing for l=0");

    // Derivative at q=0 should be 0 for even l
    double eps = 1e-12;
    check(near(interp0.derivative(0.0), 0.0, 1e-12), "d(beta)/dq = 0 at q=0 for l=0");

    // l=1 interpolator
    BetaqInterpolator interp1(r, rab, beta, 1, dq, q_max, RadialMeshType::Uniform);
    check(near(interp1(0.0), 0.0, 1e-14), "beta(0) = 0 for l=1 (j_1(0)=0)");

    // Table properties
    check(interp0.step() == dq, "step check");
    check(interp0.maxQ() == q_max, "q_max check");
    check(interp0.angularMomentum() == 0, "l check");
}

auto main() -> int {
    try {
        std::println("=== Fourier-Bessel / Beta(q) Interpolator Tests ===");
        test_beta_q_interpolator();
        std::println("\nAll Fourier-Bessel tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
