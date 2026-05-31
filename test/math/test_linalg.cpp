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


auto test_identity() -> void {
    std::println("\n=== Identity Matrix ===");
    double eps = 1e-14;

    std::vector a{1.0, 0.0, 0.0,
                  0.0, 2.0, 0.0,
                  0.0, 0.0, 3.0};
    auto res = diagonalize_jacobi(a, 3);

    check(near(res.eigenvalues[0], 1.0, eps), "λ₀ = 1");
    check(near(res.eigenvalues[1], 2.0, eps), "λ₁ = 2");
    check(near(res.eigenvalues[2], 3.0, eps), "λ₂ = 3");
}

auto test_2x2() -> void {
    std::println("\n=== 2×2 Matrix ===");
    double eps = 1e-14;

    // A = [[2, 1], [1, 2]]  →  λ = 1, 3
    std::vector a{2.0, 1.0,
                  1.0, 2.0};
    auto res = diagonalize_jacobi(a, 2);

    check(near(res.eigenvalues[0], 1.0, eps), "λ₀ = 1");
    check(near(res.eigenvalues[1], 3.0, eps), "λ₁ = 3");

    // Check A * v = λ * v
    auto& v = res.eigenvectors;
    for (int j = 0; j < 2; ++j) {
        double lhs0 = a[0] * v[0 * 2 + j] + a[2] * v[1 * 2 + j];
        double lhs1 = a[1] * v[0 * 2 + j] + a[3] * v[1 * 2 + j];
        double rhs0 = res.eigenvalues[j] * v[0 * 2 + j];
        double rhs1 = res.eigenvalues[j] * v[1 * 2 + j];
        check(near(lhs0, rhs0, eps) && near(lhs1, rhs1, eps),
              std::format("eigenvector {} satisfies A*v = λ*v", j));
    }
}

auto test_3x3() -> void {
    std::println("\n=== 3×3 Matrix ===");
    double eps = 1e-13;

    // A = [[3, 0, 0], [0, -1, 2], [0, 2, 2]]
    // Eigenvalues: 3, -2, 3
    std::vector a{3.0, 0.0, 0.0,
                  0.0, -1.0, 2.0,
                  0.0, 2.0, 2.0};
    auto res = diagonalize_jacobi(a, 3);

    check(near(res.eigenvalues[0], -2.0, eps), "λ₀ = -2");
    check(near(res.eigenvalues[1], 3.0, eps), "λ₁ = 3 (first)");
    check(near(res.eigenvalues[2], 3.0, eps), "λ₂ = 3 (second)");

    // Check A * V = V * Λ
    auto& v = res.eigenvectors;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            double lhs = 0.0;
            for (int k = 0; k < 3; ++k) {
                lhs += a[i * 3 + k] * v[k * 3 + j];
            }
            double rhs = res.eigenvalues[j] * v[i * 3 + j];
            check(near(lhs, rhs, eps),
                  std::format("A*v{} = λ*v{} at row {}", j, j, i));
        }
    }
}

auto test_random_5x5() -> void {
    std::println("\n=== Random 5×5 Matrix ===");
    double eps = 1e-12;

    int n = 5;
    std::vector a{ 4.0,  1.0,  0.0,  2.0, -1.0,
                   1.0,  3.0, -1.0,  0.0,  0.0,
                   0.0, -1.0,  5.0,  1.0,  2.0,
                   2.0,  0.0,  1.0,  2.0,  1.0,
                  -1.0,  0.0,  2.0,  1.0,  4.0};
    auto res = diagonalize_jacobi(a, n);

    // Check A * v_j = λ_j * v_j for each eigenvector
    auto& v = res.eigenvectors;
    for (int j = 0; j < n; ++j) {
        double max_err = 0.0;
        for (int i = 0; i < n; ++i) {
            double lhs = 0.0;
            for (int k = 0; k < n; ++k) {
                lhs += a[i * n + k] * v[k * n + j];
            }
            double rhs = res.eigenvalues[j] * v[i * n + j];
            max_err = std::max(max_err, std::abs(lhs - rhs));
        }
        check(max_err < eps, std::format("eigenvector {} residual < 1e-12 (max_err={:.2e})", j, max_err));
    }

    // Check orthogonality: V^T * V ≈ I
    double ortho_err = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dot = 0.0;
            for (int k = 0; k < n; ++k) dot += v[k * n + i] * v[k * n + j];
            double expected = (i == j) ? 1.0 : 0.0;
            ortho_err = std::max(ortho_err, std::abs(dot - expected));
        }
    }
    check(ortho_err < eps, std::format("V^T V ≈ I (max_err={:.2e})", ortho_err));
}

auto test_1x1() -> void {
    std::println("\n=== 1×1 Matrix ===");
    double eps = 1e-14;

    std::vector a{42.0};
    auto res = diagonalize_jacobi(a, 1);

    check(near(res.eigenvalues[0], 42.0, eps), "λ = 42");
    check(near(res.eigenvectors[0], 1.0, eps), "V = 1");
}

auto test_convergence() -> void {
    std::println("\n=== Convergence ===");

    std::vector a{2.0, 0.5,
                  0.5, 3.0};
    auto res = diagonalize_jacobi(a, 2);
    check(res.converged, "2×2 matrix converges");
}


auto main() -> int {
    try {
        std::println("=== Diagonalize Tests ===");
        test_identity();
        test_2x2();
        test_3x3();
        test_random_5x5();
        test_1x1();
        test_convergence();
        std::println("\nAll Diagonalize tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
