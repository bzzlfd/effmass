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

auto test_grad_phi() -> void {
    std::println("\n=== get_grad_phi() ===");

    double eps = 1e-14;

    // ---------- m = 0: (1/sinθ)·∂Y_l0/∂φ = 0 for any l ----------
    {
        double theta = 0.7, phi = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};
        for (int l = 0; l <= 4; ++l) {
            auto gp = eng.get_grad_phi(l, 0);
            check(near(gp[0], 0.0, eps),
                  std::format("get_grad_phi({},0) = 0", l));
        }
    }

    // ---------- m > 0: (1/sinθ)·∂Y_{lm}/∂φ = -m/sinθ · Y_{l,-m} ----------
    {
        double theta = 0.7, phi = 0.3;
        double sin_t = std::sin(theta);
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};
        RealSphericalHarmonicsData data;
        for (int l = 1; l <= 4; ++l) {
            for (int m = 1; m <= l; ++m) {
                auto gp = eng.get_grad_phi(l, m);
                const auto& y_l_neg_m = data.get(eng, l, -m);
                double expected = -static_cast<double>(m) / sin_t * y_l_neg_m[0];
                check(near(gp[0], expected, 1e-12),
                      std::format("grad_phi({},{}) = -{}/sinθ · Y({},{})", l, m, m, l, -m));
            }
        }
    }

    // ---------- m < 0: (1/sinθ)·∂Y_{lm}/∂φ = |m|/sinθ · Y_{l,|m|} ----------
    {
        double theta = 0.7, phi = 0.3;
        double sin_t = std::sin(theta);
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 4};
        RealSphericalHarmonicsData data;
        for (int l = 1; l <= 4; ++l) {
            for (int m = -l; m <= -1; ++m) {
                int m_abs = -m;
                auto gp = eng.get_grad_phi(l, m);
                const auto& y_l_pos_m = data.get(eng, l, m_abs);
                double expected = static_cast<double>(m_abs) / sin_t * y_l_pos_m[0];
                check(near(gp[0], expected, 1e-12),
                      std::format("grad_phi({},{}) = {}/sinθ · Y({},{})", l, m, m_abs, l, m_abs));
            }
        }
    }

    // ---------- Finite-difference check ----------
    {
        double theta = 0.7, phi_base = 0.3;
        double sin_t = std::sin(theta);
        double dphi = 1e-8;
        double phi_lo = phi_base - 0.5 * dphi;
        double phi_hi = phi_base + 0.5 * dphi;

        auto run = [&](double ph) -> RealSphericalHarmonicsEngine {
            return RealSphericalHarmonicsEngine{std::span<const double>(&theta, 1),
                                                std::span<const double>(&ph, 1), 4};
        };

        auto eng_lo = run(phi_lo);
        auto eng_hi = run(phi_hi);
        RealSphericalHarmonicsData data_lo, data_hi;

        for (int l = 0; l <= 4; ++l) {
            for (int m = -l; m <= l; ++m) {
                const auto& y_lo = data_lo.get(eng_lo, l, m);
                const auto& y_hi = data_hi.get(eng_hi, l, m);
                double fd = (y_hi[0] - y_lo[0]) / dphi;   // ∂Y/∂φ by central diff

                RealSphericalHarmonicsEngine eng_ref{std::span<const double>(&theta, 1),
                                                     std::span<const double>(&phi_base, 1), 4};
                double analytic = eng_ref.get_grad_phi(l, m)[0];
                // fd ≈ sinθ · analytic  →  fd/sinθ ≈ analytic
                check(near(fd, sin_t * analytic, 1e-6),
                      std::format("finite diff dY({},{})/dphi = sinθ·grad_phi", l, m));
            }
        }
    }

    // ---------- theta = 0 edge: |m|=1 finite, others zero ----------
    {
        double theta_zero = 0.0, phi = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta_zero, 1),
                                         std::span<const double>(&phi, 1), 4};

        // m = 0: always 0
        for (int l = 0; l <= 4; ++l) {
            auto gp = eng.get_grad_phi(l, 0);
            check(near(gp[0], 0.0, eps),
                  std::format("grad_phi({},0) = 0 at θ=0", l));
        }

        // l=1, m=±1: known analytic values at θ=0
        {
            double phi_val = phi;
            double expected_11 = std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::sin(phi_val);
            auto gp11 = eng.get_grad_phi(1, 1);
            check(near(gp11[0], expected_11, 1e-12),
                  "grad_phi(1,1) at θ=0 matches analytic limit");

            double expected_1m1 = -std::sqrt(3.0 / (4.0 * std::numbers::pi)) * std::cos(phi_val);
            auto gp1m1 = eng.get_grad_phi(1, -1);
            check(near(gp1m1[0], expected_1m1, 1e-12),
                  "grad_phi(1,-1) at θ=0 matches analytic limit");
        }

        // |m| >= 2: all vanish at θ=0
        for (int l = 2; l <= 4; ++l) {
            for (int m = -l; m <= l; ++m) {
                if (std::abs(m) < 2) continue;
                auto gp = eng.get_grad_phi(l, m);
                check(near(gp[0], 0.0, 1e-14),
                      std::format("grad_phi({},{}) = 0 at θ=0 (|m|>=2)", l, m));
            }
        }
    }

    // ---------- Mode=None throws ----------
    {
        double theta = 0.5, phi = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1),
                                   0, RealSphericalHarmonicsEngine::CacheMode::None};
        try {
            eng.get_grad_phi(0, 0);
            check(false, "get_grad_phi should throw in CacheMode::None");
        } catch (const std::runtime_error&) {
            check(true, "get_grad_phi throws in CacheMode::None");
        }
    }

    // ---------- l > l_max_resident throws ----------
    {
        double theta = 0.5, phi = 0.3;
        RealSphericalHarmonicsEngine eng{std::span<const double>(&theta, 1), std::span<const double>(&phi, 1), 2};
        try {
            eng.get_grad_phi(3, 0);
            check(false, "get_grad_phi should throw when l > l_max_resident");
        } catch (const std::runtime_error&) {
            check(true, "get_grad_phi throws when l > l_max_resident");
        }
    }
}

auto test_grad_phi_cache_lifecycle() -> void {
    std::println("\n=== get_grad_phi() cache lifecycle ===");

    double eps = 1e-14;
    double theta_arr[] = {0.5, 0.7};
    double phi_arr[] = {0.3, 0.6};

    // ---------- setLMax expand after get_grad_phi ----------
    {
        RealSphericalHarmonicsEngine eng{theta_arr, phi_arr, 2};
        // Fill cache for l_max=2
        auto gp20 = eng.get_grad_phi(2, 1);
        double ref = gp20[1];   // second G-vector (θ=0.7)
        check(std::isfinite(ref), "grad_phi(2,1) finite at l_max=2");

        // Expand and read — should work (cache repopulated lazily)
        eng.setLMax(4);
        auto gp31 = eng.get_grad_phi(3, 1);
        check(std::isfinite(gp31[0]), "grad_phi(3,1) works after expand");
        check(std::isfinite(gp31[1]), "grad_phi(3,1) works after expand (point 2)");
    }

    // ---------- reinit → setLMax → get_grad_phi ----------
    {
        RealSphericalHarmonicsEngine eng{theta_arr, phi_arr, 2};
        eng.get_grad_phi(2, 1);   // populate cache

        // Reinit with same theta/phi but different array
        eng.reinit(theta_arr, phi_arr);
        eng.setLMax(2);
        auto gp = eng.get_grad_phi(2, 1);
        check(std::isfinite(gp[0]), "grad_phi works after reinit + setLMax");
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
    // Pre-shrink grad_phi cache references
    auto g_ref_21 = eng.get_grad_phi(2, 1);
    auto g_ref_41 = eng.get_grad_phi(l_max, 1);
    auto g_ref_4m1 = eng.get_grad_phi(l_max, -1);

    // Shrink and verify remaining cache intact
    eng.setLMax(2);
    check(near(data.get(eng, 2, 1)[0], ref_low[0], eps),
          "get(2,1) unchanged after shrink");
    // get() falls back to compute() even after shrink
    auto y30 = data.get(eng, 3, 0);
    check(std::isfinite(y30[0]), "get(3,0) works after shrink via compute fallback");
    // grad_phi cache preserved (inner vectors for m_abs≥1 retained)
    check(near(eng.get_grad_phi(2, 1)[0], g_ref_21[0], eps),
          "grad_phi(2,1) unchanged after shrink");

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
    // grad_phi cache preserved through shrink/expand cycle
    check(near(eng.get_grad_phi(2, 1)[0], g_ref_21[0], eps),
          "grad_phi(2,1) matches pre-shrink after expand");
    check(near(eng.get_grad_phi(l_max, 1)[0], g_ref_41[0], eps),
          "grad_phi(l_max,1) matches pre-shrink after expand");
    check(near(eng.get_grad_phi(l_max, -1)[0], g_ref_4m1[0], eps),
          "grad_phi(l_max,-1) matches pre-shrink after expand");
}

auto main() -> int {
    try {
        std::println("=== Spherical Harmonics Tests ===");
        test_spherical_harmonics();
        test_spherical_harmonics_batch();
        test_compute();
        test_cache_mode_none();
        test_grad_phi();
        test_grad_phi_cache_lifecycle();
        test_get_beyond_lmax();
        test_reset_shrink_expand();
        std::println("\nAll Spherical Harmonics tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
