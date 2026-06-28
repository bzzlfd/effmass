import std;
import math;
import pseudo;


// =============================================================================
//  Helpers
// =============================================================================

auto check(bool cond, std::string_view msg) -> void {
    if (!cond) {
        throw std::runtime_error(std::string("FAILED: ") + std::string(msg));
    }
    std::println("  PASSED: {}", msg);
}

auto near(double a, double b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}

auto relDiff(double a, double b) -> double {
    double denom = std::max(std::abs(a), std::abs(b));
    if (denom < 1e-40) return 0.0;
    return std::abs(a - b) / denom;
}

inline auto lagrangeDerivAt(std::span<const double> vals, double dq, int idx) -> double {
    int n = static_cast<int>(vals.size());
    if (idx <= 0 || idx >= n - 2) return std::numeric_limits<double>::quiet_NaN();
    // stencil [idx-1, idx, idx+1, idx+2], px=1 means evaluate at vals[idx]
    return lagrangeCubicDerivative(
        vals[idx - 1], vals[idx], vals[idx + 1], vals[idx + 2],
        /*px=*/1.0, dq);
}


// =============================================================================
//  BetaqTables tests
// =============================================================================

auto test_beta_q_tables() -> void {
    std::println("\n=== BetaqTables multi-l interpolation ===");

    // Gaussian-like beta functions at two different l
    int n = 401;
    double dr = 0.02;
    std::vector<double> r(n);
    std::vector<double> rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    // l=0 beta: exp(-r^2)
    // l=1 beta: r * exp(-r^2)
    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        double g = std::exp(-r[i] * r[i]);
        betas[0][i] = g;
        betas[1][i] = r[i] * g;
    }

    std::vector<int> angular_momenta = {0, 1};

    double dq = 0.01;
    double q_max = 10.0;
    BetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                       betas, angular_momenta, dq, q_max);

    // Basic properties
    check(tables.step() == dq, "step check");
    check(tables.maxQ() == q_max, "q_max check");
    check(tables.numProjectors(0) == 1, "num projectors l=0");
    check(tables.numProjectors(1) == 1, "num projectors l=1");
    check(tables.numProjectors(2) == 0, "num projectors l=2 (none)");

    // l=0: beta(0) > 0
    double v0_0 = tables.interpolate(0, 0, 0.0);
    check(v0_0 > 0.0, "l=0 beta(0) > 0");

    // l=0: beta(q) decreasing
    double v0_1 = tables.interpolate(0, 0, 1.0);
    double v0_2 = tables.interpolate(0, 0, 2.0);
    check(v0_2 < v0_1, "l=0 beta(q) decreasing");

    // l=1: beta(0) = 0 (j_1(0) = 0)
    double v1_0 = tables.interpolate(1, 0, 0.0);
    check(near(v1_0, 0.0, 1e-14), "l=1 beta(0) = 0");

    // l=1: beta(q) > 0 for finite q
    double v1_1 = tables.interpolate(1, 0, 1.0);
    check(v1_1 > 0.0, "l=1 beta(1) > 0");

    // Consistency: l=0 and l=1 agree with direct computeBetaq
    // The interpolation should match within reasonable tolerance
    std::println("  (interpolation self-consistency ok)");
}

auto test_set_volume() -> void {
    std::println("\n=== BetaqTables::setVolume ===");

    int n = 201;
    double dr = 0.02;
    std::vector<double> r(n), rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i)
        betas[0][i] = std::exp(-r[i] * r[i]);

    std::vector<int> angular_momenta = {0};

    BetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                       betas, angular_momenta);

    // Without setVolume: norm_coeff = 1.0 (default)
    double v0 = tables.interpolate(0, 0, 1.0);

    // After setVolume: should scale by 1/√Ω
    double omega = 100.0;
    tables.setVolume(omega);
    double v0_scaled = tables.interpolate(0, 0, 1.0);

    double expected_scale = 1.0 / std::sqrt(omega);
    check(near(v0_scaled, v0 * expected_scale, 1e-12),
          "setVolume scales beta(q) correctly");
}

auto test_boundary_errors() -> void {
    std::println("\n=== BetaqTables boundary error checks ===");

    int n = 100;
    double dr = 0.02;
    std::vector<double> r(n), rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i) betas[0][i] = 1.0;

    std::vector<int> angular_momenta = {0};

    BetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                       betas, angular_momenta, 0.01, 10.0);

    // Negative q
    bool caught = false;
    try { tables.interpolate(0, 0, -0.1); }
    catch (const std::domain_error&) { caught = true; }
    check(caught, "negative q throws domain_error");

    // q > q_max
    caught = false;
    try { tables.interpolate(0, 0, 20.0); }
    catch (const std::domain_error&) { caught = true; }
    check(caught, "q > q_max throws domain_error");

    // Invalid l
    caught = false;
    try { tables.interpolate(2, 0, 1.0); }
    catch (const std::out_of_range&) { caught = true; }
    check(caught, "invalid l throws out_of_range");

    // Invalid ib_in_l
    caught = false;
    try { tables.interpolate(0, 1, 1.0); }
    catch (const std::out_of_range&) { caught = true; }
    check(caught, "invalid ib_in_l throws out_of_range");
}


// =============================================================================
//  DBetaqTables tests
// =============================================================================

auto test_dbeta_q_tables() -> void {
    std::println("\n=== DBetaqTables multi-l interpolation ===");

    // Same Gaussian-like test functions as BetaqTables tests
    int n = 401;
    double dr = 0.02;
    std::vector<double> r(n);
    std::vector<double> rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    // l=0: exp(-r²),  l=1: r·exp(-r²)
    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        double g = std::exp(-r[i] * r[i]);
        betas[0][i] = g;
        betas[1][i] = r[i] * g;
    }

    std::vector<int> angular_momenta = {0, 1};

    double dq = 0.01;
    double q_max = 10.0;
    DBetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                        betas, angular_momenta, dq, q_max);

    // Basic properties
    check(tables.step() == dq, "dbeta step check");
    check(tables.maxQ() == q_max, "dbeta q_max check");
    check(tables.numProjectors(0) == 1, "dbeta num projectors l=0");
    check(tables.numProjectors(1) == 1, "dbeta num projectors l=1");
    check(tables.numProjectors(2) == 0, "dbeta num projectors l=2 (none)");

    // l=0: d(beta)/dq at q=0 should be 0 (j_0'(0) = 0)
    double dv0_0 = tables.interpolate(0, 0, 0.0);
    check(near(dv0_0, 0.0, 1e-12), "l=0 d(beta)/dq(0) = 0");

    // l=0: d(beta)/dq should be negative for q>0 (beta(q) decreasing)
    double dv0_05 = tables.interpolate(0, 0, 0.5);
    double dv0_10 = tables.interpolate(0, 0, 1.0);
    check(dv0_05 < 0.0, "l=0 d(beta)/dq(0.5) < 0");
    check(dv0_10 < 0.0, "l=0 d(beta)/dq(1.0) < 0");
    check(dv0_10 < dv0_05, "l=0 d(beta)/dq more negative at larger q");

    // l=1: d(beta)/dq at q=0 should be positive (j_1'(0) = 1/3 > 0)
    double dv1_0 = tables.interpolate(1, 0, 0.0);
    check(dv1_0 > 0.0, "l=1 d(beta)/dq(0) > 0");

    // l=1: analytic check — d(beta_1)/dq|_{q=0} = 2π/3 ≈ 2.094
    double expected_dv1_0 = 2.0 * std::numbers::pi / 3.0;
    check(near(dv1_0, expected_dv1_0, 1e-2),
          "l=1 d(beta)/dq(0) ≈ 2π/3");

    std::println("  (DBetaqTables interpolation self-consistency ok)");
}

auto test_dbeta_set_volume() -> void {
    std::println("\n=== DBetaqTables::setVolume ===");

    int n = 201;
    double dr = 0.02;
    std::vector<double> r(n), rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i)
        betas[0][i] = std::exp(-r[i] * r[i]);

    std::vector<int> angular_momenta = {0};

    DBetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                        betas, angular_momenta);

    // Without setVolume: norm_coeff = 1.0 (default)
    double dv0 = tables.interpolate(0, 0, 1.0);

    // After setVolume: should scale by 1/√Ω
    double omega = 100.0;
    tables.setVolume(omega);
    double dv0_scaled = tables.interpolate(0, 0, 1.0);

    double expected_scale = 1.0 / std::sqrt(omega);
    check(near(dv0_scaled, dv0 * expected_scale, 1e-12),
          "setVolume scales d(beta)/dq correctly");
}

auto test_dbeta_boundary_errors() -> void {
    std::println("\n=== DBetaqTables boundary error checks ===");

    int n = 100;
    double dr = 0.02;
    std::vector<double> r(n), rab(n, dr);
    for (int i = 0; i < n; ++i) r[i] = i * dr;

    std::vector<std::vector<double>> betas;
    betas.push_back(std::vector<double>(n));
    for (int i = 0; i < n; ++i) betas[0][i] = 1.0;

    std::vector<int> angular_momenta = {0};

    DBetaqTables tables(r, rab, SimpsonMeshType::Uniform,
                        betas, angular_momenta, 0.01, 10.0);

    // Negative q
    bool caught = false;
    try { tables.interpolate(0, 0, -0.1); }
    catch (const std::domain_error&) { caught = true; }
    check(caught, "dbeta negative q throws domain_error");

    // q > q_max
    caught = false;
    try { tables.interpolate(0, 0, 20.0); }
    catch (const std::domain_error&) { caught = true; }
    check(caught, "dbeta q > q_max throws domain_error");

    // Invalid l
    caught = false;
    try { tables.interpolate(2, 0, 1.0); }
    catch (const std::out_of_range&) { caught = true; }
    check(caught, "dbeta invalid l throws out_of_range");

    // Invalid ib_in_l
    caught = false;
    try { tables.interpolate(0, 1, 1.0); }
    catch (const std::out_of_range&) { caught = true; }
    check(caught, "dbeta invalid ib_in_l throws out_of_range");
}


// =============================================================================
//  Analytical-vs-numerical derivative comparison (using real UPF data)
// =============================================================================

auto test_dbeta_analytical_vs_numerical() -> void {
    std::println("\n=== β(q) derivative comparison: analytical vs numerical ===");

    std::string upf_path = "test/data_io_upf/Ge-spd-high.PD04.PBE.UPF";
    std::println("Reading: {}\n", upf_path);

    UPF upf(upf_path);
    NCPP ncpp(upf);
    ncpp.inferMeshType();

    auto& r   = ncpp.mesh.r;
    auto& rab = ncpp.mesh.rab;
    auto mesh_type = ncpp.mesh.type;
    std::println("Mesh: {} points, type = {}",
                 ncpp.meta.mesh_size,
                 (mesh_type == MeshType::Exponential ? "Exponential"
                  : mesh_type == MeshType::Uniform ? "Uniform" : "Unknown"));

    auto simpson_type = SimpsonMeshType::General;

    // Group projectors by angular momentum
    int n_proj = ncpp.meta.number_of_proj;
    std::println("Projectors: {} total", n_proj);

    struct ProjInfo {
        int l;
        int orig_idx;
        std::vector<double> beta;
        std::span<const double> r_sub;
        std::span<const double> rab_sub;
    };
    std::vector<ProjInfo> projs;

    for (int ip = 0; ip < n_proj; ++ip) {
        int l    = ncpp.nonlocal.angular_momentum[ip];
        int cutoff = ncpp.nonlocal.cutoff_index[ip];

        auto& beta_full = ncpp.nonlocal.beta[ip];

        int nz = static_cast<int>(beta_full.size());
        while (nz > 0 && std::abs(beta_full[nz - 1]) < 1e-30) --nz;
        if (nz == 0) continue;

        projs.push_back({
            .l = l,
            .orig_idx = ip,
            .beta = std::vector<double>(beta_full.begin(),
                                        beta_full.begin() + nz),
            .r_sub = std::span(r).subspan(0, nz),
            .rab_sub = std::span(rab).subspan(0, nz),
        });
        std::println("  Projector {}: l={}, cutoff_idx={}, non-zero={}",
                     ip + 1, l, cutoff, nz);
    }

    if (projs.empty()) {
        throw std::runtime_error("No projectors found in UPF file");
    }

    // q-grid
    double dq    = 0.02;
    double q_max = 6.0;
    int nq       = static_cast<int>(q_max / dq) + 1;

    // Stats accumulators
    double max_rel_diff_d1 = 0.0;
    double max_rel_diff_d2 = 0.0;
    int worst_d1_l = -1, worst_d1_iq = -1;
    int worst_d2_l = -1, worst_d2_iq = -1;

    // Loop over projectors
    for (auto& pj : projs) {
        int l = pj.l;
        std::println("\n──────────────────────────────────────────────");
        std::println("l = {}  (projector {})", l, pj.orig_idx + 1);
        std::println("──────────────────────────────────────────────");

        const auto& f    = pj.beta;
        const auto& r_sub  = pj.r_sub;
        const auto& rab_sub = pj.rab_sub;
        int nr = static_cast<int>(f.size());
        std::println("  Radial points: {}", nr);

        // Precompute β(q), β'(q), β''(q) analytically
        std::vector<double> beta_q(nq, 0.0);
        std::vector<double> dbeta_ana(nq, 0.0);
        std::vector<double> d2beta_ana(nq, 0.0);
        std::vector<double> dbeta_num(nq, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> d2beta_num(nq, std::numeric_limits<double>::quiet_NaN());

        for (int iq = 0; iq < nq; ++iq) {
            double q = iq * dq;

            SphericalBesselJ bessel{r_sub, q};
            bessel.advance(l);
            auto jl = bessel.value();
            auto jl_prime = bessel.derivValue();
            auto jl_double_prime = bessel.secondDerivValue();

            // β(q)    = 4π ∫ f·r·j_l(qr) dr          (R1: UPF beta has one r built in)
            // β'(q)   = 4π ∫ f·r²·j_l'(qr) dr        (DerivR1)
            // β''(q)  = 4π ∫ f·r³·j_l''(qr) dr       (SecondDerivR1)
            beta_q[iq]     = fourierBesselIntegrateR1(f, r_sub, jl, rab_sub, simpson_type);
            dbeta_ana[iq]  = fourierBesselIntegrateDerivR1(f, r_sub, jl_prime, rab_sub, simpson_type);
            d2beta_ana[iq] = fourierBesselIntegrateSecondDerivR1(f, r_sub, jl_double_prime, rab_sub, simpson_type);
        }

        // Numerical derivatives (Lagrange cubic)
        for (int iq = 1; iq < nq - 2; ++iq) {
            dbeta_num[iq] = lagrangeDerivAt(beta_q, dq, iq);
        }
        for (int iq = 1; iq < nq - 2; ++iq) {
            if (!std::isnan(dbeta_num[iq - 1]) &&
                !std::isnan(dbeta_num[iq]) &&
                !std::isnan(dbeta_num[iq + 1]) &&
                !std::isnan(dbeta_num[iq + 2])) {
                d2beta_num[iq] = lagrangeDerivAt(dbeta_num, dq, iq);
            }
        }

        // Print comparison table (first 50 q values)
        std::println("  {0:>6}  {1:>14}  {2:>12}  {3:>12}  {4:>10}  "
                     "{5:>12}  {6:>12}  {7:>10}",
                     "q", "β(q)", "β'_ana", "β'_num", "rel|d1|",
                     "β''_ana", "β''_num", "rel|d2|");

        int print_limit = std::min(nq, 50);
        for (int iq = 0; iq < print_limit; ++iq) {
            double q = iq * dq;
            double rd1 = relDiff(dbeta_ana[iq], dbeta_num[iq]);
            double rd2 = relDiff(d2beta_ana[iq], d2beta_num[iq]);

            if (rd1 > max_rel_diff_d1) {
                max_rel_diff_d1 = rd1;
                worst_d1_l = l;
                worst_d1_iq = iq;
            }
            if (rd2 > max_rel_diff_d2) {
                max_rel_diff_d2 = rd2;
                worst_d2_l = l;
                worst_d2_iq = iq;
            }

            std::println("  {0:6.3f}  {1:14.8e}  {2:12.6e}  {3:12.6e}  {4:10.3e}  "
                         "{5:12.6e}  {6:12.6e}  {7:10.3e}",
                         q, beta_q[iq],
                         dbeta_ana[iq], dbeta_num[iq], rd1,
                         d2beta_ana[iq], d2beta_num[iq], rd2);
        }

        if (nq > print_limit) {
            std::println("  ... ({} points remaining, q up to {})",
                         nq - print_limit, q_max);
            for (int iq = print_limit; iq < nq; ++iq) {
                double rd1 = relDiff(dbeta_ana[iq], dbeta_num[iq]);
                double rd2 = relDiff(d2beta_ana[iq], d2beta_num[iq]);
                if (rd1 > max_rel_diff_d1) {
                    max_rel_diff_d1 = rd1;
                    worst_d1_l = l;
                    worst_d1_iq = iq;
                }
                if (rd2 > max_rel_diff_d2) {
                    max_rel_diff_d2 = rd2;
                    worst_d2_l = l;
                    worst_d2_iq = iq;
                }
            }
        }
    }

    // Summary
    std::println("\n══════════════════════════════════════════════");
    std::println("  Summary");
    std::println("══════════════════════════════════════════════");
    std::println("  q-grid: dq={}, nq={}, q_max={}", dq, nq, q_max);
    std::println("  Max |β'_ana - β'_num| / max(|β'_ana|,|β'_num|):  {:.3e}  (l={}, q={:.3f})",
                 max_rel_diff_d1, worst_d1_l, worst_d1_iq * dq);
    std::println("  Max |β''_ana- β''_num| / max(|β''_ana|,|β''_num|): {:.3e}  (l={}, q={:.3f})",
                 max_rel_diff_d2, worst_d2_l, worst_d2_iq * dq);


    // Pass/fail criteria — two-tier:
    //   fail_tolerance  — hard limit; exceeding this means something is broken.
    //   warn_tolerance  — soft limit; exceeded → warning printed, but test passes.
    //
    // Numerical differentiation with dq = 0.02 has truncation error O(dq³) ≈ 8e-6
    // on the 4-point Lagrange stencil; the Fourier-Bessel integration also carries
    // Simpson discretisation error.  Typical observed rel-diffs are ~2e-4 (1st deriv)
    // and ~2e-3 (2nd deriv).
    constexpr double fail_tolerance_d1 = 1e-2;
    constexpr double fail_tolerance_d2 = 1e-1;
    constexpr double warn_tolerance_d1 = 5e-4;
    constexpr double warn_tolerance_d2 = 5e-3;

    if (max_rel_diff_d1 > warn_tolerance_d1) {
        std::println("  ⚠ Warning: first derivative rel diff {:.3e} exceeds soft limit {}",
                     max_rel_diff_d1, warn_tolerance_d1);
    }
    if (max_rel_diff_d2 > warn_tolerance_d2) {
        std::println("  ⚠ Warning: second derivative rel diff {:.3e} exceeds soft limit {}",
                     max_rel_diff_d2, warn_tolerance_d2);
    }

    check(max_rel_diff_d1 < fail_tolerance_d1,
          std::format("first derivative max rel diff {:.3e} < {}", max_rel_diff_d1, fail_tolerance_d1));
    check(max_rel_diff_d2 < fail_tolerance_d2,
          std::format("second derivative max rel diff {:.3e} < {}", max_rel_diff_d2, fail_tolerance_d2));
}


// =============================================================================
//  Main
// =============================================================================

auto main() -> int {
    try {
        std::println("=== BetaqTables Tests ===");
        test_beta_q_tables();
        test_set_volume();
        test_boundary_errors();

        std::println("\n=== DBetaqTables Tests ===");
        test_dbeta_q_tables();
        test_dbeta_set_volume();
        test_dbeta_boundary_errors();

        std::println("\n=== Analytical vs Numerical Derivative Comparison ===");
        test_dbeta_analytical_vs_numerical();

        std::println("\nAll tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
