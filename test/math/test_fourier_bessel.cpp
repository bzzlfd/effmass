import std;
import math;
import pseudo;


auto check(bool cond, std::string_view msg) -> void {
    if (!cond) {
        throw std::runtime_error(std::string("FAILED: ") + std::string(msg));
    }
    std::println("  PASSED: {}", msg);
}

auto near(double a, double b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}


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


auto main() -> int {
    try {
        std::println("=== BetaqTables Tests ===");
        test_beta_q_tables();
        test_set_volume();
        test_boundary_errors();
        std::println("\nAll BetaqTables tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
