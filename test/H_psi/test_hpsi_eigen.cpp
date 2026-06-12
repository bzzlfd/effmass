// Test: H|ψ⟩ = E|ψ⟩ eigenstate consistency.
//
// For selected (k-point, band) pairs:
//   1. Load ψ from OUT.WG
//   2. Compute H|ψ⟩ via the Callable framework
//   3. Extract eigenvalue via Rayleigh quotient:  E = ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩
//   4. Compare with reference eigenvalue from OUT.EIGEN
//   5. Compute linear regression a·ψ + b = H|ψ⟩ for quality diagnostics
//   6. Check wavefunction normalisation
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)
// Uses test data under test/data_io-nonlocal/

import std;
import io;
import math;
import H_psi;

auto main() -> int {
    Hamiltonian h("test/data_io-nonlocal");
    h.loadFromDirectory();

    auto& gkk   = h.gkk();
    auto& wg    = h.wg();
    const auto& eigen = h.eigen();
    double volume = gkk.meta.lattice.volume();

    int nband = wg.meta.nband;

    // --- Rayeligh helper (generic span) ---
    struct Sums { std::complex<double> num; double den; };
    auto rayleigh = [](auto psi, auto hpsi) -> std::complex<double> {
        Sums s{0.0, 0.0};
        for (std::size_t i = 0; i < psi.size(); ++i) {
            s.num += std::conj(psi[i]) * hpsi[i];
            s.den += std::norm(psi[i]);
        }
        return s.num / s.den;
    };

    // --- Select test pairs: all bands at ikpt=0 and ikpt=1 ---
    struct TestPair { int ikpt; int iband; };
    std::vector<TestPair> tests;

    for (int k : {0, 1}) {
        if (k >= gkk.meta.nkpt) continue;
        for (int b = 0; b < nband; ++b)
            tests.push_back({k, b});
    }

    std::println(
        "Testing {} (ikpt, iband) pairs  |  nband={}  nkpt={}  volume={:.2f} Bohr³\n",
        tests.size(), nband, gkk.meta.nkpt, volume);

    // Column headers
    struct Result { int ikpt{}, iband{}; double e_file{}, e_rayleigh{}, delta{}, im_e{}, rel_res{}, norm_check{}; bool pass{}; };
    std::vector<Result> results;
    results.reserve(tests.size());

    for (auto [ikpt, iband] : tests) {
        auto op = h.at_k(ikpt);
        int ng = op.dim();

        auto psi = wg.loadBand(ikpt, iband).up;

        // Normalization
        double psi_norm2 = 0.0;
        for (int ig = 0; ig < ng; ++ig)
            psi_norm2 += std::norm(psi[ig]);
        double norm_check = psi_norm2 * volume;

        // H|ψ⟩
        std::vector<std::complex<double>> hpsi(static_cast<std::size_t>(ng), 0.0);
        op(psi, hpsi);

        // Rayleigh quotient
        auto e_rayleigh = rayleigh(psi, hpsi);
        double e_file = eigen[iband, ikpt];

        // 2-param linear regression: a·ψ + b = H|ψ⟩
        // a = exact Rayleigh (for eigenstate, same as one-term Rayleigh)
        // b = constant intercept, residual quality measure
        std::complex<double> sum_psi_conj_hpsi = std::conj(psi[0]) * hpsi[0];
        std::complex<double> sum_psi_conj = std::conj(psi[0]);
        std::complex<double> sum_psi = psi[0];
        std::complex<double> sum_hpsi = hpsi[0];
        double sum_psi2 = std::norm(psi[0]);
        double sum_hpsi2 = std::norm(hpsi[0]);

        for (int ig = 1; ig < ng; ++ig) {
            sum_psi_conj_hpsi += std::conj(psi[ig]) * hpsi[ig];
            sum_psi_conj += std::conj(psi[ig]);
            sum_psi += psi[ig];
            sum_hpsi += hpsi[ig];
            sum_psi2 += std::norm(psi[ig]);
            sum_hpsi2 += std::norm(hpsi[ig]);
        }

        double D = sum_psi2 * static_cast<double>(ng) - std::norm(sum_psi);
        std::complex<double> a_reg = (sum_psi_conj_hpsi * static_cast<double>(ng)
                                     - sum_psi_conj * sum_hpsi) / D;
        std::complex<double> b_reg = (sum_psi2 * sum_hpsi
                                     - sum_psi * sum_psi_conj_hpsi) / D;

        // Residual
        double sum_residual2 = 0.0;
        for (int ig = 0; ig < ng; ++ig) {
            auto r = hpsi[ig] - a_reg * psi[ig] - b_reg;
            sum_residual2 += std::norm(r);
        }
        double rel_residual = (sum_hpsi2 > 0.0)
            ? std::sqrt(sum_residual2 / sum_hpsi2) : 0.0;

        double delta = std::abs(std::real(e_rayleigh) - e_file);

        // Pass criteria (known: nonlocal pseudopotential needs fixing)
        double const tol = 2.0e-5;
        bool pass = (delta < tol) && (std::imag(e_rayleigh) < tol);

        results.push_back({ikpt, iband, e_file, std::real(e_rayleigh), delta,
                           std::imag(e_rayleigh), rel_residual, norm_check, pass});
    }

    // Summary table — printed after all H|ψ⟩ computations (and their HPSI_DEBUG
    // output) so the table is not interleaved with per-band diagnostics.
    std::println("");
    std::println("{}", std::string(100, '='));
    std::println("{:>4} {:>5}  {:>18}  {:>18}  {:>10}  {:>10}  {:>10}  {:>10}",
                 "ikpt", "iband", "E_file", "E_Rayleigh", "ΔE",
                 "Im(E)", "rel_res", "‖ψ‖²Ω");
    std::println("{}", std::string(100, '-'));
    int n_passed = 0;
    for (const auto& r : results) {
        if (r.pass) ++n_passed;
        std::println(
            "{:4d} {:5d}  {:>18.10f}  {:>18.10f}  {:>10.2e}  {:>10.2e}  {:>10.2e}  {:>10.6f}",
            r.ikpt, r.iband,
            r.e_file, r.e_rayleigh, r.delta,
            r.im_e, r.rel_res, r.norm_check);
    }
    std::println("{}", std::string(100, '-'));
    std::println("Passed {}/{} tests  (tol={:.0e} Ha)", n_passed, tests.size(), 2.0e-5);
    std::println("");

    if (n_passed == static_cast<int>(tests.size())) {
        std::println("All tests passed.");
        return 0;
    } else {
        // Note: failures are expected while the nonlocal pseudopotential
        // implementation is being debugged.
        std::println("Some tests did not meet tolerance.");
        std::println("(Expected: the nonlocal V_NL contribution needs debugging.)");
        return 1;
    }
}
