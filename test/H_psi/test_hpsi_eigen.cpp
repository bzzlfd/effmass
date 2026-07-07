// Test: H|ψ⟩ = E|ψ⟩ eigenstate consistency.
//
// For every (k-point, band) pair across all data:
//   1. Load ψ from OUT.WG
//   2. Compute H|ψ⟩ via the Callable framework
//   3. Extract eigenvalue via Rayleigh quotient:  E = ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩
//   4. Compare with reference eigenvalue from OUT.EIGEN
//   5. Check wavefunction normalisation
//
// Performance note: Callable is constructed once per k-point and reused
// across bands, avoiding redundant Ylm construction and GKK data reloads.
// Tests k-points 0 and 1 (52 pairs total).
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)
// Uses test data under test/data_scf/

import std;
import io;
import H_psi;

auto main() -> int {
    Hamiltonian h("test/data_scf");
    h.loadFromDirectory();

    auto& gkk   = h.gkk();
    auto& wg    = h.wg();
    const auto& eigen = h.eigen();
    double volume = gkk.meta.lattice.volume();
    int nband = wg.meta.nband;

    int const test_kpts[] = {0, 1};
    int n_test = static_cast<int>(std::size(test_kpts));

    std::println("H|ψ⟩ eigenstate test  |  {} bands × {} k-points = {} pairs  |  volume = {:.2f} Bohr³\n",
                 nband, n_test, nband * n_test, volume);

    // Rayleigh quotient:  E = ⟨ψ|H|ψ⟩ / ⟨ψ|ψ⟩
    struct RayleighSums { std::complex<double> num; double den; };
    auto rayleigh = [](auto psi, auto hpsi) -> std::complex<double> {
        RayleighSums s{0.0, 0.0};
        for (std::size_t i = 0; i < psi.size(); ++i) {
            s.num += std::conj(psi[i]) * hpsi[i];
            s.den += std::norm(psi[i]);
        }
        return s.num / s.den;
    };

    // Per-k-point summary
    struct KptSummary { int ikpt; int ng; int tested; int passed; double max_delta; double max_im; double min_norm; double max_norm; };
    std::vector<KptSummary> kpt_summaries;
    kpt_summaries.reserve(static_cast<std::size_t>(n_test));

    // Detailed results for global summary
    struct Detail { int ikpt, iband; double e_file, e_rayleigh, delta, im_e, norm; bool pass; };
    std::vector<Detail> details;
    details.reserve(static_cast<std::size_t>(n_test * nband));

    int total_tested = 0;
    int total_passed = 0;

    for (int ikpt : test_kpts) {
        auto op = h.at_k(ikpt);
        int ng = op.dim();

        KptSummary ks{ikpt, ng, 0, 0, 0.0, 0.0, 1e10, 0.0};

        for (int iband = 0; iband < nband; ++iband) {
            auto psi = wg.loadBand(ikpt, iband).up;

            // Wavefunction normalisation:  ⟨ψ|ψ⟩ × Ω  should ≈ 1
            double psi_norm2 = 0.0;
            for (int ig = 0; ig < ng; ++ig)
                psi_norm2 += std::norm(psi[ig]);
            double norm = psi_norm2 * volume;

            // H|ψ⟩
            std::vector<std::complex<double>> hpsi(static_cast<std::size_t>(ng), 0.0);
            op(psi, hpsi);

            // Rayleigh quotient → eigenvalue estimate
            auto e_rayleigh = rayleigh(psi, hpsi);
            double e_file = eigen[iband, ikpt];
            double delta = std::abs(std::real(e_rayleigh) - e_file);
            double im_e  = std::imag(e_rayleigh);

            // Pass criteria
            double const tol = 2.0e-5;
            bool pass = (delta < tol) && (std::abs(im_e) < tol);

            ks.tested++;
            if (pass) ks.passed++;
            ks.max_delta = std::max(ks.max_delta, delta);
            ks.max_im    = std::max(ks.max_im, std::abs(im_e));
            ks.min_norm  = std::min(ks.min_norm, norm);
            ks.max_norm  = std::max(ks.max_norm, norm);

            details.push_back({ikpt, iband, e_file, std::real(e_rayleigh), delta, im_e, norm, pass});
        }

        kpt_summaries.push_back(ks);
        total_tested += ks.tested;
        total_passed += ks.passed;
    }

    // ── Print per-k-point summary ──
    std::println("{:>5} {:>6}  {:>6} {:>6}  {:>10}  {:>10}  {:>10} {:>10}",
                 "ikpt", "ng", "tested", "passed", "max|ΔE|", "max|Im(E)|", "min‖ψ‖²Ω", "max‖ψ‖²Ω");
    std::println("{}", std::string(90, '-'));
    for (const auto& ks : kpt_summaries) {
        std::println("{:5d} {:6d}  {:6d} {:6d}  {:>10.2e}  {:>10.2e}  {:>10.6f} {:>10.6f}",
                     ks.ikpt, ks.ng, ks.tested, ks.passed,
                     ks.max_delta, ks.max_im, ks.min_norm, ks.max_norm);
    }
    std::println("{}", std::string(90, '-'));
    std::println("{:>23} {:6d}  (tol = {:.0e} Ha)", "", total_passed, 2.0e-5);
    std::println("");

    // ── If failures exist, print details ──
    int n_failed = total_tested - total_passed;
    if (n_failed > 0) {
        std::println("Failed pairs:");
        std::println("{:>4} {:>5}  {:>18}  {:>18}  {:>10}  {:>10}",
                     "ikpt", "iband", "E_file", "E_Rayleigh", "ΔE", "Im(E)");
        std::println("{}", std::string(70, '-'));
        for (const auto& d : details) {
            if (!d.pass) {
                std::println("{:4d} {:5d}  {:>18.10f}  {:>18.10f}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.e_file, d.e_rayleigh, d.delta, d.im_e);
            }
        }
        std::println("{}", std::string(70, '-'));
        return 1;
    }

    std::println("All {} tests passed.", total_tested);
    return 0;
}
