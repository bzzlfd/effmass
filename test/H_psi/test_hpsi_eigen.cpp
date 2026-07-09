// Test: H|ψ⟩ = E|ψ⟩ eigenstate consistency.
//
// For every (k-point, band) pair across all data:
//   1. Load ψ from OUT.WG
//   2. Compute H|ψ⟩ via the Callable framework
//   3. Extract eigenvalue via Rayleigh quotient:  E = ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩
//   4. Compare with reference eigenvalue from OUT.EIGEN
//   5. Check wavefunction normalisation
//   6. Residual: ‖(H - ⟨H⟩)|ψ⟩‖ / ‖ψ⟩‖  —  vector-component eigenstate check
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

    // Residual tolerance ‖(H - ⟨H⟩)|ψ⟩‖ / ‖ψ‖  (cf. WG_ERROR convention)
    double constexpr TOL_RESIDUAL_EASY       = 1.0e-4;  // WG_ERROR easy
    double constexpr TOL_RESIDUAL_DIFFICULT  = 5.0e-5;  // WG_ERROR difficult (0.5E-4)
    // Active tolerance — change to EASY if test data was produced with easy convergence
    double constexpr TOL_RESIDUAL = TOL_RESIDUAL_DIFFICULT;
    double constexpr TOL_DELTA = 2.0e-5;  // eigenvalue match tolerance (Rayleigh vs file)

    int const test_kpts[] = {0, 1};
    int const test_bands_1[] = {nband - 1};  // highest band (hardest — best coverage per call)
    int n_test = static_cast<int>(std::size(test_kpts));
    int n_test_bands = static_cast<int>(std::size(test_bands_1));

    std::println("H|ψ⟩ eigenstate test  |  {} bands × {} k-points = {} pairs  |  volume = {:.2f} Bohr³\n",
                 n_test_bands, n_test, n_test_bands * n_test, volume);

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
    struct KptSummary { int ikpt; int ng; int tested; int passed; double max_delta; double max_im; double max_residual; double min_norm; double max_norm; };
    std::vector<KptSummary> kpt_summaries;
    kpt_summaries.reserve(static_cast<std::size_t>(n_test));

    // Detailed results for global summary
    struct Detail { int ikpt, iband; double e_file, e_rayleigh, delta, im_e, norm, residual; bool pass; };
    std::vector<Detail> details;
    details.reserve(static_cast<std::size_t>(n_test * n_test_bands));

    int total_tested = 0;
    int total_passed = 0;

    for (int ikpt : test_kpts) {
        auto op = h.at_k(ikpt);
        int ng = op.dim();

        KptSummary ks{ikpt, ng, 0, 0, 0.0, 0.0, 0.0, 1e10, 0.0};

        for (int iband : test_bands_1) {
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

            // Residual: ‖(H - ⟨H⟩)|ψ⟩‖ / ‖ψ‖
            double residual = 0.0;
            auto er = std::real(e_rayleigh);
            for (int ig = 0; ig < ng; ++ig) {
                auto diff = hpsi[ig] - er * psi[ig];
                residual += std::norm(diff);
            }
            residual = std::sqrt(residual / psi_norm2);

            // Pass criteria: eigenvalue match + eigenstate quality (residual)
            bool pass = (delta < TOL_DELTA) && (std::abs(im_e) < TOL_DELTA) && (residual < TOL_RESIDUAL);

            ks.tested++;
            if (pass) ks.passed++;
            ks.max_delta = std::max(ks.max_delta, delta);
            ks.max_im    = std::max(ks.max_im, std::abs(im_e));
            ks.max_residual = std::max(ks.max_residual, residual);
            ks.min_norm  = std::min(ks.min_norm, norm);
            ks.max_norm  = std::max(ks.max_norm, norm);

            details.push_back({ikpt, iband, e_file, std::real(e_rayleigh), delta, im_e, norm, residual, pass});
        }

        kpt_summaries.push_back(ks);
        total_tested += ks.tested;
        total_passed += ks.passed;
    }

    // ── Print per-k-point summary ──
    std::println("{:>5} {:>6}  {:>6} {:>6}  {:>10}  {:>10}  {:>10}  {:>10} {:>10}",
                 "ikpt", "ng", "tested", "passed", "max|ΔE|", "max|Im(E)|", "max‖R‖", "min‖ψ‖²Ω", "max‖ψ‖²Ω");
    std::println("{}", std::string(100, '-'));
    for (const auto& ks : kpt_summaries) {
        std::println("{:5d} {:6d}  {:6d} {:6d}  {:>10.2e}  {:>10.2e}  {:>10.2e}  {:>10.6f} {:>10.6f}",
                     ks.ikpt, ks.ng, ks.tested, ks.passed,
                     ks.max_delta, ks.max_im, ks.max_residual, ks.min_norm, ks.max_norm);
    }
    std::println("{}", std::string(100, '-'));
    std::println("{:>23} {:6d}  (ΔE/Im tol = {:.0e}, ‖R‖ tol = {:.0e} Ha)", "", total_passed, TOL_DELTA, TOL_RESIDUAL);
    std::println("");

    // ── If failures exist, print details ──
    int n_failed = total_tested - total_passed;
    if (n_failed > 0) {
        std::println("Failed pairs:");
        std::println("{:>4} {:>5}  {:>18}  {:>18}  {:>10}  {:>10}  {:>10}",
                     "ikpt", "iband", "E_file", "E_Rayleigh", "ΔE", "Im(E)", "‖R‖");
        std::println("{}", std::string(82, '-'));
        for (const auto& d : details) {
            if (!d.pass) {
                std::println("{:4d} {:5d}  {:>18.10f}  {:>18.10f}  {:>10.2e}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.e_file, d.e_rayleigh, d.delta, d.im_e, d.residual);
            }
        }
        std::println("{}", std::string(82, '-'));
    }

    // ── Residual summary (SCF) ──
    double sum_residual = 0.0, max_residual = 0.0;
    for (const auto& d : details) {
        sum_residual += d.residual;
        max_residual = std::max(max_residual, d.residual);
    }
    double mean_residual = sum_residual / static_cast<double>(details.size());
    std::println("Residual ‖(H - ⟨H⟩)|ψ⟩‖ / ‖ψ‖  ―  max = {:.2e}, mean = {:.2e}", max_residual, mean_residual);

    bool scf_ok = (n_failed == 0);
    std::println("data_scf: {}/{} passed  {}\n", total_passed, total_tested,
                 scf_ok ? "✓" : "✗  (see above)");

    // =========================================================================
    //  Test 2: data_nonscf (Si, OUT1.* k-point path along [100], IN.VR)
    //
    //  Non-SCF dataset: uses an input potential (IN.VR) shared between two
    //  k-point meshes (OUT1.* and OUT2.*).  Verifies that even a non-SCF
    //  potential produces a consistent Hamiltonian.
    // =========================================================================
    Hamiltonian h2("test/data_nonscf");
    h2.loadGKK("OUT1.GKK");
    h2.loadWG("OUT1.WG");
    h2.loadVR("IN.VR");
    h2.loadATOM("atom.config");
    h2.loadNCPP("Si.SG15.PBE.UPF");
    h2.loadEIGEN("OUT1.EIGEN");
    h2.finalize({ExtendedCheck::NCPPAtomCoverage},
                static_cast<std::uint64_t>(PSPFeature::Nonlocal));

    auto& gkk2 = h2.gkk();
    auto& wg2  = h2.wg();
    const auto& eigen2 = h2.eigen();
    int nband2 = eigen2.meta.nband;
    int nscf_nkpt = eigen2.meta.nkpt;
    int const nscf_test_bands[] = {nband2 - 1};  // highest band only (hardest)
    int const nscf_test_kpts[] = {0, nscf_nkpt - 1};  // first + last
    int nscf_n_test_bands = static_cast<int>(std::size(nscf_test_bands));
    int nscf_n_test_kpts  = static_cast<int>(std::size(nscf_test_kpts));

    int nscf_total_tested = 0, nscf_total_passed = 0;

    std::println("=== data_nonscf (Si, {} bands) ===", nband2);
    std::println("Testing {} k-points × {} bands = {} pairs\n",
                 nscf_n_test_kpts, nscf_n_test_bands, nscf_n_test_kpts * nscf_n_test_bands);

    struct NscfKptSummary { int ikpt; int ng; int tested; int passed; double max_delta; double max_im; double max_residual; double min_norm; double max_norm; };
    std::vector<NscfKptSummary> nscf_kpt_summaries;
    nscf_kpt_summaries.reserve(static_cast<std::size_t>(nscf_n_test_kpts));

    struct NscfDetail { int ikpt, iband; double e_file, e_rayleigh, delta, im_e, norm, residual; bool pass; };
    std::vector<NscfDetail> nscf_details;
    nscf_details.reserve(static_cast<std::size_t>(nscf_n_test_kpts * nscf_n_test_bands));

    for (int ikpt : nscf_test_kpts) {
        auto op2 = h2.at_k(ikpt);
        int ng2 = op2.dim();

        NscfKptSummary ks{ikpt, ng2, 0, 0, 0.0, 0.0, 0.0, 1e10, 0.0};

        for (int iband : nscf_test_bands) {
            auto psi = wg2.loadBand(ikpt, iband).up;

            double psi_norm2 = 0.0;
            for (int ig = 0; ig < ng2; ++ig)
                psi_norm2 += std::norm(psi[ig]);
            double norm = psi_norm2 * gkk2.meta.lattice.volume();

            std::vector<std::complex<double>> hpsi2(static_cast<std::size_t>(ng2), 0.0);
            op2(psi, hpsi2);

            auto e_rayleigh = rayleigh(psi, hpsi2);
            double e_file = eigen2[iband, ikpt];
            double delta = std::abs(std::real(e_rayleigh) - e_file);
            double im_e  = std::imag(e_rayleigh);

            double residual = 0.0;
            auto er = std::real(e_rayleigh);
            for (int ig = 0; ig < ng2; ++ig) {
                auto diff = hpsi2[ig] - er * psi[ig];
                residual += std::norm(diff);
            }
            residual = std::sqrt(residual / psi_norm2);

            bool pass = (delta < TOL_DELTA) && (std::abs(im_e) < TOL_DELTA) && (residual < TOL_RESIDUAL);

            ks.tested++;
            if (pass) ks.passed++;
            ks.max_delta = std::max(ks.max_delta, delta);
            ks.max_im    = std::max(ks.max_im, std::abs(im_e));
            ks.max_residual = std::max(ks.max_residual, residual);
            ks.min_norm  = std::min(ks.min_norm, norm);
            ks.max_norm  = std::max(ks.max_norm, norm);

            nscf_details.push_back({ikpt, iband, e_file, std::real(e_rayleigh), delta, im_e, norm, residual, pass});
        }

        nscf_kpt_summaries.push_back(ks);
        nscf_total_tested += ks.tested;
        nscf_total_passed += ks.passed;
    }

    // ── Print data_nonscf per-k-point summary ──
    std::println("{:>5} {:>6}  {:>6} {:>6}  {:>10}  {:>10}  {:>10}  {:>10} {:>10}",
                 "ikpt", "ng", "tested", "passed", "max|ΔE|", "max|Im(E)|", "max‖R‖", "min‖ψ‖²Ω", "max‖ψ‖²Ω");
    std::println("{}", std::string(100, '-'));
    for (const auto& ks : nscf_kpt_summaries) {
        std::println("{:5d} {:6d}  {:6d} {:6d}  {:>10.2e}  {:>10.2e}  {:>10.2e}  {:>10.6f} {:>10.6f}",
                     ks.ikpt, ks.ng, ks.tested, ks.passed,
                     ks.max_delta, ks.max_im, ks.max_residual, ks.min_norm, ks.max_norm);
    }
    std::println("{}", std::string(100, '-'));
    std::println("data_nonscf: {}/{} passed (ΔE/Im tol = {:.0e}, ‖R‖ tol = {:.0e} Ha)",
                 nscf_total_passed, nscf_total_tested, TOL_DELTA, TOL_RESIDUAL);

    int nscf_n_failed = nscf_total_tested - nscf_total_passed;
    bool nscf_ok = (nscf_n_failed == 0);

    // ── Print data_nonscf failures if any ──
    if (nscf_n_failed > 0) {
        std::println("data_nonscf FAILED pairs:");
        std::println("{:>4} {:>5}  {:>18}  {:>18}  {:>10}  {:>10}  {:>10}",
                     "ikpt", "iband", "E_file", "E_Rayleigh", "ΔE", "Im(E)", "‖R‖");
        std::println("{}", std::string(82, '-'));
        for (const auto& d : nscf_details) {
            if (!d.pass) {
                std::println("{:4d} {:5d}  {:>18.10f}  {:>18.10f}  {:>10.2e}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.e_file, d.e_rayleigh, d.delta, d.im_e, d.residual);
            }
        }
        std::println("{}", std::string(82, '-'));
    }

    // =========================================================================
    //  Test 3: data_nonscf (Si, OUT2.*, same IN.VR)
    //
    //  Same non-SCF dataset, but using OUT2.* (second k-point mesh) instead
    //  of OUT1.*.  Tests whether IN.VR is consistent with OUT2 even if it
    //  isn't with OUT1.
    // =========================================================================
    Hamiltonian h3("test/data_nonscf");
    h3.loadGKK("OUT2.GKK");
    h3.loadWG("OUT2.WG");
    h3.loadVR("IN.VR");
    h3.loadATOM("atom.config");
    h3.loadNCPP("Si.SG15.PBE.UPF");
    h3.loadEIGEN("OUT2.EIGEN");
    h3.finalize({ExtendedCheck::NCPPAtomCoverage},
                static_cast<std::uint64_t>(PSPFeature::Nonlocal));

    auto& gkk3   = h3.gkk();
    auto& wg3    = h3.wg();
    const auto& eigen3 = h3.eigen();
    int nband3 = eigen3.meta.nband;
    int nkpt3  = eigen3.meta.nkpt;
    int const out2_test_kpts[] = {0, 1};               // Γ + non-Γ
    int const out2_test_bands[] = {nband3 - 1};         // highest band only
    int out2_n_test_kpts  = static_cast<int>(std::size(out2_test_kpts));
    int out2_n_test_bands = static_cast<int>(std::size(out2_test_bands));

    int out2_tested = 0, out2_passed = 0;

    std::println("\n=== data_nonscf OUT2.* (Si, {} k-pts, {} bands) ===", nkpt3, nband3);
    std::println("Testing {} k-points × {} bands = {} pairs\n",
                 out2_n_test_kpts, out2_n_test_bands, out2_n_test_kpts * out2_n_test_bands);

    struct Out2KptSummary { int ikpt; int ng; int tested; int passed; double max_delta; double max_im; double max_residual; double min_norm; double max_norm; };
    std::vector<Out2KptSummary> out2_kpt_summaries;
    out2_kpt_summaries.reserve(static_cast<std::size_t>(out2_n_test_kpts));

    struct Out2Detail { int ikpt, iband; double e_file, e_rayleigh, delta, im_e, norm, residual; bool pass; };
    std::vector<Out2Detail> out2_details;
    out2_details.reserve(static_cast<std::size_t>(out2_n_test_kpts * out2_n_test_bands));

    for (int ikpt : out2_test_kpts) {
        auto op3 = h3.at_k(ikpt);
        int ng3 = op3.dim();

        Out2KptSummary ks{ikpt, ng3, 0, 0, 0.0, 0.0, 0.0, 1e10, 0.0};

        for (int iband : out2_test_bands) {
            auto psi = wg3.loadBand(ikpt, iband).up;

            double psi_norm2 = 0.0;
            for (int ig = 0; ig < ng3; ++ig)
                psi_norm2 += std::norm(psi[ig]);
            double norm = psi_norm2 * gkk3.meta.lattice.volume();

            std::vector<std::complex<double>> hpsi3(static_cast<std::size_t>(ng3), 0.0);
            op3(psi, hpsi3);

            auto e_rayleigh = rayleigh(psi, hpsi3);
            double e_file = eigen3[iband, ikpt];
            double delta = std::abs(std::real(e_rayleigh) - e_file);
            double im_e  = std::imag(e_rayleigh);

            double residual = 0.0;
            auto er = std::real(e_rayleigh);
            for (int ig = 0; ig < ng3; ++ig) {
                auto diff = hpsi3[ig] - er * psi[ig];
                residual += std::norm(diff);
            }
            residual = std::sqrt(residual / psi_norm2);

            bool pass = (delta < TOL_DELTA) && (std::abs(im_e) < TOL_DELTA) && (residual < TOL_RESIDUAL);

            ks.tested++;
            if (pass) ks.passed++;
            ks.max_delta = std::max(ks.max_delta, delta);
            ks.max_im    = std::max(ks.max_im, std::abs(im_e));
            ks.max_residual = std::max(ks.max_residual, residual);
            ks.min_norm  = std::min(ks.min_norm, norm);
            ks.max_norm  = std::max(ks.max_norm, norm);

            out2_details.push_back({ikpt, iband, e_file, std::real(e_rayleigh), delta, im_e, norm, residual, pass});
        }

        out2_kpt_summaries.push_back(ks);
        out2_tested += ks.tested;
        out2_passed += ks.passed;
    }

    // ── Print OUT2.* per-k-point summary ──
    std::println("{:>5} {:>6}  {:>6} {:>6}  {:>10}  {:>10}  {:>10}  {:>10} {:>10}",
                 "ikpt", "ng", "tested", "passed", "max|ΔE|", "max|Im(E)|", "max‖R‖", "min‖ψ‖²Ω", "max‖ψ‖²Ω");
    std::println("{}", std::string(100, '-'));
    for (const auto& ks : out2_kpt_summaries) {
        std::println("{:5d} {:6d}  {:6d} {:6d}  {:>10.2e}  {:>10.2e}  {:>10.2e}  {:>10.6f} {:>10.6f}",
                     ks.ikpt, ks.ng, ks.tested, ks.passed,
                     ks.max_delta, ks.max_im, ks.max_residual, ks.min_norm, ks.max_norm);
    }
    std::println("{}", std::string(100, '-'));
    std::println("OUT2.*: {}/{} passed (ΔE/Im tol = {:.0e}, ‖R‖ tol = {:.0e} Ha)",
                 out2_passed, out2_tested, TOL_DELTA, TOL_RESIDUAL);

    int out2_n_failed = out2_tested - out2_passed;
    bool out2_ok = (out2_n_failed == 0);

    // ── Print OUT2.* failures if any ──
    if (out2_n_failed > 0) {
        std::println("OUT2.* FAILED pairs:");
        std::println("{:>4} {:>5}  {:>18}  {:>18}  {:>10}  {:>10}  {:>10}",
                     "ikpt", "iband", "E_file", "E_Rayleigh", "ΔE", "Im(E)", "‖R‖");
        std::println("{}", std::string(82, '-'));
        for (const auto& d : out2_details) {
            if (!d.pass) {
                std::println("{:4d} {:5d}  {:>18.10f}  {:>18.10f}  {:>10.2e}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.e_file, d.e_rayleigh, d.delta, d.im_e, d.residual);
            }
        }
        std::println("{}", std::string(82, '-'));
    }

    // ── Combined result ──
    bool all_passed = scf_ok && nscf_ok && out2_ok;
    std::println("\nCombined:  data_scf {}/{}  |  data_nonscf(OUT1) {}/{}  |  data_nonscf(OUT2) {}/{}  |  {}",
                 total_passed, total_tested,
                 nscf_total_passed, nscf_total_tested,
                 out2_passed, out2_tested,
                 all_passed ? "ALL PASSED ✓" : "FAILURES ✗");

    return all_passed ? 0 : 1;
}
