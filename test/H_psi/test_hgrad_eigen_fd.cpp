// Test: ⟨ψ_nk|∂H/∂k_α|ψ_nk⟩ = ∂ε_nk/∂k_α  (Hellmann-Feynman theorem)
//
// Compares the analytical dε/dk from Hamiltonian::Gradient::Callable against
// a central finite-difference approximation from OUT.EIGEN eigenvalues.
//
// Uses test/data_nonscf/ (OUT1.*: 7 k-points, Si diamond, 8 atoms).
//   kpts 0–4:  5 points along [100]  — standard HF vs FD validation
//   kpt 5:     (0, 0.15, 0)  along [010]  — symmetry test
//   kpt 6:     (0, 0, 0.15)  along [001]  — symmetry test
//
// For kpts 5 and 6, the analytical HF derivative in the y/z direction
// should equal FD_x at kpt 1 (symmetry-equivalent under cubic rotation).
//
// Bands tested (non-degenerate; most have non-zero slope):
//   1, 6, 13, 16, 17, 22    (0-based)
//   16 is the CBM (∂ε/∂k ≈ 0 near minimum) — included for completeness,
//   but the HF vs FD comparison is weaker since |FD_x| may be below TOL.
//
// Interior k-points (central finite-difference): 1, 2, 3  (0-based)
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)

import std;
import io;
import H_psi;

auto main() -> int {
    // =========================================================================
    // 1. Load data (step-by-step, since file names are non-standard)
    // =========================================================================
    Hamiltonian h("test/data_nonscf");
    h.loadGKK("OUT1.GKK");
    h.loadWG("OUT1.WG");
    h.loadVR("IN.VR");            // NB: local potential is named IN.VR here
    h.loadATOM("atom.config");
    h.loadNCPP("Si.SG15.PBE.UPF");
    h.loadEIGEN("OUT1.EIGEN");
    h.finalize({ExtendedCheck::NCPPAtomCoverage},
               static_cast<std::uint64_t>(PSPFeature::Nonlocal)
             | static_cast<std::uint64_t>(PSPFeature::DBetaq));

    const auto& eigen = h.eigen();
    auto& gkk  = h.gkk();          // needed for GKK data views
    auto& wg   = h.wg();           // loadBand is non-const (modifies cache)
    int nkpt  = eigen.meta.nkpt;    // 5
    int nband = eigen.meta.nband;   // 26
    double volume = gkk.meta.lattice.volume();

    // =========================================================================
    // 2. Constants
    // =========================================================================
    double constexpr TOL_FD_HF  = 2.0e-3;   // |HF_x - FD_x|  (Ha*Bohr)
    double constexpr TOL_ZERO   = 1.0e-4;   // |HF_y|, |HF_z|  (Ha*Bohr)
    double constexpr TOL_IMAG   = 1.0e-5;   // |Im(<dH/dk_x>)|  (Ha*Bohr)

    // =========================================================================
    // 3. Test bands — all non-degenerate, non-zero slope
    // =========================================================================
    // Degenerate groups (excluded): 2–5, 7–10, 11–12, 14–15, 18–21, 23–24.
    // 16 is the CBM — ∂ε/∂k ≈ 0 near the minimum, so the HF vs FD
    // comparison at kpt 2 (|FD_x| ≈ 0.0015 Ha·Bohr < TOL) cannot
    // reliably detect sign errors. Included for completeness.
    int const test_bands[] = {1, 6, 13, 16, 17, 22};

    // Interior k-points only (central FD needs both neighbours)
    int const test_kpts[] = {1, 2, 3};

    int n_bands = static_cast<int>(std::size(test_bands));
    int n_kpts  = static_cast<int>(std::size(test_kpts));

    // =========================================================================
    // 4. Result tracking
    // =========================================================================
    struct Detail {
        int ikpt, iband;
        double epsilon;         // ε_nk  (Ha)
        double delta_e;        // │E_Rayleigh − E_file│ (Ha)
        double de_x_fd;        // ∂ε/∂k_x from central FD  (Ha·Bohr)
        double de_x_hf;        // ∂ε/∂k_x from Hellmann–Feynman  (Ha·Bohr)
        double de_x_kin;       // kinetic-only contribution  (Ha·Bohr)
        double de_x_nl;        // nonlocal-only contribution  (Ha·Bohr)
        double need_nl;        // FD_x − kin_x  — what nl SHOULD be  (Ha·Bohr)
        double de_y_hf;        // ∂ε/∂k_y  (Ha·Bohr, should be ~0)
        double de_y_kin;       // kinetic y contribution  (Ha·Bohr)
        double de_z_hf;        // ∂ε/∂k_z  (Ha·Bohr, should be ~0)
        double de_z_kin;       // kinetic z contribution  (Ha·Bohr)
        double delta_x;        // │HF_x − FD_x│
        double im_x;           // Im(⟨∂H/∂k_x⟩)
        double norm_check;     // ⟨ψ|ψ⟩ × Ω  (should be ≈ 1)
        bool pass;
    };
    std::vector<Detail> details;
    details.reserve(static_cast<std::size_t>(n_kpts * n_bands));

    int total_tested = 0;
    int total_passed = 0;

    // =========================================================================
    // 5. Print header
    // =========================================================================
    std::println("Hamiltonian k-gradient (Hellmann–Feynman) test");
    std::println("  data   : test/data_nonscf/  (OUT1.*, 7 k-pts: 5 along [100] + 2 along [010]/[001])");
    std::println("  bands  :  {}  (non-degenerate; 16=CBM with |∂ε/∂k| < TOL)", "1  6  13  16  17  22");
    std::println("  k-pts  :  interior 1–3  (central FD)");
    std::println("  tol    :  │Δx│ < {:.0e}, │y/z│ < {:.0e}, │Im│ < {:.0e}  (Ha·Bohr)\n",
                 TOL_FD_HF, TOL_ZERO, TOL_IMAG);

    // =========================================================================
    // 7. K-vector diagnostic: check symmetry and coordinate system
    // =========================================================================
    std::println("--- K-vector diagnostic (GKK: ky, kz should vanish for [100] path) ---");
    for (int ikpt = 0; ikpt < nkpt; ++ikpt) {
        auto tmp_op = h.gradient().at_k(ikpt);
        int ng_k = tmp_op.dim();
        const auto& kdata = gkk.currentData();
        auto k_gkk = gkk.inferCurrent_k();

        // search for G=(0,0,0): find minimum kinetic energy
        double min_kin = 1e10;
        int min_ig = -1;
        for (int ig = 0; ig < ng_k; ++ig) {
            if (kdata.kinetic[ig] < min_kin) { min_kin = kdata.kinetic[ig]; min_ig = ig; }
        }
        std::println("  kpt {}: GKK_k=( {:.6f}, {:.6f}, {:.6f} ) frac  |  EIGEN_k=( {:.6f}, {:.6f}, {:.6f} ) Bohr⁻¹",
                     ikpt, k_gkk.x, k_gkk.y, k_gkk.z,
                     eigen.kpt_vec[ikpt].x, eigen.kpt_vec[ikpt].y, eigen.kpt_vec[ikpt].z);
        std::println("         minKin @G[{}]: {}  K=({:.4f}, {:.4f}, {:.4f})",
                     min_ig, min_kin, kdata.Kx[min_ig], kdata.Ky[min_ig], kdata.Kz[min_ig]);
        std::println("         First 8: kinetic  Kx        Ky        Kz");
        for (int ig = 0; ig < 8 && ig < ng_k; ++ig) {
            std::println("           G[{:2d}]: {:.6f}  {:.6f}  {:.6f}  {:.6f}",
                         ig, kdata.kinetic[ig], kdata.Kx[ig], kdata.Ky[ig], kdata.Kz[ig]);
        }
        // Check mirror-symmetry pairing: sum G_y²·ψ² over pairs (Gy > 0 and Gy < 0)
        // by classifying G-vectors
        double sum_pos = 0.0, sum_neg = 0.0;
        int npos = 0, nneg = 0;
        auto psi_tmp = wg.loadBand(ikpt, 1).up;
        for (int ig = 0; ig < ng_k; ++ig) {
            double gy = -kdata.Ky[ig];    // G_y (since -Ky = G+k, and k_y=0 so -Ky=G_y)
            double w  = std::norm(psi_tmp[ig]);
            if (gy > 1e-8)  { sum_pos += w * gy; npos++; }
            if (gy < -1e-8) { sum_neg += w * gy; nneg++; }
        }
        std::println("         band=1: ΣG_y>0|ψ|²·G_y = {:.6f} (n={}),  ΣG_y<0|ψ|²·G_y = {:.6f} (n={}),  net={:.6f}",
                     sum_pos, npos, sum_neg, nneg, sum_pos + sum_neg);
    }
    std::println();

    // =========================================================================
    // 8. Main test loop
    // =========================================================================
    for (int ikpt : test_kpts) {
        // Construct the gradient callable once per k-point, reuse across bands
        auto grad_op = h.gradient().at_k(ikpt);
        int ng = grad_op.dim();

        // Verify that GKK and WG agree on ng for this k-point
        int ng_wg = wg.meta.ng_tot_per_kpt[ikpt];
        if (ng != ng_wg) {
            std::println("ERROR: ng mismatch at kpt {}: grad_op.dim()={}, wg.ng={}", ikpt, ng, ng_wg);
            return 1;
        }

        for (int iband : test_bands) {
            // --- a. Load wavefunction ---
            auto psi = wg.loadBand(ikpt, iband).up;

            // --- a1. Direct kinetic sum from GKK data (cross-check) ---
            double kin_x_direct = 0.0, kin_y_direct = 0.0, kin_z_direct = 0.0;
            double psi_norm2_direct = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                double nrm2 = std::norm(psi[ig]);
                psi_norm2_direct += nrm2;
                kin_x_direct += nrm2 * (-gkk.currentData().Kx[ig]);
                kin_y_direct += nrm2 * (-gkk.currentData().Ky[ig]);
                kin_z_direct += nrm2 * (-gkk.currentData().Kz[ig]);
            }
            // Cross-check on first iteration only
            if (ikpt == test_kpts[0] && iband == test_bands[0]) {
                double cx = kin_y_direct / psi_norm2_direct;
                std::println("  [CC] kpt={}: kin_x={:.6f} kin_y={:.6f} kin_z={:.6f} | ψ_norm²={:.9f}",
                             ikpt, kin_x_direct/psi_norm2_direct, cx,
                             kin_z_direct/psi_norm2_direct, psi_norm2_direct);
            }

            // --- b. Compute ∂H/∂k_α|ψ⟩ for all three Cartesian directions ---
            std::vector<std::complex<double>> dx(static_cast<std::size_t>(ng), 0.0);
            std::vector<std::complex<double>> dy(static_cast<std::size_t>(ng), 0.0);
            std::vector<std::complex<double>> dz(static_cast<std::size_t>(ng), 0.0);
            grad_op(psi, GradHPsi{dx, dy, dz});

            // --- b1. Rayleigh check: verify H|ψ⟩ = ε|ψ⟩ for this dataset ---
            double e_file = eigen[iband, ikpt];
            double e_rayleigh = e_file;  // default: assume exact match
            if (ikpt == test_kpts[0] && iband == test_bands[0]) {
                std::vector<std::complex<double>> hpsi(static_cast<std::size_t>(ng), 0.0);
                h.at_k(ikpt)(psi, hpsi);
                std::complex<double> rayleigh_sum{0.0, 0.0};
                double rayleigh_norm = 0.0;
                for (int ig = 0; ig < ng; ++ig) {
                    rayleigh_sum += std::conj(psi[ig]) * hpsi[ig];
                    rayleigh_norm += std::norm(psi[ig]);
                }
                e_rayleigh = std::real(rayleigh_sum) / rayleigh_norm;

                std::println("  [Rl] E_file={:.8f} E_Rayleigh={:.8f} ΔE={:.2e}",
                             e_file, e_rayleigh, std::abs(e_rayleigh - e_file));
                // Check dbetaq interpolation
                const auto& ncp = h.ncpp(14);
                if (ncp.meta.l_max >= 0) {
                    const auto& dbt = h.dbetaqTables(14);
                    double q = gkk.currentData().q[10];  // some non-zero q
                    std::println("  [nl] dβ/dq l=0 ib=0: q=0 dβ={:.4e}  q={:.4f} dβ={:.4e}",
                                 dbt.interpolate(0, 0, 0.0), q, dbt.interpolate(0, 0, q));
                    std::println("  [nl] dβ/dq l=1 ib=0: q=0 dβ={:.4e}  q={:.4f} dβ={:.4e}",
                                 dbt.interpolate(1, 0, 0.0), q, dbt.interpolate(1, 0, q));
                    std::println("  [nl] β(q)    l=0 ib=0: q=0 β={:.4e}  q={:.4f} β={:.4e}",
                                 h.betaqTables(14).interpolate(0, 0, 0.0), q,
                                 h.betaqTables(14).interpolate(0, 0, q));
                }
            }

            // --- c. Hellmann–Feynman:  ∂ε/∂k_α = ⟨ψ|∂H/∂k_α|ψ⟩ / ⟨ψ|ψ⟩ ---
            double psi_norm2 = 0.0;
            std::complex<double> sum_x{0.0, 0.0}, sum_y{0.0, 0.0}, sum_z{0.0, 0.0};
            // Kinetic-only contribution (for debug) — computed directly from
            // the stored K vectors:  ∂T/∂k_α = (-K_α)  since  K = -(G+k).
            // So  ⟨ψ|∂T/∂k_α|ψ⟩ = Σ |ψ|² · (-K_α) / Σ|ψ|²  (all real).
            double kin_x = 0.0, kin_y = 0.0, kin_z = 0.0;
            for (int ig = 0; ig < ng; ++ig) {
                double nrm2 = std::norm(psi[ig]);
                psi_norm2 += nrm2;
                sum_x += std::conj(psi[ig]) * dx[ig];
                sum_y += std::conj(psi[ig]) * dy[ig];
                sum_z += std::conj(psi[ig]) * dz[ig];
                kin_x += nrm2 * (-gkk.currentData().Kx[ig]);
                kin_y += nrm2 * (-gkk.currentData().Ky[ig]);
                kin_z += nrm2 * (-gkk.currentData().Kz[ig]);
            }

            double de_x_hf = std::real(sum_x) / psi_norm2;
            double de_y_hf = std::real(sum_y) / psi_norm2;
            double de_z_hf = std::real(sum_z) / psi_norm2;
            double im_x     = std::imag(sum_x) / psi_norm2;
            double de_x_kin = kin_x / psi_norm2;
            double de_y_kin = kin_y / psi_norm2;
            double de_z_kin = kin_z / psi_norm2;
            double psi_norm_check = psi_norm2 * volume;

            // --- d. Finite-difference for x-direction ---
            // Central:  ∂ε/∂k_x  ≈  [ε(k_{i+1}) − ε(k_{i-1})] / [k_x(i+1) − k_x(i-1)]
            double e_lo = eigen[iband, ikpt - 1];
            double e_hi = eigen[iband, ikpt + 1];
            double k_lo = eigen.kpt_vec[ikpt - 1].x;  // Cartesian x-coordinate (Bohr⁻¹)
            double k_hi = eigen.kpt_vec[ikpt + 1].x;
            double de_x_fd = (e_hi - e_lo) / (k_hi - k_lo);

            double epsilon = eigen[iband, ikpt];
            double delta_x = std::abs(de_x_hf - de_x_fd);

            // --- e. Pass / fail ---
            bool pass = (delta_x < TOL_FD_HF)
                     && (std::abs(de_y_hf) < TOL_ZERO)
                     && (std::abs(de_z_hf) < TOL_ZERO)
                     && (std::abs(im_x)   < TOL_IMAG);

            ++total_tested;
            if (pass) ++total_passed;

            double de_x_nl = de_x_hf - de_x_kin;
            double need_nl = de_x_fd - de_x_kin;   // what nl_x SHOULD be for HF=FD
            double de = std::abs(e_rayleigh - e_file);
            details.push_back({ikpt, iband, epsilon, de,
                               de_x_fd, de_x_hf, de_x_kin, de_x_nl, need_nl,
                               de_y_hf, de_y_kin, de_z_hf, de_z_kin,
                               delta_x, im_x, psi_norm_check, pass});
        }
    }

    // =========================================================================
    // 7. Print results table
    // =========================================================================
    std::println("{:>4} {:>5}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {}",
                 "kpt", "band", "ε(Ha)", "FD_x", "HF_x",
                 "kin_x", "nl_x", "need_nl", "HF_y", "HF_z",
                 "│Δx│", "│Im│", "pass");
    std::println("{}", std::string(160, '-'));

    for (const auto& d : details) {
        std::println("{:4d} {:5d}  {:>10.6f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.1e}  {:>10.1e}  {}",
                     d.ikpt, d.iband, d.epsilon,
                     d.de_x_fd, d.de_x_hf, d.de_x_kin, d.de_x_nl, d.need_nl,
                     d.de_y_hf, d.de_z_hf,
                     d.delta_x, d.im_x,
                     d.pass ? "✓" : "✗");
    }

    // =========================================================================
    // 9. Symmetry test: y/z k-points vs x-direction reference
    // =========================================================================
    //
    // kpt 5 = (0, 0.15, 0) and kpt 6 = (0, 0, 0.15) are symmetry-equivalent
    // to kpt 1 = (0.15, 0, 0) under cubic rotations of the diamond lattice.
    // Therefore:
    //   ∂ε/∂k_y at kpt 5  =  ∂ε/∂k_x at kpt 1
    //   ∂ε/∂k_z at kpt 6  =  ∂ε/∂k_x at kpt 1
    //
    // We validate by comparing the analytical HF derivative at kpts 5/6
    // against the central-FD reference at kpt 1.
    // =========================================================================
    std::println("\n--- Symmetry test: y/z k-points vs kpt 1 (x-direction) ---\n");

    // Pre-compute FD_x at kpt 1 for each test band (central FD reference)
    // Uses the same formula as the main test loop.
    std::vector<double> fd_x_kpt1(static_cast<std::size_t>(n_bands));
    for (int i = 0; i < n_bands; ++i) {
        int iband = test_bands[i];
        double e_lo  = eigen[iband, 0];
        double e_hi  = eigen[iband, 2];
        double k_lo  = eigen.kpt_vec[0].x;
        double k_hi  = eigen.kpt_vec[2].x;
        fd_x_kpt1[static_cast<std::size_t>(i)] = (e_hi - e_lo) / (k_hi - k_lo);
    }

    struct SymmDetail {
        int ikpt, iband;
        double epsilon;         // ε_nk  (Ha)
        double fd_x_ref;        // FD_x at kpt 1 (reference)  (Ha·Bohr)
        double hf_dir;          // HF_y (kpt 5) or HF_z (kpt 6)  (Ha·Bohr)
        double cross1, cross2;  // perpendicular components  (Ha·Bohr)
        double delta_sym;       // |hf_dir - fd_x_ref|
        double im_dir;          // Im(⟨∂H/∂k_dir⟩)  (Ha·Bohr)
        double psi_norm_check;  // ⟨ψ|ψ⟩ × Ω
        bool pass;
    };
    std::vector<SymmDetail> symm_details;
    symm_details.reserve(static_cast<std::size_t>(2 * n_bands));
    int symm_tested = 0, symm_passed = 0;

    std::println("{:>4} {:>5}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {}",
                 "kpt", "band", "ε(Ha)", "FD_x(ref)", "HF_dir",
                 "cross1", "cross2", "│Δsym│", "│Im│", "pass");
    std::println("{}", std::string(140, '-'));

    for (int ikpt : {5, 6}) {
        auto grad_op = h.gradient().at_k(ikpt);
        int ng = grad_op.dim();

        // kpt 5 → dir=1 (y);  kpt 6 → dir=2 (z)
        int dir = (ikpt == 5) ? 1 : 2;

        for (int i = 0; i < n_bands; ++i) {
            int iband = test_bands[i];
            auto psi = wg.loadBand(ikpt, iband).up;

            std::vector<std::complex<double>> dx(static_cast<std::size_t>(ng), 0.0);
            std::vector<std::complex<double>> dy(static_cast<std::size_t>(ng), 0.0);
            std::vector<std::complex<double>> dz(static_cast<std::size_t>(ng), 0.0);
            grad_op(psi, GradHPsi{dx, dy, dz});

            double psi_norm2 = 0.0;
            std::complex<double> sum_x{0.0, 0.0}, sum_y{0.0, 0.0}, sum_z{0.0, 0.0};
            for (int ig = 0; ig < ng; ++ig) {
                double nrm2 = std::norm(psi[ig]);
                psi_norm2 += nrm2;
                sum_x += std::conj(psi[ig]) * dx[ig];
                sum_y += std::conj(psi[ig]) * dy[ig];
                sum_z += std::conj(psi[ig]) * dz[ig];
            }

            double de[3] = {
                std::real(sum_x) / psi_norm2,
                std::real(sum_y) / psi_norm2,
                std::real(sum_z) / psi_norm2
            };
            std::complex<double> dir_sum = (dir == 0) ? sum_x : ((dir == 1) ? sum_y : sum_z);
            double im_dir = std::imag(dir_sum) / psi_norm2;

            double hf_dir    = de[dir];                           // primary: along kpt direction
            double cross1    = de[(dir == 0) ? 1 : 0];            // one perpendicular component
            double cross2    = de[(dir == 2) ? 1 : 2];            // the other perpendicular
            double fd_x_ref_val = fd_x_kpt1[static_cast<std::size_t>(i)];
            double delta_sym = std::abs(hf_dir - fd_x_ref_val);
            double epsilon   = eigen[iband, ikpt];
            double psi_norm_check = psi_norm2 * volume;

            bool pass = (delta_sym < TOL_FD_HF)
                     && (std::abs(cross1) < TOL_ZERO)
                     && (std::abs(cross2) < TOL_ZERO)
                     && (std::abs(im_dir) < TOL_IMAG);

            ++symm_tested;
            if (pass) ++symm_passed;

            symm_details.push_back({ikpt, iband, epsilon, fd_x_ref_val, hf_dir,
                                    cross1, cross2, delta_sym, im_dir, psi_norm_check, pass});
        }
    }

    // Print symmetry results table
    for (const auto& d : symm_details) {
        std::println("{:4d} {:5d}  {:>10.6f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.4f}  {:>10.1e}  {:>10.1e}  {}",
                     d.ikpt, d.iband, d.epsilon,
                     d.fd_x_ref, d.hf_dir,
                     d.cross1, d.cross2,
                     d.delta_sym, d.im_dir,
                     d.pass ? "✓" : "✗");
    }

    // Symmetry summary
    double max_delta_sym = 0.0, max_cross_sym = 0.0, max_im_sym = 0.0;
    for (const auto& d : symm_details) {
        max_delta_sym  = std::max(max_delta_sym, d.delta_sym);
        max_cross_sym  = std::max({max_cross_sym, std::abs(d.cross1), std::abs(d.cross2)});
        max_im_sym     = std::max(max_im_sym, std::abs(d.im_dir));
    }

    std::println("{}", std::string(125, '-'));
    std::println("Symmetry: {}/{} passed  |  max│Δsym│ = {:.2e}, max│cross│ = {:.2e}, max│Im│ = {:.2e}",
                 symm_passed, symm_tested, max_delta_sym, max_cross_sym, max_im_sym);

    int n_symm_failed = symm_tested - symm_passed;
    if (n_symm_failed > 0) {
        std::println("\nSymmetry FAILED cases:");
        std::println("{:>4} {:>5}  {:>10}  {:>10}  {:>10}",
                     "kpt", "band", "│Δsym│", "│cross│", "│Im│");
        for (const auto& d : symm_details) {
            if (!d.pass) {
                double cross = std::max(std::abs(d.cross1), std::abs(d.cross2));
                std::println("{:4d} {:5d}  {:>10.2e}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.delta_sym, cross, d.im_dir);
            }
        }
        return 1;
    }

    // =========================================================================
    // 8. Summary
    // =========================================================================
    std::println("{}", std::string(125, '-'));
    int n_failed = total_tested - total_passed;

    double max_delta = 0.0, max_yz = 0.0, max_im = 0.0, sum_delta = 0.0;
    for (const auto& d : details) {
        max_delta = std::max(max_delta, d.delta_x);
        max_yz    = std::max({max_yz, std::abs(d.de_y_hf), std::abs(d.de_z_hf)});
        max_im    = std::max(max_im, std::abs(d.im_x));
        sum_delta += d.delta_x;
    }
    double mean_delta = sum_delta / static_cast<double>(details.size());

    std::println("{:>4} / {:>4} passed  |  max │Δx│ = {:.2e}, mean │Δx│ = {:.2e}",
                 total_passed, total_tested, max_delta, mean_delta);
    std::println("max │y/z│ = {:.2e}  |  max │Im│ = {:.2e}", max_yz, max_im);

    bool all_fd_ok = (n_failed == 0);
    if (n_failed > 0) {
        std::println("\nFAILED cases:");
        std::println("{:>4} {:>5}  {:>10}  {:>10}  {:>10}",
                     "kpt", "band", "│Δx│", "│y/z│", "│Im│");
        for (const auto& d : details) {
            if (!d.pass) {
                double yz = std::max(std::abs(d.de_y_hf), std::abs(d.de_z_hf));
                std::println("{:4d} {:5d}  {:>10.2e}  {:>10.2e}  {:>10.2e}",
                             d.ikpt, d.iband, d.delta_x, yz, d.im_x);
            }
        }
    }

    bool all_ok = all_fd_ok && (n_symm_failed == 0);
    std::println("\nCombined:  x-FD {}/{}  |  symmetry {}/{}  |  {}",
                 total_passed, total_tested,
                 symm_passed, symm_tested,
                 all_ok ? "ALL PASSED ✓" : "FAILURES ✗");

    return all_ok ? 0 : 1;
}
