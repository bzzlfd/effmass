import std;
import math;
import io;


auto check(bool cond, std::string_view msg) -> void {
    if (!cond) throw std::runtime_error(std::string("FAILED: ") + std::string(msg));
    std::println("  PASSED: {}", msg);
}


auto near(double a, double b, double eps = 1e-6) -> bool {
    return std::abs(a - b) < eps;
}


auto computeGspaceNorm(std::span<const std::complex<double>> wfc, double vol) -> double {
    double sum = 0.0;
    for (auto& v : wfc) sum += std::norm(v);
    return sum * vol;
}


auto main() -> int {
    try {
        // =========================================================================
        // 0. Open data files
        // =========================================================================
        WG   wg("test/data_io-nonlocal/OUT.WG");
        GKK  gkk("test/data_io-nonlocal/OUT.GKK");
        RHO  rho("test/data_io-nonlocal/OUT.RHO");
        EIGEN eigen("test/data_io-nonlocal/OUT.EIGEN");

        // --- grid & volume consistency ---
        check(wg.meta.n1 == gkk.meta.n1 && gkk.meta.n1 == rho.meta.n1, "n1 consistency");
        check(wg.meta.n2 == gkk.meta.n2 && gkk.meta.n2 == rho.meta.n2, "n2 consistency");
        check(wg.meta.n3 == gkk.meta.n3 && gkk.meta.n3 == rho.meta.n3, "n3 consistency");
        int n1 = rho.meta.n1, n2 = rho.meta.n2, n3 = rho.meta.n3;
        int n123 = n1 * n2 * n3;

        double vol_wg  = wg.meta.lattice.volume();
        double vol_gkk = gkk.meta.lattice.volume();
        double vol_rho = rho.lattice.volume();
        check(near(vol_wg, vol_gkk) && near(vol_gkk, vol_rho), "volume consistency");
        double vol = vol_rho;
        std::println("  grid: {} x {} x {}  volume = {:.6f} Bohr³", n1, n2, n3, vol);

        // --- metadata ---
        check(eigen.meta.nkpt == wg.meta.nkpt, "EIGEN/WG nkpt consistency");
        check(eigen.meta.nband == wg.meta.nband, "EIGEN/WG nband consistency");

        // =========================================================================
        // 1. GKK.inferCurrent_k  vs  EIGEN k-points
        // =========================================================================
        std::println("\n--- GKK.inferCurrent_k  vs  EIGEN ---");
        auto Alat = gkk.meta.lattice.A();
        constexpr double TWO_PI = 2.0 * std::numbers::pi;
        auto cartToFrac = [&](double kx, double ky, double kz) -> std::array<double, 3> {
            std::array<double, 3> f{};
            for (int d = 0; d < 3; ++d)
                f[d] = (Alat[d][0]*kx + Alat[d][1]*ky + Alat[d][2]*kz) / TWO_PI;
            return f;
        };
        for (int ikpt = 0; ikpt < wg.meta.nkpt; ++ikpt) {
            gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
            auto kv = gkk.loadKPoint(ikpt);
            auto ik = kv.kPoint;
            auto& ek = eigen.kpt_vec[ikpt];
            auto ef = cartToFrac(ek.x, ek.y, ek.z);
            double dk = std::sqrt(
                (ik[0]-ef[0])*(ik[0]-ef[0]) + (ik[1]-ef[1])*(ik[1]-ef[1]) + (ik[2]-ef[2])*(ik[2]-ef[2]));
            check(dk < 1e-12,
                  std::format("kpt={}: infer_k matches EIGEN k-point", ikpt));
        }

        // =========================================================================
        // 2. G-space normalization:  Σ|W_g|² · Ω = 1  for every band
        // =========================================================================
        std::println("\n--- G-space normalization  Σ|W_g|² · Ω == 1 ---");
        int n_norm_fail = 0;
        for (int ikpt = 0; ikpt < wg.meta.nkpt; ++ikpt) {
            gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
            auto& kvecs = gkk.loadKPoint(ikpt);
            int ng = static_cast<int>(kvecs.g_idx.size());

            for (int iband = 0; iband < wg.meta.nband; ++iband) {
                auto coeffs = wg.loadBand(ikpt, iband);
                check(static_cast<int>(coeffs.up.size()) == ng,
                      std::format("kpt={} band={}: WG/GKK G-vector count match", ikpt, iband));

                double norm_g = computeGspaceNorm(coeffs.up, vol);
                if (!near(norm_g, 1.0)) {
                    std::println("    kpt={:2d} band={:2d}: norm = {:.6f}", ikpt, iband, norm_g);
                    ++n_norm_fail;
                }
            }
        }
        if (n_norm_fail > 0)
            throw std::runtime_error(std::format("{} bands deviate from 1.0", n_norm_fail));

        // =========================================================================
        // 3. Parseval: G-space ↔ R-space  for one band
        // =========================================================================
        std::println("\n--- Parseval  G ↔ R (kpt=0 band=0) ---");
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto& kvecs0 = gkk.loadKPoint(0);
        int ng0 = static_cast<int>(kvecs0.g_idx.size());

        auto coeffs = wg.loadBand(0, 0);
        double norm_g0 = computeGspaceNorm(coeffs.up, vol);

        std::vector<std::complex<double>> grid(static_cast<std::size_t>(n123), 0.0);
        for (int ig = 0; ig < ng0; ++ig) {
            int i_idx = ((kvecs0.g_idx[ig].x % n1) + n1) % n1;
            int j_idx = ((kvecs0.g_idx[ig].y % n2) + n2) % n2;
            int k_idx = ((kvecs0.g_idx[ig].z % n3) + n3) % n3;
            grid[static_cast<std::size_t>(i_idx) * n2 * n3
               + static_cast<std::size_t>(j_idx) * n3
               + static_cast<std::size_t>(k_idx)] = coeffs.up[ig];
        }

        FFT3D fft(n1, n2, n3);
        fft(grid, G2R);

        double sum_r = 0.0;
        for (auto& v : grid) sum_r += std::norm(v);
        double norm_r = sum_r * vol / static_cast<double>(n123);

        std::println("  G-space: Σ|W_g|²·Ω              = {:.6f}", norm_g0);
        std::println("  R-space: Σ|ψ(r)|²·Ω/(n1·n2·n3) = {:.6f}", norm_r);
        check(near(norm_g0, norm_r), "G-space ↔ R-space Parseval consistency");

        // =========================================================================
        // 4. Single-band |ψ|²  vs  PWmat OUT.WG2RHO reference
        // =========================================================================

        // --- 4a.  Γ-point  (kpt=0, band=15)  →  OUT.WG2RHO_1_16 ---
        std::println("\n--- Single-band density: kpt=0 band=15  vs  OUT.WG2RHO_1_16 ---");
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto& kv_gamma = gkk.loadKPoint(0);
        int ng_gamma = static_cast<int>(kv_gamma.g_idx.size());

        auto coeffs_gamma = wg.loadBand(0, 15);
        check(static_cast<int>(coeffs_gamma.up.size()) == ng_gamma,
              "kpt=0 band=15: WG/GKK G-vector count match");

        std::vector<std::complex<double>> gr_gamma(static_cast<std::size_t>(n123), 0.0);
        for (int ig = 0; ig < ng_gamma; ++ig) {
            int i_idx = ((kv_gamma.g_idx[ig].x % n1) + n1) % n1;
            int j_idx = ((kv_gamma.g_idx[ig].y % n2) + n2) % n2;
            int k_idx = ((kv_gamma.g_idx[ig].z % n3) + n3) % n3;
            gr_gamma[static_cast<std::size_t>(i_idx) * n2 * n3
                   + static_cast<std::size_t>(j_idx) * n3
                   + static_cast<std::size_t>(k_idx)] = coeffs_gamma.up[ig];
        }
        FFT3D fft_gamma(n1, n2, n3);
        fft_gamma(gr_gamma, G2R);

        std::vector<double> rho_gamma(static_cast<std::size_t>(n123));
        for (int i = 0; i < n123; ++i) rho_gamma[i] = std::norm(gr_gamma[i]);

        RHO rho_gamma_ref("test/data_io-nonlocal/OUT.WG2RHO_1_16");
        check(rho_gamma_ref.meta.n1 == n1 && rho_gamma_ref.meta.n2 == n2 && rho_gamma_ref.meta.n3 == n3,
              "OUT.WG2RHO_1_16 grid match");
        check(rho_gamma_ref.meta.nstate == 1, "OUT.WG2RHO_1_16 nstate == 1");
        {
            double max_d = 0.0;
            for (int i = 0; i < n1; ++i)
                for (int j = 0; j < n2; ++j)
                    for (int k = 0; k < n3; ++k) {
                        auto idx = static_cast<std::size_t>(i) * n2 * n3
                                 + static_cast<std::size_t>(j) * n3 + k;
                        double d = std::abs(rho_gamma[idx] - rho_gamma_ref[i, j, k]);
                        if (d > max_d) max_d = d;
                    }
            std::println("  max |comp - ref| = {:.4e}", max_d);
            check(max_d < 1e-10, "|ψ|² matches OUT.WG2RHO_1_16");
        }

        // --- 4b.  non-Γ  (kpt=7, band=16)  →  OUT.WG2RHO_8_17 ---
        std::println("\n--- Single-band density: kpt=7 band=16  vs  OUT.WG2RHO_8_17 ---");
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto& kv_k7 = gkk.loadKPoint(7);
        int ng_k7 = static_cast<int>(kv_k7.g_idx.size());

        auto coeffs_k7 = wg.loadBand(7, 16);
        check(static_cast<int>(coeffs_k7.up.size()) == ng_k7,
              "kpt=7 band=16: WG/GKK G-vector count match");

        std::vector<std::complex<double>> gr_k7(static_cast<std::size_t>(n123), 0.0);
        for (int ig = 0; ig < ng_k7; ++ig) {
            int i_idx = ((kv_k7.g_idx[ig].x % n1) + n1) % n1;
            int j_idx = ((kv_k7.g_idx[ig].y % n2) + n2) % n2;
            int k_idx = ((kv_k7.g_idx[ig].z % n3) + n3) % n3;
            gr_k7[static_cast<std::size_t>(i_idx) * n2 * n3
                + static_cast<std::size_t>(j_idx) * n3
                + static_cast<std::size_t>(k_idx)] = coeffs_k7.up[ig];
        }
        FFT3D fft_k7(n1, n2, n3);
        fft_k7(gr_k7, G2R);

        std::vector<double> rho_k7(static_cast<std::size_t>(n123));
        for (int i = 0; i < n123; ++i) rho_k7[i] = std::norm(gr_k7[i]);

        RHO rho_k7_ref("test/data_io-nonlocal/OUT.WG2RHO_8_17");
        check(rho_k7_ref.meta.n1 == n1 && rho_k7_ref.meta.n2 == n2 && rho_k7_ref.meta.n3 == n3,
              "OUT.WG2RHO_8_17 grid match");
        check(rho_k7_ref.meta.nstate == 1, "OUT.WG2RHO_8_17 nstate == 1");
        {
            double max_d = 0.0;
            for (int i = 0; i < n1; ++i)
                for (int j = 0; j < n2; ++j)
                    for (int k = 0; k < n3; ++k) {
                        auto idx = static_cast<std::size_t>(i) * n2 * n3
                                 + static_cast<std::size_t>(j) * n3 + k;
                        double d = std::abs(rho_k7[idx] - rho_k7_ref[i, j, k]);
                        if (d > max_d) max_d = d;
                    }
            std::println("  max |comp - ref| = {:.4e}", max_d);
            check(max_d < 1e-10, "|ψ|² matches OUT.WG2RHO_8_17");
        }

        // =========================================================================
        // 5. Charge density reconstruction:  Σ occ·|ψ|²  vs  OUT.RHO
        // =========================================================================
        std::println("\n--- Charge density  Σ occ·|ψ|²  vs  OUT.RHO ---");

        OCC occ("test/data_io-nonlocal/OUT.OCC");
        check(occ.meta.nkpt == wg.meta.nkpt, "OCC/WG nkpt consistency");
        check(occ.meta.nband == wg.meta.nband, "OCC/WG nband consistency");

        std::vector<double>              rho_acc(static_cast<std::size_t>(n123), 0.0);
        std::vector<std::complex<double>> buf(static_cast<std::size_t>(n123));
        FFT3D fft_acc(n1, n2, n3);

        int nband = wg.meta.nband;
        double total_occ = 0.0;
        for (int ikpt = 0; ikpt < wg.meta.nkpt; ++ikpt) {
            gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
            auto& kv = gkk.loadKPoint(ikpt);
            int ng = static_cast<int>(kv.g_idx.size());

            for (int iband = 0; iband < nband; ++iband) {
                double occ_val = occ.occupation(iband, ikpt);
                if (occ_val == 0.0) continue;
                total_occ += occ_val;

                auto coeffs = wg.loadBand(ikpt, iband);
                check(static_cast<int>(coeffs.up.size()) == ng,
                      std::format("kpt={} band={}: WG/GKK G-vector count match", ikpt, iband));

                std::fill(buf.begin(), buf.end(), 0.0);
                for (int ig = 0; ig < ng; ++ig) {
                    int i_idx = ((kv.g_idx[ig].x % n1) + n1) % n1;
                    int j_idx = ((kv.g_idx[ig].y % n2) + n2) % n2;
                    int k_idx = ((kv.g_idx[ig].z % n3) + n3) % n3;
                    buf[static_cast<std::size_t>(i_idx) * n2 * n3
                      + static_cast<std::size_t>(j_idx) * n3
                      + static_cast<std::size_t>(k_idx)] = coeffs.up[ig];
                }
                fft_acc(buf, G2R);
                for (int ir = 0; ir < n123; ++ir)
                    rho_acc[ir] += occ_val * std::norm(buf[ir]);
            }
        }

        // --- integrate OUT.RHO for comparison ---
        double rho_int = 0.0;
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
                for (int k = 0; k < n3; ++k)
                    rho_int += rho[i, j, k];
        rho_int *= vol / static_cast<double>(n123);

        std::println("  Σ occ_val        = {:.6f}", total_occ);
        std::println("  ∫ OUT.RHO d³r    = {:.6f}", rho_int);

        // --- point-by-point error ---
        {
            double max_d = 0.0, sum_sq = 0.0, sum_abs_ref = 0.0;
            for (int i = 0; i < n1; ++i)
                for (int j = 0; j < n2; ++j)
                    for (int k = 0; k < n3; ++k) {
                        auto idx = static_cast<std::size_t>(i) * n2 * n3
                                 + static_cast<std::size_t>(j) * n3 + k;
                        double d = rho_acc[idx] - rho[i, j, k];
                        sum_sq += d * d;
                        sum_abs_ref += std::abs(rho[i, j, k]);
                        double ad = std::abs(d);
                        if (ad > max_d) max_d = ad;
                    }
            double rmse    = std::sqrt(sum_sq / static_cast<double>(n123));
            double avg_abs = sum_abs_ref / static_cast<double>(n123);
            std::println("  max  |comp - ref| = {:.4e}", max_d);
            std::println("  RMSE              = {:.4e}  ({:.2f}% of avg|ref|)",
                         rmse, rmse / avg_abs * 100.0);
            check(rmse < 1e-2, "Charge density matches OUT.RHO   RMSE < 1e-2");
        }

        std::println("\nAll tests passed!");
        return 0;

    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
