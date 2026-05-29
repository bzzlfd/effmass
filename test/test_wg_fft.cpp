import std;
import math;
import io;

#define OCC_READER_IMPORT_STD
#include "io/occ_reader.hpp"


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
        std::println("=== WG FFT Normalization Test ===");

        WG wg("test/data_io-nonlocal/OUT.WG");
        GKK gkk("test/data_io-nonlocal/OUT.GKK");
        RHO rho("test/data_io-nonlocal/OUT.RHO");

        // consistency of FFT grid dimensions and lattice volume
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

        std::println("  grid: {} x {} x {}, volume = {:.6f} Bohr³", n1, n2, n3, vol);

        // Load EIGEN for k-point comparison
        EIGEN eigen("test/data_io-nonlocal/OUT.EIGEN");
        check(eigen.meta.nkpt == wg.meta.nkpt, "EIGEN/WG nkpt consistency");
        check(eigen.meta.nband == wg.meta.mx, "EIGEN/WG nband consistency");

        // Compare infer_k with EIGEN k-points.
        // EIGEN stores k in Cartesian (Bohr^-1). Convert to fractional via A·k/(2π).
        auto Alat = gkk.meta.lattice.A();
        constexpr double TWO_PI = 2.0 * std::numbers::pi;
        auto cartToFrac = [&](double kx, double ky, double kz) -> std::array<double, 3> {
            std::array<double, 3> f{};
            for (int d = 0; d < 3; ++d)
                f[d] = (Alat[d][0]*kx + Alat[d][1]*ky + Alat[d][2]*kz) / TWO_PI;
            return f;
        };

        std::println("  infer_k vs EIGEN k-point comparison:");
        for (int ikpt = 0; ikpt < wg.meta.nkpt; ++ikpt) {
            gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
            auto kv = gkk.loadKPoint(ikpt);
            auto ik = kv.kPoint;
            auto& ek = eigen.kpt_vec[ikpt];
            auto ef = cartToFrac(ek.x, ek.y, ek.z);
            double dk = std::sqrt(
                (ik[0]-ef[0])*(ik[0]-ef[0]) + (ik[1]-ef[1])*(ik[1]-ef[1]) + (ik[2]-ef[2])*(ik[2]-ef[2]));
            std::println("    kpt={:2d}: infer_k=({:.12f},{:.12f},{:.12f})  "
                         "eigen_frac=({:.12f},{:.12f},{:.12f})  dk={:.2e}",
                         ikpt, ik[0], ik[1], ik[2], ef[0], ef[1], ef[2], dk);
            check(dk < 1e-12,
                  std::format("kpt={}: infer_k matches EIGEN k-point", ikpt));
        }

        // Normalization check: sum|W|^2 * volume should be 1.0
        std::println("  Normalization (sum|W|^2 * volume == 1.0):");
        int n_norm_fail = 0;
        for (int ikpt = 0; ikpt < wg.meta.nkpt; ++ikpt) {
            gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
            auto kvecs2 = gkk.loadKPoint(ikpt);
            int ng2 = static_cast<int>(kvecs2.iG.size());

            for (int iband = 0; iband < wg.meta.mx; ++iband) {
                auto coeffs = wg.loadBand(ikpt, iband);
                check(static_cast<int>(coeffs.up.size()) == ng2,
                      std::format("kpt={} band={}: WG/GKK G-vector count match", ikpt, iband));

                double norm_g = computeGspaceNorm(coeffs.up, vol);
                if (!near(norm_g, 1.0)) {
                    std::println("    kpt={:2d} band={:2d}: FAIL norm = {:.6f}", ikpt, iband, norm_g);
                    ++n_norm_fail;
                }
            }
        }
        if (n_norm_fail > 0) {
            throw std::runtime_error(
                std::format("Normalization: {} bands deviate from 1.0", n_norm_fail));
        }
        std::println("  Normalization: all bands = 1.0");

        // Full FFT roundtrip for band 0 (reload kpt=0 to ensure valid references)
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto kvecs0 = gkk.loadKPoint(0);
        int ng0 = static_cast<int>(kvecs0.iG.size());

        auto coeffs = wg.loadBand(0, 0);
        double norm_g0 = computeGspaceNorm(coeffs.up, vol);

        // Build full FFT grid from sparse G-vector coefficients
        std::vector<std::complex<double>> grid(static_cast<std::size_t>(n123), 0.0);

        for (int ig = 0; ig < ng0; ++ig) {
            int i = kvecs0.iG[ig];
            int j = kvecs0.jG[ig];
            int k = kvecs0.kG[ig];

            int i_idx = ((i % n1) + n1) % n1;
            int j_idx = ((j % n2) + n2) % n2;
            int k_idx = ((k % n3) + n3) % n3;

            grid[static_cast<std::size_t>(i_idx) * n2 * n3
               + static_cast<std::size_t>(j_idx) * n3
               + static_cast<std::size_t>(k_idx)] = coeffs.up[ig];
        }

        // FFT G → R (unnormalized inverse DFT)
        // Physical wavefunction: ψ(r) = IFFT[W_g]
        // By Parseval: sum_r |ψ|² = N_tot * sum_g |W_g|²
        // So: sum_r |ψ|² * vol / N_tot = sum_g |W_g|² * vol = norm_g0
        // G-space and R-space normalization should be identical.
        FFT3D fft(n1, n2, n3);
        fft(grid, G2R);

        double sum_r = 0.0;
        for (auto& v : grid) sum_r += std::norm(v);
        double norm_r = sum_r * vol / static_cast<double>(n123);

        std::println("  G-space (band 0): sum|W_g|^2 * volume            = {:.6f}", norm_g0);
        std::println("  R-space (band 0): sum|W_r|^2 * volume/(n1*n2*n3) = {:.6f} (via G2R)", norm_r);
        check(near(norm_g0, norm_r), "G-space ↔ R-space Parseval consistency");

        std::println("\nAll tests passed!");
        return 0;

    } catch (const std::exception& e) {
        std::println("Error: {}", e.what());
        return 1;
    }
}
