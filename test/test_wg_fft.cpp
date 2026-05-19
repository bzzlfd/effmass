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
        std::println("=== WG FFT Normalization Test ===");

        WG wg("test/test_io-local/OUT.WG");
        GKK gkk("test/test_io-local/OUT.GKK");
        RHO rho("test/test_io-local/OUT.RHO");

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

        // load GKK with Integer view to get G-vector indices (iG,jG,kG)
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto kvecs = gkk.loadKPoint(0);
        int ng = static_cast<int>(kvecs.iG.size());

        // Print normalization of first few bands at kpt=0
        std::println("  Normalization check (kpt=0):");
        for (int iband = 0; iband < std::min(5, wg.meta.mx); ++iband) {
            auto coeffs = wg.loadBand(0, iband);
            check(static_cast<int>(coeffs.up.size()) == ng,
                  std::format("kpt=0 band={}: WG/GKK G-vector count match", iband));
            double norm_g = computeGspaceNorm(coeffs.up, vol);
            std::println("    band {}: sum|W|^2 * volume = {:.6f}", iband, norm_g);
        }

        // Full FFT roundtrip for band 0
        auto coeffs = wg.loadBand(0, 0);
        double norm_g0 = computeGspaceNorm(coeffs.up, vol);

        // Build full FFT grid from sparse G-vector coefficients
        std::vector<std::complex<double>> grid(static_cast<std::size_t>(n123), 0.0);

        for (int ig = 0; ig < ng; ++ig) {
            int i = kvecs.iG[ig];
            int j = kvecs.jG[ig];
            int k = kvecs.kG[ig];

            int i_idx = ((i % n1) + n1) % n1;
            int j_idx = ((j % n2) + n2) % n2;
            int k_idx = ((k % n3) + n3) % n3;

            grid[static_cast<std::size_t>(i_idx) * n2 * n3
               + static_cast<std::size_t>(j_idx) * n3
               + static_cast<std::size_t>(k_idx)] = coeffs.up[ig];
        }

        // FFT G → R (backward FFT with 1/N scaling)
        // Physical wavefunction: ψ(r) = N_tot * IFFT[W_g]
        // By Parseval: sum_r |ψ|² = N_tot * sum_g |W_g|²
        // So: sum_r |ψ|² * vol / N_tot = sum_g |W_g|² * vol = norm_g0
        // G-space and R-space normalization should be identical.
        FFT3D fft(n1, n2, n3);
        fft(grid, G2R);

        double sum_r = 0.0;
        for (auto& v : grid) sum_r += std::norm(v);
        double norm_r = sum_r * static_cast<double>(n123) * vol;

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
