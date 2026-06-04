import io;
import std;


auto main() -> int {
    try {
        WG wg("test/data_io-local/OUT.WG");

        const auto& m = wg.meta;

        // Metadata checks
        if (m.n1 != 20)     throw std::runtime_error("n1 mismatch: expected 20");
        if (m.n2 != 20)     throw std::runtime_error("n2 mismatch: expected 20");
        if (m.n3 != 32)     throw std::runtime_error("n3 mismatch: expected 32");
        if (m.nband != 26)     throw std::runtime_error("nband mismatch: expected 26");
        if (m.mg_nx != 1042) throw std::runtime_error("mg_nx mismatch: expected 1042");
        if (m.nnode != 2)  throw std::runtime_error("nnode mismatch: expected 2");
        if (m.nkpt != 10)   throw std::runtime_error("nkpt mismatch: expected 10");
        if (m.is_SO != 0)   throw std::runtime_error("is_SO mismatch: expected 0");
        if (m.islda != 1)   throw std::runtime_error("islda mismatch: expected 1");
        if (std::abs(m.Ecut - 25.0) > 1e-12) {
            throw std::runtime_error("Ecut mismatch: expected 25.0");
        }

        // ng_tot_per_kpt checks
        if (m.ng_tot_per_kpt.size() != 10) {
            throw std::runtime_error("ng_tot_per_kpt size mismatch");
        }
        const std::vector<int> expected_ng = {1737, 1714, 1719, 1716, 1704,
                                              1718, 1729, 1710, 1702, 1708};
        for (int ikpt = 0; ikpt < 10; ++ikpt) {
            if (m.ng_tot_per_kpt[ikpt] != expected_ng[ikpt]) {
                throw std::runtime_error(
                    "ng_tot_per_kpt[" + std::to_string(ikpt) + "] mismatch: expected " +
                    std::to_string(expected_ng[ikpt]) + ", got " +
                    std::to_string(m.ng_tot_per_kpt[ikpt])
                );
            }
        }

        // Load band (kpt=0, band=0)
        const auto& wfc00 = wg.loadBand(0, 0);
        if (wg.current_ikpt() != 0) {
            throw std::runtime_error("current_ikpt should be 0 after loadBand(0,0)");
        }
        if (wg.current_iband() != 0) {
            throw std::runtime_error("current_iband should be 0 after loadBand(0,0)");
        }
        if (wfc00.up.size() != 1737) {
            throw std::runtime_error("wfc up size mismatch: expected 1737");
        }
        if (!wfc00.down.empty()) {
            throw std::runtime_error("wfc down should be empty when is_SO==0");
        }

        // Spot-check first few coefficients
        auto check_cplx = [](std::complex<double> a, std::complex<double> b, double tol, const std::string& msg) {
            if (std::abs(a - b) > tol) {
                throw std::runtime_error(msg + ": expected " + std::to_string(b.real()) +
                    "+" + std::to_string(b.imag()) + "i, got " +
                    std::to_string(a.real()) + "+" + std::to_string(a.imag()) + "i");
            }
        };

        check_cplx(wfc00.up[0],   std::complex<double>(5.192385287955403e-04,  5.399459972977638e-03), 1e-15, "wfc[0]");
        check_cplx(wfc00.up[1],   std::complex<double>(-6.130161928012967e-04, -6.359556224197149e-03), 1e-15, "wfc[1]");
        check_cplx(wfc00.up[2],   std::complex<double>(1.313827087869868e-04,  1.359764137305319e-03), 1e-15, "wfc[2]");
        check_cplx(wfc00.up[100], std::complex<double>(1.333363979938440e-05,  1.674602390266955e-04), 1e-15, "wfc[100]");

        // Cache test: loading same band again should return identical data
        // Note: loadBand returns a reference; subsequent calls invalidate prior references
        auto first_coeff_b0 = wfc00.up[0];
        const auto& wfc00_cached = wg.loadBand(0, 0);
        if (wfc00_cached.up.size() != 1737) {
            throw std::runtime_error("cached wfc size mismatch");
        }
        if (std::abs(wfc00_cached.up[0] - first_coeff_b0) > 1e-15) {
            throw std::runtime_error("cached wfc[0] mismatch");
        }

        // Load a different band at the same k-point
        const auto& wfc01 = wg.loadBand(0, 1);
        if (wg.current_iband() != 1) {
            throw std::runtime_error("current_iband should be 1 after loadBand(0,1)");
        }
        if (wfc01.up.size() != 1737) {
            throw std::runtime_error("band 1 wfc size mismatch");
        }
        // Band 0 and band 1 should have different coefficients
        if (std::abs(wfc01.up[0] - first_coeff_b0) < 1e-15) {
            throw std::runtime_error("band 0 and band 1 should have different coefficients");
        }

        // Load a band at a different k-point
        const auto& wfc10 = wg.loadBand(1, 0);
        if (wg.current_ikpt() != 1) {
            throw std::runtime_error("current_ikpt should be 1 after loadBand(1,0)");
        }
        if (wfc10.up.size() != 1714) {
            throw std::runtime_error("kpt 1 wfc size mismatch: expected 1714");
        }

        // Cache should still hold previously loaded entries (capacity=4 default)
        const auto& wfc00_reload = wg.loadBand(0, 0);
        if (wfc00_reload.up.size() != 1737) {
            throw std::runtime_error("kpt 0 band 0 reload size mismatch");
        }

        std::println("All WG tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
