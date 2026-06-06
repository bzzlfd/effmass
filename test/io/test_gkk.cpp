import io;
import std;


auto main() -> int {
    try {
        GKK gkk("test/data_io-local/OUT.GKK");

        const auto& m = gkk.meta;

        // Metadata checks
        if (m.n1 != 20)     throw std::runtime_error("n1 mismatch: expected 20");
        if (m.n2 != 20)     throw std::runtime_error("n2 mismatch: expected 20");
        if (m.n3 != 32)     throw std::runtime_error("n3 mismatch: expected 32");
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

        // Default view is Cartesian
        if (gkk.currentView() != KVecsView::Cartesian) {
            throw std::runtime_error("default view should be Cartesian");
        }

        // Load k-point 0
        const auto& kv0 = gkk.loadKPoint(0);
        if (gkk.current_ikpt() != 0) {
            throw std::runtime_error("current_ikpt should be 0 after loadKPoint(0)");
        }
        if (kv0.kinetic.size() != 1737) {
            throw std::runtime_error("kpt 0 kinetic size mismatch: expected 1737");
        }
        if (kv0.Kx.size() != 1737 || kv0.Ky.size() != 1737 || kv0.Kz.size() != 1737) {
            throw std::runtime_error("kpt 0 Kx/Ky/Kz size mismatch");
        }

        // Spot-check some values
        if (std::abs(kv0.kinetic[0] - 0.7529635525206351) > 1e-12) {
            throw std::runtime_error("kinetic[0] mismatch: expected ~0.7529635525206351");
        }
        if (std::abs(kv0.kinetic[1] - 2.258890654552743) > 1e-12) {
            throw std::runtime_error("kinetic[1] mismatch: expected ~2.258890654552743");
        }
        if (std::abs(kv0.kinetic[2] - 5.270744861626122) > 1e-12) {
            throw std::runtime_error("kinetic[2] mismatch: expected ~5.270744861626122");
        }
        if (std::abs(kv0.Kx[0] - 1.062753653146566) > 1e-12) {
            throw std::runtime_error("Kx[0] mismatch: expected ~1.062753653146566");
        }
        if (std::abs(kv0.Kx[1] - 2.125507306293132) > 1e-12) {
            throw std::runtime_error("Kx[1] mismatch: expected ~2.125507306293132");
        }
        if (std::abs(kv0.Ky[0] - 0.6135811093611817) > 1e-12) {
            throw std::runtime_error("Ky[0] mismatch: expected ~0.6135811093611817");
        }
        if (std::abs(kv0.Kz[0] - 0.0) > 1e-15) {
            throw std::runtime_error("Kz[0] mismatch: expected 0.0");
        }

        // Verify merge across nodes: last of node 0, first of node 1
        if (std::abs(kv0.kinetic[873] - 0.7529635525206351) > 1e-12) {
            throw std::runtime_error("kinetic[873] (last of node 0) mismatch");
        }
        if (std::abs(kv0.kinetic[874] - 1.631402507680678) > 1e-12) {
            throw std::runtime_error("kinetic[874] (first of node 1) mismatch");
        }

        // Spherical and Integer views should not be available yet
        if (!kv0.q.empty() || !kv0.theta.empty() || !kv0.phi.empty()) {
            throw std::runtime_error("Spherical data should not be available without setDataView");
        }
        if (!kv0.g_idx.empty()) {
            throw std::runtime_error("Integer data should not be available without setDataView");
        }

        // Cache test: loading same k-point again should work
        const auto& kv0_cached = gkk.loadKPoint(0);
        if (kv0_cached.kinetic.size() != 1737) {
            throw std::runtime_error("cached kpt 0 kinetic size mismatch");
        }
        if (std::abs(kv0_cached.kinetic[0] - kv0.kinetic[0]) > 1e-15) {
            throw std::runtime_error("cached data mismatch");
        }

        // Test Spherical view
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Spherical);
        const auto& kv_sph = gkk.loadKPoint(0);
        if (kv_sph.q.empty() || kv_sph.theta.empty() || kv_sph.phi.empty()) {
            throw std::runtime_error("Spherical data should be available after setDataView");
        }
        if (kv_sph.q.size() != 1737) {
            throw std::runtime_error("spherical r size mismatch");
        }
        // Verify spherical relationship: r = sqrt(Kx^2 + Ky^2 + Kz^2)
        for (int i = 0; i < 10; ++i) {
            double r_expected = std::sqrt(
                kv_sph.Kx[i] * kv_sph.Kx[i] +
                kv_sph.Ky[i] * kv_sph.Ky[i] +
                kv_sph.Kz[i] * kv_sph.Kz[i]
            );
            if (std::abs(kv_sph.q[i] - r_expected) > 1e-12) {
                throw std::runtime_error("spherical r[" + std::to_string(i) + "] mismatch");
            }
        }

        // Test Integer view
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Spherical | KVecsView::Integer);
        const auto& kv_int = gkk.loadKPoint(0);
        if (kv_int.g_idx.empty()) {
            throw std::runtime_error("Integer data should be available after setDataView");
        }
        if (kv_int.g_idx.size() != 1737) {
            throw std::runtime_error("integer g_idx size mismatch");
        }
        // Verify G-vectors are centered around origin: both signs present,
        // and each component's sum is near zero (indicating no systematic offset).
        {
            const auto& gs = kv_int.g_idx;
            int cnt_pos_x = 0, cnt_neg_x = 0;
            int cnt_pos_y = 0, cnt_neg_y = 0;
            int cnt_pos_z = 0, cnt_neg_z = 0;
            long long sum_x = 0, sum_y = 0, sum_z = 0;

            for (const auto& g : gs) {
                     if (g.x > 0) ++cnt_pos_x; else if (g.x < 0) ++cnt_neg_x;
                     if (g.y > 0) ++cnt_pos_y; else if (g.y < 0) ++cnt_neg_y;
                     if (g.z > 0) ++cnt_pos_z; else if (g.z < 0) ++cnt_neg_z;
                sum_x += g.x; sum_y += g.y; sum_z += g.z;
            }
            std::println("sum_x={}, sum_y={}, sum_z={} (ng={})", sum_x, sum_y, sum_z, gs.size());

            if (cnt_pos_x == 0 || cnt_neg_x == 0)
                throw std::runtime_error("g_idx.x must have both positive and negative values");
            if (cnt_pos_y == 0 || cnt_neg_y == 0)
                throw std::runtime_error("g_idx.y must have both positive and negative values");
            if (cnt_pos_z == 0 || cnt_neg_z == 0)
                throw std::runtime_error("g_idx.z must have both positive and negative values");

            int ng = static_cast<int>(gs.size());
            if (std::abs(sum_x) > ng)
                throw std::runtime_error("g_idx.x sum too large: "
                    + std::to_string(sum_x) + " (|sum| > ng=" + std::to_string(ng) + ")");
            if (std::abs(sum_y) > ng)
                throw std::runtime_error("g_idx.y sum too large: "
                    + std::to_string(sum_y) + " (|sum| > ng=" + std::to_string(ng) + ")");
            if (std::abs(sum_z) > ng)
                throw std::runtime_error("g_idx.z sum too large: "
                    + std::to_string(sum_z) + " (|sum| > ng=" + std::to_string(ng) + ")");
        }
        // Per-k-point metadata available in Integer view
        // kpt=0 is Gamma (0,0,0)

        // Load different k-point
        const auto& kv5 = gkk.loadKPoint(5);
        if (gkk.current_ikpt() != 5) {
            throw std::runtime_error("current_ikpt should be 5 after loadKPoint(5)");
        }
        if (kv5.kinetic.size() != 1718) {
            throw std::runtime_error("kpt 5 kinetic size mismatch: expected 1718");
        }

        std::println("All GKK tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
