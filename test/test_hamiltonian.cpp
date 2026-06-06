// Test: Hamiltonian step-by-step loading + loadFromDirectory + error handling.
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)
// Uses test data under test/data_io-nonlocal/

import std;
import io;
import utils.vector3d;
import H_psi;

auto main() -> int {
    // =====================================================================
    //  1. Default construction (base_dir = cwd)
    // =====================================================================
    {
        Hamiltonian h;
        if (h.hasGKK() || h.hasWG() || h.hasVR() ||
            h.hasATOM() || h.hasEIGEN())
        {
            std::println("FAIL: default-constructed Hamiltonian should have no data");
            return 1;
        }
        std::println("PASS: default construction, all data absent");
    }

    // =====================================================================
    //  2. Step-by-step load — nonlocal AlN data
    // =====================================================================
    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadATOM("atom.config");
        if (!h.hasATOM()) { std::println("FAIL: ATOM not loaded"); return 1; }
        if (h.hasGKK())  { std::println("FAIL: GKK should not be loaded yet"); return 1; }
        std::println("PASS: ATOM loaded, GKK absent");
    }

    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadATOM("atom.config");
        h.loadNCPPs(".");
        if (!h.hasATOM()) { std::println("FAIL: ATOM not loaded"); return 1; }
        std::println("PASS: ATOM + NCPPs loaded");

        // ATOM species: Al (Z=13), N (Z=7)
        auto ncpp_al = h.ncpp(13);
        auto ncpp_n  = h.ncpp(7);
        if (ncpp_al.meta.element != "Al") { std::println("FAIL: NCPP Al not found"); return 1; }
        if (ncpp_n.meta.element  != "N")  { std::println("FAIL: NCPP N not found"); return 1; }
        std::println("PASS: NCPP lookup by atomic number");
    }

    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadATOM ("atom.config");
        h.loadNCPPs(".");
        h.loadGKK  ("OUT.GKK");
        h.loadWG   ("OUT.WG");
        h.loadVR   ("OUT.VR");

        if (!h.hasGKK() || !h.hasWG() || !h.hasVR()) {
            std::println("FAIL: data missing after loading all required files");
            return 1;
        }

        const auto& gkk_meta = h.gkk().meta;
        const auto& wg_meta  = h.wg().meta;
        const auto& vr_meta  = h.vr().meta;

        if (gkk_meta.nkpt < 1) { std::println("FAIL: no k-points in GKK"); return 1; }
        if (wg_meta.nband < 1)    { std::println("FAIL: no bands in WG"); return 1; }
        if (vr_meta.n1 < 1)   { std::println("FAIL: VR grid empty"); return 1; }

        std::println("PASS: step-by-step load (ATOM + NCPPs + GKK + WG + VR)");
    }

    // =====================================================================
    //  3. loadFromDirectory convenience
    // =====================================================================
    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadFromDirectory();

        if (!h.hasGKK() || !h.hasWG() || !h.hasVR() ||
            !h.hasATOM() || !h.hasEIGEN())
        {
            std::println("FAIL: loadFromDirectory missing some data");
            return 1;
        }

        // Verify some cross-file consistency by accessing loaded data
        if (h.gkk().meta.nkpt != h.wg().meta.nkpt) {
            std::println("FAIL: nkpt mismatch after loadFromDirectory");
            return 1;
        }
        std::println("PASS: loadFromDirectory (GKK + WG + VR + ATOM + EIGEN)");
    }

    // =====================================================================
    //  4. Accessor throws when data not loaded
    // =====================================================================
    {
        Hamiltonian h;
        bool threw = false;
        try { (void)h.gkk(); }
        catch (const std::runtime_error&) { threw = true; }
        if (!threw) { std::println("FAIL: gkk() should throw when not loaded"); return 1; }

        threw = false;
        try { (void)h.wg(); }
        catch (const std::runtime_error&) { threw = true; }
        if (!threw) { std::println("FAIL: wg() should throw when not loaded"); return 1; }

        // ncpp with missing element
        Hamiltonian h2("test/data_io-nonlocal");
        h2.loadATOM("atom.config");
        h2.loadNCPPs(".");
        threw = false;
        try { (void)h2.ncpp(999); }
        catch (const std::out_of_range&) { threw = true; }
        if (!threw) { std::println("FAIL: ncpp(999) should throw"); return 1; }

        std::println("PASS: missing-data accessors throw correctly");
    }

    // =====================================================================
    //  5. EIGEN is truly optional — H|ψ⟩ does not require it
    // =====================================================================
    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadFromDirectory();

        // at_k should succeed even though eigen is loaded — it's optional (extra)
        // But the important test: at_k doesn't check eigen
        auto callable = h.at_k(0);
        (void)callable;
        std::println("PASS: at_k() works with EIGEN present (optional)");

        // H without EIGEN
        Hamiltonian h2("test/data_io-nonlocal");
        h2.loadATOM ("atom.config");
        h2.loadNCPPs(".");
        h2.loadGKK  ("OUT.GKK");
        h2.loadWG   ("OUT.WG");
        h2.loadVR   ("OUT.VR");
        // Note: no loadEIGEN
        auto callable2 = h2.at_k(0);
        (void)callable2;
        std::println("PASS: at_k() works without EIGEN");
    }

    // =====================================================================
    //  6. Constructor throws when base_dir does not exist
    // =====================================================================
    {
        bool threw = false;
        try { Hamiltonian h("nonexistent-directory"); }
        catch (const std::runtime_error&) { threw = true; }
        if (!threw) {
            std::println("FAIL: constructor should throw for nonexistent directory");
            return 1;
        }
        std::println("PASS: constructor rejects nonexistent directory");
    }

    // =====================================================================
    //  7. checkConsistencyExtended  —  Part 3 heavy checks
    // =====================================================================
    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadFromDirectory();

        // All checks
        h.checkConsistencyExtended();

        // Individual checks
        h.checkConsistencyExtended({ExtendedCheck::RHOReconstruct});
        h.checkConsistencyExtended({ExtendedCheck::ValenceCount});

        // Duplicate entries should be deduped
        h.checkConsistencyExtended({
            ExtendedCheck::RHOReconstruct,
            ExtendedCheck::RHOReconstruct,
        });

        std::println("PASS: checkConsistencyExtended with all data");
    }

    {
        // Partial data — checks that can't run should skip gracefully
        Hamiltonian h("test/data_io-nonlocal");
        h.loadATOM("atom.config");
        h.loadGKK("OUT.GKK");
        h.checkConsistencyExtended();   // all skip (no WG/OCC/RHO/NCPPs)
        std::println("PASS: checkConsistencyExtended skips when data missing");
    }

    // =====================================================================
    //  8. StructureFactor — cached and non-cached
    // =====================================================================
    {
        Hamiltonian h("test/data_io-nonlocal");
        h.loadFromDirectory();
        auto& gkk = h.gkk();
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        const auto& kv = gkk.loadKPoint(0);
        const auto& atom = h.atom();

        StructureFactor sf;
        double eps = 1e-12;

        // 8a. default τ = (0,0,0) → S = 1.0 for any (g,k)
        {
            vector3d<int>   g{5, -3, 2};
            vector3d<double> k{0.1, 0.2, 0.3};
            auto val = sf(g, k);
            if (std::abs(val - 1.0) > eps) {
                std::println("FAIL: default τ should give S=1.0, got ({}, {})", val.real(), val.imag());
                return 1;
            }
            std::println("PASS: default τ gives S=1.0");
        }

        // 8b. set_tau, verify with manual reference
        {
            double tau_arr[] = {0.125, 0.25, 0.375};
            sf.set_tau(tau_arr);

            vector3d<int>   g{2, -1, 3};
            vector3d<double> k{0.1, 0.2, 0.3};
            auto val = sf(g, k);

            double arg = -2.0 * std::numbers::pi
                       * ((k.x + g.x) * tau_arr[0]
                        + (k.y + g.y) * tau_arr[1]
                        + (k.z + g.z) * tau_arr[2]);
            std::complex<double> expected(std::cos(arg), std::sin(arg));
            if (std::abs(val - expected) > eps) {
                std::println("FAIL: single-G S mismatch, got ({}, {}) expected ({}, {})",
                             val.real(), val.imag(), expected.real(), expected.imag());
                return 1;
            }
            std::println("PASS: single-G operator() matches manual");
        }

        // 8c. cached: prepare + operator()(ig, k) matches per-G
        {
            const auto& g_idx = kv.g_idx;
            int ng = static_cast<int>(g_idx.size());
            sf.prepare(g_idx);

            vector3d<double> k{0.1, 0.2, 0.3};
            for (int ig = 0; ig < std::min(ng, 20); ++ig) {
                auto cached = sf(ig, k);
                auto direct = sf(g_idx[ig], k);
                if (std::abs(cached - direct) > eps) {
                    std::println("FAIL: cached ig={} mismatch", ig);
                    return 1;
                }
            }
            std::println("PASS: cached operator()(ig,k) matches direct");
        }

        // 8d. batch version matches per-G
        {
            const auto& g_idx = kv.g_idx;
            int ng = static_cast<int>(g_idx.size());
            vector3d<double> k{0.05, 0.15, 0.25};

            std::vector<std::complex<double>> batch(static_cast<std::size_t>(ng));
            sf(g_idx, k, batch);

            for (int ig = 0; ig < std::min(ng, 20); ++ig) {
                auto direct = sf(g_idx[ig], k);
                if (std::abs(batch[static_cast<std::size_t>(ig)] - direct) > eps) {
                    std::println("FAIL: batch ig={} mismatch", ig);
                    return 1;
                }
            }
            std::println("PASS: batch operator() matches per-G");
        }

        // 8e. |S| = 1.0 for all
        {
            const auto& g_idx = kv.g_idx;
            int ng = static_cast<int>(g_idx.size());
            vector3d<double> k{0.125, 0.25, 0.375};
            for (int ig = 0; ig < std::min(ng, 50); ++ig) {
                auto val = sf(g_idx[ig], k);
                double mag = std::abs(val);
                if (std::abs(mag - 1.0) > eps) {
                    std::println("FAIL: |S| != 1.0 at ig={}, got {}", ig, mag);
                    return 1;
                }
            }
            std::println("PASS: |S| = 1.0 for all checked G-vectors");
        }

        // 8f. set_frac_atomic_position alias matches set_tau
        {
            double tau_arr[] = {0.2, 0.3, 0.4};
            StructureFactor sf_a, sf_b;
            sf_a.set_tau(tau_arr);
            sf_b.set_frac_atomic_position(tau_arr);
            vector3d<int> g{1, 2, 3};
            vector3d<double> k{0.1, 0.1, 0.1};
            auto va = sf_a(g, k);
            auto vb = sf_b(g, k);
            if (std::abs(va - vb) > eps) {
                std::println("FAIL: set_frac_atomic_position != set_tau");
                return 1;
            }
            std::println("PASS: set_frac_atomic_position matches set_tau");
        }

        std::println("PASS: StructureFactor section complete");
    }

    std::println("\nAll tests passed.");
    return 0;
}
