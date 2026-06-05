// Test: Hamiltonian step-by-step loading + loadFromDirectory + error handling.
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)
// Uses test data under test/data_io-nonlocal/

import std;
import io;
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

    std::println("\nAll tests passed.");
    return 0;
}
