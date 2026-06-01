import pseudo;
import std;


auto main() -> int {
    try {
        UPF upf("test/data_io_upf/Si-PD04-FR.upf");

        const auto& h = upf.header();
        if (h.element != "Si") throw std::runtime_error("element mismatch");
        if (h.relativistic != "full") throw std::runtime_error("relativistic mismatch");
        if (!h.has_so) throw std::runtime_error("has_so should be true");
        if (h.number_of_proj != 10) throw std::runtime_error("nproj should be 10 for SOC Si");
        if (h.mesh_size != 1528) throw std::runtime_error("mesh_size mismatch");

        // Check jjj values from PP_SPIN_ORB/PP_RELBETA
        const auto& nl = upf.nonlocal();
        if (nl.jjj.size() != 10) throw std::runtime_error("jjj size mismatch");

        // Expected j values from PP_RELBETA:
        // index 1-2  (l=0): j=0.5
        // index 3-4  (l=1): j=0.5, 1.5
        // index 5-6  (l=1): j=0.5, 1.5
        // index 7-8  (l=2): j=1.5, 2.5
        // index 9-10 (l=2): j=1.5, 2.5
        const std::vector<double> expected_jjj = {0.5, 0.5, 0.5, 1.5, 0.5, 1.5, 1.5, 2.5, 1.5, 2.5};
        for (int i = 0; i < 10; ++i) {
            if (std::abs(nl.jjj[i] - expected_jjj[i]) > 1e-10) {
                throw std::runtime_error("jjj[" + std::to_string(i) + "] mismatch: expected " +
                    std::to_string(expected_jjj[i]) + ", got " + std::to_string(nl.jjj[i]));
            }
        }

        // Check jchi from PP_SPIN_ORB/PP_RELWFC
        const auto& wfc = upf.wavefunctions();
        if (wfc.jchi.size() != 3) throw std::runtime_error("jchi size mismatch");
        // PP_RELWFC: 1: jchi=0.5 (3S), 2: jchi=1.5 (3P), 3: jchi=0.5 (3P)
        if (std::abs(wfc.jchi[0] - 0.5) > 1e-10) throw std::runtime_error("jchi[0] mismatch");
        if (std::abs(wfc.jchi[1] - 1.5) > 1e-10) throw std::runtime_error("jchi[1] mismatch");
        if (std::abs(wfc.jchi[2] - 0.5) > 1e-10) throw std::runtime_error("jchi[2] mismatch");

        // Check beta data integrity (same format as non-SOC)
        if (nl.beta.size() != 10) throw std::runtime_error("beta count mismatch");
        for (int i = 0; i < 10; ++i) {
            if (nl.beta[i].size() == 0 || nl.beta[i].size() > 1528) {
                throw std::runtime_error("beta[" + std::to_string(i) + "] invalid size");
            }
            if (nl.beta[i].size() != static_cast<std::size_t>(nl.kbeta[i])) {
                throw std::runtime_error("beta[" + std::to_string(i) + "] size mismatch with kbeta");
            }
        }

        if (nl.dion.rows != 10 || nl.dion.cols != 10) {
            throw std::runtime_error("dion should be 10x10 for SOC");
        }
        if (nl.dion.size() != 100) {
            throw std::runtime_error("dion size mismatch: expected 100");
        }

        // Non-SOC file should have jjj initialized to 0.0
        UPF ncsoc("test/data_io_upf/Ge-spd-high.PD04.PBE.UPF");
        const auto& nl2 = ncsoc.nonlocal();
        if (nl2.jjj.size() != 6) throw std::runtime_error("non-SOC jjj size should be 6");
        for (int i = 0; i < 6; ++i) {
            if (nl2.jjj[i] != 0.0) {
                throw std::runtime_error("non-SOC jjj[" + std::to_string(i) + "] should be 0.0");
            }
        }
        if (ncsoc.wavefunctions().jchi.size() != 5) throw std::runtime_error("non-SOC jchi size mismatch");
        for (int i = 0; i < 5; ++i) {
            if (ncsoc.wavefunctions().jchi[i] != 0.0) {
                throw std::runtime_error("non-SOC jchi[" + std::to_string(i) + "] should be 0.0");
            }
        }

        std::println("All SOC UPF tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
