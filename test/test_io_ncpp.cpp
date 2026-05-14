import pseudo;
import std;


auto main() -> int {
    try {
        NCPPUPF upf("test/test_io_ncpp/Ge-spd-high.PD04.PBE.UPF");

        const auto& h = upf.header();

        // Header checks
        if (h.element != "Ge") {
            throw std::runtime_error("element mismatch: expected Ge, got " + h.element);
        }
        if (h.pseudo_type != "NC") {
            throw std::runtime_error("pseudo_type mismatch: expected NC, got " + h.pseudo_type);
        }
        if (h.z_valence != 22.0) {
            throw std::runtime_error("z_valence mismatch");
        }
        if (h.mesh_size != 1560) {
            throw std::runtime_error("mesh_size mismatch");
        }
        if (h.number_of_wfc != 5) {
            throw std::runtime_error("number_of_wfc mismatch");
        }
        if (h.number_of_proj != 6) {
            throw std::runtime_error("number_of_proj mismatch");
        }
        if (h.l_max != 2) {
            throw std::runtime_error("l_max mismatch");
        }
        if (!h.core_correction) {
            throw std::runtime_error("core_correction should be true for Ge");
        }

        // Mesh checks
        const auto& mesh = upf.mesh();
        if (mesh.r.size() != 1560) {
            throw std::runtime_error("mesh.r size mismatch");
        }
        if (mesh.rab.size() != 1560) {
            throw std::runtime_error("mesh.rab size mismatch");
        }

        // Local potential
        if (upf.localPotential().size() != 1560) {
            throw std::runtime_error("localPotential size mismatch");
        }

        // Nonlocal checks
        const auto& nl = upf.nonlocal();
        if (nl.beta.size() != 6) {
            throw std::runtime_error("beta count mismatch");
        }

        // Expected beta attributes for Ge-spd-high.PD04.PBE.UPF
        const std::vector<int> expected_lll = {0, 0, 1, 1, 2, 2};
        const std::vector<int> expected_kbeta = {156, 156, 156, 156, 156, 156};
        const std::vector<double> expected_rcut = {
            1.55, 1.55, 1.55, 1.55, 1.55, 1.55
        };
        // Spot-check beta values from raw UPF data (before truncation)
        const std::vector<double> expected_beta_0_0 = {1.6747178089e-09};
        const std::vector<double> expected_beta_0_127 = {-3.7097811099e-07};

        for (int i = 0; i < 6; ++i) {
            if (nl.lll[i] != expected_lll[i]) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] lll mismatch: expected " + std::to_string(expected_lll[i]) +
                    ", got " + std::to_string(nl.lll[i]));
            }
            if (nl.kbeta[i] != expected_kbeta[i]) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] kbeta mismatch: expected " + std::to_string(expected_kbeta[i]) +
                    ", got " + std::to_string(nl.kbeta[i]));
            }
            if (std::abs(nl.rcut[i] - expected_rcut[i]) > 1e-10) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] rcut mismatch");
            }
            if (nl.beta[i].size() != static_cast<std::size_t>(nl.kbeta[i])) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] size mismatch with kbeta");
            }
            if (nl.beta[i].size() == 0 || nl.beta[i].size() > 1560) {
                throw std::runtime_error("beta[" + std::to_string(i) + "] invalid size");
            }
        }

        // Spot-check: first and last non-zero elements of beta[0] must match raw data
        if (std::abs(nl.beta[0][0] - expected_beta_0_0[0]) > 1e-20) {
            throw std::runtime_error("beta[0][0] value mismatch");
        }
        if (std::abs(nl.beta[0][127] - expected_beta_0_127[0]) > 1e-20) {
            throw std::runtime_error("beta[0][127] value mismatch");
        }

        if (nl.dion.size() != 36) {
            throw std::runtime_error("dion size mismatch");
        }

        // Wavefunction checks
        const auto& wfc = upf.wavefunctions();
        const auto& rab = mesh.rab;
        if (wfc.chi.size() != 5) {
            throw std::runtime_error("chi count mismatch");
        }

        // Expected chi attributes for Ge-spd-high.PD04.PBE.UPF
        const std::vector<int> expected_lchi = {0, 1, 2, 0, 1};
        const std::vector<std::string> expected_labels = {"3S", "3P", "3D", "4S", "4P"};
        const std::vector<double> expected_oc = {2.0, 6.0, 10.0, 2.0, 2.0};
        const std::vector<int> expected_kchi = {1273, 1368, 1560, 1560, 1560};
        // Spot-check chi values from raw UPF data
        const double expected_chi_0_0 = -9.894862707e-12;
        const double expected_chi_0_1272 = -5.7978579654e-21;

        for (int i = 0; i < 5; ++i) {
            if (wfc.lchi[i] != expected_lchi[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] lchi mismatch: expected " + std::to_string(expected_lchi[i]) +
                    ", got " + std::to_string(wfc.lchi[i]));
            }
            if (wfc.labels[i] != expected_labels[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] label mismatch: expected " + expected_labels[i] +
                    ", got " + wfc.labels[i]);
            }
            if (std::abs(wfc.oc[i] - expected_oc[i]) > 1e-10) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] occupation mismatch");
            }
            if (wfc.kchi[i] != expected_kchi[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] kchi mismatch: expected " + std::to_string(expected_kchi[i]) +
                    ", got " + std::to_string(wfc.kchi[i]));
            }
            if (wfc.chi[i].size() != static_cast<std::size_t>(wfc.kchi[i])) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] size mismatch with kchi");
            }
            if (wfc.chi[i].size() == 0 || wfc.chi[i].size() > 1560) {
                throw std::runtime_error("chi[" + std::to_string(i) + "] invalid size");
            }

            // Truncation boundary: if size < mesh_size, last element must be non-zero
            if (wfc.chi[i].size() < static_cast<std::size_t>(h.mesh_size)) {
                if (wfc.chi[i].back() == 0.0) {
                    throw std::runtime_error("chi[" + std::to_string(i) +
                        "] last element is zero after truncation");
                }
            }

            // Norm check: ∫ |chi|² dr should be 1 for bound states
            double norm = 0.0;
            for (std::size_t ir = 0; ir < wfc.chi[i].size(); ++ir) {
                norm += wfc.chi[i][ir] * wfc.chi[i][ir] * rab[ir];
            }
            if (std::abs(norm - 1.0) > 1e-5) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] norm mismatch: expected 1.0, got " + std::to_string(norm));
            }
        }

        // Spot-check: first and last non-zero elements of chi[0] must match raw data
        if (std::abs(wfc.chi[0][0] - expected_chi_0_0) > 1e-20) {
            throw std::runtime_error("chi[0][0] value mismatch");
        }
        if (std::abs(wfc.chi[0][1272] - expected_chi_0_1272) > 1e-20) {
            throw std::runtime_error("chi[0][1272] value mismatch");
        }

        // Rho atom
        if (upf.rhoAtom().size() != 1560) {
            throw std::runtime_error("rhoAtom size mismatch");
        }

        // Mesh type inference
        auto mt = upf.meshType();
        if (mt == MeshType::Unknown) {
            throw std::runtime_error("meshType should not be Unknown for Ge-spd-high");
        }

        // nonlocalByL tests
        // Ge-spd-high has lll = {0, 0, 1, 1, 2, 2}
        auto nl0 = upf.nonlocalByL(0);
        if (nl0.beta.size() != 2) {
            throw std::runtime_error("nonlocalByL(0) beta count mismatch: expected 2, got " +
                std::to_string(nl0.beta.size()));
        }
        if (nl0.dion.rows != 2 || nl0.dion.cols != 2) {
            throw std::runtime_error("nonlocalByL(0) dion size mismatch: expected 2x2");
        }
        if (nl0.dion.size() != 4) {
            throw std::runtime_error("nonlocalByL(0) dion data size mismatch");
        }

        auto nl1 = upf.nonlocalByL(1);
        if (nl1.beta.size() != 2) {
            throw std::runtime_error("nonlocalByL(1) beta count mismatch: expected 2, got " +
                std::to_string(nl1.beta.size()));
        }
        if (nl1.dion.rows != 2 || nl1.dion.cols != 2) {
            throw std::runtime_error("nonlocalByL(1) dion size mismatch: expected 2x2");
        }

        auto nl2 = upf.nonlocalByL(2);
        if (nl2.beta.size() != 2) {
            throw std::runtime_error("nonlocalByL(2) beta count mismatch: expected 2, got " +
                std::to_string(nl2.beta.size()));
        }
        if (nl2.dion.rows != 2 || nl2.dion.cols != 2) {
            throw std::runtime_error("nonlocalByL(2) dion size mismatch: expected 2x2");
        }

        auto nl3 = upf.nonlocalByL(3);
        if (nl3.beta.size() != 0) {
            throw std::runtime_error("nonlocalByL(3) beta count mismatch: expected 0, got " +
                std::to_string(nl3.beta.size()));
        }
        if (nl3.dion.rows != 0 || nl3.dion.cols != 0) {
            throw std::runtime_error("nonlocalByL(3) dion size mismatch: expected 0x0");
        }

        // Consistency check: values in nonlocalByL must match original
        for (int i = 0; i < 2; ++i) {
            if (nl0.beta[i] != nl.beta[i]) {
                throw std::runtime_error("nonlocalByL(0) beta[" + std::to_string(i) +
                    "] value mismatch");
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (nl1.beta[i] != nl.beta[i + 2]) {
                throw std::runtime_error("nonlocalByL(1) beta[" + std::to_string(i) +
                    "] value mismatch");
            }
        }
        for (int i = 0; i < 2; ++i) {
            if (nl2.beta[i] != nl.beta[i + 4]) {
                throw std::runtime_error("nonlocalByL(2) beta[" + std::to_string(i) +
                    "] value mismatch");
            }
        }

        // Spot-check dion submatrix values
        if (std::abs(nl0.dion[0, 0] - nl.dion[0, 0]) > 1e-20) {
            throw std::runtime_error("nonlocalByL(0) dion[0,0] mismatch");
        }
        if (std::abs(nl0.dion[1, 1] - nl.dion[1, 1]) > 1e-20) {
            throw std::runtime_error("nonlocalByL(0) dion[1,1] mismatch");
        }
        if (std::abs(nl1.dion[0, 0] - nl.dion[2, 2]) > 1e-20) {
            throw std::runtime_error("nonlocalByL(1) dion[0,0] mismatch");
        }
        if (std::abs(nl2.dion[1, 1] - nl.dion[5, 5]) > 1e-20) {
            throw std::runtime_error("nonlocalByL(2) dion[1,1] mismatch");
        }

        std::println("All NCPPUPF tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
