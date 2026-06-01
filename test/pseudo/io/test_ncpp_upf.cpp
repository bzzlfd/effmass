import pseudo;
import std;


auto main() -> int {
    try {
        UPF upf("test/data_io_upf/Ge-spd-high.PD04.PBE.UPF");
        NCPP ncpp(upf);

        const int meshSize = ncpp.meta.mesh_size;

        // Meta checks
        if (ncpp.meta.element != "Ge") {
            throw std::runtime_error("element mismatch: expected Ge, got " + ncpp.meta.element);
        }
        if (ncpp.meta.pseudo_type != PseudoType::NC) {
            throw std::runtime_error("pseudo_type mismatch: expected NC");
        }
        if (ncpp.meta.z_valence != 22.0) {
            throw std::runtime_error("z_valence mismatch");
        }
        if (meshSize != 1560) {
            throw std::runtime_error("mesh_size mismatch");
        }
        if (ncpp.meta.number_of_wfc != 5) {
            throw std::runtime_error("number_of_wfc mismatch");
        }
        if (ncpp.meta.number_of_proj != 6) {
            throw std::runtime_error("number_of_proj mismatch");
        }
        if (ncpp.meta.l_max != 2) {
            throw std::runtime_error("l_max mismatch");
        }
        if (!ncpp.meta.core_correction) {
            throw std::runtime_error("core_correction should be true for Ge");
        }

        // Mesh checks
        if (ncpp.mesh.r.size() != 1560) {
            throw std::runtime_error("mesh.r size mismatch");
        }
        if (ncpp.mesh.rab.size() != 1560) {
            throw std::runtime_error("mesh.rab size mismatch");
        }

        // Local potential
        if (ncpp.local.size() != 1560) {
            throw std::runtime_error("local size mismatch");
        }

        // Nonlocal checks
        if (ncpp.meta.number_of_proj != 6) {
            throw std::runtime_error("beta count mismatch");
        }

        const auto& nl = ncpp.nonlocal;

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
            if (nl.angular_momentum[i] != expected_lll[i]) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] angular_momentum mismatch: expected " + std::to_string(expected_lll[i]) +
                    ", got " + std::to_string(nl.angular_momentum[i]));
            }
            if (nl.cutoff_index[i] != expected_kbeta[i]) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] cutoff_index mismatch: expected " + std::to_string(expected_kbeta[i]) +
                    ", got " + std::to_string(nl.cutoff_index[i]));
            }
            if (std::abs(nl.cutoff_radius[i] - expected_rcut[i]) > 1e-10) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] cutoff_radius mismatch");
            }
            if (nl.beta[i].size() != static_cast<std::size_t>(nl.cutoff_index[i])) {
                throw std::runtime_error("beta[" + std::to_string(i) +
                    "] size mismatch with cutoff_index");
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

        // Wavefunction checks
        if (ncpp.pseudoWfc.chi.size() != 5) {
            throw std::runtime_error("chi count mismatch");
        }

        const auto& wfc = ncpp.pseudoWfc;
        const auto& rab = ncpp.mesh.rab;

        // Expected chi attributes for Ge-spd-high.PD04.PBE.UPF
        const std::vector<int> expected_lchi = {0, 1, 2, 0, 1};
        const std::vector<std::string> expected_labels = {"3S", "3P", "3D", "4S", "4P"};
        const std::vector<double> expected_oc = {2.0, 6.0, 10.0, 2.0, 2.0};
        const std::vector<int> expected_kchi = {1273, 1368, 1560, 1560, 1560};
        // Spot-check chi values from raw UPF data
        const double expected_chi_0_0 = -9.894862707e-12;
        const double expected_chi_0_1272 = -5.7978579654e-21;

        for (int i = 0; i < 5; ++i) {
            if (wfc.angular_momentum[i] != expected_lchi[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] angular_momentum mismatch: expected " + std::to_string(expected_lchi[i]) +
                    ", got " + std::to_string(wfc.angular_momentum[i]));
            }
            if (wfc.label[i] != expected_labels[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] label mismatch: expected " + expected_labels[i] +
                    ", got " + wfc.label[i]);
            }
            if (std::abs(wfc.occupation[i] - expected_oc[i]) > 1e-10) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] occupation mismatch");
            }
            if (wfc.cutoff_index[i] != expected_kchi[i]) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] cutoff_index mismatch: expected " + std::to_string(expected_kchi[i]) +
                    ", got " + std::to_string(wfc.cutoff_index[i]));
            }
            if (wfc.chi[i].size() != static_cast<std::size_t>(wfc.cutoff_index[i])) {
                throw std::runtime_error("chi[" + std::to_string(i) +
                    "] size mismatch with cutoff_index");
            }
            if (wfc.chi[i].size() == 0 || wfc.chi[i].size() > 1560) {
                throw std::runtime_error("chi[" + std::to_string(i) + "] invalid size");
            }

            // Truncation boundary: if size < mesh_size, last element must be non-zero
            if (wfc.chi[i].size() < static_cast<std::size_t>(meshSize)) {
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
        if (ncpp.rho_atom.size() != 1560) {
            throw std::runtime_error("rho_atom size mismatch");
        }

        // Mesh type inference
        if (ncpp.inferMeshType() != MeshType::Uniform) {
            throw std::runtime_error("inferMeshType should be Uniform for Ge-spd-high");
        }

        // projectorBlock tests
        // Ge-spd-high has angular_momentum = {0, 0, 1, 1, 2, 2}
        auto nl0 = ncpp.projectorBlock(0);
        if (nl0.beta.size() != 2) {
            throw std::runtime_error("projectorBlock(0) beta count mismatch: expected 2, got " +
                std::to_string(nl0.beta.size()));
        }
        if (nl0.B.rows != 2 || nl0.B.cols != 2) {
            throw std::runtime_error("projectorBlock(0) B size mismatch: expected 2x2");
        }
        if (nl0.B.size() != 4) {
            throw std::runtime_error("projectorBlock(0) B data size mismatch");
        }

        auto nl1 = ncpp.projectorBlock(1);
        if (nl1.beta.size() != 2) {
            throw std::runtime_error("projectorBlock(1) beta count mismatch: expected 2, got " +
                std::to_string(nl1.beta.size()));
        }
        if (nl1.B.rows != 2 || nl1.B.cols != 2) {
            throw std::runtime_error("projectorBlock(1) B size mismatch: expected 2x2");
        }

        auto nl2 = ncpp.projectorBlock(2);
        if (nl2.beta.size() != 2) {
            throw std::runtime_error("projectorBlock(2) beta count mismatch: expected 2, got " +
                std::to_string(nl2.beta.size()));
        }
        if (nl2.B.rows != 2 || nl2.B.cols != 2) {
            throw std::runtime_error("projectorBlock(2) B size mismatch: expected 2x2");
        }

        auto nl3 = ncpp.projectorBlock(3);
        if (nl3.beta.size() != 0) {
            throw std::runtime_error("projectorBlock(3) beta count mismatch: expected 0, got " +
                std::to_string(nl3.beta.size()));
        }
        if (nl3.B.rows != 0 || nl3.B.cols != 0) {
            throw std::runtime_error("projectorBlock(3) B size mismatch: expected 0x0");
        }

        // Consistency check: projectorBlock values must match full nonlocal storage
        if (nl0.beta[0] != nl.beta[0]) {
            throw std::runtime_error("projectorBlock(0) beta[0] value mismatch");
        }
        if (nl0.beta[1] != nl.beta[1]) {
            throw std::runtime_error("projectorBlock(0) beta[1] value mismatch");
        }
        if (nl1.beta[0] != nl.beta[2]) {
            throw std::runtime_error("projectorBlock(1) beta[0] value mismatch");
        }
        if (nl1.beta[1] != nl.beta[3]) {
            throw std::runtime_error("projectorBlock(1) beta[1] value mismatch");
        }
        if (nl2.beta[0] != nl.beta[4]) {
            throw std::runtime_error("projectorBlock(2) beta[0] value mismatch");
        }
        if (nl2.beta[1] != nl.beta[5]) {
            throw std::runtime_error("projectorBlock(2) beta[1] value mismatch");
        }

        // Spot-check B submatrix values
        if (std::abs(nl0.B[0, 0] - nl.B[0, 0]) > 1e-20) {
            throw std::runtime_error("projectorBlock(0) B[0,0] mismatch");
        }
        if (std::abs(nl0.B[1, 1] - nl.B[1, 1]) > 1e-20) {
            throw std::runtime_error("projectorBlock(0) B[1,1] mismatch");
        }
        if (std::abs(nl1.B[0, 0] - nl.B[2, 2]) > 1e-20) {
            throw std::runtime_error("projectorBlock(1) B[0,0] mismatch");
        }
        if (std::abs(nl2.B[1, 1] - nl.B[5, 5]) > 1e-20) {
            throw std::runtime_error("projectorBlock(2) B[1,1] mismatch");
        }

        std::println("All UPF tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
