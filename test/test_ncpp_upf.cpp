import pseudo;
import std;

auto main() -> int {
    try {
        NCPPUPF upf("test/test_ncpp_upf/Ge-spd-high.PD04.PBE.UPF");

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
        for (int i = 0; i < 6; ++i) {
            if (nl.beta[i].size() != 1560) {
                throw std::runtime_error("beta[" + std::to_string(i) + "] size mismatch");
            }
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
        for (int i = 0; i < 5; ++i) {
            if (wfc.chi[i].size() != 1560) {
                throw std::runtime_error("chi[" + std::to_string(i) + "] size mismatch");
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

        // Rho atom
        if (upf.rhoAtom().size() != 1560) {
            throw std::runtime_error("rhoAtom size mismatch");
        }

        std::println("All NCPPUPF tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
