import std;
import io;

constexpr double HARTREE_TO_EV = 27.211386245988;


auto main() -> int {
    try {
        EIGEN eig("test/data_io-local/OUT.EIGEN");

        const auto& m = eig.meta;

        // Metadata checks
        if (m.islda != 1) {
            throw std::runtime_error("islda mismatch: expected 1, got " + std::to_string(m.islda));
        }
        if (m.nkpt != 10) {
            throw std::runtime_error("nkpt mismatch: expected 10, got " + std::to_string(m.nkpt));
        }
        if (m.nband != 26) {
            throw std::runtime_error("nband mismatch: expected 26, got " + std::to_string(m.nband));
        }
        if (m.nref_tot_8 != 32) {
            throw std::runtime_error("nref_tot_8 mismatch: expected 32, got " + std::to_string(m.nref_tot_8));
        }
        if (m.natom != 4) {
            throw std::runtime_error("natom mismatch: expected 4, got " + std::to_string(m.natom));
        }
        if (m.nnode != 2) {
            throw std::runtime_error("nnode mismatch: expected 2, got " + std::to_string(m.nnode));
        }
        if (m.is_SO != 0) {
            throw std::runtime_error("is_SO mismatch: expected 0, got " + std::to_string(m.is_SO));
        }

        // kpt vector checks
        if (eig.kpt_vec.size() != 10) {
            throw std::runtime_error("kpt_vec size mismatch");
        }
        if (std::abs(eig.kpt_vec[0].x - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0].x mismatch");
        }
        if (std::abs(eig.kpt_vec[0].y - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0].y mismatch");
        }
        if (std::abs(eig.kpt_vec[0].z - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0].z mismatch");
        }
        // kVec operator[] access
        if (std::abs(eig.kpt_vec[0][0] - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0][0] mismatch");
        }
        if (std::abs(eig.kpt_vec[0][1] - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0][1] mismatch");
        }
        if (std::abs(eig.kpt_vec[0][2] - 0.0) > 1e-15) {
            throw std::runtime_error("kpt_vec[0][2] mismatch");
        }

        // Weight checks
        if (eig.kpt_weight.size() != 10) {
            throw std::runtime_error("kpt_weight size mismatch");
        }
        const std::vector<double> expected_weights = {
            0.02, 0.02, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12
        };
        double sum_weight = 0.0;
        for (int ikpt = 0; ikpt < 10; ++ikpt) {
            sum_weight += eig.kpt_weight[ikpt];
            if (std::abs(eig.kpt_weight[ikpt] - expected_weights[ikpt]) > 1e-15) {
                throw std::runtime_error(
                    "kpt_weight[" + std::to_string(ikpt) + "] mismatch: expected " +
                    std::to_string(expected_weights[ikpt]) + ", got " +
                    std::to_string(eig.kpt_weight[ikpt])
                );
            }
        }
        if (std::abs(sum_weight - 1.0) > 1e-12) {
            throw std::runtime_error("sum(kpt_weight) mismatch: expected 1.0, got " + std::to_string(sum_weight));
        }

        // print_info should not throw
        eig.print_info();

        // ================================================================
        // Part 2: Cross-validate EIGEN against OUT.OCC (same k-mesh)
        // ================================================================

        EIGEN eig_nonlocal("test/data_io-nonlocal/OUT.EIGEN");
        OCC occ("test/data_io-nonlocal/OUT.OCC");

        if (occ.meta.nkpt != eig_nonlocal.meta.nkpt) {
            throw std::runtime_error(
                "nkpt mismatch: OCC=" + std::to_string(occ.meta.nkpt)
                + " EIGEN=" + std::to_string(eig_nonlocal.meta.nkpt)
            );
        }
        if (occ.meta.nband != eig_nonlocal.meta.nband) {
            throw std::runtime_error(
                "nband mismatch: OCC=" + std::to_string(occ.meta.nband)
                + " EIGEN=" + std::to_string(eig_nonlocal.meta.nband)
            );
        }

        // k-point vectors: EIGEN full-precision, OCC rounded to 4 dp
        for (int ik = 0; ik < occ.meta.nkpt; ++ik) {
            if (std::abs(eig_nonlocal.kpt_vec[ik].x - occ.kpt_vec[ik].x) > 1e-4) {
                throw std::runtime_error(
                    "k-point " + std::to_string(ik + 1) + " x mismatch"
                );
            }
            if (std::abs(eig_nonlocal.kpt_vec[ik].y - occ.kpt_vec[ik].y) > 1e-4) {
                throw std::runtime_error(
                    "k-point " + std::to_string(ik + 1) + " y mismatch"
                );
            }
            if (std::abs(eig_nonlocal.kpt_vec[ik].z - occ.kpt_vec[ik].z) > 1e-4) {
                throw std::runtime_error(
                    "k-point " + std::to_string(ik + 1) + " z mismatch"
                );
            }
        }

        // Energies: EIGEN (Hartree) → eV vs OCC (eV, rounded to 4 dp)
        for (int ik = 0; ik < occ.meta.nkpt; ++ik) {
            for (int ib = 0; ib < occ.meta.nband; ++ib) {
                double eig_ev = eig_nonlocal[ib, ik] * HARTREE_TO_EV;
                double occ_ev = occ.energy(ib, ik);
                if (std::abs(eig_ev - occ_ev) > 0.001) {
                    throw std::runtime_error(
                        "Energy mismatch at k-point " + std::to_string(ik + 1)
                        + ", band " + std::to_string(ib + 1)
                        + ": EIGEN=" + std::to_string(eig_ev)
                        + " eV, OCC=" + std::to_string(occ_ev) + " eV"
                    );
                }
            }
        }

        // ================================================================
        // Part 3: Direct OCC-specific tests
        // ================================================================

        {
            // Metadata checks via OCC
            if (occ.meta.islda != 1) {
                throw std::runtime_error("OCC islda mismatch: expected 1, got " + std::to_string(occ.meta.islda));
            }
            if (occ.meta.nkpt != 10) {
                throw std::runtime_error("OCC nkpt mismatch: expected 10, got " + std::to_string(occ.meta.nkpt));
            }
            if (occ.meta.nband != 26) {
                throw std::runtime_error("OCC nband mismatch: expected 26, got " + std::to_string(occ.meta.nband));
            }

            // kpt_vec consistency: both EIGEN and OCC use the same k-mesh
            if (occ.kpt_vec.size() != 10) {
                throw std::runtime_error("OCC kpt_vec size mismatch");
            }
            // OCC k-point coordinates are rounded to 4 dp, so tolerance is 1e-4
            for (int ik = 0; ik < occ.meta.nkpt; ++ik) {
                if (std::abs(eig_nonlocal.kpt_vec[ik].x - occ.kpt_vec[ik].x) > 1e-4) {
                    throw std::runtime_error("OCC kpt_vec[" + std::to_string(ik) + "].x mismatch");
                }
                if (std::abs(eig_nonlocal.kpt_vec[ik].y - occ.kpt_vec[ik].y) > 1e-4) {
                    throw std::runtime_error("OCC kpt_vec[" + std::to_string(ik) + "].y mismatch");
                }
                if (std::abs(eig_nonlocal.kpt_vec[ik].z - occ.kpt_vec[ik].z) > 1e-4) {
                    throw std::runtime_error("OCC kpt_vec[" + std::to_string(ik) + "].z mismatch");
                }
            }

            // Occupations should be non-negative
            for (int ik = 0; ik < occ.meta.nkpt; ++ik) {
                for (int ib = 0; ib < occ.meta.nband; ++ib) {
                    double occ_val = occ.occupation(ib, ik);
                    if (occ_val < 0.0) {
                        throw std::runtime_error(
                            "OCC negative occupation at kpt=" + std::to_string(ik + 1)
                            + " band=" + std::to_string(ib + 1)
                        );
                    }
                }
            }

            // print_info should not throw
            occ.print_info();
        }

        // ================================================================
        // Part 4: Spin=2 OCC test
        // ================================================================
        {
            OCC occ_spin2("test/data_io-nonlocal/OUT.OCC_SPIN2");

            if (occ_spin2.meta.islda != 2) {
                throw std::runtime_error(
                    "OCC spin2 islda mismatch: expected 2, got "
                    + std::to_string(occ_spin2.meta.islda)
                );
            }
            if (occ_spin2.meta.nkpt != occ.meta.nkpt) {
                throw std::runtime_error(
                    "OCC spin2 nkpt mismatch: expected " + std::to_string(occ.meta.nkpt)
                    + ", got " + std::to_string(occ_spin2.meta.nkpt)
                );
            }
            if (occ_spin2.meta.nband != 31) {
                throw std::runtime_error(
                    "OCC spin2 nband mismatch: expected 31, got "
                    + std::to_string(occ_spin2.meta.nband)
                );
            }

            // kpt_vec: same k-mesh as non-spin OCC
            if (occ_spin2.kpt_vec.size() != occ.meta.nkpt) {
                throw std::runtime_error("OCC spin2 kpt_vec size mismatch");
            }
            for (int ik = 0; ik < occ_spin2.meta.nkpt; ++ik) {
                if (std::abs(occ_spin2.kpt_vec[ik].x - occ.kpt_vec[ik].x) > 1e-4
                    || std::abs(occ_spin2.kpt_vec[ik].y - occ.kpt_vec[ik].y) > 1e-4
                    || std::abs(occ_spin2.kpt_vec[ik].z - occ.kpt_vec[ik].z) > 1e-4)
                {
                    throw std::runtime_error(
                        "OCC spin2 kpt_vec[" + std::to_string(ik) + "] mismatch"
                    );
                }
            }

            // Occupations non-negative for both spins
            for (int is = 0; is < 2; ++is) {
                for (int ik = 0; ik < occ_spin2.meta.nkpt; ++ik) {
                    for (int ib = 0; ib < occ_spin2.meta.nband; ++ib) {
                        double occ_val = occ_spin2.occupation(ib, ik, is);
                        if (occ_val < 0.0) {
                            throw std::runtime_error(
                                "OCC spin2 negative occupation at spin="
                                + std::to_string(is) + " kpt=" + std::to_string(ik + 1)
                            );
                        }
                    }
                }
            }

            // Spin 0 and spin 1 data should differ somewhere
            bool found_difference = false;
            for (int ik = 0; ik < occ_spin2.meta.nkpt && !found_difference; ++ik) {
                for (int ib = 0; ib < occ_spin2.meta.nband && !found_difference; ++ib) {
                    double e0 = occ_spin2.energy(ib, ik, 0);
                    double e1 = occ_spin2.energy(ib, ik, 1);
                    if (std::abs(e0 - e1) > 1e-6) found_difference = true;
                }
            }
            if (!found_difference) {
                throw std::runtime_error("OCC spin2: spin 0 and 1 energies are identical");
            }

            occ_spin2.print_info();
        }

        std::println("All EIGEN / OCC tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
