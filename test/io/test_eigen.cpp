import std;
import io;

#define OCC_READER_IMPORT_STD
#include "occ_reader.hpp"

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
        auto occ = parseOCC("test/data_io-nonlocal/OUT.OCC");

        if (occ.nkpt != eig_nonlocal.meta.nkpt) {
            throw std::runtime_error(
                "nkpt mismatch: OCC=" + std::to_string(occ.nkpt)
                + " EIGEN=" + std::to_string(eig_nonlocal.meta.nkpt)
            );
        }
        if (occ.nband != eig_nonlocal.meta.nband) {
            throw std::runtime_error(
                "nband mismatch: OCC=" + std::to_string(occ.nband)
                + " EIGEN=" + std::to_string(eig_nonlocal.meta.nband)
            );
        }

        // k-point vectors: EIGEN full-precision, OCC rounded to 4 dp
        for (int ik = 0; ik < occ.nkpt; ++ik) {
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
        for (int ik = 0; ik < occ.nkpt; ++ik) {
            for (int ib = 0; ib < occ.nband; ++ib) {
                double eig_ev = eig_nonlocal[ib, ik] * HARTREE_TO_EV;
                double occ_ev = occ.energy(ik, ib);
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

        std::println("All EIGEN tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
