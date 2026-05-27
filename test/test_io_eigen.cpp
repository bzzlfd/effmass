import io;
import std;


auto main() -> int {
    try {
        EIGEN eig("test/test_io-local/OUT.EIGEN");

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
        if (std::abs(eig.kpt_vec[9].x - 0.6376521918879398) > 1e-12) {
            throw std::runtime_error("kpt_vec[9].x mismatch");
        }
        if (std::abs(eig.kpt_vec[9].y - 0.12271622187223635) > 1e-12) {
            throw std::runtime_error("kpt_vec[9].y mismatch");
        }
        if (std::abs(eig.kpt_vec[9].z - 0.3313681780059838) > 1e-12) {
            throw std::runtime_error("kpt_vec[9].z mismatch");
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

        // Eigenvalue shape via metadata
        if (m.nband != 26 || m.nkpt != 10 || m.islda != 1) {
            throw std::runtime_error(
                "eigenvalue dims mismatch: got (" +
                std::to_string(m.nband) + ", " + std::to_string(m.nkpt) + ", " + std::to_string(m.islda) + ")"
            );
        }

        // Eigenvalue values: [iband, ikpt] 2-arg access (islda == 1)
        if (std::abs(eig[0, 0] - (-103.8315071578049)) > 1e-10) {
            throw std::runtime_error("eig[0, 0] mismatch");
        }
        if (std::abs(eig[1, 0] - (-103.8242383349260)) > 1e-10) {
            throw std::runtime_error("eig[1, 0] mismatch");
        }
        if (std::abs(eig[2, 0] - (-30.16721646660804)) > 1e-10) {
            throw std::runtime_error("eig[2, 0] mismatch");
        }
        if (std::abs(eig[23, 9] - 22.83226372800616) > 1e-10) {
            throw std::runtime_error("eig[23, 9] mismatch");
        }
        if (std::abs(eig[24, 9] - 23.41453607942439) > 1e-10) {
            throw std::runtime_error("eig[24, 9] mismatch");
        }
        if (std::abs(eig[25, 9] - 23.41453901951371) > 1e-10) {
            throw std::runtime_error("eig[25, 9] mismatch");
        }

        // 3-arg access [iband, ikpt, ispin]
        if (std::abs(eig[0, 0, 0] - (-103.8315071578049)) > 1e-10) {
            throw std::runtime_error("eig[0, 0, 0] mismatch");
        }
        if (std::abs(eig[25, 9, 0] - 23.41453901951371) > 1e-10) {
            throw std::runtime_error("eig[25, 9, 0] mismatch");
        }

        // print_info should not throw
        eig.print_info();

        std::println("All EIGEN tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
