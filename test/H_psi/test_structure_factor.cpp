// Test: StructureFactor — None and Separable modes.
//
// WORKING_DIRECTORY: ${CMAKE_SOURCE_DIR} (project root)
// Uses test data under test/data_scf/

import std;
import io;
import utils.vector3d;
import H_psi;

auto main() -> int {
    // =====================================================================
    //  StructureFactor — None and Separable modes
    // =====================================================================
    {
        Hamiltonian h("test/data_scf");
        h.loadFromDirectory();
        auto& gkk = h.gkk();
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        const auto& kv = gkk.loadKPoint(0);
        double eps = 1e-12;

        auto n1 = h.vr().meta.n1;
        auto n2 = h.vr().meta.n2;
        auto n3 = h.vr().meta.n3;

        // 1. None mode, τ=(0,0,0) → S=1.0 for any (g,k)
        {
            StructureFactor sf({0.0, 0.0, 0.0}, 0, 0, 0, StructureFactor::CacheMode::None);
            vector3d<int>   g{5, -3, 2};
            vector3d<double> k{0.1, 0.2, 0.3};
            auto val = sf(g, k);
            if (std::abs(val - 1.0) > eps) {
                std::println("FAIL: τ=0 should give S=1.0, got ({}, {})", val.real(), val.imag());
                return 1;
            }
            std::println("PASS: τ=0 gives S=1.0 (None mode)");
        }

        // 2. reset_tau, verify with manual reference
        {
            StructureFactor sf({0.0, 0.0, 0.0}, 0, 0, 0, StructureFactor::CacheMode::None);
            sf.reset_tau({0.125, 0.25, 0.375});

            vector3d<int>   g{2, -1, 3};
            vector3d<double> k{0.1, 0.2, 0.3};
            auto val = sf(g, k);

            double arg = -2.0 * std::numbers::pi
                       * ((k.x + g.x) * 0.125
                        + (k.y + g.y) * 0.25
                        + (k.z + g.z) * 0.375);
            std::complex<double> expected(std::cos(arg), std::sin(arg));
            if (std::abs(val - expected) > eps) {
                std::println("FAIL: single-G S mismatch, got ({}, {}) expected ({}, {})",
                             val.real(), val.imag(), expected.real(), expected.imag());
                return 1;
            }
            std::println("PASS: single-G operator() matches manual (None mode)");
        }

        // 3. Separable mode matches None mode for GKK g_idx
        {
            const auto& g_idx = kv.g_idx;
            int ng = static_cast<int>(g_idx.size());

            vector3d<double> tau{0.125, 0.25, 0.375};
            StructureFactor sf_sep(tau, n1, n2, n3, StructureFactor::CacheMode::Separable);
            StructureFactor sf_none(tau, 0, 0, 0, StructureFactor::CacheMode::None);

            vector3d<double> k{0.1, 0.2, 0.3};
            for (int ig = 0; ig < std::min(ng, 20); ++ig) {
                auto cached = sf_sep(g_idx[ig], k);
                auto direct = sf_none(g_idx[ig], k);
                if (std::abs(cached - direct) > eps) {
                    std::println("FAIL: Separable ig={} mismatch with None, "
                                 "g=({},{},{}), sep=({},{}), ref=({},{})",
                                 ig, g_idx[ig].x, g_idx[ig].y, g_idx[ig].z,
                                 cached.real(), cached.imag(),
                                 direct.real(), direct.imag());
                    return 1;
                }
            }
            std::println("PASS: Separable mode matches None mode");
        }

        // 4. |S| = 1.0 for all (Separable mode)
        {
            const auto& g_idx = kv.g_idx;
            int ng = static_cast<int>(g_idx.size());
            vector3d<double> tau{0.125, 0.25, 0.375};
            StructureFactor sf(tau, n1, n2, n3, StructureFactor::CacheMode::Separable);
            vector3d<double> k{0.125, 0.25, 0.375};
            for (int ig = 0; ig < std::min(ng, 50); ++ig) {
                auto val = sf(g_idx[ig], k);
                double mag = std::abs(val);
                if (std::abs(mag - 1.0) > eps) {
                    std::println("FAIL: |S| != 1.0 at ig={}, got {}", ig, mag);
                    return 1;
                }
            }
            std::println("PASS: |S| = 1.0 for all checked G-vectors (Separable mode)");
        }

        // 5. reset_frac_atomic_position alias matches reset_tau
        {
            StructureFactor sf_a({0.0, 0.0, 0.0}, 0, 0, 0, StructureFactor::CacheMode::None);
            StructureFactor sf_b({0.0, 0.0, 0.0}, 0, 0, 0, StructureFactor::CacheMode::None);
            sf_a.reset_tau({0.2, 0.3, 0.4});
            sf_b.reset_frac_atomic_position({0.2, 0.3, 0.4});
            vector3d<int>   g{1, 2, 3};
            vector3d<double> k{0.1, 0.1, 0.1};
            auto va = sf_a(g, k);
            auto vb = sf_b(g, k);
            if (std::abs(va - vb) > eps) {
                std::println("FAIL: reset_frac_atomic_position != reset_tau");
                return 1;
            }
            std::println("PASS: reset_frac_atomic_position matches reset_tau");
        }

        // 6. reset_tau in Separable mode rebuilds cache correctly
        {
            const auto& g_idx = kv.g_idx;
            StructureFactor sf({0.0, 0.0, 0.0}, n1, n2, n3, StructureFactor::CacheMode::Separable);
            sf.reset_tau({0.125, 0.25, 0.375});

            StructureFactor sf_ref({0.125, 0.25, 0.375}, n1, n2, n3, StructureFactor::CacheMode::Separable);

            vector3d<double> k{0.1, 0.2, 0.3};
            for (int ig = 0; ig < std::min(static_cast<int>(g_idx.size()), 20); ++ig) {
                auto va = sf(g_idx[ig], k);
                auto vb = sf_ref(g_idx[ig], k);
                if (std::abs(va - vb) > eps) {
                    std::println("FAIL: reset_tau in Separable mode mismatch at ig={}", ig);
                    return 1;
                }
            }
            std::println("PASS: reset_tau rebuilds cache in Separable mode");
        }

        std::println("PASS: StructureFactor section complete");
    }

    std::println("\nAll tests passed.");
    return 0;
}
