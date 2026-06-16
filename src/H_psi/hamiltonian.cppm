module;

export module H_psi.hamiltonian;

import std;
import io;
import pseudo;
import math;

// ===================================================================
//  ExtendedCheck  —  select which heavyweight checks to run in
//                    checkConsistencyExtended()
// ===================================================================
export {
    enum class ExtendedCheck : std::uint64_t {
        RHOReconstruct   = 1ull << 0,  // Σ occ·|WG|² → FFT → RHO'  vs  file RHO
        ValenceCount     = 1ull << 1,  // Σ(NCPP.z_valence × count)  ≈  ∫RHO d³r
        NCPPAtomCoverage = 1ull << 2,  // every ATOM species has a matching NCPP
    };
}

// ===================================================================
//  KDir — Cartesian direction specifier for gradient/hessian
// ===================================================================
export {
    enum class KDir : int {
        X = 0,
        Y = 1,
        Z = 2,
    };
}

// ===========================================================================
//  Hamiltonian  —  owns input data files for H|ψ⟩ / gradient / hessian
//
//  Data is loaded incrementally via loadXxx() methods.  Each load triggers
//  cross-file consistency checks that log what was verified.
//
//  Always-needed:  GKK, ATOM, NCPPs
//  Optional:       VR, RHO
//                  EIGEN, WG
//
//  k-point coordinates always come from GKK (inferCurrent_k).
//
//  Three groups (peers at the Hamiltonian level):
//    1.  Callable          —  H|ψ⟩             H.at_k()
//    2.  Gradient                              H.gradient
//        2.1.  Callable    —  ∂H/∂k_α|ψ⟩       H.gradient.at_k_a()
//    3.  Hessian                               H.hessian
//        3.1.  Callable    —  ∂²H/∂k_α∂k_β|ψ⟩  H.hessian.at_k_ab()
//
//  Implementations live in sibling .cpp files (module implementation units).
// ===========================================================================
export {
    class Hamiltonian {
    public:
        /// Construct with a base directory (default: current working directory).
        /// All relative paths in loadXxx() are resolved against this directory.
        explicit Hamiltonian(std::filesystem::path base_dir = ".");

        // --- step-by-step loading (relative paths resolved against base_dir_) ---
        auto loadGKK(const std::string& path) -> void;
        auto loadWG(const std::string& path) -> void;
        auto loadVR(const std::string& path) -> void;
        auto loadRHO(const std::string& path) -> void;
        auto loadATOM(const std::string& path) -> void;
        auto loadEIGEN(const std::string& path) -> void;
        auto loadOCC(const std::string& path) -> void;
        auto loadNCPP(const std::string& path) -> void;

        /// Post-load initialisation.  Prepares NCPP data (sortByL,
        /// diagonalizeNonlocal), constructs BetaqTables from NCPP data and
        /// GKK cell volume, then runs the selected Part‑3 consistency checks.
        /// Must be called after all loadXxx() and before
        /// at_k() / gradient() / hessian().
        auto finalize(std::initializer_list<ExtendedCheck> checks
                      = {ExtendedCheck::NCPPAtomCoverage},
                      bool enable_psp_nonlocal = true) -> void;

        /// Convenience: load all standard files from the base directory.
        auto loadFromDirectory() -> void;

        // --- queries ---
        auto hasGKK()   const -> bool { return gkk_.has_value(); }
        auto hasWG()    const -> bool { return wg_.has_value(); }
        auto hasVR()    const -> bool { return vr_.has_value(); }
        auto hasRHO()   const -> bool { return rho_.has_value(); }
        auto hasATOM()  const -> bool { return atom_.has_value(); }
        auto hasEIGEN() const -> bool { return eigen_.has_value(); }
        auto hasOCC()   const -> bool { return occ_.has_value(); }

        // --- owned data access (throws if not loaded) ---
        auto gkk()    const -> GKK&;
        auto wg()     const -> WG&;
        auto vr()     const -> const VR&;
        auto rho()    const -> const RHO&;
        auto atom()   const -> const ATOM&;
        auto eigen()  const -> const EIGEN&;
        auto occ()    const -> const OCC&;
        auto ncpp(int atomic_number) const -> const NCPP&;
        auto betaqTables(int atomic_number) const -> const BetaqTables&;

        // -------------------------------------------------------------------
        //  1.  Callable  —  H|ψ⟩ at fixed k
        //      H.at_k(ikpt) → Callable
        //      callable(psi, hpsi)  or  callable(n_bands, psi, hpsi)
        //      callable.set_ikpt(ikpt)  —  rebind to a different k-point
        // -------------------------------------------------------------------
        class Callable {
        public:
            Callable() = default;

            void operator()(std::span<const std::complex<double>> psi,
                            std::span<std::complex<double>> hpsi) const;

            void operator()(int n_bands,
                            std::span<const std::complex<double>> psi,
                            std::span<std::complex<double>> hpsi) const;

            auto dim() const -> int;

            /// Rebind this Callable to a different k-point.
            /// Reloads k-point data and reconstructs Ylm.
            auto set_ikpt(int ikpt) -> void;

        private:
            friend class Hamiltonian;
            Callable(const Hamiltonian* parent, int ikpt);

            const Hamiltonian* parent_{nullptr};
            int ikpt_{0};
            int ng_{0};
            int n1_{0}, n2_{0}, n3_{0};             // FFT grid dimensions

            // Ylm — cached Y_lm(theta, phi) for bound k-point.
            // Pre-allocated to parent_->ng_max_ via Ylm::reserveNg() so that
            // set_ikpt() can switch k-points without heap reallocation.
            mutable std::optional<RealSphericalHarmonics> ylm_;
        };
        auto at_k(int ikpt) const -> Callable { return Callable(this, ikpt); }


        // ===================================================================
        //  2.  Gradient  —  ∂H/∂k_α|ψ⟩ at fixed k
        //      H.gradient.at_k_a(ikpt, dir) → Callable
        //      callable(psi, out)
        // ===================================================================
        class Gradient {
        public:
            Gradient() = default;

            class Callable {
            public:
                Callable() = default;

                void operator()(std::span<const std::complex<double>> psi,
                                std::span<std::complex<double>> out) const;

                auto dim() const -> int;

            private:
                friend class Gradient;
                Callable(const Hamiltonian* parent, int ikpt, KDir dir);

                const Hamiltonian* parent_{nullptr};
                int ikpt_{0};
                KDir dir_{KDir::X};
            };
            auto at_k_a(int ikpt, KDir dir) const -> Callable { return Callable(parent_, ikpt, dir); }

        private:
            friend class Hamiltonian;
            explicit Gradient(const Hamiltonian* parent);

            const Hamiltonian* parent_{nullptr};
        };
        auto gradient() const -> Gradient { return Gradient(this); }


        // ===================================================================
        //  3.  Hessian  —  ∂²H/∂k_α∂k_β|ψ⟩ at fixed k
        //      H.hessian.at_k_ab(ikpt, d1, d2) → Callable
        //      callable(psi, out)
        // ===================================================================
        class Hessian {
        public:
            Hessian() = default;

            class Callable {
            public:
                Callable() = default;

                void operator()(std::span<const std::complex<double>> psi,
                                std::span<std::complex<double>> out) const;

                auto dim() const -> int;

            private:
                friend class Hessian;
                Callable(const Hamiltonian* parent, int ikpt,
                         KDir d1, KDir d2);

                const Hamiltonian* parent_{nullptr};
                int ikpt_{0};
                KDir d1_{KDir::X};
                KDir d2_{KDir::X};
            };
            auto at_k_ab(int ikpt, KDir d1, KDir d2) const -> Callable { return Callable(parent_, ikpt, d1, d2); }

        private:
            friend class Hamiltonian;
            explicit Hessian(const Hamiltonian* parent);

            const Hamiltonian* parent_{nullptr};
        };
        auto hessian() const -> Hessian { return Hessian(this); }

        // -- Part 3 (explicit user invocation) --
        auto checkConsistencyExtended() -> void;
        auto checkConsistencyExtended(std::initializer_list<ExtendedCheck> checks) -> void;

    private:
        std::filesystem::path base_dir_;
        mutable std::optional<GKK>   gkk_;
        mutable std::optional<WG>    wg_;
        std::optional<VR>    vr_;
        std::optional<RHO>   rho_;
        std::optional<ATOM>  atom_;
        std::optional<EIGEN> eigen_;
        std::optional<OCC>   occ_;
        struct ElementData { NCPP ncpp; std::optional<BetaqTables> betaq_tables; };
        std::vector<ElementData> elements_;

        /// Bitmask to run each Part 2 file-pair integrity check exactly once.
        /// Reset by loadXxx() when a file is loaded/reloaded so affected pairs
        /// are re-checked.
        std::uint64_t part2_done_{0};

        // -------------------------------------------------------------------
        //  Part 1  —  canonical physical quantities (quantity-centered check)
        //  First-loaded file sets the value; subsequent files must match.
        //  These are available for H|ψ⟩ (Callable / gradient / hessian).
        // -------------------------------------------------------------------
        struct FFTGrid { int n1, n2, n3; };

        std::optional<Lattice>                   canonical_lattice_;
        std::optional<int>                       canonical_nkpt_;
        std::optional<int>                       canonical_nband_;
        std::optional<FFTGrid>                   canonical_fft_grid_;
        std::optional<int>                       canonical_is_SO_;
        std::optional<int>                       canonical_islda_;
        std::optional<double>                    canonical_Ecut_;
        std::optional<int>                       canonical_natom_;
        std::optional<std::vector<std::array<double, 3>>> canonical_kpt_vec_;

        // -- Derived quantities (computed in finalize) --
        int ng_max_{};          // maximum number of G-vectors across all k-points (from GKK)
        int l_max_{};           // maximum angular momentum across all NCPP elements
        bool enable_psp_nonlocal_{true};  // whether nonlocal PSP was prepared in finalize()
        std::optional<FFT3D> fft_;       // FFT plan (constructed in finalize)

        /// Resolve a user-provided path against base_dir_.
        /// Absolute paths are returned as-is; relative paths are prefixed with base_dir_.
        auto resolve(const std::string& path) const -> std::string;

        auto checkConsistency() -> void;
        auto checkPart1() -> void;
        auto checkPart2() -> void;
    };
}
