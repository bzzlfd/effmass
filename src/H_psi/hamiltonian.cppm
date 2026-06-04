module;

export module H_psi.hamiltonian;

import std;
import io;
import pseudo;

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
//  Always-needed:  GKK, WG, ATOM, NCPPs
//  Optional:       VR, RHO
//                  EIGEN
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
        Hamiltonian() = default;

        // --- step-by-step loading ---
        auto loadGKK(const std::string& path) -> void;
        auto loadWG(const std::string& path) -> void;
        auto loadVR(const std::string& path) -> void;
        auto loadRHO(const std::string& path) -> void;
        auto loadATOM(const std::string& path) -> void;
        auto loadEIGEN(const std::string& path) -> void;
        auto loadNCPPs(const std::string& directory) -> void;

        /// Convenience: load everything from a working directory.
        auto loadFromDirectory(const std::string& directory) -> void;

        // --- queries ---
        auto hasGKK()   const -> bool { return gkk_.has_value(); }
        auto hasWG()    const -> bool { return wg_.has_value(); }
        auto hasVR()    const -> bool { return vr_.has_value(); }
        auto hasRHO()   const -> bool { return rho_.has_value(); }
        auto hasATOM()  const -> bool { return atom_.has_value(); }
        auto hasEIGEN() const -> bool { return eigen_.has_value(); }

        // --- owned data access (throws if not loaded) ---
        auto gkk()   const -> const GKK&;
        auto wg()    const -> const WG&;
        auto vr()    const -> const VR&;
        auto rho()   const -> const RHO&;
        auto atom()  const -> const ATOM&;
        auto eigen() const -> const EIGEN&;
        auto ncpp(int atomic_number) const -> const NCPP&;

        // -------------------------------------------------------------------
        //  1.  Callable  —  H|ψ⟩ at fixed k
        //      H.at_k(ikpt) → Callable
        //      callable(psi, hpsi)  or  callable(n_bands, psi, hpsi)
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

        private:
            friend class Hamiltonian;
            Callable(const Hamiltonian* parent, int ikpt);

            const Hamiltonian* parent_{nullptr};
            int ikpt_{0};
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

    private:
        std::optional<GKK>   gkk_;
        std::optional<WG>    wg_;
        std::optional<VR>    vr_;
        std::optional<RHO>   rho_;
        std::optional<ATOM>  atom_;
        std::optional<EIGEN> eigen_;
        std::vector<NCPP>    ncpps_;

        auto checkConsistency() const -> void;
    };
}
