module;

export module H_psi.hamiltonian;

import std;

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
//  Hamiltonian  —  shared data container that owns nothing
//
//  Three groups (peers at the Hamiltonian level):
//    1.  Callable          —  H|ψ⟩             H.at_k()
//    2.  Gradient                              H.gradient
//        2.1.  Callable    —  ∂H/∂k_α|ψ⟩       H.gradient.at_k_a()
//    3.  Hessian                               H.hessian
//        3.1.  Callable    —  ∂²H/∂k_α∂k_β|ψ⟩  H.hessian.at_k_ab()
// ===========================================================================
export {
    class Hamiltonian {
    public:
        Hamiltonian() = default;

        class Callable;
        auto at_k(int ikpt) const -> Callable { return Callable(this, ikpt); }

        // -------------------------------------------------------------------
        //  1.  Callable  —  H|ψ⟩ at fixed k
        //      H.at_k(ikpt) → Callable
        //      callable(psi, hpsi)  or  callable(n_bands, psi, hpsi)
        // --------------------------------------------------------------------
        class Callable {
        public:
            Callable() = default;

            void operator()(std::span<const std::complex<double>> psi,
                            std::span<std::complex<double>> hpsi) const {}

            void operator()(int n_bands,
                            std::span<const std::complex<double>> psi,
                            std::span<std::complex<double>> hpsi) const {}

            auto dim() const -> int { return 0; }

        private:
            friend class Hamiltonian;
            Callable(const Hamiltonian* parent, int ikpt)
                : parent_(parent), ikpt_(ikpt) {}

            const Hamiltonian* parent_{nullptr};
            int ikpt_{0};
        };

        // ===================================================================
        //  2.  Gradient  —  ∂H/∂k_α|ψ⟩ at fixed k
        //      H.gradient.at_k_a(ikpt, dir) → Callable
        //      callable(psi, out)
        // ===================================================================
        class Gradient {
        public:
            Gradient() = default;

            class Callable;
            auto at_k_a(int ikpt, KDir dir) const -> Callable { return Callable(parent_, ikpt, dir); }
            
            // --------
            // Callable  (nested functor)  
            // --------
            class Callable {
            public:
                Callable() = default;

                void operator()(std::span<const std::complex<double>> psi,
                                std::span<std::complex<double>> out) const {}

                auto dim() const -> int { return 0; }

            private:
                friend class Gradient;
                Callable(const Hamiltonian* parent, int ikpt, KDir dir)
                    : parent_(parent), ikpt_(ikpt), dir_(dir) {}

                const Hamiltonian* parent_{nullptr};
                int ikpt_{0};
                KDir dir_{KDir::X};
            };

        private:
            friend class Hamiltonian;
            explicit Gradient(const Hamiltonian* parent) : parent_(parent) {}
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

            class Callable;
            auto at_k_ab(int ikpt, KDir d1, KDir d2) const -> Callable { return Callable(parent_, ikpt, d1, d2); }
            
            // --------
            // Callable  (nested functor)  
            // --------
            class Callable {
            public:
                Callable() = default;

                void operator()(std::span<const std::complex<double>> psi,
                                std::span<std::complex<double>> out) const {}

                auto dim() const -> int { return 0; }

            private:
                friend class Hessian;
                Callable(const Hamiltonian* parent, int ikpt,
                         KDir d1, KDir d2)
                    : parent_(parent), ikpt_(ikpt), d1_(d1), d2_(d2) {}

                const Hamiltonian* parent_{nullptr};
                int ikpt_{0};
                KDir d1_{KDir::X};
                KDir d2_{KDir::X};
            };

        private:
            friend class Hamiltonian;
            explicit Hessian(const Hamiltonian* parent) : parent_(parent) {}
            const Hamiltonian* parent_{nullptr};
        };
        auto hessian() const -> Hessian { return Hessian(this); }

    private:
        // Data references — filled when real types are available
        // const GKK& gkk_;
        // const FFT3D& fft_;
        // std::span<const double> vlocal_;
    };
}
