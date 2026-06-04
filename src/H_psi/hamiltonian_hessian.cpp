module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Hessian  —  ∂²H/∂k_α∂k_β|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Hessian::Hessian(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Hessian::Callable::Callable(
    const Hamiltonian* parent, int ikpt, KDir d1, KDir d2)
    : parent_(parent), ikpt_(ikpt), d1_(d1), d2_(d2) {}

void Hamiltonian::Hessian::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> out) const {}

auto Hamiltonian::Hessian::Callable::dim() const -> int {
    return 0;
}
