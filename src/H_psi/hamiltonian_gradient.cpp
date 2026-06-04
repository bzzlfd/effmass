module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Gradient  —  ∂H/∂k_α|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Gradient::Gradient(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Gradient::Callable::Callable(
    const Hamiltonian* parent, int ikpt, KDir dir)
    : parent_(parent), ikpt_(ikpt), dir_(dir) {}

void Hamiltonian::Gradient::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> out) const {}

auto Hamiltonian::Gradient::Callable::dim() const -> int {
    return 0;
}
