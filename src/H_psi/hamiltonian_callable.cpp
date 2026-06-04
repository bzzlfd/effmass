module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Callable  —  H|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Callable::Callable(const Hamiltonian* parent, int ikpt)
    : parent_(parent), ikpt_(ikpt) {}

void Hamiltonian::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> hpsi) const {}

void Hamiltonian::Callable::operator()(
    int n_bands,
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> hpsi) const {}

auto Hamiltonian::Callable::dim() const -> int {
    return 0;
}
