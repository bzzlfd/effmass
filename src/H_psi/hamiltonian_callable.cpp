module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Callable  —  H|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Callable::Callable(const Hamiltonian* parent, int ikpt)
    : parent_(parent), ikpt_(ikpt)
{
    if (!parent_->gkk_)  throw std::runtime_error("H|ψ⟩ requires GKK");
    if (!parent_->wg_)   throw std::runtime_error("H|ψ⟩ requires WG");
    if (!parent_->vr_)   throw std::runtime_error("H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("H|ψ⟩ requires ATOM");
    if (parent_->ncpps_.empty()) throw std::runtime_error("H|ψ⟩ requires NCPPs");
}

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
