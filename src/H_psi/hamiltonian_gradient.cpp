module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Gradient  —  ∂H/∂k_α|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Gradient::Gradient(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Gradient::Callable::Callable(
    const Hamiltonian* parent, int ikpt, KDir dir)
    : parent_(parent), ikpt_(ikpt), dir_(dir)
{
    if (!parent_->gkk_)  throw std::runtime_error("gradient H|ψ⟩ requires GKK");
    if (!parent_->wg_)   throw std::runtime_error("gradient H|ψ⟩ requires WG");
    if (!parent_->vr_)   throw std::runtime_error("gradient H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("gradient H|ψ⟩ requires ATOM");
    if (parent_->ncpps_.empty()) throw std::runtime_error("gradient H|ψ⟩ requires NCPPs");
}

void Hamiltonian::Gradient::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> out) const {}

auto Hamiltonian::Gradient::Callable::dim() const -> int {
    return 0;
}
