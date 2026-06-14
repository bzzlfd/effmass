module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Hessian  —  ∂²H/∂k_α∂k_β|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Hessian::Hessian(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Hessian::Callable::Callable(
    const Hamiltonian* parent, int ikpt, KDir d1, KDir d2)
    : parent_(parent), ikpt_(ikpt), d1_(d1), d2_(d2)
{
    if (!parent_->gkk_)  throw std::runtime_error("hessian H|ψ⟩ requires GKK");
    if (!parent_->wg_)   throw std::runtime_error("hessian H|ψ⟩ requires WG");
    if (!parent_->vr_)   throw std::runtime_error("hessian H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("hessian H|ψ⟩ requires ATOM");
    if (parent_->elements_.empty()) throw std::runtime_error("hessian H|ψ⟩ requires NCPPs");
}

void Hamiltonian::Hessian::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> out) const {}

auto Hamiltonian::Hessian::Callable::dim() const -> int {
    return 0;
}
