module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Callable  —  H|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Callable::Callable(const Hamiltonian* parent, int ikpt)
    : parent_(parent), ikpt_(ikpt)
{
    if (!parent_->gkk_)  throw std::runtime_error("H|ψ⟩ requires GKK");
    if (!parent_->vr_)   throw std::runtime_error("H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("H|ψ⟩ requires ATOM");
    if (parent_->ncpps_.empty()) throw std::runtime_error("H|ψ⟩ requires NCPPs");

    const auto& kv = parent_->gkk().loadKPoint(ikpt_);
    parent_->gkk().validateKineticConsistency();
    ng_ = static_cast<int>(kv.kinetic.size());
}

// -----------------------------------------------------------------------------
//  dim  —  number of G-vectors for the bound k-point
// -----------------------------------------------------------------------------
auto Hamiltonian::Callable::dim() const -> int {
    return ng_;
}

// -----------------------------------------------------------------------------
//  operator()  —  H|ψ⟩ at fixed k
//
//  The caller is responsible for zeroing hpsi before calling, or for ensuring
//  that other terms (local potential, non-local pseudopotential) have already
//  been accumulated.
// -----------------------------------------------------------------------------
void Hamiltonian::Callable::operator()(
    std::span<const std::complex<double>> psi,
    std::span<std::complex<double>> hpsi) const
{
    const auto& kv = parent_->gkk().loadKPoint(ikpt_);

    // --- kinetic-energy contribution ------------------------------------------
    //  T|ψ⟩[ig]  +=  kinetic[ig] * ψ[ig]
    //  kinetic[ig] = |G+k|²/2, read from OUT.GKK.  Since T is diagonal in
    //  G-space, each plane-wave coefficient is multiplied independently.
    // ---------------------------------------------------------------------------
    for (int ig = 0; ig < ng_; ++ig) {
        hpsi[ig] += kv.kinetic[ig] * psi[ig];
    }
}
