module H_psi.hamiltonian;

// =============================================================================
//  Hamiltonian::Hessian  —  ∂²H/∂k_α∂k_β|ψ⟩ at fixed k
// =============================================================================

Hamiltonian::Hessian::Hessian(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Hessian::Callable::Callable(
    const Hamiltonian* parent, int ikpt)
    : parent_(parent), ikpt_(ikpt)
{
    if (!parent_->gkk_)  throw std::runtime_error("hessian H|ψ⟩ requires GKK");
    if (!parent_->wg_)   throw std::runtime_error("hessian H|ψ⟩ requires WG");
    if (!parent_->vr_)   throw std::runtime_error("hessian H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("hessian H|ψ⟩ requires ATOM");

    // Validate NCPP, BetaqTables, DBetaqTables, and D2BetaqTables for all
    // element types.  Catches missing/broken element data at construction time,
    // so operator() can use bare access without try/catch noise.
    if (parent_->psp_features_ & static_cast<std::uint64_t>(PSPFeature::Nonlocal)) {
        constexpr auto DBETAQ  = static_cast<std::uint64_t>(PSPFeature::DBetaq);
        constexpr auto D2BETAQ = static_cast<std::uint64_t>(PSPFeature::D2Betaq);
        for (auto&& type : parent_->atom().eachType()) {
            try {
                const auto& ncp = parent_->ncpp(type.z);
                if (ncp.meta.l_max >= 0) {
                    parent_->betaqTables(type.z);
                    if (parent_->psp_features_ & DBETAQ)
                        parent_->dbetaqTables(type.z);
                    if (parent_->psp_features_ & D2BETAQ)
                        parent_->d2betaqTables(type.z);
                }
            } catch (const std::exception&) {
                std::throw_with_nested(std::runtime_error(
                    std::format("hess H|ψ⟩ (ikpt={}, element={})",
                        ikpt_, ATOM::elementName(type.z))));
            }
        }
    }

    set_ikpt(ikpt);
}

void Hamiltonian::Hessian::Callable::operator()(
    std::span<const std::complex<double>> psi,
    HessHPsi out) const {}

auto Hamiltonian::Hessian::Callable::dim() const -> int {
    return ng_;
}

// -----------------------------------------------------------------------------
//  set_ikpt  —  rebind to a different k-point
//
//  Reloads k-point data so that dim() and future operator() reflect the
//  correct number of G-vectors.
// -----------------------------------------------------------------------------
auto Hamiltonian::Hessian::Callable::set_ikpt(int ikpt) -> void {
    ikpt_ = ikpt;

    parent_->gkk().setDataView(KVecsView::Cartesian | KVecsView::Integer | KVecsView::Spherical);
    const auto& kv = parent_->gkk().loadKPoint(ikpt_);
    parent_->gkk().validateKineticConsistency();

    ng_ = static_cast<int>(kv.kinetic.size());
}
