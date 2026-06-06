module H_psi.hamiltonian;

import math;

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

    // Enable Integer view alongside Cartesian — needed for FFT grid placement
    // in the local-potential term.  Keep it active for the lifetime of the
    // Hamiltonian (the memory overhead is modest).
    parent_->gkk().setDataView(KVecsView::Cartesian | KVecsView::Integer);
    const auto& kv = parent_->gkk().loadKPoint(ikpt_);
    parent_->gkk().validateKineticConsistency();

    ng_ = static_cast<int>(kv.kinetic.size());
    n1_ = parent_->gkk().meta.n1;
    n2_ = parent_->gkk().meta.n2;
    n3_ = parent_->gkk().meta.n3;
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

    // --- local-potential contribution ------------------------------------------
    //  ψ(r) = Σ ψ[G] exp(iG·r)
    //  V_loc|ψ⟩ = FFT[ V_loc(r) · ψ(r) ]
    //
    //  1. Place ψ[G] onto the FFT grid at the correct G-vector positions
    //  2. FFT G2R  (plane-wave → real-space grid)
    //  3. Multiply by VR(r) in real space
    //  4. FFT R2G  (real-space grid → plane-wave coefficients)
    //  5. Extract the relevant G-vector coefficients into hpsi
    // ---------------------------------------------------------------------------
    auto n123 = static_cast<std::size_t>(n1_) * n2_ * n3_;
    std::vector<std::complex<double>> grid(n123, 0.0);

    for (int ig = 0; ig < ng_; ++ig) {
        auto i = ((kv.g_idx[ig].x % n1_) + n1_) % n1_;
        auto j = ((kv.g_idx[ig].y % n2_) + n2_) % n2_;
        auto k = ((kv.g_idx[ig].z % n3_) + n3_) % n3_;
        auto idx = static_cast<std::size_t>(i) * n2_ * n3_
                 + static_cast<std::size_t>(j) * n3_
                 + static_cast<std::size_t>(k);
        grid[idx] = psi[ig];
    }

    FFT3D fft(n1_, n2_, n3_);
    fft(grid, G2R);

    const auto& vr = parent_->vr();
    for (int i = 0; i < n1_; ++i) {
        for (int j = 0; j < n2_; ++j) {
            auto base = static_cast<std::size_t>(i) * n2_ * n3_
                      + static_cast<std::size_t>(j) * n3_;
            for (int k = 0; k < n3_; ++k) {
                grid[base + k] *= vr[i, j, k];
            }
        }
    }

    fft(grid, R2G);

    for (int ig = 0; ig < ng_; ++ig) {
        auto i = ((kv.g_idx[ig].x % n1_) + n1_) % n1_;
        auto j = ((kv.g_idx[ig].y % n2_) + n2_) % n2_;
        auto k = ((kv.g_idx[ig].z % n3_) + n3_) % n3_;
        auto idx = static_cast<std::size_t>(i) * n2_ * n3_
                 + static_cast<std::size_t>(j) * n3_
                 + static_cast<std::size_t>(k);
        hpsi[ig] += grid[idx];
    }
}
