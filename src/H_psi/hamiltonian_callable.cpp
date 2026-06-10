module H_psi.hamiltonian;

import math;
import H_psi.structure_factor;

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

    // Enable Integer view alongside Cartesian and Spherical — needed for
    // FFT grid placement (Integer), Y_lm evaluation (Spherical theta/phi),
    // and beta_q interpolation (Spherical q = |K|).
    parent_->gkk().setDataView(KVecsView::Cartesian | KVecsView::Integer | KVecsView::Spherical);
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

    // --- nonlocal pseudopotential contribution ---------------------------------
    //  V_NL = Σ_{a} Σ_{lm} Σ_{ij} D_{ij} |β_i Y_lm S_a⟩ ⟨β_j Y_lm S_a|
    //
    //  After diagonalizeNonlocal() (called during NCPP loading), D_ij =
    //  λ_i δ_{ij} is diagonal, so each projector contributes independently:
    //
    //    H|ψ⟩ += λ_i · ⟨p_i|ψ⟩ · |p_i⟩
    //
    //  where |p_i⟩ = β_i(q) · Y_lm(θ,φ) · S(g,k,τ_a)  per G-vector.
    //
    //  Y_lm and S(g,k,τ) are independent of β_i, so we precompute S(g,k,τ)
    //  once per atom and re-evaluate β_i(q) per projector.
    // ---------------------------------------------------------------------------
    // Check whether any pseudopotential has nonlocal projectors.
    // If none, the nonlocal contribution is zero.
    int global_max_l = -1;
    for (const auto& ncpp : parent_->ncpps_) {
        if (ncpp.meta.l_max > global_max_l) global_max_l = ncpp.meta.l_max;
    }
    if (global_max_l < 0) return;

    RealSphericalHarmonics ylm(kv.theta, kv.phi, global_max_l);

    const auto& atom = parent_->atom();
    double q_max = std::sqrt(2.0 * parent_->gkk().meta.Ecut) + 2.0;
    double norm_coeff = 1.0 / std::sqrt(parent_->gkk().meta.lattice.volume());

    for (auto&& type : atom.eachType()) {
        const auto& ncp = parent_->ncpp(type.z);
        int l_max = ncp.meta.l_max;
        if (l_max < 0) continue;

        // Convert radial mesh type for BetaqInterpolator
        auto mesh_type = (ncp.mesh.type == MeshType::Uniform)
                         ? RadialMeshType::Uniform
                         : RadialMeshType::General;

        // Collect all projector blocks for this atom type so we don't
        // recompute them for every atom of the same element.
        struct BlockInfo { int l; int nb; ProjectorBlock block; };
        std::vector<BlockInfo> blocks;
        for (int l = 0; l <= l_max; ++l) {
            auto block = ncp.projectorBlock(l);
            if (block.B.rows > 0) {
                blocks.push_back({l, block.B.rows, std::move(block)});
            }
        }
        if (blocks.empty()) continue;

        // Loop over atoms of this type
        for (auto&& ad : atom.eachTypeAtom(type.ityp)) {
            auto tau = ad.coord;

            // --- Structure factor: S(g,k) = exp(-i·2π·(k+g)·τ) ---
            StructureFactor sf(tau, n1_, n2_, n3_);

            // Precompute S(g,k) for all G-vectors (independent of l,m,ib)
            std::vector<std::complex<double>> sf_vals(
                static_cast<std::size_t>(ng_));
            for (int ig = 0; ig < ng_; ++ig) {
                sf_vals[ig] = sf(kv.g_idx[ig], kv.kPoint);
            }

            for (const auto& bi : blocks) {
                // --- Beta(q) interpolation for all projectors in this block ---
                // Create one interpolator for the first projector, then
                // reset_beta for subsequent ones to reuse table memory.
                auto first_beta = bi.block.beta[0];
                BetaqInterpolator betaq_interp(
                    ncp.mesh.r, ncp.mesh.rab, first_beta,
                    bi.l, 0.01, q_max, mesh_type, norm_coeff);

                for (int ib = 0; ib < bi.nb; ++ib) {
                    if (ib > 0) {
                        betaq_interp.reset_beta(
                            ncp.mesh.r, ncp.mesh.rab, bi.block.beta[ib]);
                    }

                    double lambda = bi.block.B[ib, ib];
                    if (std::abs(lambda) < 1e-15) continue;

                    // Precompute β(q) for all G-vectors
                    std::vector<double> beta_q(static_cast<std::size_t>(ng_));
                    for (int ig = 0; ig < ng_; ++ig) {
                        beta_q[ig] = betaq_interp(kv.q[ig]);
                    }

                    for (int m = -bi.l; m <= bi.l; ++m) {
                        const auto& ylm_lm = ylm.get(bi.l, m);

                        // First pass:  inner = ⟨projector|ψ⟩
                        std::complex<double> inner = 0.0;
                        for (int ig = 0; ig < ng_; ++ig) {
                            auto p = std::complex<double>(
                                beta_q[ig] * ylm_lm[ig], 0.0) * sf_vals[ig];
                            inner += std::conj(p) * psi[ig];
                        }

                        // Second pass:  hpsi += λ · inner · |projector⟩
                        for (int ig = 0; ig < ng_; ++ig) {
                            auto p = std::complex<double>(
                                beta_q[ig] * ylm_lm[ig], 0.0) * sf_vals[ig];
                            hpsi[ig] += lambda * inner * p;
                        }
                    }
                }
            }
        }
    }
}
