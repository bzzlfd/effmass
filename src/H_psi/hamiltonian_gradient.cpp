module H_psi.hamiltonian;

import math;
import H_psi.structure_factor;

// =============================================================================
//  Hamiltonian::Gradient  —  ∂H/∂k_α|ψ⟩ at fixed k
//
//  The Callable computes all three Cartesian gradient directions (X/Y/Z) in
//  one batch, because the nonlocal pseudopotential term shares spherical-
//  gradient computation (g_q, g_θ, g_φ) across all three directions.
// =============================================================================

Hamiltonian::Gradient::Gradient(const Hamiltonian* parent)
    : parent_(parent) {}

//  ------- nested Callable -------

Hamiltonian::Gradient::Callable::Callable(
    const Hamiltonian* parent, int ikpt)
    : parent_(parent), ikpt_(ikpt)
{
    if (!parent_->gkk_)  throw std::runtime_error("gradient H|ψ⟩ requires GKK");
    if (!parent_->wg_)   throw std::runtime_error("gradient H|ψ⟩ requires WG");
    if (!parent_->vr_)   throw std::runtime_error("gradient H|ψ⟩ requires VR");
    if (!parent_->atom_) throw std::runtime_error("gradient H|ψ⟩ requires ATOM");

    // FFT grid dimensions are k-point-independent.
    n1_ = parent_->canonical_fft_grid_->n1;
    n2_ = parent_->canonical_fft_grid_->n2;
    n3_ = parent_->canonical_fft_grid_->n3;

    // Validate NCPP, BetaqTables, and DBetaqTables for all element types.
    // Catches missing/broken element data at construction time, so operator()
    // can use bare access without try/catch noise.
    if (parent_->psp_features_ & static_cast<std::uint64_t>(PSPFeature::Nonlocal)) {
        for (auto&& type : parent_->atom().eachType()) {
            try {
                const auto& ncp = parent_->ncpp(type.z);
                if (ncp.meta.l_max >= 0) {
                    parent_->betaqTables(type.z);
                    parent_->dbetaqTables(type.z);
                }
            } catch (const std::exception&) {
                std::throw_with_nested(std::runtime_error(
                    std::format("grad H|ψ⟩ (ikpt={}, element={})",
                        ikpt_, ATOM::elementName(type.z))));
            }
        }
    }

    // Delegate k-point-dependent setup (loadKPoint, trig arrays, Ylm) to set_ikpt.
    set_ikpt(ikpt);
}


// -----------------------------------------------------------------------------
//  dim  —  number of G-vectors for the bound k-point
// -----------------------------------------------------------------------------
auto Hamiltonian::Gradient::Callable::dim() const -> int {
    return ng_;
}

// -----------------------------------------------------------------------------
//  set_ikpt  —  rebind to a different k-point
//
//  Reloads k-point data, recomputes trig arrays for in-loop Cartesian
//  rotation, and (re)initialises the Ylm engine when nonlocal projectors
//  are present.
// -----------------------------------------------------------------------------
auto Hamiltonian::Gradient::Callable::set_ikpt(int ikpt) -> void {
    ikpt_ = ikpt;

    parent_->gkk().setDataView(KVecsView::Cartesian | KVecsView::Integer | KVecsView::Spherical);
    const auto& kv = parent_->gkk().loadKPoint(ikpt_);
    parent_->gkk().validateKineticConsistency();

    ng_ = static_cast<int>(kv.kinetic.size());

    // --- Ylm engine (only when nonlocal projectors exist) --------------------
    if ((parent_->psp_features_ & static_cast<std::uint64_t>(PSPFeature::Nonlocal))
        && parent_->l_max_ >= 0)
    {
        if (!(parent_->psp_features_ & static_cast<std::uint64_t>(PSPFeature::DBetaq)))
            throw std::runtime_error("gradient H|ψ⟩ requires PSPFeature::DBetaq for nonlocal contribution");

        if (engine_) {
            engine_->reinit(kv.theta, kv.phi);
            engine_->setLMax(parent_->l_max_);
        } else {
            engine_.emplace(kv.theta, kv.phi, parent_->l_max_);
            engine_->reserveNg(parent_->ng_max_);
            ylm_data_.reserveNg(parent_->ng_max_);
            ylm_grad_theta_data_.reserveNg(parent_->ng_max_);
            ylm_grad_phi_data_.reserveNg(parent_->ng_max_);
        }
    }
}


// -----------------------------------------------------------------------------
//  operator()  —  ∂H/∂k_α|ψ⟩  at fixed k  (all three directions)
//
//  The Hamiltonian k-gradient has two non-zero contributions:
//
//    1. Kinetic:  ∂T/∂k_α = (k+G)_α = K_α   (diagonal in G-space; GKK stores K = (G+k))
//    2. Nonlocal pseudopotential:
//         ∂(V_NL)/∂k_α  involves differentiating  β(q)·Y_lm(θ,φ)·S(g,k)
//         via the product rule  ∂(⟨p|ψ⟩·|p⟩) = ⟨∂p|ψ⟩·|p⟩ + ⟨p|ψ⟩·|∂p⟩.
//
//  The local potential V_loc(r) is k-independent — no contribution.
//
//  Only directions with non-empty spans in `out` are computed; all three
//  share the direction-independent spherical gradient (g_q, g_θ, g_φ).
// -----------------------------------------------------------------------------
void Hamiltonian::Gradient::Callable::operator()(
    std::span<const std::complex<double>> psi,
    GradHPsi out) const
{
    // Guard: ensure no external code changed GKK's loaded k-point since
    // construction (which would make our precomputed trig arrays and Ylm
    // engine inconsistent with the current GKK internal data).
    if (parent_->gkk().current_ikpt() != ikpt_)
        throw std::runtime_error(
            "Gradient::Callable: GKK k-point changed externally — expected "
            + std::to_string(ikpt_) + ", got "
            + std::to_string(parent_->gkk().current_ikpt()));

    const auto& kv = parent_->gkk().currentData();

    bool need_x = !out.x.empty();
    bool need_y = !out.y.empty();
    bool need_z = !out.z.empty();

    // --- 1. kinetic-energy contribution --------------------------------------------
    //  ∂/∂k_α ( |G+k|²/2 ) = (G+k)_α = K_α
    //  (Cartesian coords are  K = (G+k).)
    // --------------------------------------------------------------------------------
    if (need_x) for (int ig = 0; ig < ng_; ++ig) out.x[ig] += kv.Kx[ig] * psi[ig];
    if (need_y) for (int ig = 0; ig < ng_; ++ig) out.y[ig] += kv.Ky[ig] * psi[ig];
    if (need_z) for (int ig = 0; ig < ng_; ++ig) out.z[ig] += kv.Kz[ig] * psi[ig];

    // --- 2. nonlocal pseudopotential contribution -----------------------------------
    //  After diagonalizeNonlocal() (called in finalize), D_ij = λ_i δ_ij, so:
    //
    //    V_NL|ψ⟩ = Σ_{a} Σ_{lm} Σ_{i} λ_i · ⟨p_i|ψ⟩ · |p_i⟩
    //
    //  where  |p_i⟩ = β_i(q) · Y_lm(θ,φ) · S(g,k).
    //
    //  The k-gradient of this term involves the product rule:
    //
    //    ∂ (V_NL|ψ⟩) / ∂k_α
    //      = Σ λ_i · ⟨∂p_i/∂k_α|ψ⟩ · |p_i⟩  +  Σ λ_i · ⟨p_i|ψ⟩ · |∂p_i/∂k_α⟩
    //
    //  Spherical gradient components (direction-independent, shared across α):
    //    grad_q     = β′(q)           · Y_lm
    //    grad_θ     = (β(q)/q)        · ∂Y_lm/∂θ
    //    grad_φ     = (β(q)/q)        · (1/sinθ)·∂Y_lm/∂φ
    //
    //  Cartesian projection (direction-dependent):
    //    ∇_α =  R_α_q · grad_q  +  R_α_θ · grad_θ  +  R_α_φ · grad_φ
    //
    //  ∂S/∂k_α terms cancel between the two product-rule branches (required
    //  by τ-translation invariance), so the projector derivative reduces to
    //  the spherical-gradient derivative alone:
    //
    //    ∂p_i/∂k_α  =  ∇_α[β_i(q)·Y_lm(θ,φ)] · S(g,k)
    // --------------------------------------------------------------------------------

    if (!(parent_->psp_features_ & static_cast<std::uint64_t>(PSPFeature::Nonlocal))) return;

    // Trig views from the Ylm engine — avoids redundant storage and sin/cos
    // evaluation (RealSphericalHarmonicsEngine precomputes these internally).
    const auto& sin_theta = engine_->sinTheta();
    const auto& cos_theta = engine_->cosTheta();
    const auto& sin_phi   = engine_->sinPhi();
    const auto& cos_phi   = engine_->cosPhi();

    const auto& atom = parent_->atom();

    for (auto&& type : atom.eachType()) {
        const auto& ncp = parent_->ncpp(type.z);
        int l_max = ncp.meta.l_max;
        if (l_max < 0) continue;

        const auto& betaq_tables  = parent_->betaqTables(type.z);
        const auto& dbetaq_tables = parent_->dbetaqTables(type.z);

        // Precompute β(q) and β′(q) per (l, ib) — τ-independent, shared across atoms.
        struct BetaGradCache { int l; int ib; double lambda; double dbeta_at_0; std::vector<double> beta_q; std::vector<double> dbeta_q; };
        std::vector<BetaGradCache> beta_cache;
        for (int l = 0; l <= l_max; ++l) {
            auto pb = ncp.projectorBlock(l);
            if (pb.B.rows <= 0) continue;

            for (int ib = 0; ib < pb.B.rows; ++ib) {
                double lambda = pb.B[ib, ib];
                if (std::abs(lambda) < 1e-15) continue;

                double dbeta_at_0 = dbetaq_tables.interpolate(l, ib, 0.0);

                std::vector<double> beta_q(static_cast<std::size_t>(ng_));
                std::vector<double> dbeta_q(static_cast<std::size_t>(ng_));
                for (int ig = 0; ig < ng_; ++ig) {
                    double q = kv.q[ig];
                    beta_q[ig]  = betaq_tables.interpolate(l, ib, q);
                    dbeta_q[ig] = dbetaq_tables.interpolate(l, ib, q);
                }

                beta_cache.push_back({l, ib, lambda, dbeta_at_0, std::move(beta_q), std::move(dbeta_q)});
            }
        }
        if (beta_cache.empty()) continue;

        // Loop over atoms of this type
        for (auto&& ad : atom.eachTypeAtom(type.ityp)) {
            auto tau = ad.coord;

            // --- Structure factor ---
            StructureFactor sf(tau, n1_, n2_, n3_);

            std::vector<std::complex<double>> sf_vals(static_cast<std::size_t>(ng_));
            for (int ig = 0; ig < ng_; ++ig) {
                sf_vals[ig] = sf(kv.g_idx[ig], kv.kPoint);
            }

            for (const auto& bc : beta_cache) {
                for (int m = -bc.l; m <= bc.l; ++m) {
                    const auto& ylm         = ylm_data_.get(*engine_, bc.l, m);
                    const auto& ylm_grad_theta  = ylm_grad_theta_data_.get(*engine_, bc.l, m);
                    const auto& ylm_grad_phi    = ylm_grad_phi_data_.get(*engine_, bc.l, m);

                    // ------------------------------------------------------------------
                    //  Per-G-vector shared computation (inlined by compiler at the
                    //  two call sites below):
                    //
                    //    |p⟩      = β·Y_lm · S
                    //    |∂p/∂k_α⟩ = ∇_α[β·Y_lm] · S          (per requested dir.)
                    //
                    // ------------------------------------------------------------------
                    struct ProjG { std::complex<double> p, dp_x, dp_y, dp_z; };
                    auto proj_at = [&](int ig) -> ProjG {
                        double q = kv.q[ig];
                        double beta_ylm = bc.beta_q[ig] * ylm[ig];
                        double dbeta_q = bc.dbeta_q[ig];

                        // β/q  with q→0 limit:  lim_{q→0} β(q)/q = β′(0)
                        double beta_over_q = (q > 1e-12) ? bc.beta_q[ig] / q : bc.dbeta_at_0;

                        // Spherical gradient components of β·Y_lm
                        double grad_q     = dbeta_q * ylm[ig];
                        double grad_theta = beta_over_q * ylm_grad_theta[ig];
                        double grad_phi   = beta_over_q * ylm_grad_phi[ig];

                        // Trig lookups for spherical→Cartesian rotation
                        double s_th = sin_theta[ig], c_th = cos_theta[ig];
                        double s_ph = sin_phi[ig],   c_ph = cos_phi[ig];

                        auto S = sf_vals[ig];
                        ProjG ret;
                        ret.p = std::complex<double>(beta_ylm, 0.0) * S;
                        if (need_x) ret.dp_x = std::complex<double>(
                                        s_th * c_ph * grad_q + c_th * c_ph * grad_theta - s_ph * grad_phi, 0.0) * S;
                        if (need_y) ret.dp_y = std::complex<double>(
                                        s_th * s_ph * grad_q + c_th * s_ph * grad_theta + c_ph * grad_phi, 0.0) * S;
                        if (need_z) ret.dp_z = std::complex<double>(
                                        c_th * grad_q - s_th * grad_theta, 0.0) * S;
                        return ret;
                    };

                    // Pass 1: inner products  ⟨p|ψ⟩  and  ⟨∂p/∂α|ψ⟩
                    std::complex<double> inner = 0.0;
                    std::complex<double> inner_dp_x = 0.0, inner_dp_y = 0.0, inner_dp_z = 0.0;
                    for (int ig = 0; ig < ng_; ++ig) {
                        auto [p, dp_x, dp_y, dp_z] = proj_at(ig);
                        inner += std::conj(p) * psi[ig];
                        if (need_x) inner_dp_x += std::conj(dp_x) * psi[ig];
                        if (need_y) inner_dp_y += std::conj(dp_y) * psi[ig];
                        if (need_z) inner_dp_z += std::conj(dp_z) * psi[ig];
                    }

                    // Pass 2:  out_α += λ · (⟨∂p|ψ⟩·|p⟩ + ⟨p|ψ⟩·|∂p_α⟩)
                    //         (product-rule structure visible in each expression)
                    for (int ig = 0; ig < ng_; ++ig) {
                        auto [p, dp_x, dp_y, dp_z] = proj_at(ig);
                        if (need_x) out.x[ig] += bc.lambda * (inner_dp_x * p + inner * dp_x);
                        if (need_y) out.y[ig] += bc.lambda * (inner_dp_y * p + inner * dp_y);
                        if (need_z) out.z[ig] += bc.lambda * (inner_dp_z * p + inner * dp_z);
                    }
                }
            }
        }
    }
}
