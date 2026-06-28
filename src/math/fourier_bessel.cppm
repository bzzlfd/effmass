export module math.fourier_bessel;

// All functions return the real radial part of the Fourier-Bessel transform.
// The (-i)^l phase factor from the plane-wave expansion is intentionally
// omitted — apply it at the other complex layer.

import std;
import math.numerical;
import math.sph_bessel;

export {
    // --- High-level: complete Fourier-Bessel transform --------------------------
    //  beta(q) = 4π ∫ f(r)·r²·j_l(qr)·dr
    //  Single (l,q): self-contained, manages Bessel internally.
    auto fourierBesselTransform(
        std::span<const double> r, std::span<const double> rab,
        std::span<const double> f, int l, double q,
        SimpsonMeshType mesh_type) -> double;

    //  Batch: single l, multiple q (reuses Bessel buffer).
    auto fourierBesselTransform(
        std::span<const double> r, std::span<const double> rab,
        std::span<const double> f, int l,
        std::span<const double> qs, std::span<double> out,
        SimpsonMeshType mesh_type) -> void;

    // --- Low-level: integrand + simpson, j_l(qr) provided by caller ------------
    //  Each pair has a General form and an R1 form.
    //
    //    General                R1
    //    ──────────────────────────────────────────────────
    //    Integrate       r²·jl  IntegrateR1        r¹·jl
    //    IntegrateDeriv  r³·jl' IntegrateDerivR1   r²·jl'
    //    IntegrateSecondDeriv r⁴·jl'' IntegrateSecondDerivR1 r³·jl''
    //
    //  The R1 variant has one less explicit r factor — for callers whose f(r)
    //  already includes one factor of r (the UPF beta convention: β_UPF = r·β_phys).

    //  General — standard form:  4π ∫ f(r)·r²·j_l(qr)·dr.
    auto fourierBesselIntegrate(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  R1 — standard form:  4π ∫ f(r)·r·j_l(qr)·dr.
    auto fourierBesselIntegrateR1(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  General — first-derivative form:  4π ∫ f(r)·r³·j_l'(qr)·dr.
    auto fourierBesselIntegrateDeriv(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl_deriv,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  R1 — first-derivative form:  4π ∫ f(r)·r²·j_l'(qr)·dr.
    auto fourierBesselIntegrateDerivR1(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl_deriv,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  General — second-derivative form:  4π ∫ f(r)·r⁴·j_l''(qr)·dr.
    auto fourierBesselIntegrateSecondDeriv(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl_second_deriv,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  R1 — second-derivative form:  4π ∫ f(r)·r³·j_l''(qr)·dr.
    auto fourierBesselIntegrateSecondDerivR1(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl_second_deriv,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;
}


// =============================================================================
//  Implementation
// =============================================================================

namespace {

constexpr double FOUR_PI = 4.0 * std::numbers::pi;

} // anonymous namespace


// -----------------------------------------------------------------------------
//  fourierBesselTransform  —  single (l,q)
// -----------------------------------------------------------------------------

auto fourierBesselTransform(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> f,
    int l,
    double q,
    SimpsonMeshType mesh_type
) -> double {
    double out = 0.0;
    fourierBesselTransform(r, rab, f, l,
                           std::span<const double>(&q, 1),
                           std::span<double>(&out, 1), mesh_type);
    return out;
}


// -----------------------------------------------------------------------------
//  fourierBesselTransform  —  single l, batch q
// -----------------------------------------------------------------------------

auto fourierBesselTransform(
    std::span<const double> r,
    std::span<const double> rab,
    std::span<const double> f,
    int l,
    std::span<const double> qs,
    std::span<double> out,
    SimpsonMeshType mesh_type
) -> void {
    int nq = static_cast<int>(qs.size());
    if (static_cast<std::size_t>(nq) != out.size()) {
        throw std::invalid_argument("fourierBesselTransform: qs and out size mismatch");
    }
    if (nq == 0) return;

    int n = static_cast<int>(f.size());
    if (n == 0) {
        std::fill(out.begin(), out.end(), 0.0);
        return;
    }

    int nr = static_cast<int>(r.size());
    int nrab = static_cast<int>(rab.size());
    int nmax = std::min({n, nr, nrab});
    while (nmax > 0 && f[static_cast<std::size_t>(nmax) - 1] == 0.0) {
        --nmax;
    }
    if (nmax == 0) {
        std::fill(out.begin(), out.end(), 0.0);
        return;
    }

    auto r_sub = r.subspan(0, static_cast<std::size_t>(nmax));
    auto rab_sub = rab.subspan(0, static_cast<std::size_t>(nmax));
    auto f_sub = f.subspan(0, static_cast<std::size_t>(nmax));

    std::vector<double> integrand(nmax);
    SphericalBesselJ bessel{r_sub, qs.front()};

    for (int iq = 0; iq < nq; ++iq) {
        bessel.reset(qs[iq]);
        bessel.advance(l);
        auto jl = bessel.value();

        for (int i = 0; i < nmax; ++i)
            integrand[i] = f_sub[i] * r_sub[i] * jl[i];

        out[iq] = FOUR_PI * simpson(integrand, rab_sub, mesh_type);
    }
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrate  —  General standard form, r²·j_l(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrate(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrate: r and jl must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * r[i] * jl[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrateR1  —  R1 standard form, r·j_l(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrateR1(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrateR1: r and jl must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * jl[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrateDeriv  —  General first-derivative, r³·j_l'(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrateDeriv(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl_deriv,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl_deriv.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrateDeriv: r and jl_deriv must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * r[i] * r[i] * jl_deriv[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrateDerivR1  —  R1 first-derivative, r²·j_l'(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrateDerivR1(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl_deriv,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl_deriv.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrateDerivR1: r and jl_deriv must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * r[i] * jl_deriv[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrateSecondDeriv  —  General second-derivative, r⁴·j_l''(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrateSecondDeriv(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl_second_deriv,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl_second_deriv.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrateSecondDeriv: r and jl_second_deriv must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * r[i] * r[i] * r[i] * jl_second_deriv[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}


// -----------------------------------------------------------------------------
//  fourierBesselIntegrateSecondDerivR1  —  R1 second-derivative, r³·j_l''(qr)
// -----------------------------------------------------------------------------

auto fourierBesselIntegrateSecondDerivR1(
    std::span<const double> f,
    std::span<const double> r,
    std::span<const double> jl_second_deriv,
    std::span<const double> rab,
    SimpsonMeshType mesh_type
) -> double {
    int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    if (static_cast<int>(r.size()) < n || static_cast<int>(jl_second_deriv.size()) < n) {
        throw std::invalid_argument(
            "fourierBesselIntegrateSecondDerivR1: r and jl_second_deriv must be at least as large as f");
    }

    std::vector<double> integrand(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        integrand[i] = f[i] * r[i] * r[i] * r[i] * jl_second_deriv[i];

    return FOUR_PI * simpson(integrand, rab, mesh_type);
}
