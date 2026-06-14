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
    //  Fourier-Bessel integral with single radial power: 4π ∫ f(r)·r·j_l(qr)·dr.
    //  For callers whose f(r) already includes one factor of r (e.g. UPF beta
    //  projectors are stored as r·β(r)).  Callers needing the standard r² form
    //  should absorb the extra r into f before calling.
    auto fourierBesselIntegrateR1(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl,
        std::span<const double> rab,
        SimpsonMeshType mesh_type) -> double;

    //  Standard form:  4π ∫ f(r)·r²·j_l(qr)·dr.
    auto fourierBesselIntegrate(
        std::span<const double> f,
        std::span<const double> r,
        std::span<const double> jl,
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
//  fourierBesselIntegrateR1  —  integrand + simpson, j_l(qr) provided by caller
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
//  fourierBesselIntegrate  —  standard form, j_l(qr) provided by caller
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
