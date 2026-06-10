module;

#include <cassert>

export module H_psi.structure_factor;

import std;
import io;
import utils.vector3d;

export {
    class StructureFactor;
}


/// Structure factor S(g,k) = exp(-i·2π·(k+g)·τ) for the nonlocal pseudopotential.
///
/// Two modes:
///   None      — operator() computes directly via complex exp (no caching).
///   Separable — operator() uses precomputed 1D phase tables of size n1×n2×n3,
///               exploiting the identity
///                 exp(-i·2π·(k+g)·τ) = exp(-i·2π·k·τ)
///                                     × exp(-i·2π·g_x·τ_x)
///                                     × exp(-i·2π·g_y·τ_y)
///                                     × exp(-i·2π·g_z·τ_z)
///
/// In Separable mode the cache stores:
///   p_x[i] = exp(-i·2π·(G_x)·τ_x)   G_x = i - n1/2 - B,  B = CACHE_BUFFER
///   (p_y, p_z similarly)
/// so that S(g,k) = exp(-i·2π·k·τ) × p_x[g.x + n1/2 + B] × p_y[g.y + n2/2 + B] × p_z[g.z + n3/2 + B]
/// with G-vector components from the centered integer convention (e.g. GKK integer view).
class StructureFactor {
public:
    enum class CacheMode { None, Separable };

    /// Construct a StructureFactor.
    /// @param tau   atomic position in fractional coordinates, each component
    ///              in (-0.5, 0.5] (the convention from GKK's inferkPoint)
    /// @param n1,n2,n3   FFT grid dimensions (used in Separable mode only)
    /// @param mode  None or Separable (default: Separable)
    explicit StructureFactor(
        vector3d<double> tau,
        int n1,
        int n2,
        int n3,
        CacheMode mode = CacheMode::Separable
    );

    /// S(g, k) = exp(-i·2π·(k+g)·τ)
    auto operator()(vector3d<int> g, vector3d<double> k) const -> std::complex<double>;

    /// Reset atomic position τ (fractional coordinates).
    /// In Separable mode this rebuilds the phase caches.
    auto reset_tau(vector3d<double> tau) -> void;

    /// Alias for reset_tau — makes the meaning of τ explicit to readers.
    auto reset_frac_atomic_position(vector3d<double> fc) -> void;

private:
    static constexpr int CACHE_BUFFER = 4;   // extra entries on each side of p_?_
    CacheMode mode_;
    vector3d<double> tau_{};
    int n1_{}, n2_{}, n3_{};
    std::vector<std::complex<double>> p_x_, p_y_, p_z_;

    /// Map G integer component (centered convention) to cache index.
    /// Valid range: g ∈ [-n/2 - CACHE_BUFFER, n/2 - 1 + CACHE_BUFFER].
    static auto cacheIndex(int g, int n) -> std::size_t {
        int idx = g + n / 2 + CACHE_BUFFER;
        assert(idx >= 0);
        return static_cast<std::size_t>(idx);
    }

    auto buildCache() -> void;
};


// =============================================================================
//  Implementation
// =============================================================================

StructureFactor::StructureFactor(
    vector3d<double> tau,
    int n1, int n2, int n3,
    CacheMode mode
) : mode_(mode), tau_(tau), n1_(n1), n2_(n2), n3_(n3) {
    if (mode_ == CacheMode::Separable) {
        buildCache();
    }
}

auto StructureFactor::operator()(
    vector3d<int> g,
    vector3d<double> k) const -> std::complex<double>
{
    if (mode_ == CacheMode::None) {
        double arg = -2.0 * std::numbers::pi
                   * ((k.x + g.x) * tau_.x
                    + (k.y + g.y) * tau_.y
                    + (k.z + g.z) * tau_.z);
        return std::exp(std::complex<double>(0, arg));
    }

    // Separable mode:
    //   exp(-i·2π·k·τ) × p_x[g.x + n1/2 + B] × p_y[g.y + n2/2 + B] × p_z[g.z + n3/2 + B]
    //   where B = CACHE_BUFFER.
    //
    // Precondition: g components stay within the buffered range,
    //   g ∈ [-n/2 - B, n/2 - 1 + B], guaranteed by the caller.
    double kd = -2.0 * std::numbers::pi
              * (k.x * tau_.x + k.y * tau_.y + k.z * tau_.z);
    auto kp = std::exp(std::complex<double>(0, kd));

    auto ix = cacheIndex(g.x, n1_);
    auto iy = cacheIndex(g.y, n2_);
    auto iz = cacheIndex(g.z, n3_);

    assert(ix < p_x_.size());
    assert(iy < p_y_.size());
    assert(iz < p_z_.size());

    return kp * p_x_[ix] * p_y_[iy] * p_z_[iz];
}

auto StructureFactor::reset_tau(vector3d<double> tau) -> void {
    tau_ = tau;
    if (mode_ == CacheMode::Separable) {
        buildCache();
    }
}

auto StructureFactor::reset_frac_atomic_position(vector3d<double> fc) -> void {
    reset_tau(fc);
}

auto StructureFactor::buildCache() -> void {
    auto buf = static_cast<int>(CACHE_BUFFER);
    p_x_.resize(static_cast<std::size_t>(n1_) + CACHE_BUFFER * 2);
    p_y_.resize(static_cast<std::size_t>(n2_) + CACHE_BUFFER * 2);
    p_z_.resize(static_cast<std::size_t>(n3_) + CACHE_BUFFER * 2);

    int shift_x = n1_ / 2;
    int shift_y = n2_ / 2;
    int shift_z = n3_ / 2;

    // Iterate over the full G range (including buffer) and use cacheIndex
    // directly, so that build and lookup share the identical mapping.
    int g_min = -shift_x - buf;
    int g_max = n1_ - shift_x + buf;   // exclusive
    for (int g = g_min; g < g_max; ++g)
        p_x_[cacheIndex(g, n1_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * g * tau_.x));
    for (int g = -shift_y - buf; g < n2_ - shift_y + buf; ++g)
        p_y_[cacheIndex(g, n2_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * g * tau_.y));
    for (int g = -shift_z - buf; g < n3_ - shift_z + buf; ++g)
        p_z_[cacheIndex(g, n3_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * g * tau_.z));
}
