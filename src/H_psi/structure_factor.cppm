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
///   p_x[i] = exp(-i·2π·(i - n1/2)·τ_x)   i = 0 … n1-1
///   (p_y, p_z similarly)
/// so that S(g,k) = exp(-i·2π·k·τ) × p_x[g.x + n1/2] × p_y[g.y + n2/2] × p_z[g.z + n3/2]
/// with G-vector components from the centered integer convention (e.g. GKK integer view).
class StructureFactor {
public:
    enum class CacheMode { None, Separable };

    /// Construct a StructureFactor.
    /// @param tau   atomic position in fractional coordinates
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
    CacheMode mode_;
    vector3d<double> tau_{};
    int n1_{}, n2_{}, n3_{};
    std::vector<std::complex<double>> p_x_, p_y_, p_z_;

    /// Map G integer component (centered convention, [-n/2, n/2-1]) to cache index [0, n-1].
    static auto cacheIndex(int g, int n) -> std::size_t {
        return static_cast<std::size_t>(g + n / 2);
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
    //   exp(-i·2π·k·τ) × p_x[g.x + n1/2] × p_y[g.y + n2/2] × p_z[g.z + n3/2]
    //
    // Precondition: FFT grid guarantees g components in [-n/2, n/2-1], so
    //   g.x + n1/2 ∈ [0, n1-1]  (and similarly for y, z).
    // See H_psi/structure_factor.cppm design notes for proof.
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
    p_x_.resize(static_cast<std::size_t>(n1_));
    p_y_.resize(static_cast<std::size_t>(n2_));
    p_z_.resize(static_cast<std::size_t>(n3_));

    int shift_x = n1_ / 2;
    int shift_y = n2_ / 2;
    int shift_z = n3_ / 2;

    // p_x[i] = exp(-i·2π·(i - n1/2)·τ_x)
    // Uses the SAME n_/2 (integer) as cacheIndex(), so that lookup and build
    // refer to the same centered-convention G:  index = g + n/2 ↔ G = i - n/2.
    for (int i = 0; i < n1_; ++i)
        p_x_[static_cast<std::size_t>(i)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * (i - shift_x) * tau_.x));
    for (int i = 0; i < n2_; ++i)
        p_y_[static_cast<std::size_t>(i)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * (i - shift_y) * tau_.y));
    for (int i = 0; i < n3_; ++i)
        p_z_[static_cast<std::size_t>(i)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * (i - shift_z) * tau_.z));
}
