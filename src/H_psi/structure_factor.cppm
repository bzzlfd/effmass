module;

export module H_psi.structure_factor;

import std;
import io;
import utils.vector3d;

export {

/// Structure factor for the nonlocal pseudopotential.
///
///   S(G) = exp(-i(k+G) · τ)
///
/// τ is in fractional coordinates (set via set_tau).
/// k is in fractional coordinates (passed to operator()).
/// G is the plane-wave integer triplet (g_x, g_y, g_z).
///
/// Three access paths:
///   - Non-cached:  operator()(g, k)      — 1 complex exp per call
///   - Cached:      prepare(gs) then
///                  operator()(ig, k)      — 1 complex exp + 3 table lookups
///   - Batch:       operator()(gs, k, out) — N lookups, no persistent cache
class StructureFactor {
public:
    StructureFactor() = default;

    /// Set atomic position τ in fractional coordinates.
    auto set_tau(std::span<const double, 3> frac_coord) -> void;
    auto set_frac_atomic_position(std::span<const double, 3> frac_coord) -> void;

    // ---- non-cached (direct complex exp) ----

    /// S(g, k) = exp(-i(k+G)·τ)
    auto operator()(const vector3d<int>& g,
                    const vector3d<double>& k) const -> std::complex<double>;

    // ---- cached (after prepare) ----

    /// Precompute component-wise phase tables for a set of G-vectors.
    /// Subsequent operator()(ig, k) will use the table.
    auto prepare(std::span<const vector3d<int>> gs) -> void;

    /// S(ig, k) using precomputed tables.  Requires prepare() first.
    auto operator()(int ig, const vector3d<double>& k) const -> std::complex<double>;

    // ---- batch (self-contained) ----

    /// fill out[i] = S(gs[i], k); builds internal tables on the fly.
    auto operator()(std::span<const vector3d<int>> gs,
                    const vector3d<double>& k,
                    std::span<std::complex<double>> out) const -> void;

private:
    vector3d<double> tau_{};

    // cache state (valid after prepare())
    std::vector<vector3d<int>>   g_idx_;
    std::vector<std::complex<double>> p_x_, p_y_, p_z_;
    int gx_offset_{0}, gy_offset_{0}, gz_offset_{0};
};

} // export


// =============================================================================
//  Implementation
// =============================================================================

auto StructureFactor::set_tau(std::span<const double, 3> fc) -> void {
    tau_.x = fc[0];
    tau_.y = fc[1];
    tau_.z = fc[2];
}

auto StructureFactor::set_frac_atomic_position(std::span<const double, 3> fc) -> void {
    set_tau(fc);
}

auto StructureFactor::operator()(
    const vector3d<int>& g,
    const vector3d<double>& k) const -> std::complex<double>
{
    double arg = -2.0 * std::numbers::pi
               * ((k.x + g.x) * tau_.x
                + (k.y + g.y) * tau_.y
                + (k.z + g.z) * tau_.z);
    return std::exp(std::complex<double>(0, arg));
}

auto StructureFactor::operator()(
    int ig,
    const vector3d<double>& k) const -> std::complex<double>
{
    // k-phase (one complex exp).
    double kd = -2.0 * std::numbers::pi
              * (k.x * tau_.x + k.y * tau_.y + k.z * tau_.z);
    auto kp = std::exp(std::complex<double>(0, kd));

    // G-phase: three table lookups.
    const auto& g = g_idx_[static_cast<std::size_t>(ig)];
    std::size_t ix = static_cast<std::size_t>(g.x + gx_offset_);
    std::size_t iy = static_cast<std::size_t>(g.y + gy_offset_);
    std::size_t iz = static_cast<std::size_t>(g.z + gz_offset_);
    return kp * p_x_[ix] * p_y_[iy] * p_z_[iz];
}

auto StructureFactor::prepare(std::span<const vector3d<int>> gs) -> void {
    auto n = gs.size();
    g_idx_.assign(gs.begin(), gs.end());

    if (n == 0) {
        p_x_.clear(); p_y_.clear(); p_z_.clear();
        gx_offset_ = gy_offset_ = gz_offset_ = 0;
        return;
    }

    // G-component range.
    int gx_min = gs[0].x, gx_max = gs[0].x;
    int gy_min = gs[0].y, gy_max = gs[0].y;
    int gz_min = gs[0].z, gz_max = gs[0].z;
    for (std::size_t i = 1; i < n; ++i) {
        if      (gs[i].x < gx_min) gx_min = gs[i].x;
        else if (gs[i].x > gx_max) gx_max = gs[i].x;
        if      (gs[i].y < gy_min) gy_min = gs[i].y;
        else if (gs[i].y > gy_max) gy_max = gs[i].y;
        if      (gs[i].z < gz_min) gz_min = gs[i].z;
        else if (gs[i].z > gz_max) gz_max = gs[i].z;
    }

    gx_offset_ = -gx_min;
    gy_offset_ = -gy_min;
    gz_offset_ = -gz_min;

    int gx_size = gx_max - gx_min + 1;
    int gy_size = gy_max - gy_min + 1;
    int gz_size = gz_max - gz_min + 1;

    p_x_.resize(static_cast<std::size_t>(gx_size));
    p_y_.resize(static_cast<std::size_t>(gy_size));
    p_z_.resize(static_cast<std::size_t>(gz_size));

    for (int gx = gx_min; gx <= gx_max; ++gx)
        p_x_[static_cast<std::size_t>(gx + gx_offset_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gx * tau_.x));
    for (int gy = gy_min; gy <= gy_max; ++gy)
        p_y_[static_cast<std::size_t>(gy + gy_offset_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gy * tau_.y));
    for (int gz = gz_min; gz <= gz_max; ++gz)
        p_z_[static_cast<std::size_t>(gz + gz_offset_)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gz * tau_.z));
}

void StructureFactor::operator()(
    std::span<const vector3d<int>> gs,
    const vector3d<double>& k,
    std::span<std::complex<double>> out) const
{
    auto n = gs.size();

    if (n == 0) return;

    // k-phase (one complex exp for the whole batch).
    double kd = -2.0 * std::numbers::pi
              * (k.x * tau_.x + k.y * tau_.y + k.z * tau_.z);
    auto kp = std::exp(std::complex<double>(0, kd));

    // G-component range.
    int gx_min = gs[0].x, gx_max = gs[0].x;
    int gy_min = gs[0].y, gy_max = gs[0].y;
    int gz_min = gs[0].z, gz_max = gs[0].z;
    for (std::size_t i = 1; i < n; ++i) {
        if      (gs[i].x < gx_min) gx_min = gs[i].x;
        else if (gs[i].x > gx_max) gx_max = gs[i].x;
        if      (gs[i].y < gy_min) gy_min = gs[i].y;
        else if (gs[i].y > gy_max) gy_max = gs[i].y;
        if      (gs[i].z < gz_min) gz_min = gs[i].z;
        else if (gs[i].z > gz_max) gz_max = gs[i].z;
    }

    int gx_size = gx_max - gx_min + 1;
    int gy_size = gy_max - gy_min + 1;
    int gz_size = gz_max - gz_min + 1;
    int gx_off = -gx_min;
    int gy_off = -gy_min;
    int gz_off = -gz_min;

    // Build on-the-fly phase tables.
    std::vector<std::complex<double>> p_x(static_cast<std::size_t>(gx_size));
    std::vector<std::complex<double>> p_y(static_cast<std::size_t>(gy_size));
    std::vector<std::complex<double>> p_z(static_cast<std::size_t>(gz_size));

    for (int gx = gx_min; gx <= gx_max; ++gx)
        p_x[static_cast<std::size_t>(gx + gx_off)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gx * tau_.x));
    for (int gy = gy_min; gy <= gy_max; ++gy)
        p_y[static_cast<std::size_t>(gy + gy_off)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gy * tau_.y));
    for (int gz = gz_min; gz <= gz_max; ++gz)
        p_z[static_cast<std::size_t>(gz + gz_off)]
            = std::exp(std::complex<double>(0, -2.0 * std::numbers::pi * gz * tau_.z));

    // Fill output using three table lookups per entry.
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t ix = static_cast<std::size_t>(gs[i].x + gx_off);
        std::size_t iy = static_cast<std::size_t>(gs[i].y + gy_off);
        std::size_t iz = static_cast<std::size_t>(gs[i].z + gz_off);
        out[i] = kp * p_x[ix] * p_y[iy] * p_z[iz];
    }
}
