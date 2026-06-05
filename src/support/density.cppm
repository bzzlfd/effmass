export module support.density;

import std;
import math;
import io;


export {

/// Result structure for compareDensity().
struct RHOCompareResult {
    double max_diff;
    double rmse;
    double avg_ref;
};

/// Build charge density from occupied wavefunctions:
///   ρ(r) = Σ_{k,b} occ_{k,b} · |ψ_{k,b}(r)|²
///
/// Returns a flat array indexed as [i][j][k] = i*n2*n3 + j*n3 + k,
/// directly comparable to RHO file data (same units).
///
/// Only handles non-spin-polarized, non-SO cases.
/// GKK and WG are passed as non-const because their internal cache/view
/// state is modified during the multi-pass reconstruction.
auto buildDensity(GKK& gkk, WG& wg, const OCC& occ) -> std::vector<double>;

/// Volume integral of charge density:  ∫ ρ(r) d³r
auto integrateDensity(const RHO& rho) -> double;

/// Compare reconstructed density against file RHO (point-by-point).
auto compareDensity(const std::vector<double>& reconstructed, const RHO& rho_ref) -> RHOCompareResult;

} // export


// ===========================================================================
//  Implementation
// ===========================================================================

auto buildDensity(GKK& gkk, WG& wg, const OCC& occ) -> std::vector<double> {
    if (gkk.meta.islda != 1 || gkk.meta.is_SO != 0) {
        throw std::runtime_error(
            "buildDensity: only non-spin-polarized, non-SO data is supported");
    }

    int n1  = gkk.meta.n1;
    int n2  = gkk.meta.n2;
    int n3  = gkk.meta.n3;
    auto n123 = static_cast<std::size_t>(n1) * n2 * n3;
    int nband = wg.meta.nband;
    int nkpt  = gkk.meta.nkpt;

    std::vector<double>               rho_acc(n123, 0.0);
    std::vector<std::complex<double>> buf(n123);
    FFT3D fft(n1, n2, n3);

    for (int ikpt = 0; ikpt < nkpt; ++ikpt) {
        gkk.setDataView(KVecsView::Cartesian | KVecsView::Integer);
        auto& kv = gkk.loadKPoint(ikpt);
        int ng = static_cast<int>(kv.g_idx.size());

        for (int iband = 0; iband < nband; ++iband) {
            double occ_val = occ.occupation(iband, ikpt);
            if (occ_val == 0.0) continue;

            auto coeffs = wg.loadBand(ikpt, iband);

            if (static_cast<int>(coeffs.up.size()) != ng) {
                throw std::runtime_error(
                    "buildDensity: WG coefficient count (" +
                    std::to_string(coeffs.up.size()) +
                    ") != GKK G-vector count (" + std::to_string(ng) +
                    ") at kpt=" + std::to_string(ikpt) +
                    " band=" + std::to_string(iband));
            }

            std::fill(buf.begin(), buf.end(), 0.0);
            for (int ig = 0; ig < ng; ++ig) {
                int i_idx = ((kv.g_idx[ig].x % n1) + n1) % n1;
                int j_idx = ((kv.g_idx[ig].y % n2) + n2) % n2;
                int k_idx = ((kv.g_idx[ig].z % n3) + n3) % n3;
                buf[static_cast<std::size_t>(i_idx) * n2 * n3
                  + static_cast<std::size_t>(j_idx) * n3
                  + static_cast<std::size_t>(k_idx)] = coeffs.up[ig];
            }
            fft(buf, G2R);
            for (std::size_t ir = 0; ir < n123; ++ir)
                rho_acc[ir] += occ_val * std::norm(buf[ir]);
        }
    }

    return rho_acc;
}


auto integrateDensity(const RHO& rho) -> double {
    int n1  = rho.meta.n1;
    int n2  = rho.meta.n2;
    int n3  = rho.meta.n3;
    auto n123 = static_cast<std::size_t>(n1) * n2 * n3;
    double vol = rho.lattice.volume();

    double sum = 0.0;
    for (int i = 0; i < n1; ++i)
        for (int j = 0; j < n2; ++j)
            for (int k = 0; k < n3; ++k)
                sum += rho[i, j, k];

    return sum * vol / static_cast<double>(n123);
}


auto compareDensity(const std::vector<double>& reconstructed, const RHO& rho_ref) -> RHOCompareResult {
    int n1  = rho_ref.meta.n1;
    int n2  = rho_ref.meta.n2;
    int n3  = rho_ref.meta.n3;
    auto n123 = static_cast<std::size_t>(n1) * n2 * n3;

    double max_d = 0.0, sum_sq = 0.0, sum_abs_ref = 0.0;
    for (int i = 0; i < n1; ++i)
        for (int j = 0; j < n2; ++j)
            for (int k = 0; k < n3; ++k) {
                auto idx = static_cast<std::size_t>(i) * n2 * n3
                         + static_cast<std::size_t>(j) * n3 + k;
                double d = reconstructed[idx] - rho_ref[i, j, k];
                sum_sq += d * d;
                sum_abs_ref += std::abs(rho_ref[i, j, k]);
                double ad = std::abs(d);
                if (ad > max_d) max_d = ad;
            }

    return {
        .max_diff = max_d,
        .rmse     = std::sqrt(sum_sq / static_cast<double>(n123)),
        .avg_ref  = sum_abs_ref / static_cast<double>(n123),
    };
}
