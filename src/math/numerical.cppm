export module math.numerical;

import std;

export {
    enum class SimpsonMeshType : int;
    auto simpson(std::span<const double> f, std::span<const double> rab,
                 SimpsonMeshType mesh_type) -> double;
    auto lagrangeCubic(double f0, double f1, double f2, double f3, double px) -> double;
    auto lagrangeCubicDerivative(double f0, double f1, double f2, double f3,
                                 double px, double dq) -> double;
}


// =============================================================================
//  Enums
// =============================================================================

enum class SimpsonMeshType : int {
    Uniform,
    General,
};


// =============================================================================
//  Simpson integration
// =============================================================================

// Simpson integration on the index i (dr/di encoded in rab).
// For a Uniform mesh, rab(i) is constant (dr) and the formula is simplified.
// Handles even number of points without dropping the last point by combining
// Simpson 3/8 rule for the first 4 points and Simpson 1/3 for the remainder.
auto simpson(std::span<const double> f, std::span<const double> rab,
             SimpsonMeshType mesh_type) -> double {
    int n = static_cast<int>(f.size());
    int nrab = static_cast<int>(rab.size());
    if (n == 0) return 0.0;
    if (n == 1) return f[0] * rab[0];
    if (nrab < n) {
        throw std::invalid_argument("simpson: rab size must be >= f size");
    }

    if (mesh_type == SimpsonMeshType::Uniform) {
        double dr = rab[0];
        if (n % 2 == 1) {
            double sum = f[0] + f[n - 1];
            for (int i = 1; i < n - 1; ++i) {
                sum += (i % 2 == 0 ? 2.0 : 4.0) * f[i];
            }
            return dr * sum / 3.0;
        } else {
            if (n == 2) {
                return dr * 0.5 * (f[0] + f[1]);
            }
            double sum38 = f[0] + 3.0 * f[1] + 3.0 * f[2] + f[3];
            double sum13 = f[3] + f[n - 1];
            for (int i = 4; i < n - 1; ++i) {
                sum13 += (i % 2 == 0 ? 4.0 : 2.0) * f[i];
            }
            return dr * (3.0 / 8.0 * sum38 + sum13 / 3.0);
        }
    } else {
        // General mesh: g(i) = f(i) * rab(i), integrate over uniform index i.
        if (n % 2 == 1) {
            double sum = f[0] * rab[0] + f[n - 1] * rab[n - 1];
            for (int i = 1; i < n - 1; ++i) {
                sum += (i % 2 == 0 ? 2.0 : 4.0) * f[i] * rab[i];
            }
            return sum / 3.0;
        } else {
            if (n == 2) {
                return 0.5 * (f[0] * rab[0] + f[1] * rab[1]);
            }
            double g0 = f[0] * rab[0];
            double g1 = f[1] * rab[1];
            double g2 = f[2] * rab[2];
            double g3 = f[3] * rab[3];
            double sum38 = g0 + 3.0 * g1 + 3.0 * g2 + g3;
            double sum13 = g3 + f[n - 1] * rab[n - 1];
            for (int i = 4; i < n - 1; ++i) {
                sum13 += (i % 2 == 0 ? 4.0 : 2.0) * f[i] * rab[i];
            }
            return 3.0 / 8.0 * sum38 + sum13 / 3.0;
        }
    }
}


// =============================================================================
//  Lagrange 4-point cubic interpolation (same as QE interp_beta / interp_dbeta)
// =============================================================================

auto lagrangeCubic(
    double f0, double f1, double f2, double f3,
    double px
) -> double {
    double ux = 1.0 - px;
    double vx = 2.0 - px;
    double wx = 3.0 - px;
    return f0 * ux * vx * wx / 6.0 +
           f1 * px * vx * wx / 2.0 -
           f2 * px * ux * wx / 2.0 +
           f3 * px * ux * vx / 6.0;
}

auto lagrangeCubicDerivative(
    double f0, double f1, double f2, double f3,
    double px, double dq
) -> double {
    double ux = 1.0 - px;
    double vx = 2.0 - px;
    double wx = 3.0 - px;
    return (f0 * (-vx*wx - ux*wx - ux*vx) / 6.0 +
            f1 * ( vx*wx - px*wx - px*vx) / 2.0 -
            f2 * ( ux*wx - px*wx - px*ux) / 2.0 +
            f3 * ( ux*vx - px*vx - px*ux) / 6.0) / dq;
}
