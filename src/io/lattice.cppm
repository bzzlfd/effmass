module;
#include "../physical_constants.hpp"

export module io.lattice;

import std;
import utils.array2d;


export {
    class Lattice;
        enum class LengthUnit : int;
}


enum class LengthUnit {
    Bohr,
    Angstrom,
};


class Lattice {
public:
    Lattice() = default;
    explicit Lattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit);
    explicit Lattice(std::span<const double, 9> flat, LengthUnit unit);
    Lattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit);

    auto setLattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit) -> void;
    auto setLattice(std::span<const double, 9> flat, LengthUnit unit) -> void;
    auto setLattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit) -> void;

    auto A(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;
    auto B(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;
    auto volume() const -> double;

private:
    array2d<double, 3, 3> A_{};

    auto setLatticeFromFlat(std::span<const double, 9> flat, LengthUnit unit) -> void;
    auto computeReciprocalLattice() const -> std::array<std::array<double, 3>, 3>;

    static auto flattenInitList(std::initializer_list<std::initializer_list<double>> A)
        -> std::array<double, 9>;
};


auto Lattice::flattenInitList(std::initializer_list<std::initializer_list<double>> A)
    -> std::array<double, 9>
{
    if (A.size() != 3) {
        throw std::runtime_error("Lattice: expected 3 rows");
    }
    std::array<double, 9> result{};
    int n = 0;
    for (const auto& row : A) {
        if (row.size() != 3) {
            throw std::runtime_error("Lattice: expected 3 columns per row");
        }
        int c = 0;
        for (double val : row) {
            result[n * 3 + c] = val;
            ++c;
        }
        ++n;
    }
    return result;
}


Lattice::Lattice(std::span<const double, 9> flat, LengthUnit unit)
{
    setLattice(flat, unit);
}


Lattice::Lattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit)
{
    setLattice(A, unit);
}


Lattice::Lattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit)
{
    setLattice(A, unit);
}


auto Lattice::setLattice(std::span<const double, 9> flat, LengthUnit unit) -> void
{
    setLatticeFromFlat(flat, unit);
}


auto Lattice::setLattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit) -> void
{
    std::array<double, 9> flat{};
    for (int n = 0; n < 3; ++n) {
        for (int c = 0; c < 3; ++c) {
            flat[n * 3 + c] = A[n][c];
        }
    }
    setLatticeFromFlat(flat, unit);
}


auto Lattice::setLattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit) -> void
{
    setLatticeFromFlat(flattenInitList(A), unit);
}


auto Lattice::setLatticeFromFlat(std::span<const double, 9> flat, LengthUnit unit) -> void
{
    for (int i = 0; i < 9; ++i) {
        A_.data[i] = flat[i];
    }
    if (unit == LengthUnit::Angstrom) {
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                A_[n, c] /= BOHR_RADIUS_ANGSTROM;
            }
        }
    }
    computeReciprocalLattice(); // eager validation
}


auto Lattice::A(LengthUnit unit) const -> std::array<std::array<double, 3>, 3>
{
    std::array<std::array<double, 3>, 3> result{};
    for (int n = 0; n < 3; ++n) {
        for (int c = 0; c < 3; ++c) {
            result[n][c] = A_[n, c];
        }
    }
    if (unit == LengthUnit::Angstrom) {
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                result[n][c] *= BOHR_RADIUS_ANGSTROM;
            }
        }
    }
    return result;
}


auto Lattice::B(LengthUnit unit) const -> std::array<std::array<double, 3>, 3>
{
    auto B_mat = computeReciprocalLattice();

    switch (unit) {
    case LengthUnit::Bohr:
        return B_mat;
    case LengthUnit::Angstrom: {
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                B_mat[n][c] /= BOHR_RADIUS_ANGSTROM;
            }
        }
        return B_mat;
    }
    default:
        throw std::runtime_error("Lattice::B: unknown LengthUnit");
    }
}


auto Lattice::volume() const -> double
{
    auto cross = [](const double a[3], const double b[3]) -> std::array<double, 3> {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };
    };

    auto dot = [](const double a[3], const double b[3]) -> double {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    const double* a1 = &A_[0, 0];
    const double* a2 = &A_[1, 0];
    const double* a3 = &A_[2, 0];

    auto a2xa3 = cross(a2, a3);
    double V = dot(a1, a2xa3.data());
    return std::abs(V);
}


auto Lattice::computeReciprocalLattice() const -> std::array<std::array<double, 3>, 3>
{
    constexpr double TWO_PI = 2.0 * std::numbers::pi;

    auto cross = [](const double a[3], const double b[3]) -> std::array<double, 3> {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };
    };

    auto dot = [](const double a[3], const double b[3]) -> double {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    const double* a1 = &A_[0, 0];
    const double* a2 = &A_[1, 0];
    const double* a3 = &A_[2, 0];

    auto a2xa3 = cross(a2, a3);
    double V = dot(a1, a2xa3.data());
    if (std::abs(V) < 1e-12) {
        throw std::runtime_error("Lattice::computeReciprocalLattice: zero cell volume");
    }
    double invV = TWO_PI / V;

    auto a3xa1 = cross(a3, a1);
    auto a1xa2 = cross(a1, a2);

    std::array<std::array<double, 3>, 3> B_mat{};
    for (int c = 0; c < 3; ++c) {
        B_mat[0][c] = a2xa3[c] * invV;
        B_mat[1][c] = a3xa1[c] * invV;
        B_mat[2][c] = a1xa2[c] * invV;
    }
    return B_mat;
}
