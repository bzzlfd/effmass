module;
#include "../physical_constants.hpp"

export module io.lattice;

import std;


export enum class LengthUnit {
    Bohr,
    Angstrom,
};


export class Lattice {
public:
    // Lattice() = default;  // Lattice is not yet mutable, so a default constructor is nonsensical.
    explicit Lattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit);
    explicit Lattice(std::span<const double, 9> flat, LengthUnit unit);
    Lattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit);

    auto A(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;
    auto B(LengthUnit unit = LengthUnit::Bohr) const -> std::array<std::array<double, 3>, 3>;

private:
    std::array<std::array<double, 3>, 3> A_{};
    std::array<std::array<double, 3>, 3> B_{};

    auto computeReciprocalLattice() -> void;

    static auto toArray(std::initializer_list<std::initializer_list<double>> A)
        -> std::array<std::array<double, 3>, 3>;
};


auto Lattice::toArray(std::initializer_list<std::initializer_list<double>> A)
    -> std::array<std::array<double, 3>, 3>
{
    if (A.size() != 3) {
        throw std::runtime_error("Lattice: expected 3 rows");
    }
    std::array<std::array<double, 3>, 3> result{};
    int n = 0;
    for (const auto& row : A) {
        if (row.size() != 3) {
            throw std::runtime_error("Lattice: expected 3 columns per row");
        }
        int c = 0;
        for (double val : row) {
            result[n][c] = val;
            ++c;
        }
        ++n;
    }
    return result;
}


Lattice::Lattice(std::span<const double, 9> flat, LengthUnit unit)
{
    for (int n = 0; n < 3; ++n) {
        for (int c = 0; c < 3; ++c) {
            A_[n][c] = flat[n * 3 + c];
        }
    }
    if (unit == LengthUnit::Angstrom) {
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                A_[n][c] /= BOHR_RADIUS_ANGSTROM;
            }
        }
    }
    computeReciprocalLattice();
}


Lattice::Lattice(std::initializer_list<std::initializer_list<double>> A, LengthUnit unit)
    : Lattice(toArray(A), unit)
{
}


auto Lattice::A(LengthUnit unit) const -> std::array<std::array<double, 3>, 3> {
    switch (unit) {
    case LengthUnit::Bohr:
        return A_;
    case LengthUnit::Angstrom: {
        auto result = A_;
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                result[n][c] *= BOHR_RADIUS_ANGSTROM;
            }
        }
        return result;
    }
    default:
        throw std::runtime_error("Lattice::A: unknown LengthUnit");
    }
}


auto Lattice::B(LengthUnit unit) const -> std::array<std::array<double, 3>, 3> {
    switch (unit) {
    case LengthUnit::Bohr:
        return B_;
    case LengthUnit::Angstrom: {
        auto result = B_;
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                result[n][c] /= BOHR_RADIUS_ANGSTROM;
            }
        }
        return result;
    }
    default:
        throw std::runtime_error("Lattice::B: unknown LengthUnit");
    }
}


Lattice::Lattice(const std::array<std::array<double, 3>, 3>& A, LengthUnit unit)
    : A_(A)
{
    if (unit == LengthUnit::Angstrom) {
        for (int n = 0; n < 3; ++n) {
            for (int c = 0; c < 3; ++c) {
                A_[n][c] /= BOHR_RADIUS_ANGSTROM;
            }
        }
    }
    computeReciprocalLattice();
}


auto Lattice::computeReciprocalLattice() -> void {
    // A_[n][c] = real-space lattice vector a_n (component c)
    // Reciprocal lattice: b1 = 2π * (a2 × a3) / V, b2 = 2π * (a3 × a1) / V, b3 = 2π * (a1 × a2) / V
    // where V = a1 · (a2 × a3)
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

    const double* a1 = A_[0].data();
    const double* a2 = A_[1].data();
    const double* a3 = A_[2].data();

    auto a2xa3 = cross(a2, a3);
    double V = dot(a1, a2xa3.data());
    if (std::abs(V) < 1e-12) {
        throw std::runtime_error("Lattice::computeReciprocalLattice: zero cell volume");
    }
    double invV = TWO_PI / V;

    auto a3xa1 = cross(a3, a1);
    auto a1xa2 = cross(a1, a2);

    for (int c = 0; c < 3; ++c) {
        B_[0][c] = a2xa3[c] * invV;
        B_[1][c] = a3xa1[c] * invV;
        B_[2][c] = a1xa2[c] * invV;
    }
}
