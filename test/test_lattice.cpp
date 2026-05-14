import io;
import std;


// Gauss-Jordan elimination on augmented matrix [A | I] to compute A^{-1},
// then return its transpose: (A^{-1})^T.
auto gaussJordanInverseTranspose(const std::array<std::array<double, 3>, 3>& A)
    -> std::array<std::array<double, 3>, 3>
{
    double aug[3][6];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            aug[i][j] = A[i][j];
        }
        for (int j = 0; j < 3; ++j) {
            aug[i][3 + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        double maxVal = std::abs(aug[col][col]);
        for (int row = col + 1; row < 3; ++row) {
            if (std::abs(aug[row][col]) > maxVal) {
                maxVal = std::abs(aug[row][col]);
                pivot = row;
            }
        }
        if (maxVal < 1e-15) {
            throw std::runtime_error("gaussJordanInverseTranspose: singular matrix");
        }
        if (pivot != col) {
            for (int j = 0; j < 6; ++j) {
                std::swap(aug[col][j], aug[pivot][j]);
            }
        }

        double pivotVal = aug[col][col];
        for (int j = 0; j < 6; ++j) {
            aug[col][j] /= pivotVal;
        }

        for (int row = 0; row < 3; ++row) {
            if (row == col) continue;
            double factor = aug[row][col];
            if (std::abs(factor) > 0.0) {
                for (int j = 0; j < 6; ++j) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }
    }

    std::array<std::array<double, 3>, 3> A_inv{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A_inv[i][j] = aug[i][3 + j];
        }
    }

    std::array<std::array<double, 3>, 3> result{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = A_inv[j][i];
        }
    }
    return result;
}


auto checkLattice(
    const std::array<std::array<double, 3>, 3>& A_input,
    LengthUnit unit,
    const std::string& name
) -> void
{
    Lattice lattice(A_input, unit);

    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    constexpr double tol = 1e-12;

    auto B_ref = gaussJordanInverseTranspose(lattice.A());
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            B_ref[i][j] *= TWO_PI;
        }
    }

    const auto& B_comp = lattice.B();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double diff = std::abs(B_comp[i][j] - B_ref[i][j]);
            if (diff > tol) {
                throw std::runtime_error(
                    name + ": B[" + std::to_string(i) + "][" + std::to_string(j) +
                    "] mismatch: expected " + std::to_string(B_ref[i][j]) +
                    ", got " + std::to_string(B_comp[i][j]) +
                    ", diff = " + std::to_string(diff)
                );
            }
        }
    }
}


auto main() -> int {
    try {
        // 1. Simple cubic lattice (using std::array)
        std::array<std::array<double, 3>, 3> cubic = {{
            {5.0, 0.0, 0.0},
            {0.0, 5.0, 0.0},
            {0.0, 0.0, 5.0}
        }};
        checkLattice(cubic, LengthUnit::Bohr, "simple cubic (Bohr, array)");

        // 1b. Simple cubic lattice (using initializer list)
        checkLattice(
            {{{5.0, 0.0, 0.0}, {0.0, 5.0, 0.0}, {0.0, 0.0, 5.0}}},
            LengthUnit::Bohr,
            "simple cubic (Bohr, init list)"
        );

        // 2. FCC lattice (non-orthogonal, using initializer list)
        double a = 2.7;
        checkLattice(
            {{{0.0, a, a}, {a, 0.0, a}, {a, a, 0.0}}},
            LengthUnit::Bohr,
            "FCC (Bohr)"
        );

        // 3. Simple cubic in Angstrom (using initializer list)
        // 2.645886... Angstrom = 5.0 Bohr (since 1 Bohr = 0.52917721067 Å)
        checkLattice(
            {{{2.64588605335, 0.0, 0.0},
              {0.0, 2.64588605335, 0.0},
              {0.0, 0.0, 2.64588605335}}},
            LengthUnit::Angstrom,
            "simple cubic (Angstrom)"
        );

        // 4. Simple cubic from std::vector<double> (using std::span)
        std::vector<double> flat = {5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0};
        Lattice lattice_from_vec(std::span<const double, 9>(flat.data(), 9), LengthUnit::Bohr);
        const auto& B_vec = lattice_from_vec.B();
        const auto& B_ref = lattice_from_vec.B(); // same object, just re-verify
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (std::abs(B_vec[i][j] - B_ref[i][j]) > 1e-12) {
                    throw std::runtime_error("std::vector constructor mismatch");
                }
            }
        }

        // Verify Angstrom conversion explicitly (using initializer list)
        Lattice lattice_ang({{2.64588605335, 0.0, 0.0},
                             {0.0, 2.64588605335, 0.0},
                             {0.0, 0.0, 2.64588605335}},
                            LengthUnit::Angstrom);
        const auto& A_ang = lattice_ang.A();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double expected = (i == j) ? 5.0 : 0.0;
                if (std::abs(A_ang[i][j] - expected) > 1e-12) {
                    throw std::runtime_error("Angstrom conversion mismatch");
                }
            }
        }

        std::println("All Lattice tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
