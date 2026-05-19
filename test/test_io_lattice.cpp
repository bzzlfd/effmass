import io;
import std;


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


auto verifyB(const Lattice& lattice, const std::string& name) -> void
{
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    constexpr double tol = 1e-12;

    auto B_comp = lattice.B();
    auto B_ref = gaussJordanInverseTranspose(lattice.A());
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            B_ref[i][j] *= TWO_PI;
        }
    }

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


auto checkLattice(
    const std::array<std::array<double, 3>, 3>& A_input,
    LengthUnit unit,
    const std::string& name
) -> void
{
    Lattice lattice(A_input, unit);
    verifyB(lattice, name);
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
        verifyB(lattice_from_vec, "std::span constructor");

        // Verify Angstrom conversion explicitly (using initializer list)
        Lattice lattice_ang({{2.64588605335, 0.0, 0.0},
                             {0.0, 2.64588605335, 0.0},
                             {0.0, 0.0, 2.64588605335}},
                            LengthUnit::Angstrom);
        const auto A_ang = lattice_ang.A();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double expected = (i == j) ? 5.0 : 0.0;
                if (std::abs(A_ang[i][j] - expected) > 1e-12) {
                    throw std::runtime_error("Angstrom conversion mismatch");
                }
            }
        }

        // 5. Default constructor + setLattice (std::array)
        Lattice lattice_def;
        lattice_def.setLattice(cubic, LengthUnit::Bohr);
        const auto A_def = lattice_def.A();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (std::abs(A_def[i][j] - cubic[i][j]) > 1e-12) {
                    throw std::runtime_error("default ctor + setLattice(array): A mismatch");
                }
            }
        }
        verifyB(lattice_def, "default ctor + setLattice(array)");

        // 6. default ctor + setLattice (std::span)
        Lattice lattice_def_span;
        lattice_def_span.setLattice(std::span<const double, 9>(flat.data(), 9), LengthUnit::Bohr);
        verifyB(lattice_def_span, "default ctor + setLattice(span)");

        // 7. default ctor + setLattice (initializer_list)
        Lattice lattice_def_il;
        lattice_def_il.setLattice(
            {{{5.0, 0.0, 0.0}, {0.0, 5.0, 0.0}, {0.0, 0.0, 5.0}}},
            LengthUnit::Bohr
        );
        verifyB(lattice_def_il, "default ctor + setLattice(init list)");

        // 8. Modify existing lattice with setLattice
        Lattice lattice_mod(cubic, LengthUnit::Bohr);
        std::array<std::array<double, 3>, 3> fcc = {{
            {0.0, a, a},
            {a, 0.0, a},
            {a, a, 0.0}
        }};
        lattice_mod.setLattice(fcc, LengthUnit::Bohr);
        const auto A_mod = lattice_mod.A();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (std::abs(A_mod[i][j] - fcc[i][j]) > 1e-12) {
                    throw std::runtime_error("setLattice modify: A mismatch");
                }
            }
        }
        verifyB(lattice_mod, "setLattice modify");

        // 9. setLattice with Angstrom unit
        Lattice lattice_set_ang;
        lattice_set_ang.setLattice(
            {{{2.64588605335, 0.0, 0.0},
              {0.0, 2.64588605335, 0.0},
              {0.0, 0.0, 2.64588605335}}},
            LengthUnit::Angstrom
        );
        const auto A_set_ang = lattice_set_ang.A();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double expected = (i == j) ? 5.0 : 0.0;
                if (std::abs(A_set_ang[i][j] - expected) > 1e-12) {
                    throw std::runtime_error("setLattice Angstrom conversion: A mismatch");
                }
            }
        }
        verifyB(lattice_set_ang, "setLattice Angstrom");

        std::println("All Lattice tests passed!");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "Test failed: {}", e.what());
        return 1;
    }
}
