import io;
import std;


// --- helpers ---

static auto check(bool cond, std::string_view msg) -> void {
    if (!cond) throw std::runtime_error(std::string(msg));
}

static auto close(double a, double b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}


auto main() -> int {
    try {
        auto test = [](std::string_view label, const std::string& path) {
            std::println("=== {} ===", label);
            ATOM atom(path);

            check(atom.natom == 4, "natom != 4");
            std::println("  natom = {} [OK]", atom.natom);

            // Lattice vectors (Angstrom)
            auto A = atom.lattice.A(LengthUnit::Angstrom);

            check(close(A[0][0],  1.56429407), "a1x mismatch");
            check(close(A[0][1], -2.70943680), "a1y mismatch");
            check(close(A[0][2],  0.00000000), "a1z mismatch");

            check(close(A[1][0],  1.56429407), "a2x mismatch");
            check(close(A[1][1],  2.70943680), "a2y mismatch");
            check(close(A[1][2],  0.00000000), "a2z mismatch");

            check(close(A[2][0],  0.00000000), "a3x mismatch");
            check(close(A[2][1],  0.00000000), "a3y mismatch");
            check(close(A[2][2],  5.01695500), "a3z mismatch");
            std::println("  lattice vectors [OK]");

            // Atom species
            check(atom.species[0] == 13, "species[0] != 13");
            check(atom.species[1] == 13, "species[1] != 13");
            check(atom.species[2] ==  7, "species[2] != 7");
            check(atom.species[3] ==  7, "species[3] != 7");

            // Fractional coordinates
            check(close(atom.x[0], 0.66666667) &&
                  close(atom.y[0], 0.33333333) &&
                  close(atom.z[0], 0.49928700), "atom[0] position mismatch");
            check(close(atom.x[1], 0.33333333) &&
                  close(atom.y[1], 0.66666667) &&
                  close(atom.z[1], 0.99928700), "atom[1] position mismatch");
            check(close(atom.x[2], 0.66666667) &&
                  close(atom.y[2], 0.33333333) &&
                  close(atom.z[2], 0.88071300), "atom[2] position mismatch");
            check(close(atom.x[3], 0.33333333) &&
                  close(atom.y[3], 0.66666667) &&
                  close(atom.z[3], 0.38071300), "atom[3] position mismatch");
            std::println("  atoms [OK]");

            // Move semantics
            ATOM moved = std::move(atom);
            check(moved.natom == 4, "natom lost after move");
            check(atom.natom == 0, "source natom not zeroed after move");
            check(atom.species.empty(), "source species not empty after move");
            check(atom.x.empty(), "source x not empty after move");
            std::println("  move semantics [OK]");

            // print_info()
            moved.print_info();
            std::println("  print_info() [OK]");
        };

        test("test ATOM (nonlocal)", "test/data_io-nonlocal/atom.config");
        test("test ATOM (local)",    "test/data_io-local/atom.config");

        std::println("\nall tests passed");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "FAILED: {}", e.what());
        return 1;
    }
}
