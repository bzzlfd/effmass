import io;
import std;


// --- helpers ---

static auto check(bool cond, std::string_view msg) -> void {
    if (!cond) throw std::runtime_error(std::string(msg));
}

static auto close(double a, double b, double eps = 1e-12) -> bool {
    return std::abs(a - b) < eps;
}


// --- element conversion tests ---

static auto testElementConversion() -> void {
    std::println("=== element name/number conversion ===");

    // Z → name
    check(ATOM::elementName(1)  == "H",   "Z=1 should be H");
    check(ATOM::elementName(6)  == "C",   "Z=6 should be C");
    check(ATOM::elementName(7)  == "N",   "Z=7 should be N");
    check(ATOM::elementName(13) == "Al",  "Z=13 should be Al");
    check(ATOM::elementName(79) == "Au",  "Z=79 should be Au");
    check(ATOM::elementName(112) == "Cn", "Z=112 should be Cn");
    std::println("  elementName [OK]");

    // name → Z
    check(ATOM::atomicNumber("H")  == 1,  "\"H\" → Z=1");
    check(ATOM::atomicNumber("C")  == 6,  "\"C\" → Z=6");
    check(ATOM::atomicNumber("N")  == 7,  "\"N\" → Z=7");
    check(ATOM::atomicNumber("Al") == 13, "\"Al\" → Z=13");
    check(ATOM::atomicNumber("Au") == 79, "\"Au\" → Z=79");
    check(ATOM::atomicNumber("Cn") == 112,"\"Cn\" → Z=112");
    std::println("  atomicNumber [OK]");

    // round-trip
    for (int z = 1; z <= 112; ++z) {
        auto name = ATOM::elementName(z);
        int z2 = ATOM::atomicNumber(name);
        check(z == z2, std::format("round-trip failed for Z={}, name={}", z, name));
    }
    std::println("  round-trip Z→name→Z for 1..112 [OK]");

    // out-of-range
    bool caught_oor = false;
    try { ATOM::elementName(0); } catch (const std::out_of_range&) { caught_oor = true; }
    check(caught_oor, "elementName(0) should throw out_of_range");

    caught_oor = false;
    try { ATOM::elementName(113); } catch (const std::out_of_range&) { caught_oor = true; }
    check(caught_oor, "elementName(113) should throw out_of_range");

    std::println("  out-of-range [OK]");

    // unknown symbol
    bool caught_inv = false;
    try { ATOM::atomicNumber("Zz"); } catch (const std::invalid_argument&) { caught_inv = true; }
    check(caught_inv, "atomicNumber(\"Zz\") should throw invalid_argument");
    std::println("  unknown symbol [OK]");
}


auto main() -> int {
    try {
        testElementConversion();

        auto test = [](std::string_view label, const std::string& path) {
            std::println("=== {} ===", label);
            ATOM atom(path);

            // --- parsing ---
            check(atom.natom == 4, "natom != 4");
            std::println("  natom = {} [OK]", atom.natom);

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

            // --- species analysis ---
            check(atom.ntyp == 2, "ntyp != 2");
            check(atom.zval[0] ==  7, "zval[0] != 7");
            check(atom.zval[1] == 13, "zval[1] != 13");
            check(atom.type_count[0] == 2, "type_count[0] != 2");
            check(atom.type_count[1] == 2, "type_count[1] != 2");
            std::println("  ntyp = {}, zval = [{}, {}], counts = [{}, {}] [OK]",
                         atom.ntyp, atom.zval[0], atom.zval[1],
                         atom.type_count[0], atom.type_count[1]);

            // atom_type: original order: [13, 13, 7, 7] → type [1, 1, 0, 0]
            check(atom.atom_type[0] == 1, "atom_type[0] != 1");
            check(atom.atom_type[1] == 1, "atom_type[1] != 1");
            check(atom.atom_type[2] == 0, "atom_type[2] != 0");
            check(atom.atom_type[3] == 0, "atom_type[3] != 0");

            // sorted_idx: groups by type → [type0=7, type0=7, type1=13, type1=13]
            //             → original indices: [2, 3, 0, 1]
            check(atom.sorted_idx[0] == 2, "sorted_idx[0] != 2");
            check(atom.sorted_idx[1] == 3, "sorted_idx[1] != 3");
            check(atom.sorted_idx[2] == 0, "sorted_idx[2] != 0");
            check(atom.sorted_idx[3] == 1, "sorted_idx[3] != 1");
            std::println("  sorted_idx = [{}, {}, {}, {}] [OK]",
                         atom.sorted_idx[0], atom.sorted_idx[1],
                         atom.sorted_idx[2], atom.sorted_idx[3]);

            // sorted_idx is a valid permutation: each index 0..natom-1 appears exactly once
            {
                std::vector<bool> seen(atom.natom, false);
                for (int k = 0; k < atom.natom; ++k) {
                    int idx = atom.sorted_idx[k];
                    check(idx >= 0 && idx < atom.natom, "sorted_idx: out of range");
                    check(!seen[idx], "sorted_idx: duplicate");
                    seen[idx] = true;
                }
                for (int i = 0; i < atom.natom; ++i) {
                    check(seen[i], std::format("sorted_idx: index {} missing", i));
                }
            }

            // sorted_idx groups atoms by type
            for (int k = 1; k < atom.natom; ++k) {
                int ta = atom.atom_type[atom.sorted_idx[k - 1]];
                int tb = atom.atom_type[atom.sorted_idx[k]];
                check(ta <= tb, std::format("sorted_idx: type out of order at k={}", k));
            }
            std::println("  sorted_idx groups by type [OK]");

            // --- move semantics ---
            ATOM moved = std::move(atom);
            check(moved.natom == 4, "natom lost after move");
            check(moved.ntyp  == 2, "ntyp lost after move");
            check(moved.zval.size() == 2, "zval lost after move");
            check(atom.natom == 0, "source natom not zeroed after move");
            check(atom.ntyp  == 0, "source ntyp not zeroed after move");
            check(atom.zval.empty(), "source zval not empty after move");
            std::println("  move semantics [OK]");

            // --- print_info ---
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
