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

            // type info via eachType()
            {
                int seen = 0;
                int expected_z[] = {7, 13};
                int expected_count[] = {2, 2};
                for (auto&& t : atom.eachType()) {
                    check(t.z == expected_z[seen], "eachType: Z mismatch");
                    check(t.count == expected_count[seen], "eachType: count mismatch");
                    ++seen;
                }
                check(seen == 2, "eachType: expected 2 types");
                std::println("  ntyp = {}, eachType = [Z={}(x{}), Z={}(x{})] [OK]",
                             atom.ntyp,
                             expected_z[0], expected_count[0],
                             expected_z[1], expected_count[1]);
            }

            // iteration views below verify the grouping indirectly.

            // --- iteration views ---
            {
                // eachType: should yield [type0: Z=7 count=2, type1: Z=13 count=2] in order
                int seen_types = 0;
                for (auto&& t : atom.eachType()) {
                    check(seen_types < 2, "eachType: too many iterations");
                    if (seen_types == 0) {
                        check(t.z == 7, "eachType[0]: z != 7");
                        check(t.count == 2, "eachType[0]: count != 2");
                        check(t.ityp == 0, "eachType[0]: ityp != 0");
                    } else {
                        check(t.z == 13, "eachType[1]: z != 13");
                        check(t.count == 2, "eachType[1]: count != 2");
                        check(t.ityp == 1, "eachType[1]: ityp != 1");
                    }
                    ++seen_types;
                }
                check(seen_types == 2, "eachType: expected 2 types");
                std::println("  eachType [OK]");

                // eachType with custom order: Al then N
                {
                    int custom_order[] = {1, 0};
                    int seen = 0;
                    int expected_z[] = {13, 7};
                    int expected_ityp[] = {1, 0};
                    for (auto&& t : atom.eachType(custom_order)) {
                        check(t.z == expected_z[seen], "eachType(order): Z mismatch");
                        check(t.ityp == expected_ityp[seen], "eachType(order): ityp mismatch");
                        ++seen;
                    }
                    check(seen == 2, "eachType(order): expected 2 types");
                    std::println("  eachType(custom order) [OK]");
                }
            }

            {
                // eachTypeAtom(0) = type N (Z=7): 2 atoms
                int count0 = 0;
                for (auto&& a : atom.eachTypeAtom(0)) {
                    check(a.species == 7, std::format("eachTypeAtom(0)[{}]: species != 7", count0));
                    ++count0;
                }
                check(count0 == 2, "eachTypeAtom(0): expected 2 atoms");

                // eachTypeAtom(1) = type Al (Z=13): 2 atoms
                int count1 = 0;
                for (auto&& a : atom.eachTypeAtom(1)) {
                    check(a.species == 13, std::format("eachTypeAtom(1)[{}]: species != 13", count1));
                    ++count1;
                }
                check(count1 == 2, "eachTypeAtom(1): expected 2 atoms");
                std::println("  eachTypeAtom [OK]");
            }

            {
                // eachSpecie: original order [13, 13, 7, 7]
                int idx = 0;
                int expected[] = {13, 13, 7, 7};
                for (auto&& a : atom.eachSpecie()) {
                    check(a.species == expected[idx],
                          std::format("eachSpecie[{}]: species {} != expected {}", idx, a.species, expected[idx]));
                    ++idx;
                }
                check(idx == 4, "eachSpecie: expected 4 atoms");
                std::println("  eachSpecie [OK]");
            }

            // --- move semantics ---
            ATOM moved = std::move(atom);
            check(moved.natom == 4, "natom lost after move");
            check(moved.ntyp  == 2, "ntyp lost after move");
            std::println("  move semantics [OK]");

            // --- print_info ---
            moved.print_info();
            std::println("  print_info() [OK]");
        };

        test("test ATOM (nonlocal)", "test/data_scf/atom.config");

        std::println("\nall tests passed");
        return 0;
    } catch (const std::exception& e) {
        std::println(std::cerr, "FAILED: {}", e.what());
        return 1;
    }
}
