module;

export module io.ATOM;

import std;
import io.lattice;
import utils.vector3d;

// --- forward declarations for archived reference classes ---
namespace archived {
    class TypeView;
    class TypeAtomView;
}


export {
    class ATOM;
}


export class ATOM {
public:
    explicit ATOM(const std::string& filename);

    // Rule of zero: no FILE* or other resource handle

    auto print_info() const -> void;
    // element name ↔ atomic number conversion (Z = 1..112)
    static auto elementName(int atomic_number) -> std::string_view;
    static auto atomicNumber(std::string_view name) -> int;

    // --- member variables ---

    // parsed data (as-read order, unsorted)
    int natom{};
    Lattice lattice;
    // Note: separate vectors (species / x / y / z) rather than vector<Atom>.
    // Tradeoffs: flat arrays let callers pass individual buffers, and are
    // compact when only species are needed; sorting requires an external
    // permutation (typeSort_order_) instead of a trivial struct-sort.
    //
    // Positions are stored as vector3d<double> (fractional coordinates).
    std::vector<int> species;
    std::vector<vector3d<double>> frac_coord;

    // species analysis (computed during construction)
    int ntyp{};                     // number of distinct species types


    // --- iteration views ---
    struct TypeDatum { int z; int count; int ityp; };
    struct AtomDatum { int species; vector3d<double> coord; };

    // Iterate over each species type in ascending Z order.
    // Each step yields TypeDatum{z, count, ityp} for one type.
    // Usage: for (auto&& t : atom.eachType())
    class TypeView {
        friend class ATOM;
        const ATOM* a_;
        explicit TypeView(const ATOM* a) : a_(a) {}
      public:
        class Iterator {
            friend class TypeView;
            const ATOM* a_;
            int it_;
            Iterator(const ATOM* a, int it) : a_(a), it_(it) {}
          public:
            auto operator*() const -> TypeDatum { return {a_->typeGroup_z_[it_], a_->typeGroup_count_[it_], it_}; }
            auto operator++() -> Iterator& { ++it_; return *this; }
            bool operator!=(const Iterator& o) const { return it_ != o.it_; }
        };
        auto begin() const { return Iterator(a_, 0); }
        auto end() const { return Iterator(a_, a_->ntyp); }
    };

    // Iterate over species types in a caller-specified order.
    // `order` is a sequence of type indices (0..ntyp-1) in the desired
    // iteration order.  Each step yields TypeDatum{z, count, ityp}.
    // Usage: for (auto&& t : atom.eachType({2, 0, 1}))
    class TypeOrderView {
        friend class ATOM;
        const ATOM* a_;
        std::span<const int> order_;
        TypeOrderView(const ATOM* a, std::span<const int> order) : a_(a), order_(order) {}
      public:
        class Iterator {
            friend class TypeOrderView;
            const ATOM* a_;
            const int* order_;
            int pos_;
            Iterator(const ATOM* a, const int* order, int pos) : a_(a), order_(order), pos_(pos) {}
          public:
            auto operator*() const -> TypeDatum {
                int it = order_[pos_];
                return {a_->typeGroup_z_[it], a_->typeGroup_count_[it], it};
            }
            auto operator++() -> Iterator& { ++pos_; return *this; }
            bool operator!=(const Iterator& o) const { return pos_ != o.pos_; }
        };
        auto begin() const { return Iterator(a_, order_.data(), 0); }
        auto end() const { return Iterator(a_, order_.data(), static_cast<int>(order_.size())); }
    };

    // Iterate over all atoms of a given type, in type-grouped order
    // (i.e. the order produced by typeSort_order_).
    // Each step yields AtomDatum{species, x, y, z}.
    // Usage: for (auto&& a : atom.eachTypeAtom(ityp))
    class TypeAtomView {
        friend class ATOM;
        const ATOM* a_;
        int start_;
        int count_;
        TypeAtomView(const ATOM* a, int ityp) : a_(a) {
            start_ = 0;
            for (int t = 0; t < ityp; ++t) start_ += a->typeGroup_count_[t];
            count_ = a->typeGroup_count_[ityp];
        }
    public:
        class Iterator {
            friend class TypeAtomView;
            const ATOM* a_;
            int start_;
            int k_;
            Iterator(const ATOM* a, int start, int k) : a_(a), start_(start), k_(k) {}
        public:
            auto operator*() const -> AtomDatum {
                int orig = a_->typeSort_order_[start_ + k_];
                return {a_->species[orig], a_->frac_coord[orig]};
            }
            auto operator++() -> Iterator& { ++k_; return *this; }
            bool operator!=(const Iterator& o) const { return k_ != o.k_; }
        };
        auto begin() const { return Iterator(a_, start_, 0); }
        auto end() const { return Iterator(a_, start_, count_); }
    };

    // Iterate over all atoms in the original file-read order.
    // Each step yields AtomDatum{species, x, y, z}.
    // Usage: for (auto&& a : atom.eachSpecie())
    class SpecieView {
        friend class ATOM;
        const ATOM* a_;
        explicit SpecieView(const ATOM* a) : a_(a) {}
    public:
        class Iterator {
            friend class SpecieView;
            const ATOM* a_;
            int i_;
            Iterator(const ATOM* a, int i) : a_(a), i_(i) {}
        public:
            auto operator*() const -> AtomDatum {
                return {a_->species[i_], a_->frac_coord[i_]};
            }
            auto operator++() -> Iterator& { ++i_; return *this; }
            bool operator!=(const Iterator& o) const { return i_ != o.i_; }
        };
        auto begin() const { return Iterator(a_, 0); }
        auto end() const { return Iterator(a_, a_->natom); }
    };

    auto eachType() const -> TypeView { return TypeView(this); }
    auto eachType(std::span<const int> order) const -> TypeOrderView { return TypeOrderView(this, order); }
    auto eachTypeAtom(int ityp) const -> TypeAtomView { return TypeAtomView(this, ityp); }
    auto eachSpecie() const -> SpecieView { return SpecieView(this); }

private:
    friend class archived::TypeView;
    friend class archived::TypeAtomView;

    std::vector<int> typeGroup_z_;            // typeGroup_z_[it] = atomic number Z, sorted ascending, length ntyp
    std::vector<int> typeGroup_count_;        // typeGroup_count_[it] = how many atoms of that type, length ntyp
    std::vector<int> atom_type_idx_;          // atom_type_idx_[ia] = which type (0..ntyp-1) atom ia belongs to

    // typeSort_order_[k] = original index of k-th atom when grouped by type
    std::vector<int> typeSort_order_;

    auto analyzeSpecies() -> void;

};


// --- module-internal helpers (anonymous namespace) ---

namespace {

auto skipWhitespace(char*& p) -> void {
    while (*p == ' ' || *p == '\t') ++p;
}


auto trimTrailing(char* p) -> void {
    auto* end = p + std::strlen(p) - 1;
    while (end >= p && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '\r')) {
        *end = '\0';
        --end;
    }
}


auto parseAtomConfigFile(std::FILE* fp, int& natom, Lattice& lattice,
                         std::vector<int>& species,
                         std::vector<vector3d<double>>& coords) -> void
{
    char line[1024];

    // first line: natom
    if (!std::fgets(line, sizeof(line), fp)) {
        throw std::runtime_error("atom.config: unexpected EOF reading natom");
    }
    natom = static_cast<int>(std::strtol(line, nullptr, 10));

    // card loop: freely-ordered cards
    bool lattice_found = false;
    bool position_found = false;

    while (std::fgets(line, sizeof(line), fp)) {
        char* p = line;
        skipWhitespace(p);
        trimTrailing(p);
        if (*p == '\0') continue;  // skip blank lines

        if (std::strcmp(p, "LATTICE") == 0) {
            std::array<std::array<double, 3>, 3> A{};
            for (int i = 0; i < 3; ++i) {
                if (!std::fgets(line, sizeof(line), fp)) {
                    throw std::runtime_error(
                        "atom.config: unexpected EOF reading lattice vector "
                        + std::to_string(i)
                    );
                }
                char* ep;
                A[i][0] = std::strtod(line, &ep);
                A[i][1] = std::strtod(ep, &ep);
                A[i][2] = std::strtod(ep, nullptr);
            }
            lattice = Lattice(A, LengthUnit::Angstrom);
            lattice_found = true;
        }
        else if (std::strcmp(p, "POSITION") == 0) {
            species.resize(natom);
            coords.resize(natom);
            for (int i = 0; i < natom; ++i) {
                if (!std::fgets(line, sizeof(line), fp)) {
                    throw std::runtime_error(
                        "atom.config: unexpected EOF reading atom "
                        + std::to_string(i)
                    );
                }
                char* ep;
                species[i] = static_cast<int>(std::strtol(line, &ep, 10));
                double fx = std::strtod(ep, &ep);
                double fy = std::strtod(ep, &ep);
                double fz = std::strtod(ep, nullptr);
                coords[i] = {fx, fy, fz};
            }
            position_found = true;
        }
        // else: unknown card → skip (next fgets reads on)
    }

    if (!lattice_found) {
        throw std::runtime_error("atom.config: missing LATTICE card");
    }
    if (!position_found) {
        throw std::runtime_error("atom.config: missing POSITION card");
    }
}

} // anonymous namespace


// --- ATOM implementation ---

ATOM::ATOM(const std::string& filename) {
    auto* fp = std::fopen(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    try {
        parseAtomConfigFile(fp, natom, lattice, species, frac_coord);
    } catch (...) {
        std::fclose(fp);
        throw;
    }
    std::fclose(fp);
    analyzeSpecies();
}


auto ATOM::analyzeSpecies() -> void {
    // 1. sort unique species → type_z
    auto sorted_species = species;
    std::ranges::sort(sorted_species);
    typeGroup_z_.clear();
    for (int z : sorted_species) {
        if (typeGroup_z_.empty() || typeGroup_z_.back() != z) typeGroup_z_.push_back(z);
    }
    ntyp = static_cast<int>(typeGroup_z_.size());

    // 2. count atoms per type
    typeGroup_count_.assign(ntyp, 0);
    atom_type_idx_.resize(natom);
    for (int i = 0; i < natom; ++i) {
        auto it = std::ranges::lower_bound(typeGroup_z_, species[i]);
        int ityp = static_cast<int>(it - typeGroup_z_.begin());
        atom_type_idx_[i] = ityp;
        ++typeGroup_count_[ityp];
    }

    // 3. permutation that groups atoms by type (stable within type)
    typeSort_order_.resize(natom);
    std::iota(typeSort_order_.begin(), typeSort_order_.end(), 0);
    std::ranges::stable_sort(typeSort_order_, [&](int a, int b) {
        return atom_type_idx_[a] < atom_type_idx_[b];
    });
}


// --- element symbols, Z = 1..112 (matching gen_element_name_number.f90) ---

namespace {

constexpr std::string_view element_symbols[] = {
    "",   // index 0 unused; 1-based below
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn",
};

constexpr int max_element = static_cast<int>(std::size(element_symbols)) - 1;

} // anonymous namespace


auto ATOM::elementName(int atomic_number) -> std::string_view {
    if (atomic_number < 1 || atomic_number > max_element) {
        throw std::out_of_range(
            "atomic number not between 1 and " + std::to_string(max_element)
        );
    }
    return element_symbols[atomic_number];
}


auto ATOM::atomicNumber(std::string_view name) -> int {
    for (int z = 1; z <= max_element; ++z) {
        if (element_symbols[z] == name) return z;
    }
    throw std::invalid_argument(
        "unknown element symbol: " + std::string(name)
    );
}


auto ATOM::print_info() const -> void {
    std::println("ATOM: natom = {}, ntyp = {}", natom, ntyp);
    auto A_ang = lattice.A(LengthUnit::Angstrom);
    std::println("  Lattice vectors (Angstrom):");
    for (int i = 0; i < 3; ++i) {
        std::println("    {:14.8f}{:14.8f}{:14.8f}",
                     A_ang[i][0], A_ang[i][1], A_ang[i][2]);
    }
    for (int it = 0; it < ntyp; ++it) {
        std::println("  type {}: Z = {}, count = {}", it, typeGroup_z_[it], typeGroup_count_[it]);
    }
    int nshow = (std::cmp_less(natom, 5)) ? natom : 5;
    for (int i = 0; i < nshow; ++i) {
        std::println("  atom[{}]: species={}, type={},  frac=({:.8f}, {:.8f}, {:.8f})",
                     i, species[i], atom_type_idx_[i],
                     frac_coord[i].x, frac_coord[i].y, frac_coord[i].z);
    }
    if (std::cmp_greater(natom, 5)) {
        std::println("  ... and {} more atoms", natom - 5);
    }
}


// --- archived: alternative "view-as-iterator" iteration pattern (module-internal, for reference) ---
//
// Compared to the main ATOM iteration views (eachType / eachTypeAtom / eachSpecie):
//   - View class IS the iterator: all 5 range-for operations on one class
//   - No separate Iterator nested class
//   - Uses EndTag sentinel so end() returns a non-dereferenceable type
//   - Tradeoff: less type separation vs. simpler structure
//   - Tradeoff: `end()` sentinel has no data (one bool compare instead of field-vs-field)
//   - This approach is viable whenever the view only needs forward iteration

namespace archived {

struct EndTag {};

class TypeView {
    const ATOM* a_ = nullptr;
    int it_ = 0;
    friend ATOM;
    TypeView(const ATOM* a, int it) : a_(a), it_(it) {}
public:
    auto operator*() const -> ATOM::TypeDatum { return {a_->typeGroup_z_[it_], a_->typeGroup_count_[it_], it_}; }
    auto operator++() -> TypeView& { ++it_; return *this; }
    bool operator!=(EndTag) const { return it_ < a_->ntyp; }

    auto begin() const -> TypeView { return *this; }
    auto end() const -> EndTag { return {}; }
};


class TypeAtomView {
    const ATOM* a_ = nullptr;
    int start_ = 0;
    int k_ = 0;
    int count_ = 0;
    friend ATOM;
    TypeAtomView(const ATOM* a, int ityp) : a_(a) {
        start_ = 0;
        for (int t = 0; t < ityp; ++t) start_ += a->typeGroup_count_[t];
        count_ = a->typeGroup_count_[ityp];
    }
public:
    auto operator*() const -> ATOM::AtomDatum {
        int orig = a_->typeSort_order_[start_ + k_];
        return {a_->species[orig], a_->frac_coord[orig]};
    }
    auto operator++() -> TypeAtomView& { ++k_; return *this; }
    bool operator!=(EndTag) const { return k_ < count_; }

    auto begin() const -> TypeAtomView { return *this; }
    auto end() const -> EndTag { return {}; }
};


class SpecieView {
    const ATOM* a_ = nullptr;
    int i_ = 0;
    friend ATOM;
    SpecieView(const ATOM* a, int i) : a_(a), i_(i) {}
public:
    auto operator*() const -> ATOM::AtomDatum {
        return {a_->species[i_], a_->frac_coord[i_]};
    }
    auto operator++() -> SpecieView& { ++i_; return *this; }
    bool operator!=(EndTag) const { return i_ < a_->natom; }

    auto begin() const -> SpecieView { return *this; }
    auto end() const -> EndTag { return {}; }
};

} // namespace archived
