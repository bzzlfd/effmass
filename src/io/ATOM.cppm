module;

export module io.ATOM;

import std;
import io.lattice;


export {
    class ATOM;
}


export class ATOM {
public:
    explicit ATOM(const std::string& filename);
    ~ATOM();

    ATOM(const ATOM&) = delete;
    auto operator=(const ATOM&) -> ATOM& = delete;

    ATOM(ATOM&&) noexcept;
    auto operator=(ATOM&&) noexcept -> ATOM&;

    auto print_info() const -> void;

    int natom{};
    Lattice lattice;
    std::vector<int> species;
    std::vector<double> x, y, z;   // fractional coordinates

private:
    std::FILE* fp_ = nullptr;
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
                         std::vector<double>& x,
                         std::vector<double>& y,
                         std::vector<double>& z) -> void
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
            x.resize(natom);
            y.resize(natom);
            z.resize(natom);
            for (int i = 0; i < natom; ++i) {
                if (!std::fgets(line, sizeof(line), fp)) {
                    throw std::runtime_error(
                        "atom.config: unexpected EOF reading atom "
                        + std::to_string(i)
                    );
                }
                char* ep;
                species[i] = static_cast<int>(std::strtol(line, &ep, 10));
                x[i] = std::strtod(ep, &ep);
                y[i] = std::strtod(ep, &ep);
                z[i] = std::strtod(ep, nullptr);
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
    fp_ = std::fopen(filename.c_str(), "r");
    if (!fp_) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    parseAtomConfigFile(fp_, natom, lattice, species, x, y, z);
}


ATOM::~ATOM() {
    if (fp_) std::fclose(fp_);
}


ATOM::ATOM(ATOM&& other) noexcept
    : natom(std::exchange(other.natom, 0))
    , lattice(std::move(other.lattice))
    , species(std::move(other.species))
    , x(std::move(other.x))
    , y(std::move(other.y))
    , z(std::move(other.z))
    , fp_(std::exchange(other.fp_, nullptr))
{}


auto ATOM::operator=(ATOM&& other) noexcept -> ATOM& {
    if (this != &other) {
        if (fp_) std::fclose(fp_);
        natom    = std::exchange(other.natom, 0);
        lattice  = std::move(other.lattice);
        species  = std::move(other.species);
        x        = std::move(other.x);
        y        = std::move(other.y);
        z        = std::move(other.z);
        fp_      = std::exchange(other.fp_, nullptr);
    }
    return *this;
}


auto ATOM::print_info() const -> void {
    std::println("ATOM: natom = {}", natom);
    auto A_ang = lattice.A(LengthUnit::Angstrom);
    std::println("  Lattice vectors (Angstrom):");
    for (int i = 0; i < 3; ++i) {
        std::println("    {:14.8f}{:14.8f}{:14.8f}",
                     A_ang[i][0], A_ang[i][1], A_ang[i][2]);
    }
    int nshow = (std::cmp_less(natom, 5)) ? natom : 5;
    for (int i = 0; i < nshow; ++i) {
        std::println("  atom[{}]: species={},  frac=({:.8f}, {:.8f}, {:.8f})",
                     i, species[i], x[i], y[i], z[i]);
    }
    if (std::cmp_greater(natom, 5)) {
        std::println("  ... and {} more atoms", natom - 5);
    }
}


// --- archived: alternative eager-close implementation (module-internal, for reference) ---

namespace archived {

struct SimpleATOM {
    int natom{};
    Lattice lattice;
    std::vector<int> species;
    std::vector<double> x, y, z;

    explicit SimpleATOM(const std::string& filename) {
        auto* fp = std::fopen(filename.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        try {
            parseAtomConfigFile(fp, natom, lattice, species, x, y, z);
        } catch (...) {
            std::fclose(fp);
            throw;
        }
        std::fclose(fp);
    }
};

} // namespace archived
