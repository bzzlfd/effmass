module;
#include "pugixml.hpp"

export module pseudo.io.ncpp_upf;

import std;


export {
    class NCPPUPF;
        struct NCPPUPFHeader;
        struct NCPPUPFMesh;
            enum class MeshType : int;
        struct NCPPUPFNonlocal;
            struct NCPPUPFNonlocalByL;
        struct NCPPUPFWavefunction;
}


// === Simple dense matrix (flat storage) ===
struct Matrix {
    std::vector<double> data;
    int rows = 0;
    int cols = 0;

    auto operator[](int i, int j) const -> const double& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }
    auto operator[](int i, int j) -> double& {
        return data[static_cast<std::size_t>(i) * cols + j];
    }

    [[nodiscard]] auto size() const -> std::size_t { return data.size(); }
    [[nodiscard]] auto empty() const -> bool { return data.empty(); }
};


// === NCPPUPF Header structure ===
struct NCPPUPFHeader {
    std::string generated;
    std::string author;
    std::string date;
    std::string comment;
    std::string element;
    std::string pseudo_type;
    std::string relativistic;
    bool is_ultrasoft = false;
    bool is_paw = false;
    bool is_coulomb = false;
    bool has_so = false;
    bool has_wfc = false;
    bool has_gipaw = false;
    bool core_correction = false;
    std::string functional;
    double z_valence = 0.0;
    double total_psenergy = 0.0;
    double wfc_cutoff = 0.0;
    double rho_cutoff = 0.0;
    int l_max = 0;
    int l_local = 0;
    int mesh_size = 0;
    int number_of_wfc = 0;
    int number_of_proj = 0;
};


// === Mesh data ===
struct NCPPUPFMesh {
    std::vector<double> r;
    std::vector<double> rab;
};


// === Nonlocal data (beta projectors + D_ij) ===
struct NCPPUPFNonlocal {
    std::vector<std::vector<double>> beta;  // [nbeta][mesh]
    std::vector<int> lll;                    // [nbeta] angular momentum
    std::vector<int> kbeta;                  // [nbeta] cutoff radius index
    std::vector<double> rcut;                // [nbeta] cutoff radius
    Matrix dion;                             // D_ij matrix
};


// === Mesh type enumeration ===
enum class MeshType {
    Uniform,
    Exponential,
    Unknown
};


// === Nonlocal data filtered by angular momentum l ===
struct NCPPUPFNonlocalByL {
    std::vector<std::vector<double>> beta;  // [n_beta_l][mesh]
    Matrix dion;                             // D_ij submatrix for this l
};


// === Wavefunction data (pseudo atomic orbitals) ===
struct NCPPUPFWavefunction {
    std::vector<std::vector<double>> chi;   // [nwfc][effective_mesh]
    std::vector<int> kchi;                   // [nwfc] effective length (truncated trailing zeros)
    std::vector<int> lchi;                   // [nwfc] angular momentum
    std::vector<double> oc;                  // [nwfc] occupation
    std::vector<std::string> labels;         // [nwfc] label
};


// === NCPPUPF class ===
class NCPPUPF {
public:
    explicit NCPPUPF(const std::string& filename);

    auto header() const -> const NCPPUPFHeader& { return header_; }
    auto mesh() const -> const NCPPUPFMesh& { return mesh_; }
    auto localPotential() const -> std::span<const double> { return vloc_; }
    auto nonlocal() const -> const NCPPUPFNonlocal& { return nonlocal_; }
    auto wavefunctions() const -> const NCPPUPFWavefunction& { return wfc_; }
    auto rhoAtom() const -> std::span<const double> { return rho_at_; }

    auto meshType() const -> MeshType;
    auto nonlocalByL(int l) const -> NCPPUPFNonlocalByL;

private:
    NCPPUPFHeader header_;
    NCPPUPFMesh mesh_;
    std::vector<double> vloc_;
    NCPPUPFNonlocal nonlocal_;
    NCPPUPFWavefunction wfc_;
    std::vector<double> rho_at_;

    auto readHeader(const pugi::xml_node& root) -> void;
    auto readMesh(const pugi::xml_node& root) -> void;
    auto readLocalPotential(const pugi::xml_node& root) -> void;
    auto readNonlocal(const pugi::xml_node& root) -> void;
    auto readWavefunctions(const pugi::xml_node& root) -> void;
    auto readRhoAtom(const pugi::xml_node& root) -> void;
};


// --- helper: trim whitespace ---
static auto trimWhitespace(const std::string& s) -> std::string {
    auto start = s.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

// --- helper: parse text content to vector of doubles ---
static auto parseTextToVector(const std::string& text) -> std::vector<double> {
    std::vector<double> data;
    std::istringstream iss(text);
    std::string token;
    while (iss >> token) {
        if (token == "...") continue;
        try {
            std::size_t pos = 0;
            double val = std::stod(token, &pos);
            if (pos == token.size()) {
                data.push_back(val);
            }
        } catch (...) {
            // skip non-numeric tokens
        }
    }
    return data;
}

// --- helper: read string attribute ---
static auto getAttrString(const pugi::xml_node& node, const char* name) -> std::string {
    pugi::xml_attribute attr = node.attribute(name);
    if (!attr) {
        throw std::runtime_error(std::string("UPF: missing attribute '") + name + "'");
    }
    return trimWhitespace(attr.value());
}

// --- helper: read bool attribute (T/F) ---
static auto getAttrBool(const pugi::xml_node& node, const char* name) -> bool {
    std::string value = getAttrString(node, name);
    if (value == "T" || value == "t" || value == "true" || value == "TRUE") return true;
    if (value == "F" || value == "f" || value == "false" || value == "FALSE") return false;
    throw std::runtime_error(std::string("UPF: invalid boolean attribute '") + name + "' = '" + value + "'");
}

// --- helper: read double attribute ---
static auto getAttrDouble(const pugi::xml_node& node, const char* name) -> double {
    std::string value = getAttrString(node, name);
    try {
        return std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("UPF: invalid double attribute '") + name + "' = '" + value + "'");
    }
}

// --- helper: read int attribute ---
static auto getAttrInt(const pugi::xml_node& node, const char* name) -> int {
    std::string value = getAttrString(node, name);
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("UPF: invalid int attribute '") + name + "' = '" + value + "'");
    }
}

// Implementation of NCPPUPF

NCPPUPF::NCPPUPF(const std::string& filename) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if (!result) {
        throw std::runtime_error("UPF: failed to parse file: " + filename + " (" + result.description() + ")");
    }

    pugi::xml_node root = doc.child("UPF");
    if (!root) {
        throw std::runtime_error("UPF: missing root element <UPF> in file: " + filename);
    }

    readHeader(root);
    readMesh(root);
    readLocalPotential(root);
    readNonlocal(root);
    readWavefunctions(root);
    readRhoAtom(root);
}


auto NCPPUPF::readHeader(const pugi::xml_node& root) -> void {
    pugi::xml_node pp_header = root.child("PP_HEADER");
    if (!pp_header) {
        throw std::runtime_error("UPF: missing <PP_HEADER>");
    }

    header_.generated       = getAttrString(pp_header, "generated");
    header_.author          = getAttrString(pp_header, "author");
    header_.date            = getAttrString(pp_header, "date");
    header_.comment         = getAttrString(pp_header, "comment");
    header_.element         = getAttrString(pp_header, "element");
    header_.pseudo_type     = getAttrString(pp_header, "pseudo_type");
    header_.relativistic    = getAttrString(pp_header, "relativistic");
    header_.is_ultrasoft    = getAttrBool(pp_header, "is_ultrasoft");
    header_.is_paw          = getAttrBool(pp_header, "is_paw");
    header_.is_coulomb      = getAttrBool(pp_header, "is_coulomb");
    header_.has_so          = getAttrBool(pp_header, "has_so");
    header_.has_wfc         = getAttrBool(pp_header, "has_wfc");
    header_.has_gipaw       = getAttrBool(pp_header, "has_gipaw");
    header_.core_correction = getAttrBool(pp_header, "core_correction");
    header_.functional      = getAttrString(pp_header, "functional");
    header_.z_valence       = getAttrDouble(pp_header, "z_valence");
    header_.total_psenergy  = getAttrDouble(pp_header, "total_psenergy");
    header_.wfc_cutoff      = getAttrDouble(pp_header, "wfc_cutoff");
    header_.rho_cutoff      = getAttrDouble(pp_header, "rho_cutoff");
    header_.l_max           = getAttrInt(pp_header, "l_max");
    header_.l_local         = getAttrInt(pp_header, "l_local");
    header_.mesh_size       = getAttrInt(pp_header, "mesh_size");
    header_.number_of_wfc   = getAttrInt(pp_header, "number_of_wfc");
    header_.number_of_proj  = getAttrInt(pp_header, "number_of_proj");
}


auto NCPPUPF::readMesh(const pugi::xml_node& root) -> void {
    pugi::xml_node pp_mesh = root.child("PP_MESH");
    if (!pp_mesh) {
        throw std::runtime_error("UPF: missing <PP_MESH>");
    }

    pugi::xml_node r_node = pp_mesh.child("PP_R");
    if (!r_node) throw std::runtime_error("UPF: missing <PP_R>");
    mesh_.r = parseTextToVector(r_node.child_value());
    if (mesh_.r.size() != static_cast<std::size_t>(header_.mesh_size)) {
        throw std::runtime_error("UPF: <PP_R> size mismatch");
    }

    pugi::xml_node rab_node = pp_mesh.child("PP_RAB");
    if (!rab_node) throw std::runtime_error("UPF: missing <PP_RAB>");
    mesh_.rab = parseTextToVector(rab_node.child_value());
    if (mesh_.rab.size() != static_cast<std::size_t>(header_.mesh_size)) {
        throw std::runtime_error("UPF: <PP_RAB> size mismatch");
    }
}


auto NCPPUPF::readLocalPotential(const pugi::xml_node& root) -> void {
    if (header_.is_coulomb) {
        // Coulomb pseudopotential has no local potential
        return;
    }

    pugi::xml_node vloc_node = root.child("PP_LOCAL");
    if (!vloc_node) throw std::runtime_error("UPF: missing <PP_LOCAL>");
    vloc_ = parseTextToVector(vloc_node.child_value());
    if (vloc_.size() != static_cast<std::size_t>(header_.mesh_size)) {
        throw std::runtime_error("UPF: <PP_LOCAL> size mismatch");
    }
}


auto NCPPUPF::readNonlocal(const pugi::xml_node& root) -> void {
    pugi::xml_node pp_nonlocal = root.child("PP_NONLOCAL");
    if (!pp_nonlocal) {
        throw std::runtime_error("UPF: missing <PP_NONLOCAL>");
    }

    const int nb = header_.number_of_proj;
    const int mesh = header_.mesh_size;

    nonlocal_.beta.resize(nb, std::vector<double>(mesh));
    nonlocal_.lll.resize(nb);
    nonlocal_.kbeta.resize(nb);
    nonlocal_.rcut.resize(nb);

    for (int i = 0; i < nb; ++i) {
        std::string tag = "PP_BETA." + std::to_string(i + 1);
        pugi::xml_node beta_node = pp_nonlocal.child(tag.c_str());
        if (!beta_node) {
            throw std::runtime_error("UPF: missing <" + tag + ">");
        }

        nonlocal_.lll[i]   = getAttrInt(beta_node, "angular_momentum");
        nonlocal_.kbeta[i] = getAttrInt(beta_node, "cutoff_radius_index");
        nonlocal_.rcut[i]  = getAttrDouble(beta_node, "cutoff_radius");

        auto data = parseTextToVector(beta_node.child_value());
        if (data.size() != static_cast<std::size_t>(mesh)) {
            throw std::runtime_error("UPF: <" + tag + "> size mismatch");
        }
        // Truncate trailing zeros using cutoff_radius_index (1-based in UPF)
        int cutoff = nonlocal_.kbeta[i];
        if (cutoff < 0 || cutoff > mesh) {
            throw std::runtime_error("UPF: <" + tag + "> invalid cutoff_radius_index");
        }
        data.resize(static_cast<std::size_t>(cutoff));
        nonlocal_.beta[i] = std::move(data);
    }

    pugi::xml_node dij_node = pp_nonlocal.child("PP_DIJ");
    if (!dij_node) throw std::runtime_error("UPF: missing <PP_DIJ>");
    nonlocal_.dion.data = parseTextToVector(dij_node.child_value());
    nonlocal_.dion.rows = nb;
    nonlocal_.dion.cols = nb;
    if (nonlocal_.dion.size() != static_cast<std::size_t>(nb * nb)) {
        throw std::runtime_error("UPF: <PP_DIJ> size mismatch");
    }
}


auto NCPPUPF::readWavefunctions(const pugi::xml_node& root) -> void {
    pugi::xml_node pp_pswfc = root.child("PP_PSWFC");
    if (!pp_pswfc) {
        throw std::runtime_error("UPF: missing <PP_PSWFC>");
    }

    const int nw = header_.number_of_wfc;
    const int mesh = header_.mesh_size;

    wfc_.chi.resize(nw, std::vector<double>(mesh));
    wfc_.kchi.resize(nw);
    wfc_.lchi.resize(nw);
    wfc_.oc.resize(nw);
    wfc_.labels.resize(nw);

    for (int i = 0; i < nw; ++i) {
        std::string tag = "PP_CHI." + std::to_string(i + 1);
        pugi::xml_node chi_node = pp_pswfc.child(tag.c_str());
        if (!chi_node) {
            throw std::runtime_error("UPF: missing <" + tag + ">");
        }

        wfc_.lchi[i]   = getAttrInt(chi_node, "l");
        wfc_.oc[i]     = getAttrDouble(chi_node, "occupation");
        wfc_.labels[i] = getAttrString(chi_node, "label");

        auto data = parseTextToVector(chi_node.child_value());
        if (data.size() != static_cast<std::size_t>(mesh)) {
            throw std::runtime_error("UPF: <" + tag + "> size mismatch");
        }
        // Truncate trailing zeros
        int cutoff = mesh;
        while (cutoff > 0 && data[static_cast<std::size_t>(cutoff) - 1] == 0.0) {
            --cutoff;
        }
        data.resize(static_cast<std::size_t>(cutoff));
        wfc_.kchi[i] = cutoff;
        wfc_.chi[i] = std::move(data);
    }
}


auto NCPPUPF::readRhoAtom(const pugi::xml_node& root) -> void {
    pugi::xml_node rho_node = root.child("PP_RHOATOM");
    if (!rho_node) throw std::runtime_error("UPF: missing <PP_RHOATOM>");
    rho_at_ = parseTextToVector(rho_node.child_value());
    if (rho_at_.size() != static_cast<std::size_t>(header_.mesh_size)) {
        throw std::runtime_error("UPF: <PP_RHOATOM> size mismatch");
    }
}


auto NCPPUPF::meshType() const -> MeshType {
    const auto& r = mesh_.r;
    const auto& rab = mesh_.rab;
    const int n = static_cast<int>(r.size());
    if (n < 2) return MeshType::Unknown;

    // Check uniform: r[i] - r[i-1] ≈ constant
    double drSum = 0.0;
    for (int i = 1; i < n; ++i) {
        drSum += r[i] - r[i - 1];
    }
    double drMean = drSum / (n - 1);
    bool isUniform = true;
    for (int i = 1; i < n; ++i) {
        if (std::abs((r[i] - r[i - 1]) - drMean) > 1e-6 * std::max(std::abs(drMean), 1e-10)) {
            isUniform = false;
            break;
        }
    }
    if (isUniform) return MeshType::Uniform;

    // Check exponential: rab[i] / r[i] ≈ constant (skip r == 0)
    int firstNonZero = 0;
    while (firstNonZero < n && r[firstNonZero] == 0.0) ++firstNonZero;
    if (firstNonZero >= n - 1) return MeshType::Unknown;

    double ratioSum = 0.0;
    int ratioCount = 0;
    for (int i = firstNonZero; i < n; ++i) {
        if (r[i] != 0.0) {
            ratioSum += rab[i] / r[i];
            ++ratioCount;
        }
    }
    if (ratioCount < 2) return MeshType::Unknown;
    double ratioMean = ratioSum / ratioCount;

    bool isExponential = true;
    for (int i = firstNonZero; i < n; ++i) {
        if (r[i] != 0.0) {
            if (std::abs(rab[i] / r[i] - ratioMean) > 1e-6 * std::max(std::abs(ratioMean), 1e-10)) {
                isExponential = false;
                break;
            }
        }
    }
    if (isExponential) return MeshType::Exponential;

    return MeshType::Unknown;
}


auto NCPPUPF::nonlocalByL(int l) const -> NCPPUPFNonlocalByL {
    std::vector<int> indices;
    const int nb = header_.number_of_proj;
    for (int i = 0; i < nb; ++i) {
        if (nonlocal_.lll[i] == l) {
            indices.push_back(i);
        }
    }

    NCPPUPFNonlocalByL result;
    const int nbl = static_cast<int>(indices.size());
    result.beta.resize(nbl);
    for (int i = 0; i < nbl; ++i) {
        result.beta[i] = nonlocal_.beta[indices[i]];
    }

    result.dion.rows = nbl;
    result.dion.cols = nbl;
    result.dion.data.resize(static_cast<std::size_t>(nbl * nbl));
    for (int i = 0; i < nbl; ++i) {
        for (int j = 0; j < nbl; ++j) {
            result.dion[i, j] = nonlocal_.dion[indices[i], indices[j]];
        }
    }

    return result;
}
