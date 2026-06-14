module;
#include "pugixml.hpp"

export module pseudo.io.upf;

import std;
import utils.array2d;


export {
    class UPF;
        struct UPFHeader;
        struct UPFMesh;
        struct UPFNonlocal;
        struct UPFWavefunction;
        struct UPFSOC;
}


// === UPF Header structure ===
struct UPFHeader {
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
struct UPFMesh {
    std::vector<double> r;
    std::vector<double> rab;
};


// === Nonlocal data (beta projectors + D_ij) ===
struct UPFNonlocal {
    std::vector<std::vector<double>> beta;
    std::vector<int> lll;
    std::vector<int> kbeta;
    std::vector<double> rcut;
    array2d<double> dion;
};


// === Wavefunction data (pseudo atomic orbitals) ===
struct UPFWavefunction {
    std::vector<std::vector<double>> chi;
    std::vector<int> lchi;
    std::vector<double> oc;
    std::vector<std::string> labels;
};


// === SOC extension data (only populated when has_so == true) ===
struct UPFSOC {
    enum class Format { pp_spin_orb, ps_library };

    Format format = Format::pp_spin_orb;
    std::vector<double> jjj;   // [nbeta] total angular momentum
    std::vector<double> jchi;  // [nwfc]  total angular momentum
};


// === UPF class ===
class UPF {
public:
    explicit UPF(const std::string& filename);

    auto header() const -> const UPFHeader& { return header_; }
    auto mesh() const -> const UPFMesh& { return mesh_; }
    auto localPotential() const -> std::span<const double> { return vloc_; }
    auto nonlocal() const -> const UPFNonlocal& { return nonlocal_; }
    auto wavefunctions() const -> const UPFWavefunction& { return wfc_; }
    auto rhoAtom() const -> std::span<const double> { return rho_at_; }
    auto socData() const -> const UPFSOC* { return header_.has_so ? &soc_ : nullptr; }
    auto sourceFile() const -> const std::string& { return filename_; }

private:
    UPFHeader header_;
    UPFMesh mesh_;
    std::vector<double> vloc_;
    UPFNonlocal nonlocal_;
    UPFWavefunction wfc_;
    std::vector<double> rho_at_;
    UPFSOC soc_;
    std::string filename_;

    auto readSpinOrbit(const pugi::xml_node& root) -> void;
    auto readHeader(const pugi::xml_node& root) -> void;
    auto readMesh(const pugi::xml_node& root) -> void;
    auto readLocalPotential(const pugi::xml_node& root) -> void;
    auto readNonlocal(const pugi::xml_node& root) -> void;
    auto readWavefunctions(const pugi::xml_node& root) -> void;
    auto readRhoAtom(const pugi::xml_node& root) -> void;
};


namespace {
    auto trimWhitespace(const std::string& s) -> std::string {
        auto start = s.find_first_not_of(" \t\n\r");
        if (start == std::string::npos) return "";
        auto end = s.find_last_not_of(" \t\n\r");
        return s.substr(start, end - start + 1);
    }

    auto parseTextToVector(const std::string& text) -> std::vector<double> {
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

    auto getAttrString(const pugi::xml_node& node, const char* name) -> std::string {
        pugi::xml_attribute attr = node.attribute(name);
        if (!attr) {
            throw std::runtime_error(std::string("UPF: missing attribute '") + name + "'");
        }
        return trimWhitespace(attr.value());
    }

    auto getAttrBool(const pugi::xml_node& node, const char* name) -> bool {
        std::string value = getAttrString(node, name);
        if (value == "T" || value == "t" || value == "true" || value == "TRUE") return true;
        if (value == "F" || value == "f" || value == "false" || value == "FALSE") return false;
        throw std::runtime_error(std::string("UPF: invalid boolean attribute '") + name + "' = '" + value + "'");
    }

    auto getAttrDouble(const pugi::xml_node& node, const char* name) -> double {
        std::string value = getAttrString(node, name);
        try {
            return std::stod(value);
        } catch (const std::exception&) {
            throw std::runtime_error(std::string("UPF: invalid double attribute '") + name + "' = '" + value + "'");
        }
    }

    auto getAttrInt(const pugi::xml_node& node, const char* name) -> int {
        std::string value = getAttrString(node, name);
        try {
            return std::stoi(value);
        } catch (const std::exception&) {
            throw std::runtime_error(std::string("UPF: invalid int attribute '") + name + "' = '" + value + "'");
        }
    }
}

// Implementation of UPF

UPF::UPF(const std::string& filename) : filename_(filename) {
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

    if (header_.is_ultrasoft) {
        throw std::runtime_error("UPF: USPP reader not yet implemented");
    }
    if (header_.is_paw) {
        throw std::runtime_error("UPF: PAW reader not yet implemented");
    }

    readMesh(root);
    readLocalPotential(root);
    readNonlocal(root);
    readWavefunctions(root);
    readRhoAtom(root);

    if (header_.has_so) {
        readSpinOrbit(root);
    }
}


auto UPF::readHeader(const pugi::xml_node& root) -> void {
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
    if (pp_header.attribute("wfc_cutoff")) {
        header_.wfc_cutoff = getAttrDouble(pp_header, "wfc_cutoff");
    }
    if (pp_header.attribute("rho_cutoff")) {
        header_.rho_cutoff = getAttrDouble(pp_header, "rho_cutoff");
    }
    header_.l_max           = getAttrInt(pp_header, "l_max");
    header_.l_local         = getAttrInt(pp_header, "l_local");
    header_.mesh_size       = getAttrInt(pp_header, "mesh_size");
    header_.number_of_wfc   = getAttrInt(pp_header, "number_of_wfc");
    header_.number_of_proj  = getAttrInt(pp_header, "number_of_proj");
}


auto UPF::readMesh(const pugi::xml_node& root) -> void {
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


auto UPF::readLocalPotential(const pugi::xml_node& root) -> void {
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


auto UPF::readNonlocal(const pugi::xml_node& root) -> void {
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


auto UPF::readWavefunctions(const pugi::xml_node& root) -> void {
    pugi::xml_node pp_pswfc = root.child("PP_PSWFC");
    if (!pp_pswfc) {
        throw std::runtime_error("UPF: missing <PP_PSWFC>");
    }

    const int nw = header_.number_of_wfc;
    const int mesh = header_.mesh_size;

    wfc_.chi.resize(nw, std::vector<double>(mesh));
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
        wfc_.chi[i] = std::move(data);
    }
}


auto UPF::readRhoAtom(const pugi::xml_node& root) -> void {
    pugi::xml_node rho_node = root.child("PP_RHOATOM");
    if (!rho_node) throw std::runtime_error("UPF: missing <PP_RHOATOM>");
    rho_at_ = parseTextToVector(rho_node.child_value());
    if (rho_at_.size() != static_cast<std::size_t>(header_.mesh_size)) {
        throw std::runtime_error("UPF: <PP_RHOATOM> size mismatch");
    }
}


auto UPF::readSpinOrbit(const pugi::xml_node& root) -> void {
    const int nb = header_.number_of_proj;
    const int nw = header_.number_of_wfc;

    soc_.jjj .resize(nb);
    soc_.jchi.resize(nw);

    pugi::xml_node so = root.child("PP_SPIN_ORB");

    if (so) {
        // PseudoDojo format: jjj/jchi in PP_SPIN_ORB/PP_RELBETA.* and PP_RELWFC.*
        soc_.format = UPFSOC::Format::pp_spin_orb;

        for (int i = 0; i < nb; ++i) {
            std::string tag = "PP_RELBETA." + std::to_string(i + 1);
            pugi::xml_node rb = so.child(tag.c_str());
            if (!rb) {
                throw std::runtime_error("UPF: missing <" + tag + "> in PP_SPIN_ORB");
            }
            soc_.jjj[i] = getAttrDouble(rb, "jjj");
        }

        for (int i = 0; i < nw; ++i) {
            std::string tag = "PP_RELWFC." + std::to_string(i + 1);
            pugi::xml_node rw = so.child(tag.c_str());
            if (!rw) {
                throw std::runtime_error("UPF: missing <" + tag + "> in PP_SPIN_ORB");
            }
            soc_.jchi[i] = getAttrDouble(rw, "jchi");
        }
    } else {
        // PSlibrary format: jjj on PP_BETA.* / jchi on PP_CHI.*
        soc_.format = UPFSOC::Format::ps_library;

        pugi::xml_node pp_nonlocal = root.child("PP_NONLOCAL");
        if (!pp_nonlocal) {
            throw std::runtime_error("UPF: missing <PP_NONLOCAL> for SOC pseudopotential");
        }
        for (int i = 0; i < nb; ++i) {
            std::string tag = "PP_BETA." + std::to_string(i + 1);
            pugi::xml_node beta_node = pp_nonlocal.child(tag.c_str());
            if (!beta_node) {
                throw std::runtime_error("UPF: missing <" + tag + ">");
            }
            if (!beta_node.attribute("jjj")) {
                throw std::runtime_error("UPF: missing jjj attribute on <" + tag + ">");
            }
            soc_.jjj[i] = getAttrDouble(beta_node, "jjj");
        }

        pugi::xml_node pp_pswfc = root.child("PP_PSWFC");
        if (!pp_pswfc) {
            throw std::runtime_error("UPF: missing <PP_PSWFC> for SOC pseudopotential");
        }
        for (int i = 0; i < nw; ++i) {
            std::string tag = "PP_CHI." + std::to_string(i + 1);
            pugi::xml_node chi_node = pp_pswfc.child(tag.c_str());
            if (!chi_node) {
                throw std::runtime_error("UPF: missing <" + tag + ">");
            }
            if (!chi_node.attribute("jchi")) {
                throw std::runtime_error("UPF: missing jchi attribute on <" + tag + ">");
            }
            soc_.jchi[i] = getAttrDouble(chi_node, "jchi");
        }
    }
}
