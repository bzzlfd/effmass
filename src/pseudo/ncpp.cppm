module;
#include "../physical_constants.hpp"

export module pseudo.ncpp;

import std;
import utils.array2d;
import pseudo.io.upf;

export {
    class NCPP;
    enum class MeshType : int;
    enum class Relativistic : int;
    enum class PseudoType : int;
    struct ProjectorBlock;
}


enum class MeshType {
    Uniform,
    Exponential,
    Unknown
};


enum class Relativistic {
    None,
    Scalar,
    Full
};


enum class PseudoType {
    NC,
    SL,
    US,
    PAW,
    Coulomb
};


struct ProjectorBlock {
    std::vector<std::vector<double>> beta;
    array2d<double> B;
};


class NCPP {
public:
    struct Meta {
        std::string element;
        std::string source_file;
        double z_valence = 0.0;
        Relativistic relativistic = Relativistic::None;
        std::string functional;
        PseudoType pseudo_type = PseudoType::NC;

        int l_max = 0;
        int mesh_size = 0;
        int number_of_wfc = 0;
        int number_of_proj = 0;
        bool has_so = false;
        bool core_correction = false;
        double total_psenergy = 0.0;
        double suggested_ecutwfc_min = 0.0;
        double suggested_ecutrho_min = 0.0;
    };

    struct RadialMesh {
        std::vector<double> r;
        std::vector<double> rab;
        MeshType type = MeshType::Unknown;
    };

    struct VNonlocal {
        std::vector<std::vector<double>> beta;
        std::vector<int> angular_momentum;
        std::vector<int> cutoff_index;
        std::vector<double> cutoff_radius;
        array2d<double> B;
        std::vector<double> jjj;  // total angular momentum j (SOC only)
    };

    struct PseudoWfc {
        std::vector<std::vector<double>> chi;
        std::vector<int> angular_momentum;
        std::vector<int> cutoff_index;
        std::vector<double> cutoff_radius;
        std::vector<double> occupation;
        std::vector<std::string> label;
        std::vector<double> jchi;  // total angular momentum j (SOC only)
    };

    Meta meta;
    RadialMesh mesh;
    std::vector<double> local;
    VNonlocal nonlocal;
    PseudoWfc pseudoWfc;
    std::vector<double> rho_atom;

    // === 别名系统 (Alias) 参见 docs/proposal/005.md ===
    struct Alias {
    private:
        NCPP& p;
        friend class NCPP;
        explicit Alias(NCPP& pp) : p(pp) {}

    public:
        // D/B 矩阵 (UPF: dion)
        auto dion()       -> array2d<double>&       { return p.nonlocal.B; }
        auto dion() const -> const array2d<double>& { return p.nonlocal.B; }

        // 局域势 (UPF: vloc)
        auto vloc()       -> std::vector<double>&       { return p.local; }
        auto vloc() const -> const std::vector<double>& { return p.local; }

        // 投影子 (UPF: beta)
        auto beta()       -> std::vector<std::vector<double>>&       { return p.nonlocal.beta; }
        auto beta() const -> const std::vector<std::vector<double>>& { return p.nonlocal.beta; }
        auto lll()        -> std::vector<int>&       { return p.nonlocal.angular_momentum; }
        auto lll() const  -> const std::vector<int>& { return p.nonlocal.angular_momentum; }
        auto jjj()        -> std::vector<double>&       { return p.nonlocal.jjj; }
        auto jjj() const  -> const std::vector<double>& { return p.nonlocal.jjj; }
        auto kbeta()       -> std::vector<int>&       { return p.nonlocal.cutoff_index; }
        auto kbeta() const -> const std::vector<int>& { return p.nonlocal.cutoff_index; }
        auto rcut()       -> std::vector<double>&       { return p.nonlocal.cutoff_radius; }
        auto rcut() const -> const std::vector<double>& { return p.nonlocal.cutoff_radius; }

        // 赝波函数 (UPF: chi / wfc)
        auto chi()       -> std::vector<std::vector<double>>&       { return p.pseudoWfc.chi; }
        auto chi() const -> const std::vector<std::vector<double>>& { return p.pseudoWfc.chi; }
        auto wfc()       -> std::vector<std::vector<double>>&       { return p.pseudoWfc.chi; }
        auto wfc() const -> const std::vector<std::vector<double>>& { return p.pseudoWfc.chi; }
        auto lchi()       -> std::vector<int>&       { return p.pseudoWfc.angular_momentum; }
        auto lchi() const -> const std::vector<int>& { return p.pseudoWfc.angular_momentum; }
        auto jchi()       -> std::vector<double>&       { return p.pseudoWfc.jchi; }
        auto jchi() const -> const std::vector<double>& { return p.pseudoWfc.jchi; }
        auto oc()       -> std::vector<double>&       { return p.pseudoWfc.occupation; }
        auto oc() const -> const std::vector<double>& { return p.pseudoWfc.occupation; }

        // 文献/通用别名
        auto D()           -> array2d<double>&       { return p.nonlocal.B; }
        auto D() const     -> const array2d<double>& { return p.nonlocal.B; }
        auto projectors()       -> std::vector<std::vector<double>>&       { return p.nonlocal.beta; }
        auto projectors() const -> const std::vector<std::vector<double>>& { return p.nonlocal.beta; }
        auto radialGrid()       -> std::vector<double>&       { return p.mesh.r; }
        auto radialGrid() const -> const std::vector<double>& { return p.mesh.r; }
    };

    Alias alias{*this};

    explicit NCPP(const UPF& upf);

    auto inferMeshType() -> MeshType;
    auto projectorBlock(int l) const -> ProjectorBlock;

private:
    auto checkSemilocal(const UPF& upf) -> void;
    auto checkConsistency() -> void;
};


NCPP::NCPP(const UPF& upf) {
    checkSemilocal(upf);

    const auto& h = upf.header();

    if (h.pseudo_type != "NC") {
        throw std::runtime_error("NCPP: expected NC pseudopotential, got " + h.pseudo_type);
    }

    const auto& nl = upf.nonlocal();
    const auto& wf = upf.wavefunctions();

    meta.element     = h.element;
    meta.source_file = upf.sourceFile();
    meta.z_valence   = h.z_valence;
    meta.functional = h.functional;

    if (h.relativistic == "no")         meta.relativistic = Relativistic::None;
    else if (h.relativistic == "scalar")   meta.relativistic = Relativistic::Scalar;
    else if (h.relativistic == "full")     meta.relativistic = Relativistic::Full;

    meta.l_max            = h.l_max;
    meta.mesh_size        = h.mesh_size;
    meta.number_of_wfc    = h.number_of_wfc;
    meta.number_of_proj   = h.number_of_proj;
    meta.has_so           = h.has_so;
    meta.core_correction  = h.core_correction;
    meta.total_psenergy   = convertEnergy(h.total_psenergy, EnergyUnit::Rydberg, EnergyUnit::Hartree);
    meta.suggested_ecutwfc_min = convertEnergy(h.wfc_cutoff, EnergyUnit::Rydberg, EnergyUnit::Hartree);
    meta.suggested_ecutrho_min = convertEnergy(h.rho_cutoff, EnergyUnit::Rydberg, EnergyUnit::Hartree);

    mesh.r   = upf.mesh().r;
    mesh.rab = upf.mesh().rab;

    const auto vloc_src = upf.localPotential();
    local.assign(vloc_src.begin(), vloc_src.end());
    for (auto& v : local) v = convertEnergy(v, EnergyUnit::Rydberg, EnergyUnit::Hartree);

    nonlocal.beta             = nl.beta;
    nonlocal.angular_momentum = nl.lll;
    nonlocal.cutoff_index     = nl.kbeta;
    nonlocal.cutoff_radius    = nl.rcut;
    nonlocal.B                = nl.dion;
    for (auto& d : nonlocal.B.data) d = convertEnergy(d, EnergyUnit::Rydberg, EnergyUnit::Hartree);

    // Truncate beta to effective length (cutoff_radius_index)
    for (int i = 0; i < meta.number_of_proj; ++i) {
        int cutoff = nonlocal.cutoff_index[i];
        nonlocal.beta[i].resize(static_cast<std::size_t>(cutoff));
    }

    pseudoWfc.angular_momentum = wf.lchi;
    pseudoWfc.occupation       = wf.oc;
    pseudoWfc.label            = wf.labels;
    pseudoWfc.chi              = wf.chi;
    pseudoWfc.cutoff_index.resize(meta.number_of_wfc);
    pseudoWfc.cutoff_radius.resize(meta.number_of_wfc);

    // Truncate chi trailing zeros
    for (int i = 0; i < meta.number_of_wfc; ++i) {
        int cutoff = static_cast<int>(pseudoWfc.chi[i].size());
        while (cutoff > 0 && pseudoWfc.chi[i][static_cast<std::size_t>(cutoff) - 1] == 0.0) {
            --cutoff;
        }
        pseudoWfc.chi[i].resize(static_cast<std::size_t>(cutoff));
        pseudoWfc.cutoff_index[i] = cutoff;
        pseudoWfc.cutoff_radius[i] = mesh.r[static_cast<std::size_t>(cutoff) - 1];
    }

    const auto rho_src = upf.rhoAtom();
    rho_atom.assign(rho_src.begin(), rho_src.end());

    if (h.has_so) {
        const auto* soc = upf.socData();
        nonlocal.jjj = soc->jjj;
        pseudoWfc.jchi = soc->jchi;
    }

    checkConsistency();
}


auto NCPP::checkSemilocal(const UPF& upf) -> void {
    if (upf.header().pseudo_type == "SL") {
        throw std::runtime_error(
            "NCPP: semilocal NC pseudopotential not yet supported");
    }
}

auto NCPP::checkConsistency() -> void {
    int ms = meta.mesh_size;
    int nb = meta.number_of_proj;
    int nw = meta.number_of_wfc;

    if (ms <= 0)
        throw std::runtime_error("NCPP: mesh_size (" + std::to_string(ms) + ") must be positive");

    if (nb < 0)
        throw std::runtime_error("NCPP: number_of_proj (" + std::to_string(nb) + ") must be non-negative");
    if (nw < 0)
        throw std::runtime_error("NCPP: number_of_wfc (" + std::to_string(nw) + ") must be non-negative");

    if (static_cast<int>(mesh.r.size()) != ms)
        throw std::runtime_error("NCPP: mesh.r size (" + std::to_string(mesh.r.size())
            + ") != mesh_size (" + std::to_string(ms) + ")");
    if (static_cast<int>(mesh.rab.size()) != ms)
        throw std::runtime_error("NCPP: mesh.rab size != mesh_size");
    if (static_cast<int>(local.size()) != ms)
        throw std::runtime_error("NCPP: local potential size != mesh_size");
    if (static_cast<int>(rho_atom.size()) != ms)
        throw std::runtime_error("NCPP: rho_atom size != mesh_size");

    if (nb > 0) {
        int bsz = static_cast<int>(nonlocal.beta.size());
        if (bsz != nb)
            throw std::runtime_error("NCPP: beta.size() (" + std::to_string(bsz)
                + ") != number_of_proj (" + std::to_string(nb) + ")");
        if (static_cast<int>(nonlocal.angular_momentum.size()) != nb)
            throw std::runtime_error("NCPP: angular_momentum.size() != number_of_proj");
        if (static_cast<int>(nonlocal.cutoff_index.size()) != nb)
            throw std::runtime_error("NCPP: cutoff_index.size() != number_of_proj");
        if (static_cast<int>(nonlocal.cutoff_radius.size()) != nb)
            throw std::runtime_error("NCPP: cutoff_radius.size() != number_of_proj");
        if (nonlocal.B.rows != nb || nonlocal.B.cols != nb)
            throw std::runtime_error("NCPP: B matrix (" + std::to_string(nonlocal.B.rows)
                + "x" + std::to_string(nonlocal.B.cols) + ") != " + std::to_string(nb) + "x" + std::to_string(nb));

        for (int i = 0; i < nb; ++i) {
            int ci = nonlocal.cutoff_index[i];
            if (ci < 0 || ci > ms)
                throw std::runtime_error("NCPP: cutoff_index[" + std::to_string(i) + "] = "
                    + std::to_string(ci) + " out of range [0, " + std::to_string(ms) + "]");
            if (static_cast<int>(nonlocal.beta[i].size()) != ci)
                throw std::runtime_error("NCPP: beta[" + std::to_string(i) + "] size ("
                    + std::to_string(nonlocal.beta[i].size()) + ") != cutoff_index (" + std::to_string(ci) + ")");
        }
    }

    if (nw > 0) {
        if (static_cast<int>(pseudoWfc.chi.size()) != nw)
            throw std::runtime_error("NCPP: chi.size() != number_of_wfc");
        if (static_cast<int>(pseudoWfc.angular_momentum.size()) != nw)
            throw std::runtime_error("NCPP: angular_momentum(lchi).size() != number_of_wfc");
        if (static_cast<int>(pseudoWfc.occupation.size()) != nw)
            throw std::runtime_error("NCPP: occupation.size() != number_of_wfc");
        if (static_cast<int>(pseudoWfc.label.size()) != nw)
            throw std::runtime_error("NCPP: label.size() != number_of_wfc");
    }

    if (meta.has_so) {
        if (static_cast<int>(nonlocal.jjj.size()) != nb)
            throw std::runtime_error("NCPP: jjj.size() (" + std::to_string(nonlocal.jjj.size())
                + ") != number_of_proj (" + std::to_string(nb) + ")");
        if (static_cast<int>(pseudoWfc.jchi.size()) != nw)
            throw std::runtime_error("NCPP: jchi.size() (" + std::to_string(pseudoWfc.jchi.size())
                + ") != number_of_wfc (" + std::to_string(nw) + ")");
    }
}


auto NCPP::inferMeshType() -> MeshType {
    const int n = static_cast<int>(mesh.r.size());
    if (n < 2) {
        mesh.type = MeshType::Unknown;
        return mesh.type;
    }

    double drSum = 0.0;
    for (int i = 1; i < n; ++i) {
        drSum += mesh.r[i] - mesh.r[i - 1];
    }
    double drMean = drSum / (n - 1);
    bool isUniform = true;
    for (int i = 1; i < n; ++i) {
        if (std::abs((mesh.r[i] - mesh.r[i - 1]) - drMean) > 1e-6 * std::max(std::abs(drMean), 1e-10)) {
            isUniform = false;
            break;
        }
    }
    if (isUniform) {
        mesh.type = MeshType::Uniform;
        return mesh.type;
    }

    int firstNonZero = 0;
    while (firstNonZero < n && mesh.r[firstNonZero] == 0.0) ++firstNonZero;
    if (firstNonZero >= n - 1) {
        mesh.type = MeshType::Unknown;
        return mesh.type;
    }

    double ratioSum = 0.0;
    int ratioCount = 0;
    for (int i = firstNonZero; i < n; ++i) {
        if (mesh.r[i] != 0.0) {
            ratioSum += mesh.rab[i] / mesh.r[i];
            ++ratioCount;
        }
    }
    if (ratioCount < 2) {
        mesh.type = MeshType::Unknown;
        return mesh.type;
    }
    double ratioMean = ratioSum / ratioCount;

    bool isExponential = true;
    for (int i = firstNonZero; i < n; ++i) {
        if (mesh.r[i] != 0.0) {
            if (std::abs(mesh.rab[i] / mesh.r[i] - ratioMean) > 1e-6 * std::max(std::abs(ratioMean), 1e-10)) {
                isExponential = false;
                break;
            }
        }
    }
    if (isExponential) {
        mesh.type = MeshType::Exponential;
        return mesh.type;
    }

    mesh.type = MeshType::Unknown;
    return mesh.type;
}


auto NCPP::projectorBlock(int l) const -> ProjectorBlock {
    const int nb = meta.number_of_proj;

    std::vector<int> indices;
    for (int i = 0; i < nb; ++i) {
        if (nonlocal.angular_momentum[i] == l) {
            indices.push_back(i);
        }
    }

    ProjectorBlock result;
    const int nbl = static_cast<int>(indices.size());
    result.beta.resize(nbl);
    for (int i = 0; i < nbl; ++i) {
        result.beta[i] = nonlocal.beta[indices[i]];
    }

    result.B.rows = nbl;
    result.B.cols = nbl;
    result.B.data.resize(static_cast<std::size_t>(nbl * nbl));
    for (int i = 0; i < nbl; ++i) {
        for (int j = 0; j < nbl; ++j) {
            result.B[i, j] = nonlocal.B[indices[i], indices[j]];
        }
    }

    return result;
}


