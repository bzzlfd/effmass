export module pseudo.ncpp;

import std;
import utils.matrix;
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
    DenseMatrix<double> B;
};


class NCPP {
public:
    struct Meta {
        std::string element;
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
        double suggested_ecutwfc_min = 0.0;  // Ry
        double suggested_ecutrho_min = 0.0;  // Ry
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
        DenseMatrix<double> B;
    };

    struct PseudoWfc {
        std::vector<std::vector<double>> chi;
        std::vector<int> angular_momentum;
        std::vector<int> cutoff_index;
        std::vector<double> cutoff_radius;
        std::vector<double> occupation;
        std::vector<std::string> label;
    };

    Meta meta;
    RadialMesh mesh;
    std::vector<double> local;
    VNonlocal nonlocal;
    PseudoWfc pseudoWfc;
    std::vector<double> rho_atom;

    explicit NCPP(const UPF& upf);

    auto inferMeshType() -> MeshType;
    auto projectorBlock(int l) const -> ProjectorBlock;

private:
    auto checkSemilocal(const UPF& upf) -> void;
};


NCPP::NCPP(const UPF& upf) {
    checkSemilocal(upf);

    const auto& h = upf.header();

    meta.element    = h.element;
    meta.z_valence  = h.z_valence;
    meta.functional = h.functional;

    if (h.relativistic == "no")         meta.relativistic = Relativistic::None;
    else if (h.relativistic == "scalar")   meta.relativistic = Relativistic::Scalar;
    else if (h.relativistic == "full")     meta.relativistic = Relativistic::Full;

    if (h.pseudo_type != "NC") {
        throw std::runtime_error("NCPP: expected NC pseudopotential, got " + h.pseudo_type);
    }

    meta.l_max            = h.l_max;
    meta.mesh_size        = h.mesh_size;
    meta.number_of_wfc    = h.number_of_wfc;
    meta.number_of_proj   = h.number_of_proj;
    meta.has_so           = h.has_so;
    meta.core_correction  = h.core_correction;
    meta.total_psenergy   = h.total_psenergy;
    meta.suggested_ecutwfc_min = h.wfc_cutoff;
    meta.suggested_ecutrho_min = h.rho_cutoff;

    mesh.r   = upf.mesh().r;
    mesh.rab = upf.mesh().rab;

    const auto vloc_src = upf.localPotential();
    local.assign(vloc_src.begin(), vloc_src.end());

    const auto& nl = upf.nonlocal();
    nonlocal.beta             = nl.beta;
    nonlocal.angular_momentum = nl.lll;
    nonlocal.cutoff_index     = nl.kbeta;
    nonlocal.cutoff_radius    = nl.rcut;
    nonlocal.B                = nl.dion;

    // Truncate beta to effective length (cutoff_radius_index)
    for (int i = 0; i < meta.number_of_proj; ++i) {
        int cutoff = nonlocal.cutoff_index[i];
        nonlocal.beta[i].resize(static_cast<std::size_t>(cutoff));
    }

    const auto& wf = upf.wavefunctions();
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
}


auto NCPP::checkSemilocal(const UPF& upf) -> void {
    if (upf.header().pseudo_type == "SL") {
        throw std::runtime_error(
            "NCPP: semilocal NC pseudopotential not yet supported");
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
