export module pseudo.ncpp;

import std;
import utils.array2d;
import math.linalg;
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

    struct Advance {
    private:
        NCPP& p;
        friend class NCPP;
        explicit Advance(NCPP& pp) : p(pp) {}
    public:
        auto diagonalizeNonlocal() -> void;
        auto upf2upfSO() -> void;
    };

    Advance advance{*this};

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

    meta.element    = h.element;
    meta.z_valence  = h.z_valence;
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
    meta.total_psenergy   = h.total_psenergy;
    meta.suggested_ecutwfc_min = h.wfc_cutoff;
    meta.suggested_ecutrho_min = h.rho_cutoff;

    mesh.r   = upf.mesh().r;
    mesh.rab = upf.mesh().rab;

    const auto vloc_src = upf.localPotential();
    local.assign(vloc_src.begin(), vloc_src.end());

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


auto NCPP::Advance::diagonalizeNonlocal() -> void {
    const int n = p.meta.number_of_proj;
    if (n == 0) return;

    // Diagonalize B = U * Lambda * U^T
    auto result = diagonalize_jacobi(p.nonlocal.B.data, n);
    if (!result.converged) {
        throw std::runtime_error("NCPP::Advance::diagonalizeNonlocal: Jacobi diagonalization did not converge");
    }

    // B <- Lambda (diagonal)
    std::ranges::fill(p.nonlocal.B.data, 0.0);
    for (int i = 0; i < n; ++i) {
        p.nonlocal.B[i, i] = result.eigenvalues[static_cast<std::size_t>(i)];
    }

    // beta_new[k][r] = sum_i beta[i][r] * U[i][k]
    // eigenvectors row-major, row k <-> eigenvector k -> U[i][k] = eigvec[k, i]
    int max_cutoff = 0;
    for (const auto& b : p.nonlocal.beta) {
        max_cutoff = std::max(max_cutoff, static_cast<int>(b.size()));
    }

    auto beta_old = std::move(p.nonlocal.beta);
    p.nonlocal.beta.resize(n);
    for (int k = 0; k < n; ++k) {
        p.nonlocal.beta[k].assign(max_cutoff, 0.0);
        for (int i = 0; i < n; ++i) {
            double u_ik = result.eigenvectors[k][i];
            const auto& beta_i = beta_old[i];
            for (int r = 0; r < static_cast<int>(beta_i.size()); ++r) {
                p.nonlocal.beta[k][r] += beta_i[r] * u_ik;
            }
        }
    }

    // Update cutoff_index and cutoff_radius to max cutoff
    std::ranges::fill(p.nonlocal.cutoff_index, max_cutoff);
    if (max_cutoff > 0 && !p.mesh.r.empty()) {
        std::ranges::fill(p.nonlocal.cutoff_radius, p.mesh.r[static_cast<std::size_t>(max_cutoff - 1)]);
    }
}


auto NCPP::Advance::upf2upfSO() -> void {
    if (!p.meta.has_so) {
        throw std::runtime_error("NCPP::Advance::upf2upfSO: has_so must be true");
    }
    if (p.meta.l_max < 1) {
        throw std::runtime_error("NCPP::Advance::upf2upfSO: l_max must be >= 1");
    }
    if (p.meta.l_max > 3) {
        throw std::runtime_error("NCPP::Advance::upf2upfSO: l_max > 3 not yet supported");
    }

    const int n = p.meta.mesh_size;
    const int nbeta = p.meta.number_of_proj;
    const int nwave = p.meta.number_of_wfc;
    const int lmax = p.meta.l_max;
    const auto& r = p.mesh.r;

    // Integration weights: dr[i] = (r[i+1] - r[i-1]) / 2,  i = 1..n-2
    std::vector<double> dr(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        dr[i] = (r[i + 1] - r[i - 1]) / 2.0;
    }

    // Pad beta to full mesh
    std::vector<std::vector<double>> beta_tmp(nbeta, std::vector<double>(n, 0.0));
    for (int ll = 0; ll < nbeta; ++ll) {
        const auto& src = p.nonlocal.beta[ll];
        std::copy(src.begin(), src.end(), beta_tmp[ll].begin());
    }

    // Pad chi to full mesh
    std::vector<std::vector<double>> wave_tmp(nwave, std::vector<double>(n, 0.0));
    for (int ll = 0; ll < nwave; ++ll) {
        const auto& src = p.pseudoWfc.chi[ll];
        std::copy(src.begin(), src.end(), wave_tmp[ll].begin());
    }

    // DD diagonal access
    auto dd_diag = [&](int ll) -> double {
        if (ll < 0 || ll >= nbeta) return 0.0;
        return p.nonlocal.B.data[static_cast<std::size_t>(ll) * nbeta + ll];
    };

    // ------------------------------------------------------------
    // 1. Map projectors to (L, J, m) grouped structure
    // ------------------------------------------------------------
    std::vector<int> map_L(nbeta), map_J(nbeta);
    for (int ll = 0; ll < nbeta; ++ll) {
        map_L[ll] = p.nonlocal.angular_momentum[ll];
        map_J[ll] = (p.nonlocal.jjj[ll] - map_L[ll] < 0) ? 0 : 1;
    }

    // beta_arr[L][J][m][i] -- 0-based L=0..lmax, J=0,1, m=0,1
    using Beta4D = std::vector<std::vector<std::vector<std::vector<double>>>>;
    Beta4D beta_arr(lmax + 1,
        std::vector<std::vector<std::vector<double>>>(2,
            std::vector<std::vector<double>>(2,
                std::vector<double>(n, 0.0))));
    std::vector<std::vector<int>> ind(lmax + 1, std::vector<int>(2, 0));
    std::vector<std::vector<std::vector<int>>> ll_map(lmax + 1,
        std::vector<std::vector<int>>(2, std::vector<int>(2, -1)));

    for (int ll = 0; ll < nbeta; ++ll) {
        int L = map_L[ll];
        int J = map_J[ll];
        int m = ind[L][J];       // 0-based m, before increment
        ll_map[L][J][m] = ll;
        for (int i = 0; i < n; ++i) {
            beta_arr[L][J][m][i] = beta_tmp[ll][i];
        }
        ind[L][J]++;
    }

    // Duplicate: if only 1 projector for (L,J), copy m=0 -> m=1
    for (int L = 0; L <= lmax; ++L) {
        for (int J = 0; J <= 1; ++J) {
            if (ind[L][J] == 1) {
                for (int i = 0; i < n; ++i) {
                    beta_arr[L][J][1][i] = beta_arr[L][J][0][i];
                }
                ind[L][J] = 2;
            }
        }
    }

    // ------------------------------------------------------------
    // 2. Map wavefunctions to (L, J) grouped structure
    // ------------------------------------------------------------
    std::vector<std::vector<std::vector<double>>> wave_arr(lmax + 1,
        std::vector<std::vector<double>>(2, std::vector<double>(n, 0.0)));
    std::vector<std::vector<double>> occ(lmax + 1, std::vector<double>(2, 0.0));
    std::vector<int> map_Lw(nwave), map_Jw(nwave);
    int lchi_max = 0;
    for (int ll = 0; ll < nwave; ++ll) {
        int Lw = p.pseudoWfc.angular_momentum[ll];
        int Jw = (p.pseudoWfc.jchi[ll] - Lw < 0) ? 0 : 1;
        map_Lw[ll] = Lw;
        map_Jw[ll] = Jw;
        lchi_max = std::max(lchi_max, Lw);
        for (int i = 0; i < n; ++i) {
            wave_arr[Lw][Jw][i] = wave_tmp[ll][i];
        }
        occ[Lw][Jw] = p.pseudoWfc.occupation[ll];
    }

    // ------------------------------------------------------------
    // 3. Compute beta0 and betaSO for each L = 1..lmax
    // ------------------------------------------------------------
    std::vector<std::vector<std::vector<double>>> beta0_scalar(lmax + 1,
        std::vector<std::vector<double>>(2, std::vector<double>(n, 0.0)));
    std::vector<std::vector<double>> beta_SO(lmax + 1, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> DD_0(lmax + 1, std::vector<double>(2, 0.0));
    std::vector<double> DD_SO(lmax + 1, 0.0);

    for (int L = 1; L <= lmax; ++L) {
        // Overlap sums
        double sum1 = 0.0, sum2 = 0.0, sumw = 0.0, sumw1 = 0.0, sumw2 = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            sum1  += beta_arr[L][0][0][i] * beta_arr[L][1][0][i] * dr[i];
            sum2  += beta_arr[L][0][1][i] * beta_arr[L][1][1][i] * dr[i];
            sumw  += wave_arr[L][0][i]     * wave_arr[L][1][i]     * dr[i];
            sumw1 += wave_arr[L][0][i]     * wave_arr[L][0][i]     * dr[i];
            sumw2 += wave_arr[L][1][i]     * wave_arr[L][1][i]     * dr[i];
        }

        int isign1 = (sum1 > 0) ? 1 : -1;
        int isign2 = (sum2 > 0) ? 1 : -1;
        int isignw = (sumw > 0) ? 1 : -1;

        // Mixing coefficients
        double x2 = std::sqrt(L + 1.0) / (std::sqrt(L + 1.0) + std::sqrt(L * 1.0));
        double x1 = 1.0 - x2;

        // Integrals with test functions
        double s1_1 = 0.0, s2_1 = 0.0, s1_2 = 0.0, s2_2 = 0.0;
        double s1_1T = 0.0, s2_1T = 0.0, s1_2T = 0.0, s2_2T = 0.0;
        int lw_ref = (L <= lchi_max) ? L : lchi_max;

        for (int i = 1; i < n - 1; ++i) {
            double waveT1 = x2 * beta_arr[L][1][0][i] + x1 * static_cast<double>(isign1) * beta_arr[L][0][0][i];
            double waveT2 = x2 * beta_arr[L][1][1][i] + x1 * static_cast<double>(isign2) * beta_arr[L][0][1][i];
            double waveT  = x2 * wave_arr[lw_ref][1][i] + x1 * static_cast<double>(isignw) * wave_arr[lw_ref][0][i];

            s1_1 += beta_arr[L][0][0][i] * waveT * dr[i];
            s2_1 += beta_arr[L][1][0][i] * waveT * dr[i];
            s1_2 += beta_arr[L][0][1][i] * waveT * dr[i];
            s2_2 += beta_arr[L][1][1][i] * waveT * dr[i];

            s1_1T += beta_arr[L][0][0][i] * waveT1 * dr[i];
            s2_1T += beta_arr[L][1][0][i] * waveT1 * dr[i];
            s1_2T += beta_arr[L][0][1][i] * waveT2 * dr[i];
            s2_2T += beta_arr[L][1][1][i] * waveT2 * dr[i];
        }

        // Scale by D diagonal elements
        s1_1 *= dd_diag(ll_map[L][0][0]);
        s2_1 *= dd_diag(ll_map[L][1][0]);
        s1_2 *= dd_diag(ll_map[L][0][1]);
        s2_2 *= dd_diag(ll_map[L][1][1]);
        s1_1T *= dd_diag(ll_map[L][0][0]);
        s2_1T *= dd_diag(ll_map[L][1][0]);
        s1_2T *= dd_diag(ll_map[L][0][1]);
        s2_2T *= dd_diag(ll_map[L][1][1]);

        // Construct beta0 (scalar-relativistic projectors)
        double Lfac = static_cast<double>(L) / (2.0 * L + 1.0);
        double L1fac = (L + 1.0) / (2.0 * L + 1.0);
        for (int i = 0; i < n; ++i) {
            beta0_scalar[L][0][i] = s1_1T * Lfac * beta_arr[L][0][0][i] + L1fac * s2_1T * beta_arr[L][1][0][i];
            beta0_scalar[L][1][i] = s1_2T * Lfac * beta_arr[L][0][1][i] + L1fac * s2_2T * beta_arr[L][1][1][i];
        }

        // Construct betaSO
        double so_fac = 2.0 / (2.0 * L + 1.0);
        for (int i = 0; i < n; ++i) {
            beta_SO[L][i] = so_fac * (s2_1 * beta_arr[L][1][0][i] - s1_1 * beta_arr[L][0][0][i]
                                     + s2_2 * beta_arr[L][1][1][i] - s1_2 * beta_arr[L][0][1][i]);
        }

        // Compute integrals for normalization AND DD matrix elements
        // (both use RAW unnormalized beta0/betaSO, matching Fortran order)
        double sum1t = 0.0, sum2t = 0.0, sumSOt = 0.0;
        double rsum1 = 0.0, rsum2 = 0.0, rsumSO = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            double wT1 = x2 * beta_arr[L][1][0][i] + x1 * static_cast<double>(isign1) * beta_arr[L][0][0][i];
            double wT2 = x2 * beta_arr[L][1][1][i] + x1 * static_cast<double>(isign2) * beta_arr[L][0][1][i];
            double wT  = x2 * wave_arr[lw_ref][1][i] + x1 * static_cast<double>(isignw) * wave_arr[lw_ref][0][i];

            rsum1  += beta0_scalar[L][0][i] * wT1 * dr[i];
            rsum2  += beta0_scalar[L][1][i] * wT2 * dr[i];
            rsumSO += beta_SO[L][i]          * wT  * dr[i];

            sum1t  += beta0_scalar[L][0][i] * beta0_scalar[L][0][i] * dr[i];
            sum2t  += beta0_scalar[L][1][i] * beta0_scalar[L][1][i] * dr[i];
            sumSOt += beta_SO[L][i]          * beta_SO[L][i]          * dr[i];
        }

        // Normalize
        if (sum1t > 1e-30) {
            double norm = 1.0 / std::sqrt(sum1t);
            for (int i = 0; i < n; ++i) beta0_scalar[L][0][i] *= norm;
        }
        if (sum2t > 1e-30) {
            double norm = 1.0 / std::sqrt(sum2t);
            for (int i = 0; i < n; ++i) beta0_scalar[L][1][i] *= norm;
        }
        if (sumSOt > 1e-30) {
            double norm = 1.0 / std::sqrt(sumSOt);
            for (int i = 0; i < n; ++i) beta_SO[L][i] *= norm;
        }

        // DD matrix elements (using rsum computed from RAW beta0/betaSO)
        if (std::abs(rsum1) > 1e-30)  DD_0[L][0] = sum1t / rsum1;
        if (std::abs(rsum2) > 1e-30)  DD_0[L][1] = sum2t / rsum2;
        if (std::abs(rsumSO) > 1e-30) DD_SO[L]   = sumSOt / rsumSO;
    }

    // ------------------------------------------------------------
    // 4. Assemble output projectors
    // ------------------------------------------------------------
    int nbeta_new = 2 + 3 * lmax;   // 2xL=0 + 2xscalar(L) + SO(L) for each L=1..lmax

    std::vector<std::vector<double>> new_beta(nbeta_new, std::vector<double>(n, 0.0));
    std::vector<int> new_lll(nbeta_new, 0);
    array2d<double> new_D;
    new_D.rows = nbeta_new;
    new_D.cols = nbeta_new;
    new_D.data.assign(static_cast<std::size_t>(nbeta_new * nbeta_new), 0.0);

    // L=0: pass through
    int p_idx = 0;
    for (int m = 0; m < 2 && m < ind[0][1]; ++m) {
        int ll_orig = ll_map[0][1][m];   // L=0 only has J=1 (J=L+1/2)
        if (ll_orig >= 0) {
            for (int i = 0; i < n; ++i) {
                new_beta[p_idx][i] = beta_tmp[ll_orig][i];
            }
            new_D[p_idx, p_idx] = dd_diag(ll_orig);
        }
        new_lll[p_idx] = 0;
        p_idx++;
    }

    // L=1..lmax: beta0 scalar (m=0, m=1) then betaSO
    for (int L = 1; L <= lmax; ++L) {
        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta0_scalar[L][0][i];
        }
        new_D[p_idx, p_idx] = DD_0[L][0];
        new_lll[p_idx] = L;
        p_idx++;

        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta0_scalar[L][1][i];
        }
        new_D[p_idx, p_idx] = DD_0[L][1];
        new_lll[p_idx] = L;
        p_idx++;
    }

    for (int L = 1; L <= lmax; ++L) {
        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta_SO[L][i];
        }
        new_D[p_idx, p_idx] = DD_SO[L];
        new_lll[p_idx] = L * 10;    // 10, 20, 30 for L=1,2,3
        p_idx++;
    }

    // ------------------------------------------------------------
    // 5. Assemble output wavefunctions
    // ------------------------------------------------------------
    std::vector<int> l0_wfc_idx;
    for (int ll = 0; ll < nwave; ++ll) {
        if (map_Lw[ll] == 0) l0_wfc_idx.push_back(ll);
    }
    int nwfc_new = static_cast<int>(l0_wfc_idx.size());

    std::vector<std::vector<double>> new_chi;
    std::vector<int> new_lchi;
    std::vector<double> new_oc;

    for (int jj : l0_wfc_idx) {
        new_chi.push_back(wave_tmp[jj]);
        new_lchi.push_back(0);
        new_oc.push_back(p.pseudoWfc.occupation[jj]);
    }

    for (int L = 1; L <= lmax; ++L) {
        double sumw = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            sumw += wave_arr[L][0][i] * wave_arr[L][1][i] * dr[i];
        }
        int signw = (sumw > 0) ? 1 : -1;

        double x2 = std::sqrt(L + 1.0) / (std::sqrt(L + 1.0) + std::sqrt(L * 1.0));
        double x1 = 1.0 - x2;

        double sumt = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            double wT = x2 * wave_arr[L][1][i] + x1 * static_cast<double>(signw) * wave_arr[L][0][i];
            sumt += wT * wT * dr[i];
        }

        if (sumt > 1e-30) {
            std::vector<double> chi_L(n);
            double norm = 1.0 / std::sqrt(sumt);
            for (int i = 0; i < n; ++i) {
                chi_L[i] = (x2 * wave_arr[L][1][i] + x1 * static_cast<double>(signw) * wave_arr[L][0][i]) * norm;
            }
            new_chi.push_back(std::move(chi_L));
            new_lchi.push_back(L);
            new_oc.push_back(occ[L][0] + occ[L][1]);
            nwfc_new++;
        }
    }

    // ------------------------------------------------------------
    // 6. Update NCPP data
    // ------------------------------------------------------------
    // Truncate trailing zeros from new projectors
    for (auto& b : new_beta) {
        int cut = n;
        while (cut > 0 && std::abs(b[static_cast<std::size_t>(cut) - 1]) < 1e-30) --cut;
        b.resize(cut);
    }

    // Truncate trailing zeros from new wavefunctions
    for (auto& c : new_chi) {
        int cut = n;
        while (cut > 0 && std::abs(c[static_cast<std::size_t>(cut) - 1]) < 1e-30) --cut;
        c.resize(cut);
    }

    p.meta.has_so = false;
    p.meta.number_of_proj = nbeta_new;
    p.meta.number_of_wfc = nwfc_new;

    p.nonlocal.beta = std::move(new_beta);
    p.nonlocal.angular_momentum = std::move(new_lll);
    p.nonlocal.jjj.clear();
    p.nonlocal.B = std::move(new_D);
    p.nonlocal.cutoff_index.assign(nbeta_new, n);
    p.nonlocal.cutoff_radius.assign(nbeta_new,
        n > 0 ? p.mesh.r.back() : 0.0);

    p.pseudoWfc.chi = std::move(new_chi);
    p.pseudoWfc.angular_momentum = std::move(new_lchi);
    p.pseudoWfc.jchi.clear();
    p.pseudoWfc.occupation = std::move(new_oc);
}
