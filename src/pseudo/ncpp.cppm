export module pseudo.ncpp;

import std;
import utils.matrix;
import pseudo.io.ncpp_upf;

export {
    class NCPP;
    enum class MeshType : int;
    struct NonlocalByL;
}


enum class MeshType {
    Uniform,
    Exponential,
    Unknown
};


struct NonlocalByL {
    std::vector<std::vector<double>> beta;
    DenseMatrix<double> dion;
};


class NCPP {
public:
    explicit NCPP(const NCPPUPF& upf);

    // Access raw UPF data
    auto header() const -> const NCPPUPFHeader& { return upf_.header(); }
    auto mesh() const -> const NCPPUPFMesh& { return upf_.mesh(); }
    auto localPotential() const -> std::span<const double> { return upf_.localPotential(); }
    auto nonlocal() const -> const NCPPUPFNonlocal& { return upf_.nonlocal(); }
    auto wavefunctions() const -> const NCPPUPFWavefunction& { return upf_.wavefunctions(); }
    auto rhoAtom() const -> std::span<const double> { return upf_.rhoAtom(); }
    auto upfData() const -> const NCPPUPF& { return upf_; }

    // Higher-level interfaces
    auto meshType() const -> MeshType;
    auto nonlocalByL(int l) const -> NonlocalByL;

private:
    NCPPUPF upf_;
};


NCPP::NCPP(const NCPPUPF& upf) : upf_(upf) {}


auto NCPP::meshType() const -> MeshType {
    const auto& r = upf_.mesh().r;
    const auto& rab = upf_.mesh().rab;
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


auto NCPP::nonlocalByL(int l) const -> NonlocalByL {
    const auto& nl = upf_.nonlocal();
    const int nb = upf_.header().number_of_proj;

    std::vector<int> indices;
    for (int i = 0; i < nb; ++i) {
        if (nl.lll[i] == l) {
            indices.push_back(i);
        }
    }

    NonlocalByL result;
    const int nbl = static_cast<int>(indices.size());
    result.beta.resize(nbl);
    for (int i = 0; i < nbl; ++i) {
        result.beta[i] = nl.beta[indices[i]];
    }

    result.dion.rows = nbl;
    result.dion.cols = nbl;
    result.dion.data.resize(static_cast<std::size_t>(nbl * nbl));
    for (int i = 0; i < nbl; ++i) {
        for (int j = 0; j < nbl; ++j) {
            result.dion[i, j] = nl.dion[indices[i], indices[j]];
        }
    }

    return result;
}
