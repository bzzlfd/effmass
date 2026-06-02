module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.GKK;

import io.lattice;
import utils.array2d;
import utils.vector3d;
import std;


export {
    class GKK;
        struct GKKMetadata;
        enum class KVecsView : unsigned int;
            constexpr auto operator|(KVecsView a, KVecsView b) -> KVecsView;
            constexpr auto operator&(KVecsView a, KVecsView b) -> KVecsView;
            constexpr auto operator~(KVecsView a) -> KVecsView;
            constexpr auto hasView(KVecsView flags, KVecsView view) -> bool;
        struct KVecs;
}


// GKK file metadata structure
struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};


// Bitmask controlling which representations are computed and exposed in KVecs
enum class KVecsView : unsigned int {
    Cartesian = 1 << 0,  // kinetic, Kx, Ky, Kz
    Spherical = 1 << 1,  // q, theta, phi
    Integer   = 1 << 2,  // iG, jG, kG, kPoint, reciprocalLattice
};


constexpr auto operator|(KVecsView a, KVecsView b) -> KVecsView {
    return static_cast<KVecsView>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}
constexpr auto operator&(KVecsView a, KVecsView b) -> KVecsView {
    return static_cast<KVecsView>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
}
constexpr auto operator~(KVecsView a) -> KVecsView {
    return static_cast<KVecsView>(~static_cast<unsigned int>(a));
}
constexpr auto hasView(KVecsView flags, KVecsView view) -> bool {
    return (static_cast<unsigned int>(flags) & static_cast<unsigned int>(view)) != 0;
}


// k-point G-vector data view - non-owning spans to contiguous memory
struct KVecs {
    // Cartesian: data represents -K = -(G+k), negated at load time vs PWmat file convention
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, -(G_x+k_x), -(G_y+k_y), -(G_z+k_z)

    // Spherical representation of -K
    std::span<const double> q, theta, phi;        // |K|, polar angle [0,π], azimuthal angle [-π,π]

    // Integer indices of G vector in reciprocal lattice basis
    std::span<const int>    iG, jG, kG;           // G = iG*b1 + jG*b2 + kG*b3
    // Per-k-point metadata (valid whenever Integer view is enabled)
    vector3d<double>                    kPoint{};            // fractional coordinate of current k-point
    array2d<double, 3, 3>               reciprocalLattice{}; // reciprocal lattice vectors: row n = b_n
};


// GKK class - abstraction for OUT.GKK file
class GKK {
public:
    GKKMetadata meta;                             // metadata

    explicit GKK(const std::string& filename);  // open file and read metadata
    ~GKK();
    GKK(const GKK&) = delete;                   // disable copy
    auto operator=(const GKK&) -> GKK& = delete;
    GKK(GKK&& other) noexcept;                  // enable move
    auto operator=(GKK&& other) noexcept -> GKK&;

    // Control which representations are computed and exposed in currentData()
    auto setDataView(KVecsView view) -> void;
    auto currentView() const -> KVecsView { return desired_views_; }

    auto loadKPoint(int ikpt) -> const KVecs&;    // load k-point data (with cache)
    auto current_ikpt() const -> int { return current_ikpt_; }     // current k-point index
    auto currentData() const -> const KVecs& { return current_data_; }  // current data (filtered by currentView)
    auto inferCurrent_k() const -> vector3d<double>;  // infer k fractional coord from G-k data
    auto validateKineticConsistency() const -> void;        // validate kinetic = 0.5*|K|^2 for loaded k-point

private:
    auto readRecordLength() -> int;                     // read record length marker
    auto checkRecordLength(int expected) -> void;       // verify record length marker
    auto readRecord(void* dst, const char* context) -> void; // read full record
    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void; // read nbytes, skip rest

    auto readMetadata() -> void;                        // read file metadata
    auto readNgtotnod(int record_len) -> void;          // read ngtotnod array
    auto skipRecord() -> void;                            // skip one Fortran record by reading its length markers
    auto computeOffsets() -> void;                        // compute file offset per k-point
    auto seekToKPoint(int ikpt) -> void;                // seek to k-point data

    auto updateDataSpans(std::size_t ng) -> void;       // update current_data_ spans according to desired_views_
    auto computeSpherical(std::size_t ng) -> void;      // compute q, theta, phi from Kx, Ky, Kz
    auto computeIntegerIndices(std::size_t ng) -> void; // compute iG, jG, kG from Kx, Ky, Kz and inferred k

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
    std::vector<std::vector<int>> ngtotnod_;    // G-vector count per k-point per node
    std::vector<long> kpt_data_offsets_;        // file offset per k-point
    int current_ikpt_ = -1;                      // currently cached k-point

    KVecsView desired_views_ = KVecsView::Cartesian; // desired representations
    KVecsView ready_views_ = KVecsView{};  // which representations have been computed for the cached k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::size_t max_ng_ = 0;  // maximum possible G-vectors per k-point

    // buffers: working arrays (contiguous)
    std::vector<double> kinetic_, Kx_, Ky_, Kz_;
    std::vector<double> q_, theta_, phi_;
    std::vector<int>    iG_, jG_, kG_;
    KVecs current_data_;  // data view
};


GKK::GKK(const std::string& filename)
    : filename_(filename)
    , fp_(nullptr)
    , current_ikpt_(-1)
{
    fp_ = std::fopen(filename.c_str(), "rb");
    if (!fp_) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    readMetadata();

    // preallocate working buffers (maximum possible size) for Cartesian only
    max_ng_ = static_cast<std::size_t>(meta.mg_nx) * meta.nnodes;
    kinetic_.resize(max_ng_);
    Kx_.resize(max_ng_);
    Ky_.resize(max_ng_);
    Kz_.resize(max_ng_);

    // initialize empty data view
    current_data_ = {};

    // compute file offset for each k-point
    computeOffsets();
}


GKK::~GKK() {
    if (fp_) {
        std::fclose(fp_);
    }
}


GKK::GKK(GKK&& other) noexcept
    : filename_(std::move(other.filename_))
    , fp_(other.fp_)
    , meta(std::move(other.meta))
    , ngtotnod_(std::move(other.ngtotnod_))
    , kpt_data_offsets_(std::move(other.kpt_data_offsets_))
    , current_ikpt_(other.current_ikpt_)
    , desired_views_(other.desired_views_)
    , ready_views_(other.ready_views_)
    , max_ng_(other.max_ng_)
    , kinetic_(std::move(other.kinetic_))
    , Kx_(std::move(other.Kx_))
    , Ky_(std::move(other.Ky_))
    , Kz_(std::move(other.Kz_))
    , q_(std::move(other.q_))
    , theta_(std::move(other.theta_))
    , phi_(std::move(other.phi_))
    , iG_(std::move(other.iG_))
    , jG_(std::move(other.jG_))
    , kG_(std::move(other.kG_))
    , current_data_(other.current_data_)
{
    other.fp_ = nullptr;
    other.current_ikpt_ = -1;
    // update current_data_ spans to point to our own buffers
    if (!current_data_.kinetic.empty()) {
        const auto ng = current_data_.kinetic.size();
        current_data_.kinetic = std::span<const double>(kinetic_.data(), ng);
        current_data_.Kx      = std::span<const double>(Kx_.data(), ng);
        current_data_.Ky      = std::span<const double>(Ky_.data(), ng);
        current_data_.Kz      = std::span<const double>(Kz_.data(), ng);
        if (!current_data_.q.empty()) {
            current_data_.q     = std::span<const double>(q_.data(), ng);
            current_data_.theta = std::span<const double>(theta_.data(), ng);
            current_data_.phi   = std::span<const double>(phi_.data(), ng);
        }
        if (!current_data_.iG.empty()) {
            current_data_.iG    = std::span<const int>(iG_.data(), ng);
            current_data_.jG    = std::span<const int>(jG_.data(), ng);
            current_data_.kG    = std::span<const int>(kG_.data(), ng);
        }
    }
}


auto GKK::operator=(GKK&& other) noexcept -> GKK& {
    if (this != &other) {
        if (fp_) std::fclose(fp_);

        filename_ = std::move(other.filename_);
        fp_ = other.fp_;
        meta = std::move(other.meta);
        ngtotnod_ = std::move(other.ngtotnod_);
        kpt_data_offsets_ = std::move(other.kpt_data_offsets_);
        current_ikpt_ = other.current_ikpt_;
        desired_views_ = other.desired_views_;
        ready_views_ = other.ready_views_;
        max_ng_ = other.max_ng_;
        kinetic_ = std::move(other.kinetic_);
        Kx_ = std::move(other.Kx_);
        Ky_ = std::move(other.Ky_);
        Kz_ = std::move(other.Kz_);
        q_ = std::move(other.q_);
        theta_ = std::move(other.theta_);
        phi_ = std::move(other.phi_);
        iG_ = std::move(other.iG_);
        jG_ = std::move(other.jG_);
        kG_ = std::move(other.kG_);
        current_data_ = other.current_data_;

        other.fp_ = nullptr;
        other.current_ikpt_ = -1;

        // update spans
        if (!current_data_.kinetic.empty()) {
            const auto ng = current_data_.kinetic.size();
            current_data_.kinetic = std::span<const double>(kinetic_.data(), ng);
            current_data_.Kx      = std::span<const double>(Kx_.data(), ng);
            current_data_.Ky      = std::span<const double>(Ky_.data(), ng);
            current_data_.Kz      = std::span<const double>(Kz_.data(), ng);
            if (!current_data_.q.empty()) {
                current_data_.q     = std::span<const double>(q_.data(), ng);
                current_data_.theta = std::span<const double>(theta_.data(), ng);
                current_data_.phi   = std::span<const double>(phi_.data(), ng);
            }
            if (!current_data_.iG.empty()) {
                current_data_.iG    = std::span<const int>(iG_.data(), ng);
                current_data_.jG    = std::span<const int>(jG_.data(), ng);
                current_data_.kG    = std::span<const int>(kG_.data(), ng);
            }
        }
    }
    return *this;
}


auto GKK::readRecordLength() -> int {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    return length;
}


auto GKK::checkRecordLength(int expected_length) -> void {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    if (length != expected_length) {
        throw std::runtime_error("Record length mismatch");
    }
}


auto GKK::readRecord(void* dst, const char* context) -> void {
    int len = readRecordLength();
    if (std::fread(dst, 1, len, fp_) != static_cast<std::size_t>(len)) {
        throw std::runtime_error(std::string(context) + ": read failed");
    }
    checkRecordLength(len);
}

auto GKK::readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
    int len = readRecordLength();
    if (static_cast<int>(nbytes) > len) {
        throw std::runtime_error(std::string(context) + ": nbytes exceeds record length");
    }
    if (nbytes > 0) {
        if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
            throw std::runtime_error(std::string(context) + ": read failed");
        }
    }
    std::size_t remaining = static_cast<std::size_t>(len) - nbytes;
    if (remaining > 0) {
        if (std::fseek(fp_, static_cast<long>(remaining), SEEK_CUR) != 0) {
            throw std::runtime_error(std::string(context) + ": seek failed");
        }
    }
    checkRecordLength(len);
}


auto GKK::readMetadata() -> void {
    // Record 1: n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda
    int header[8];
    readRecord(header, "header");
    meta.n1 = header[0];
    meta.n2 = header[1];
    meta.n3 = header[2];
    meta.mg_nx = header[3];
    meta.nnodes = header[4];
    meta.nkpt = header[5];
    meta.is_SO = header[6];
    meta.islda = header[7];

    // handle spin-orbit coupling
    if (meta.is_SO == 1) {
        meta.mg_nx /= 2;
    }

    // Record 2: Ecut
    readRecord(&meta.Ecut, "Ecut");

    // Record 3: AL(3,3) - note Fortran is column-major
    double al_flat[9];
    readRecord(al_flat, "AL");
    meta.lattice.setLattice(al_flat, LengthUnit::Angstrom);

    // Record 4: nnodes, ngtotnod
    int len = readRecordLength();
    readNgtotnod(len);
    checkRecordLength(len);
}


auto GKK::readNgtotnod(int record_len) -> void {
    int nnodes_check;
    if (std::fread(&nnodes_check, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read nnodes");
    }
    if (nnodes_check != meta.nnodes) {
        throw std::runtime_error("nnodes mismatch");
    }

    // read G-vector count per k-point per node
    ngtotnod_.resize(meta.nkpt, std::vector<int>(meta.nnodes));
    meta.ng_tot_per_kpt.resize(meta.nkpt, 0);

    for (int ikpt = 0; ikpt < meta.nkpt; ++ikpt) {
        int ng_total = 0;
        for (int n = 0; n < meta.nnodes; ++n) {
            int ng;
            if (std::fread(&ng, sizeof(int), 1, fp_) != 1) {
                throw std::runtime_error("Failed to read ngtotnod");
            }
            if (meta.is_SO == 1) {
                ng /= 2;
            }
            ngtotnod_[ikpt][n] = ng;
            ng_total += ng;
        }
        meta.ng_tot_per_kpt[ikpt] = ng_total;
    }
}


auto GKK::skipRecord() -> void {
    int len = readRecordLength();
    if (std::fseek(fp_, len, SEEK_CUR) != 0) {
        throw std::runtime_error("Failed to skip record data");
    }
    checkRecordLength(len);
}


auto GKK::computeOffsets() -> void {
    // record the starting file offset for each k-point's data
    kpt_data_offsets_.resize(meta.nkpt);

    for (int ikpt = 0; ikpt < meta.nkpt; ++ikpt) {
        kpt_data_offsets_[ikpt] = std::ftell(fp_);

        // skip all data for this k-point by walking through records
        // This correctly handles compiler-dependent record padding (alignment),
        // because we trust the length markers in the file rather than computing
        // the offset from ng.
        for (int n = 0; n < meta.nnodes; ++n) {
            // 4 arrays (gkk, gkk_x, gkk_y, gkk_z), each in its own record
            for (int i = 0; i < 4; ++i) {
                skipRecord();
            }
        }
    }
}


auto GKK::seekToKPoint(int ikpt) -> void {
    if (ikpt < 0 || ikpt >= meta.nkpt) {
        throw std::out_of_range("Invalid k-point index");
    }

    if (std::fseek(fp_, kpt_data_offsets_[ikpt], SEEK_SET) != 0) {
        throw std::runtime_error("Failed to seek to k-point");
    }
}


auto GKK::updateDataSpans(std::size_t ng) -> void {
    // Cartesian is always loaded when loadKPoint succeeds
    current_data_.kinetic = std::span<const double>(kinetic_.data(), ng);
    current_data_.Kx      = std::span<const double>(Kx_.data(), ng);
    current_data_.Ky      = std::span<const double>(Ky_.data(), ng);
    current_data_.Kz      = std::span<const double>(Kz_.data(), ng);

    if (hasView(desired_views_, KVecsView::Spherical) && hasView(ready_views_, KVecsView::Spherical)) {
        current_data_.q     = std::span<const double>(q_.data(), ng);
        current_data_.theta = std::span<const double>(theta_.data(), ng);
        current_data_.phi   = std::span<const double>(phi_.data(), ng);
    } else {
        current_data_.q = current_data_.theta = current_data_.phi = {};
    }

    if (hasView(desired_views_, KVecsView::Integer) && hasView(ready_views_, KVecsView::Integer)) {
        current_data_.iG = std::span<const int>(iG_.data(), ng);
        current_data_.jG = std::span<const int>(jG_.data(), ng);
        current_data_.kG = std::span<const int>(kG_.data(), ng);
    } else {
        current_data_.iG = current_data_.jG = current_data_.kG = {};
    }
}


auto GKK::computeSpherical(std::size_t ng) -> void {
    // Handles kinetic=0 (|K|=0, i.e. Kx=Ky=Kz=0, the Gamma-point G=G-k vector):
    // q=0 → theta=0 via guard below, phi=0 via atan2(0,0)=0. No division-by-zero.
    for (std::size_t ig = 0; ig < ng; ++ig) {
        const double kx = Kx_[ig];
        const double ky = Ky_[ig];
        const double kz = Kz_[ig];
        const double q  = std::sqrt(kx * kx + ky * ky + kz * kz);
        q_[ig]     = q;
        theta_[ig] = (q > 0.0) ? std::acos(kz / q) : 0.0;
        phi_[ig]   = std::atan2(ky, kx);
    }
}


auto GKK::computeIntegerIndices(std::size_t ng) -> void {
    auto k_frac = inferCurrent_k();
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    auto A_mat = meta.lattice.A();
    for (std::size_t ig = 0; ig < ng; ++ig) {
        const double Kx = - Kx_[ig];
        const double Ky = - Ky_[ig];
        const double Kz = - Kz_[ig];

        // c[n] = A[n][0] * Kx + A[n][1] * Ky + A[n][2] * Kz
        // (iG,jG,kG) = round(c / (2π) + k_frac)
        // note: Data is read as -K = -(G+k), so A·(-K)/(2π) = -(iG + k_frac).
        //       iG = round(-A·(-K)/(2π) - k_frac) = round(iG + k_frac - k_frac) = iG.
        double cx = (A_mat[0][0] * Kx + A_mat[0][1] * Ky + A_mat[0][2] * Kz);
        double cy = (A_mat[1][0] * Kx + A_mat[1][1] * Ky + A_mat[1][2] * Kz);
        double cz = (A_mat[2][0] * Kx + A_mat[2][1] * Ky + A_mat[2][2] * Kz);

        iG_[ig] = static_cast<int>(std::lround(cx / TWO_PI - k_frac[0]));
        jG_[ig] = static_cast<int>(std::lround(cy / TWO_PI - k_frac[1]));
        kG_[ig] = static_cast<int>(std::lround(cz / TWO_PI - k_frac[2]));
    }
}


auto GKK::setDataView(KVecsView view) -> void {
    if (view == desired_views_) {
        return;
    }

    // If new view does not need Spherical, release buffers and clear cached state
    if (!hasView(view, KVecsView::Spherical)) {
        q_.clear();
        q_.shrink_to_fit();
        theta_.clear();
        theta_.shrink_to_fit();
        phi_.clear();
        phi_.shrink_to_fit();
        current_data_.q = current_data_.theta = current_data_.phi = {};
        ready_views_ = ready_views_ & ~KVecsView::Spherical;
    } else if (q_.empty()) {
        q_.resize(max_ng_);
        theta_.resize(max_ng_);
        phi_.resize(max_ng_);
    }

    // If new view does not need Integer, release buffers and clear cached state
    if (!hasView(view, KVecsView::Integer)) {
        iG_.clear();
        iG_.shrink_to_fit();
        jG_.clear();
        jG_.shrink_to_fit();
        kG_.clear();
        kG_.shrink_to_fit();
        current_data_.iG = current_data_.jG = current_data_.kG = {};
        ready_views_ = ready_views_ & ~KVecsView::Integer;
    } else if (iG_.empty()) {
        iG_.resize(max_ng_);
        jG_.resize(max_ng_);
        kG_.resize(max_ng_);
    }

    desired_views_ = view;
}


auto GKK::loadKPoint(int ikpt) -> const KVecs& {
    // exact cache hit: same k-point and all desired representations already computed
    if (ikpt == current_ikpt_ && ready_views_ == desired_views_) {
        return current_data_;
    }

    if (ikpt < 0 || ikpt >= meta.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }

    bool same_kpt = (ikpt == current_ikpt_);
    std::size_t total_pos = 0;

    if (!same_kpt) {
        seekToKPoint(ikpt);

        // read all nodes for this k-point and merge into contiguous buffers
        // Fortran records are mg_nx-sized; readRecord(dst, nbytes, ...) reads nbytes
        // and automatically skips the remaining (mg_nx - ng) elements.
        for (int inode = 0; inode < meta.nnodes; ++inode) {
            int ng = ngtotnod_[ikpt][inode];
            if (ng == 0) continue;

            auto nbytes = static_cast<std::size_t>(ng) * sizeof(double);
            readRecord(kinetic_.data() + total_pos, nbytes, "gkk");
            readRecord(Kx_.data() + total_pos,      nbytes, "gkk_x");
            readRecord(Ky_.data() + total_pos,      nbytes, "gkk_y");
            readRecord(Kz_.data() + total_pos,      nbytes, "gkk_z");

            total_pos += static_cast<std::size_t>(ng);
        }

        current_ikpt_ = ikpt;
        ready_views_  = KVecsView::Cartesian;  // at minimum Cartesian is now valid
    } else {
        // same k-point, just need to compute additional representations
        total_pos = current_data_.kinetic.size();
    }

    // compute additional representations as needed
    auto missing = static_cast<unsigned int>(desired_views_) & ~static_cast<unsigned int>(ready_views_);

    if (total_pos > 0) {
        updateDataSpans(total_pos);  // Cartesian spans needed by computeSpherical / computeIntegerIndices

        if ((missing & static_cast<unsigned int>(KVecsView::Spherical)) != 0) {
            computeSpherical(total_pos);
            ready_views_ = ready_views_ | KVecsView::Spherical;
        }
        if ((missing & static_cast<unsigned int>(KVecsView::Integer)) != 0) {
            computeIntegerIndices(total_pos);
            ready_views_ = ready_views_ | KVecsView::Integer;
        }
    }

    // populate per-k-point metadata for Integer view
    if (hasView(desired_views_, KVecsView::Integer)) {
        auto B = meta.lattice.B();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                current_data_.reciprocalLattice[i, j] = B[i][j];
        if (total_pos > 0) {
            current_data_.kPoint = inferCurrent_k();
        } else {
            current_data_.kPoint = {0.0, 0.0, 0.0};
        }
    } else {
        current_data_.kPoint = {};
        current_data_.reciprocalLattice = {};
    }

    // re-wire spans (now includes Spherical/Integer if just computed)
    updateDataSpans(total_pos);
    return current_data_;
}


auto GKK::inferCurrent_k() const -> vector3d<double> {
    if (current_ikpt_ < 0) {
        throw std::runtime_error("inferCurrent_k: no k-point loaded");
    }

    const auto& data = current_data_;
    if (data.Kx.empty()) {
        throw std::runtime_error("inferCurrent_k: current k-point has no G-vectors");
    }

    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    const std::size_t ng = data.Kx.size();

    vector3d<double> k_frac{};

    auto A_mat = meta.lattice.A();

    for (int dim = 0; dim < 3; ++dim) {
        double sum_cos = 0.0;
        double sum_sin = 0.0;

        for (std::size_t ig = 0; ig < ng; ++ig) {
            // Data is interpreted as -K = -(G+k) in Bohr^-1.
            const double Kx = - data.Kx[ig];
            const double Ky = - data.Ky[ig];
            const double Kz = - data.Kz[ig];

            // A·(K) / (2π) = A·(G+k) / (2π) = iG + k_frac  (no leading minus needed)
            const double c = (A_mat[dim][0] * Kx + A_mat[dim][1] * Ky + A_mat[dim][2] * Kz) / TWO_PI;

            // d = k_frac mod 1 (fractional part of iG + k_frac is just k_frac)
            const double d = c - std::floor(c);   // [0, 1)
            const double theta = d * TWO_PI;
            sum_cos += std::cos(theta);
            sum_sin += std::sin(theta);
        }

        double avg_theta = std::atan2(sum_sin, sum_cos);
        if (avg_theta < 0.0) {
            avg_theta += TWO_PI;
        }
        const double d_avg = avg_theta / TWO_PI;   // [0, 1)

        // wrap d_avg into (-0.5, 0.5]
        double k = d_avg;
        if (k > 0.5) k -= 1.0;
        k_frac[dim] = k;
    }

    return k_frac;
}


auto GKK::validateKineticConsistency() const -> void {
    const auto& data = current_data_;
    const auto ng = data.kinetic.size();
    constexpr double eps = 1e-10;

    for (std::size_t i = 0; i < ng; ++i) {
        const double kx = data.Kx[i];
        const double ky = data.Ky[i];
        const double kz = data.Kz[i];
        const double expected = 0.5 * (kx * kx + ky * ky + kz * kz);
        const double diff = data.kinetic[i] - expected;
        if (diff > eps * (1.0 + expected) || diff < -eps * (1.0 + expected)) {
            throw std::runtime_error(
                "GKK kinetic-energy mismatch at k-point " + std::to_string(current_ikpt_) +
                ", G-vector " + std::to_string(i) +
                ": kinetic=" + std::to_string(data.kinetic[i]) +
                " != 0.5*(Kx^2+Ky^2+Kz^2)=" + std::to_string(expected));
        }
    }
}
