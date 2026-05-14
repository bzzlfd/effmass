module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.GKK;

import io.lattice;
import std;


// GKK file metadata structure
export struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};


// Bitmask controlling which representations are computed and exposed in KVecs
export enum class KVecsView : unsigned int {
    Cartesian = 1 << 0,  // kinetic, Kx, Ky, Kz
    Spherical = 1 << 1,  // r, theta, phi
    Integer   = 1 << 2,  // iG, jG, kG, kPoint, reciprocalLattice
};


export constexpr auto operator|(KVecsView a, KVecsView b) -> KVecsView {
    return static_cast<KVecsView>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}
export constexpr auto operator&(KVecsView a, KVecsView b) -> KVecsView {
    return static_cast<KVecsView>(static_cast<unsigned int>(a) & static_cast<unsigned int>(b));
}
export constexpr auto operator~(KVecsView a) -> KVecsView {
    return static_cast<KVecsView>(~static_cast<unsigned int>(a));
}
export constexpr auto hasView(KVecsView flags, KVecsView view) -> bool {
    return (static_cast<unsigned int>(flags) & static_cast<unsigned int>(view)) != 0;
}


// k-point G-vector data view - non-owning spans to contiguous memory
export struct KVecs {
    // Cartesian representation: K = G - k
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, G_x-k_x, G_y-k_y, G_z-k_z
    
    // Spherical representation of K = G - k
    std::span<const double> r, theta, phi;        // |K|, polar angle [0,π], azimuthal angle [-π,π]

    // Integer indices of G vector in reciprocal lattice basis
    std::span<const int>    iG, jG, kG;           // G = iG*b1 + jG*b2 + kG*b3
    // Per-k-point metadata (valid whenever Integer view is enabled)
    std::array<double, 3>               kPoint{};            // fractional coordinate of current k-point
    std::array<std::array<double, 3>, 3> reciprocalLattice{}; // reciprocal lattice vectors: row n = b_n
};


// GKK class - abstraction for OUT.GKK file
export class GKK {
public:
    explicit GKK(const std::string& filename);  // open file and read metadata
    ~GKK();
    GKK(const GKK&) = delete;                   // disable copy
    auto operator=(const GKK&) -> GKK& = delete;
    GKK(GKK&& other) noexcept;                  // enable move
    auto operator=(GKK&& other) noexcept -> GKK&;

    auto metadata() const -> const GKKMetadata& { return meta_; }  // get metadata

    // Control which representations are computed and exposed in currentData()
    auto setDataView(KVecsView view) -> void;
    auto currentView() const -> KVecsView { return desired_views_; }

    auto loadKPoint(int ikpt) -> const KVecs&;    // load k-point data (with cache)
    auto current_ikpt() const -> int { return current_ikpt_; }     // current k-point index
    auto currentData() const -> const KVecs& { return current_data_; }  // current data (filtered by currentView)
    auto inferCurrent_k() const -> std::array<double, 3>;  // infer k fractional coord from G-k data

private:
    auto readRecordLength() -> int;                     // read record length marker
    auto checkRecordLength(int expected) -> void;       // verify record length marker
    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void; // read full record data

    auto readMetadata() -> void;                        // read file metadata
    auto readNgtotnod(int record_len) -> void;          // read ngtotnod array
    auto skipRecord() -> void;                            // skip one Fortran record by reading its length markers
    auto computeOffsets() -> void;                        // compute file offset per k-point
    auto seekToKPoint(int ikpt) -> void;                // seek to k-point data

    auto updateDataSpans(std::size_t ng) -> void;       // update current_data_ spans according to desired_views_
    auto computeSpherical(std::size_t ng) -> void;      // compute r, theta, phi from Kx, Ky, Kz
    auto computeIntegerIndices(std::size_t ng) -> void; // compute iG, jG, kG from Kx, Ky, Kz and inferred k

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
    GKKMetadata meta_;                          // metadata
    std::vector<std::vector<int>> ngtotnod_;    // G-vector count per k-point per node
    std::vector<long> kpt_data_offsets_;        // file offset per k-point
    int current_ikpt_ = -1;                      // currently cached k-point

    KVecsView desired_views_ = KVecsView::Cartesian; // desired representations
    KVecsView ready_views_ = KVecsView{};  // which representations have been computed for the cached k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::size_t max_ng_ = 0;  // maximum possible G-vectors per k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::vector<double> kinetic_buf_, Kx_buf_, Ky_buf_, Kz_buf_;
    std::vector<double> r_buf_, theta_buf_, phi_buf_;
    std::vector<int>    iG_buf_, jG_buf_, kG_buf_;
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
    max_ng_ = static_cast<std::size_t>(meta_.mg_nx) * meta_.nnodes;
    kinetic_buf_.resize(max_ng_);
    Kx_buf_.resize(max_ng_);
    Ky_buf_.resize(max_ng_);
    Kz_buf_.resize(max_ng_);

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
    , meta_(std::move(other.meta_))
    , ngtotnod_(std::move(other.ngtotnod_))
    , kpt_data_offsets_(std::move(other.kpt_data_offsets_))
    , current_ikpt_(other.current_ikpt_)
    , desired_views_(other.desired_views_)
    , ready_views_(other.ready_views_)
    , max_ng_(other.max_ng_)
    , kinetic_buf_(std::move(other.kinetic_buf_))
    , Kx_buf_(std::move(other.Kx_buf_))
    , Ky_buf_(std::move(other.Ky_buf_))
    , Kz_buf_(std::move(other.Kz_buf_))
    , r_buf_(std::move(other.r_buf_))
    , theta_buf_(std::move(other.theta_buf_))
    , phi_buf_(std::move(other.phi_buf_))
    , iG_buf_(std::move(other.iG_buf_))
    , jG_buf_(std::move(other.jG_buf_))
    , kG_buf_(std::move(other.kG_buf_))
    , current_data_(other.current_data_)
{
    other.fp_ = nullptr;
    other.current_ikpt_ = -1;
    // update current_data_ spans to point to our own buffers
    if (!current_data_.kinetic.empty()) {
        const auto ng = current_data_.kinetic.size();
        current_data_.kinetic = std::span<const double>(kinetic_buf_.data(), ng);
        current_data_.Kx      = std::span<const double>(Kx_buf_.data(), ng);
        current_data_.Ky      = std::span<const double>(Ky_buf_.data(), ng);
        current_data_.Kz      = std::span<const double>(Kz_buf_.data(), ng);
        if (!current_data_.r.empty()) {
            current_data_.r     = std::span<const double>(r_buf_.data(), ng);
            current_data_.theta = std::span<const double>(theta_buf_.data(), ng);
            current_data_.phi   = std::span<const double>(phi_buf_.data(), ng);
        }
        if (!current_data_.iG.empty()) {
            current_data_.iG    = std::span<const int>(iG_buf_.data(), ng);
            current_data_.jG    = std::span<const int>(jG_buf_.data(), ng);
            current_data_.kG    = std::span<const int>(kG_buf_.data(), ng);
        }
    }
}


auto GKK::operator=(GKK&& other) noexcept -> GKK& {
    if (this != &other) {
        if (fp_) std::fclose(fp_);

        filename_ = std::move(other.filename_);
        fp_ = other.fp_;
        meta_ = std::move(other.meta_);
        ngtotnod_ = std::move(other.ngtotnod_);
        kpt_data_offsets_ = std::move(other.kpt_data_offsets_);
        current_ikpt_ = other.current_ikpt_;
        desired_views_ = other.desired_views_;
        ready_views_ = other.ready_views_;
        max_ng_ = other.max_ng_;
        kinetic_buf_ = std::move(other.kinetic_buf_);
        Kx_buf_ = std::move(other.Kx_buf_);
        Ky_buf_ = std::move(other.Ky_buf_);
        Kz_buf_ = std::move(other.Kz_buf_);
        r_buf_ = std::move(other.r_buf_);
        theta_buf_ = std::move(other.theta_buf_);
        phi_buf_ = std::move(other.phi_buf_);
        iG_buf_ = std::move(other.iG_buf_);
        jG_buf_ = std::move(other.jG_buf_);
        kG_buf_ = std::move(other.kG_buf_);
        current_data_ = other.current_data_;

        other.fp_ = nullptr;
        other.current_ikpt_ = -1;

        // update spans
        if (!current_data_.kinetic.empty()) {
            const auto ng = current_data_.kinetic.size();
            current_data_.kinetic = std::span<const double>(kinetic_buf_.data(), ng);
            current_data_.Kx      = std::span<const double>(Kx_buf_.data(), ng);
            current_data_.Ky      = std::span<const double>(Ky_buf_.data(), ng);
            current_data_.Kz      = std::span<const double>(Kz_buf_.data(), ng);
            if (!current_data_.r.empty()) {
                current_data_.r     = std::span<const double>(r_buf_.data(), ng);
                current_data_.theta = std::span<const double>(theta_buf_.data(), ng);
                current_data_.phi   = std::span<const double>(phi_buf_.data(), ng);
            }
            if (!current_data_.iG.empty()) {
                current_data_.iG    = std::span<const int>(iG_buf_.data(), ng);
                current_data_.jG    = std::span<const int>(jG_buf_.data(), ng);
                current_data_.kG    = std::span<const int>(kG_buf_.data(), ng);
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


auto GKK::readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
    int len = readRecordLength();
    if (len != static_cast<int>(nbytes)) {
        throw std::runtime_error(std::string(context) + ": record size mismatch");
    }
    if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
        throw std::runtime_error(std::string(context) + ": read failed");
    }
    checkRecordLength(len);
}


auto GKK::readMetadata() -> void {
    // Record 1: n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda
    int header[8];
    readRecord(header, sizeof(header), "header");
    meta_.n1 = header[0];
    meta_.n2 = header[1];
    meta_.n3 = header[2];
    meta_.mg_nx = header[3];
    meta_.nnodes = header[4];
    meta_.nkpt = header[5];
    meta_.is_SO = header[6];
    meta_.islda = header[7];

    // handle spin-orbit coupling
    if (meta_.is_SO == 1) {
        meta_.mg_nx /= 2;
    }

    // Record 2: Ecut
    readRecord(&meta_.Ecut, sizeof(double), "Ecut");

    // Record 3: AL(3,3) - note Fortran is column-major
    double al_flat[9];
    readRecord(al_flat, sizeof(al_flat), "AL");
    meta_.lattice = Lattice(al_flat, LengthUnit::Angstrom);

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
    if (nnodes_check != meta_.nnodes) {
        throw std::runtime_error("nnodes mismatch");
    }

    // read G-vector count per k-point per node
    ngtotnod_.resize(meta_.nkpt, std::vector<int>(meta_.nnodes));
    meta_.ng_tot_per_kpt.resize(meta_.nkpt, 0);

    for (int ikpt = 0; ikpt < meta_.nkpt; ++ikpt) {
        int ng_total = 0;
        for (int n = 0; n < meta_.nnodes; ++n) {
            int ng;
            if (std::fread(&ng, sizeof(int), 1, fp_) != 1) {
                throw std::runtime_error("Failed to read ngtotnod");
            }
            if (meta_.is_SO == 1) {
                ng /= 2;
            }
            ngtotnod_[ikpt][n] = ng;
            ng_total += ng;
        }
        meta_.ng_tot_per_kpt[ikpt] = ng_total;
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
    kpt_data_offsets_.resize(meta_.nkpt);

    for (int ikpt = 0; ikpt < meta_.nkpt; ++ikpt) {
        kpt_data_offsets_[ikpt] = std::ftell(fp_);

        // skip all data for this k-point by walking through records
        // This correctly handles compiler-dependent record padding (alignment),
        // because we trust the length markers in the file rather than computing
        // the offset from ng.
        for (int n = 0; n < meta_.nnodes; ++n) {
            // 4 arrays (gkk, gkk_x, gkk_y, gkk_z), each in its own record
            for (int i = 0; i < 4; ++i) {
                skipRecord();
            }
        }
    }
}


auto GKK::seekToKPoint(int ikpt) -> void {
    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index");
    }

    if (std::fseek(fp_, kpt_data_offsets_[ikpt], SEEK_SET) != 0) {
        throw std::runtime_error("Failed to seek to k-point");
    }
}


auto GKK::updateDataSpans(std::size_t ng) -> void {
    // Cartesian is always loaded when loadKPoint succeeds
    current_data_.kinetic = std::span<const double>(kinetic_buf_.data(), ng);
    current_data_.Kx      = std::span<const double>(Kx_buf_.data(), ng);
    current_data_.Ky      = std::span<const double>(Ky_buf_.data(), ng);
    current_data_.Kz      = std::span<const double>(Kz_buf_.data(), ng);

    if (hasView(desired_views_, KVecsView::Spherical) && hasView(ready_views_, KVecsView::Spherical)) {
        current_data_.r     = std::span<const double>(r_buf_.data(), ng);
        current_data_.theta = std::span<const double>(theta_buf_.data(), ng);
        current_data_.phi   = std::span<const double>(phi_buf_.data(), ng);
    } else {
        current_data_.r = current_data_.theta = current_data_.phi = {};
    }

    if (hasView(desired_views_, KVecsView::Integer) && hasView(ready_views_, KVecsView::Integer)) {
        current_data_.iG = std::span<const int>(iG_buf_.data(), ng);
        current_data_.jG = std::span<const int>(jG_buf_.data(), ng);
        current_data_.kG = std::span<const int>(kG_buf_.data(), ng);
    } else {
        current_data_.iG = current_data_.jG = current_data_.kG = {};
    }
}


auto GKK::computeSpherical(std::size_t ng) -> void {
    for (std::size_t ig = 0; ig < ng; ++ig) {
        const double kx = Kx_buf_[ig];
        const double ky = Ky_buf_[ig];
        const double kz = Kz_buf_[ig];
        const double r  = std::sqrt(kx * kx + ky * ky + kz * kz);
        r_buf_[ig]     = r;
        theta_buf_[ig] = (r > 0.0) ? std::acos(kz / r) : 0.0;
        phi_buf_[ig]   = std::atan2(ky, kx);
    }
}


auto GKK::computeIntegerIndices(std::size_t ng) -> void {
    auto k_frac = inferCurrent_k();
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    for (std::size_t ig = 0; ig < ng; ++ig) {
        const double kx = Kx_buf_[ig];
        const double ky = Ky_buf_[ig];
        const double kz = Kz_buf_[ig];

        // c[n] = A[n][0] * Kx + A[n][1] * Ky + A[n][2] * Kz
        // (iG,jG,kG) = round(c / (2π) + k_frac)
        double cx = meta_.lattice.A()[0][0] * kx + meta_.lattice.A()[0][1] * ky + meta_.lattice.A()[0][2] * kz;
        double cy = meta_.lattice.A()[1][0] * kx + meta_.lattice.A()[1][1] * ky + meta_.lattice.A()[1][2] * kz;
        double cz = meta_.lattice.A()[2][0] * kx + meta_.lattice.A()[2][1] * ky + meta_.lattice.A()[2][2] * kz;

        iG_buf_[ig] = static_cast<int>(std::lround(cx / TWO_PI + k_frac[0]));
        jG_buf_[ig] = static_cast<int>(std::lround(cy / TWO_PI + k_frac[1]));
        kG_buf_[ig] = static_cast<int>(std::lround(cz / TWO_PI + k_frac[2]));
    }
}


auto GKK::setDataView(KVecsView view) -> void {
    if (view == desired_views_) {
        return;
    }

    // If new view does not need Spherical, release buffers and clear cached state
    if (!hasView(view, KVecsView::Spherical)) {
        r_buf_.clear();
        r_buf_.shrink_to_fit();
        theta_buf_.clear();
        theta_buf_.shrink_to_fit();
        phi_buf_.clear();
        phi_buf_.shrink_to_fit();
        current_data_.r = current_data_.theta = current_data_.phi = {};
        ready_views_ = ready_views_ & ~KVecsView::Spherical;
    } else if (r_buf_.empty()) {
        r_buf_.resize(max_ng_);
        theta_buf_.resize(max_ng_);
        phi_buf_.resize(max_ng_);
    }

    // If new view does not need Integer, release buffers and clear cached state
    if (!hasView(view, KVecsView::Integer)) {
        iG_buf_.clear();
        iG_buf_.shrink_to_fit();
        jG_buf_.clear();
        jG_buf_.shrink_to_fit();
        kG_buf_.clear();
        kG_buf_.shrink_to_fit();
        current_data_.iG = current_data_.jG = current_data_.kG = {};
        ready_views_ = ready_views_ & ~KVecsView::Integer;
    } else if (iG_buf_.empty()) {
        iG_buf_.resize(max_ng_);
        jG_buf_.resize(max_ng_);
        kG_buf_.resize(max_ng_);
    }

    desired_views_ = view;
}


auto GKK::loadKPoint(int ikpt) -> const KVecs& {
    // exact cache hit: same k-point and all desired representations already computed
    if (ikpt == current_ikpt_ && ready_views_ == desired_views_) {
        return current_data_;
    }

    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }

    bool same_kpt = (ikpt == current_ikpt_);
    std::size_t total_pos = 0;

    if (!same_kpt) {
        seekToKPoint(ikpt);

        // read all nodes for this k-point and merge into contiguous buffers
        for (int inode = 0; inode < meta_.nnodes; ++inode) {
            int ng = ngtotnod_[ikpt][inode];
            if (ng == 0) continue;

            readRecord(kinetic_buf_.data() + total_pos, ng * sizeof(double), "gkk");
            readRecord(Kx_buf_.data() + total_pos,     ng * sizeof(double), "gkk_x");
            readRecord(Ky_buf_.data() + total_pos,     ng * sizeof(double), "gkk_y");
            readRecord(Kz_buf_.data() + total_pos,     ng * sizeof(double), "gkk_z");

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
        current_data_.reciprocalLattice = meta_.lattice.B();
        if (total_pos > 0) {
            current_data_.kPoint = inferCurrent_k();
        } else {
            current_data_.kPoint = {0.0, 0.0, 0.0};
        }
    } else {
        current_data_.kPoint = {};
        current_data_.reciprocalLattice = {};
    }

    updateDataSpans(total_pos);
    return current_data_;
}


auto GKK::inferCurrent_k() const -> std::array<double, 3> {
    if (current_ikpt_ < 0) {
        throw std::runtime_error("inferCurrent_k: no k-point loaded");
    }

    const auto& data = current_data_;
    if (data.Kx.empty()) {
        throw std::runtime_error("inferCurrent_k: current k-point has no G-vectors");
    }

    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    const std::size_t ng = data.Kx.size();

    std::array<double, 3> k_frac = {0.0, 0.0, 0.0};

    for (int dim = 0; dim < 3; ++dim) {
        double sum_cos = 0.0;
        double sum_sin = 0.0;

        for (std::size_t ig = 0; ig < ng; ++ig) {
            // v = G - k (Cartesian, in Bohr^-1, consistent with AL in Bohr)
            const double Kx = data.Kx[ig];
            const double Ky = data.Ky[ig];
            const double Kz = data.Kz[ig];

            // c = a_dim · v / (2π), fractional coordinate in reciprocal basis
            const double c = (meta_.lattice.A()[dim][0] * Kx + meta_.lattice.A()[dim][1] * Ky + meta_.lattice.A()[dim][2] * Kz) / TWO_PI;

            // c = n - k_frac, thus c mod 1 = (-k_frac) mod 1
            // Use circular mean on [0,1) to robustly handle wrap-around at 0/1 boundary
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

        // k = wrap(-d_avg) into (-0.5, 0.5]
        double k = -d_avg;
        k = std::fmod(k, 1.0);
        if (k <= -0.5) k += 1.0;
        if (k > 0.5)   k -= 1.0;
        k_frac[dim] = k;
    }

    return k_frac;
}
