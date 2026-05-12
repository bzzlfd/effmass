module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.GKK;

import std;

// GKK file metadata structure
export struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    double AL[3][3];          // lattice vectors (Bohr)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};

// k-point G-vector data view - non-owning spans to contiguous memory
export struct KVecs {
    std::span<const double> kinetic, Kx, Ky, Kz;  // |G+k|²/2, G_x-k_x, G_y-k_y, G_z-k_z
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
    auto loadKPoint(int ikpt) -> const KVecs&;    // load k-point data (with cache)
    auto current_ikpt() const -> int { return current_ikpt_; }     // current k-point index
    auto currentData() const -> const KVecs& { return current_data_; }  // current data
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

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
    GKKMetadata meta_;                          // metadata
    std::vector<std::vector<int>> ngtotnod_;    // G-vector count per k-point per node
    std::vector<long> kpt_data_offsets_;        // file offset per k-point
    int current_ikpt_ = -1;                      // currently cached k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::vector<double> kinetic_buf_, Kx_buf_, Ky_buf_, Kz_buf_;
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

    // preallocate working buffers (maximum possible size)
    std::size_t max_ng = static_cast<std::size_t>(meta_.mg_nx) * meta_.nnodes;
    kinetic_buf_.resize(max_ng);
    Kx_buf_.resize(max_ng);
    Ky_buf_.resize(max_ng);
    Kz_buf_.resize(max_ng);

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
    , kinetic_buf_(std::move(other.kinetic_buf_))
    , Kx_buf_(std::move(other.Kx_buf_))
    , Ky_buf_(std::move(other.Ky_buf_))
    , Kz_buf_(std::move(other.Kz_buf_))
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
        kinetic_buf_ = std::move(other.kinetic_buf_);
        Kx_buf_ = std::move(other.Kx_buf_);
        Ky_buf_ = std::move(other.Ky_buf_);
        Kz_buf_ = std::move(other.Kz_buf_);
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
    // Fortran column-major to C++ row-major, AL in file is Angstrom, convert to Bohr
    // Convention B: AL[n][c] where n=vector(a1,a2,a3), c=component(x,y,z), contiguous in memory
    for (int n = 0; n < 3; ++n) {
        for (int c = 0; c < 3; ++c) {
            meta_.AL[n][c] = al_flat[n * 3 + c] / BOHR_RADIUS_ANGSTROM;
        }
    }

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

auto GKK::loadKPoint(int ikpt) -> const KVecs& {
    // check if already in buffer
    if (ikpt == current_ikpt_) {
        return current_data_;
    }

    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }

    seekToKPoint(ikpt);

    // read all nodes for this k-point and merge into contiguous buffers
    std::size_t total_pos = 0;

    for (int inode = 0; inode < meta_.nnodes; ++inode) {
        int ng = ngtotnod_[ikpt][inode];
        if (ng == 0) continue;

        readRecord(kinetic_buf_.data() + total_pos, ng * sizeof(double), "gkk");
        readRecord(Kx_buf_.data() + total_pos, ng * sizeof(double), "gkk_x");
        readRecord(Ky_buf_.data() + total_pos, ng * sizeof(double), "gkk_y");
        readRecord(Kz_buf_.data() + total_pos, ng * sizeof(double), "gkk_z");

        total_pos += static_cast<std::size_t>(ng);
    }

    // update current data view
    current_data_.kinetic = std::span<const double>(kinetic_buf_.data(), total_pos);
    current_data_.Kx      = std::span<const double>(Kx_buf_.data(), total_pos);
    current_data_.Ky      = std::span<const double>(Ky_buf_.data(), total_pos);
    current_data_.Kz      = std::span<const double>(Kz_buf_.data(), total_pos);

    current_ikpt_ = ikpt;
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
            const double c = (meta_.AL[dim][0] * Kx + meta_.AL[dim][1] * Ky + meta_.AL[dim][2] * Kz) / TWO_PI;

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
