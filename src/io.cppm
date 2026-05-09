module;
#include <cstdio>

export module io;

import std;

constexpr double BOHR_RADIUS_ANGSTROM = 0.52917721067;  // 1 Bohr = 0.52917721067 Å

// GKK file metadata structure
export struct GKKMetadata {
    int n1, n2, n3, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    double AL[3][3];          // lattice vectors (Bohr)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};

// k-point G-vector data view - non-owning spans to contiguous memory
export struct KPointGVecs {
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
    auto loadKPoint(int ikpt) -> const KPointGVecs&;    // load k-point data (with cache)
    auto current_ikpt() const -> int { return current_ikpt_; }     // current k-point index
    auto currentData() const -> const KPointGVecs& { return current_data_; }  // current data
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
    KPointGVecs current_data_;  // data view
};

// WG file metadata structure
export struct WGMetadata {
    int n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / bands / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    double AL[3][3];          // lattice vectors (Bohr)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};

// Wavefunction coefficient view for a single k-point and band
export struct WGCoeffs {
    std::span<const std::complex<double>> up;     // spin-up component (always valid)
    std::span<const std::complex<double>> down;   // spin-down component (valid only when is_SO == 1)
};

// WG class - abstraction for OUT.WG file
export class WG {
public:
    explicit WG(const std::string& filename, std::size_t cache_capacity = 4);  // open file and read metadata
    ~WG();
    WG(const WG&) = delete;                    // disable copy
    auto operator=(const WG&) -> WG& = delete;
    WG(WG&& other) noexcept;                   // enable move
    auto operator=(WG&& other) noexcept -> WG&;

    auto metadata() const -> const WGMetadata& { return meta_; }   // get metadata
    auto loadBand(int ikpt, int iband) -> const WGCoeffs&;        // load band data (with cache)
    auto current_ikpt() const -> int { return current_ikpt_; }     // current k-point index
    auto current_iband() const -> int { return current_iband_; }      // current band index
    auto currentData() const -> const WGCoeffs& { return current_data_view_; }  // current data

private:
    struct CacheEntry {
        int ikpt = -1;
        int iband = -1;
        std::vector<std::complex<double>> up;
        std::vector<std::complex<double>> down;
        WGCoeffs view;
    };

    auto readRecordLength() -> int;                     // read record length marker
    auto checkRecordLength(int expected) -> void;       // verify record length marker
    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void; // read full record data

    auto readMetadata() -> void;                        // read file metadata
    auto readNgtotnod(int record_len) -> void;          // read ngtotnod array
    auto skipRecord() -> void;                          // skip one Fortran record
    auto computeOffsets() -> void;                      // compute file offset per (k-point, band)
    auto seekToBand(int ikpt, int iband) -> void;       // seek to band data

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
    WGMetadata meta_;                           // metadata
    std::vector<std::vector<int>> ngtotnod_;    // G-vector count per k-point per node
    std::vector<long> band_offsets_;            // file offset per band per k-point

    // multi-band cache
    std::vector<std::unique_ptr<CacheEntry>> cache_;  // cache entries (stable addresses via unique_ptr)
    std::size_t cache_capacity_;                      // max cached bands
    std::size_t cache_next_slot_ = 0;                 // next slot to replace in FIFO order

    // track most recent loadBand
    int current_ikpt_ = -1;
    int current_iband_ = -1;
    WGCoeffs current_data_view_;

    std::vector<std::complex<double>> tmp_buf_;  // temporary buffer for is_SO record splitting
};

// Implementation of GKK

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

auto GKK::loadKPoint(int ikpt) -> const KPointGVecs& {
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

// Implementation of WG

WG::WG(const std::string& filename, std::size_t cache_capacity)
    : filename_(filename)
    , fp_(nullptr)
    , cache_capacity_(cache_capacity)
    , cache_next_slot_(0)
    , current_ikpt_(-1)
    , current_iband_(-1)
{
    fp_ = std::fopen(filename.c_str(), "rb");
    if (!fp_) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    readMetadata();

    current_data_view_ = {};

    computeOffsets();
}

WG::~WG() {
    if (fp_) {
        std::fclose(fp_);
    }
}

WG::WG(WG&& other) noexcept
    : filename_(std::move(other.filename_))
    , fp_(other.fp_)
    , meta_(std::move(other.meta_))
    , ngtotnod_(std::move(other.ngtotnod_))
    , band_offsets_(std::move(other.band_offsets_))
    , cache_(std::move(other.cache_))
    , cache_capacity_(other.cache_capacity_)
    , cache_next_slot_(other.cache_next_slot_)
    , current_ikpt_(other.current_ikpt_)
    , current_iband_(other.current_iband_)
    , current_data_view_(other.current_data_view_)
    , tmp_buf_(std::move(other.tmp_buf_))
{
    other.fp_ = nullptr;
    other.current_ikpt_ = -1;
    other.current_iband_ = -1;
    other.current_data_view_ = {};
}

auto WG::operator=(WG&& other) noexcept -> WG& {
    if (this != &other) {
        if (fp_) std::fclose(fp_);

        filename_ = std::move(other.filename_);
        fp_ = other.fp_;
        meta_ = std::move(other.meta_);
        ngtotnod_ = std::move(other.ngtotnod_);
        band_offsets_ = std::move(other.band_offsets_);
        cache_ = std::move(other.cache_);
        cache_capacity_ = other.cache_capacity_;
        cache_next_slot_ = other.cache_next_slot_;
        current_ikpt_ = other.current_ikpt_;
        current_iband_ = other.current_iband_;
        current_data_view_ = other.current_data_view_;
        tmp_buf_ = std::move(other.tmp_buf_);

        other.fp_ = nullptr;
        other.current_ikpt_ = -1;
        other.current_iband_ = -1;
        other.current_data_view_ = {};
    }
    return *this;
}

auto WG::readRecordLength() -> int {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    return length;
}

auto WG::checkRecordLength(int expected_length) -> void {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    if (length != expected_length) {
        throw std::runtime_error("Record length mismatch");
    }
}

auto WG::readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
    int len = readRecordLength();
    if (len != static_cast<int>(nbytes)) {
        throw std::runtime_error(std::string(context) + ": record size mismatch");
    }
    if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
        throw std::runtime_error(std::string(context) + ": read failed");
    }
    checkRecordLength(len);
}

auto WG::readMetadata() -> void {
    // Record 1: n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda (9 ints)
    int header[9];
    readRecord(header, sizeof(header), "header");
    meta_.n1 = header[0];
    meta_.n2 = header[1];
    meta_.n3 = header[2];
    meta_.mx = header[3];
    meta_.mg_nx = header[4];
    meta_.nnodes = header[5];
    meta_.nkpt = header[6];
    meta_.is_SO = header[7];
    meta_.islda = header[8];

    if (meta_.is_SO == 1) {
        meta_.mg_nx /= 2;
    }

    // Record 2: Ecut
    readRecord(&meta_.Ecut, sizeof(double), "Ecut");

    // Record 3: AL(3,3) - Fortran column-major to C++ row-major, Å to Bohr
    double al_flat[9];
    readRecord(al_flat, sizeof(al_flat), "AL");
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

auto WG::readNgtotnod(int record_len) -> void {
    int nnodes_check;
    if (std::fread(&nnodes_check, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read nnodes");
    }
    if (nnodes_check != meta_.nnodes) {
        throw std::runtime_error("nnodes mismatch");
    }

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

auto WG::skipRecord() -> void {
    int len = readRecordLength();
    if (std::fseek(fp_, len, SEEK_CUR) != 0) {
        throw std::runtime_error("Failed to skip record data");
    }
    checkRecordLength(len);
}

auto WG::computeOffsets() -> void {
    // record the starting file offset for each (k-point, band) pair
    band_offsets_.resize(static_cast<std::size_t>(meta_.nkpt) * meta_.mx);

    for (int ikpt = 0; ikpt < meta_.nkpt; ++ikpt) {
        for (int b = 0; b < meta_.mx; ++b) {
            band_offsets_[static_cast<std::size_t>(ikpt) * meta_.mx + b] = std::ftell(fp_);
            for (int n = 0; n < meta_.nnodes; ++n) {
                skipRecord();
            }
        }
    }
}

auto WG::seekToBand(int ikpt, int iband) -> void {
    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index");
    }
    if (iband < 0 || iband >= meta_.mx) {
        throw std::out_of_range("Invalid band index");
    }

    long offset = band_offsets_[static_cast<std::size_t>(ikpt) * meta_.mx + iband];
    if (std::fseek(fp_, offset, SEEK_SET) != 0) {
        throw std::runtime_error("Failed to seek to band");
    }
}

auto WG::loadBand(int ikpt, int iband) -> const WGCoeffs& {
    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }
    if (iband < 0 || iband >= meta_.mx) {
        throw std::out_of_range("Invalid band index: " + std::to_string(iband));
    }

    // check cache
    for (const auto& entry : cache_) {
        if (entry->ikpt == ikpt && entry->iband == iband) {
            current_ikpt_ = ikpt;
            current_iband_ = iband;
            current_data_view_ = entry->view;
            return current_data_view_;
        }
    }

    // cache miss: read from file
    seekToBand(ikpt, iband);

    CacheEntry* entry = nullptr;
    if (cache_.size() < cache_capacity_) {
        cache_.emplace_back(std::make_unique<CacheEntry>());
        entry = cache_.back().get();
    } else {
        entry = cache_[cache_next_slot_].get();
        cache_next_slot_ = (cache_next_slot_ + 1) % cache_capacity_;
    }

    entry->ikpt = ikpt;
    entry->iband = iband;

    const std::size_t total_ng = meta_.ng_tot_per_kpt[ikpt];
    entry->up.resize(total_ng);
    if (meta_.is_SO == 1) {
        entry->down.resize(total_ng);
    }

    std::size_t pos = 0;
    for (int inode = 0; inode < meta_.nnodes; ++inode) {
        int ng = ngtotnod_[ikpt][inode];
        if (ng == 0) continue;

        if (meta_.is_SO == 1) {
            tmp_buf_.resize(static_cast<std::size_t>(ng) * 2);
            readRecord(tmp_buf_.data(), tmp_buf_.size() * sizeof(std::complex<double>), "wg");
            for (int ig = 0; ig < ng; ++ig) {
                entry->up[pos + ig] = tmp_buf_[ig];
                entry->down[pos + ig] = tmp_buf_[ng + ig];
            }
        } else {
            readRecord(entry->up.data() + pos, ng * sizeof(std::complex<double>), "wg");
        }
        pos += static_cast<std::size_t>(ng);
    }

    entry->view.up = std::span<const std::complex<double>>(entry->up.data(), entry->up.size());
    if (meta_.is_SO == 1) {
        entry->view.down = std::span<const std::complex<double>>(entry->down.data(), entry->down.size());
    } else {
        entry->view.down = {};
    }

    current_ikpt_ = ikpt;
    current_iband_ = iband;
    current_data_view_ = entry->view;
    return current_data_view_;
}
