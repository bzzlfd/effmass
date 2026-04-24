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
    std::span<const double> g, gx, gy, gz;  // |G+k|²/2, Gx, Gy, Gz
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
    auto currentKPoint() const -> int { return current_kpt_; }     // current k-point index
    auto currentData() const -> const KPointGVecs& { return current_data_; }  // current data

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
    int current_kpt_ = -1;                      // currently cached k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::vector<double> g_buf_, gx_buf_, gy_buf_, gz_buf_;
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
    explicit WG(const std::string& filename);  // open file and read metadata
    ~WG();
    WG(const WG&) = delete;                    // disable copy
    auto operator=(const WG&) -> WG& = delete;
    WG(WG&& other) noexcept;                   // enable move
    auto operator=(WG&& other) noexcept -> WG&;

    auto metadata() const -> const WGMetadata& { return meta_; }   // get metadata
    auto loadBand(int ikpt, int iband) -> const WGCoeffs&;        // load band data (with cache)
    auto currentKPoint() const -> int { return current_kpt_; }     // current k-point index
    auto currentBand() const -> int { return current_band_; }      // current band index
    auto currentData() const -> const WGCoeffs& { return current_data_; }  // current data

private:
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
    int current_kpt_ = -1;                      // currently cached k-point
    int current_band_ = -1;                     // currently cached band

    std::vector<std::complex<double>> up_buf_, down_buf_;
    WGCoeffs current_data_;
};

// Implementation of GKK

GKK::GKK(const std::string& filename)
    : filename_(filename)
    , fp_(nullptr)
    , current_kpt_(-1)
{
    fp_ = std::fopen(filename.c_str(), "rb");
    if (!fp_) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    readMetadata();

    // preallocate working buffers (maximum possible size)
    std::size_t max_ng = static_cast<std::size_t>(meta_.mg_nx) * meta_.nnodes;
    g_buf_.resize(max_ng);
    gx_buf_.resize(max_ng);
    gy_buf_.resize(max_ng);
    gz_buf_.resize(max_ng);

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
    , current_kpt_(other.current_kpt_)
    , g_buf_(std::move(other.g_buf_))
    , gx_buf_(std::move(other.gx_buf_))
    , gy_buf_(std::move(other.gy_buf_))
    , gz_buf_(std::move(other.gz_buf_))
    , current_data_(other.current_data_)
{
    other.fp_ = nullptr;
    other.current_kpt_ = -1;
    // update current_data_ spans to point to our own buffers
    if (!current_data_.g.empty()) {
        const auto ng = current_data_.g.size();
        current_data_.g  = std::span<const double>(g_buf_.data(), ng);
        current_data_.gx = std::span<const double>(gx_buf_.data(), ng);
        current_data_.gy = std::span<const double>(gy_buf_.data(), ng);
        current_data_.gz = std::span<const double>(gz_buf_.data(), ng);
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
        current_kpt_ = other.current_kpt_;
        g_buf_ = std::move(other.g_buf_);
        gx_buf_ = std::move(other.gx_buf_);
        gy_buf_ = std::move(other.gy_buf_);
        gz_buf_ = std::move(other.gz_buf_);
        current_data_ = other.current_data_;

        other.fp_ = nullptr;
        other.current_kpt_ = -1;

        // update spans
        if (!current_data_.g.empty()) {
            const auto ng = current_data_.g.size();
            current_data_.g  = std::span<const double>(g_buf_.data(), ng);
            current_data_.gx = std::span<const double>(gx_buf_.data(), ng);
            current_data_.gy = std::span<const double>(gy_buf_.data(), ng);
            current_data_.gz = std::span<const double>(gz_buf_.data(), ng);
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

    for (int k = 0; k < meta_.nkpt; ++k) {
        int ng_total = 0;
        for (int n = 0; n < meta_.nnodes; ++n) {
            int ng;
            if (std::fread(&ng, sizeof(int), 1, fp_) != 1) {
                throw std::runtime_error("Failed to read ngtotnod");
            }
            if (meta_.is_SO == 1) {
                ng /= 2;
            }
            ngtotnod_[k][n] = ng;
            ng_total += ng;
        }
        meta_.ng_tot_per_kpt[k] = ng_total;
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

    for (int k = 0; k < meta_.nkpt; ++k) {
        kpt_data_offsets_[k] = std::ftell(fp_);

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
    if (ikpt == current_kpt_) {
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

        readRecord(g_buf_.data() + total_pos, ng * sizeof(double), "gkk");
        readRecord(gx_buf_.data() + total_pos, ng * sizeof(double), "gkk_x");
        readRecord(gy_buf_.data() + total_pos, ng * sizeof(double), "gkk_y");
        readRecord(gz_buf_.data() + total_pos, ng * sizeof(double), "gkk_z");

        total_pos += static_cast<std::size_t>(ng);
    }

    // update current data view
    current_data_.g  = std::span<const double>(g_buf_.data(), total_pos);
    current_data_.gx = std::span<const double>(gx_buf_.data(), total_pos);
    current_data_.gy = std::span<const double>(gy_buf_.data(), total_pos);
    current_data_.gz = std::span<const double>(gz_buf_.data(), total_pos);

    current_kpt_ = ikpt;
    return current_data_;
}

// Implementation of WG

WG::WG(const std::string& filename)
    : filename_(filename)
    , fp_(nullptr)
    , current_kpt_(-1)
    , current_band_(-1)
{
    fp_ = std::fopen(filename.c_str(), "rb");
    if (!fp_) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    readMetadata();

    // preallocate working buffers
    std::size_t max_ng = static_cast<std::size_t>(meta_.mg_nx) * meta_.nnodes;
    up_buf_.resize(max_ng);
    if (meta_.is_SO == 1) {
        down_buf_.resize(max_ng);
    }

    current_data_ = {};

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
    , current_kpt_(other.current_kpt_)
    , current_band_(other.current_band_)
    , up_buf_(std::move(other.up_buf_))
    , down_buf_(std::move(other.down_buf_))
    , current_data_(other.current_data_)
{
    other.fp_ = nullptr;
    other.current_kpt_ = -1;
    other.current_band_ = -1;
    if (!current_data_.up.empty()) {
        const auto ng = current_data_.up.size();
        current_data_.up = std::span<const std::complex<double>>(up_buf_.data(), ng);
        if (meta_.is_SO == 1) {
            current_data_.down = std::span<const std::complex<double>>(down_buf_.data(), ng);
        }
    }
}

auto WG::operator=(WG&& other) noexcept -> WG& {
    if (this != &other) {
        if (fp_) std::fclose(fp_);

        filename_ = std::move(other.filename_);
        fp_ = other.fp_;
        meta_ = std::move(other.meta_);
        ngtotnod_ = std::move(other.ngtotnod_);
        band_offsets_ = std::move(other.band_offsets_);
        current_kpt_ = other.current_kpt_;
        current_band_ = other.current_band_;
        up_buf_ = std::move(other.up_buf_);
        down_buf_ = std::move(other.down_buf_);
        current_data_ = other.current_data_;

        other.fp_ = nullptr;
        other.current_kpt_ = -1;
        other.current_band_ = -1;

        if (!current_data_.up.empty()) {
            const auto ng = current_data_.up.size();
            current_data_.up = std::span<const std::complex<double>>(up_buf_.data(), ng);
            if (meta_.is_SO == 1) {
                current_data_.down = std::span<const std::complex<double>>(down_buf_.data(), ng);
            }
        }
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

    for (int k = 0; k < meta_.nkpt; ++k) {
        int ng_total = 0;
        for (int n = 0; n < meta_.nnodes; ++n) {
            int ng;
            if (std::fread(&ng, sizeof(int), 1, fp_) != 1) {
                throw std::runtime_error("Failed to read ngtotnod");
            }
            if (meta_.is_SO == 1) {
                ng /= 2;
            }
            ngtotnod_[k][n] = ng;
            ng_total += ng;
        }
        meta_.ng_tot_per_kpt[k] = ng_total;
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

    for (int k = 0; k < meta_.nkpt; ++k) {
        for (int b = 0; b < meta_.mx; ++b) {
            band_offsets_[static_cast<std::size_t>(k) * meta_.mx + b] = std::ftell(fp_);
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
    if (ikpt == current_kpt_ && iband == current_band_) {
        return current_data_;
    }

    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }
    if (iband < 0 || iband >= meta_.mx) {
        throw std::out_of_range("Invalid band index: " + std::to_string(iband));
    }

    seekToBand(ikpt, iband);

    std::size_t total_pos = 0;

    for (int inode = 0; inode < meta_.nnodes; ++inode) {
        int ng = ngtotnod_[ikpt][inode];
        if (ng == 0) continue;

        if (meta_.is_SO == 1) {
            // record contains both spin-up and spin-down components
            std::vector<std::complex<double>> tmp(static_cast<std::size_t>(ng) * 2);
            readRecord(tmp.data(), tmp.size() * sizeof(std::complex<double>), "wg");
            for (int ig = 0; ig < ng; ++ig) {
                up_buf_[total_pos + ig] = tmp[ig];
                down_buf_[total_pos + ig] = tmp[ng + ig];
            }
        } else {
            readRecord(up_buf_.data() + total_pos, ng * sizeof(std::complex<double>), "wg");
        }

        total_pos += static_cast<std::size_t>(ng);
    }

    current_data_.up = std::span<const std::complex<double>>(up_buf_.data(), total_pos);
    if (meta_.is_SO == 1) {
        current_data_.down = std::span<const std::complex<double>>(down_buf_.data(), total_pos);
    } else {
        current_data_.down = {};
    }

    current_kpt_ = ikpt;
    current_band_ = iband;
    return current_data_;
}
