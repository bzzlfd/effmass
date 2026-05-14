module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.WG;

import io.lattice;
import std;

// WG file metadata structure
export struct WGMetadata {
    int n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / bands / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
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
    meta_.lattice = Lattice(al_flat, LengthUnit::Angstrom);

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
