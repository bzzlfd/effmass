module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.WG;

import io.lattice;
import std;


export {
    class WG;
        struct WGMetadata;
        struct WGCoeffs;
}


// WG file metadata structure
struct WGMetadata {
    int n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda;  // FFT grid / bands / record length / node / k-point / spin
    double Ecut;              // cutoff energy
    Lattice lattice;          // lattice vectors (Bohr) and reciprocal lattice (Bohr^-1)
    std::vector<int> ng_tot_per_kpt;  // total G-vectors per k-point
};


// Wavefunction coefficient view for a single k-point and band
struct WGCoeffs {
    std::span<const std::complex<double>> up;     // spin-up component (always valid)
    std::span<const std::complex<double>> down;   // spin-down component (valid only when is_SO == 1)
};


// WG class - abstraction for OUT.WG file
class WG {
public:
    WGMetadata meta;                              // metadata

    explicit WG(const std::string& filename, std::size_t cache_capacity = 4);  // open file and read metadata
    ~WG();
    WG(const WG&) = delete;                    // disable copy
    auto operator=(const WG&) -> WG& = delete;
    WG(WG&& other) noexcept;                   // enable move
    auto operator=(WG&& other) noexcept -> WG&;

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
    auto readRecord(void* dst, const char* context) -> void; // read full record
    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void; // read nbytes, skip rest
    auto readRecordSOC(std::complex<double>* dst_up, std::complex<double>* dst_down,
                       int ng, const char* context) -> void; // read SOC record: up[0..ng) + down[0..ng)

    auto readMetadata() -> void;                        // read file metadata
    auto readNgtotnod(int record_len) -> void;          // read ngtotnod array
    auto skipRecord() -> void;                          // skip one Fortran record
    auto computeOffsets() -> void;                      // compute file offset per (k-point, band)
    auto seekToBand(int ikpt, int iband) -> void;       // seek to band data

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
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
    , meta(std::move(other.meta))
    , ngtotnod_(std::move(other.ngtotnod_))
    , band_offsets_(std::move(other.band_offsets_))
    , cache_(std::move(other.cache_))
    , cache_capacity_(other.cache_capacity_)
    , cache_next_slot_(other.cache_next_slot_)
    , current_ikpt_(other.current_ikpt_)
    , current_iband_(other.current_iband_)
    , current_data_view_(other.current_data_view_)
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
        meta = std::move(other.meta);
        ngtotnod_ = std::move(other.ngtotnod_);
        band_offsets_ = std::move(other.band_offsets_);
        cache_ = std::move(other.cache_);
        cache_capacity_ = other.cache_capacity_;
        cache_next_slot_ = other.cache_next_slot_;
        current_ikpt_ = other.current_ikpt_;
        current_iband_ = other.current_iband_;
        current_data_view_ = other.current_data_view_;

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


auto WG::readRecord(void* dst, const char* context) -> void {
    int len = readRecordLength();
    if (std::fread(dst, 1, len, fp_) != static_cast<std::size_t>(len)) {
        throw std::runtime_error(std::string(context) + ": read failed");
    }
    checkRecordLength(len);
}

auto WG::readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
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

auto WG::readRecordSOC(std::complex<double>* dst_up, std::complex<double>* dst_down,
                       int ng, const char* context) -> void {
    int len = readRecordLength();
    std::size_t cplx_size = sizeof(std::complex<double>);
    std::size_t half = static_cast<std::size_t>(len) / (2 * cplx_size);

    auto read_cplx = [this, cplx_size](std::complex<double>* dst, int n, const char* ctx) {
        auto nbytes = static_cast<std::size_t>(n) * cplx_size;
        if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
            throw std::runtime_error(std::string(ctx) + ": read failed");
        }
    };

    auto skip_cplx = [this, cplx_size](std::size_t n) {
        if (n > 0) {
            std::fseek(fp_, static_cast<long>(n) * cplx_size, SEEK_CUR);
        }
    };

    read_cplx(dst_up, ng, context);
    skip_cplx(half - static_cast<std::size_t>(ng));
    read_cplx(dst_down, ng, context);
    skip_cplx(half - static_cast<std::size_t>(ng));

    checkRecordLength(len);
}


auto WG::readMetadata() -> void {
    // Record 1: n1, n2, n3, mx, mg_nx, nnodes, nkpt, is_SO, islda (9 ints)
    int header[9];
    readRecord(header, "header");
    meta.n1 = header[0];
    meta.n2 = header[1];
    meta.n3 = header[2];
    meta.mx = header[3];
    meta.mg_nx = header[4];
    meta.nnodes = header[5];
    meta.nkpt = header[6];
    meta.is_SO = header[7];
    meta.islda = header[8];

    if (meta.is_SO == 1) {
        meta.mg_nx /= 2;
    }

    // Record 2: Ecut
    readRecord(&meta.Ecut, "Ecut");

    // Record 3: AL(3,3) - Fortran column-major to C++ row-major, Å to Bohr
    double al_flat[9];
    readRecord(al_flat, "AL");
    meta.lattice.setLattice(al_flat, LengthUnit::Angstrom);

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
    if (nnodes_check != meta.nnodes) {
        throw std::runtime_error("nnodes mismatch");
    }

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


auto WG::skipRecord() -> void {
    int len = readRecordLength();
    if (std::fseek(fp_, len, SEEK_CUR) != 0) {
        throw std::runtime_error("Failed to skip record data");
    }
    checkRecordLength(len);
}


auto WG::computeOffsets() -> void {
    // record the starting file offset for each (k-point, band) pair
    // File layout per kpt: node=0 band=0..mx-1, node=1 band=0..mx-1, ...
    band_offsets_.resize(static_cast<std::size_t>(meta.nkpt) * meta.mx);

    for (int ikpt = 0; ikpt < meta.nkpt; ++ikpt) {
        for (int b = 0; b < meta.mx; ++b) {
            band_offsets_[static_cast<std::size_t>(ikpt) * meta.mx + b] = std::ftell(fp_);
            for (int n = 0; n < meta.nnodes; ++n) {
                skipRecord();
            }
        }
    }
}


auto WG::seekToBand(int ikpt, int iband) -> void {
    if (ikpt < 0 || ikpt >= meta.nkpt) {
        throw std::out_of_range("Invalid k-point index");
    }
    if (iband < 0 || iband >= meta.mx) {
        throw std::out_of_range("Invalid band index");
    }

    long offset = band_offsets_[static_cast<std::size_t>(ikpt) * meta.mx + iband];
    if (std::fseek(fp_, offset, SEEK_SET) != 0) {
        throw std::runtime_error("Failed to seek to band");
    }
}


auto WG::loadBand(int ikpt, int iband) -> const WGCoeffs& {
    if (ikpt < 0 || ikpt >= meta.nkpt) {
        throw std::out_of_range("Invalid k-point index: " + std::to_string(ikpt));
    }
    if (iband < 0 || iband >= meta.mx) {
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

    const std::size_t total_ng = meta.ng_tot_per_kpt[ikpt];
    entry->up.resize(total_ng);
    if (meta.is_SO == 1) {
        entry->down.resize(total_ng);
    }

    std::size_t pos = 0;
    for (int inode = 0; inode < meta.nnodes; ++inode) {
        int ng = ngtotnod_[ikpt][inode];
        if (ng == 0) continue;

        if (meta.is_SO == 1) {
            readRecordSOC(entry->up.data() + pos, entry->down.data() + pos, ng, "wg");
        } else {
            auto nbytes = static_cast<std::size_t>(ng) * sizeof(std::complex<double>);
            readRecord(entry->up.data() + pos, nbytes, "wg");
        }

        pos += static_cast<std::size_t>(ng);
    }

    entry->view.up = std::span<const std::complex<double>>(entry->up.data(), entry->up.size());
    if (meta.is_SO == 1) {
        entry->view.down = std::span<const std::complex<double>>(entry->down.data(), entry->down.size());
    } else {
        entry->view.down = {};
    }

    current_ikpt_ = ikpt;
    current_iband_ = iband;
    current_data_view_ = entry->view;
    return current_data_view_;
}
