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

// k-point G-vector data view - pointers to contiguous memory
export struct KPointGVecs {
    std::size_t ng;           // number of G-vectors
    const double *g, *gx, *gy, *gz;  // |G+k|²/2, Gx, Gy, Gz
};

// GKK class - abstraction for OUT.GKK file
export class GKK {
public:
    explicit GKK(const std::string& filename);  // open file and read metadata
    ~GKK();
    GKK(const GKK&) = delete;                   // disable copy
    GKK& operator=(const GKK&) = delete;
    GKK(GKK&& other) noexcept;                  // enable move
    GKK& operator=(GKK&& other) noexcept;

    const GKKMetadata& metadata() const { return meta_; }  // get metadata
    const KPointGVecs& loadKPoint(int ikpt);    // load k-point data (with cache)
    int currentKPoint() const { return current_kpt_; }     // current k-point index
    const KPointGVecs& currentData() const { return current_data_; }  // current data

private:
    int readRecordLength();                     // read record length marker
    void checkRecordLength(int expected);       // verify record length marker
    void readRecord(void* dst, std::size_t nbytes, const char* context); // read full record data

    void readMetadata();                        // read file metadata
    void readNgtotnod(int record_len);          // read ngtotnod array
    void skipRecord();                            // skip one Fortran record by reading its length markers
    void computeOffsets();                        // compute file offset per k-point
    void seekToKPoint(int ikpt);                // seek to k-point data

    std::string filename_;                      // file name
    std::FILE* fp_;                             // file handle
    GKKMetadata meta_;                          // metadata
    std::vector<std::vector<int>> ngtotnod_;    // G-vector count per k-point per node
    std::vector<long> kpt_data_offsets_;        // file offset per k-point
    int current_kpt_ = -1;                      // currently cached k-point

    // buffers: working arrays (contiguous) + file read buffer (reused)
    std::vector<double> g_buf_, gx_buf_, gy_buf_, gz_buf_;
    KPointGVecs current_data_{0, nullptr, nullptr, nullptr, nullptr};  // data view
};

// Implementation

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

    // initialize data view pointers
    current_data_ = {0, nullptr, nullptr, nullptr, nullptr};

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
    // update current_data_ pointers to point to our own buffers
    if (current_data_.ng > 0) {
        current_data_.g = g_buf_.data();
        current_data_.gx = gx_buf_.data();
        current_data_.gy = gy_buf_.data();
        current_data_.gz = gz_buf_.data();
    }
}

GKK& GKK::operator=(GKK&& other) noexcept {
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

        // update pointers
        if (current_data_.ng > 0) {
            current_data_.g = g_buf_.data();
            current_data_.gx = gx_buf_.data();
            current_data_.gy = gy_buf_.data();
            current_data_.gz = gz_buf_.data();
        }
    }
    return *this;
}

int GKK::readRecordLength() {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    return length;
}

void GKK::checkRecordLength(int expected_length) {
    int length;
    if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
        throw std::runtime_error("Failed to read record length marker");
    }
    if (length != expected_length) {
        throw std::runtime_error("Record length mismatch");
    }
}

void GKK::readRecord(void* dst, std::size_t nbytes, const char* context) {
    int len = readRecordLength();
    if (len != static_cast<int>(nbytes)) {
        throw std::runtime_error(std::string(context) + ": record size mismatch");
    }
    if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
        throw std::runtime_error(std::string(context) + ": read failed");
    }
    checkRecordLength(len);
}

void GKK::readMetadata() {
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

void GKK::readNgtotnod(int record_len) {
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

void GKK::skipRecord() {
    int len = readRecordLength();
    if (std::fseek(fp_, len, SEEK_CUR) != 0) {
        throw std::runtime_error("Failed to skip record data");
    }
    checkRecordLength(len);
}

void GKK::computeOffsets() {
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

void GKK::seekToKPoint(int ikpt) {
    if (ikpt < 0 || ikpt >= meta_.nkpt) {
        throw std::out_of_range("Invalid k-point index");
    }

    if (std::fseek(fp_, kpt_data_offsets_[ikpt], SEEK_SET) != 0) {
        throw std::runtime_error("Failed to seek to k-point");
    }
}

const KPointGVecs& GKK::loadKPoint(int ikpt) {
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
    current_data_.ng = total_pos;
    current_data_.g = g_buf_.data();
    current_data_.gx = gx_buf_.data();
    current_data_.gy = gy_buf_.data();
    current_data_.gz = gz_buf_.data();

    current_kpt_ = ikpt;
    return current_data_;
}
