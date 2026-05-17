module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.EIGEN;

import std;


export struct EIGENMetadata {
    int islda;
    int nkpt;
    int nband;
    int nref_tot_8;
    int natom;
    int nnode;
    int is_SO;
};


export struct KPointVec {
    double x{}, y{}, z{};

    auto operator[](int i) -> double& {
        if (i < 0 || i > 2) {
            throw std::out_of_range("KPointVec index out of range");
        }
        return i == 0 ? x : (i == 1 ? y : z);
    }

    auto operator[](int i) const -> double {
        if (i < 0 || i > 2) {
            throw std::out_of_range("KPointVec index out of range");
        }
        return i == 0 ? x : (i == 1 ? y : z);
    }
};


export class Array3d {
public:
    Array3d() = default;
    Array3d(int nband, int nkpt, int islda)
        : nband_(nband)
        , nkpt_(nkpt)
        , islda_(islda)
        , data_(static_cast<std::size_t>(nband) * nkpt * islda)
    {}

    auto dims() const -> std::array<int, 3> {
        return {nband_, nkpt_, islda_};
    }

    auto empty() const -> bool { return data_.empty(); }
    auto size() const -> std::size_t { return data_.size(); }
    auto data() -> double* { return data_.data(); }
    auto data() const -> const double* { return data_.data(); }

    // ======================================================================== 
    // implement of eigenvalue[iband, ikpt, ispin] (when ISPIN=1,2)
    // ======================================================================== 
    auto operator[](int iband, int ikpt, int ispin) -> double& {
        return data_[index(iband, ikpt, ispin)];
    }
    auto operator[](int iband, int ikpt, int ispin) const -> double {
        return data_[index(iband, ikpt, ispin)];
    }

    // implement of eigenvalue[iband][ikpt] (when ISPIN=1)
    auto operator[](int iband, int ikpt) -> double& {
        if (islda_ != 1) {
            throw std::runtime_error(
                "Array3d: islda=" + std::to_string(islda_) + " != 1, use [iband, ikpt, ispin]"
            );
        }
        return data_[index(iband, ikpt, 0)];
    }
    auto operator[](int iband, int ikpt) const -> double {
        if (islda_ != 1) {
            throw std::runtime_error(
                "Array3d: islda=" + std::to_string(islda_) + " != 1, use [iband, ikpt, ispin]"
            );
        }
        return data_[index(iband, ikpt, 0)];
    }

    // ======================================================================== 
    // implement of Array3d[ispin][ikpt][iband] access via proxy classes
    // ======================================================================== 
    class Slice2d {
    public:
        Slice2d(double* d, int nband, int nkpt)
            : d_(d), nband_(nband), nkpt_(nkpt) {}
        auto operator[](int ikpt) -> std::span<double> {
            return {d_ + ikpt * nband_, static_cast<std::size_t>(nband_)};
        }
    private:
        double* d_;
        int nband_, nkpt_;
    };

    class ConstSlice2d {
    public:
        ConstSlice2d(const double* d, int nband, int nkpt)
            : d_(d), nband_(nband), nkpt_(nkpt) {}
        auto operator[](int ikpt) -> std::span<const double> {
            return {d_ + ikpt * nband_, static_cast<std::size_t>(nband_)};
        }
    private:
        const double* d_;
        int nband_, nkpt_;
    };

    auto operator[](int ispin) -> Slice2d {
        return {data_.data() + ispin * nkpt_ * nband_, nband_, nkpt_};
    }

    auto operator[](int ispin) const -> ConstSlice2d {
        return {data_.data() + ispin * nkpt_ * nband_, nband_, nkpt_};
    }

private:
    auto index(int iband, int ikpt, int ispin) const -> std::size_t {
        return static_cast<std::size_t>(ispin) * nkpt_ * nband_
             + static_cast<std::size_t>(ikpt) * nband_
             + static_cast<std::size_t>(iband);
    }

    int nband_ = 0, nkpt_ = 0, islda_ = 0;
    std::vector<double> data_;
};


export class EIGEN {
public:
    explicit EIGEN(const std::string& filename)
        : meta{}
        , fp_(nullptr)
    {
        fp_ = std::fopen(filename.c_str(), "rb");
        if (!fp_) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        readMetadata();

        kpt_weight.resize(meta.nkpt);
        kpt_vec.resize(meta.nkpt);
        eigenvalue = Array3d(meta.nband, meta.nkpt, meta.islda);

        readData();
    }

    ~EIGEN() {
        if (fp_) {
            std::fclose(fp_);
        }
    }

    EIGEN(const EIGEN&) = delete;
    auto operator=(const EIGEN&) -> EIGEN& = delete;

    EIGEN(EIGEN&& other) noexcept
        : meta(std::move(other.meta))
        , kpt_weight(std::move(other.kpt_weight))
        , kpt_vec(std::move(other.kpt_vec))
        , eigenvalue(std::move(other.eigenvalue))
        , fp_(other.fp_)
    {
        other.fp_ = nullptr;
    }

    auto operator=(EIGEN&& other) noexcept -> EIGEN& {
        if (this != &other) {
            if (fp_) std::fclose(fp_);
            meta = std::move(other.meta);
            kpt_weight = std::move(other.kpt_weight);
            kpt_vec = std::move(other.kpt_vec);
            eigenvalue = std::move(other.eigenvalue);
            fp_ = other.fp_;
            other.fp_ = nullptr;
        }
        return *this;
    }

    auto print_info() const -> void {
        double sum_weight = 0.0;
        for (const auto& w : kpt_weight) {
            sum_weight += w;
        }
        std::println("EIGEN metadata:");
        std::println("  islda      = {}", meta.islda);
        std::println("  nkpt       = {}", meta.nkpt);
        std::println("  nband      = {}", meta.nband);
        std::println("  nref_tot_8 = {}", meta.nref_tot_8);
        std::println("  natom      = {}", meta.natom);
        std::println("  nnode      = {}", meta.nnode);
        std::println("  is_SO      = {}", meta.is_SO);
        std::println("  sum(weight_kpt) = {}", sum_weight);
    }

    EIGENMetadata meta;
    std::vector<double> kpt_weight;
    std::vector<KPointVec> kpt_vec;
    Array3d eigenvalue;

private:
    // Fortran record layout (sequential, no padding): int, int, double x4
    // C++ layout: int(4) + int(4) + double(8) + double(8) + double(8) + double(8) = 40 bytes
    // No padding: all ints precede the first double, maintaining natural 8-byte alignment.
    // sizeof(KptRecord) == 2*sizeof(int) + 4*sizeof(double) verified by static_assert below.
    struct KptRecord {
        int islda_tmp;
        int ikpt_tmp;
        double weight;
        double akx;
        double aky;
        double akz;
    };
    static_assert(sizeof(KptRecord) == 2 * sizeof(int) + 4 * sizeof(double),
                  "KptRecord layout mismatch");

    auto readRecordLength() -> int {
        int length;
        if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
            throw std::runtime_error("Failed to read record length marker");
        }
        return length;
    }

    auto checkRecordLength(int expected) -> void {
        int length;
        if (std::fread(&length, sizeof(int), 1, fp_) != 1) {
            throw std::runtime_error("Failed to read record length marker");
        }
        if (length != expected) {
            throw std::runtime_error(
                "Record length mismatch: expected " + std::to_string(expected) +
                ", got " + std::to_string(length)
            );
        }
    }

    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
        int len = readRecordLength();
        if (static_cast<std::size_t>(len) != nbytes) {
            throw std::runtime_error(
                std::string(context) + ": record size mismatch: expected " +
                std::to_string(nbytes) + ", got " + std::to_string(len)
            );
        }
        if (std::fread(dst, 1, nbytes, fp_) != nbytes) {
            throw std::runtime_error(std::string(context) + ": read failed");
        }
        checkRecordLength(len);
    }

    auto readMetadata() -> void {
        int len = readRecordLength();

        if (len == static_cast<int>(7 * sizeof(int))) {
            int header[7];
            if (std::fread(header, sizeof(int), 7, fp_) != 7) {
                throw std::runtime_error("Failed to read EIGEN header (7 ints)");
            }
            meta.islda     = header[0];
            meta.nkpt      = header[1];
            meta.nband     = header[2];
            meta.nref_tot_8 = header[3];
            meta.natom     = header[4];
            meta.nnode     = header[5];
            meta.is_SO     = header[6];
        } else if (len == static_cast<int>(6 * sizeof(int))) {
            int header[6];
            if (std::fread(header, sizeof(int), 6, fp_) != 6) {
                throw std::runtime_error("Failed to read EIGEN header (6 ints)");
            }
            meta.islda     = header[0];
            meta.nkpt      = header[1];
            meta.nband     = header[2];
            meta.nref_tot_8 = header[3];
            meta.natom     = header[4];
            meta.nnode     = header[5];
            meta.is_SO     = 0;
        } else {
            throw std::runtime_error(
                "Unexpected EIGEN header record length: " + std::to_string(len) + " bytes"
            );
        }

        checkRecordLength(len);
    }

    auto readData() -> void {
        for (int ispin = 0; ispin < meta.islda; ++ispin) {
            for (int ikpt = 0; ikpt < meta.nkpt; ++ikpt) {
                KptRecord kpt_rec;
                readRecord(&kpt_rec, sizeof(KptRecord), "kpt metadata");

                if (kpt_rec.islda_tmp != ispin + 1) {
                    throw std::runtime_error(
                        "EIGEN islda_tmp mismatch: expected " +
                        std::to_string(ispin + 1) + ", got " +
                        std::to_string(kpt_rec.islda_tmp)
                    );
                }
                if (kpt_rec.ikpt_tmp != ikpt + 1) {
                    throw std::runtime_error(
                        "EIGEN ikpt_tmp mismatch: expected " +
                        std::to_string(ikpt + 1) + ", got " +
                        std::to_string(kpt_rec.ikpt_tmp)
                    );
                }

                kpt_weight[ikpt] = kpt_rec.weight;
                kpt_vec[ikpt] = {kpt_rec.akx, kpt_rec.aky, kpt_rec.akz};

                double* eig_slice = eigenvalue.data()
                    + (static_cast<std::size_t>(ispin) * meta.nkpt + ikpt) * meta.nband;
                readRecord(eig_slice, meta.nband * sizeof(double), "eigenvalues");
            }
        }
    }

    std::FILE* fp_ = nullptr;
};
