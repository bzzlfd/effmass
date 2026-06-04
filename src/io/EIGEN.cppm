module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.EIGEN;

import std;
import utils.vector3d;


export {
    class EIGEN;
        struct EIGENMetadata;
}

using kVec = vector3d<double>;


struct EIGENMetadata {
    int islda;
    int nkpt;
    int nband;
    int nref_tot_8;
    int natom;
    int nnode;
    int is_SO;
};


class EIGEN {
public:
    explicit EIGEN(const std::string& filename)
        : meta{}
    {
        auto* fp = std::fopen(filename.c_str(), "rb");
        if (!fp) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        try {
            readMetadata(fp);
        } catch (...) {
            std::fclose(fp);
            throw;
        }

        kpt_weight.resize(meta.nkpt);
        kpt_vec.resize(meta.nkpt);
        eigenvalue_data_.resize(static_cast<std::size_t>(meta.nband) * meta.nkpt * meta.islda);

        try {
            readData(fp);
        } catch (...) {
            std::fclose(fp);
            throw;
        }

        std::fclose(fp);
    }

    EIGEN(const EIGEN&) = default;
    auto operator=(const EIGEN&) -> EIGEN& = default;
    EIGEN(EIGEN&&) noexcept = default;
    auto operator=(EIGEN&&) noexcept -> EIGEN& = default;

    auto operator[](int iband, int ikpt) -> double& {
        if (meta.islda != 1) {
            throw std::runtime_error(
                "EIGEN: islda=" + std::to_string(meta.islda) + " != 1, use [iband, ikpt, ispin]"
            );
        }
        return eigenvalue_data_[eigIndex(iband, ikpt, 0)];
    }
    auto operator[](int iband, int ikpt) const -> double {
        if (meta.islda != 1) {
            throw std::runtime_error(
                "EIGEN: islda=" + std::to_string(meta.islda) + " != 1, use [iband, ikpt, ispin]"
            );
        }
        return eigenvalue_data_[eigIndex(iband, ikpt, 0)];
    }

    auto operator[](int iband, int ikpt, int ispin) -> double& {
        return eigenvalue_data_[eigIndex(iband, ikpt, ispin)];
    }
    auto operator[](int iband, int ikpt, int ispin) const -> double {
        return eigenvalue_data_[eigIndex(iband, ikpt, ispin)];
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
    std::vector<kVec> kpt_vec;

private:
    // Fortran record layout (sequential, no padding): int, int, double x4
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

    auto eigIndex(int iband, int ikpt, int ispin) const -> std::size_t {
        return static_cast<std::size_t>(ispin) * meta.nkpt * meta.nband
             + static_cast<std::size_t>(ikpt) * meta.nband
             + static_cast<std::size_t>(iband);
    }

    auto readRecordLength(std::FILE* fp) -> int {
        int length;
        if (std::fread(&length, sizeof(int), 1, fp) != 1) {
            throw std::runtime_error("Failed to read record length marker");
        }
        return length;
    }

    auto checkRecordLength(std::FILE* fp, int expected) -> void {
        int length;
        if (std::fread(&length, sizeof(int), 1, fp) != 1) {
            throw std::runtime_error("Failed to read record length marker");
        }
        if (length != expected) {
            throw std::runtime_error(
                "Record length mismatch: expected " + std::to_string(expected) +
                ", got " + std::to_string(length)
            );
        }
    }

    auto readRecord(std::FILE* fp, void* dst, std::size_t nbytes, const char* context) -> void {
        int len = readRecordLength(fp);
        if (static_cast<std::size_t>(len) != nbytes) {
            throw std::runtime_error(
                std::string(context) + ": record size mismatch: expected " +
                std::to_string(nbytes) + ", got " + std::to_string(len)
            );
        }
        if (std::fread(dst, 1, nbytes, fp) != nbytes) {
            throw std::runtime_error(std::string(context) + ": read failed");
        }
        checkRecordLength(fp, len);
    }

    auto readMetadata(std::FILE* fp) -> void {
        int len = readRecordLength(fp);

        if (len == static_cast<int>(7 * sizeof(int))) {
            int header[7];
            if (std::fread(header, sizeof(int), 7, fp) != 7) {
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
            if (std::fread(header, sizeof(int), 6, fp) != 6) {
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

        checkRecordLength(fp, len);
    }

    auto readData(std::FILE* fp) -> void {
        for (int ispin = 0; ispin < meta.islda; ++ispin) {
            for (int ikpt = 0; ikpt < meta.nkpt; ++ikpt) {
                KptRecord kpt_rec;
                readRecord(fp, &kpt_rec, sizeof(KptRecord), "kpt metadata");

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

                double* eig_slice = eigenvalue_data_.data()
                    + (static_cast<std::size_t>(ispin) * meta.nkpt + ikpt) * meta.nband;
                readRecord(fp, eig_slice, meta.nband * sizeof(double), "eigenvalues");
                for (int iband = 0; iband < meta.nband; ++iband) {
                    eig_slice[iband] /= HARTREE_TO_EV;  // file stores eV → Hartree
                }
            }
        }
    }

    std::vector<double> eigenvalue_data_;
};
