module;
#include <cstdio>
#include "../physical_constants.hpp"

export module io.RHO;

import io.lattice;
import std;


export struct RHOMetadata {
    int n1, n2, n3;
    int nnodes;
    int nstate;
};


export class RHO {
public:
    explicit RHO(const std::string& filename)
        : fp_(nullptr)
    {
        fp_ = std::fopen(filename.c_str(), "rb");
        if (!fp_) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        readMetadata();
        readData();
    }

    ~RHO() {
        if (fp_) {
            std::fclose(fp_);
        }
    }

    RHO(const RHO&) = delete;
    auto operator=(const RHO&) -> RHO& = delete;

    RHO(RHO&& other) noexcept
        : meta(other.meta)
        , lattice(std::move(other.lattice))
        , data_(std::move(other.data_))
        , fp_(other.fp_)
    {
        other.fp_ = nullptr;
    }

    auto operator=(RHO&& other) noexcept -> RHO& {
        if (this != &other) {
            if (fp_) std::fclose(fp_);
            meta = other.meta;
            lattice = std::move(other.lattice);
            data_ = std::move(other.data_);
            fp_ = other.fp_;
            other.fp_ = nullptr;
        }
        return *this;
    }

    // [i, j, k] access state 0, no runtime check, direct index
    auto operator[](int i, int j, int k) -> double& {
        return data_[static_cast<std::size_t>(i) * meta.n2 * meta.n3
                   + static_cast<std::size_t>(j) * meta.n3
                   + static_cast<std::size_t>(k)];
    }

    auto operator[](int i, int j, int k) const -> double {
        return data_[static_cast<std::size_t>(i) * meta.n2 * meta.n3
                   + static_cast<std::size_t>(j) * meta.n3
                   + static_cast<std::size_t>(k)];
    }

    // [i, j, k, state] explicit state access
    auto operator[](int i, int j, int k, int state) -> double& {
        return data_[(static_cast<std::size_t>(state) * meta.n1 + i) * meta.n2 * meta.n3
                   + static_cast<std::size_t>(j) * meta.n3
                   + static_cast<std::size_t>(k)];
    }

    auto operator[](int i, int j, int k, int state) const -> double {
        return data_[(static_cast<std::size_t>(state) * meta.n1 + i) * meta.n2 * meta.n3
                   + static_cast<std::size_t>(j) * meta.n3
                   + static_cast<std::size_t>(k)];
    }

    auto print_info() const -> void {
        std::println("RHO/VR metadata:");
        std::println("  n1     = {}", meta.n1);
        std::println("  n2     = {}", meta.n2);
        std::println("  n3     = {}", meta.n3);
        std::println("  nnodes = {}", meta.nnodes);
        std::println("  nstate = {}", meta.nstate);
        std::println("  nr     = {}", meta.n1 * meta.n2 * meta.n3);
    }

    RHOMetadata meta;
    Lattice lattice;

private:
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

    auto readRecord(void* dst, const char* context) -> void {
        int len = readRecordLength();
        if (std::fread(dst, 1, len, fp_) != static_cast<std::size_t>(len)) {
            throw std::runtime_error(std::string(context) + ": read failed");
        }
        checkRecordLength(len);
    }

    auto readRecord(void* dst, std::size_t nbytes, const char* context) -> void {
        int len = readRecordLength();
        if (static_cast<int>(nbytes) > len) {
            throw std::runtime_error(
                std::string(context) + ": nbytes exceeds record length: expected ≤ " +
                std::to_string(len) + ", got " + std::to_string(nbytes)
            );
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

    auto readMetadata() -> void {
        int len = readRecordLength();

        if (len == static_cast<int>(5 * sizeof(int))) {
            int header[5];
            if (std::fread(header, sizeof(int), 5, fp_) != 5) {
                throw std::runtime_error("Failed to read RHO header (5 ints)");
            }
            meta.n1      = header[0];
            meta.n2      = header[1];
            meta.n3      = header[2];
            meta.nnodes  = header[3];
            meta.nstate  = header[4];
        } else if (len == static_cast<int>(4 * sizeof(int))) {
            int header[4];
            if (std::fread(header, sizeof(int), 4, fp_) != 4) {
                throw std::runtime_error("Failed to read RHO header (4 ints)");
            }
            meta.n1      = header[0];
            meta.n2      = header[1];
            meta.n3      = header[2];
            meta.nnodes  = header[3];
            meta.nstate  = 1;
        } else {
            throw std::runtime_error(
                "Unexpected RHO header record length: " + std::to_string(len) + " bytes"
            );
        }

        checkRecordLength(len);

        // Read lattice vectors (Angstrom in file, converted to Bohr internally)
        double al_flat[9];
        readRecord(al_flat, "AL");
        lattice.setLattice(al_flat, LengthUnit::Angstrom);
    }

    auto readData() -> void {
        int nr = meta.n1 * meta.n2 * meta.n3;
        int nr_n = nr / meta.nnodes;
        data_.resize(static_cast<std::size_t>(meta.nstate) * nr);

        auto nbytes = static_cast<std::size_t>(nr_n) * sizeof(double);
        for (int istate = 0; istate < meta.nstate; ++istate) {
            for (int inode = 0; inode < meta.nnodes; ++inode) {
                double* slice = data_.data()
                    + (static_cast<std::size_t>(istate) * nr + inode * nr_n);
                readRecord(slice, nbytes, "grid data");
            }
        }
    }

    std::FILE* fp_ = nullptr;
    std::vector<double> data_;
};
