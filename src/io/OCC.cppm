module;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

export module io.OCC;

import std;
import utils.vector3d;


export {
    class OCC;
        struct OCCMetadata;
}

using kVec = vector3d<double>;


struct OCCMetadata {
    int islda{};
    int nkpt{};
    int nband{};
};


class OCC {
public:
    // Open file and read all data (text format, eager loading)
    explicit OCC(const std::string& filename);

    // Rule of zero: no FILE* or other resource handle
    OCC() = default;
    OCC(const OCC&) = default;
    auto operator=(const OCC&) -> OCC& = default;
    OCC(OCC&&) noexcept = default;
    auto operator=(OCC&&) noexcept -> OCC& = default;

    auto energy(int iband, int ikpt) const -> double {
        return energies_[idx(iband, ikpt, 0)];
    }

    auto energy(int iband, int ikpt, int ispin) const -> double {
        return energies_[idx(iband, ikpt, ispin)];
    }

    auto occupation(int iband, int ikpt) const -> double {
        return occupations_[idx(iband, ikpt, 0)];
    }

    auto occupation(int iband, int ikpt, int ispin) const -> double {
        return occupations_[idx(iband, ikpt, ispin)];
    }

    auto print_info() const -> void {
        std::println("OCC: islda = {}, nkpt = {}, nband = {}",
                     meta.islda, meta.nkpt, meta.nband);
        std::println("  k-points:");
        int nshow = (meta.nkpt < 5) ? meta.nkpt : 5;
        for (int i = 0; i < nshow; ++i) {
            std::println("    {:3}: ({:.4f}, {:.4f}, {:.4f})",
                         i + 1, kpt_vec[i].x, kpt_vec[i].y, kpt_vec[i].z);
        }
        if (meta.nkpt > 5) {
            std::println("    ... and {} more", meta.nkpt - 5);
        }
    }

    OCCMetadata meta;
    std::vector<kVec> kpt_vec;

private:
    auto idx(int iband, int ikpt, int ispin) const -> std::size_t {
        return (static_cast<std::size_t>(ispin) * meta.nkpt + ikpt)
             * meta.nband + iband;
    }

    std::vector<double> energies_;
    std::vector<double> occupations_;
};


OCC::OCC(const std::string& filename) {
    auto* fp = std::fopen(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    try {
        char line[4096];
        int current_spin = 0;    // 0-based: 0=spin1, 1=spin2
        int current_kpt = -1;
        int bands_this_kpt = 0;
        int nband = 0;
        bool spin_mode = false;

        while (std::fgets(line, sizeof(line), fp)) {
            char* p = line;
            while (*p == ' ' || *p == '\t') ++p;
            if (*p == '\0' || *p == '\n') continue;

            // ==========  SPIN N  ==========
            if (*p == '=') {
                spin_mode = true;

                // Finalize previous spin's last k-point
                if (current_kpt >= 0) {
                    if (nband == 0) {
                        nband = bands_this_kpt;
                    } else if (bands_this_kpt != nband) {
                        throw std::runtime_error(
                            "OCC: inconsistent nband in SPIN "
                            + std::to_string(current_spin + 1)
                        );
                    }
                }

                char* sp = std::strstr(p, "SPIN");
                if (!sp) {
                    throw std::runtime_error("OCC: malformed SPIN header");
                }
                sp += 4;
                while (*sp == ' ') ++sp;
                long spin_id = std::strtol(sp, nullptr, 10);
                if (spin_id < 1 || spin_id > 2) {
                    throw std::runtime_error(
                        "OCC: invalid spin number: " + std::to_string(spin_id)
                    );
                }
                current_spin = static_cast<int>(spin_id) - 1;
                current_kpt = -1;
                bands_this_kpt = 0;
                continue;
            }

            // KPOINTS  N:  x  y  z
            if (std::strncmp(p, "KPOINTS", 7) == 0) {
                // Finalize previous k-point
                if (current_kpt >= 0) {
                    if (nband == 0) {
                        nband = bands_this_kpt;
                    } else if (bands_this_kpt != nband) {
                        throw std::runtime_error(
                            "OCC: inconsistent nband at k-point "
                            + std::to_string(current_kpt + 2)
                        );
                    }
                }

                p += 7;
                while (*p == ' ' || *p == '\t') ++p;

                char* end;
                /*long kpt_id =*/ std::strtol(p, &end, 10);
                if (end == p) {
                    throw std::runtime_error("OCC: failed to parse k-point index");
                }
                p = end;
                while (*p == ' ' || *p == '\t' || *p == ':') ++p;

                double kx = std::strtod(p, &end); p = end;
                double ky = std::strtod(p, &end); p = end;
                double kz = std::strtod(p, &end);

                ++current_kpt;
                if (current_spin == 0) {
                    kpt_vec.emplace_back(kx, ky, kz);
                } else {
                    // Verify second spin shares the same k-mesh
                    if (current_kpt >= static_cast<int>(kpt_vec.size())
                        || std::abs(kpt_vec[current_kpt].x - kx) > 1e-4
                        || std::abs(kpt_vec[current_kpt].y - ky) > 1e-4
                        || std::abs(kpt_vec[current_kpt].z - kz) > 1e-4)
                    {
                        throw std::runtime_error(
                            "OCC: SPIN 2 k-point " + std::to_string(current_kpt + 1)
                            + " does not match SPIN 1"
                        );
                    }
                }
                bands_this_kpt = 0;
                continue;
            }

            // Data line:  iband  energy  occupation
            if (*p >= '0' && *p <= '9') {
                if (current_kpt < 0) {
                    throw std::runtime_error("OCC: data before any KPOINTS header");
                }

                char* end;
                /*long iband =*/ std::strtol(p, &end, 10);
                if (end == p) continue;

                double energy = std::strtod(end, &end);
                double occ_val = std::strtod(end, nullptr);

                ++bands_this_kpt;
                energies_.push_back(energy);
                occupations_.push_back(occ_val);
            }
        }

        // Finalize last k-point
        if (current_kpt >= 0) {
            if (nband == 0) {
                nband = bands_this_kpt;
            } else if (bands_this_kpt != nband) {
                throw std::runtime_error("OCC: inconsistent nband at last k-point");
            }
        }

        meta.islda = spin_mode ? (current_spin + 1) : 1;
        meta.nkpt = static_cast<int>(kpt_vec.size());
        meta.nband = nband;

        if (meta.nkpt == 0) throw std::runtime_error("OCC: no k-points found");
        if (meta.nband == 0) throw std::runtime_error("OCC: no band data found");

        auto expected = static_cast<std::size_t>(meta.islda) * meta.nkpt * meta.nband;
        if (energies_.size() != expected || occupations_.size() != expected) {
            throw std::runtime_error(
                "OCC: data size mismatch: expected " + std::to_string(expected)
                + " entries, got " + std::to_string(energies_.size())
            );
        }
    } catch (...) {
        std::fclose(fp);
        throw;
    }

    std::fclose(fp);
}
