#pragma once

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>


// Lightweight header-only reader for OUT.OCC text format.
// Shared by tests — if an io.OCC module is added later, this can serve as
// a reference or be replaced by importing that module.


struct kPoint {
    double x{}, y{}, z{};
};


struct OCCData {
    int nkpt{};
    int nband{};  // bands per k-point (inferred from file)

    std::vector<kPoint> kpt_vec;
    // Flat: [ikpt * nband + iband]
    std::vector<double> energies;      // eV
    std::vector<double> occupations;

    auto energy(int ikpt, int iband) const -> double {
        return energies[idx(ikpt, iband)];
    }

    auto occupation(int ikpt, int iband) const -> double {
        return occupations[idx(ikpt, iband)];
    }

private:
    auto idx(int ikpt, int iband) const -> std::size_t {
        return static_cast<std::size_t>(ikpt) * nband + iband;
    }
};


inline auto parseOCC(const std::string& filename) -> OCCData {
    auto* fp = std::fopen(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    OCCData occ;

    try {
        char line[256];
        int current_kpt = -1;
        int bands_this_kpt = 0;

        while (std::fgets(line, sizeof(line), fp)) {
            char* p = line;
            while (*p == ' ' || *p == '\t') ++p;

            if (std::strncmp(p, "KPOINTS", 7) == 0) {
                // Finalize previous k-point's band count
                if (current_kpt >= 0) {
                    if (occ.nband == 0) {
                        occ.nband = bands_this_kpt;
                    } else if (bands_this_kpt != occ.nband) {
                        throw std::runtime_error(
                            "OCC: inconsistent nband at k-point "
                            + std::to_string(current_kpt + 1)
                        );
                    }
                }

                // Parse: KPOINTS  N:  x  y  z
                p += 7;
                while (*p == ' ' || *p == '\t') ++p;

                char* end;
                long kpt_id = std::strtol(p, &end, 10);
                if (end == p || kpt_id < 1) {
                    throw std::runtime_error("OCC: failed to parse k-point index");
                }
                p = end;
                while (*p == ' ' || *p == '\t' || *p == ':') ++p;

                double kx = std::strtod(p, &end); p = end;
                double ky = std::strtod(p, &end); p = end;
                double kz = std::strtod(p, &end);

                ++current_kpt;
                occ.kpt_vec.push_back({kx, ky, kz});
                bands_this_kpt = 0;
            }
            else if (*p >= '0' && *p <= '9') {
                if (current_kpt < 0) {
                    throw std::runtime_error("OCC: data before any KPOINTS header");
                }

                char* end;
                /*long iband =*/ std::strtol(p, &end, 10);
                if (end == p) continue;

                double energy = std::strtod(end, &end);
                double occ_val = std::strtod(end, nullptr);

                ++bands_this_kpt;
                occ.energies.push_back(energy);
                occ.occupations.push_back(occ_val);
            }
            // else: header line or blank → skip
        }

        // Finalize last k-point
        if (current_kpt >= 0) {
            if (occ.nband == 0) {
                occ.nband = bands_this_kpt;
            } else if (bands_this_kpt != occ.nband) {
                throw std::runtime_error("OCC: inconsistent nband at last k-point");
            }
        }

        occ.nkpt = current_kpt + 1;
        if (occ.nkpt == 0) throw std::runtime_error("OCC: no k-points found");
        if (occ.nband == 0) throw std::runtime_error("OCC: no band data found");

        if (static_cast<std::size_t>(occ.nkpt * occ.nband) != occ.energies.size()) {
            throw std::runtime_error(
                "OCC: data size mismatch: expected "
                + std::to_string(occ.nkpt * occ.nband)
                + " entries, got " + std::to_string(occ.energies.size())
            );
        }
    } catch (...) {
        std::fclose(fp);
        throw;
    }

    std::fclose(fp);
    return occ;
}
