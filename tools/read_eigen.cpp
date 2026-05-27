// Standalone OUT.EIGEN reader — no module/C++23 dependency.
// Compile: c++ -std=c++17 -o read_eigen read_eigen.cpp
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr double HARTREE_TO_EV = 27.211386245988;  // CODATA 2018

struct EIGENMeta {
    int islda, nkpt, nband, nref_tot_8, natom, nnode, is_SO;
};

struct KPoint {
    double kx, ky, kz, weight;
};

static auto readI32(std::FILE* fp, const char* ctx) -> int {
    std::int32_t v;
    if (std::fread(&v, sizeof(v), 1, fp) != 1)
        throw std::runtime_error(std::string("read ") + ctx + " failed");
    return static_cast<int>(v);
}

static auto readRecord(std::FILE* fp, void* dst, int nbytes, const char* ctx) -> void {
    int len = readI32(fp, ctx);
    if (len != nbytes)
        throw std::runtime_error(std::string(ctx) + ": size mismatch (" +
                                 std::to_string(len) + " != " + std::to_string(nbytes) + ")");
    if (std::fread(dst, 1, static_cast<std::size_t>(nbytes), fp) != static_cast<std::size_t>(nbytes))
        throw std::runtime_error(std::string("read ") + ctx + " data failed");
    int trail = readI32(fp, ctx);
    if (trail != len)
        throw std::runtime_error(std::string(ctx) + ": trailer mismatch");
}

auto main(int argc, char* argv[]) -> int {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <OUT.EIGEN>\n";
        return 1;
    }

    std::FILE* fp = std::fopen(argv[1], "rb");
    if (!fp) {
        std::cerr << "Cannot open: " << argv[1] << "\n";
        return 1;
    }

    try {
        // --- header ---
        int hdr_len = readI32(fp, "header");
        EIGENMeta m{};
        if (hdr_len == static_cast<int>(7 * sizeof(int))) {
            int hdr[7];
            if (std::fread(hdr, sizeof(int), 7, fp) != 7)
                throw std::runtime_error("read header data failed");
            m.islda = hdr[0]; m.nkpt = hdr[1]; m.nband = hdr[2];
            m.nref_tot_8 = hdr[3]; m.natom = hdr[4]; m.nnode = hdr[5];
            m.is_SO = hdr[6];
        } else if (hdr_len == static_cast<int>(6 * sizeof(int))) {
            int hdr[6];
            if (std::fread(hdr, sizeof(int), 6, fp) != 6)
                throw std::runtime_error("read header data failed");
            m.islda = hdr[0]; m.nkpt = hdr[1]; m.nband = hdr[2];
            m.nref_tot_8 = hdr[3]; m.natom = hdr[4]; m.nnode = hdr[5];
            m.is_SO = 0;
        } else {
            throw std::runtime_error("unexpected header size: " + std::to_string(hdr_len));
        }
        {   // trailing marker
            int t = readI32(fp, "header");
            if (t != hdr_len) throw std::runtime_error("header trailer mismatch");
        }

        // --- read k-points + eigenvalues ---
        std::vector<KPoint> kpts(static_cast<std::size_t>(m.nkpt));
        // eigenvalues[ispin][ikpt * nband + iband]
        std::vector<std::vector<double>> ev(m.islda);
        for (int s = 0; s < m.islda; ++s)
            ev[s].resize(static_cast<std::size_t>(m.nkpt) * m.nband);

        for (int s = 0; s < m.islda; ++s) {
            for (int ik = 0; ik < m.nkpt; ++ik) {
                // KptRecord: 2 int + 4 double = 40 bytes
                struct { int islda_tmp, ikpt_tmp; double w, kx, ky, kz; } rec;
                readRecord(fp, &rec, static_cast<int>(sizeof(rec)), "kpt");
                if (rec.islda_tmp != s + 1 || rec.ikpt_tmp != ik + 1)
                    throw std::runtime_error("kpt index mismatch");

                auto& kp = kpts[static_cast<std::size_t>(ik)];
                kp.weight = rec.w;
                kp.kx = rec.kx;
                kp.ky = rec.ky;
                kp.kz = rec.kz;

                int nbytes = m.nband * static_cast<int>(sizeof(double));
                readRecord(fp, ev[s].data() + static_cast<std::size_t>(ik) * m.nband,
                           nbytes, "ev");
            }
        }
        std::fclose(fp);
        fp = nullptr;

        // --- print info ---
        double sum_w = 0;
        for (const auto& k : kpts) sum_w += k.weight;
        std::cout << "File: " << argv[1] << "\n"
                  << "  islda=" << m.islda << "  nkpt=" << m.nkpt
                  << "  nband=" << m.nband << "  nref_tot_8=" << m.nref_tot_8
                  << "  natom=" << m.natom << "  nnode=" << m.nnode
                  << "  is_SO=" << m.is_SO << "\n"
                  << "  sum(weight)=" << sum_w << "\n\n"
                  << "  ikpt  weight        kx             ky             kz\n";
        for (int ik = 0; ik < m.nkpt; ++ik) {
            auto& k = kpts[static_cast<std::size_t>(ik)];
            std::printf("  %4d  %8.6f  %14.12f  %14.12f  %14.12f\n",
                        ik, k.weight, k.kx, k.ky, k.kz);
        }

        // --- interactive query ---
        while (true) {
            std::cout << "\nikpt [0," << m.nkpt << ")  iband [0," << m.nband
                      << ")  [ispin]  (-1 to quit): ";
            int ik, ib, is = 0;
            std::cin >> ik;
            if (ik < 0) break;
            std::cin >> ib;

            // allow optional ispin
            if (std::cin.peek() != '\n') std::cin >> is;

            if (ik < 0 || ik >= m.nkpt || ib < 0 || ib >= m.nband ||
                is < 0 || is >= m.islda) {
                std::cerr << "Invalid indices.\n";
                continue;
            }

            double ev_h = ev[static_cast<std::size_t>(is)]
                            [static_cast<std::size_t>(ik) * m.nband + ib];
            double ev_ev = ev_h / HARTREE_TO_EV;
            std::printf("  eig[%d,%d,%d] = %14.8f Hartree = %12.8f eV\n",
                        ib, ik, is, ev_h, ev_ev);
        }
        return 0;
    } catch (const std::exception& e) {
        if (fp) std::fclose(fp);
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
