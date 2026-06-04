module H_psi.hamiltonian;

import std;
import utils.logger;

// ===========================================================================
//  Constructor
// ===========================================================================

Hamiltonian::Hamiltonian(std::filesystem::path base_dir) : base_dir_(std::move(base_dir)) {
    if (!std::filesystem::is_directory(base_dir_)) {
        throw std::runtime_error(
            "Hamiltonian: base directory does not exist or is not a directory: "
            + base_dir_.string());
    }
    base_dir_ = std::filesystem::absolute(base_dir_);
}


// ===========================================================================
//  Path resolution — relative paths are prefixed with base_dir_
// ===========================================================================

auto Hamiltonian::resolve(const std::string& path) const -> std::string {
    auto p = std::filesystem::path(path);
    if (p.is_absolute()) return path;
    return (base_dir_ / p).string();
}

namespace {

// Compare two 3×3 lattice matrices with tolerance.
auto latticesMatch(const Lattice& a, const Lattice& b) -> bool {
    auto A_a = a.A();  // Hartree atomic units (Bohr)
    auto A_b = b.A();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (std::abs(A_a[i][j] - A_b[i][j]) > 1e-10) return false;
    return true;
}

} // anonymous namespace


// ===========================================================================
//  Step-by-step loading
// ===========================================================================

auto Hamiltonian::loadGKK(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading GKK: {}", full);
    gkk_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadWG(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading WG: {}", full);
    wg_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadVR(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading VR: {}", full);
    vr_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadRHO(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading RHO: {}", full);
    rho_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadATOM(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading ATOM: {}", full);
    atom_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadEIGEN(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading EIGEN: {}", full);
    eigen_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadOCC(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading OCC: {}", full);
    occ_.emplace(full);
    checkConsistency();
}

auto Hamiltonian::loadNCPPs(const std::string& directory) -> void {
    auto full = resolve(directory);
    std::println("[Hamiltonian] loading NCPPs from: {}", full);
    auto upf_dir = std::filesystem::path(full);
    if (!std::filesystem::is_directory(upf_dir)) {
        throw std::runtime_error(
            "Hamiltonian: directory not found: " + upf_dir.string());
    }

    for (const auto& entry : std::filesystem::directory_iterator(upf_dir)) {
        auto ext = entry.path().extension();
        if (ext == ".UPF" || ext == ".upf") {
            ncpps_.emplace_back(UPF(entry.path().string()));
            std::println("  loaded: {}", entry.path().filename().string());
        }
    }
    std::println("  total NCPP objects: {}", ncpps_.size());
    checkConsistency();
}


// ===========================================================================
//  Convenience
// ===========================================================================

auto Hamiltonian::loadFromDirectory() -> void {
    // Order: load the ones with fewer cross-file dependencies first so that
    // each checkConsistency() call gets to verify progressively more pairs.
    loadATOM("atom.config");

    // Try UPF/ subdirectory first, fall back to the base directory itself.
    auto upf_subdir = base_dir_ / "UPF";
    if (std::filesystem::is_directory(upf_subdir)) {
        loadNCPPs("UPF");
    } else {
        loadNCPPs(".");
    }

    loadGKK("OUT.GKK");
    loadWG("OUT.WG");
    loadVR("OUT.VR");
    loadRHO("OUT.RHO");
    loadEIGEN("OUT.EIGEN");
    loadOCC("OUT.OCC");
}


// ===========================================================================
//  Data accessors
// ===========================================================================

auto Hamiltonian::gkk() const -> const GKK& {
    if (!gkk_) throw std::runtime_error("Hamiltonian: GKK not loaded");
    return *gkk_;
}

auto Hamiltonian::wg() const -> const WG& {
    if (!wg_) throw std::runtime_error("Hamiltonian: WG not loaded");
    return *wg_;
}

auto Hamiltonian::vr() const -> const VR& {
    if (!vr_) throw std::runtime_error("Hamiltonian: VR not loaded");
    return *vr_;
}

auto Hamiltonian::rho() const -> const RHO& {
    if (!rho_) throw std::runtime_error("Hamiltonian: RHO not loaded");
    return *rho_;
}

auto Hamiltonian::atom() const -> const ATOM& {
    if (!atom_) throw std::runtime_error("Hamiltonian: ATOM not loaded");
    return *atom_;
}

auto Hamiltonian::eigen() const -> const EIGEN& {
    if (!eigen_) throw std::runtime_error("Hamiltonian: EIGEN not loaded");
    return *eigen_;
}

auto Hamiltonian::occ() const -> const OCC& {
    if (!occ_) throw std::runtime_error("Hamiltonian: OCC not loaded");
    return *occ_;
}

auto Hamiltonian::ncpp(int atomic_number) const -> const NCPP& {
    std::string_view want = ATOM::elementName(atomic_number);
    for (const auto& p : ncpps_) {
        if (p.meta.element == want) return p;
    }
    throw std::out_of_range(
        "Hamiltonian: no NCPP for element " + std::string(want) +
        " (Z=" + std::to_string(atomic_number) + ')');
}


// ===========================================================================
//  Consistency check  —  runs after every load step
// ===========================================================================

auto Hamiltonian::checkConsistency() -> void {
    using namespace std;
    auto& log = Logger::instance();

    // -- helper: log a debug-level skip message ---------------------------------
    auto debug = [&](string_view msg) {
        log.log(LogLevel::Debug, "  — {}: skipped (data not loaded)", msg);
    };

    // -- helper: check int equality, log at info-level --------------------------
    auto checkInt = [&](int a, int b, string_view label) -> void {
        bool ok = (a == b);
        log.log(LogLevel::Info, "  {} {}: {} vs {}", ok ? "✓" : "✗", label, a, b);
        if (!ok) throw runtime_error(
            format("Hamiltonian consistency: {} mismatch ({} vs {})",
                   label, a, b));
    };

    // -- helper: check lattice match -------------------------------------------
    auto checkLattice = [&](const Lattice& a, const Lattice& b,
                            string_view label) -> void {
        bool ok = latticesMatch(a, b);
        log.log(LogLevel::Info, "  {} lattice [{}]", ok ? "✓" : "✗", label);
        if (!ok) throw runtime_error(
            format("Hamiltonian consistency: lattice mismatch [{}]", label));
    };

    log.log(LogLevel::Info, "[Hamiltonian] consistency check");

    // ----  pair: GKK  x  WG  ---------------------------------------------------
    if (gkk_ && wg_ && !checked_[Pair_GKK_WG]) {
        const auto& g = gkk_->meta;
        const auto& w = wg_->meta;
        checkInt(g.nkpt,  w.nkpt,  "nkpt  [GKK vs WG]");
        checkInt(g.n1,    w.n1,    "n1    [GKK vs WG]");
        checkInt(g.n2,    w.n2,    "n2    [GKK vs WG]");
        checkInt(g.n3,    w.n3,    "n3    [GKK vs WG]");
        checkInt(g.mg_nx, w.mg_nx, "mg_nx [GKK vs WG]");
        checkInt(g.nnode, w.nnode, "nnode [GKK vs WG]");
        checkInt(g.is_SO, w.is_SO, "is_SO [GKK vs WG]");
        checkInt(g.islda, w.islda, "islda [GKK vs WG]");

        bool ecut_ok = (abs(g.Ecut - w.Ecut) < 1e-12);
        log.log(LogLevel::Info, "  {} Ecut  [GKK vs WG]: {} vs {}",
                ecut_ok ? "✓" : "✗", g.Ecut, w.Ecut);
        if (!ecut_ok) throw runtime_error(
            "Hamiltonian consistency: Ecut mismatch");

        bool ng_ok = (g.ng_tot_per_kpt == w.ng_tot_per_kpt);
        log.log(LogLevel::Info, "  {} ng_tot_per_kpt [GKK vs WG]", ng_ok ? "✓" : "✗");
        if (!ng_ok) throw runtime_error(
            "Hamiltonian consistency: ng_tot_per_kpt mismatch");

        checkLattice(g.lattice, w.lattice, "GKK vs WG");
        checked_.set(Pair_GKK_WG);
    } else {
        debug("GKK vs WG");
    }

    // ----  pair: GKK  x  VR  ---------------------------------------------------
    if (gkk_ && vr_ && !checked_[Pair_GKK_VR]) {
        const auto& g = gkk_->meta;
        const auto& v = vr_->meta;
        checkInt(g.n1,     v.n1,     "n1    [GKK vs VR]");
        checkInt(g.n2,     v.n2,     "n2    [GKK vs VR]");
        checkInt(g.n3,     v.n3,     "n3    [GKK vs VR]");
        checkLattice(gkk_->meta.lattice, vr_->lattice, "GKK vs VR");
        checked_.set(Pair_GKK_VR);
    } else {
        debug("GKK vs VR");
    }

    // ----  pair: GKK  x  ATOM  ------------------------------------------------
    if (gkk_ && atom_ && !checked_[Pair_GKK_ATOM]) {
        checkLattice(gkk_->meta.lattice, atom_->lattice, "GKK vs ATOM");
        checked_.set(Pair_GKK_ATOM);
    } else {
        debug("GKK vs ATOM");
    }

    // ----  pair: GKK  x  EIGEN  -----------------------------------------------
    if (gkk_ && eigen_ && !checked_[Pair_GKK_EIGEN]) {
        const auto& g = gkk_->meta;
        const auto& e = eigen_->meta;
        checkInt(g.nkpt,  e.nkpt,  "nkpt  [GKK vs EIGEN]");
        checkInt(g.is_SO, e.is_SO, "is_SO [GKK vs EIGEN]");
        checkInt(g.islda, e.islda, "islda [GKK vs EIGEN]");
        checkInt(g.nnode, e.nnode, "nnode [GKK vs EIGEN]");

        // Compare k-point vectors, allowing periodic translation
        {
            auto saved_view = gkk_->currentView();
            gkk_->setDataView(saved_view | KVecsView::Integer);
            bool kpts_ok = true;
            int bad_ik = -1;
            double bad_gkk_x{}, bad_gkk_y{}, bad_gkk_z{};
            double bad_eig_x{}, bad_eig_y{}, bad_eig_z{};
            constexpr double TWO_PI = 2.0 * std::numbers::pi;
            auto A = g.lattice.A();  // direct lattice (Bohr)
            for (int ik = 0; ik < g.nkpt; ++ik) {
                const auto& kv = gkk_->loadKPoint(ik);
                const auto& eig_k = eigen_->kpt_vec[ik];
                // EIGEN stores k-points in Cartesian (Bohr^-1).
                // GKK infers fractional coordinates.  Convert EIGEN → fractional
                // via  kf[i] = (A[i] · k_cart) / (2π).
                double eig_kf_x = (A[0][0] * eig_k.x + A[0][1] * eig_k.y + A[0][2] * eig_k.z) / TWO_PI;
                double eig_kf_y = (A[1][0] * eig_k.x + A[1][1] * eig_k.y + A[1][2] * eig_k.z) / TWO_PI;
                double eig_kf_z = (A[2][0] * eig_k.x + A[2][1] * eig_k.y + A[2][2] * eig_k.z) / TWO_PI;
                double dx = kv.kPoint.x - eig_kf_x;
                double dy = kv.kPoint.y - eig_kf_y;
                double dz = kv.kPoint.z - eig_kf_z;
                constexpr double eps = 1e-8;
                if (abs(dx - round(dx)) > eps ||
                    abs(dy - round(dy)) > eps ||
                    abs(dz - round(dz)) > eps) {
                    kpts_ok = false;
                    bad_ik = ik;
                    bad_gkk_x = kv.kPoint.x; bad_gkk_y = kv.kPoint.y; bad_gkk_z = kv.kPoint.z;
                    bad_eig_x = eig_k.x;      bad_eig_y = eig_k.y;      bad_eig_z = eig_k.z;
                    break;
                }
            }
            gkk_->setDataView(saved_view);
            if (!kpts_ok) throw runtime_error(format(
                "Hamiltonian consistency: k-point[{}] mismatch (GKK vs EIGEN): "
                "GKK_frac({},{},{}) EIGEN_cart({},{},{})",
                bad_ik,
                bad_gkk_x, bad_gkk_y, bad_gkk_z,
                bad_eig_x, bad_eig_y, bad_eig_z));
            log.log(LogLevel::Info, "  ✓ k-point vectors [GKK vs EIGEN]");
        }

        checked_.set(Pair_GKK_EIGEN);
    } else {
        debug("GKK vs EIGEN");
    }

    // ----  pair: EIGEN  x  ATOM  ---------------------------------------------
    if (eigen_ && atom_ && !checked_[Pair_EIGEN_ATOM]) {
        checkInt(eigen_->meta.natom, atom_->natom, "natom [EIGEN vs ATOM]");
        checked_.set(Pair_EIGEN_ATOM);
    } else {
        debug("EIGEN vs ATOM");
    }

    // ----  pair: WG  x  EIGEN  ------------------------------------------------
    if (wg_ && eigen_ && !checked_[Pair_WG_EIGEN]) {
        checkInt(wg_->meta.nband, eigen_->meta.nband, "nband [WG vs EIGEN]");
        checkInt(wg_->meta.nkpt,  eigen_->meta.nkpt,  "nkpt  [WG vs EIGEN]");
        checkInt(wg_->meta.is_SO, eigen_->meta.is_SO, "is_SO [WG vs EIGEN]");
        checkInt(wg_->meta.islda, eigen_->meta.islda, "islda [WG vs EIGEN]");
        checked_.set(Pair_WG_EIGEN);
    } else {
        debug("WG vs EIGEN");
    }

    // ----  pair: VR  x  ATOM  -------------------------------------------------
    if (vr_ && atom_ && !checked_[Pair_VR_ATOM]) {
        checkLattice(vr_->lattice, atom_->lattice, "VR vs ATOM");
        checked_.set(Pair_VR_ATOM);
    } else {
        debug("VR vs ATOM");
    }

    // ----  pair: VR  x  WG  ---------------------------------------------------
    if (vr_ && wg_ && !checked_[Pair_VR_WG]) {
        checkInt(vr_->meta.n1, wg_->meta.n1, "n1 [VR vs WG]");
        checkInt(vr_->meta.n2, wg_->meta.n2, "n2 [VR vs WG]");
        checkInt(vr_->meta.n3, wg_->meta.n3, "n3 [VR vs WG]");
        checked_.set(Pair_VR_WG);
    } else {
        debug("VR vs WG");
    }

    // ----  pair: RHO  x  VR  -------------------------------------------------
    if (rho_ && vr_ && !checked_[Pair_RHO_VR]) {
        const auto& r = rho_->meta;
        const auto& v = vr_->meta;
        checkInt(r.n1,     v.n1,     "n1    [RHO vs VR]");
        checkInt(r.n2,     v.n2,     "n2    [RHO vs VR]");
        checkInt(r.n3,     v.n3,     "n3    [RHO vs VR]");
        checkInt(r.nstate, v.nstate, "nstate[RHO vs VR]");
        checkLattice(rho_->lattice, vr_->lattice, "RHO vs VR");
        checked_.set(Pair_RHO_VR);
    } else {
        debug("RHO vs VR");
    }

    // ----  pair: RHO  x  ATOM  -----------------------------------------------
    if (rho_ && atom_ && !checked_[Pair_RHO_ATOM]) {
        checkLattice(rho_->lattice, atom_->lattice, "RHO vs ATOM");
        checked_.set(Pair_RHO_ATOM);
    } else {
        debug("RHO vs ATOM");
    }

    // ----  pair: RHO  x  GKK  ------------------------------------------------
    if (rho_ && gkk_ && !checked_[Pair_RHO_GKK]) {
        const auto& r = rho_->meta;
        const auto& g = gkk_->meta;
        checkInt(r.n1,     g.n1,     "n1    [RHO vs GKK]");
        checkInt(r.n2,     g.n2,     "n2    [RHO vs GKK]");
        checkInt(r.n3,     g.n3,     "n3    [RHO vs GKK]");
        checkInt(r.nnode, g.nnode, "nnode [RHO vs GKK]");
        checkLattice(rho_->lattice, g.lattice, "RHO vs GKK");
        checked_.set(Pair_RHO_GKK);
    } else {
        debug("RHO vs GKK");
    }

    // ----  pair: RHO  x  WG  -------------------------------------------------
    if (rho_ && wg_ && !checked_[Pair_RHO_WG]) {
        checkInt(rho_->meta.n1, wg_->meta.n1, "n1 [RHO vs WG]");
        checkInt(rho_->meta.n2, wg_->meta.n2, "n2 [RHO vs WG]");
        checkInt(rho_->meta.n3, wg_->meta.n3, "n3 [RHO vs WG]");
        checked_.set(Pair_RHO_WG);
    } else {
        debug("RHO vs WG");
    }

    // ----  NCPPs  x  ATOM  ----------------------------------------------------
    if (!ncpps_.empty() && atom_ && !checked_[Pair_NCPP_ATOM]) {
        // Only check that every ATOM species has a matching NCPP.
        // Extra NCPPs are allowed (over-loading).
        for (int it = 0; it < atom_->ntyp; ++it) {
            int         z        = atom_->zvals[it];
            string_view expected = ATOM::elementName(z);
            bool found = false;
            for (const auto& p : ncpps_) {
                if (p.meta.element == expected) { found = true; break; }
            }
            if (found) {
                log.log(LogLevel::Info, "  ✓ element {} [NCPP vs ATOM]", expected);
            } else {
                log.log(LogLevel::Info, "  ✗ element {} [NCPP vs ATOM]: no matching UPF", expected);
                throw runtime_error(
                    format("Hamiltonian consistency: no NCPP for element {} (Z={})",
                           expected, z));
            }
        }
        checked_.set(Pair_NCPP_ATOM);
    } else {
        debug("NCPP vs ATOM");
    }

    // ----  pair: EIGEN  x  VR  -------------------------------------------------
    if (eigen_ && vr_ && !checked_[Pair_EIGEN_VR]) {
        checkInt(eigen_->meta.nnode, vr_->meta.nnode, "nnode [EIGEN vs VR]");
        if (atom_) {
            checkLattice(vr_->lattice, atom_->lattice, "EIGEN vs VR");
        }
        checked_.set(Pair_EIGEN_VR);
    } else {
        debug("EIGEN vs VR");
    }

    // ----  pair: EIGEN  x  RHO  ------------------------------------------------
    if (eigen_ && rho_ && !checked_[Pair_EIGEN_RHO]) {
        checkInt(eigen_->meta.nnode, rho_->meta.nnode, "nnode [EIGEN vs RHO]");
        if (atom_) {
            checkLattice(rho_->lattice, atom_->lattice, "EIGEN vs RHO");
        }
        checked_.set(Pair_EIGEN_RHO);
    } else {
        debug("EIGEN vs RHO");
    }

    // ----  pair: OCC  x  EIGEN  ------------------------------------------------
    if (occ_ && eigen_ && !checked_[Pair_OCC_EIGEN]) {
        const auto& o = occ_->meta;
        const auto& e = eigen_->meta;
        checkInt(o.islda, e.islda, "islda [OCC vs EIGEN]");
        checkInt(o.nkpt,  e.nkpt,  "nkpt  [OCC vs EIGEN]");
        checkInt(o.nband, e.nband, "nband [OCC vs EIGEN]");

        // Compare k-point vectors directly (OCC text format has ~4 dp precision)
        bool kpts_ok = true;
        int bad_ik = -1;
        constexpr double eps = 1e-4;
        for (int ik = 0; ik < o.nkpt; ++ik) {
            auto dx = std::abs(occ_->kpt_vec[ik].x - eigen_->kpt_vec[ik].x);
            auto dy = std::abs(occ_->kpt_vec[ik].y - eigen_->kpt_vec[ik].y);
            auto dz = std::abs(occ_->kpt_vec[ik].z - eigen_->kpt_vec[ik].z);
            if (dx > eps || dy > eps || dz > eps) {
                kpts_ok = false;
                bad_ik = ik;
                break;
            }
        }
        if (!kpts_ok) throw runtime_error(format(
            "Hamiltonian consistency: k-point[{}] mismatch (OCC vs EIGEN)",
            bad_ik));
        log.log(LogLevel::Info, "  ✓ k-point vectors [OCC vs EIGEN]");

        checked_.set(Pair_OCC_EIGEN);
    } else {
        debug("OCC vs EIGEN");
    }

    // ----  pair: OCC  x  WG  --------------------------------------------------
    if (occ_ && wg_ && !checked_[Pair_OCC_WG]) {
        checkInt(occ_->meta.nkpt,  wg_->meta.nkpt,  "nkpt  [OCC vs WG]");
        checkInt(occ_->meta.nband, wg_->meta.nband, "nband [OCC vs WG]");
        checked_.set(Pair_OCC_WG);
    } else {
        debug("OCC vs WG");
    }

    // ----  pair: OCC  x  GKK  -------------------------------------------------
    if (occ_ && gkk_ && !checked_[Pair_OCC_GKK]) {
        checkInt(occ_->meta.nkpt, gkk_->meta.nkpt, "nkpt [OCC vs GKK]");
        checked_.set(Pair_OCC_GKK);
    } else {
        debug("OCC vs GKK");
    }

    log.log(LogLevel::Info, "[Hamiltonian] consistency check complete");
}
