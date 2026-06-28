module H_psi.hamiltonian;

import std;
import utils.logger;
import support.density;

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

namespace {
    enum : std::uint64_t {
        PART2_GKK_WG    = 1ull << 0,   // mg_nx, ng_tot_per_kpt, nnode
        PART2_GKK_EIGEN = 1ull << 1,   // nnode
        PART2_VR_RHO    = 1ull << 2,   // nstate
    };
}


// ===========================================================================
//  Step-by-step loading
// ===========================================================================

auto Hamiltonian::loadGKK(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading GKK: {}", full);
    gkk_.emplace(full);
    part2_done_ &= ~(PART2_GKK_WG | PART2_GKK_EIGEN);
    checkConsistency();
}

auto Hamiltonian::loadWG(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading WG: {}", full);
    wg_.emplace(full);
    part2_done_ &= ~PART2_GKK_WG;
    checkConsistency();
}

auto Hamiltonian::loadVR(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading VR: {}", full);
    vr_.emplace(full);
    part2_done_ &= ~PART2_VR_RHO;
    checkConsistency();
}

auto Hamiltonian::loadRHO(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading RHO: {}", full);
    rho_.emplace(full);
    part2_done_ &= ~PART2_VR_RHO;
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
    part2_done_ &= ~PART2_GKK_EIGEN;
    checkConsistency();
}

auto Hamiltonian::loadOCC(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading OCC: {}", full);
    occ_.emplace(full);
    // OCC does not participate in any Part 2 pair check.
    checkConsistency();
}

auto Hamiltonian::loadNCPP(const std::string& path) -> void {
    auto full = resolve(path);
    std::println("[Hamiltonian] loading NCPP: {}", full);

    auto ncpp = NCPP(UPF(full));
    auto name = ncpp.meta.element;

    elements_.push_back({std::move(ncpp), std::nullopt});
    std::println("  loaded: {} ({} element(s) total)", name, elements_.size());

    checkConsistency();
}


// ===========================================================================
//  Convenience
// ===========================================================================

auto Hamiltonian::loadFromDirectory() -> void {
    // Order: load the ones with fewer cross-file dependencies first so that
    // each checkConsistency() call gets to verify progressively more pairs.
    loadATOM("atom.config");

    // Load each UPF file from UPF/ subdirectory, falling back to base_dir_.
    auto upf_dir = base_dir_ / "UPF";
    if (!std::filesystem::is_directory(upf_dir))
        upf_dir = base_dir_;
    for (const auto& entry : std::filesystem::directory_iterator(upf_dir)) {
        auto ext = entry.path().extension();
        if (ext == ".UPF" || ext == ".upf")
            loadNCPP(std::filesystem::relative(entry.path(), base_dir_).string());
    }

    loadGKK("OUT.GKK");
    loadWG("OUT.WG");
    loadVR("OUT.VR");
    loadRHO("OUT.RHO");
    loadEIGEN("OUT.EIGEN");
    loadOCC("OUT.OCC");
    finalize();
}


// ===========================================================================
//  Data accessors
// ===========================================================================

auto Hamiltonian::gkk() const -> GKK& {
    if (!gkk_) throw std::runtime_error("Hamiltonian: GKK not loaded");
    return *gkk_;
}

auto Hamiltonian::wg() const -> WG& {
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



// -------------------------------------------------------------------------
//  finalize  —  post-load initialisation
// -------------------------------------------------------------------------

auto Hamiltonian::finalize(std::initializer_list<ExtendedCheck> checks,
                           std::uint64_t psp_features) -> void {
    if (!canonical_lattice_) throw std::runtime_error("finalize: no lattice loaded");
    if (elements_.empty()) return;

    // -------------------------------------------------------------------
    //  PSP features  —  validate and store
    // -------------------------------------------------------------------

    constexpr auto NONLOCAL = static_cast<std::uint64_t>(PSPFeature::Nonlocal);
    constexpr auto DBETAQ   = static_cast<std::uint64_t>(PSPFeature::DBetaq);
    constexpr auto D2BETAQ  = static_cast<std::uint64_t>(PSPFeature::D2Betaq);

    if ((psp_features & D2BETAQ) && !(psp_features & DBETAQ)) {
        throw std::invalid_argument(
            "finalize: PSPFeature::D2Betaq requires PSPFeature::DBetaq");
    }
    if ((psp_features & DBETAQ) && !(psp_features & NONLOCAL)) {
        throw std::invalid_argument(
            "finalize: PSPFeature::DBetaq requires PSPFeature::Nonlocal");
    }
    psp_features_ = psp_features;

    // -------------------------------------------------------------------
    //  BetaqTables / DBetaqTables / D2BetaqTables  —  NCPP initialisation
    // -------------------------------------------------------------------

    double omega = 0.0;
    if (psp_features & NONLOCAL) {
        std::println("[Hamiltonian] finalize: constructing BetaqTables...");
        omega = canonical_lattice_->volume();

        auto buildTables = [&](ElementData& elem) -> void {
            sortByL(elem.ncpp);
            diagonalizeNonlocal(elem.ncpp);

            const auto& ncpp = elem.ncpp;
            auto mesh_type = (ncpp.mesh.type == MeshType::Uniform)
                             ? SimpsonMeshType::Uniform : SimpsonMeshType::General;

            // BetaqTables — always built when Nonlocal is set
            elem.betaq_tables.emplace(
                ncpp.mesh.r, ncpp.mesh.rab, mesh_type,
                ncpp.nonlocal.beta, ncpp.nonlocal.angular_momentum);
            elem.betaq_tables->setVolume(omega);

            // DBetaqTables — built only when DBetaq flag is set
            if (psp_features & DBETAQ) {
                elem.dbetaq_tables.emplace(
                    ncpp.mesh.r, ncpp.mesh.rab, mesh_type,
                    ncpp.nonlocal.beta, ncpp.nonlocal.angular_momentum);
                elem.dbetaq_tables->setVolume(omega);
            }

            // D2BetaqTables — built only when D2Betaq flag is set
            if (psp_features & D2BETAQ) {
                elem.d2betaq_tables.emplace(
                    ncpp.mesh.r, ncpp.mesh.rab, mesh_type,
                    ncpp.nonlocal.beta, ncpp.nonlocal.angular_momentum);
                elem.d2betaq_tables->setVolume(omega);
            }
        };

        // Build tables only for elements that actually appear in the
        // crystal (atom.config).  When no ATOM file is loaded, conservatively
        // build for every loaded NCPP instead.
        if (atom_) {
            for (auto&& t : atom_->eachType()) {
                std::string_view want = ATOM::elementName(t.z);
                auto it = std::ranges::find_if(elements_, [want](const auto& e) {
                    return e.ncpp.meta.element == want;
                });
                if (it != elements_.end())
                    buildTables(*it);
            }
        } else {
            for (auto& elem : elements_)
                buildTables(elem);
        }
        if (psp_features & DBETAQ) {
            std::println("[Hamiltonian] finalize: DBetaqTables constructed");
        }
        if (psp_features & D2BETAQ) {
            std::println("[Hamiltonian] finalize: D2BetaqTables constructed");
        }
    }

    // -------------------------------------------------------------------
    //  Derived quantities
    // -------------------------------------------------------------------

    // l_max_ — maximum angular momentum across all elements.
    l_max_ = 0;
    for (const auto& elem : elements_) {
        if (elem.ncpp.meta.l_max > l_max_) l_max_ = elem.ncpp.meta.l_max;
    }

    // ng_max_ — maximum G-vectors across k-points (from GKK).
    ng_max_ = 0;
    if (gkk_) {
        const auto& per_kpt = gkk_->meta.ng_tot_per_kpt;
        if (!per_kpt.empty()) {
            ng_max_ = *std::ranges::max_element(per_kpt);
        }
    }

    // -------------------------------------------------------------------
    //  FFT plan  —  construct from canonical FFT grid dimensions
    // -------------------------------------------------------------------

    if (canonical_fft_grid_) {
        fft_.emplace(canonical_fft_grid_->n1, canonical_fft_grid_->n2, canonical_fft_grid_->n3);
    }

    // -------------------------------------------------------------------
    //  Consistency checks
    // -------------------------------------------------------------------

    checkConsistencyExtended(checks);

    // -------------------------------------------------------------------
    //  Done
    // -------------------------------------------------------------------

    std::println("  done ({} element(s), omega={:.6f}, ng_max={}, l_max={})",
                 elements_.size(), omega, ng_max_, l_max_);
}

auto Hamiltonian::ncpp(int atomic_number) const -> const NCPP& {
    std::string_view want = ATOM::elementName(atomic_number);
    for (const auto& elem : elements_) {
        if (elem.ncpp.meta.element == want) return elem.ncpp;
    }
    throw std::out_of_range(
        "Hamiltonian: no NCPP for element " + std::string(want) +
        " (Z=" + std::to_string(atomic_number) + ')');
}

auto Hamiltonian::betaqTables(int atomic_number) const -> const BetaqTables& {
    std::string_view want = ATOM::elementName(atomic_number);
    for (const auto& elem : elements_) {
        if (elem.ncpp.meta.element == want) {
            if (!elem.betaq_tables) {
                throw std::runtime_error(
                    "Hamiltonian: BetaqTables not initialised for element "
                    + std::string(want) + " (did you call finalize()?)");
            }
            return *elem.betaq_tables;
        }
    }
    throw std::out_of_range(
        "Hamiltonian: no BetaqTables for element " + std::string(want) +
        " (Z=" + std::to_string(atomic_number) + ')');
}

auto Hamiltonian::dbetaqTables(int atomic_number) const -> const DBetaqTables& {
    std::string_view want = ATOM::elementName(atomic_number);
    for (const auto& elem : elements_) {
        if (elem.ncpp.meta.element == want) {
            if (!elem.dbetaq_tables) {
                throw std::runtime_error(
                    "Hamiltonian: DBetaqTables not initialised for element "
                    + std::string(want)
                    + " (did you enable PSPFeature::DBetaq in finalize()?)");
            }
            return *elem.dbetaq_tables;
        }
    }
    throw std::out_of_range(
        "Hamiltonian: no DBetaqTables for element " + std::string(want) +
        " (Z=" + std::to_string(atomic_number) + ')');
}

auto Hamiltonian::d2betaqTables(int atomic_number) const -> const D2BetaqTables& {
    std::string_view want = ATOM::elementName(atomic_number);
    for (const auto& elem : elements_) {
        if (elem.ncpp.meta.element == want) {
            if (!elem.d2betaq_tables) {
                throw std::runtime_error(
                    "Hamiltonian: D2BetaqTables not initialised for element "
                    + std::string(want)
                    + " (did you enable PSPFeature::D2Betaq in finalize()?)");
            }
            return *elem.d2betaq_tables;
        }
    }
    throw std::out_of_range(
        "Hamiltonian: no D2BetaqTables for element " + std::string(want) +
        " (Z=" + std::to_string(atomic_number) + ')');
}


// ===========================================================================
//  Consistency check  —  runs after every load step
//  Three-part structure:
//    Part 1 — canonical physical quantities (first-load-sets-canonical)
//    Part 2 — file‑to‑file integrity checks
//    Part 3 — heavyweight computational verification (separate method)
// ===========================================================================

auto Hamiltonian::checkConsistency() -> void {
    Logger::instance().log(LogLevel::Info, "[Hamiltonian] consistency check");
    checkPart1();
    checkPart2();
    Logger::instance().log(LogLevel::Info, "[Hamiltonian] consistency check complete");
}


// ===========================================================================
//  Part 1  —  canonical physical quantities
//  Each quantity follows "first load sets the canonical, subsequent loads
//  must match".  The canonical values are stored in Hamiltonian members
//  and are available for H|ψ⟩ (Callable / gradient / hessian).
// ===========================================================================

auto Hamiltonian::checkPart1() -> void {
    using namespace std;
    auto& log = Logger::instance();

    // -- helper: check int equality, log at info-level --------------------------
    auto checkInt = [&](int a, int b, string_view label) -> void {
        bool ok = (a == b);
        log.log(LogLevel::Info, "  {} {}: {} vs {}", ok ? "✓" : "✗", label, a, b);
        if (!ok) throw runtime_error(
            format("Hamiltonian consistency: {} mismatch ({} vs {})", label, a, b));
    };

    // -- helper: check lattice match -------------------------------------------
    auto checkLattice = [&](const Lattice& a, const Lattice& b,
                            string_view label) -> void {
        bool ok = latticesMatch(a, b);
        log.log(LogLevel::Info, "  {} lattice [{}]", ok ? "✓" : "✗", label);
        if (!ok) throw runtime_error(
            format("Hamiltonian consistency: lattice mismatch [{}]", label));
    };

    // ------------------------------------------------------------------
    //  lattice   (GKK, WG, VR, RHO, ATOM)
    // ------------------------------------------------------------------
    if (!canonical_lattice_) {
        if      (gkk_)  canonical_lattice_ = gkk_->meta.lattice;
        else if (wg_)   canonical_lattice_ = wg_->meta.lattice;
        else if (vr_)   canonical_lattice_ = vr_->lattice;
        else if (rho_)  canonical_lattice_ = rho_->lattice;
        else if (atom_) canonical_lattice_ = atom_->lattice;
        if (canonical_lattice_)
            log.log(LogLevel::Info, "  ✓ canonical [lattice] set");
    }
    if (canonical_lattice_) {
        if (gkk_)  checkLattice(gkk_->meta.lattice, *canonical_lattice_, "lattice [GKK]");
        if (wg_)   checkLattice(wg_->meta.lattice,  *canonical_lattice_, "lattice [WG]");
        if (vr_)   checkLattice(vr_->lattice,       *canonical_lattice_, "lattice [VR]");
        if (rho_)  checkLattice(rho_->lattice,      *canonical_lattice_, "lattice [RHO]");
        if (atom_) checkLattice(atom_->lattice,     *canonical_lattice_, "lattice [ATOM]");
    }

    // ------------------------------------------------------------------
    //  nkpt   (GKK, WG, EIGEN, OCC)
    // ------------------------------------------------------------------
    if (!canonical_nkpt_) {
        if      (gkk_)   canonical_nkpt_ = gkk_->meta.nkpt;
        else if (wg_)    canonical_nkpt_ = wg_->meta.nkpt;
        else if (eigen_) canonical_nkpt_ = eigen_->meta.nkpt;
        else if (occ_)   canonical_nkpt_ = occ_->meta.nkpt;
        if (canonical_nkpt_)
            log.log(LogLevel::Info, "  ✓ canonical [nkpt] = {}", *canonical_nkpt_);
    }
    if (canonical_nkpt_) {
        if (gkk_)   checkInt(gkk_->meta.nkpt,   *canonical_nkpt_, "nkpt [GKK]");
        if (wg_)    checkInt(wg_->meta.nkpt,    *canonical_nkpt_, "nkpt [WG]");
        if (eigen_) checkInt(eigen_->meta.nkpt, *canonical_nkpt_, "nkpt [EIGEN]");
        if (occ_)   checkInt(occ_->meta.nkpt,   *canonical_nkpt_, "nkpt [OCC]");
    }

    // ------------------------------------------------------------------
    //  nband   (WG, EIGEN, OCC)
    // ------------------------------------------------------------------
    if (!canonical_nband_) {
        if      (wg_)    canonical_nband_ = wg_->meta.nband;
        else if (eigen_) canonical_nband_ = eigen_->meta.nband;
        else if (occ_)   canonical_nband_ = occ_->meta.nband;
        if (canonical_nband_)
            log.log(LogLevel::Info, "  ✓ canonical [nband] = {}", *canonical_nband_);
    }
    if (canonical_nband_) {
        if (wg_)    checkInt(wg_->meta.nband,    *canonical_nband_, "nband [WG]");
        if (eigen_) checkInt(eigen_->meta.nband, *canonical_nband_, "nband [EIGEN]");
        if (occ_)   checkInt(occ_->meta.nband,   *canonical_nband_, "nband [OCC]");
    }

    // ------------------------------------------------------------------
    //  FFT grid n1,n2,n3   (GKK, WG, VR, RHO)
    // ------------------------------------------------------------------
    if (!canonical_fft_grid_) {
        if      (gkk_) canonical_fft_grid_ = FFTGrid{gkk_->meta.n1, gkk_->meta.n2, gkk_->meta.n3};
        else if (wg_)  canonical_fft_grid_ = FFTGrid{wg_->meta.n1,  wg_->meta.n2,  wg_->meta.n3};
        else if (vr_)  canonical_fft_grid_ = FFTGrid{vr_->meta.n1,  vr_->meta.n2,  vr_->meta.n3};
        else if (rho_) canonical_fft_grid_ = FFTGrid{rho_->meta.n1, rho_->meta.n2, rho_->meta.n3};
        if (canonical_fft_grid_)
            log.log(LogLevel::Info, "  ✓ canonical [n1,n2,n3] = {},{},{}",
                    canonical_fft_grid_->n1, canonical_fft_grid_->n2, canonical_fft_grid_->n3);
    }
    if (canonical_fft_grid_) {
        auto checkGrid = [&](int n1, int n2, int n3, string_view tag) {
            checkInt(n1, canonical_fft_grid_->n1, format("n1 [{}]", tag));
            checkInt(n2, canonical_fft_grid_->n2, format("n2 [{}]", tag));
            checkInt(n3, canonical_fft_grid_->n3, format("n3 [{}]", tag));
        };
        if (gkk_)  checkGrid(gkk_->meta.n1,  gkk_->meta.n2,  gkk_->meta.n3,  "GKK");
        if (wg_)   checkGrid(wg_->meta.n1,   wg_->meta.n2,   wg_->meta.n3,   "WG");
        if (vr_)   checkGrid(vr_->meta.n1,   vr_->meta.n2,   vr_->meta.n3,   "VR");
        if (rho_)  checkGrid(rho_->meta.n1,  rho_->meta.n2,  rho_->meta.n3,  "RHO");
    }

    // ------------------------------------------------------------------
    //  is_SO   (GKK, WG, EIGEN)
    // ------------------------------------------------------------------
    if (!canonical_is_SO_) {
        if      (gkk_)   canonical_is_SO_ = gkk_->meta.is_SO;
        else if (wg_)    canonical_is_SO_ = wg_->meta.is_SO;
        else if (eigen_) canonical_is_SO_ = eigen_->meta.is_SO;
        if (canonical_is_SO_)
            log.log(LogLevel::Info, "  ✓ canonical [is_SO] = {}", *canonical_is_SO_);
    }
    if (canonical_is_SO_) {
        if (gkk_)   checkInt(gkk_->meta.is_SO,   *canonical_is_SO_, "is_SO [GKK]");
        if (wg_)    checkInt(wg_->meta.is_SO,    *canonical_is_SO_, "is_SO [WG]");
        if (eigen_) checkInt(eigen_->meta.is_SO, *canonical_is_SO_, "is_SO [EIGEN]");
    }

    // ------------------------------------------------------------------
    //  islda   (GKK, WG, EIGEN, OCC)
    // ------------------------------------------------------------------
    if (!canonical_islda_) {
        if      (gkk_)   canonical_islda_ = gkk_->meta.islda;
        else if (wg_)    canonical_islda_ = wg_->meta.islda;
        else if (eigen_) canonical_islda_ = eigen_->meta.islda;
        else if (occ_)   canonical_islda_ = occ_->meta.islda;
        if (canonical_islda_)
            log.log(LogLevel::Info, "  ✓ canonical [islda] = {}", *canonical_islda_);
    }
    if (canonical_islda_) {
        if (gkk_)   checkInt(gkk_->meta.islda,   *canonical_islda_, "islda [GKK]");
        if (wg_)    checkInt(wg_->meta.islda,    *canonical_islda_, "islda [WG]");
        if (eigen_) checkInt(eigen_->meta.islda, *canonical_islda_, "islda [EIGEN]");
        if (occ_)   checkInt(occ_->meta.islda,   *canonical_islda_, "islda [OCC]");
    }

    // ------------------------------------------------------------------
    //  Ecut   (GKK, WG)
    // ------------------------------------------------------------------
    if (!canonical_Ecut_) {
        if      (gkk_) canonical_Ecut_ = gkk_->meta.Ecut;
        else if (wg_)  canonical_Ecut_ = wg_->meta.Ecut;
        if (canonical_Ecut_)
            log.log(LogLevel::Info, "  ✓ canonical [Ecut] = {}", *canonical_Ecut_);
    }
    if (canonical_Ecut_) {
        if (gkk_) {
            bool ok = (abs(gkk_->meta.Ecut - *canonical_Ecut_) < 1e-12);
            log.log(LogLevel::Info, "  {} Ecut [GKK]: {}", ok ? "✓" : "✗", gkk_->meta.Ecut);
            if (!ok) throw runtime_error("Hamiltonian consistency: Ecut [GKK] mismatch");
        }
        if (wg_) {
            bool ok = (abs(wg_->meta.Ecut - *canonical_Ecut_) < 1e-12);
            log.log(LogLevel::Info, "  {} Ecut [WG]: {}", ok ? "✓" : "✗", wg_->meta.Ecut);
            if (!ok) throw runtime_error("Hamiltonian consistency: Ecut [WG] mismatch");
        }
    }

    // ------------------------------------------------------------------
    //  natom   (EIGEN, ATOM)
    // ------------------------------------------------------------------
    if (!canonical_natom_) {
        if      (eigen_) canonical_natom_ = eigen_->meta.natom;
        else if (atom_)  canonical_natom_ = atom_->natom;
        if (canonical_natom_)
            log.log(LogLevel::Info, "  ✓ canonical [natom] = {}", *canonical_natom_);
    }
    if (canonical_natom_) {
        if (eigen_) checkInt(eigen_->meta.natom, *canonical_natom_, "natom [EIGEN]");
        if (atom_)  checkInt(atom_->natom,        *canonical_natom_, "natom [ATOM]");
    }

    // ------------------------------------------------------------------
    //  kpt_vec   (fractional coordinates;  GKK, EIGEN, OCC)
    //
    //  Setting the canonical requires either GKK (loadKPoint →
    //  fractional) or EIGEN / OCC (Cartesian → fractional, needs
    //  canonical_lattice_ for the conversion matrix).
    //
    //  Once set, each loaded source is compared against the canonical
    //  every time checkPart1 runs.  Comparison allows periodic integer
    //  translation (k  and  k+G  represent the same physical state).
    // ------------------------------------------------------------------
    if (!canonical_kpt_vec_ && canonical_nkpt_) {
        if (gkk_) {
            auto saved = gkk_->currentView();
            gkk_->setDataView(saved | KVecsView::Integer);
            vector<array<double,3>> kpts(*canonical_nkpt_);
            for (int ik = 0; ik < *canonical_nkpt_; ++ik) {
                auto& kv = gkk_->loadKPoint(ik);
                kpts[ik] = {kv.kPoint.x, kv.kPoint.y, kv.kPoint.z};
            }
            gkk_->setDataView(saved);
            canonical_kpt_vec_ = std::move(kpts);
            log.log(LogLevel::Info, "  ✓ canonical [kpt_vec] set from GKK ({} pts)", *canonical_nkpt_);
        } else if (canonical_lattice_ && eigen_) {
            auto A = canonical_lattice_->A();
            constexpr double TWO_PI = 2.0 * numbers::pi;
            vector<array<double,3>> kpts(eigen_->meta.nkpt);
            for (int ik = 0; ik < eigen_->meta.nkpt; ++ik) {
                auto& v = eigen_->kpt_vec[ik];
                kpts[ik] = {
                    (A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z) / TWO_PI,
                    (A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z) / TWO_PI,
                    (A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z) / TWO_PI
                };
            }
            canonical_kpt_vec_ = std::move(kpts);
            log.log(LogLevel::Info, "  ✓ canonical [kpt_vec] set from EIGEN ({} pts)", eigen_->meta.nkpt);
        } else if (canonical_lattice_ && occ_) {
            auto A = canonical_lattice_->A();
            constexpr double TWO_PI = 2.0 * numbers::pi;
            vector<array<double,3>> kpts(occ_->meta.nkpt);
            for (int ik = 0; ik < occ_->meta.nkpt; ++ik) {
                auto& v = occ_->kpt_vec[ik];
                kpts[ik] = {
                    (A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z) / TWO_PI,
                    (A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z) / TWO_PI,
                    (A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z) / TWO_PI
                };
            }
            canonical_kpt_vec_ = std::move(kpts);
            log.log(LogLevel::Info, "  ✓ canonical [kpt_vec] set from OCC ({} pts)", occ_->meta.nkpt);
        }
    }
    if (canonical_kpt_vec_) {
        constexpr double TWO_PI = 2.0 * numbers::pi;

        // Compare GKK k-points (loadKPoint → fractional) against canonical
        if (gkk_ && canonical_lattice_) {
            auto saved = gkk_->currentView();
            gkk_->setDataView(saved | KVecsView::Integer);
            bool ok = true;
            int bad_ik = -1;
            constexpr double eps = 1e-8;
            for (int ik = 0; ik < static_cast<int>(canonical_kpt_vec_->size()) && ok; ++ik) {
                auto& kv = gkk_->loadKPoint(ik);
                double dx = kv.kPoint.x - (*canonical_kpt_vec_)[ik][0];
                double dy = kv.kPoint.y - (*canonical_kpt_vec_)[ik][1];
                double dz = kv.kPoint.z - (*canonical_kpt_vec_)[ik][2];
                if (abs(dx - round(dx)) > eps || abs(dy - round(dy)) > eps || abs(dz - round(dz)) > eps) {
                    ok = false; bad_ik = ik; break;
                }
            }
            gkk_->setDataView(saved);
            if (!ok) throw runtime_error(
                format("Hamiltonian consistency: kpt_vec [GKK] mismatch at k-point {}", bad_ik));
            log.log(LogLevel::Info, "  ✓ kpt_vec [GKK]");
        }

        // Compare EIGEN k-points (Cartesian → fractional) against canonical
        if (eigen_) {
            auto A = canonical_lattice_->A();
            bool ok = true;
            int bad_ik = -1;
            constexpr double eps = 1e-8;
            for (int ik = 0; ik < eigen_->meta.nkpt && ok; ++ik) {
                auto& v = eigen_->kpt_vec[ik];
                double fx = (A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z) / TWO_PI;
                double fy = (A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z) / TWO_PI;
                double fz = (A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z) / TWO_PI;
                double dx = fx - (*canonical_kpt_vec_)[ik][0];
                double dy = fy - (*canonical_kpt_vec_)[ik][1];
                double dz = fz - (*canonical_kpt_vec_)[ik][2];
                if (abs(dx - round(dx)) > eps || abs(dy - round(dy)) > eps || abs(dz - round(dz)) > eps) {
                    ok = false; bad_ik = ik; break;
                }
            }
            if (!ok) throw runtime_error(
                format("Hamiltonian consistency: kpt_vec [EIGEN] mismatch at k-point {}", bad_ik));
            log.log(LogLevel::Info, "  ✓ kpt_vec [EIGEN]");
        }

        // Compare OCC k-points (Cartesian → fractional; looser tolerance for text format)
        if (occ_) {
            auto A = canonical_lattice_->A();
            bool ok = true;
            int bad_ik = -1;
            constexpr double eps = 1e-4;
            for (int ik = 0; ik < occ_->meta.nkpt && ok; ++ik) {
                auto& v = occ_->kpt_vec[ik];
                double fx = (A[0][0]*v.x + A[0][1]*v.y + A[0][2]*v.z) / TWO_PI;
                double fy = (A[1][0]*v.x + A[1][1]*v.y + A[1][2]*v.z) / TWO_PI;
                double fz = (A[2][0]*v.x + A[2][1]*v.y + A[2][2]*v.z) / TWO_PI;
                double dx = fx - (*canonical_kpt_vec_)[ik][0];
                double dy = fy - (*canonical_kpt_vec_)[ik][1];
                double dz = fz - (*canonical_kpt_vec_)[ik][2];
                if (abs(dx - round(dx)) > eps || abs(dy - round(dy)) > eps || abs(dz - round(dz)) > eps) {
                    ok = false; bad_ik = ik; break;
                }
            }
            if (!ok) throw runtime_error(
                format("Hamiltonian consistency: kpt_vec [OCC] mismatch at k-point {}", bad_ik));
            log.log(LogLevel::Info, "  ✓ kpt_vec [OCC]");
        }
    }

    // ------------------------------------------------------------------
    //  NCPP internal  —  no duplicate element symbols
    // ------------------------------------------------------------------
    if (elements_.size() >= 2) {
        auto& last = elements_.back().ncpp.meta;
        for (size_t i = 0; i < elements_.size() - 1; ++i) {
            if (elements_[i].ncpp.meta.element == last.element) {
                int z = ATOM::atomicNumber(last.element);
                throw std::runtime_error(
                    std::format("Hamiltonian consistency: duplicate NCPP for"
                                " element {} (Z={})\n  file 1: {}\n  file 2: {}",
                                last.element, z,
                                elements_[i].ncpp.meta.source_file,
                                last.source_file));
            }
        }
    }
}


// ===========================================================================
//  Part 2  —  file‑to‑file integrity checks
//  Pure data‑integrity checks that are not needed for H|ψ⟩ computation.
//  Each file‑pair check runs exactly once (first time both files are
//  present).  The `part2_done_` bitmask tracks which pairs have been
//  checked; loadXxx() resets the relevant bits when data is replaced.
// ===========================================================================

auto Hamiltonian::checkPart2() -> void {
    using namespace std;
    auto& log = Logger::instance();

    auto checkInt = [&](int a, int b, string_view label) -> void {
        bool ok = (a == b);
        log.log(LogLevel::Info, "  {} {}: {} vs {}", ok ? "✓" : "✗", label, a, b);
        if (!ok) throw runtime_error(
            format("Hamiltonian consistency: {} mismatch ({} vs {})", label, a, b));
    };

    // ------------------------------------------------------------------
    //  mg_nx, ng_tot_per_kpt, nnode  —  GKK vs WG
    // ------------------------------------------------------------------
    if (gkk_ && wg_ && !(part2_done_ & PART2_GKK_WG)) {
        part2_done_ |= PART2_GKK_WG;

        checkInt(gkk_->meta.mg_nx, wg_->meta.mg_nx, "mg_nx [GKK vs WG]");

        bool ok = (gkk_->meta.ng_tot_per_kpt == wg_->meta.ng_tot_per_kpt);
        log.log(LogLevel::Info, "  {} ng_tot_per_kpt [GKK vs WG]", ok ? "✓" : "✗");
        if (!ok) throw runtime_error("Hamiltonian consistency: ng_tot_per_kpt mismatch");

        checkInt(gkk_->meta.nnode, wg_->meta.nnode, "nnode [GKK vs WG]");
    }

    // ------------------------------------------------------------------
    //  nnode  —  GKK vs EIGEN
    // ------------------------------------------------------------------
    if (gkk_ && eigen_ && !(part2_done_ & PART2_GKK_EIGEN)) {
        part2_done_ |= PART2_GKK_EIGEN;
        checkInt(gkk_->meta.nnode, eigen_->meta.nnode, "nnode [GKK vs EIGEN]");
    }

    // ------------------------------------------------------------------
    //  nstate  —  VR vs RHO
    // ------------------------------------------------------------------
    if (vr_ && rho_ && !(part2_done_ & PART2_VR_RHO)) {
        part2_done_ |= PART2_VR_RHO;
        checkInt(vr_->meta.nstate, rho_->meta.nstate, "nstate [VR vs RHO]");
    }

}


// ===========================================================================
//  Part 3  —  advanced consistency checks (computationally heavy)
//
//  Not called automatically — the user invokes this explicitly.
//  Each check independently handles missing data by logging a skip message.
//
//  Internal dedup:  std::initializer_list → std::uint64_t bitmask
//  (supports up to 64 distinct check types).
// ===========================================================================

auto Hamiltonian::checkConsistencyExtended() -> void {
    checkConsistencyExtended({
        ExtendedCheck::RHOReconstruct,
        ExtendedCheck::ValenceCount,
    });
}

auto Hamiltonian::checkConsistencyExtended(std::initializer_list<ExtendedCheck> checks) -> void {
    auto& log = Logger::instance();
    log.log(LogLevel::Info, "[Hamiltonian] extended consistency check");

    auto mask = 0ull;
    for (auto c : checks)
        mask |= static_cast<std::uint64_t>(c);

    using EC = ExtendedCheck;
    // ------------------------------------------------------------------
    //  1.  RHOReconstruct  —  Σ occ·|WG|² → FFT → compare vs file RHO
    // ------------------------------------------------------------------
    if (mask & static_cast<std::uint64_t>(EC::RHOReconstruct)) {
        log.log(LogLevel::Info, "  [RHOReconstruct]");
        if (!gkk_ || !wg_ || !occ_ || !rho_) {
            log.log(LogLevel::Info, "    skipped — missing data (need GKK+WG+OCC+RHO)");
        } else {
            auto reconstructed = buildDensity(*gkk_, *wg_, *occ_);
            auto result = compareDensity(reconstructed, *rho_);

            log.log(LogLevel::Info, "    max diff = {:.4e}", result.max_diff);
            log.log(LogLevel::Info, "    RMSE     = {:.4e}  ({:.2f}% of avg|ref|)",
                    result.rmse, result.rmse / result.avg_ref * 100.0);

            if (result.rmse > 1e-2) {
                throw std::runtime_error(std::format(
                    "Hamiltonian checkConsistencyExtended: "
                    "RHO reconstruction RMSE ({:.4e}) exceeds 1e-2", result.rmse));
            }
            log.log(LogLevel::Info, "    ✓");
        }
    }

    // ------------------------------------------------------------------
    //  2.  ValenceCount  —  Σ(NCPP.z_valence × count)  ≈  ∫RHO d³r
    // ------------------------------------------------------------------
    if (mask & static_cast<std::uint64_t>(EC::ValenceCount)) {
        log.log(LogLevel::Info, "  [ValenceCount]");
        if (!atom_ || !rho_ || elements_.empty()) {
            log.log(LogLevel::Info, "    skipped — missing data (need ATOM+RHO+NCPPs)");
        } else {
            double total_valence = 0.0;
            for (auto&& t : atom_->eachType()) {
                std::string_view name = ATOM::elementName(t.z);
                for (const auto& elem : elements_) {
                    if (elem.ncpp.meta.element == name) {
                        total_valence += elem.ncpp.meta.z_valence * t.count;
                        break;
                    }
                }
            }

            double rho_int = integrateDensity(*rho_);

            log.log(LogLevel::Info, "    Σ(NCPP.z_valence × count) = {:.1f}", total_valence);
            log.log(LogLevel::Info, "    ∫RHO d³r                   = {:.6f}", rho_int);

            double diff = std::abs(rho_int - total_valence);
            if (diff > 1.0) {
                throw std::runtime_error(std::format(
                    "Hamiltonian checkConsistencyExtended: "
                    "valence electron count mismatch: "
                    "NCPP sum={:.1f}, ∫RHO={:.6f}, diff={:.4f}",
                    total_valence, rho_int, diff));
            }
            log.log(LogLevel::Info, "    ✓");
        }
    }

    // ------------------------------------------------------------------
    //  3.  NCPPAtomCoverage  —  every ATOM species has a matching NCPP
    // ------------------------------------------------------------------
    if (mask & static_cast<std::uint64_t>(EC::NCPPAtomCoverage)) {
        log.log(LogLevel::Info, "  [NCPPAtomCoverage]");
        if (!atom_ || elements_.empty()) {
            log.log(LogLevel::Info, "    skipped — missing data (need ATOM+NCPPs)");
        } else {
            for (auto&& t : atom_->eachType()) {
                std::string_view expected = ATOM::elementName(t.z);
                bool found = false;
                for (const auto& elem : elements_) {
                    if (elem.ncpp.meta.element == expected) { found = true; break; }
                }
                if (found)
                    log.log(LogLevel::Info, "    ✓ element {} covered", expected);
                else
                    throw std::runtime_error(std::format(
                        "Hamiltonian checkConsistencyExtended: "
                        "no NCPP for element {} (Z={})", expected, t.z));
            }
            log.log(LogLevel::Info, "    ✓");
        }
    }

    log.log(LogLevel::Info, "[Hamiltonian] extended consistency check complete");
}
