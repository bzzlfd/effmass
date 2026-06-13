export module pseudo.ncpp.advance;

import std;
import pseudo.ncpp;
import math.linalg;


export auto sortByL(NCPP& ncpp) -> void {
    const int n = ncpp.meta.number_of_proj;
    if (n <= 1) return;

    std::vector<int> idx(static_cast<std::size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);

    bool has_jjj = !ncpp.nonlocal.jjj.empty();

    if (has_jjj) {
        std::ranges::stable_sort(idx, [&](int a, int b) {
            int la = ncpp.nonlocal.angular_momentum[a];
            int lb = ncpp.nonlocal.angular_momentum[b];
            if (la != lb) return la < lb;
            return ncpp.nonlocal.jjj[a] < ncpp.nonlocal.jjj[b];
        });
    } else {
        std::ranges::stable_sort(idx, [&](int a, int b) {
            return ncpp.nonlocal.angular_momentum[a] < ncpp.nonlocal.angular_momentum[b];
        });
    }

    // idx[i] == i for all i → permutation is identity, the data was already
    // ordered by l and stable_sort was a no-op. Skip reordering below.
    bool is_identity = true;
    for (int i = 0; i < n; ++i) {
        if (idx[i] != i) { is_identity = false; break; }
    }
    if (is_identity) return;

    // Apply permutation to per-projector vectors (move semantics for efficiency)
    auto reorder = [&]<typename T>(std::vector<T>& vec) -> void {
        auto copy = std::move(vec);
        vec.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i)
            vec[i] = std::move(copy[static_cast<std::size_t>(idx[i])]);
    };

    reorder(ncpp.nonlocal.beta);
    reorder(ncpp.nonlocal.angular_momentum);
    reorder(ncpp.nonlocal.cutoff_index);
    reorder(ncpp.nonlocal.cutoff_radius);
    if (has_jjj) reorder(ncpp.nonlocal.jjj);

    // Permute B matrix symmetrically: B_new[i,j] = B_old[idx[i], idx[j]]
    auto B_old = ncpp.nonlocal.B.data;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ncpp.nonlocal.B[i, j] = B_old[static_cast<std::size_t>(idx[i] * n + idx[j])];
        }
    }
}


export auto diagonalizeNonlocal(NCPP& ncpp) -> void {
    const int n = ncpp.meta.number_of_proj;
    if (n == 0) return;

    // B matrix is block-diagonal in angular momentum l: projectors with
    // different l couple to different spherical harmonics and must NOT mix.
    // Diagonalize each l-block independently.

    // 1. Group projector indices by angular momentum
    int l_max = ncpp.meta.l_max;
    std::vector<std::vector<int>> l_indices(l_max + 1);
    for (int i = 0; i < n; ++i) {
        int l = ncpp.nonlocal.angular_momentum[i];
        if (l < 0 || l > l_max) {
            throw std::runtime_error(
                "diagonalizeNonlocal: projector " + std::to_string(i)
                + " has angular_momentum=" + std::to_string(l)
                + ", expected [0, " + std::to_string(l_max) + "]");
        }
        l_indices[l].push_back(i);
    }

    // 2. Save old beta and determine max cutoff
    auto beta_old = std::move(ncpp.nonlocal.beta);
    int max_cutoff = 1;
    for (const auto& b : beta_old) {
        // b can be empty (zero-size); skip to avoid zero-size allocation
        int sz = static_cast<int>(b.size());
        if (sz > max_cutoff) max_cutoff = sz;
    }

    // 3. Save original B and reset (will fill eigenvalues per block)
    auto B_orig = std::move(ncpp.nonlocal.B.data);
    ncpp.nonlocal.B.data.assign(static_cast<std::size_t>(n * n), 0.0);

    // 4. Prepare new beta storage (all zeros, max cutoff size)
    ncpp.nonlocal.beta.resize(n);
    for (auto& b : ncpp.nonlocal.beta) {
        b.assign(static_cast<std::size_t>(max_cutoff), 0.0);
    }

    // 5. Process each l-block independently
    for (int l = 0; l <= l_max; ++l) {
        const auto& idx = l_indices[l];
        const int nbl = static_cast<int>(idx.size());
        if (nbl == 0) continue;

        // 5a. Extract B_l submatrix from saved original B
        std::vector<double> B_l(static_cast<std::size_t>(nbl * nbl));
        for (int i = 0; i < nbl; ++i) {
            for (int j = 0; j < nbl; ++j) {
                B_l[static_cast<std::size_t>(i * nbl + j)] = B_orig[static_cast<std::size_t>(idx[i] * n + idx[j])];
            }
        }

        // 5b. Diagonalize this l-block
        auto result = diagonalize_jacobi(B_l, nbl);
        if (!result.converged) {
            throw std::runtime_error(
                "diagonalizeNonlocal: Jacobi diagonalization "
                "did not converge for l=" + std::to_string(l));
        }

        // 5c. Store eigenvalues into diagonal of full B
        for (int i = 0; i < nbl; ++i) {
            ncpp.nonlocal.B[idx[i], idx[i]] = result.eigenvalues[static_cast<std::size_t>(i)];
        }

        // 5d. Rotate beta projectors within this l-block
        //     beta_new[idx[k]][r] = sum_i beta_old[idx[i]][r] * U[i][k]
        //     U[i][k] = eigvec[k, i] (eigenvectors row-major)
        for (int k = 0; k < nbl; ++k) {
            auto& beta_k = ncpp.nonlocal.beta[idx[k]];
            for (int i = 0; i < nbl; ++i) {
                double u_ik = result.eigenvectors[k][i];
                const auto& beta_i = beta_old[idx[i]];
                for (int r = 0; r < static_cast<int>(beta_i.size()); ++r) {
                    beta_k[r] += beta_i[r] * u_ik;
                }
            }
        }
    }

    // 6. Update cutoff_index and cutoff_radius to max cutoff
    std::ranges::fill(ncpp.nonlocal.cutoff_index, max_cutoff);
    if (max_cutoff > 0 && !ncpp.mesh.r.empty()) {
        std::ranges::fill(ncpp.nonlocal.cutoff_radius, ncpp.mesh.r[static_cast<std::size_t>(max_cutoff - 1)]);
    }
}


export auto upf2upfSO(NCPP& ncpp) -> void {
    if (!ncpp.meta.has_so) {
        throw std::runtime_error("upf2upfSO: has_so must be true");
    }
    if (ncpp.meta.l_max < 1) {
        throw std::runtime_error("upf2upfSO: l_max must be >= 1");
    }
    if (ncpp.meta.l_max > 3) {
        throw std::runtime_error("upf2upfSO: l_max > 3 not yet supported");
    }

    const int n = ncpp.meta.mesh_size;
    const int nbeta = ncpp.meta.number_of_proj;
    const int nwave = ncpp.meta.number_of_wfc;
    const int lmax = ncpp.meta.l_max;
    const auto& r = ncpp.mesh.r;

    // Integration weights: dr[i] = (r[i+1] - r[i-1]) / 2,  i = 1..n-2
    std::vector<double> dr(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        dr[i] = (r[i + 1] - r[i - 1]) / 2.0;
    }

    // Pad beta to full mesh
    std::vector<std::vector<double>> beta_tmp(nbeta, std::vector<double>(n, 0.0));
    for (int ll = 0; ll < nbeta; ++ll) {
        const auto& src = ncpp.nonlocal.beta[ll];
        std::copy(src.begin(), src.end(), beta_tmp[ll].begin());
    }

    // Pad chi to full mesh
    std::vector<std::vector<double>> wave_tmp(nwave, std::vector<double>(n, 0.0));
    for (int ll = 0; ll < nwave; ++ll) {
        const auto& src = ncpp.pseudoWfc.chi[ll];
        std::copy(src.begin(), src.end(), wave_tmp[ll].begin());
    }

    // DD diagonal access
    auto dd_diag = [&](int ll) -> double {
        if (ll < 0 || ll >= nbeta) return 0.0;
        return ncpp.nonlocal.B.data[static_cast<std::size_t>(ll) * nbeta + ll];
    };

    // ------------------------------------------------------------
    // 1. Map projectors to (L, J, m) grouped structure
    // ------------------------------------------------------------
    std::vector<int> map_L(nbeta), map_J(nbeta);
    for (int ll = 0; ll < nbeta; ++ll) {
        map_L[ll] = ncpp.nonlocal.angular_momentum[ll];
        map_J[ll] = (ncpp.nonlocal.jjj[ll] - map_L[ll] < 0) ? 0 : 1;
    }

    // beta_arr[L][J][m][i] -- 0-based L=0..lmax, J=0,1, m=0,1
    using Beta4D = std::vector<std::vector<std::vector<std::vector<double>>>>;
    Beta4D beta_arr(lmax + 1,
        std::vector<std::vector<std::vector<double>>>(2,
            std::vector<std::vector<double>>(2,
                std::vector<double>(n, 0.0))));
    std::vector<std::vector<int>> ind(lmax + 1, std::vector<int>(2, 0));
    std::vector<std::vector<std::vector<int>>> ll_map(lmax + 1,
        std::vector<std::vector<int>>(2, std::vector<int>(2, -1)));

    for (int ll = 0; ll < nbeta; ++ll) {
        int L = map_L[ll];
        int J = map_J[ll];
        int m = ind[L][J];       // 0-based m, before increment
        ll_map[L][J][m] = ll;
        for (int i = 0; i < n; ++i) {
            beta_arr[L][J][m][i] = beta_tmp[ll][i];
        }
        ind[L][J]++;
    }

    // Duplicate: if only 1 projector for (L,J), copy m=0 -> m=1
    for (int L = 0; L <= lmax; ++L) {
        for (int J = 0; J <= 1; ++J) {
            if (ind[L][J] == 1) {
                for (int i = 0; i < n; ++i) {
                    beta_arr[L][J][1][i] = beta_arr[L][J][0][i];
                }
                ind[L][J] = 2;
            }
        }
    }

    // ------------------------------------------------------------
    // 2. Map wavefunctions to (L, J) grouped structure
    // ------------------------------------------------------------
    std::vector<std::vector<std::vector<double>>> wave_arr(lmax + 1,
        std::vector<std::vector<double>>(2, std::vector<double>(n, 0.0)));
    std::vector<std::vector<double>> occ(lmax + 1, std::vector<double>(2, 0.0));
    std::vector<int> map_Lw(nwave), map_Jw(nwave);
    int lchi_max = 0;
    for (int ll = 0; ll < nwave; ++ll) {
        int Lw = ncpp.pseudoWfc.angular_momentum[ll];
        int Jw = (ncpp.pseudoWfc.jchi[ll] - Lw < 0) ? 0 : 1;
        map_Lw[ll] = Lw;
        map_Jw[ll] = Jw;
        lchi_max = std::max(lchi_max, Lw);
        for (int i = 0; i < n; ++i) {
            wave_arr[Lw][Jw][i] = wave_tmp[ll][i];
        }
        occ[Lw][Jw] = ncpp.pseudoWfc.occupation[ll];
    }

    // ------------------------------------------------------------
    // 3. Compute beta0 and betaSO for each L = 1..lmax
    // ------------------------------------------------------------
    std::vector<std::vector<std::vector<double>>> beta0_scalar(lmax + 1,
        std::vector<std::vector<double>>(2, std::vector<double>(n, 0.0)));
    std::vector<std::vector<double>> beta_SO(lmax + 1, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> DD_0(lmax + 1, std::vector<double>(2, 0.0));
    std::vector<double> DD_SO(lmax + 1, 0.0);

    for (int L = 1; L <= lmax; ++L) {
        // Overlap sums
        double sum1 = 0.0, sum2 = 0.0, sumw = 0.0, sumw1 = 0.0, sumw2 = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            sum1  += beta_arr[L][0][0][i] * beta_arr[L][1][0][i] * dr[i];
            sum2  += beta_arr[L][0][1][i] * beta_arr[L][1][1][i] * dr[i];
            sumw  += wave_arr[L][0][i]     * wave_arr[L][1][i]     * dr[i];
            sumw1 += wave_arr[L][0][i]     * wave_arr[L][0][i]     * dr[i];
            sumw2 += wave_arr[L][1][i]     * wave_arr[L][1][i]     * dr[i];
        }

        int isign1 = (sum1 > 0) ? 1 : -1;
        int isign2 = (sum2 > 0) ? 1 : -1;
        int isignw = (sumw > 0) ? 1 : -1;

        // Mixing coefficients
        double x2 = std::sqrt(L + 1.0) / (std::sqrt(L + 1.0) + std::sqrt(L * 1.0));
        double x1 = 1.0 - x2;

        // Integrals with test functions
        double s1_1 = 0.0, s2_1 = 0.0, s1_2 = 0.0, s2_2 = 0.0;
        double s1_1T = 0.0, s2_1T = 0.0, s1_2T = 0.0, s2_2T = 0.0;
        int lw_ref = (L <= lchi_max) ? L : lchi_max;

        for (int i = 1; i < n - 1; ++i) {
            double waveT1 = x2 * beta_arr[L][1][0][i] + x1 * static_cast<double>(isign1) * beta_arr[L][0][0][i];
            double waveT2 = x2 * beta_arr[L][1][1][i] + x1 * static_cast<double>(isign2) * beta_arr[L][0][1][i];
            double waveT  = x2 * wave_arr[lw_ref][1][i] + x1 * static_cast<double>(isignw) * wave_arr[lw_ref][0][i];

            s1_1 += beta_arr[L][0][0][i] * waveT * dr[i];
            s2_1 += beta_arr[L][1][0][i] * waveT * dr[i];
            s1_2 += beta_arr[L][0][1][i] * waveT * dr[i];
            s2_2 += beta_arr[L][1][1][i] * waveT * dr[i];

            s1_1T += beta_arr[L][0][0][i] * waveT1 * dr[i];
            s2_1T += beta_arr[L][1][0][i] * waveT1 * dr[i];
            s1_2T += beta_arr[L][0][1][i] * waveT2 * dr[i];
            s2_2T += beta_arr[L][1][1][i] * waveT2 * dr[i];
        }

        // Scale by D diagonal elements
        s1_1 *= dd_diag(ll_map[L][0][0]);
        s2_1 *= dd_diag(ll_map[L][1][0]);
        s1_2 *= dd_diag(ll_map[L][0][1]);
        s2_2 *= dd_diag(ll_map[L][1][1]);
        s1_1T *= dd_diag(ll_map[L][0][0]);
        s2_1T *= dd_diag(ll_map[L][1][0]);
        s1_2T *= dd_diag(ll_map[L][0][1]);
        s2_2T *= dd_diag(ll_map[L][1][1]);

        // Construct beta0 (scalar-relativistic projectors)
        double Lfac = static_cast<double>(L) / (2.0 * L + 1.0);
        double L1fac = (L + 1.0) / (2.0 * L + 1.0);
        for (int i = 0; i < n; ++i) {
            beta0_scalar[L][0][i] = s1_1T * Lfac * beta_arr[L][0][0][i] + L1fac * s2_1T * beta_arr[L][1][0][i];
            beta0_scalar[L][1][i] = s1_2T * Lfac * beta_arr[L][0][1][i] + L1fac * s2_2T * beta_arr[L][1][1][i];
        }

        // Construct betaSO
        double so_fac = 2.0 / (2.0 * L + 1.0);
        for (int i = 0; i < n; ++i) {
            beta_SO[L][i] = so_fac * (s2_1 * beta_arr[L][1][0][i] - s1_1 * beta_arr[L][0][0][i]
                                     + s2_2 * beta_arr[L][1][1][i] - s1_2 * beta_arr[L][0][1][i]);
        }

        // Compute integrals for normalization AND DD matrix elements
        // (both use RAW unnormalized beta0/betaSO, matching Fortran order)
        double sum1t = 0.0, sum2t = 0.0, sumSOt = 0.0;
        double rsum1 = 0.0, rsum2 = 0.0, rsumSO = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            double wT1 = x2 * beta_arr[L][1][0][i] + x1 * static_cast<double>(isign1) * beta_arr[L][0][0][i];
            double wT2 = x2 * beta_arr[L][1][1][i] + x1 * static_cast<double>(isign2) * beta_arr[L][0][1][i];
            double wT  = x2 * wave_arr[lw_ref][1][i] + x1 * static_cast<double>(isignw) * wave_arr[lw_ref][0][i];

            rsum1  += beta0_scalar[L][0][i] * wT1 * dr[i];
            rsum2  += beta0_scalar[L][1][i] * wT2 * dr[i];
            rsumSO += beta_SO[L][i]          * wT  * dr[i];

            sum1t  += beta0_scalar[L][0][i] * beta0_scalar[L][0][i] * dr[i];
            sum2t  += beta0_scalar[L][1][i] * beta0_scalar[L][1][i] * dr[i];
            sumSOt += beta_SO[L][i]          * beta_SO[L][i]          * dr[i];
        }

        // Normalize
        if (sum1t > 1e-30) {
            double norm = 1.0 / std::sqrt(sum1t);
            for (int i = 0; i < n; ++i) beta0_scalar[L][0][i] *= norm;
        }
        if (sum2t > 1e-30) {
            double norm = 1.0 / std::sqrt(sum2t);
            for (int i = 0; i < n; ++i) beta0_scalar[L][1][i] *= norm;
        }
        if (sumSOt > 1e-30) {
            double norm = 1.0 / std::sqrt(sumSOt);
            for (int i = 0; i < n; ++i) beta_SO[L][i] *= norm;
        }

        // DD matrix elements (using rsum computed from RAW beta0/betaSO)
        if (std::abs(rsum1) > 1e-30)  DD_0[L][0] = sum1t / rsum1;
        if (std::abs(rsum2) > 1e-30)  DD_0[L][1] = sum2t / rsum2;
        if (std::abs(rsumSO) > 1e-30) DD_SO[L]   = sumSOt / rsumSO;
    }

    // ------------------------------------------------------------
    // 4. Assemble output projectors
    // ------------------------------------------------------------
    int nbeta_new = 2 + 3 * lmax;   // 2xL=0 + 2xscalar(L) + SO(L) for each L=1..lmax

    std::vector<std::vector<double>> new_beta(nbeta_new, std::vector<double>(n, 0.0));
    std::vector<int> new_lll(nbeta_new, 0);
    array2d<double> new_D;
    new_D.rows = nbeta_new;
    new_D.cols = nbeta_new;
    new_D.data.assign(static_cast<std::size_t>(nbeta_new * nbeta_new), 0.0);

    // L=0: pass through
    int p_idx = 0;
    for (int m = 0; m < 2 && m < ind[0][1]; ++m) {
        int ll_orig = ll_map[0][1][m];   // L=0 only has J=1 (J=L+1/2)
        if (ll_orig >= 0) {
            for (int i = 0; i < n; ++i) {
                new_beta[p_idx][i] = beta_tmp[ll_orig][i];
            }
            new_D[p_idx, p_idx] = dd_diag(ll_orig);
        }
        new_lll[p_idx] = 0;
        p_idx++;
    }

    // L=1..lmax: beta0 scalar (m=0, m=1) then betaSO
    for (int L = 1; L <= lmax; ++L) {
        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta0_scalar[L][0][i];
        }
        new_D[p_idx, p_idx] = DD_0[L][0];
        new_lll[p_idx] = L;
        p_idx++;

        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta0_scalar[L][1][i];
        }
        new_D[p_idx, p_idx] = DD_0[L][1];
        new_lll[p_idx] = L;
        p_idx++;
    }

    for (int L = 1; L <= lmax; ++L) {
        for (int i = 0; i < n; ++i) {
            new_beta[p_idx][i] = beta_SO[L][i];
        }
        new_D[p_idx, p_idx] = DD_SO[L];
        new_lll[p_idx] = L * 10;    // 10, 20, 30 for L=1,2,3
        p_idx++;
    }

    // ------------------------------------------------------------
    // 5. Assemble output wavefunctions
    // ------------------------------------------------------------
    std::vector<int> l0_wfc_idx;
    for (int ll = 0; ll < nwave; ++ll) {
        if (map_Lw[ll] == 0) l0_wfc_idx.push_back(ll);
    }
    int nwfc_new = static_cast<int>(l0_wfc_idx.size());

    std::vector<std::vector<double>> new_chi;
    std::vector<int> new_lchi;
    std::vector<double> new_oc;

    for (int jj : l0_wfc_idx) {
        new_chi.push_back(wave_tmp[jj]);
        new_lchi.push_back(0);
        new_oc.push_back(ncpp.pseudoWfc.occupation[jj]);
    }

    for (int L = 1; L <= lmax; ++L) {
        double sumw = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            sumw += wave_arr[L][0][i] * wave_arr[L][1][i] * dr[i];
        }
        int signw = (sumw > 0) ? 1 : -1;

        double x2 = std::sqrt(L + 1.0) / (std::sqrt(L + 1.0) + std::sqrt(L * 1.0));
        double x1 = 1.0 - x2;

        double sumt = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            double wT = x2 * wave_arr[L][1][i] + x1 * static_cast<double>(signw) * wave_arr[L][0][i];
            sumt += wT * wT * dr[i];
        }

        if (sumt > 1e-30) {
            std::vector<double> chi_L(n);
            double norm = 1.0 / std::sqrt(sumt);
            for (int i = 0; i < n; ++i) {
                chi_L[i] = (x2 * wave_arr[L][1][i] + x1 * static_cast<double>(signw) * wave_arr[L][0][i]) * norm;
            }
            new_chi.push_back(std::move(chi_L));
            new_lchi.push_back(L);
            new_oc.push_back(occ[L][0] + occ[L][1]);
            nwfc_new++;
        }
    }

    // ------------------------------------------------------------
    // 6. Update NCPP data
    // ------------------------------------------------------------
    // Truncate trailing zeros from new projectors
    for (auto& b : new_beta) {
        int cut = n;
        while (cut > 0 && std::abs(b[static_cast<std::size_t>(cut) - 1]) < 1e-30) --cut;
        b.resize(cut);
    }

    // Truncate trailing zeros from new wavefunctions
    for (auto& c : new_chi) {
        int cut = n;
        while (cut > 0 && std::abs(c[static_cast<std::size_t>(cut) - 1]) < 1e-30) --cut;
        c.resize(cut);
    }

    ncpp.meta.has_so = false;
    ncpp.meta.number_of_proj = nbeta_new;
    ncpp.meta.number_of_wfc = nwfc_new;

    ncpp.nonlocal.beta = std::move(new_beta);
    ncpp.nonlocal.angular_momentum = std::move(new_lll);
    ncpp.nonlocal.jjj.clear();
    ncpp.nonlocal.B = std::move(new_D);
    ncpp.nonlocal.cutoff_index.assign(nbeta_new, n);
    ncpp.nonlocal.cutoff_radius.assign(nbeta_new,
        n > 0 ? ncpp.mesh.r.back() : 0.0);

    ncpp.pseudoWfc.chi = std::move(new_chi);
    ncpp.pseudoWfc.angular_momentum = std::move(new_lchi);
    ncpp.pseudoWfc.jchi.clear();
    ncpp.pseudoWfc.occupation = std::move(new_oc);
}
