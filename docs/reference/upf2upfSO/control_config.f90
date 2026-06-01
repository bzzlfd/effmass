    module control_config
    !
    implicit none
    !
    character(len=200) :: job, atom_file, &
                          in_wg_file(2), &
                          out_wg_file(2), &
                          in_rho_file(2), &
                          out_rho_file(2), &
                          in_vr_file(2), &
                          out_vr_file(2), &
                          in_vext_file, &
                          out_dens_file, &
                          out_force_file, &
                          sym_file, kpt_file
    character(len=4) :: ACCURACY

    !
    integer :: num_group_k, num_nodes_b, num_group_b, &
               IQIJ_PD, &
               n1, n2, n3, &
               n1s, n2s, n3s, &
               n1L, n2L, n3L, &
               nk1, nk2, nk3, sk1, sk2, sk3, kp_inv,&
               spin, &
               igga, &
               coulomb_flag, &
               in_wg, out_wg, &
               in_rho, out_rho, &
               in_vr, out_vr, &
               in_vext, &
               flag_dens, out_dens, &
               kpt_dens(2), spin_dens(2), wg_dens(2), &
               out_force, &
               in_symm, in_kpt, &
               iQijL0_GS, isbf, &
               num_iter0, num_line0, scf_iter0, cg_bad0, &
               do_dos,&
               rmtheod1, num_mov, mv_cont, &
               num_iter1, num_line1, scf_iter1, cg_bad1, &
               local, &
               md_mth, scale_temp_VVMD, step_temp_VVMD, mx, &
               irmask_ldau, ldaU, ldau_flag, istress_cal, &
               use_hse, USE_HSE_THISETOT, USE_SBFFT_HSE, &
               HYB_GGA_TYPE, HSE_NQ1, HSE_NQ2, HSE_NQ3, &
               HSE_N1, HSE_N2, HSE_N3, iupdate_sxp_dn, NUM_MOV_HSE, &
               IRELAX_HSELDA, IS_RCUT_TYPE_SET, NUM_MOV_LDA


    integer, allocatable :: &
               cg_mth0(:), type_fermi0(:), &
               cg_mth1(:), type_fermi1(:), &
               ldau_l(:), is_LDAU_l(:), LDAU_nwfc(:), nref_type(:), &
               nref_wfc_type(:),  LDAU_occ(:)

    !
    real(kind=8) :: ecut, ecut2, ecut2L, &
                    xgga, &
                    smth, &
                    xcoul(3), &
                    totnel, &
                    ug_thred, e_thred, &
                    w_cg, w_scf, &
                    force_thred, &
                    rcut, &
                    mdstep, dtime, temper_start, temper_end, &
                    scfmth, &
                    scaler_ldau, &
                    london_s6, london_rcut, tol_stress, &
                    HSE_OMEGA, HSE_ALPHA, FAC_HSELDA, TOLRHO, &
                    MAX_RCUT

    real(kind=8), allocatable :: &
                    fermi_de0(:), &
                    mx_mtn0(:), scf_mth0(:),&
                    fermi_de1(:), &
                    mx_mth1(:), scf_mth1(:), &
                    Hubbard_U(:), &
                    LONDON_C6(:), &
                    rcut_of_type(:)
    !
    logical :: pwscf_out, out_kpt, out_symm, &
        has_london, has_vdw
    !
    end module control_config
