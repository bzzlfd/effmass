!=====================================================!

subroutine read_upf_v2 (upfname, upf, iflag)
    use default_element
    use iotk_module
    use upf_module
    use double_precision
    

    implicit none
    integer :: iflag  ! =0, no change for upf%rho_atc,upf%vloc, =1: change them
    real(dp) :: pi
    type(atompp) :: upf
    CHARACTER(len=iotk_namlenx) :: root
    CHARACTER(len=iotk_attlenx) :: attr
    integer :: ierr_
    integer :: ioupf = 10000
    character(len=*) :: upfname
    logical :: found
    real, allocatable :: zeros(:)
    integer :: nb, mb, lm, ln, nw, l, nmb, i, tmpint
    !
    call initial_elements ()
    !
    !write ( 6, * ) trim(upfname)
    call iotk_open_read ( ioupf, trim(upfname), root=root, ierr=ierr_ )
    if ( ierr_ /= 0 .or. index(root,"UPF") == 0 ) then
	write ( 6, * ) "upfv2: something wrong with the file "//trim(upfname)
	call iotk_close_read ( ioupf )
	stop
    end if

    call iotk_scan_begin ( ioupf, "PP_INFO", found=found )
    if ( found ) then
	call iotk_scan_begin ( ioupf, "PP_INPUTFILE", found=found )
	if ( found ) then
	    call iotk_scan_end ( ioupf, "PP_INPUTFILE" )
	end if
	call iotk_scan_end ( ioupf, "PP_INFO" )
    end if
    !write (6, * ) "PP_INFO"

    call iotk_scan_empty ( ioupf, "PP_HEADER", attr=attr, found=found )
    if ( found ) then
	call iotk_scan_attr ( attr, "generated", upf%generated, default=" " )
	call iotk_scan_attr ( attr, "author", upf%author, default="anonymous" )
	call iotk_scan_attr ( attr, "date", upf%date, default=" " )
	call iotk_scan_attr ( attr, "comment", upf%comment, default=" " )
	call iotk_scan_attr ( attr, "element", upf%psd, default=" " )
	upf%psd = adjustl(upf%psd)
	if ( upf%psd == " " ) then
	    write ( 6, * ) "please add the 'element' in the tag 'PP_HEADER' "
	    stop
	else
	    do i = 1, 108
		if ( upf%psd == elements(i)%label ) then
		    upf%num = elements(i)%atnum
		    exit
		end if
	    end do
	end if
	call iotk_scan_attr ( attr, "pseudo_type", upf%typ )
	call iotk_scan_attr ( attr, "relativistic", upf%rel )
	call iotk_scan_attr ( attr, "is_ultrasoft", upf%tvanp )
	call iotk_scan_attr ( attr, "has_so", upf%has_so, default=.false. )
	call iotk_scan_attr ( attr, "has_wfc", upf%has_wfc, default=.false. )
	call iotk_scan_attr ( attr, "core_correction", upf%nlcc )
	call iotk_scan_attr ( attr, "functional", upf%dft )
	call iotk_scan_attr ( attr, "z_valence", upf%zp )
	call iotk_scan_attr ( attr, "total_psenergy", upf%etotps, default=0._dp )
	if ( .not. upf%tvanp ) then
	    call iotk_scan_attr ( attr, "wfc_cutoff", upf%ecutwfc, default=elements(upf%num)%ecutwfc(1) )
	    if ( abs(upf%ecutwfc) < 0.00001 ) upf%ecutwfc = elements(upf%num)%ecutwfc(1)
	    call iotk_scan_attr ( attr, "rho_cutoff", upf%ecutrho, default=2*elements(upf%num)%ecutwfc(1) )
	    if ( abs(upf%ecutrho) < 0.00001 ) upf%ecutrho = 2*elements(upf%num)%ecutwfc(1)
	else
	    call iotk_scan_attr ( attr, "wfc_cutoff", upf%ecutwfc, default=elements(upf%num)%ecutwfc(2) )
	    if ( abs(upf%ecutwfc) < 0.00001 ) upf%ecutwfc = elements(upf%num)%ecutwfc(2)
	    call iotk_scan_attr ( attr, "rho_cutoff", upf%ecutrho, default=4*elements(upf%num)%ecutwfc(2) )
	    if ( abs(upf%ecutrho) < 0.00001 ) upf%ecutrho = 2*elements(upf%num)%ecutwfc(2)
	end if
	call iotk_scan_attr ( attr, "l_max", upf%lmax, default=0 )
	call iotk_scan_attr ( attr, "l_max_rho", upf%lmax_rho, default=2*upf%lmax )
	call iotk_scan_attr ( attr, "l_local", upf%lloc, default=0 )
	call iotk_scan_attr ( attr, "mesh_size", upf%mesh )
	call iotk_scan_attr ( attr, "number_of_wfc", upf%nwfc )
	call iotk_scan_attr ( attr, "number_of_proj", upf%nbeta )
    else
	write ( 6, * ) "Not found the tag: ""PP_HEADER"" "
	stop
    end if
    upf%mass = elements(upf%num)%mass
    allocate ( zeros(upf%mesh) )
    zeros = 0.0_dp
    !write (6, * ) "PP_HEADER"
	
    call iotk_scan_begin ( ioupf, "PP_MESH", attr=attr, found=found )
    if ( found ) then
	call iotk_scan_attr ( attr, "dx", upf%dx, default=0._dp )
	call iotk_scan_attr ( attr, "mesh", upf%mesh, default=upf%mesh )
	call iotk_scan_attr ( attr, "xmin", upf%xmin, default=0._dp )
	call iotk_scan_attr ( attr, "rmax", upf%rmax, default=0._dp )
	call iotk_scan_attr ( attr, "zmesh", upf%zmesh, default=0._dp )
	allocate ( upf%r(upf%mesh), upf%rab(upf%mesh) )
	call iotk_scan_dat ( ioupf, "PP_R", upf%r(1:upf%mesh) )
	call iotk_scan_dat ( ioupf, "PP_RAB", upf%rab(1:upf%mesh) )
    call iotk_scan_end ( ioupf, "PP_MESH" )
    else
	write ( 6, * ) "Not found the tag: ""PP_MESH"" "
	stop
    end if
    !write (6, * ) "PP_MESH"
    
    allocate ( upf%rho_atc(upf%mesh) )
    if ( upf%nlcc ) then
	call iotk_scan_dat ( ioupf, "PP_NLCC", upf%rho_atc )
    else
	upf%rho_atc(1:upf%mesh) = 0._dp
    end if
    !write (6, * ) "PP_NLCC"

    allocate ( upf%vloc(upf%mesh) )
    call iotk_scan_dat ( ioupf, "PP_LOCAL", upf%vloc )
    !write (6, * ) "PP_LOCAL"
    if (upf%nbeta == 0) then
        upf%nqf = 0
        upf%nqlc= 0
        upf%qqq_eps= -1._dp
        upf%kkbeta = 0
        ALLOCATE( upf%kbeta(1),upf%lll(1),upf%beta(upf%mesh,1),&
                  upf%dion(1,1),upf%rinner(1),upf%qqq(1,1),&
                  upf%qfunc(upf%mesh,1),upf%qfcoef(1,1,1,1),&
                  upf%rcut(1),upf%rcutus(1),upf%els_beta(1))
    else
 
    call iotk_scan_begin ( ioupf, "PP_NONLOCAL", found=found )
    if ( found ) then
    allocate ( upf%kbeta(upf%nbeta), upf%lll(upf%nbeta), upf%beta(upf%mesh,upf%nbeta), &
	upf%dion(upf%nbeta,upf%nbeta), upf%rcut(upf%nbeta), upf%rcutus(upf%nbeta), &
	upf%els_beta(upf%nbeta) )
    do nb = 1, upf%nbeta
	call iotk_scan_dat( ioupf, "PP_BETA"//iotk_index( nb ), &
	    upf%beta(:,nb), attr=attr )
	    call iotk_scan_attr ( attr, "label", upf%els_beta(nb), default="Xn" )
	    call iotk_scan_attr ( attr, "angular_momentum", upf%lll(nb) )
	    call iotk_scan_attr ( attr, "cutoff_radius_index", upf%kbeta(nb), default=upf%mesh )
	    call iotk_scan_attr ( attr, "cutoff_radius", upf%rcut(nb), default=0._dp )
	    call iotk_scan_attr ( attr, "ultrasoft_cutoff_radius", upf%rcutus(nb), default=0._dp )
	    if ( upf%rcutus(nb)==0._dp ) call iotk_scan_attr ( attr, &
	    "norm_conserving_radius", upf%rcutus(nb), default=0._dp )
    end do
    !write (6, * ) "PP_BETA"

    call iotk_scan_dat (ioupf, "PP_DIJ", upf%dion, attr=attr )
    !write (6, * ) "PP_DIJ"
    if ( upf%tvanp ) then
	call iotk_scan_begin ( ioupf, "PP_AUGMENTATION", attr=attr )
	    call iotk_scan_attr ( attr, "q_with_l", upf%q_with_l )
	    call iotk_scan_attr ( attr, "nqf", upf%nqf )
	    call iotk_scan_attr ( attr, "nqlc", upf%nqlc, default=2*upf%lmax+1 )
	    call iotk_scan_attr ( attr, "augmentation_epsilon", upf%qqq_eps, default =&
	-1._dp )
	allocate ( upf%rinner(upf%nqlc) )
	allocate ( upf%qqq(upf%nbeta,upf%nbeta) )
	if ( upf%q_with_l ) then
	    allocate ( upf%qfuncl(upf%mesh,upf%nbeta*(upf%nbeta+1)/2,0:2*upf%lmax) )
	    upf%qfuncl = 0._dp
	else
	    allocate ( upf%qfunc(upf%mesh,upf%nbeta*(upf%nbeta+1)/2) )
	end if
	call iotk_scan_dat ( ioupf, "PP_Q", upf%qqq )
	if ( upf%nqf <= 0 ) then
	    upf%rinner = 0._dp
	    allocate ( upf%qfcoef(1,1,1,1) )
	    upf%qfcoef = 0._dp
	else
	    allocate ( upf%qfcoef(max(upf%nqf,1), upf%nqlc, upf%nbeta, upf%nbeta) )
	    call iotk_scan_dat ( ioupf, "PP_QFCOEF", upf%qfcoef, attr=attr )
	    call iotk_scan_dat ( ioupf, "PP_RINNER", upf%rinner, attr=attr )
	end if
	!write (6, * ) "PP_QFCOEF", "PP_RINNER"
	if ( upf%tvanp ) then
	    do nb = 1, upf%nbeta
		ln = upf%lll(nb)
		do mb = nb, upf%nbeta
		    lm = upf%lll(mb)
		    nmb = mb * (mb - 1) / 2 + nb
		    if ( upf%q_with_l ) then
			do l = abs(ln-lm), ln+lm, 2
                  CALL iotk_scan_dat(ioupf, 'PP_QIJL'//iotk_index((/nb,mb,l/)),&
                                    upf%qfuncl(:,nmb,l), found=found )
			if (.not. found) upf%qfuncl(:,nmb,l) = zeros
			end do
		    else
               CALL iotk_scan_dat(ioupf, 'PP_QIJ'//iotk_index((/nb,mb/)),&
                                 upf%qfunc(:,nmb), found=found)
			     if (.not.found) upf%qfunc(:,nmb) = zeros
		    end if
		end do
	    end do
	end if
	call iotk_scan_end ( ioupf, "PP_AUGMENTATION" )
    end if
    upf%kkbeta = maxval( upf%kbeta(1:upf%nbeta) )
    call iotk_scan_end ( ioupf, "PP_NONLOCAL" )
    else
	write ( 6, * ) "Not found the tag: ""PP_NONLOCAL"" "
	stop
    end if
    end if
    !write (6, * ) "PP_NONLOCAL"

    call iotk_scan_begin ( ioupf, "PP_PSWFC", found=found )
    if ( found ) then
    allocate ( upf%chi(upf%mesh,upf%nwfc) )
    allocate ( upf%els(upf%nwfc), upf%oc(upf%nwfc), upf%lchi(upf%nwfc), upf%nchi(upf%nwfc), &
	upf%rcut_chi(upf%nwfc), upf%rcutus_chi(upf%nwfc), upf%epseu(upf%nwfc) )
    do nw = 1, upf%nwfc
	call iotk_scan_dat ( ioupf, "PP_CHI"//iotk_index(nw), &
	    upf%chi(:,nw), attr=attr )
	call iotk_scan_attr ( attr, "label", upf%els(nw), default="Xn" )
	call iotk_scan_attr ( attr, "l", upf%lchi(nw) )
	call iotk_scan_attr ( attr, "occupation", upf%oc(nw) )
	call iotk_scan_attr ( attr, "n", upf%nchi(nw), default=upf%lchi(nw)-1)
	call iotk_scan_attr ( attr, "pseudo_energy", upf%epseu(nw), default=0._dp )
	call iotk_scan_attr ( attr, "cutoff_radius", upf%rcut_chi(nw),&
	    default=0._dp)
	call iotk_scan_attr ( attr, "ultrasoft_cutoff_radius", upf%rcutus_chi(nw),&
	    default=0._dp )
    end do
    call iotk_scan_end ( ioupf, "PP_PSWFC" )
    else
	write ( 6, * ) "Not found the tag: ""PP_PSWFC"" "
    end if
    !write (6, * ) "PP_PSWFC"

    if ( upf%has_wfc ) then
    call iotk_scan_begin ( ioupf, "PP_FULL_WFC" )
    allocate ( upf%aewfc(upf%mesh,upf%nbeta) )
    do nb = 1, upf%nbeta
	call iotk_scan_dat ( ioupf, "PP_AEWFC"//iotk_index(nb), &
	    upf%aewfc(:,nb), attr=attr )
    end do
    allocate ( upf%pswfc(upf%mesh,upf%nbeta) )
    do nb = 1, upf%nbeta
	call iotk_scan_dat ( ioupf, "PP_PSWFC"//iotk_index(nb), &
	    upf%pswfc(:,nb), attr=attr )
    end do
    call iotk_scan_end ( ioupf, "PP_FULL_WFC" )
    end if
    !write (6, * ) "PP_FULL_WFC"
    allocate ( upf%rho_at(upf%mesh) )
    call iotk_scan_dat ( ioupf, "PP_RHOATOM", upf%rho_at )
    !write (6, * ) "PP_RHOATOM"
    if ( upf%has_so ) then
    call iotk_scan_begin ( ioupf, "PP_SPIN_ORB" )
    allocate ( upf%nn(upf%nwfc), upf%jchi(upf%nwfc) )
    do nw = 1, upf%nwfc
	call iotk_scan_empty ( ioupf, "PP_RELWFC"//iotk_index(nw), &
	    attr=attr )
	    call iotk_scan_attr ( attr, "nn", upf%nn(nw) )
	    call iotk_scan_attr ( attr, "jchi", upf%jchi(nw) )
    end do
    allocate ( upf%jjj(upf%nbeta) )
    do nb = 1, upf%nbeta
	call iotk_scan_empty ( ioupf, "PP_RELBETA"//iotk_index(nb), &
	    attr=attr )
	    call iotk_scan_attr ( attr, "lll", upf%lll(nb) )
	    call iotk_scan_attr ( attr, "jjj", upf%jjj(nb) )
    end do
    call iotk_scan_end ( ioupf, "PP_SPIN_ORB" )
    end if
    !write (6, * ) "PP_SPIN_ORB"
    
    call iotk_close_read ( ioupf )

    !call debug_read ()
    ! convert upf to uspp unit
    !

  !  cccccccccc, this has been changed!

 
    if(iflag.eq.1) then      ! this is a very bad idea, to change it here, LWW 
    pi = atan(1.) * 4.0_dp
    do i = 2, upf%mesh
	upf%vloc(i) = upf%vloc(i)*upf%r(i)
    end do
    upf%vloc(1) = 0.0
    if ( upf%nlcc ) then
	upf%rho_atc(1) = 0.0
	do i = 2, upf%mesh
	    upf%rho_atc(i) = upf%rho_atc(i) * 4.0 * pi * upf%r(i) **2
	end do
    end if
    endif

end subroutine read_upf_v2

!=====================================================!
!=====================================================!
!subroutine debug_read ()
!
!    use iotk_module 
!    use double_precision
!    use upf_module
!
!    implicit none
!    integer :: nb, mb, lm, ln, nw, l, nmb, i
!    CHARACTER(len=iotk_namlenx) :: root
!    CHARACTER(len=iotk_attlenx) :: attr
!
!    call iotk_open_write ( 110, "have_read", root="UPF" )
!    call iotk_write_begin ( 110, "PP_INFO" )
!    call iotk_write_end ( 110, "PP_INFO" )
!    
!	call iotk_write_attr ( attr, "generate", trim(generated), first=.true. )
!	call iotk_write_attr ( attr, "author", trim(author), newline=.true. )
!	call iotk_write_attr ( attr, "date", trim(date), newline=.true. )
!	call iotk_write_attr ( attr, "comment", trim(comment), newline=.true. )
!	call iotk_write_attr ( attr, "element", trim(psd), newline=.true. )
!        call iotk_write_attr ( attr, "pseudo_type", trim(typ), newline=.true. )
!	call iotk_write_attr ( attr, "relativistic", rel, newline=.true. )
!        call iotk_write_attr ( attr, "is_ultrasoft", tvanp, newline=.true. )
!	call iotk_write_attr ( attr, "has_so", has_so, newline=.true. )
!        call iotk_write_attr ( attr, "has_wfc", has_wfc, newline=.true. )
!	call iotk_write_attr ( attr, "core_correction", nlcc, newline=.true. )
!	call iotk_write_attr ( attr, "functional", trim(dft_buffer), newline=.true. )
!	call iotk_write_attr ( attr, "z_valence", zp, newline=.true. )
!	call iotk_write_attr ( attr, "total_psenergy", etotps, newline=.true. )
!	call iotk_write_attr ( attr, "wfc_cutoff", ecutwfc, newline=.true. )
!	call iotk_write_attr ( attr, "rho_cutoff", ecutrho, newline=.true. )
!	call iotk_write_attr ( attr, "l_max", lmax, newline=.true. )
!	call iotk_write_attr ( attr, "l_max_rho", lmax_rho, newline=.true. )
!	call iotk_write_attr ( attr, "l_local", lloc, newline=.true. )
!	call iotk_write_attr ( attr, "mesh_size", mesh, newline=.true. )
!	call iotk_write_attr ( attr, "number_of_wfc", nwfc, newline=.true. )
!	call iotk_write_attr ( attr, "number_of_proj", nbeta, newline=.true. )
!	call iotk_write_attr ( attr, "mass", mass, newline=.true. )
!    call iotk_write_empty ( 110, "PP_HEADER", attr=attr )
!
!	call iotk_write_attr ( attr, "dx", dx, first=.true. )
!	call iotk_write_attr ( attr, "mesh", mesh )
!	call iotk_write_attr ( attr, "xmin", xmin )
!	call iotk_write_attr ( attr, "rmax", rmax )
!	call iotk_write_attr ( attr, "zmesh", zmesh )
!    call iotk_write_begin ( 110, "PP_MESH", attr=attr )
!        call iotk_write_dat ( 110, "PP_R", r(1:mesh), columns=4 )
!	call iotk_write_dat ( 110, "PP_RAB", rab(1:mesh), columns=4 )
!    call iotk_write_end ( 110, "PP_MESH" )
!
!    if ( nlcc ) then
!	call iotk_write_dat ( 110, "PP_NLCC", rho_atc, columns=4 )
!    end if
!
!    call iotk_write_dat ( 110, "PP_LOCAL", vloc, columns=4 )
!    
!    call iotk_write_begin ( 110, "PP_NONLOCAL" )
!    do nb = 1, nbeta
!	    call iotk_write_attr ( attr, "index", nb, first=.true. )
!	    call iotk_write_attr ( attr, "label", els_beta(nb) )
!	    call iotk_write_attr ( attr, "angular_momentum", lll(nb) )
!	    call iotk_write_attr ( attr, "cutoff_radius_index", kbeta(nb) )
!	    call iotk_write_attr ( attr, "cutoff_radius", rcut(nb) )
!	    call iotk_write_attr ( attr, "ultrasoft_cutoff_radius", rcutus(nb) )
!	    if ( rcutus(nb) == 0._dp ) call iotk_write_attr ( attr, "norm_conserving_radius", rcutus(nb) )
!	call iotk_write_dat ( 110, "PP_BETA"//iotk_index(nb), beta(:,nb), attr=attr, columns=4 )
!    end do
!
!    call iotk_write_dat ( 110, "PP_DIJ", dion, columns=4 )
!    if ( tvanp ) then
!	call iotk_write_attr ( attr, "q_with_l", q_with_l, first=.true. )
!	call iotk_write_attr ( attr, "nqf", nqf )
!	call iotk_write_attr ( attr, "nqlc", nqlc )
!	call iotk_write_begin ( 110, "PP_AUGMENTATION", attr=attr )
!	call iotk_write_dat ( 110, "PP_Q", qqq, columns=4 )
!	if ( nqf > 0 ) then
!	    call iotk_write_comment ( 110, " polinomial expansion of Q_ij at small radius " )
!	    call iotk_write_dat ( 110, "PP_QFCOEF", qfcoef, attr=attr, columns=4 )
!	    call iotk_write_dat ( 110, "PP_RINNER", rinner, attr=attr, columns=4 )
!	end if
!	do nb = 1, nbeta
!	    ln = lll(nb)
!	    do mb = nb, nbeta
!		lm = lll(mb)
!		nmb = mb * (mb-1) / 2 + nb
!		if ( q_with_l ) then
!		    do l = abs(ln-lm), ln+lm, 2
!			call iotk_write_attr ( attr, "first_index", nb, first=.true. )
!			call iotk_write_attr ( attr, "second_index", mb)
!			call iotk_write_attr ( attr, "composite_index", nmb )
!			call iotk_write_attr ( attr, "angular_momentum", l )
!			call iotk_write_dat ( 110, "PP_QIJL"//iotk_index((/nb,mb,l/)), qfuncl(:,nmb,l), attr=attr, columns=4 )
!		    end do
!		else
!		    call iotk_write_attr ( attr, "first_index", nb, first=.true. )
!		    call iotk_write_attr ( attr, "second_index", mb)
!		    call iotk_write_attr ( attr, "composite_index", nmb )
!		    call iotk_write_dat ( 110, "PP_QIJ"//iotk_index((/nb,mb/)), qfunc(:,nmb), attr=attr, columns=4 )
!		end if
!	    end do
!	end do
!	call iotk_write_end ( 110, "PP_AUGMENTATION" )
!    end if
!    call iotk_write_end ( 110, "PP_NONLOCAL" )
!    
!    call iotk_write_begin ( 110, "PP_PSWFC" )
!    do nw = 1, nwfc
!	call iotk_write_attr ( attr, "index", nw, first=.true. )
!	call iotk_write_attr ( attr, "label", els(nw) )
!	call iotk_write_attr ( attr, "l", lchi(nw) )
!	call iotk_write_attr ( attr, "occupation", oc(nw) )
!	call iotk_write_attr ( attr, "n", nchi(nw) )
!	call iotk_write_attr ( attr, "pseudo_energy", epseu(nw) )
!	call iotk_write_attr ( attr, "cutoff_radius", rcut_chi(nw) )
!	call iotk_write_attr ( attr, "ultrasoft_cutoff_radius", rcutus_chi(nw) )
!	call iotk_write_dat ( 110, "PP_CHI"//iotk_index(nw), chi(:,nw), columns=4, attr=attr )
!    end do
!    call iotk_write_end ( 110, "PP_PSWFC" )
!    
!    if ( has_wfc ) then
!	CALL iotk_write_attr(attr, 'number_of_wfc', nbeta, first=.true.)
!	CALL iotk_write_begin(110, 'PP_FULL_WFC', attr=attr)
!	! All-electron wavefunctions corresponding to beta functions
!	DO nb = 1, nbeta
!	    CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
!	    CALL iotk_write_attr(attr, 'label',      els_beta(nb))
!	    CALL iotk_write_attr(attr, 'l',          lll(nb))
!	    CALL iotk_write_dat(110, 'PP_AEWFC'//iotk_index(nb), &
!           aewfc(:,nb), columns=4, attr=attr)
!	ENDDO
!	DO nb = 1,nbeta
!	    CALL iotk_write_attr(attr, 'index',      nb, first=.true.)
!	    CALL iotk_write_attr(attr, 'label',      els_beta(nb))
!	    CALL iotk_write_attr(attr, 'l',          lll(nb))
!	    CALL iotk_write_dat(110, 'PP_PSWFC'//iotk_index(nb), &
!           pswfc(:,nb), columns=4, attr=attr)
!	ENDDO
!	! Finalize
!	CALL iotk_write_end(110, 'PP_FULL_WFC')
!    end if
!
!    CALL iotk_write_dat(110, 'PP_RHOATOM', rho_at, columns=4)
!
!    if ( has_so ) then
!	CALL iotk_write_begin(110, 'PP_SPIN_ORB')
!    !
!	DO nw = 1,nwfc
!	    CALL iotk_write_attr(attr, 'index', nw, first=.true.)
!	    CALL iotk_write_attr(attr, 'els',   els(nw))
!	    CALL iotk_write_attr(attr, 'nn',    nn(nw))
!	    CALL iotk_write_attr(attr, 'lchi',  lchi(nw))
!	    CALL iotk_write_attr(attr, 'jchi',  jchi(nw))
!	    CALL iotk_write_attr(attr, 'oc',    oc(nw))
!	    CALL iotk_write_empty(110, 'PP_RELWFC'//iotk_index(nw),&
!           attr=attr)
!	ENDDO
!   !
!	DO nb = 1,nbeta
!	    CALL iotk_write_attr(attr, 'index', nb, first=.true.)
!	    CALL iotk_write_attr(attr, 'lll',   lll(nb))
!	    CALL iotk_write_attr(attr, 'jjj',   jjj(nb))
!	    CALL iotk_write_empty(110, 'PP_RELBETA'//iotk_index(nb),&
!           attr=attr)
!	ENDDO
!   !
!	CALL iotk_write_end(110, 'PP_SPIN_ORB')
!    end if
!
!    call iotk_close_write ( 110 )  
!
!end subroutine debug_read 
