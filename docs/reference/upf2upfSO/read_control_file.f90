    subroutine read_control_file (filename)
    !
    use pseudo_potential_upf !upf
    use file_unit, only: control_file_unit, long_control_file_unit, atom_file_unit, io_stdout
    use control_config
    use atom_base
    use pse_pot_config
    use xc_f90_lib_m
    use xc
    !
    implicit none
    !
    logical :: alive = .false.
    logical :: readit = .false., setit = .false., right_logical
    character(len=*) :: filename
    !
    integer :: stat, i, j, ierr, ia, i_plus, iflag_hse
    !
    character(len=200) :: right,temp_right, tmp_char, message, f_xatom, &
    temp_char_xc, temp_char
    real*8, parameter :: Hartree_ev=27.211396d0
    !
    real(kind=8) :: fackpt, dd11, dd22, dd33, ALX, ALY, ALZ
    REAL(KIND=8), PARAMETER :: A_AU_1 = 0.529177249
    real(kind=8) :: standard_al, factor
    integer :: irmethod_1, Nimage, itype_string_1, itype_xatom, imv_cont,&
    inumber_block_psi,iwrite, write2memory_7
    real(kind=8) :: tolforce, ak_string_1, Etotneb1_1, Etotneb2_1, dtstart
    integer :: tmp_int1, tmp_int2, tmp_int3
    character(len=30) mixing, vdw_method
    character(len=10) converge
    !
    !*******************************
    !
    inquire (file = trim(adjustl(filename)), exist = alive)
    if (alive) then
        open (unit = control_file_unit, file = trim(adjustl(filename)), action = "read")
        open (unit = long_control_file_unit, file = "etot.input.long", action = "write")
    else
        call error_message (trim(filename)//" doesn't exist!")
    endif
    !
    !*******************************
    !
    num_group_b = 1
    read (control_file_unit, *, iostat = stat) num_nodes_b, num_group_k
    if (stat /= 0) then
        call error_message ("the first line in "//trim(filename)//" should be 'num_nodes_b, num_group_k' ")
    else
        write (long_control_file_unit, *) num_nodes_b, num_group_k
    endif
    !
    !*******************************
    !
    call read_key_words (control_file_unit, "IN.ATOM", LEN("IN.ATOM"), right, readit)
    if (readit) then
        read (right, "(A200)") atom_file
    else
        atom_file = "atom.config"
        write (io_stdout, *)"Not set ""IN.ATOM"" in ""etot.input"", &
                             read the default atomic position file: ""atom.config"""
    endif
    write (long_control_file_unit, *) "    IN.ATOM    = "//trim(atom_file)
    !
    inquire (file = trim(atom_file), exist = alive)
    if (alive) then
        open (unit = atom_file_unit, file = trim(atom_file), action = "read")
    else
        call error_message (""""//trim(atom_file)//""""//" doesn't exist!")
    endif
    !
    !*******************************
    !
    call read_atom_file ()
    call cal_ityp ()
    !
    !*******************************
    !
    call read_key_words (control_file_unit, "JOB", len("JOB"), right, readit)
    if (readit) then
        read (right, "(A200)") tmp_char
	    call transform_to_upper (tmp_char, job)
	    job = adjustl(job)
	    select case (trim(job))
        case ("SCF", "NONSCF", "DOS", "NEB")
	    !call check_atom_mov ()
        case ("RELAX", "MD")
	    call check_atom_mov ()
        case default
	    call error_message &
	    ("In ""etot.input"", JOB should be SCF/NONSCF/MD/RELAX/NEB/DOS")
	    end select
        write (long_control_file_unit, *) &
	"    JOB    = ", trim(job)
    else
	    call error_message &
	("In ""etot.input"", JOB should be SCF/NONSCF/MD/RELAX/NEB/DOS")
    endif    
    !
    CALL read_key_words ( control_file_unit, 'RELAX_DETAIL', LEN('RELAX_DETAIL'),&
                              right, readit )
    if(readit) then
        READ ( right, *, iostat = stat) rmtheod1, num_mov, force_thred, istress_cal, tol_stress
        dtstart = 1.d0
        mv_cont = 0
        if (stat /= 0) then
            read (right, *) rmtheod1, num_mov, force_thred
            istress_cal = 0
            tol_stress = 0.0d0
        end if
        tolforce = tolforce
        num_mov = num_mov*1.001
        mv_cont = mv_cont*1.001
    else
        rmtheod1 = 1
        num_mov = 100
        force_thred = 0.01
        istress_cal = 0
        tol_stress = 0.0d0
    endif
    write(tmp_char,2005) rmtheod1, num_mov, force_thred, istress_cal, tol_stress
2005   format(i1,1x,i5,1x,E15.6,1x,i3,1x,E15.6)
    WRITE (long_control_file_unit, *) '    RELAX_DETAIL    = ', TRIM(ADJUSTL(tmp_char))
    !
    CALL read_key_words ( control_file_unit, 'NEB_DETAIL', LEN('NEB_DETAIL'),&
     &                           right, readit )
    dtstart = 1.d0
    imv_cont = 0
    if( readit) then
       READ (right,*) irmethod_1,num_mov,tolforce,Nimage,ak_string_1,&
           & itype_string_1,Etotneb1_1,Etotneb2_1, itype_xatom, f_xatom
       num_mov = num_mov*1.001
       imv_cont = imv_cont*1.001
       Etotneb1_1=Etotneb1_1/Hartree_ev
       Etotneb2_1=Etotneb2_1/Hartree_ev
       ak_string_1 = ak_string_1/Hartree_eV* A_AU_1
       write (long_control_file_unit, 1155) irmethod_1,num_mov,tolforce,Nimage,&
     & ak_string_1,itype_string_1,Etotneb1_1,Etotneb2_1,itype_xatom, trim(f_xatom)
1155   format("     NEB_DETAIL =",&
     &        i3,x,i5,x,f12.5,x,i3,x,f12.5,x,i3,x,2f12.5,x,i3,x,A)
    else
        if (job=='NEB') then
            message="for JOB=NEB, need NEB_DETAIL line"
            call error_message(message)
        end if
    endif
    !
    !*******************************
    !
    call read_key_words (control_file_unit, 'CONVERGE', len('CONVERGE'), right, readit)
    if (readit) then
        temp_right = trim(adjustl(right))
        call transform_to_upper (temp_right, right)
        read (right, '(a10)') converge
    else
        converge = 'DIFFICULT'
    end if
    if (converge /= 'EASY' .and. converge /= 'DIFFICULT') then
        write (6, *) 'CONVERGE is set to a wrong value.'
        write (6, *) 'change it to default value "DIFFICULT".'
        converge = 'DIFFICULT'
    end if
    write (long_control_file_unit, '(a)') '     CONVERGE  = '//trim(converge)
    !
    !
    !
    CALL read_key_words(control_file_unit,'ACCURACY',LEN('ACCURACY'),right,readit)
    IF(readit) then
        temp_RIGHT = adjustl(trim(right))
        call transform_to_upper(temp_right, right)
        read (right, '(A4)')  ACCURACY
    ELSE
        ACCURACY = "NORM"
    ENDIF
    if(ACCURACY .ne. "HIGH" .and. ACCURACY .ne. "NORM") then
        write(*,*) "ACCURACY  is set to a wrong value."
        write(*,*) "change it to default value 'NORM'."
        ACCURACY = "NORM"
    endif
    write (long_control_file_unit, '(a)') '     ACCURACY  = '//trim(ACCURACY)
    !
    !*******************************
    !
    ALLOCATE (pp(sum_atom_type))
    ALLOCATE (atom_pp(sum_atom))
    ALLOCATE ( upfpsp(sum_atom_type) )
      !
    DO i = 1, sum_atom_type
        WRITE (tmp_char, *) i
        tmp_char = 'IN.PSP'//TRIM(ADJUSTL(tmp_char))
        CALL read_key_words ( control_file_unit, TRIM(tmp_char), LEN_TRIM(tmp_char), right, readit)
        IF (readit) THEN
            right = ADJUSTL(right)
            READ (right, '(A200)') pp(i)%file_name
            pp(i)%file_name = ADJUSTL(pp(i)%file_name)
            inquire (file = trim(pp(i)%file_name), exist = alive)
            if (.not. alive) then
                write (6, *) trim(pp(i)%file_name)//' not exist!'
                stop
            end if
            IF (INDEX(pp(i)%file_name, 'vwr') == 1 .and. &
		.not.(index(pp(i)%file_name,"UPF") > 0 .or. index(pp(i)%file_name,"upf") > 0)) THEN
               CALL read_vwr_head(TRIM(pp(i)%file_name), pp(i)%ecut, pp(i)%mass, pp(i)%number, pp(i)%zval)
	    ELSEIF (INDEX(pp(i)%file_name, 'uspp') == 1 .and. &
		    .not.(index(pp(i)%file_name,"UPF") > 0 .or. index(pp(i)%file_name,"upf") > 0)) THEN
            !ELSEIF (INDEX(pp(i)%file_name, 'uspp') == 1) THEN
               CALL read_uspp_head(TRIM(pp(i)%file_name), pp(i)%ecut, pp(i)%mass, pp(i)%number, pp(i)%zval)
            ELSE
      ! iflag=1, means change vloc,rho_atc !
          call read_upf_v2 (trim(pp(i)%file_name), upfpsp(i), 1 )
	       pp(i)%ecut = upfpsp(i)%ecutwfc
	       pp(i)%mass = upfpsp(i)%mass
	       pp(i)%number = upfpsp(i)%num
	       pp(i)%zval = upfpsp(i)%zp
            ENDIF
            CALL gen_element_name (pp(i)%number, pp(i)%element)
         ELSE
	    call error_message &
	    ("In reading "//trim(tmp_char)//": no such item")
         ENDIF
         write(long_control_file_unit, "(A)") '     '//TRIM(tmp_char)//'    = '//TRIM(right)
      ENDDO
      !
      !  from pp info to atom pp info
      !
      loopi: DO i =  1, sum_atom
         loopj: DO j = 1, sum_atom_type
            IF (pp(j)%number /= atom(i)%number) THEN
               IF ( j == sum_atom_type ) THEN
	          call error_message &
		  ("The atomic number from pseudopotential is different from "//trim(atom_file))
               ENDIF
            ELSE
               atom_pp(i) = pp(j)
               EXIT loopj
            ENDIF
         ENDDO loopj
         CALL gen_element_name (atom_pp(i)%number, atom_pp(i)%element)
      ENDDO loopi
      !
      !-----------------------------------------------------------------
      !
      !IF (INDEX(pp(1)%file_name, 'uspp') >= 1) THEN
      IF ( upfpsp(1)%tvanp .or. index(pp(1)%file_name,'uspp') > 0 ) then
          CALL read_key_words( control_file_unit, "QIJ_PD",LEN('QIJ_PD'),right,readit)
          IF(readit ) then
            read (right, *)  IQIJ_PD 
          ELSE
            IQIJ_PD = 0
          ENDIF
          write(long_control_file_unit,1923) IQIJ_PD
1923      format("     QIJ_PD    = ",i2)
          !QIJ_PD_1 = IQIJ_PD
      ENDIF
      !
      !*********************************
      !
      allocate (LDAU_l(sum_atom_type), is_LDAU_l(sum_atom_type), LDAU_nwfc(sum_atom_type), &
          nref_type(sum_atom_type), nref_wfc_type(sum_atom_type), Hubbard_U(sum_atom_type), &
          LDAU_occ(sum_atom_type))
      LDAU_nwfc = 0
      LDAU = 0
      
      do i = 1, sum_atom_type
          !temp_char = 'LDAU_PSP'//char(48+i)
          write(temp_char,*) i
          temp_char="LDAU_PSP"//ADJUSTL(trim(temp_char))
          temp_char = adjustl(temp_char)
          call read_key_words (control_file_unit, trim(temp_char), len_trim(temp_char), right, readit)
          if (readit) then
              LDAU = 1
              is_LDAU_l(i) = 1
              read (right, *) LDAU_l(i), Hubbard_U(i)
              Hubbard_U(i) = Hubbard_U(i)/Hartree_eV
              write (long_control_file_unit, *) '   '//trim(temp_char)//'  = ', &
                  LDAU_l(i), Hubbard_U(i)*Hartree_ev, ' eV'
              if ( LDAU_l(i) < 0) then
                  Hubbard_U(i) = 1.0e-10
                  is_LDAU_l(i) = 0
                  LDAU = 0
                  LDAU_l(i) = -1
              end if
          else
              Hubbard_U(i) = 1.0e-10
              is_LDAU_l(i) = 0
              LDAU_l(i) = -1
              write (long_control_file_unit, *) '    '//trim(temp_char)//' = ',LDAU_l(i), Hubbard_U(i)*Hartree_eV, ' eV'
          end if
      end do
      !
      !*********************************
      !
      call read_key_words (control_file_unit, 'RMASK_LDAU', LEN('RMASK_LDAU'), &
          right, readit)
      if (readit) then
          read (right, *) right_logical, scaler_ldau
          if (right_logical) then
              irmask_ldau = 1
              write (long_control_file_unit, *) '    RMASK_LDAU    = ', right_logical, scaler_ldau
          else
              irmask_ldau = 0
              write (long_control_file_unit, *) '    RMASK_LDAU    = ', right_logical, scaler_ldau
          end if
      else
          irmask_ldau = 0
          scaler_ldau = 1.0
          write (long_control_file_unit, *) '    RMASK_LDAU    = ', .false., 1.0
      end if
      !
      !*********************************
      !
      do i = 1, sum_atom_type
          nref_type(i) = 0
          if (job /= 'DOS') then
              do j = 1, upfpsp(i)%nbeta
                  nref_type(i) = nref_type(i) + upfpsp(i)%lll(j)*2 + 1
              end do
              if (is_LDAU_l(i) == 1) then
                  if (upfpsp(i)%nwfc > 0) then
                      ldau_flag = 0
                      do j = 1, upfpsp(i)%nwfc
                          if (LDAU_l(i) == upfpsp(i)%lchi(j)) then
                              ldau_flag = 1
                              LDAU_nwfc(i) = LDAU_nwfc(i) + 1
                              nref_wfc_type(i) = nref_wfc_type(i) + upfpsp(i)%lchi(j)*2 + 1
                              LDAU_occ(i) = upfpsp(i)%oc(j)
                          end if
                      end do
                      if (ldau_flag == 0) then
                          message = 'For LDA+U, pseudopotential file does not &
         contain the needed atomic wave functions for l='//char(48+ldau_l(i))
                          message = upfpsp(i)%psd//trim(message)
                          write (6, *) trim(message)
                          stop
                      end if
                  else
                      message = 'For LDA+U, pseudopotential file does not &
         contain the needed atomic wave functions for l='//char(48+ldau_l(i))
                      message = upfpsp(i)%psd//trim(message)
                      write (6, *) trim(message)
                      stop
                  end if
              end if
          else
              if (upfpsp(i)%nwfc > 0) then
                  do j = 1, upfpsp(i)%nwfc
                      nref_type(i) = nref_type(i) + upfpsp(i)%lchi(j)*2 + 1
                  end do
              else
                  nref_type(i) = 9
              end if
          end if
      end do
      !
      !*********************************
      !
      call read_key_words (control_file_unit, 'MIX', len('MIX'), right, readit)
      if (readit) then
          read (right, '(a200)') temp_char
          call transform_to_upper(temp_char, mixing)
          if (index(mixing, 'CHA') > 0) then
              potential_mix = .false.
          elseif (index(mixing, 'POT') > 0) then
              potential_mix = .true.
          else
              write (6, *) 'MIX should be CHARGE/POTENTIAL'
              stop
          end if
      else
          potential_mix = .false.
          mixing = 'CHARGE'
      end if
      write (long_control_file_unit, *) '    MIX      = ', trim(mixing)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "Ecut", LEN('Ecut'), right, readit )
      if(readit) then
         read ( right, * ) Ecut
      else
         Ecut=-1.d0
         do i = 1, sum_atom_type
            IF (ecut < pp(i)%ecut) Ecut = pp(i)%ecut  
         enddo
         if (ACCURACY=='HIGH') then
             if (upfpsp(1)%tvanp) then
                 Ecut = Ecut * 1.1
             else
                 Ecut = Ecut * 1.2
             end if
         end if
      endif
      WRITE (tmp_char, '(F8.3)') Ecut
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    ECUT    = '//TRIM(tmp_char)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "Ecut2", LEN('Ecut2'), right, readit)
      if(readit) then
         READ ( right, * ) Ecut2
      else
         Ecut2 = 2*Ecut
         if (ACCURACY=='HIGH') Ecut2 = 4.0* Ecut
      endif
      WRITE (tmp_char, '(F8.3)') Ecut2
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    ECUT2     = '//TRIM(tmp_char)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "Ecut2L", LEN('Ecut2L'), right, readit )
      if(readit) then
         READ ( right, * ) Ecut2L
      else
         if(INDEX(pp(1)%file_name,'vwr') > 1 .or. (.not. upfpsp(1)%tvanp)) then
            Ecut2L = Ecut2
         else
            Ecut2L = 4*Ecut2
         endif
      endif
      WRITE (tmp_char, '(F8.3)') Ecut2L
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    ECUT2L    = '//TRIM(tmp_char)
      !
      !*********************************
      !  change Ecut, Ecut2, Ecut2L to A.U
      !
      Ecut = Ecut / 2.
      Ecut2 = Ecut2 / 2.
      Ecut2L = Ecut2L / 2.
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "N123", LEN('N123'), right, readit)
      if(readit) then
         READ ( right, * ) n1, n2, n3
      else
         fackpt=dsqrt(2.d0*Ecut2)/(4*datan(1.d0))
         dd11=fackpt*dsqrt(latt_vec_bohr(1,1)**2+latt_vec_bohr(2,1)**2+latt_vec_bohr(3,1)**2)
         dd22=fackpt*dsqrt(latt_vec_bohr(1,2)**2+latt_vec_bohr(2,2)**2+latt_vec_bohr(3,2)**2)
         dd33=fackpt*dsqrt(latt_vec_bohr(1,3)**2+latt_vec_bohr(2,3)**2+latt_vec_bohr(3,3)**2)
         call get_n123(n1,dd11) 
         call get_n123(n2,dd22)
         call get_n123(n3,dd33)
      endif
      WRITE (tmp_char, '(3I4)') n1, n2, n3
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    N123    = '//TRIM(tmp_char)
      !
      !-----------------------------------------------------------------
      !
      ALX = sqrt( latt_vec_bohr(1,1)**2 + latt_vec_bohr(2,1)**2 + latt_vec_bohr(3,1)**2) * A_AU_1
      ALY = sqrt( latt_vec_bohr(1,2)**2 + latt_vec_bohr(2,2)**2 + latt_vec_bohr(3,2)**2) * A_AU_1
      ALZ = sqrt( latt_vec_bohr(1,3)**2 + latt_vec_bohr(2,3)**2 + latt_vec_bohr(3,3)**2) * A_AU_1
      
       standard_al = 6.0
       factor = 0.80
       N1S = N1
       N2S = N2
       N3S = N3
      
       if((2*standard_al).le.(factor*ALX)) N1S = (standard_al*N1) / ALX
       if((2*standard_al).le.(factor*ALY)) N2S = (standard_al*N2) / ALY
       if((2*standard_al).le.(factor*ALZ)) N3S = (standard_al*N3) / ALZ
      
       call get_n456(N1S)
       call get_N456(N2S)
       call get_N456(N3S)
      
       if(N1S .gt. N1) N1S = N1
       if(N2S .gt. N2) N2S = N2
       if(N3S .gt. N3) N3S = N3
      
      CALL read_key_words( control_file_unit, 'NS123', LEN('NS123'), right, readit)
      if(readit) then
         READ ( right, *) n1s, n2s, n3s
      endif
      WRITE (tmp_char, '(3I4)') n1s, n2s, n3s
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    NS123    = '//TRIM(tmp_char)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "N123L", LEN('N123L'), right, readit)
      if(readit) then  ! there is a input from etot.input
         READ ( right, * ) n1L, n2L, n3L
      else       ! there is no input from etot.input, use default
         !if(INDEX(pp(1)%file_name, 'vwr').eq.1) then
	 if ( .not. upfpsp(1)%tvanp .or. (index(pp(1)%file_name,'vwr')>0)) then
           n1L=n1
           n2L=n2
           n3L=n3
         else     ! ultrasoft
           n1L=2*n1
           n2L=2*n2
           n3L=2*n3
         endif
      endif
      WRITE (tmp_char, '(3I4)') n1L, n2L, n3L
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    N123L    = '//TRIM(tmp_char)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "SPIN", LEN('SPIN'), right, readit )
      if(readit) then
         READ ( right, * ) spin
      else
         spin=1
      endif
      WRITE (tmp_char, '(I4)') spin
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    SPIN    = '//TRIM(tmp_char)
      !
      !-----------------------------------------------------------------
      !
!      CALL read_key_words ( control_file_unit, "XCFUNCTIONAL", LEN('XCFUNCTIONAL'), right, readit )
!      if(readit) then
!         !READ ( right, * ) xgga
!	 read ( right, * ) xctype(1)
!         if ( xctype(1) < 100 ) then
!         xgga = xctype(1)
!         igga=xgga*1.00001
!         xgga=xgga-igga       ! xgga will be used as a LDA mixing parameter for small density region
!         else
!	     read ( right, * ) xctype(1:3)
!         endif
!      else
!         igga=0
!         xgga=0.d0
!      endif
!      if ( xgga < 100 ) then
!      WRITE (tmp_char, * ) xctype(1)
!      else
!      write (tmp_char, * ) xctype(1:3)
!      endif
!      tmp_char = ADJUSTL(tmp_char)
!      write (long_control_file_unit, *) '    XCFUNCTIONAL    = '//TRIM(tmp_char)
      !
      !-----------------------------------------------------------------
      !
      iflag_hse = 0
      use_hse = 0
      use_hse_thisEtot = 0
      use_sbfft_hse = 1
      CALL read_key_words ( control_file_unit,"XCFUNCTIONAL",LEN('XCFUNCTIONAL'),right,readit )
       if (readit) then
           read ( right, '(A200)' ) temp_char_xc
           temp_char_xc = ADJUSTL(temp_char_xc)
           call capital(temp_char_xc, len(temp_char_xc))
           tmp_char = temp_char_xc
           if ( index(temp_char_xc,"LDA") == 1 ) then
                xctype(1) = 100
                xctype(2) = 1
                xctype(3) = 9
           elseif ( index(temp_char_xc,"SX") == 1 ) then
                hyb_gga_type = 0
                call read_key_words (control_file_unit, 'HSE_OMEGA', len('HSE_OMEGA'), right, readit)
                if (readit) then
                    read (right, *) hse_omega
                else
                    hse_omega = 0.1058
                end if
                write (long_control_file_unit, *) '    HSE_OMEGA = ', hse_omega
                call read_key_words (control_file_unit, 'HSE_ALPHA', len('HSE_ALPHA'), right, readit)
                if (readit) then
                    read (right, *) hse_alpha
                else
                    hse_alpha = 0.25
                end if
                write (long_control_file_unit, *) '    HSE_ALPHA = ', hse_alpha
                call read_key_words (control_file_unit, 'NQ123', len('NQ123'), right, readit)
                if (readit) then
                    read (right, *) hse_nq1, hse_nq2, hse_nq3
                else
                    call read_key_words(control_file_unit, 'MP_N123', len('MP_N123'), right, readit)
                    if (readit) then
                        read (right, *) hse_nq1, hse_nq2, hse_nq3
                    else
                        hse_nq1 = 1
                        hse_nq2 = 1
                        hse_nq3 = 1
                    end if
                end if
                write (long_control_file_unit, *) '    NQ123 = ',hse_nq1, hse_nq2, hse_nq3 
                use_hse = 1
                use_hse_thisEtot = 1
                if (use_sbfft_hse == 1) then
                    fackpt=dsqrt(2.0d0*Ecut)/(4.0d0*datan(1.0d0))
                    dd11=fackpt*dsqrt(latt_vec_bohr(1,1)**2+latt_vec_bohr(2,1)**2+latt_vec_bohr(3,1)**2)
                    dd22=fackpt*dsqrt(latt_vec_bohr(1,2)**2+latt_vec_bohr(2,2)**2+latt_vec_bohr(3,2)**2)
                    dd33=fackpt*dsqrt(latt_vec_bohr(1,3)**2+latt_vec_bohr(2,3)**2+latt_vec_bohr(3,3)**2)
                    call get_n123(hse_n1,dd11+2) 
                    call get_n123(hse_n2,dd22+2) 
                    call get_n123(hse_n3,dd33+2) 
                    if(hse_n1 .gt. n1)  hse_n1 = n1
                    if(hse_n2 .gt. n1)  hse_n2 = n2
                    if(hse_n3 .gt. n1)  hse_n3 = n3
                    CALL read_key_words ( control_file_unit, "P123", LEN('P123'), right,readit)
                    if(readit) then
                        READ ( right, * ) hse_n1,hse_n2, hse_n3
                    endif
                    if (ACCURACY == 'HIGH') then
                        hse_n1 = n1
                        hse_n2 = n2
                        hse_n3 = n3
                    end if
                else
                    hse_n1 = n1
                    hse_n2 = n2
                    hse_n3 = n3
                end if
                write (long_control_file_unit, *) '    P123 = ', hse_n1,hse_n2, hse_n3
            elseif ( index(temp_char_xc,"HSE") == 1 ) then
                hyb_gga_type = 0
                call read_key_words (control_file_unit, 'HSE_OMEGA', len('HSE_OMEGA'), right, readit)
                if (readit) then
                    read (right, *) hse_omega
                else
                    hse_omega = 0.1058
                end if
                write (long_control_file_unit, *) '    HSE_OMEGA = ', hse_omega
                call read_key_words (control_file_unit, 'HSE_ALPHA', len('HSE_ALPHA'), right, readit)
                if (readit) then
                    read (right, *) hse_alpha
                else
                    hse_alpha = 0.25
                end if
                write (long_control_file_unit, *) '    HSE_ALPHA = ', hse_alpha
                call read_key_words( control_file_unit, 'NQ123', len('NQ123'), right, readit)
                if (readit) then
                    read (right, *) hse_nq1, hse_nq2, hse_nq3
                else
                    call read_key_words(control_file_unit, 'MP_N123', len('MP_N123'), right, readit)
                    if (readit) then
                        read (right, *) hse_nq1, hse_nq2, hse_nq3
                    else
                        hse_nq1 = 1
                        hse_nq2 = 1
                        hse_nq3 = 1
                    end if
                end if
                write (long_control_file_unit, *) '    NQ123 = ',hse_nq1, hse_nq2, hse_nq3 
                use_hse = 1
                use_hse_thisEtot = 1
                if (use_sbfft_hse == 1) then
                    fackpt=dsqrt(2.0d0*Ecut)/(4.0d0*datan(1.0d0))
                    dd11=fackpt*dsqrt(latt_vec_bohr(1,1)**2+latt_vec_bohr(2,1)**2+latt_vec_bohr(3,1)**2)
                    dd22=fackpt*dsqrt(latt_vec_bohr(1,2)**2+latt_vec_bohr(2,2)**2+latt_vec_bohr(3,2)**2)
                    dd33=fackpt*dsqrt(latt_vec_bohr(1,3)**2+latt_vec_bohr(2,3)**2+latt_vec_bohr(3,3)**2)
                    call get_n123(hse_n1,dd11+2) 
                    call get_n123(hse_n2,dd22+2) 
                    call get_n123(hse_n3,dd33+2) 
                    if(hse_n1 .gt. n1)  hse_n1 = n1
                    if(hse_n2 .gt. n1)  hse_n2 = n2
                    if(hse_n3 .gt. n1)  hse_n3 = n3
                    CALL read_key_words ( control_file_unit, "P123", LEN('P123'), right,readit)
                    if(readit) then
                        READ ( right, * ) hse_n1,hse_n2, hse_n3
                    endif
                    if (ACCURACY == 'HIGH') then
                        hse_n1 = n1
                        hse_n2 = n2
                        hse_n3 = n3
                    end if
                else
                    hse_n1 = n1
                    hse_n2 = n2
                    hse_n3 = n3
                end if
                write (long_control_file_unit, *) '    P123 = ', hse_n1,hse_n2, hse_n3
                xctype(1) = 100
                xctype(2) = 1
                xctype(3) = 9
                iflag_hse = 1
            elseif ( index(temp_char_xc,"BP") == 1 ) then
                xctype(1) = 100
                xctype(2) = 106
                xctype(3) = 132
            elseif ( index(temp_char_xc,"PW91") == 1 ) then
                xctype(1) = 100
                xctype(2) = 109
                xctype(3) = 134
            elseif ( index(temp_char_xc,"PBESOL") == 1 ) then
                xctype(1) = 100
                xctype(2) = 116
                xctype(3) = 133
            elseif ( index(temp_char_xc,"PBE") == 1 ) then
                xctype(1) = 100
                xctype(2) = 101
                xctype(3) = 130
            elseif ( index(temp_char_xc,"M06L") == 1 ) then
                xctype(1) = 100
                xctype(2) = 203
                xctype(3) = 233
            elseif ( index(temp_char_xc,"SOGGA") == 1 ) then
                xctype(1) = 100
                xctype(2) = 151
                xctype(3) = 152
            elseif ( index(temp_char_xc,"Q2D") == 1 ) then
                xctype(1) = 100
                xctype(2) = 47
                xctype(3) = 48
            elseif ( index(temp_char_xc,"PW86PBE") == 1 ) then
                xctype(1) = 100
                xctype(2) = 108
                xctype(3) = 130
            elseif ( index(temp_char_xc,"WC") == 1 ) then
                xctype(1) = 100
                xctype(2) = 118
                xctype(3) = 130
            elseif ( index(temp_char_xc,"TPSS") == 1 ) then
                xctype(1) = 100
                xctype(2) = 202
                xctype(3) = 231
            elseif ( index(temp_char_xc,"LIBXC") > 0 ) then
                xctype(1) = 100
                i_plus = index(temp_char_xc,"+")
                if ( i_plus < 1 ) then
                    read ( temp_char_xc(7:i_plus-1), * ) xctype(2)
                    xctype(3) = xctype(2)
                else
                    read ( temp_char_xc(7:i_plus-1), * ) xctype(2)
                    temp_char_xc = temp_char_xc(i_plus+1:)
                    temp_char_xc = adjustl(temp_char_xc)
                    read ( temp_char_xc(7:), * ) xctype(3)
                end if
            elseif ( index(temp_char_xc,"XC") > 0 ) then
                xctype(1) = 100
                i_plus = index(temp_char_xc,"+")
                !if ( i_plus < 1 ) then
                !    xctype(2) = libxc_id(temp_char_xc)
                !    xctype(3) = 0
                !else
                !    xctype(2) = libxc_id(temp_char_xc(:i_plus-1))
                !    temp_char_xc = temp_char_xc(i_plus+1:)
                !    temp_char_xc = adjustl(temp_char_xc)
                !    xctype(3) = libxc_id(trim(temp_char_xc))
                !end if
                if ( i_plus < 1 ) then
                    xctype(2) = xc_f90_functional_get_number(temp_char_xc(4:))
                    xctype(3) = xctype(2)
                else
                    xctype(2) = xc_f90_functional_get_number(temp_char_xc(4:i_plus-1))
                    temp_char_xc = temp_char_xc(i_plus+1:)
                    temp_char_xc = adjustl(temp_char_xc)
                    xctype(3) = xc_f90_functional_get_number(temp_char_xc(4:))
                end if
            else
                write ( 6, * ) 'XCFUNCTIONAL should be "LDA"/"GGA"'
                write ( 6, * ) 'For libxc: "XC_XXX" or "LIBXC_YYY" '
                stop
            end if
       else
           tmp_char = "LDA"
           xctype(1) = 100
           xctype(2) = 1
           xctype(3) = 9
           igga=0
           xgga=0.d0
       end if
       write (long_control_file_unit, *) "    XCFUNCTIONAL = ", trim(tmp_char)
       igga=0
       xgga=0.d0
       
       if ( xctype(1) == 100 ) call libxc_check(xctype(2),xctype(3),spin)
       !
       CALL read_key_words ( control_file_unit, "HSE_dN", LEN('HSE_dN'), right, readit)
       if(readit) then
           READ ( right, * ) iupdate_sxp_dn
       else
           iupdate_sxp_dn=3
       endif
       write(long_control_file_unit,1124) iupdate_sxp_dn
1124   format("     HSE_dN = ",i5)
       !
       CALL read_key_words ( control_file_unit, "RELAX_HSEM", LEN('RELAX_HSEM'), right, readit)
       if(readit) then
            READ ( right, * ) num_mov_hse,num_mov_LDA,fac_HSELDA
            if( num_mov_hse.gt.0) then
                irelax_hseLDA=1
            else
                irelax_hseLDA=0
            endif
       else
           irelax_hseLDA=1
           num_mov_hse=30
           num_mov_LDA=30
           fac_HSELDA=0.05
       endif
       write(long_control_file_unit,1125) num_mov_hse,num_mov_LDA,fac_HSELDA
1125     format("     RELAX_HSEM = ",2x,i4,1x,i4,1x,E12.5)
       !
       !
       !if(igga.ne.0.and.igga.ne.1.and.inode.eq.1) then
       !write(6,*) "XcFunctional must be 0(no gga) or 1(gga), stop",igga
       !call mpi_abort(MPI_COMM_WORLD,ierr)
       !endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !
      !-----------------------------------------------------------------
      !

       CALL read_key_words ( control_file_unit, "VDW", LEN("VDW"), right, readit )
       if ( readit) then
           read ( right, '(A30)' ) vdw_method
           call capital(vdw_method, len(vdw_method))
           if ( trim(adjustl(vdw_method)) == 'DFT-D' ) then
               has_vdw = .true.
               has_london = .true.
                   write(long_control_file_unit,*) "    VDW      = ", trim(vdw_method)
           else
               has_vdw = .false.
               has_london = .false.
           end if
           if ( has_london ) then
               call  read_key_words ( control_file_unit, "LONDON_S6", len('LONDON_S6'), right, readit )
               if ( readit ) then
                    read ( right , * ) london_s6
                        write(long_control_file_unit,*) "    LONDON_S6      = ", london_s6
                    end if
                    ALLOCATE(LONDON_C6(sum_atom_type))
               do ia = 1, sum_atom_type
                   write ( temp_char, * ) ia
                   temp_char = adjustl(temp_char)
                   tmp_char = 'LONDON_C6('//trim(temp_char)//')'
                   call read_key_words ( control_file_unit, trim(tmp_char), len_trim(tmp_char), right, readit)
                   if ( readit) then
                       read ( right, * ) london_c6(ia)
                           write(long_control_file_unit,*) '   '//trim(tmp_char)//"      = ", london_c6(ia)
                   end if
               end do
               call read_key_words ( control_file_unit, 'LONDON_RCUT', len('LONDON_RCUT'), right, readit )
               if ( readit ) then
                   read ( right , * ) london_rcut
                       write(long_control_file_unit,*) "    LONDON_RCUT      = ", london_rcut
               end if
           end if
       else
           has_vdw = .false.
           write(long_control_file_unit,*) "    VDW      =  NONE"
       end if
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, "SMTH", LEN('SMTH'), right, readit )
      if(readit) then
         READ ( right, * ) Smth
      else
         Smth=1.d0
      endif
      WRITE (tmp_char, '(F8.5)') Smth
      tmp_char = ADJUSTL(tmp_char)
      write (long_control_file_unit, *) '    SMTH    = '//TRIM(tmp_char)
      !
      !*********************************
      !
      CALL read_key_words ( control_file_unit, "COULOMB", LEN('COULOMB'), right, readit )
      if(readit) then
         READ (right, *) coulomb_flag
         SELECT CASE (coulomb_flag)
         CASE (1)
            READ (right, *) coulomb_flag, xcoul(1), xcoul(2), xcoul(3)
            WRITE (long_control_file_unit, "('     COULOMB    = ', i4, 2x, 3(f10.5,1x))")  &
               coulomb_flag, xcoul(1), xcoul(2), xcoul(3)
         CASE (11)
            READ (right, *) coulomb_flag, xcoul(1)
            WRITE (long_control_file_unit, "('     COULOMB    = ', i4, 2x, f10.5)") coulomb_flag,xcoul(1)
         CASE (12)
            READ (right, *) coulomb_flag, xcoul(2)
            WRITE (long_control_file_unit, "('     COULOMB    = ', i4, 2x, f10.5)") coulomb_flag, xcoul(2)
         CASE (13)
            READ (right, *) coulomb_flag, xcoul(3)
            WRITE (long_control_file_unit, "('     COULOMB    = ', i4, 2x, f10.5)") coulomb_flag, xcoul(3)
         CASE DEFAULT
            coulomb_flag = 0
            WRITE (long_control_file_unit, "('     COULOMB    = ', i2)") coulomb_flag
         END SELECT
         coulomb_flag = coulomb_flag * 1.001
      else
         coulomb_flag = 0 
         WRITE (long_control_file_unit, "('     COULOMB    = ', i2)") coulomb_flag
      endif
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, "PWSCF", LEN('PWSCF'), right, readit )
      IF ( readit ) THEN
         READ ( right, * ) pwscf_out
      ELSE
         pwscf_out = .FALSE.
      ENDIF
      WRITE (tmp_char, *) pwscf_out
      tmp_char = ADJUSTL(tmp_char)
      WRITE (long_control_file_unit, *) '    PWSCF    = '//TRIM(ADJUSTL(tmp_char)) 
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'IN.WG', LEN('IN.WG'), right, readit )
      if(readit) then
         read ( right, * ) right_logical
         if ( right_logical ) then
            in_wg = 1
            in_wg_file(1)='IN.WG'
               INQUIRE (FILE = trim(in_wg_file(1)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_wg_file(1))//" doesn't exist")
            if(spin.eq.2) then
               in_wg_file(2)=trim(in_wg_file(1))//"_2"
               INQUIRE (FILE = trim(in_wg_file(2)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_wg_file(2))//" doesn't exist!")
            endif
            write(long_control_file_unit,*) "    IN.WG    = ", right_logical
         else
            in_wg = 0
            write(long_control_file_unit,*) "    IN.WG    = ", right_logical
         endif
      else
         if(INDEX(job,"DOS") >= 1) then
             in_wg = 1
             in_wg_file(1)='IN.WG'
             INQUIRE (FILE = trim(in_wg_file(1)), EXIST = alive)
             IF (.NOT. alive) call error_message (trim(in_wg_file(1))//" doesn't exist")
             if(spin.eq.2) then
                 in_wg_file(2)=trim(in_wg_file(1))//"_2"
                 INQUIRE (FILE = trim(in_wg_file(2)), EXIST = alive)
                 IF (.NOT. alive) call error_message (trim(in_wg_file(2))//" doesn't exist!")
             endif
             write(long_control_file_unit,*) "    IN.WG    = ", .TRUE.

         else
             in_wg = 0
             write(long_control_file_unit,*) "    IN.WG    = ", .FALSE.
         endif
      endif
      !
      !-----------------------------------------------------------------
      !
      IF (pwscf_out.OR.(INDEX(job,"SCF") >= 1).OR.(INDEX(job,"NONSCF")>=1)) THEN   !pwscf
         out_wg = 1
         out_wg_file(1) = 'OUT.WG'
         if(spin.eq.2) then
            out_wg_file(2)=trim(out_wg_file(1))//"_2"
         endif
         right_logical = .TRUE.
         WRITE (tmp_char, *) right_logical
         write(long_control_file_unit,*) "    OUT.WG    = "//TRIM(ADJUSTL(tmp_char))
      ELSE   !pwscf
         CALL read_key_words ( control_file_unit, 'OUT.WG', LEN('OUT.WG'), right, readit )
         if(readit) then
            read ( right, * ) right_logical
            IF ( right_logical ) THEN
               out_wg = 1
               out_wg_file(1) = 'OUT.WG'
               if(spin.eq.2) then
                  out_wg_file(2)=trim(out_wg_file(1))//"_2"
               endif
               write(long_control_file_unit,*) "    OUT.WG    = ", right_logical
            ELSE
               out_wg = 0
               write(long_control_file_unit,*) "    OUT.WG    = ", right_logical
            ENDIF
         else
            if (INDEX(job,"SCF")>=1) then
               out_wg = 1
               out_wg_file(1) = 'OUT.WG'
               if(spin.eq.2) then
                   out_wg_file(2)=trim(out_wg_file(1))//"_2"
               endif
               write(long_control_file_unit,*) "    OUT.WG    = ", .TRUE.
            else
               out_wg = 0
               write(long_control_file_unit,*) "    OUT.WG    = ", .FALSE.
            endif
         endif   
      ENDIF  ! pwscf
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'IN.RHO', LEN('IN.RHO'), right, readit )
      if(readit) then 
         read ( right, * ) right_logical
         IF ( right_logical ) THEN
            in_rho = 1
            in_rho_file(1)='IN.RHO'
               INQUIRE (FILE = trim(in_rho_file(1)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_rho_file(1))//" doesn't exist!")
            if(spin.eq.2) then
               in_rho_file(2)=trim(in_rho_file(1))//"_2"
               INQUIRE (FILE = trim(in_rho_file(2)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_rho_file(2)//" doesn't exist!"))
            endif
            write(long_control_file_unit,*) "    IN.RHO    = ", right_logical
         ELSE
            in_rho = 0
            write(long_control_file_unit,*) "    IN.RHO    = ", right_logical
         ENDIF
      else
         in_rho = 0
         write(long_control_file_unit,*) "    IN.RHO    = ", .FALSE.
      endif
      !
      !-----------------------------------------------------------------
      !
      IF (pwscf_out.OR.(INDEX(job,"SCF")>=1).OR.(INDEX(job,"NONSCF")>=1)) THEN   !pwscf
         out_rho = 1
         out_rho_file(1) = 'OUT.RHO'
         if(spin.eq.2) then
            out_rho_file(2)=trim(out_rho_file(1))//"_2"
         endif
         right_logical = .TRUE.
         WRITE (tmp_char, *) right_logical
         write(long_control_file_unit,*) "    OUT.RHO    = "//TRIM(ADJUSTL(tmp_char))
      ELSE   !pwscf
         CALL read_key_words ( control_file_unit, 'OUT.RHO', LEN('OUT.RHO'), right, readit )
         if(readit) then
            read ( right, * ) right_logical
            IF ( right_logical ) THEN
               out_rho = 1
               out_rho_file(1) = 'OUT.RHO'
               if(spin.eq.2) then
                   out_rho_file(2)=trim(out_rho_file(1))//"_2"
               endif
               write(long_control_file_unit,*) "    OUT.RHO    = ", right_logical
            ELSE
               out_rho = 0 
               write(long_control_file_unit,*) "    OUT.RHO    = ", right_logical
            ENDIF
         else
            if(INDEX(job,"SCF")>=1) then
               out_rho = 1
               out_rho_file(1) = 'OUT.RHO'
               if(spin.eq.2) then
                   out_rho_file(2)=trim(out_rho_file(1))//"_2"
               endif
               write(long_control_file_unit,*) "    OUT.RHO    = ", .TRUE.
            else
                out_rho = 0
                write(long_control_file_unit,*) "    OUT.RHO    = ", .FALSE.
            endif
         endif
      ENDIF   !  pwscf
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit , 'IN.VR', LEN('IN.VR'), right, readit )
      if(readit) then
         read ( right, * ) right_logical
         IF ( right_logical ) THEN
            in_vr = 1
            in_vr_file(1) = 'IN.VR'
               INQUIRE (FILE = trim(in_vr_file(1)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_vr_file(1))//" doesn't exist!")
            if(spin.eq.2) then
               in_vr_file(2)=trim(in_vr_file(1))//"_2"
               INQUIRE (FILE = trim(in_vr_file(2)), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(in_vr_file(2))//" doesn't exist!")
            endif
            write(long_control_file_unit,*) "    IN.VR    = ", right_logical
         ELSE
            in_vr = 0 
            write(long_control_file_unit,*) "    IN.VR    = ", right_logical
         ENDIF
      else
         in_vr=0
         write(long_control_file_unit,*) "    IN.VR    = ", .FALSE.
      endif
      if (job == 'NONSCF' .and. in_vr == 0) call error_message('in NONSCF calculation, please set IN.VR = T in "etot.input!"')
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'OUT.VR', LEN('OUT.VR'), right, readit )
      if(readit) then
         read ( right, * ) right_logical
         IF ( right_logical ) THEN
            out_vr = 2
            out_vr_file(1) = 'OUT.VR'
            if(spin.eq.2) then
               out_vr_file(2)=trim(out_vr_file(1))//"_2"
            endif
            write(long_control_file_unit,*) "    OUT.VR    = ", right_logical
         ELSE
            out_vr = 0
            write(long_control_file_unit,*) "    OUT.VR    = ", right_logical
         ENDIF
      else
           if(INDEX(job,"SCF")>=1) then
               out_vr = 1
               out_vr_file(1) = 'OUT.VR'
               if(spin.eq.2) then
                   out_vr_file(2)=trim(out_vr_file(1))//"_2"
               endif
               write(long_control_file_unit,*) "    OUT.VR     = ", .TRUE.
           else
                out_vr=0
                write(long_control_file_unit,*) "    OUT.VR    = ", .FALSE.
           endif
      endif
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'IN.VEXT', LEN('IN.VEXT'), right, readit )
      if(readit) then
         READ ( right, * ) right_logical
         IF ( right_logical ) THEN
            in_vext = 1
            in_vext_file = 'IN.VEXT'
               INQUIRE (FILE = trim(in_vext_file), EXIST = alive) 
               IF (.NOT. alive) call error_message (trim(in_vext_file)//" doesn't exist!")
            write(long_control_file_unit,*) "    IN.VEXT    = ", right_logical
         ELSE
            in_vext = 0
            write(long_control_file_unit,*) "    IN.VEXT    = ", right_logical
         ENDIF
      else
         in_vext =0
         write(long_control_file_unit,*) "    IN.VEXT    = ", .FALSE.
      endif
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'OUT.RHO_SP', LEN('OUT.RHO_SP'), right, readit )
      if(readit) then
         READ (right, *) flag_dens
         SELECT CASE (flag_dens)
         CASE (0)
            out_dens=0
            WRITE (tmp_char, *) flag_dens
            WRITE (long_control_file_unit, *)  '    OUT.RHO_SP    = '//TRIM(ADJUSTL(tmp_char))
         CASE (1, 11, 2, 22)    
            READ ( right, * ) flag_dens, kpt_dens(1), kpt_dens(2), &
     &      spin_dens(1), spin_dens(2), wg_dens(1), wg_dens(2) 
            out_dens=1
            out_dens_file = "OUT.RHO_SP"
            write(tmp_char, 2002) &
     &               flag_dens, kpt_dens(1),kpt_dens(2),&
     &                spin_dens(1),spin_dens(2),&
     &                wg_dens(1),wg_dens(2)
2002   format(i3,3(i5,1x,i5,4x)) 
            WRITE (long_control_file_unit, *) '    OUT.RHO_SP    = '//TRIM(ADJUSTL(tmp_char))
         CASE default
            call error_message ("OUT.RHO_SP should be 0/1/11/2/22?")
         END SELECT
      else
         out_dens=0
         WRITE (tmp_char, *) 0
         WRITE (long_control_file_unit, *) '    OUT.RHO_SP    = '//TRIM(ADJUSTL(tmp_char))
      endif
      !
      !-----------------------------------------------------------------
      !
      IF ( TRIM( JOB ) .eq. "RELAX" .OR. TRIM(JOB) .eq. "MD" ) THEN
         OUT_FORCE = 1
         out_force_file = "OUT.FORCE"
         write ( long_control_file_unit, * ) "    OUT.FORCE    = ", .TRUE.
      ELSE
         CALL read_key_words ( control_file_unit, 'OUT.FORCE', LEN('OUT.FORCE'), right, readit )
         if ( readit ) then
            READ ( right, * ) right_logical 
            IF ( right_logical ) THEN
               OUT_FORCE = 1
               out_force_file = 'OUT.FORCE'
               write(long_control_file_unit,*) "    OUT.FORCE    = ", right_logical
            ELSE
               OUT_FORCE = 0
               write(long_control_file_unit,*) "    OUT.FORCE    = ", right_logical
            ENDIF
         else
            OUT_FORCE = 0
            write(long_control_file_unit,*) "    OUT.FORCE    = ", .FALSE.
         endif
      ENDIF
      !
      !-----------------------------------------------------------------
      !
      IF ( TRIM( JOB ) .eq. "RELAX") THEN
          if(istress_cal .eq. 1) then
              write(long_control_file_unit,*) "    OUT.STRESS = ", .TRUE.
          else
              write(long_control_file_unit,*) "    OUT.STRESS = ", .FALSE.
          endif
      ELSEIF ( TRIM( JOB ) .ne. "SCF") THEN
          write(long_control_file_unit,*) "    OUT.STRESS = ", .FALSE.
      ENDIF
      
      IF ( TRIM( JOB ) .eq. "SCF") THEN
           CALL read_key_words ( control_file_unit, 'OUT.STRESS',LEN('OUT.STRESS'), right, readit)
           if (readit) then
               READ ( right, * ) right_logical 
               IF ( right_logical ) THEN
                   istress_cal = 1
                   write(long_control_file_unit,*) "    OUT.STRESS = ", right_logical
               ELSE
                   istress_cal = 0
                   write(long_control_file_unit,*) "    OUT.STRESS = ", right_logical
               ENDIF
           else
               istress_cal = 0
               write(long_control_file_unit,*) "    OUT.STRESS = ", .FALSE.
           endif
      ENDIF
      !
      !-----------------------------------------------------------------
      !
      out_kpt = .False.
      out_symm = .False.
      CALL read_key_words ( control_file_unit, "MP_N123", LEN('MP_N123'), right, readit )
      IF ( readit ) THEN
         READ ( right, *) nk1, nk2, nk3, sk1, sk2, sk3
         WRITE (tmp_char, '(6i4)') nk1, nk2, nk3, sk1, sk2, sk3
         WRITE (long_control_file_unit, *) '    MP_N123     = '//TRIM(ADJUSTL(tmp_char))
         WRITE (tmp_char, *) .TRUE.
         write(long_control_file_unit,*) "    IN.SYMM    = "//TRIM(ADJUSTL(tmp_char))
         WRITE (tmp_char, *) .TRUE.
         write(long_control_file_unit,*) "    IN.KPT    = ", TRIM(ADJUSTL(tmp_char))
         out_kpt = .true.
         out_symm = .true.
      ELSE
         CALL read_key_words ( control_file_unit, 'IN.SYMM', LEN('IN.SYMM'), right, readit )
         if(readit) then
            READ ( right, * ) right_logical
            IF ( right_logical ) THEN
               in_symm = 1
               sym_file = 'IN.SYMM'
               INQUIRE (FILE = trim(sym_file), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(sym_file)//" doesn't exist!")
               WRITE (tmp_char, *) right_logical
               write(long_control_file_unit,*) "    IN.SYMM    = "//TRIM(ADJUSTL(tmp_char))
            ELSE
               in_symm = 0
               WRITE (tmp_char, *) right_logical
               write(long_control_file_unit,*) "    IN.SYMM    = "//TRIM(ADJUSTL(tmp_char))
            ENDIF
         else
            in_symm = 0
            WRITE (tmp_char, *) .FALSE.
            write(long_control_file_unit,*) "    IN.SYMM    = "//TRIM(ADJUSTL(tmp_char))
         endif
         !
         !
         CALL read_key_words ( control_file_unit, 'IN.KPT', LEN('IN.KPT'), right, readit )
         if(readit) then
            READ ( right, * ) right_logical 
            IF ( right_logical ) THEN
               in_kpt = 1
               kpt_file = 'IN.KPT'
               INQUIRE (FILE = trim(kpt_file), EXIST = alive)
               IF (.NOT. alive) call error_message (trim(kpt_file)//" doesn't exist")
               WRITE (tmp_char, *) right_logical
               write(long_control_file_unit,*) "    IN.KPT    = ", TRIM(ADJUSTL(tmp_char))
            ELSE
               in_kpt= 0
               WRITE (tmp_char, *) right_logical
               write(long_control_file_unit,*) "    IN.KPT    = ", TRIM(ADJUSTL(tmp_char))
            ENDIF
         else
            in_kpt= 0
            WRITE (tmp_char, *) .FALSE.
            write(long_control_file_unit,*) "    IN.KPT    = ", TRIM(ADJUSTL(tmp_char))
         endif
         !
      ENDIF
      !
      !-----------------------------------------------------------------
      !
      totNel = 0.0
      DO i = 1, sum_atom
         totNel = totNel + atom_pp(i)%zval
      ENDDO
      !
      CALL read_key_words ( control_file_unit, 'NUM_ELECTRON', LEN('NUM_ELECTRON'), right, readit )
      if(readit) then
         READ ( right, * ) totNel     ! otherwise, use the totNel calculated above
      endif
      WRITE (tmp_char, '(f12.5)') totNel
      write(long_control_file_unit,*) "    NUM_ELECTRON    = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'NUM_BAND', LEN('NUM_BAND'), right, readit )
      if(readit) then
        read ( right, * ) mx
      else
        mx = 1.05*totNel/2+10
      endif
      WRITE (tmp_char, *) mx
      write(long_control_file_unit,*) "    NUM_BAND    = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'WG_ERROR', LEN('WG_ERROR'), right, readit )
      if(readit) then
        READ ( right, * ) ug_thred
      else
        ug_thred = 1.D-4
        if (CONVERGE=='DIFFICULT') ug_thred = ug_thred*0.5
      endif
      WRITE (tmp_char, '(E12.5)') ug_thred
      write(long_control_file_unit,*) "    WG_ERROR    = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'RHO_ERROR', LEN('RHO_ERROR'), right, readit )
      if(readit) then
        READ ( right, * ) tolrho
      else
        tolrho = 0.4D-3
        if (CONVERGE=='DIFFICULT') tolrho = tolrho * 0.5
      endif
      WRITE (tmp_char, '(E12.5)') tolrho
      write(long_control_file_unit,*) "    RHO_ERROR = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'E_ERROR', LEN('E_ERROR'), right, readit )
      if(readit) then
        READ ( right, * ) e_thred
      else
        e_thred = 1.D-7*totNel* Hartree_eV
        if (CONVERGE=='DIFFICULT') e_thred = e_thred*0.01
      endif
      WRITE (tmp_char, '(E12.5)') e_thred
      write(long_control_file_unit,*) "    E_ERROR = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'W_CG', LEN('W_CG'), right, readit )
      if(readit) then
        READ ( right, * ) w_cg
      else
        w_cg = 0.07d0    ! this is a test result
        if (CONVERGE=='DIFFICULT') w_cg = 0.0
      endif
      WRITE (tmp_char, '(E12.5)') w_cg
      write(long_control_file_unit,*) "    W_CG = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'W_SCF', LEN('W_SCF'), right, readit )
      if(readit) then
        READ ( right, * ) w_scf
      else
        if (trim(JOB).eq. "MD") then
           w_scf = 0.02d0       ! this is a test result
        elseif (trim(JOB).eq."RELAX") then
           w_scf = 0.003d0       ! this is a test result
        endif
        if (CONVERGE=='DIFFICULT') w_scf = w_scf * 0.05
      endif
      WRITE (tmp_char, '(E12.5)') w_scf
      write(long_control_file_unit,*) "    W_SCF = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'QijL0_GS', LEN('QijL0_GS'), right, readit )
      if(readit) then
        READ ( right, * ) iQijL0_GS
      else
        iQijL0_GS = 1     ! Gspace implementation, Real-space: 2
      endif
      WRITE (tmp_char, *) iQijL0_GS
      write(long_control_file_unit,*) "    QijL0_GS = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words( control_file_unit, 'SBFFT', LEN('SBFFT'), right, readit )
      if ( readit ) then
        READ( right, * ) isbf
      else
        isbf = 1
      endif
      WRITE (tmp_char, *) isbf
      write(long_control_file_unit,*) "    SBFFT = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      CALL read_key_words ( control_file_unit, 'SCF_ITER0', LEN('SCF_ITER0'), right, readit )
      if ( readit  ) then
        READ ( right, * ) num_iter0, num_line0
        scf_iter0=1
        cg_bad0=5
      else
        scf_iter0=0
        if (CONVERGE=='EASY') then
            num_iter0 = 40
        else
            num_iter0 = 100
        end if
        num_line0 = 6
        cg_bad0 = 5
        if ( trim ( JOB ) .eq. "NONSCF" ) num_iter0=10
      endif
      WRITE (tmp_char, *) num_iter0, num_line0
      write(long_control_file_unit,*) "    SCF_ITER0    = "//TRIM(ADJUSTL(tmp_char))
      !
      !-----------------------------------------------------------------
      !
      ALLOCATE (fermi_de0(num_iter0), cg_mth0(num_iter0), type_fermi0(num_iter0),&
             mx_mtn0(num_iter0), scf_mth0(num_iter0))
      if(scf_iter0.eq.1) then
        CALL read_key_words ( control_file_unit, 'ALGORITHM0', LEN('ALGORITHM0'), right, readit )
        if(.not. readit) then
           message="Not set ALGORITHM0 in ""etot.input""! "
           call error_message (message)
        else
           READ ( right, * ) cg_mth0(1), scfmth, fermi_de0(1), type_fermi0(1)
           cg_mth0(1) = cg_mth0(1)*1.0001
           scf_mth0(1)=scfmth+1.D-6
           mx_mtn0(1)=scfmth-scf_mth0(1)
           fermi_de0(1)=fermi_de0(1)/27.211396d0
           do i=2,num_iter0
               read(control_file_unit,*,iostat=ierr) cg_mth0(i),scfmth, fermi_de0(i),type_fermi0(i)
               if(ierr.ne.0) then
                   message="Need explicit num_iter0 ALGORITHM0 lines for scf_iter0 input xx"//char(48+i)
                   call error_message(message)
               endif
               scf_mth0(i)=scfmth+1.D-6
               mx_mtn0(i)=scfmth-scf_mth0(i)
               fermi_de0(i)=fermi_de0(i)/27.211396d0
           enddo
        endif
      else     ! scf_iter0
        do i=1,num_iter0
           cg_mth0(i)=3
           fermi_de0(i)=0.1d0/27.211396d0
           scf_mth0(i)=1
           mx_mtn0(i)=0.d0
           type_fermi0(i)=1
        enddo
        scf_mth0(1)=0
        scf_mth0(2)=0
        scf_mth0(3)=0
        if(trim(JOB).eq."NONSCF") then
           do i=1,num_iter0
               scf_mth0(i)=0
           enddo
        endif
      endif
2003       format("     ALGORITHM0    = ",i4,2x,f8.4,2x,f10.5,2x,i3) 
2004       format("                     ",i4,2x,f8.4,2x,f10.5,2x,i3) 
      write(long_control_file_unit,2003)  cg_mth0(1),scf_mth0(1)+mx_mtn0(1), &
                   fermi_de0(1)*27.211396d0,type_fermi0(1)
      do i=2,num_iter0
        write(long_control_file_unit,2004)  cg_mth0(i),scf_mth0(i)+mx_mtn0(i), &
                       fermi_de0(i)*27.211396d0,type_fermi0(i)
      enddo
      !
      !
      if(trim(JOB).eq."DOS") then
        do_dos=1
      else
        do_dos=0
      endif
      if(do_dos.eq.1.and.in_wg.eq.0) then
        message="  IN.WG = T, when JOB = DOS"
        call error_message (message)
      endif
      !
      !
      !
      !
      CALL read_key_words ( control_file_unit, 'SCF_ITER1', LEN('SCF_ITER1'), &
                         right, readit ) 
      if(readit) then
        READ ( right, * ) num_iter1, num_line1
        scf_iter1=1
        cg_bad1=5
      else
        scf_iter1=0
        if (CONVERGE=='EASY') then
            num_iter1 = 40
        else
            num_iter1 = 50
        end if
        num_line1 = 4
        cg_bad1 = 5
      endif
      WRITE (tmp_char, *) num_iter1,num_line1
      write(long_control_file_unit,*) "    SCF_ITER1    = "//TRIM(ADJUSTL(tmp_char))
      !
      !    
      ALLOCATE (fermi_de1(num_iter1), cg_mth1(num_iter1), type_fermi1(num_iter1),&
                mx_mth1(num_iter1), scf_mth1(num_iter1))
      if(scf_iter1.eq.1) then
        CALL read_key_words ( control_file_unit, 'ALGORITHM1', LEN('ALGORITHM1'), right, readit )
        if(.not.readit) then
            message= "    ALGORITHM1 = xxx, xxx, xxx, xxx?"
            call error_message (message)
        else
            READ ( right, * ) cg_mth1(1), scfmth,  &
                              fermi_de1(1), type_fermi1(1)
            cg_mth1(1) = cg_mth1(1)*1.001
            ! 
            scf_mth1(1)=scfmth+1.D-6
            mx_mth1(1)=scfmth-scf_mth1(1)
            fermi_de1(1)=fermi_de1(1)/27.211396d0
            do i=2,num_iter1
                read(control_file_unit,*,iostat=ierr) cg_mth1(i),scfmth,&
                                           fermi_de1(i),type_fermi1(i)
                if(ierr.ne.0) then
                message= &
                "Need explicit num_iter1 ALGORITHM1 lines for scf_iter1 input xx"&
                    //char(48+i)
                call error_message (message)    
                endif
                scf_mth1(i)=scfmth+1.D-6
                mx_mth1(i)=scfmth-scf_mth1(i)
                fermi_de1(i)=fermi_de1(i)/27.211396d0
            enddo
            write(long_control_file_unit,2013)  cg_mth1(1),scf_mth1(1)+mx_mth1(1),&
                            fermi_de1(1)*27.211396d0,type_fermi1(1)
           do i=2,num_iter1
               write(long_control_file_unit,2014)  cg_mth1(i),scf_mth1(i)+mx_mth1(i), &
                               fermi_de1(i)*27.211396d0,type_fermi1(i)
           enddo
       endif
     else     ! scf_iter1
       do i=1,num_iter1
           cg_mth1(i)=3
           fermi_de1(i)=0.1d0/27.211396d0
           scf_mth1(i)=1
           mx_mth1(i)=0.d0
           type_fermi1(i)=1
       enddo
2013   format("     ALGORITHM1    = ",i4,2x,f8.4,2x,f10.5,2x,i3) 
2014   format("                     ",i4,2x,f8.4,2x,f10.5,2x,i3) 
       write(long_control_file_unit,2013)  cg_mth1(1),scf_mth1(1)+mx_mth1(1),&
                       fermi_de1(1)*27.211396d0,type_fermi1(1)
       do i=2,num_iter1
           write(long_control_file_unit,2014)  cg_mth1(i),scf_mth1(i)+mx_mth1(i), &
                           fermi_de1(i)*27.211396d0,type_fermi1(i)
       enddo
     endif
     !
     !
     CALL read_key_words ( control_file_unit, 'NONLOCAL', LEN('NONLOCAL'), right, readit )
     if(readit) then
        READ ( right, * ) local
     else
        local = 2
     endif
     WRITE (tmp_char, *) local
     write(long_control_file_unit,*) "    NONLOCAL    = "//TRIM(ADJUSTL(tmp_char))
     !
     !
     max_rcut=0.d0
     is_rcut_type_set=0
     ALLOCATE(rcut_of_type(sum_atom_type))
     do i=1,sum_atom_type
         !temp_char = "IN.PSP_RCUT"//char(48+i)
         write(temp_char,*) i
         temp_char="IN.PSP_RCUT"//ADJUSTL(trim(temp_char))
         temp_char = ADJUSTL(temp_char)
         CALL read_key_words ( control_file_unit, TRIM(temp_char),len_trim(temp_char), right, readit)
         if(readit) then
             READ ( right, * ) rcut_of_type(i)
             is_rcut_type_set=1
             write(long_control_file_unit,'(a,f12.5)') '     '//TRIM(temp_char)//" = ",rcut_of_type(i)
         else
             rcut_of_type(i)=3.2d0
         endif
         if(rcut_of_type(i)>max_rcut) then
             max_rcut=rcut_of_type(i)
         endif
     enddo
     !
     !
     CALL read_key_words ( control_file_unit, 'RCUT', LEN('RCUT'), right, readit )
     if(readit) then
        READ ( right, * ) RCUT
     else
        RCUT=3.2
     endif
     if(is_rcut_type_set .eq. 1) then
         rcut=max_rcut
     else
         rcut_of_type(:)=rcut
     endif
     WRITE (tmp_char, "(f12.5)") RCUT
     write(long_control_file_unit,'(a,f12.5)') "     RCUT    = "//TRIM(ADJUSTL(tmp_char))
     do i = 1, sum_atom_type
         !temp_char = "IN.PSP_RCUT"//char(48+i)
         write(temp_char,*) i
         temp_char="IN.PSP_RCUT"//ADJUSTL(trim(temp_char))
         temp_char = ADJUSTL(temp_char)
         write(long_control_file_unit,'(a,f12.5)') '     '//TRIM(temp_char)//" = ",rcut_of_type(i)
     end do
     !
     !
     if(trim(JOB).eq."MD") then
       CALL read_key_words ( control_file_unit, 'MD_DETAIL', LEN('MD_DETAIL'), right, readit )
       if(readit) then
           READ ( right, * ) md_mth, MDstep, dtime, temper_start, temper_end
           md_mth=md_mth*1.001
           MDstep = MDstep*1.001
           write(long_control_file_unit,2006) md_mth,MDstep,dtime,temper_start,temper_end
2006       format("     MD_DETAIL    = ",i4,2x,f9.1,2x,f10.5,2x,f10.5,1x,f10.5)
       else
           message="  MD_DETAIL = xxx, xxx, xxx, xxx, xxx, when JOB = MD"
           call error_message (message)
       endif
     else
       md_mth = 0
       MDstep=0
     endif
     !
     !
     if(md_mth.eq.1.or.md_mth.eq.11) then
       CALL read_key_words ( control_file_unit, 'NSCALE_VVMD', LEN('NSCALE_VVMD'), &
                             right, readit )
       if(readit) then
           READ ( right, * ) step_temp_VVMD
       else
           step_temp_VVMD = 100
       endif
     else
       step_temp_VVMD = 100
     endif
     write (long_control_file_unit, *) '    NSCALE_VVMD = ', step_temp_VVMD
      !       
       CALL read_key_words ( control_file_unit, 'NUM_BLOCKED_PSI', &
     &                       LEN('NUM_BLOCKED_PSI'), right, readit )
       if( readit) then
           READ ( right, * ) inumber_block_psi
           if(inumber_block_psi .lt. 1) inumber_block_psi=1
       else
           inumber_block_psi= 1     ! Gspace implementation, Real-space: 2
       endif
       write(long_control_file_unit,*) "    NUM_BLOCKED_PSI = ", inumber_block_psi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       CALL read_key_words ( control_file_unit, 'WRITE2MEMORY', &
     &                       LEN('WRITE2MEMORY'), right, readit )
       if(readit) then
           READ ( right, * ) iwrite
       else
           iwrite = 1     ! Gspace implementation, Real-space: 2
       endif
       write2memory_7 = iwrite
       write(long_control_file_unit,*) "    WRITE2MEMORY = ", iwrite

      RETURN     
      !
      !*********************************
      !
    end subroutine read_control_file


      !-----------------------------------------------------------------
      !--
      !-----------------------------------------------------------------
      subroutine get_n123(n_t,dd)
      implicit double precision(a-h,o-z)
      real*8 dd
      integer n_t
      integer nnn(2000)   ! FFT dimension table
        nnn(1) =1
        nnn(2) =2
        nnn(3) =3
        nnn(4) =4
        nnn(5) =5
        nnn(6) =6
        nnn(7) =7
        nnn(8) =8
        nnn(9) =9
        nnn(10) =10
        nnn(11) =12
        nnn(12) =14
        nnn(13) =15
        nnn(14) =16
        nnn(15) =18
        nnn(16) =20
        nnn(17) =21
        nnn(18) =24
        nnn(19) =25
        nnn(20) =27
        nnn(21) =28
        nnn(22) =30
        nnn(23) =32
        nnn(24) =35
        nnn(25) =36
        nnn(26) =40
        nnn(27) =42
        nnn(28) =45
        nnn(29) =48
        nnn(30) =49
        nnn(31) =50
        nnn(32) =54
        nnn(33) =56
        nnn(34) =60
        nnn(35) =63
        nnn(36) =64
        nnn(37) =70
        nnn(38) =72
        nnn(39) =75
        nnn(40) =80
        nnn(41) =81
        nnn(42) =84
        nnn(43) =90
        nnn(44) =96
        nnn(45) =98
        nnn(46) =100
        nnn(47) =105
        nnn(48) =108
        nnn(49) =112
        nnn(50) =120
        nnn(51) =125
        nnn(52) =126
        nnn(53) =128
        nnn(54) =135
        nnn(55) =140
        nnn(56) =144
        nnn(57) =147
        nnn(58) =150
        nnn(59) =160
        nnn(60) =162
        nnn(61) =168
        nnn(62) =175
        nnn(63) =180
        nnn(64) =189
        nnn(65) =192
        nnn(66) =196
        nnn(67) =200
        nnn(68) =210
        nnn(69) =216
        nnn(70) =224
        nnn(71) =225
        nnn(72) =240
        nnn(73) =243
        nnn(74) =245
        nnn(75) =250
        nnn(76) =252
        nnn(77) =256
        nnn(78) =270
        nnn(79) =280
        nnn(80) =288
        nnn(81) =294
        nnn(82) =300
        nnn(83) =315
        nnn(84) =320
        nnn(85) =324
        nnn(86) =336
        nnn(87) =343
        nnn(88) =350
        nnn(89) =360
        nnn(90) =375
        nnn(91) =378
        nnn(92) =384
        nnn(93) =392
        nnn(94) =400
        nnn(95) =405
        nnn(96) =420
        nnn(97) =432
        nnn(98) =441
        nnn(99) =448
        nnn(100) =450
        nnn(101) =480
        nnn(102) =486
        nnn(103) =490
        nnn(104) =500
        nnn(105) =504
        nnn(106) =512
        nnn(107) =525
        nnn(108) =540
        nnn(109) =560
        nnn(110) =567
        nnn(111) =576
        nnn(112) =588
        nnn(113) =600
        nnn(114) =625
        nnn(115) =630
        nnn(116) =640
        nnn(117) =648
        nnn(118) =672
        nnn(119) =675
        nnn(120) =686
        nnn(121) =700
        nnn(122) =720
        nnn(123) =729
        nnn(124) =735
        nnn(125) =750
        nnn(126) =756
        nnn(127) =768
        nnn(128) =784
        nnn(129) =800
        nnn(130) =810
        nnn(131) =840
        nnn(132) =864
        nnn(133) =875
        nnn(134) =882
        nnn(135) =896
        nnn(136) =900
        nnn(137) =945
        nnn(138) =960
        nnn(139) =972
        nnn(140) =980
        nnn(141) =1000
        nnn(142) =1008
        nnn(143) =1024
        nnn(144) =1029
        nnn(145) =1050
        nnn(146) =1080
        nnn(147) =1120
        nnn(148) =1125
        nnn(149) =1134
        nnn(150) =1152
        nnn(151) =1176
        nnn(152) =1200
        nnn(153) =1215
        nnn(154) =1225
        nnn(155) =1250
        nnn(156) =1260
        nnn(157) =1280
        nnn(158) =1296
        nnn(159) =1323
        nnn(160) =1344
        nnn(161) =1350
        nnn(162) =1372
        nnn(163) =1400
        nnn(164) =1440
        nnn(165) =1458
        nnn(166) =1470
        nnn(167) =1500
        nnn(168) =1512
        nnn(169) =1536
        nnn(170) =1568
        nnn(171) =1575
        nnn(172) =1600
        nnn(173) =1620
        nnn(174) =1680
        nnn(175) =1701
        nnn(176) =1715
        nnn(177) =1728
        nnn(178) =1750
        nnn(179) =1764
        nnn(180) =1792
        nnn(181) =1800
        nnn(182) =1875
        nnn(183) =1890
        nnn(184) =1920
        nnn(185) =1944
        nnn(186) =1960
        nnn(187) =2000
        nnn(188) =2016
        nnn(189) =2025
        nnn(190) =2048
        nnn(191) =2058
        nnn(192) =2100
        nnn(193) =2160
        nnn(194) =2187
        nnn(195) =2205
        nnn(196) =2240
        nnn(197) =2250
        nnn(198) =2268
        nnn(199) =2304
        nnn(200) =2352
        nnn(201) =2400
        nnn(202) =2401
        nnn(203) =2430
        nnn(204) =2450
        nnn(205) =2500
        nnn(206) =2520
        nnn(207) =2560
        nnn(208) =2592
        nnn(209) =2625
        nnn(210) =2646
        nnn(211) =2688
        nnn(212) =2700
        nnn(213) =2744
        nnn(214) =2800
        nnn(215) =2835
        nnn(216) =2880
        nnn(217) =2916
        nnn(218) =2940
        nnn(219) =3000
        nnn(220) =3024
        nnn(221) =3072
        nnn(222) =3087
        nnn(223) =3125
        nnn(224) =3136
        nnn(225) =3150
        nnn(226) =3200
        nnn(227) =3240
        nnn(228) =3360
        nnn(229) =3375
        nnn(230) =3402
        nnn(231) =3430
        nnn(232) =3456
        nnn(233) =3500
        nnn(234) =3528
        nnn(235) =3584
        nnn(236) =3600
        nnn(237) =3645
        nnn(238) =3675
        nnn(239) =3750
        nnn(240) =3780
        nnn(241) =3840
        nnn(242) =3888
        nnn(243) =3920
        nnn(244) =3969
        nnn(245) =4000
        nnn(246) =4032
        nnn(247) =4050
        nnn(248) =4096
        nnn(249) =4116
        nnn(250) =4200
        nnn(251) =4320
        nnn(252) =4374
        nnn(253) =4375
        nnn(254) =4410
        nnn(255) =4480
        nnn(256) =4500
        nnn(257) =4536
        nnn(258) =4608
        nnn(259) =4704
        nnn(260) =4725
        nnn(261) =4800
        nnn(262) =4802
        nnn(263) =4860
        nnn(264) =4900
        nnn(265) =5000
        nnn(266) =5040
        nnn(267) =5103
        nnn(268) =5120
        nnn(269) =5145
        nnn(270) =5184
        nnn(271) =5250
        nnn(272) =5292
        nnn(273) =5376
        nnn(274) =5400
        nnn(275) =5488
        nnn(276) =5600
        nnn(277) =5625
        nnn(278) =5670
        nnn(279) =5760
        nnn(280) =5832
        nnn(281) =5880
        nnn(282) =6000
        nnn(283) =6048
        nnn(284) =6075
        nnn(285) =6125
        nnn(286) =6144
        nnn(287) =6174
        nnn(288) =6250
        nnn(289) =6272
        nnn(290) =6300
        nnn(291) =6400
        nnn(292) =6480
        nnn(293) =6561
        nnn(294) =6615
        nnn(295) =6720
        nnn(296) =6750
        nnn(297) =6804
        nnn(298) =6860
        nnn(299) =6912
        nnn(300) =7000
        nnn(301) =7056
        nnn(302) =7168
        nnn(303) =7200
        nnn(304) =7203
        nnn(305) =7290
        nnn(306) =7350
        nnn(307) =7500
        nnn(308) =7560
        nnn(309) =7680
        nnn(310) =7776
        nnn(311) =7840
        nnn(312) =7875
        nnn(313) =7938
        nnn(314) =8000
        nnn(315) =8064
        nnn(316) =8100
        nnn(317) =8192
        nnn(318) =8232
        nnn(319) =8400
        nnn(320) =8505
        nnn(321) =8575
        nnn(322) =8640
        nnn(323) =8748
        nnn(324) =8750
        nnn(325) =8820
        nnn(326) =8960
        nnn(327) =9000
        nnn(328) =9072
        nnn(329) =9216
        nnn(330) =9261
        nnn(331) =9375
        nnn(332) =9408
        nnn(333) =9450
        nnn(334) =9600
        nnn(335) =9604
        nnn(336) =9720
        nnn(337) =9800
        nnn(338) =10000
        nnn(339) =10080
        nnn(340) =10125
        nnn(341) =10206
        nnn(342) =10240
        nnn(343) =10290
        nnn(344) =10368
        nnn(345) =10500
        nnn(346) =10584
        nnn(347) =10752
        nnn(348) =10800
        nnn(349) =10935
        nnn(350) =10976
        nnn(351) =11025
        nnn(352) =11200
        nnn(353) =11250
        nnn(354) =11340
        nnn(355) =11520
        nnn(356) =11664
        nnn(357) =11760
        nnn(358) =11907
        nnn(359) =12000
        nnn(360) =12005
        nnn(361) =12096
        nnn(362) =12150
        nnn(363) =12250
        nnn(364) =12288
        nnn(365) =12348
        nnn(366) =12500
        nnn(367) =12544
        nnn(368) =12600
        nnn(369) =12800
        nnn(370) =12960
        nnn(371) =13122
        nnn(372) =13125
        nnn(373) =13230
        nnn(374) =13440
        nnn(375) =13500
        nnn(376) =13608
        nnn(377) =13720
        nnn(378) =13824
        nnn(379) =14000
        nnn(380) =14112
        nnn(381) =14175
        nnn(382) =14336
        nnn(383) =14400
        nnn(384) =14406
        nnn(385) =14580
        nnn(386) =14700
        nnn(387) =15000
        nnn(388) =15120
        nnn(389) =15309
        nnn(390) =15360
        nnn(391) =15435
        nnn(392) =15552
        nnn(393) =15625
        nnn(394) =15680
        nnn(395) =15750
        nnn(396) =15876
        nnn(397) =16000
        nnn(398) =16128
        nnn(399) =16200
        nnn(400) =16384
        nnn(401) =16464
        nnn(402) =16800
        nnn(403) =16807
        nnn(404) =16875
        nnn(405) =17010
        nnn(406) =17150
        nnn(407) =17280
        nnn(408) =17496
        nnn(409) =17500
        nnn(410) =17640
        nnn(411) =17920
        nnn(412) =18000
        nnn(413) =18144
        nnn(414) =18225
        nnn(415) =18375
        nnn(416) =18432
        nnn(417) =18522
        nnn(418) =18750
        nnn(419) =18816
        nnn(420) =18900
        nnn(421) =19200
        nnn(422) =19208
        nnn(423) =19440
        nnn(424) =19600
        nnn(425) =19683
        nnn(426) =19845
        nnn(427) =20000
        nnn(428) =20160
        nnn(429) =20250
        nnn(430) =20412
        nnn(431) =20480
        nnn(432) =20580
        nnn(433) =20736
        nnn(434) =21000
        nnn(435) =21168
        nnn(436) =21504
        nnn(437) =21600
        nnn(438) =21609
        nnn(439) =21870
        nnn(440) =21875
        nnn(441) =21952
        nnn(442) =22050
        nnn(443) =22400
        nnn(444) =22500
        nnn(445) =22680
        nnn(446) =23040
        nnn(447) =23328
        nnn(448) =23520
        nnn(449) =23625
        nnn(450) =23814
        nnn(451) =24000
        nnn(452) =24010
        nnn(453) =24192
        nnn(454) =24300
        nnn(455) =24500
        nnn(456) =24576
        nnn(457) =24696
        nnn(458) =25000
        nnn(459) =25088
        nnn(460) =25200
        nnn(461) =25515
        nnn(462) =25600
        nnn(463) =25725
        nnn(464) =25920
        nnn(465) =26244
        nnn(466) =26250
        nnn(467) =26460
        nnn(468) =26880
        nnn(469) =27000
        nnn(470) =27216
        nnn(471) =27440
        nnn(472) =27648
        nnn(473) =27783
        nnn(474) =28000
        nnn(475) =28125
        nnn(476) =28224
        nnn(477) =28350
        nnn(478) =28672
        nnn(479) =28800
        nnn(480) =28812
        nnn(481) =29160
        nnn(482) =29400
        nnn(483) =30000
        nnn(484) =30240
        nnn(485) =30375
        nnn(486) =30618
        nnn(487) =30625
        nnn(488) =30720
        nnn(489) =30870
        nnn(490) =31104
        nnn(491) =31250
        nnn(492) =31360
        nnn(493) =31500
        nnn(494) =31752
        nnn(495) =32000
        nnn(496) =32256
        nnn(497) =32400
        nnn(498) =32768
        nnn(499) =32805
        nnn(500) =32928
        nnn(501) =33075
        nnn(502) =33600
        nnn(503) =33614
        nnn(504) =33750
        nnn(505) =34020
        nnn(506) =34300
        nnn(507) =34560
        nnn(508) =34992
        nnn(509) =35000
        nnn(510) =35280
        nnn(511) =35721
        nnn(512) =35840
        nnn(513) =36000
        nnn(514) =36015
        nnn(515) =36288
        nnn(516) =36450
        nnn(517) =36750
        nnn(518) =36864
        nnn(519) =37044
        nnn(520) =37500
        nnn(521) =37632
        nnn(522) =37800
        nnn(523) =38400
        nnn(524) =38416
        nnn(525) =38880
        nnn(526) =39200
        nnn(527) =39366
        nnn(528) =39375
        nnn(529) =39690
        nnn(530) =40000
        nnn(531) =40320
        nnn(532) =40500
        nnn(533) =40824
        nnn(534) =40960
        nnn(535) =41160
        nnn(536) =41472
        nnn(537) =42000
        nnn(538) =42336
        nnn(539) =42525
        nnn(540) =42875
        nnn(541) =43008
        nnn(542) =43200
        nnn(543) =43218
        nnn(544) =43740
        nnn(545) =43750
        nnn(546) =43904
        nnn(547) =44100
        nnn(548) =44800
        nnn(549) =45000
        nnn(550) =45360
        nnn(551) =45927
        nnn(552) =46080
        nnn(553) =46305
        nnn(554) =46656
        nnn(555) =46875
        nnn(556) =47040
        nnn(557) =47250
        nnn(558) =47628
        nnn(559) =48000
        nnn(560) =48020
        nnn(561) =48384
        nnn(562) =48600
        nnn(563) =49000
        nnn(564) =49152
        nnn(565) =49392
        nnn(566) =50000
        nnn(567) =50176
        nnn(568) =50400
        nnn(569) =50421
        nnn(570) =50625
        nnn(571) =51030
        nnn(572) =51200
        nnn(573) =51450
        nnn(574) =51840
        nnn(575) =52488
        nnn(576) =52500
        nnn(577) =52920
        nnn(578) =53760
        nnn(579) =54000
        nnn(580) =54432
        nnn(581) =54675
        nnn(582) =54880
        nnn(583) =55125
        nnn(584) =55296
        nnn(585) =55566
        nnn(586) =56000
        nnn(587) =56250
        nnn(588) =56448
        nnn(589) =56700
        nnn(590) =57344
        nnn(591) =57600
        nnn(592) =57624
        nnn(593) =58320
        nnn(594) =58800
        nnn(595) =59049
        nnn(596) =59535
        nnn(597) =60000
        nnn(598) =60025
        nnn(599) =60480
        nnn(600) =60750
        nnn(601) =61236
        nnn(602) =61250
        nnn(603) =61440
        nnn(604) =61740
        nnn(605) =62208
        nnn(606) =62500
        nnn(607) =62720
        nnn(608) =63000
        nnn(609) =63504
        nnn(610) =64000
        nnn(611) =64512
        nnn(612) =64800
        nnn(613) =64827
        nnn(614) =65536
      ntot_tmp=614       !  to be completed later
      do i=2,ntot_tmp
      if(nnn(i-1).le.dd.and.dd.lt.nnn(i)) then
      n_t=nnn(i)

      goto 122
      endif
      enddo
      write(6,*) "dd out of range for n",dd,nnn(ntot)
      stop
122   continue
      return
      end subroutine get_n123
 

      subroutine get_n456(n_t)
      implicit double precision(a-h,o-z)
      integer n_t
      integer nnn(2000)   ! FFT dimension table
        nnn(1) =1
        nnn(2) =2
        nnn(3) =3
        nnn(4) =4
        nnn(5) =5
        nnn(6) =6
        nnn(7) =7
        nnn(8) =8
        nnn(9) =9
        nnn(10) =10
        nnn(11) =12
        nnn(12) =14
        nnn(13) =15
        nnn(14) =16
        nnn(15) =18
        nnn(16) =20
        nnn(17) =21
        nnn(18) =24
        nnn(19) =25
        nnn(20) =27
        nnn(21) =28
        nnn(22) =30
        nnn(23) =32
        nnn(24) =35
        nnn(25) =36
        nnn(26) =40
        nnn(27) =42
        nnn(28) =45
        nnn(29) =48
        nnn(30) =49
        nnn(31) =50
        nnn(32) =54
        nnn(33) =56
        nnn(34) =60
        nnn(35) =63
        nnn(36) =64
        nnn(37) =70
        nnn(38) =72
        nnn(39) =75
        nnn(40) =80
        nnn(41) =81
        nnn(42) =84
        nnn(43) =90
        nnn(44) =96
        nnn(45) =98
        nnn(46) =100
      ntot_tmp=46 !  to be completed later

      do i=2,ntot_tmp
        if(nnn(i-1).le.n_t.and.n_t.lt.nnn(i)) then
          n_t=nnn(i)
          goto 123
        endif
      enddo
123   continue
      return
      
      end subroutine get_n456


