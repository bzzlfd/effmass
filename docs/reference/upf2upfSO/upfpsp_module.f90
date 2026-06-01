!=====================================================!
!=====================================================!

       module double_precision
        implicit none
        integer, parameter :: dp = selected_real_kind(14,200)
       end module double_precision
!=====================================================!
!=====================================================!

       module element_type
        use double_precision
        implicit none

        type element
            CHARACTER(LEN=2) :: label      ! Element label
            REAL(DP) :: ecutwfc(2)           ! suggested cut-off for wfc
            real(dp) :: mass                           ! suggested mass
            integer :: atnum               ! atomic number
        end type
    
       end module element_type
!=====================================================!
!=====================================================!
       module default_element
        use element_type
        implicit none
        public
        type (element) :: elements(108) ! total 108 species atom
        contains
        subroutine initial_elements ()
            use double_precision
            implicit none
            integer :: i
            real(dp) :: mass(108), ecutwfc(108,2)
            character(len=2) :: label(108)
            ! mass from NIST SP 966 ( March 2013 )
            data mass &
     &      / 1.008, 4.002602, 6.94, 9.012182, 10.81, 12.011, 14.007, &
     &        15.999, 18.9984032, 20.197, &
     &        22.98976928, 24.3050, 26.9815386, 28.085, 30.973762, &
     &        32.06, 35.45, 39.948, 39.0983, 40.078, &
     &        44.955912, 47.867, 50.9415, 51.9961, 54.938045,&
     &        55.845, 58.933195, 58.6934, 63.546, 65.38, &
     &        69.723, 72.63, 74.92160, 78.96, 79.904, 83.798, 85.4678,&
     &        87.62, 88.90585, 91.224, &
     &        92.90638, 95.96, 98, 101.07, 102.9055, 106.42, 107.8682,&
     &       112.411, 114.818, 118.710, &
     &       121.760, 127.60, 126.90447, 131.293, 132.9054519, 137.327,&
     &       138.90547, 140.116, 140.90765, 144.242, &
     &       145, 150.36, 151.964, 157.25, 158.92535, 162.5, 164.93032,&
     &       167.259, 168.93421, 173.054, &
     &       174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23,&
     &       192.217, 195.084, 196.966569, 200.59, &
     &       204.38, 207.2, 208.98040, 209, 210, 222, 223, 226,&
     &       227, 232.03806, &
     &       231.03588, 238.02891, 237, 244, 243, 247, 247, 251,&
     &       252, 257, &
     &       258, 259, 262, 265, 268, 271, 270, 277 /
        ! element name
            data label  &
     &    / "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", &
     &    "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", &
     &    "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
     &    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", &
     &    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
     &    "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
     &    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",&
     &    "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",&
     &    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",&
     &    "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",&
     &    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs" /
            ! ecutwfc(1) NCPP | ecutwfc(2) USPP
            ! the value "10.00_dp" not test, please note! 
            data ecutwfc(1:108,1) &
            &  / 50., 50., 50., 50., 50., 50., 50., 50., 50., 50.,&
     &  50., 50., 50., 50., 50., 50., 50., 50., 50., 50., &
     &  50., 50., 50., 50., 50., 50., 50., 50., 50., 50., &
     &  50., 50., 50., 50., 50., 50., 50., 50., 50., 50., &
     &  50., 50., 50., 50., 50., 50., 50., 50., 50., 50., &
     &  50., 50., 50., 50., 50., 50., 50.00, 50.00, 50.00, 50.00,&
     &  50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00,&
     &  50.00,&
     &  50.00, 50., 50., 50., 50., 50., 50., 50., 50., 50.,&
     &  50., 50., 50., 50., 50., 50., 50., 50., 50., 50.00,&
     &  50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00,&
     &  50.00,&
     &  50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00, 50.00 / ! Ry
            data ecutwfc(1:108,2) &
     & / 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00,&
     &   20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00,&
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00,&
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00,&
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, &
     & 20.00,&
     & 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00, 20.00 / ! Ry
            do i = 1, 108
                elements(i)%atnum = i
                elements(i)%mass = mass(i)
                elements(i)%label = label(i)
                elements(i)%ecutwfc(1) = ecutwfc(i,1)
                elements(i)%ecutwfc(2) = ecutwfc(i,2)
            end do
         end subroutine initial_elements

        end module default_element

!=====================================================!
!=====================================================!

        module upf_module
          use double_precision   
          implicit none
          type atompp
          CHARACTER(LEN=80):: generated=' '! generator software
          CHARACTER(LEN=80):: author=' '   ! pseudopotential's author
          CHARACTER(LEN=80):: date=' '     ! generation date
          CHARACTER(LEN=80):: comment=' '  ! author's comment
          CHARACTER(LEN=2) :: psd=' '      ! Element label
          CHARACTER(LEN=20) :: typ=' '     ! Pseudo type ( NC or US or PAW)
          CHARACTER(len=6) :: rel=' '      ! relativistic: {no|scalar|full}
          LOGICAL :: tvanp              ! .true. if Ultrasoft
          LOGICAL :: tcoulombp          ! .true. if Coulomb 1/r potential
          LOGICAL :: nlcc               ! Non linear core corrections
          CHARACTER(LEN=25) :: dft      ! Exch-Corr type
          REAL(DP) :: zp                ! z valence
          REAL(DP) :: etotps            ! total energy
          REAL(DP) :: ecutwfc           ! suggested cut-off for wfc
          REAL(DP) :: ecutrho           ! suggested cut-off for rho
          real(dp) :: mass                           ! suggested mass
          integer :: num
          !
          CHARACTER(len=11) :: nv       ! UPF file three-digit version i.e. 2.0.0
          INTEGER :: lmax               ! maximum l component in beta
          INTEGER :: lmax_rho           ! max l component in charge (should be 2*lmax)
          REAL(DP), allocatable :: vnl(:,:,:) ! vnl(i,l,s) = V(r_i)_{ls}
          ! only for single-channel NC PP
          ! Wavefunctions and projectors
          INTEGER :: nwfc               ! number of atomic wavefunctions
          INTEGER :: nbeta              ! number of projectors
          INTEGER,  allocatable :: kbeta(:) ! kbeta(nbeta) see below
          INTEGER :: kkbeta             ! kkbeta=max(kbeta(:))
          !  kbeta<=mesh is the number of grid points for each beta function
          !              beta(r,nb) = 0 for r > r(kbeta(nb))
          ! kkbeta<=mesh is the largest of such number so that for all beta
          !              beta(r,nb) = 0 for r > r(kkbeta)
          !
          INTEGER,  allocatable :: lll(:)     ! lll(nbeta) l of each projector
          REAL(DP), allocatable :: beta(:,:)  ! beta(mesh,nbeta) projectors
          !
          CHARACTER(LEN=2), allocatable :: els(:)  ! els(nwfc) label of wfc
          CHARACTER(LEN=2), allocatable :: els_beta(:)  ! els(nbeta) label of beta
          INTEGER, allocatable  :: nchi(:)    ! lchi(nwfc) value of pseudo-n for wavefcts
          INTEGER, allocatable  :: lchi(:)    ! lchi(nwfc) value of l for wavefcts
          REAL(DP), allocatable :: oc(:)      ! oc(nwfc) occupancies for wavefcts
          REAL(DP), allocatable :: epseu(:)   ! pseudo one-particle energy (nwfc)
          REAL(DP), allocatable :: rcut_chi(:)! rcut_chi(nwfc) cutoff inner radius
          REAL(DP), allocatable :: rcutus_chi(:)! rcutus_chi(nwfc) ultrasoft outer radius
          ! Chi and rho_at are only used for initial density and initial wfcs:
          REAL(DP), allocatable :: chi(:,:)   ! chi(mesh,nwfc) atomic wavefcts
          REAL(DP), allocatable :: rho_at(:)  ! rho_at(mesh) atomic charge
          ! Minimal radial grid:
          INTEGER :: mesh               ! number of points in the radial mesh
          REAL(DP) :: xmin              ! the minimum x of the linear mesh
          REAL(DP) :: rmax              ! the maximum radius of the mesh
          REAL(DP) :: zmesh             ! the nuclear charge used for mesh
          REAL(DP) :: dx                ! the deltax of the linear mesh
          REAL(DP), allocatable :: r(:)     ! r(mesh)  radial grid
          REAL(DP), allocatable :: rab(:)   ! rab(mesh) dr(x)/dx (x=linear grid)
          ! Pseudized core charge
          REAL(DP), allocatable :: rho_atc(:) ! rho_atc(mesh) atomic core charge
          ! Local potential
          INTEGER :: lloc                 ! L of channel used to generate local potential
          ! (if < 0 it was generated by smoothing AE potential)
          REAL(DP) :: rcloc               ! vloc = v_ae for r > rcloc
          REAL(DP), allocatable :: vloc(:)    ! vloc(mesh) local atomic potential
          !
          REAL(DP), allocatable :: dion(:,:)  ! dion(nbeta,nbeta) atomic D_{mu,nu}
          ! Augmentation
          LOGICAL :: q_with_l              ! if .true. qfunc is pseudized in
          ! different ways for different l
          INTEGER :: nqf                  ! number of Q coefficients
          INTEGER :: nqlc                 ! number of angular momenta in Q
          REAL(DP):: qqq_eps              ! qfunc is null if its norm is .lt. qqq_eps
          REAL(DP), allocatable :: rinner(:)  ! rinner(0:2*lmax) r_L
          REAL(DP), allocatable :: qqq(:,:)   ! qqq(nbeta,nbeta) q_{mu,nu}
          ! Augmentation without L dependecy
          REAL(DP), allocatable :: qfunc(:,:) ! qfunc(mesh,nbeta*(nbeta+1)/2)
          ! Q_{mu,nu}(|r|) function for |r|> r_L
          ! Augmentation depending on L (optional, compulsory for PAW)
          REAL(DP), allocatable :: qfuncl(:,:,:)!  qfuncl(mesh,nbeta*(nbeta+1)/2,l)
          ! Q_{mu,nu}(|r|) function for |r|> r_L
          ! Analitycal coeffs cor small r expansion of qfunc (Vanderbilt's code)
          REAL(DP), allocatable :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          ! coefficients for Q for |r|<r_L
          ! All electron and pseudo wavefunction, pswfc differ from chi as they are
          ! one for each beta, not just some choosen for initial conditions
          LOGICAL           :: has_wfc    ! if true, UPF contain AE and PS wfc for each beta
          REAL(DP), allocatable :: aewfc(:,:) ! wfc(mesh,nbeta) all-electron wfc
          REAL(DP), allocatable :: pswfc(:,:) ! wfc(mesh,nbeta) pseudo wfc

          LOGICAL :: has_so             ! if .true. includes spin-orbit
          INTEGER, allocatable :: nn(:)     ! nn(nwfc) quantum number of wfc
          REAL(DP), allocatable :: rcut(:)  ! cut-off radius(nbeta)
          REAL(DP), allocatable :: rcutus(:)! ultrasoft cut-off radius (nbeta)
          REAL(DP), allocatable :: jchi(:)  ! jchi(nwfc) j=l+1/2 or l-1/2 of wfc
          REAL(DP), allocatable :: jjj(:)   ! jjj(nbeta) j=l+1/2 or l-1/2 of beta
          end type atompp 
       end module upf_module

!=====================================================!
!=====================================================!
       module pseudo_potential_upf
        use upf_module
        implicit none
        type(atompp), allocatable :: upfpsp(:)

       end module pseudo_potential_upf

!=====================================================!
!
