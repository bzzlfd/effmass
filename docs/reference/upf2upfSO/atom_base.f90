   MODULE atom_base
   !
   IMPLICIT NONE
   PUBLIC
   SAVE
   !
   INTEGER :: sum_atom_type, &  ! the total number of atom type
              sum_atom     ! the total number of atom
   REAL(KIND=8) :: latt_vec_ang(3,3)  ! lattice vector in the unit of angstrom (10^-10m)
   REAL(KIND=8) :: latt_vec_bohr(3,3)  ! lattice vector in the unit of bohr (0.52917720859E-10m)
   REAL(KIND=8) :: latt_vec(3,3)  ! lattice vector, unit alat
   REAL(KIND=8) :: alat    ! the parameter to uniform the lattice vector, bohr
   REAL(KIND=8) :: rec_vec(3,3)  ! the reciprocal lattice vector, unit 2pi / alat
   REAL(KIND=8) :: omega   ! the volume of the lattice unit: alat^3
   REAL(KIND=8) :: omega_bohr   ! the volume of the lattice unit: bohr^3
   INTEGER, ALLOCATABLE :: ityp(:)
   !
   TYPE :: atom_info
     INTEGER :: number ! the atomic number
     CHARACTER(LEN=2) :: label  ! the atomic name
     REAL(KIND=8) :: position(3) ! atomic position in the relative coordinates of the lattice vector
     real(kind=8) :: pos_bohr(3) ! atomic position in bohr
     real(kind=8) :: pos_ang(3)  ! atomic position in angstrom
     INTEGER :: mov(3) !  the componets i of force for the atom is multipled by if_mov(i, num_atom) which is 0 or 1. These parameters are used in MD or Optimization.
   END TYPE
   !
   TYPE(atom_info), ALLOCATABLE :: atom(:)
   !
   END MODULE atom_base
