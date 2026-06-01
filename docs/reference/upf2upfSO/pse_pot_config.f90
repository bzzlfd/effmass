   MODULE pse_pot_config
   !
   IMPLICIT NONE
   PUBLIC
   SAVE
   !
   TYPE pseudo_potential
      REAL(KIND=8) :: ecut, & ! cutoff energy of wavefunction
                      mass, & ! mass of the atom
                      zval    ! valence of the atom
      INTEGER :: number  ! atomic number, it is real type, not integer!!!
      CHARACTER(LEN=200) :: file_name ! pseudo potential file name
      CHARACTER(LEN=2) :: element
   END TYPE
   !
   TYPE(pseudo_potential), ALLOCATABLE :: pp(:)
   TYPE(pseudo_potential), ALLOCATABLE :: atom_pp(:)
   !
   END MODULE pse_pot_config
