      !-----------------------------------------------------------------
      !--
      !-----------------------------------------------------------------
      SUBROUTINE gen_element_name( element_number, element_name )
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: element_number
      INTEGER :: i
      CHARACTER(len=2), INTENT(OUT) :: element_name
      CHARACTER(len=2) :: element_names(112)
      !
      DATA  element_names &
        /                                                             &
          "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
          "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
          "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
          "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
          "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
          "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
          "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
          "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
          "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
          "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
          "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", &
          "Rg", "Ch" &
        /
      !
      IF( element_number > 112 .OR. element_number < 1 ) &
          call error_message ("atomic number not between 1 and 112")
      DO i = 1, 112
         IF( i == element_number ) THEN
            element_name = element_names(i)
         ENDIF
      ENDDO
      !
      END SUBROUTINE gen_element_name

      !-----------------------------------------------------------------
      !--
      !-----------------------------------------------------------------
      SUBROUTINE gen_element_num ( element_name, element_number )
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(OUT) :: element_number
      INTEGER :: i
      CHARACTER(len=2), INTENT(IN) :: element_name
      CHARACTER(len=2) :: element_names(112)
      logical :: dropin = .false.
      !
      DATA  element_names &
        /                                                             &
          "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
          "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
          "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
          "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
          "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
          "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
          "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
          "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
          "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
          "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
          "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", &
          "Rg", "Ch" &
        /
      !
      DO i = 1, 112
         IF( element_name == element_names(i)) THEN
           element_number = i
	   dropin = .true.
	   exit
         ENDIF
      ENDDO
      if (.not. dropin) then
          call error_message &
      ("The atomic label should like this: H, He, Li, Be, ... ... ")
      endif
      !
      END SUBROUTINE gen_element_num
      !-----------------------------------------------------------------
      !--
      !-----------------------------------------------------------------

