   MODULE constant
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
   REAL(KIND=8), PARAMETER :: bohr_radius = 0.52917720859E-10  ! m
   INTEGER, PARAMETER :: npk = 40000  ! max number of kpoints
   real(kind=8), parameter :: pi = 4*datan(1.0_dp)
   REAL(DP), PARAMETER :: tpi    = 2.0_DP * pi
   REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
   !
   END MODULE constant
