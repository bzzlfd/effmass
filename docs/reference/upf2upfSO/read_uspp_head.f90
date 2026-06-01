   !--------------------------------------------------------------------
   !--
   !--------------------------------------------------------------------
      SUBROUTINE read_uspp_head (pp_file_name, ecut, amass, zps_, zvps)
      !
      IMPLICIT NONE
      !
      LOGICAL :: alive
      CHARACTER(LEN=*), INTENT(IN) :: pp_file_name
      REAL(KIND=8), INTENT(OUT) :: ecut
      CHARACTER(LEN=200) :: file_name
      INTEGER :: us_file_io, i
      INTEGER :: iver(3), idmy(3)
      CHARACTER(LEN=20) :: title
      REAL(KIND=8) :: zps, zvps, exftps, etot
      INTEGER :: nvalps, mesh, amass, ierror
      INTEGER :: zps_
      !
      !
      us_file_io = 30
      file_name = pp_file_name
      file_name = ADJUSTL(file_name)
      !
      !
      INQUIRE ( FILE = TRIM(file_name), EXIST = alive )
      IF ( .NOT. alive ) THEN
          call error_message (""""//trim(file_name)//""" doesn't exist!")
      ELSE 
          OPEN ( UNIT = us_file_io, FILE = TRIM(file_name), IOSTAT = ierror )
          REWIND ( us_file_io )
          READ ( us_file_io, * ) (iver(i),i=1,3),(idmy(i),i=1,3)
          !
          READ ( us_file_io, * ) title, zps, zvps, exftps, nvalps,&
                                 mesh, etot, ecut, amass
          !
          CLOSE ( us_file_io )
      ENDIF
      zps_ = ANINT(zps)
      !
      !
      END SUBROUTINE read_uspp_head

