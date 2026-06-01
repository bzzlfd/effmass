   !--------------------------------------------------------------------
   !--
   !--------------------------------------------------------------------
      SUBROUTINE read_vwr_head ( pp_file_name, ecut, amass, iiatom_, zatom )
      !
      IMPLICIT NONE
      !
      LOGICAL :: alive
      CHARACTER(LEN=*), INTENT(IN) :: pp_file_name
      REAL(KIND=8), INTENT(OUT) :: ecut
      CHARACTER(LEN=200) :: file_name
      INTEGER :: vwr_file_io, i
      INTEGER :: nrr, ic, iloc, iso, ierror
      REAL(KIND=8) :: zatom, occ_s, occ_p, occ_d, amass, iiatom
      INTEGER :: iiatom_
      !
      vwr_file_io = 40
      !
      file_name = pp_file_name
      file_name = ADJUSTL(file_name)
      !
      INQUIRE ( FILE = TRIM(file_name), EXIST = alive )
      IF ( .NOT. alive ) THEN
          call error_message (""""//trim(file_name)//""" doesn't exist!")
      ELSE 
          OPEN ( UNIT = vwr_file_io, FILE = TRIM(file_name), IOSTAT = ierror )
          REWIND ( vwr_file_io )
          READ ( vwr_file_io, * ) nrr, ic, iiatom, zatom, iloc, occ_s, &
                                 occ_p, occ_d, iso, ecut, amass
      ENDIF
      !
      CLOSE ( vwr_file_io )
      iiatom_ = ANINT(iiatom)
      !
      !
      END SUBROUTINE read_vwr_head
   !--------------------------------------------------------------------
   !--
   !--------------------------------------------------------------------

