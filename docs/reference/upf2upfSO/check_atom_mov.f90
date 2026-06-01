      SUBROUTINE check_atom_mov ()
      !
      USE control_config, ONLY : job, atom_file
      USE atom_base, ONLY : atom, sum_atom
      !
      IMPLICIT NONE
      !
      INTEGER :: i, j
      LOGICAL :: tmp_logic = .FALSE.
      !
      IF (TRIM(job) == "MD") THEN
         DO i = 1, sum_atom
            DO j = 1, 3
               IF (atom(i)%mov(j) == 1) tmp_logic = .TRUE.
            ENDDO
         ENDDO
         IF (.NOT. tmp_logic) THEN
	     call error_message &
	     ("For JOB = "//trim(job)//&
	     " some of the mov x,y,z should be 1 in "//trim(atom_file))
         ENDIF        
      ENDIF
      END SUBROUTINE check_atom_mov
 
