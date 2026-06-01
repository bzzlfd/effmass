   SUBROUTINE cal_ityp ()
   !
   USE atom_base, ONLY : atom, sum_atom, ityp
   !
   IMPLICIT NONE
   !
   INTEGER :: tmp_int, i, j
   character(len=200) :: tmp_char
   !
   ALLOCATE(ityp(sum_atom))
   tmp_int = atom(1)%number
   j = 1
   DO i = 1, sum_atom
      IF (tmp_int == atom(i)%number) THEN
         ityp(i) = j
      ELSE
         j = j + 1
         ityp(i) = j
      ENDIF 
      tmp_int = atom(i)%number
   ENDDO  
   !
   do i = 1, sum_atom
       write (tmp_char, *) ityp(i)
   enddo
   !
   END SUBROUTINE cal_ityp
