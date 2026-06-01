    subroutine read_atom_file ()
    !
    use control_config, only : atom_file
    use file_unit, only : atom_file_unit
    use atom_base
    use constant
    !
    implicit none
    !
    logical :: alive, scanit
    !
    integer :: i, j, tmp_int, stat
    !
    character(len=200) :: tmp_char
    character(len=200) :: tmp_char_arr(10000)
    !
    integer, allocatable :: tmp_int_(:)
    !
    rewind (atom_file_unit)
    read (atom_file_unit, *, iostat = stat) sum_atom
    if (stat /= 0) then
        call error_message &
	("not read the total number of atoms correctly in file: "//trim(atom_file))
    else
        allocate (atom(sum_atom))
    endif
    !
    call scan_key_words (atom_file_unit, "LATTICE", len("LATTICE"), scanit)
    if (scanit) then
        do i = 1, 3
            read (atom_file_unit, *, iostat = stat) latt_vec_ang(1:3,i)
	        if (stat /= 0) then
	            call error_message &
		("Not read the lattice vectors correctly in file: "//trim(atom_file))
	        endif
	    enddo 
    else
        call error_message &
	("Not find the key words: ""LATTICE"" in file "//""""//trim(atom_file)//"""")
    endif
    do i = 1, 3
        write (tmp_char_arr(i), "(3f15.10)") latt_vec_ang(1:3,i)
    enddo    
    latt_vec_bohr = latt_vec_ang * 1.0E-10 / bohr_radius
    !
    call scan_key_words (atom_file_unit, "POSITION", len("POSITION"), scanit)
    if (scanit) then
      DO i = 1, sum_atom
         READ (atom_file_unit, *, IOSTAT = stat) &
               atom(i)%number, atom(i)%position(1:3), atom(i)%mov(1:3)
         IF (stat == 0) THEN
            CALL gen_element_name (atom(i)%number, atom(i)%label)
	    write (tmp_char_arr(i), "(A2,I3,3f10.5,3I3)") atom(i)%label, atom(i)%number, &
	                                      atom(i)%position(1:3), atom(i)%mov(1:3)
         ELSEIF (stat > 0) THEN
            call error_message &
	    ("Not read the atomic positions or if_mov correctly in file "//""""//trim(atom_file)//"""")
         ELSE
            IF (i <= sum_atom) THEN
	      write (tmp_char, *) i-1
              call error_message &
	      ('Miss atom positions, only find '//trim(adjustl(tmp_char))//' atoms')
            ENDIF
         ENDIF
	 !
	 atom(i)%pos_ang(1) = atom(i)%position(1)*latt_vec_ang(1,1) + &
	                      atom(i)%position(2)*latt_vec_ang(1,2) + &
			      atom(i)%position(3)*latt_vec_ang(1,3)
	 atom(i)%pos_ang(2) = atom(i)%position(1)*latt_vec_ang(2,1) + &
	                      atom(i)%position(2)*latt_vec_ang(2,2) + &
			      atom(i)%position(3)*latt_vec_ang(2,3)
	 atom(i)%pos_ang(3) = atom(i)%position(1)*latt_vec_ang(3,1) + &
	                      atom(i)%position(2)*latt_vec_ang(3,2) + &
			      atom(i)%position(3)*latt_vec_ang(3,3)
	 atom(i)%pos_bohr(1:3) = atom(i)%pos_ang(1:3) * 1.0E-10 / bohr_radius
	 !
	 !write ( 6, * ) atom(i)%pos_ang(1:3)
         !
      ENDDO
      !
    else
        call error_message &
	("Not find the key words: ""POSITION"" in file "//""""//trim(atom_file)//"""")
    endif  
    !
      !
      !  calculate the total number of atom type
      !
      allocate ( tmp_int_(sum_atom) )
      do i = 1, sum_atom
	  tmp_int_(i) = atom(i)%number
      end do
      !
      sum_atom_type = 1
      iloop:do i = 1, sum_atom
            if ( tmp_int_(i) == 0 ) then
		cycle iloop
	    else
		tmp_int = tmp_int_(i)
	        jloop:do j = i, sum_atom
		    if ( tmp_int == atom(j)%number ) tmp_int_(j) = 0
		end do jloop
		if ( sum(tmp_int_) /= 0 ) sum_atom_type = sum_atom_type + 1
	    end if
      end do iloop
      !write ( 6, * ) "sum_atom_type", sum_atom_type
      !stop
    end subroutine read_atom_file
