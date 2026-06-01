    module file_unit
    !
    implicit none
    !
    integer :: control_file_unit = 100, &
               long_control_file_unit = 101, &
	       atom_file_unit = 102, &
           kpoint_file_unit = 103, &
           symm_file_unit = 104, &
           pwscf_file_unit = 200, &
	   in_pwscf_unit = 201
    integer :: io_stdout = 6

    end module file_unit
