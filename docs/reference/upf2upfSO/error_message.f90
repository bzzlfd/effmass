    subroutine error_message(info_error)
    !
    use file_unit, only : io_stdout
    implicit none
    !
    character(len=*) :: info_error
    !
    
    write (io_stdout, *)
    write (io_stdout, "(A)") "**************************************************"
    write (io_stdout, "(4x,A)") "ERROR !"
    write (io_stdout, "(4x,A)") trim(info_error)
    write (io_stdout, "(4x,A)") "STOP !"
    write (io_stdout, "(A)") "**************************************************"
    write (io_stdout, *)
    stop
    !
    end subroutine error_message

