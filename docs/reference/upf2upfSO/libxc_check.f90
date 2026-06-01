subroutine libxc_check(id1, id2, spin)
    use xc_f90_lib_m
    use modlibxc
    use xc, only : ismgga
    
    implicit none
    
    integer :: id1, id2
    integer :: major, minor
    integer :: spin
    integer :: i
    type(xc_f90_pointer_t) :: str
    integer :: name

    funcs(1)%id = id1
    funcs(2)%id = id2
    do i = 1, 2
	if ( funcs(i)%id < 1 ) then
	    write ( 6, * ) "Not right libxc functional"
	    stop
	end if
    end do

    ! version check
    call xc_f90_version(major, minor)
    if ( major < 2 ) then
	write ( 6, * ) "the libxc version should > 2.1.x"
	stop
    else
	if ( minor < 1 ) then
	    write ( 6, * ) "the libxc version should > 2.1.x"
	    stop
	end if
    end if

    ! initialize all the parameter
    ismgga = .false.
    do i = 1, 2
	funcs(i)%nspin = spin
	funcs(i)%family = xc_f90_family_from_id( funcs(i)%id )
	if ( funcs(i)%family == XC_FAMILY_MGGA ) ismgga = .true.
	call xc_f90_functional_get_name( funcs(i)%id, funcs(i)%name )
	call xc_f90_func_init &
	    ( funcs(i)%p, funcs(i)%info, funcs(i)%id, funcs(i)%nspin )
	!
	if ( xc_f90_info_kind(funcs(i)%info) == XC_EXCHANGE ) funcs(i)%kind = 0
	if ( xc_f90_info_kind(funcs(i)%info) == XC_CORRELATION) funcs(i)%kind = 1
	if ( xc_f90_info_kind(funcs(i)%info) == XC_EXCHANGE_CORRELATION ) funcs(i)%kind = 2
	if ( xc_f90_info_kind(funcs(i)%info) == XC_KINETIC ) funcs(i)%kind = 3
	!
	name = 0 ! number of the reference, must be 0 in the first call
	call xc_f90_info_refs( funcs(i)%info, name, str, funcs(i)%ref )
	!
    end do
    !
    if ( .false. ) then
    write ( 6, * ) "ismgga:", ismgga
    write ( 6, * ) "*******************************************************"
	do i = 1, 2
	    write ( 6, * ) "id:", funcs(i)%id
	    write ( 6, * ) "family:", funcs(i)%family
	    write ( 6, * ) "spin:", funcs(i)%nspin
	    write ( 6, * ) "kind:", funcs(i)%kind
	    write ( 6, * ) "name:", trim(funcs(i)%name)
	    write ( 6, * ) "ref:", trim(funcs(i)%ref)
	end do
    write ( 6, * ) "*******************************************************"
    end if
    ! ensure contain LDA/GGA/MGGA
    do i = 1, 2
	if ( .not. ( &
	    funcs(i)%family == XC_FAMILY_LDA .or. &
	    funcs(i)%family == XC_FAMILY_GGA .or. &
	    funcs(i)%family == XC_FAMILY_MGGA ) ) then
	    write ( 6, * ) 'Now PWmat only support libxc with FAMILY: "LDA/GGA/MGGA"'
	    stop
	end if
    end do
    ! ensure contain exchange and correlation energy 
    if ( .not. ( &
	( funcs(1)%kind == 0 .and. funcs(2)%kind == 1 ) .or. &
	( funcs(2)%kind == 0 .and. funcs(1)%kind == 1 ) .or. &
	( funcs(1)%kind == 2 .and. funcs(2)%kind < 3 ) .or. &
	( funcs(2)%kind == 2 .and. funcs(1)%kind < 3 )  ) ) then
	write ( 6, * ) 'libxc should be set as following:'
	write ( 6, * ) '1: "_C_ + _X_"'
	write ( 6, * ) '2: "_X_ + _C_"'
	write ( 6, * ) '3: "_XC_"'
	stop
    end if
    ! ensure contain energy and potential
    if ( .not. &
        ( &
	( iand(xc_f90_info_flags(funcs(1)%info),XC_FLAGS_HAVE_EXC) /= 0 .and. &
	  iand(xc_f90_info_flags(funcs(2)%info),XC_FLAGS_HAVE_VXC) /= 0 ) .or. &
	( iand(xc_f90_info_flags(funcs(2)%info),XC_FLAGS_HAVE_EXC) /= 0 .and. &
	  iand(xc_f90_info_flags(funcs(1)%info),XC_FLAGS_HAVE_VXC) /= 0 ) .or. &
	( iand(xc_f90_info_flags(funcs(1)%info),XC_FLAGS_HAVE_EXC) /= 0 .and. &
	  iand(xc_f90_info_flags(funcs(1)%info),XC_FLAGS_HAVE_VXC) /= 0 ) .or. &
	( iand(xc_f90_info_flags(funcs(2)%info),XC_FLAGS_HAVE_EXC) /= 0 .and. &
	  iand(xc_f90_info_flags(funcs(2)%info),XC_FLAGS_HAVE_VXC) /= 0 ) &
        ) &
       ) then
	write ( 6, * ) 'unsupport libxc, please change xc functional'
	stop
    end if
    !
    do i = 1, 2
	call xc_f90_func_end ( funcs(i)%p )
    end do
    !
end subroutine libxc_check

