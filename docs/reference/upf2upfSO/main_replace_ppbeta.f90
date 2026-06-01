program main_replace_ppbeta
    use upf_module
    use pseudo_potential_upf ! upfpsp(:)


    character(len=1024) :: upfname
    character(len=1024) :: upfname_src
    character(len=2) :: elename
    integer :: ii
    integer :: min_mesh

    !
    write(*,*) "input filename of pseudopotential which PP_BETA will be replaced:"
    read(5,*) upfname
    write(*,*) "input filename of pseudopotential which PP_BETA will be sourced:"
    read(5,*) upfname_src
    allocate(upfpsp(2))
    call read_upf_v2(upfname, upfpsp(1),0)
    call read_upf_v2(upfname_src, upfpsp(2),0)

    if(upfpsp(1)%lmax .ne. upfpsp(2)%lmax) then
        write(*,*) "l_max not same in pseudopotential files"
        stop
    endif
    !if(upfpsp(1).mesh .gt. upfpsp(2).mesh) then
    !    write(*,*) "mesh_size not same in pseudopotential files "
    !    stop
    !endif
    min_mesh=minval(upfpsp(1:2)%mesh)
    do k=1,upfpsp(2)%nbeta
        upfpsp(1)%beta(:,k)=0.d0
        do ii=1,min_mesh
            upfpsp(1)%beta(ii,k) = upfpsp(2)%beta(ii,k)
        enddo
        upfpsp(1)%dion(k,k) = upfpsp(2)%dion(k,k)
    enddo

    write(*,*)
    write(*,*) "write new psps in:"
    upfname=adjustl(trim(upfname))//".REP"
    write(*,*) adjustl(trim(upfname))
    open(100,file=upfname)
    call write_upf_v2(100,upfpsp(1))
    close(100)

    deallocate(upfpsp)

end program main_replace_ppbeta
