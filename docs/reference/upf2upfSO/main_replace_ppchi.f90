program main_replace_ppchi
    use upf_module
    use pseudo_potential_upf ! upfpsp(:)


    character(len=1024) :: upfname
    character(len=1024) :: upfname_src
    character(len=2) :: elename
    integer :: ii
    integer :: min_mesh

    !
    write(*,*) "input filename of pseudopotential which PP_CHI will be replaced:"
    read(5,*) upfname
    write(*,*) "input filename of pseudopotential which PP_CHI will be sourced:"
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
    upfpsp(1)%nwfc=upfpsp(2)%nwfc
    !
    deallocate(upfpsp(1)%chi)
    deallocate(upfpsp(1)%els)
    deallocate(upfpsp(1)%oc)
    deallocate(upfpsp(1)%lchi)
    deallocate(upfpsp(1)%nchi)
    deallocate(upfpsp(1)%rcut_chi)
    deallocate(upfpsp(1)%rcutus_chi)
    deallocate(upfpsp(1)%epseu)

    allocate(upfpsp(1)%chi(upfpsp(1)%mesh,upfpsp(2)%nwfc))
    allocate(upfpsp(1)%els(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%oc(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%lchi(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%nchi(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%rcut_chi(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%rcutus_chi(upfpsp(2)%nwfc))
    allocate(upfpsp(1)%epseu(upfpsp(2)%nwfc))


    do k=1,upfpsp(2)%nwfc
        upfpsp(1)%chi(:,k)=0.d0
        do ii=1,min_mesh
            upfpsp(1)%chi(ii,k) = upfpsp(2)%chi(ii,k)
        enddo
        upfpsp(1)%els(k) = upfpsp(2)%els(k)
        upfpsp(1)%oc(k) = upfpsp(2)%oc(k)
        upfpsp(1)%lchi(k) = upfpsp(2)%lchi(k)
        upfpsp(1)%nchi(k) = upfpsp(2)%nchi(k)
        upfpsp(1)%rcut_chi(k) = upfpsp(2)%rcut_chi(k)
        upfpsp(1)%rcutus_chi(k) = upfpsp(2)%rcutus_chi(k)
        upfpsp(1)%epseu(k) = upfpsp(2)%epseu(k)
    enddo

    write(*,*)
    write(*,*) "write new psps in:"
    upfname=adjustl(trim(upfname))//".REP"
    write(*,*) adjustl(trim(upfname))
    open(100,file=upfname)
    call write_upf_v2(100,upfpsp(1))
    close(100)

    deallocate(upfpsp)

end program main_replace_ppchi
