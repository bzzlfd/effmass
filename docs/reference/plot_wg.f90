program plot_wg
    implicit double precision (a-h,o-z)
    real*8, parameter :: A_AU_1=0.52917721067d0
    integer nstep,nstate
    real*8 AL(3,3),AL_t(3,3),ALI(3,3)
    real*8,allocatable,dimension(:)::E1
    complex*16,allocatable,dimension(:)::cc_int,cphase
    complex*16 cc

    integer ngtotnod_9(100,8000),ngtotnod_9_t(100,8000)
    real*8, allocatable,dimension(:) :: gkk_n_tmp,gkk_n_xtmp, gkk_n_ytmp,gkk_n_ztmp
    real*8, allocatable,dimension(:) :: gkk,gkk_x, gkk_y,gkk_z
    complex*16,allocatable,dimension(:) :: ug_n_tmp,ug_one,ug_two
    complex*16,allocatable,dimension(:,:) :: ug_tot
    real*8,allocatable,dimension(:,:,:) :: wave_r,wave_r_dn
    real*8,allocatable,dimension(:,:,:) :: wave_i,wave_i_dn
    integer,allocatable,dimension(:,:,:) :: index
    real*8,allocatable,dimension(:) :: wrk
    real*8,allocatable,dimension(:,:) :: xatom,xyz
    integer,allocatable,dimension(:) ::  iz
    character*20 filename,filename2
    real*8,allocatable,dimension(:,:) :: ak
    real*8,allocatable,dimension(:) :: akx,aky,akz
    integer:: n1,n2,n3


    open(72,file="OUT.KPT")
    rewind(72)
    read(72,*) nkpt
    allocate(ak(3,nkpt))
    allocate(akx(nkpt))
    allocate(aky(nkpt))
    allocate(akz(nkpt))
    read(72,*) iflag
    if(iflag.ne.2) then
      write(6,*) "iflag.ne.2,in OUT.KPT"
      stop
    endif
    do i1=1,nkpt
    read(72,*) ak(1,i1),ak(2,i1),ak(3,i1)   
    enddo
    close(72)




    open(11,file="OUT.GKK",form="unformatted")
    rewind(11)
    read(11) n1,n2,n3,mg_nx,nnodes,nkpt,is_SO,islda
    read(11) Ecut
    read(11) AL
    AL=AL/A_AU_1
    read(11) nnodes, ((ngtotnod_9(inode,kpt),inode=1,nnodes),kpt=1,nkpt)

    if(is_SO.eq.1) then
    ngtotnod_9=ngtotnod_9/2
    endif

    write(6,*) "there are ", nkpt, " kpoints"
    write(6,*) "input the ikpt to plot wg"
    read(5,*) ikpt_p

    allocate(gkk_n_tmp(mg_nx))
    allocate(gkk_n_xtmp(mg_nx))
    allocate(gkk_n_ytmp(mg_nx))
    allocate(gkk_n_ztmp(mg_nx))

    allocate(gkk(mg_nx*nnodes))
    allocate(gkk_x(mg_nx*nnodes))
    allocate(gkk_y(mg_nx*nnodes))
    allocate(gkk_z(mg_nx*nnodes))

    
    call get_ALI(AL,ALI)
    pi=4*datan(1.d0)

    do i1=1,nkpt
    akx(i1)=2*pi*(ALI(1,1)*ak(1,i1)+ALI(1,2)*ak(2,i1)+ALI(1,3)*ak(3,i1))
    aky(i1)=2*pi*(ALI(2,1)*ak(1,i1)+ALI(2,2)*ak(2,i1)+ALI(2,3)*ak(3,i1))
    akz(i1)=2*pi*(ALI(3,1)*ak(1,i1)+ALI(3,2)*ak(2,i1)+ALI(3,3)*ak(3,i1))
    enddo


    do kpt=1,ikpt_p  
        num=0
        do inode=1,nnodes
            read(11) gkk_n_tmp
            read(11) gkk_n_xtmp
            read(11) gkk_n_ytmp
            read(11) gkk_n_ztmp
            do ig=1,ngtotnod_9(inode,kpt)
                gkk(num+ig)=gkk_n_tmp(ig)
                gkk_x(num+ig)=gkk_n_xtmp(ig)+akx(kpt)
                gkk_y(num+ig)=gkk_n_ytmp(ig)+aky(kpt)
                gkk_z(num+ig)=gkk_n_ztmp(ig)+akz(kpt)
            enddo
            num=num+ngtotnod_9(inode,kpt)
        enddo
        ng_tot=num
    enddo
    close(11)

    write(6,*) "input the name of OUT.WG file"
    read(5,*) filename
    write(6,*) "input the name of xatom.config file"
    read(5,*) filename2
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccco
    open(11,file=filename,form="unformatted")
    rewind(11)
    read(11) n1_t,n2_t,n3_t,mx,mg_nx,nnodes_t,nkpt_t,is_SO,islda
    read(11) Ecut_t
    read(11) AL_t
    read(11) nnodes, ((ngtotnod_9_t(inode,kpt),inode=1,nnodes),kpt=1,nkpt)

    if(is_SO.eq.1) then
    ngtotnod_9_t=ngtotnod_9_t/2
    endif

    AL_t=AL_t/A_AU_1
    sum=0.d0
    do i=1,3
        do j=1,3
            sum=sum+abs(AL(i,j)-AL_t(i,j))
        enddo
    enddo
    if(sum.gt.0.01) then
        write(6,*) "AL in OUT.GKK and ",filename, "not the same, stop"
        stop
    endif

    sum=0.d0
    do kpt=1,nkpt
        do inode=1,nnodes
            sum=sum+abs(ngtotnod_9(inode,kpt)-ngtotnod_9_t(inode,kpt))
        enddo
    enddo
    if(sum.gt.0.01) then
        write(6,*) "ngtotnod in OUT.GG and ",filename,"not thesame,stop"
        stop
    endif

    if(n1.ne.n1_t.or.n2.ne.n2_t.or.n3.ne.n3_t.or. nnodes.ne.nnodes_t.or.nkpt.ne.nkpt_t.or. abs(Ecut-Ecut_t).gt.0.001) then
        write(6,*) "Param in OUT.GKK and ",filename, "not the same"
        write(6,*) n1,n2,n3,nnodes,nkpt,Ecut
        write(6,*) n1_t,n2_t,n3_t,nnodes_t,nkpt_t,Ecut_t
        write(6,*) "stop"
        stop
    endif

    write(6,*) "there are ",mx, "wavefunction, input ind im to plot" 
    read(5,*) im_p

    allocate(ug_n_tmp(mg_nx))
    allocate(ug_one(ng_tot))
    if(is_SO.eq.1) then
      allocate(ug_two(ng_tot))  ! the down component
    endif

    do kpt=1,ikpt_p  
        if(kpt.lt.ikpt_p) then
            do inode=1,nnodes
                do im=1,mx
                    read(11) 
                enddo
            enddo
        else
            num=0
            do inode=1,nnodes
                do im=1,mx
                    if(im.ne.im_p) then
                        read(11) 
                    else
                        read(11) ug_n_tmp
                        do ig=1,ngtotnod_9(inode,kpt)
                            ug_one(num+ig)=ug_n_tmp(ig)
                        enddo

                        if(is_SO.eq.1) then
                        do ig=1,ngtotnod_9(inode,kpt)
                            ug_two(num+ig)=ug_n_tmp(ngtotnod_9(inode,kpt)+ig)
                        enddo
                        endif
                        num=num+ngtotnod_9(inode,kpt)

                    endif
                enddo
            enddo
        endif
    end do
    close(11)

    vol=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+ AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+ AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
    vol=dabs(vol)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccc
    pi=4*datan(1.d0)

    allocate(wave_r(n1,n2,n3))
    allocate(wave_i(n1,n2,n3))
    wave_r=0.d0
    wave_i=0.d0
    if(is_SO.eq.1) then
    allocate(wave_r_dn(n1,n2,n3))
    allocate(wave_i_dn(n1,n2,n3))
    wave_r_dn=0.d0
    wave_i_dn=0.d0
    endif


    allocate(index(n1,n2,n3))
    index=0 
    do ig=1,ng_tot
        aa1=AL(1,1)*gkk_x(ig)+AL(2,1)*gkk_y(ig)+AL(3,1)*gkk_z(ig)
        aa2=AL(1,2)*gkk_x(ig)+AL(2,2)*gkk_y(ig)+AL(3,2)*gkk_z(ig)
        aa3=AL(1,3)*gkk_x(ig)+AL(2,3)*gkk_y(ig)+AL(3,3)*gkk_z(ig)
        aa1=aa1/(2*pi)+n1*2+0.1
        aa2=aa2/(2*pi)+n2*2+0.1
        aa3=aa3/(2*pi)+n3*2+0.1
        i1=aa1
        i2=aa2
        i3=aa3
        i1=mod(i1,n1)+1
        i2=mod(i2,n2)+1
        i3=mod(i3,n3)+1
        if(index(i1,i2,i3).gt.0) then
            write(6,*) "something wrong,two ig,samei1,i2,i3", aa1,aa2,aa3,i1,i2,i3
            stop
        endif
        index(i1,i2,i3)=1
        wave_r(i1,i2,i3)=real(ug_one(ig))
        !c           wave_i(i1,i2,i3)=-aimag(ug_one(ig))
        wave_i(i1,i2,i3)=aimag(ug_one(ig))   ! need to test more, but I believe this is correct
        if(is_SO.eq.1) then
        wave_r_dn(i1,i2,i3)=real(ug_two(ig))
        wave_i_dn(i1,i2,i3)=aimag(ug_two(ig))   ! need to test more, but I believe this is correct
        endif
    enddo

    sum=0.d0
    do k=1,n3
        do j=1,n2
            do i=1,n1
                sum=sum+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                if(is_SO.eq.1) then
                sum=sum+wave_r_dn(i,j,k)**2+wave_i_dn(i,j,k)**2
                endif
            enddo
        enddo
    enddo
    sum=sum*vol
    write(6,*) "charge(G-space)=",sum

    lwrk=6*(n1+n2+n3)+15
    allocate(wrk(lwrk))
    call cfft(n1,n2,n3,wave_r,wave_i,wrk,lwrk,1)
    !cccccc To use cfft, need to compile and link: cfft.f and cfftd.f in the directory
    if(is_SO.eq.1) then
    call cfft(n1,n2,n3,wave_r_dn,wave_i_dn,wrk,lwrk,1)
    endif

    sum=0.d0
    do k=1,n3
        do j=1,n2
            do i=1,n1
                sum=sum+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                if(is_SO.eq.1) then
                sum=sum+wave_r_dn(i,j,k)**2+wave_i_dn(i,j,k)**2
                endif
            enddo
        enddo
    enddo
    sum=sum*vol/(n1*n2*n3)
    write(6,*) "charge(R-space)=",sum
    !cccccccccccccc  Now, we have the wave function in real space, 
    !cccccccccccccc we can do whatever we wanted for this wave function
    !cccccccccccccc Here, we will output it for plotting
    AL=AL*A_AU_1
    open(10,file=filename2)
    rewind(10)
    read(10,*) natom
    allocate(xatom(3,natom))
    allocate(xyz(3,natom))
    allocate(iz(natom))
    read(10,*)
    read(10,*) AL_t(1,1),AL_t(2,1),AL_t(3,1)
    read(10,*) AL_t(1,2),AL_t(2,2),AL_t(3,2)
    read(10,*) AL_t(1,3),AL_t(2,3),AL_t(3,3)
    sum=0.d0
    do i=1,3
        do j=1,3
            sum=sum+abs(AL(i,j)-AL_t(i,j))
        enddo
    enddo
    if(sum.gt.0.01) then
        write(6,*) "AL in xatom.config and OUT.KGG not same, stop"
        stop
    endif

    read(10,*) 
    do iat=1,natom
        read(10,*) iz(iat),xatom(1,iat),xatom(2,iat),xatom(3,iat)
    enddo
    close(10)

    do iat=1,natom
        do j=1,3
            xatom(j,iat)=mod(xatom(j,iat)+2.d0,1.d0)
        enddo
    enddo

    do iat=1,natom
        do i=1,3
            xyz(i,iat)=AL(i,1)*xatom(1,iat)+AL(i,2)*xatom(2,iat)+ AL(i,3)*xatom(3,iat)
        enddo
    enddo

    open(12,file="PSI.xsf")
    rewind(12)
    write(12,'(x,a)') 'CRYSTAL'
    write(12,'(x,a)') 'PRIMVEC'
    do i = 1, 3
        write (12, '(3f20.10)') AL(1:3,i)
    end do
    write(12,'(x,a)') 'PRIMCOORD'
    write(12, *) natom, 1
    do iat=1,natom
        write(12, '(x,i4,x,3(f20.10,x))') iz(iat),xyz(1:3,iat)
    end do
    write(12,1001)"BEGIN_BLOCK_DATAGRID_3D"
    write(12,1001)"XSF_FILE"
    write(12,1001)" BEGIN_DATAGRID_3D_XSF_FILE"
    write(12,98) n1+1,n2+1,n3+1
    write(12,200) 0.0000,0.00000,0.00000
    write(12,200) AL(1,1),AL(2,1),AL(3,1)
    write(12,200) AL(1,2),AL(2,2),AL(3,2)
    write(12,200) AL(1,3),AL(2,3),AL(3,3)
    write(12,*)
    if(is_SO.eq.0) then
        write(12,100) (((wave_r(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2+&
            wave_i(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2,i=1,n1+1),j=1,n2+1),k=1,n3+1)
    else
    write(12,100)    &
        (((wave_r(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2+&
        wave_i(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2+&
        wave_r_dn(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2+&
        wave_i_dn(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1)**2,i=1,n1+1),j=1,n2+1),k=1,n3+1)
    endif

    write(12,1001)"END_DATAGRID_3D"
    write(12,1001)"END_BLOCK_DATAGRID_3D"
    close(12)
    98     format(3(i5,2x))
    200    format(3(1x,E13.5,1x))
    !100    format(6(F12.6))
    100    format(6(E14.7,1x))
    1001   format(a28)
    write(6,*) "PSI^2 is written in PSI.xsf"
    
    if(is_SO.eq.0) then
        call write_wr_xsf('REAL_PSI.xsf', wave_r)
        write(6,*) "real part of PSI is written in REAL_PSI.xsf"
        call write_wr_xsf('IMAG_PSI.xsf',wave_i)
        write(6,*) "imaginary part of PSI is written in IMAG_PSI.xsf"
    else
        call write_wr_xsf('REAL_PSI_UP.xsf', wave_r)
        write(6,*) "real part of PSI's up component is written in REAL_PSI_UP.xsf"
        call write_wr_xsf('IMAG_PSI_UP.xsf',wave_i)
        write(6,*) "imaginary part of PSI's up component is written in IMAG_PSI_UP.xsf"
        call write_wr_xsf('REAL_PSI_DOWN.xsf', wave_r_dn)
        write(6,*) "real part of PSI's down component is written in REAL_PSI_DOWN.xsf"
        call write_wr_xsf('IMAG_PSI_DOWN.xsf',wave_i_dn)
        write(6,*) "imaginary part of PSI's down component is written in IMAG_PSI_DOWN.xsf"
    endif

contains
 subroutine write_wr_xsf(filename,wr)
     character(len=*) :: filename
     real*8 :: wr(n1,n2,n3)
     open(12,file=adjustl(trim(filename)))
     rewind(12)
     write(12,'(x,a)') 'CRYSTAL'
     write(12,'(x,a)') 'PRIMVEC'
     do i = 1, 3
         write (12, '(3f20.10)') AL(1:3,i)
     end do
     write(12,'(x,a)') 'PRIMCOORD'
     write(12, *) natom, 1
     do iat=1,natom
         write(12, '(x,i4,x,3(f20.10,x))') iz(iat),xyz(1:3,iat)
     end do
     write(12,1001)"BEGIN_BLOCK_DATAGRID_3D"
     write(12,1001)"XSF_FILE"
     write(12,1001)" BEGIN_DATAGRID_3D_XSF_FILE"
     write(12,98) n1+1,n2+1,n3+1
     write(12,200) 0.0000,0.00000,0.00000
     write(12,200) AL(1,1),AL(2,1),AL(3,1)
     write(12,200) AL(1,2),AL(2,2),AL(3,2)
     write(12,200) AL(1,3),AL(2,3),AL(3,3)
     write(12,*)
     write(12,100) (((wr(mod(i-1,n1)+1,mod(j-1,n2)+1,mod(k-1,n3)+1),i=1,n1+1),j=1,n2+1),k=1,n3+1)

     write(12,1001)"END_DATAGRID_3D"
     write(12,1001)"END_BLOCK_DATAGRID_3D"
     close(12)
    98     format(3(i5,2x))
    200    format(3(1x,E13.5,1x))
    !100    format(6(F12.6))
    100    format(6(E14.7,1x))
    1001   format(a28)
 end subroutine write_wr_xsf
end program plot_wg
