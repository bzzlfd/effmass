program convert_wg2rho
    implicit double precision (a-h,o-z)
    real*8, parameter :: A_AU_1=0.52917721067d0
    integer nstep,nstate
    real*8 AL(3,3),AL_t(3,3)
    real*8,allocatable,dimension(:)::E1
    complex*16,allocatable,dimension(:)::cc_int,cphase
    complex*16 cc

    integer,allocatable,dimension(:,:) :: ngtotnod_9,ngtotnod_9_t
    real*8, allocatable,dimension(:) :: gkk_n_tmp,gkk_n_xtmp, gkk_n_ytmp,gkk_n_ztmp
    real*8, allocatable,dimension(:) :: gkk,gkk_x, gkk_y,gkk_z
    complex*16,allocatable,dimension(:) :: ug_n_tmp,ug_one
    complex*16,allocatable,dimension(:,:) :: ug_tot
    real*8,allocatable,dimension(:,:,:) :: wave_r
    real*8,allocatable,dimension(:,:,:) :: wave_i
    real*8,allocatable,dimension(:,:,:) :: rho
    real*8,allocatable,dimension(:) :: vr_tmp
    integer,allocatable,dimension(:,:,:) :: index
    real*8,allocatable,dimension(:) :: wrk
    real*8,allocatable,dimension(:,:) :: xatom,xyz
    integer,allocatable,dimension(:) ::  iz
    character*20 filename,filename2

    write(6,*) "input the number of wave functions to construct rho"
    read(5,*) num_wave

    do 2000 iiwave=1,num_wave

    write(6,*) "input the prefactor for ", iiwave, "wavefunction"

    read(5,*) wave_factor

    write(6,*) "the following is the inform for", iiwave, "th wavefunction"


    open(11,file="OUT.GKK",form="unformatted")
    rewind(11)
    read(11) n1,n2,n3,mg_nx,nnodes,nkpt,is_SO,islda
    read(11) Ecut
    read(11) AL
    AL=AL/A_AU_1
    allocate(ngtotnod_9(nnodes,nkpt),ngtotnod_9_t(nnodes,nkpt))
    read(11) nnodes, ((ngtotnod_9(inode,kpt),inode=1,nnodes),kpt=1,nkpt)

    write(6,*) "there are ", nkpt, " kpoints"
    write(6,*) "input the ikpt to plot wg"
    read(5,*) ikpt_p

    if(.not.allocated(gkk_n_tmp)) then
    allocate(gkk_n_tmp(mg_nx))
    allocate(gkk_n_xtmp(mg_nx))
    allocate(gkk_n_ytmp(mg_nx))
    allocate(gkk_n_ztmp(mg_nx))

    allocate(gkk(mg_nx*nnodes))
    allocate(gkk_x(mg_nx*nnodes))
    allocate(gkk_y(mg_nx*nnodes))
    allocate(gkk_z(mg_nx*nnodes))
    endif

    do kpt=1,ikpt_p  
        num=0
        do inode=1,nnodes
            read(11) gkk_n_tmp
            read(11) gkk_n_xtmp
            read(11) gkk_n_ytmp
            read(11) gkk_n_ztmp
            if(is_SO==0) then
                do ig=1,ngtotnod_9(inode,kpt)
                    gkk(num+ig)=gkk_n_tmp(ig)
                    gkk_x(num+ig)=gkk_n_xtmp(ig)
                    gkk_y(num+ig)=gkk_n_ytmp(ig)
                    gkk_z(num+ig)=gkk_n_ztmp(ig)
                enddo
                num=num+ngtotnod_9(inode,kpt)
            else
                !just read half of the gkk,gkk_x,gkk_y,gkk_z, 
                !the down part are the same with the up
                do ig=1,ngtotnod_9(inode,kpt)/2
                    gkk(num+ig)=gkk_n_tmp(ig)
                    gkk_x(num+ig)=gkk_n_xtmp(ig)
                    gkk_y(num+ig)=gkk_n_ytmp(ig)
                    gkk_z(num+ig)=gkk_n_ztmp(ig)
                enddo
                num=num+ngtotnod_9(inode,kpt)/2
                !
            endif
        enddo
        if(is_SO==0) then
            ng_tot=num
        else
            ng_tot=num*2
        endif
    enddo
    close(11)

    write(6,*) "input the name of WG file"
    read(5,*) filename
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccco
    open(11,file=filename,form="unformatted")
    rewind(11)
    read(11) n1_t,n2_t,n3_t,mx,mg_nx,nnodes_t,nkpt_t,is_SO,islda
    read(11) Ecut_t
    read(11) AL_t
    read(11) nnodes, ((ngtotnod_9_t(inode,kpt),inode=1,nnodes),kpt=1,nkpt)

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

    !if(.not.allocated(ug_n_tmp)) then
    allocate(ug_n_tmp(mg_nx))
    allocate(ug_one(ng_tot))
    !endif

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
                        if(is_SO==0) then
                            do ig=1,ngtotnod_9(inode,kpt)
                                ug_one(num+ig)=ug_n_tmp(ig)
                            enddo
                            num=num+ngtotnod_9(inode,kpt)
                        else
                            !ug_one(1:ng_tot/2) -- up
                            !ug_one(ng_tot/2+1:ng_tot) -- dn
                            do ig=1,ngtotnod_9(inode,kpt)/2
                                ug_one(num+ig)=ug_n_tmp(ig)
                                ug_one(num+ig+ng_tot/2)=ug_n_tmp(ig+ngtotnod_9(inode,kpt)/2)
                            enddo
                            num=num+ngtotnod_9(inode,kpt)/2
                        endif
                    endif
                enddo
            enddo
        endif
    end do
    close(11)

    deallocate(ug_n_tmp)

    vol=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+ AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+ AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
    vol=dabs(vol)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccc
    pi=4*datan(1.d0)
    sum_chg=0.d0
    sum_chr=0.d0
    if(is_SO==0) then
        if(.not.allocated(wave_r)) then
            allocate(wave_r(n1,n2,n3))
            allocate(wave_i(n1,n2,n3))
            allocate(index(n1,n2,n3))
        endif
        wave_r=0.d0
        wave_i=0.d0
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
        enddo

        do k=1,n3
            do j=1,n2
                do i=1,n1
                    sum_chg=sum_chg+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                enddo
            enddo
        enddo

        lwrk=6*(n1+n2+n3)+15

        if(.not.allocated(wrk)) then
            allocate(wrk(lwrk))
        endif
        call cfft(n1,n2,n3,wave_r,wave_i,wrk,lwrk,1)
        !cccccc To use cfft, need to compile and link: cfft.f and cfftd.f in the directory

        if(.not.allocated(rho)) then
            allocate(rho(n1,n2,n3))
            rho=0.d0
        endif

        do k=1,n3
            do j=1,n2
                do i=1,n1
                    sum_chr=sum_chr+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                    rho(i,j,k)=rho(i,j,k)+wave_factor*(wave_r(i,j,k)**2+wave_i(i,j,k)**2)
                enddo
            enddo
        enddo
    else
        do ispinor=1,2
            if(.not.allocated(wave_r)) then
                allocate(wave_r(n1,n2,n3))
                allocate(wave_i(n1,n2,n3))
                allocate(index(n1,n2,n3))
            endif
            wave_r=0.d0
            wave_i=0.d0
            index=0 
            do ig=1,ng_tot/2
                aa1=al(1,1)*gkk_x(ig)+al(2,1)*gkk_y(ig)+al(3,1)*gkk_z(ig)
                aa2=al(1,2)*gkk_x(ig)+al(2,2)*gkk_y(ig)+al(3,2)*gkk_z(ig)
                aa3=al(1,3)*gkk_x(ig)+al(2,3)*gkk_y(ig)+al(3,3)*gkk_z(ig)
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
                    write(6,*) "something wrong,two ig,same i1,i2,i3", aa1,aa2,aa3,i1,i2,i3
                    stop
                endif
                index(i1,i2,i3)=1
                wave_r(i1,i2,i3)=real(ug_one(ig+(ispinor-1)*ng_tot/2))
                !c           wave_i(i1,i2,i3)=-aimag(ug_one(ig))
                wave_i(i1,i2,i3)=aimag(ug_one(ig+(ispinor-1)*ng_tot/2))   ! need to test more, but i believe this is correct
            enddo

            do k=1,n3
                do j=1,n2
                    do i=1,n1
                        sum_chg=sum_chg+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                    enddo
                enddo
            enddo

            lwrk=6*(n1+n2+n3)+15

            if(.not.allocated(wrk)) then
                allocate(wrk(lwrk))
            endif
            call cfft(n1,n2,n3,wave_r,wave_i,wrk,lwrk,1)
            !cccccc to use cfft, need to compile and link: cfft.f and cfftd.f in the directory

            if(.not.allocated(rho)) then
                allocate(rho(n1,n2,n3))
                rho=0.d0
            endif

            do k=1,n3
                do j=1,n2
                    do i=1,n1
                        sum_chr=sum_chr+wave_r(i,j,k)**2+wave_i(i,j,k)**2
                        rho(i,j,k)=rho(i,j,k)+wave_factor*(wave_r(i,j,k)**2+wave_i(i,j,k)**2)
                    enddo
                enddo
            enddo
        enddo
    endif
    sum_chg=sum_chg*vol
    sum_chr=sum_chr*vol/(n1*n2*n3)
    write(6,*) "charge(g-space)=",sum_chg
    write(6,*) "charge(R-space)=",sum_chr
    !cccccccccccccc  Now, we have the wave function in real space, 
    !cccccccccccccc we can do whatever we wanted for this wave function
    !cccccccccccccc Here, we will output it for plotting
    deallocate(ug_one)
2000 continue
     

    AL=AL*A_AU_1

    open(12,file="OUT.WG2RHO",form="unformatted")
    rewind(12)
    nstate=1
    write(12) n1,n2,n3,nnodes,nstate
    write(12) AL

    nr=n1*n2*n3
    nr_n=nr/nnodes
    allocate(vr_tmp(nr_n))

    do iread=1,nnodes
    do ii=1,nr_n
    jj=ii+(iread-1)*nr_n
    i=(jj-1)/(n2*n3)+1
    j=(jj-1-(i-1)*n2*n3)/n3+1
    k=jj-(i-1)*n2*n3-(j-1)*n3

    vr_tmp(ii)=rho(i,j,k)
    enddo
    write(12) vr_tmp
    enddo

    close(12)

    write(6,*) "the constructed rho is in OUT.WG2RHO" 

    stop
end program convert_wg2rho
