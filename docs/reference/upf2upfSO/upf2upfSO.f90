    !! 
    !!
    module fhi
	!
	type angular_comp
	    real(kind=8), allocatable :: pot(:), wfc(:), grid(:)
	    real(kind=8) :: amesh
	    integer :: nmesh
	    integer :: lcomp
	end type angular_comp
	!
	real(kind=8) :: zval
	integer :: lmax
	logical :: nlcc_
	real(kind=8), allocatable :: rho_atc(:)
	!
	type (angular_comp), allocatable :: comp(:)
	!
	integer :: zatom, lloc
    end module fhi
    !! 
    !! 
    !!
    !! 
    !! 
    !!
        subroutine upf_to_upfSO (upf)
	!
	use upf_module
!	use fhi
!	use vwr
	!
	implicit double precision (a-h,o-z)
	type(atompp) :: upf
	real(kind=8), allocatable :: aux(:)
	real(kind=8) :: vll,sum1,sum2
	character(len=2) :: label
	character(len=2), external :: atom_name
	integer :: i, l, ir, iv

       real*8 r(3000),beta_tmp(3000,14)
       real*8 dr(3000),vloc(3000)
       real*8 wave(3000,0:3,2),wave_tmp(3000,8)
       real*8 oc_tmp(8),occ(0:3,2)
 !ccc   wave(i,L,j=1,2)
        real*8 beta(3000,0:3,2,2)  !(i,L,j,m)
        real*8 beta0(3000,0:3,2)   ! (i,L,m)
        real*8 betaSO(3000,0:3)    !  (i,L)
        real*8 DD_0(0:3,2)   ! (L,m)
        real*8 DD_SO(0:3)    !  L
        real*8 DD(0:30,0:30), DD_new(30,30)
        integer map_J(20),map_L(20),map_Jw(20),map_Lw(20)
        integer ind(0:3,1:2),ll_map(0:3,1:2,2)

        !init zeros
        r=0.d0
        beta_tmp=0.d0
        dr=0.d0
        vloc=0.d0
        wave=0.d0
        wave_tmp=0.d0
        oc_tmp=0.d0
        occ=0.d0
        beta=0.d0
        beta0=0.d0
        betaSO=0.d0
        DD_0=0.d0
        DD_SO=0.d0
        DD=0.d0
        DD_new=0.d0
        map_J=0
        map_L=0
        map_Jw=0
        map_Lw=0
        ind=0
        ll_map=0
        
	!

       ! write(6,*) upf%nwfc
	upf%nv = "2.0.1"
	upf%generated = "from upf To upfSO"
	upf%typ = "NC"
        if(upf%has_so.ne..true.) then
        write(6,*) "the input must be a norm conserv. PP with SO coupling "
        stop
        endif
        upf%has_so = .false.

        nwave=upf%nwfc
        nbeta=upf%nbeta
        n=upf%mesh
    
        !if(upf%lmax.eq.1.and.nbeta.ne.6) then
        !write(6,*) "lmax.eq.1,nbeta.ne.6,stop",upf%lmax,nbeta
        !stop
        !endif

        !if(upf%lmax.eq.2.and.nbeta.ne.10) then
        !write(6,*) "lmax.eq.2,nbeta.ne.10,stop",upf%lmax,nbeta
        !stop
        !endif

        !if(upf%lmax.eq.3.and.nbeta.ne.14) then
        !write(6,*) "lmax.eq.3,nbeta.ne.14,stop",upf%lmax,nbeta
        !stop
        !endif
        if(upf%lmax>3) then
            write(6,*) "upf%lmax should .le. 3"
            stop
        endif

        if(nwave.eq.0) then
            write(6,*) "number_of_wfc in pseudopotential file must great than 0"
            stop
        endif

        r(1:n)=upf%r(1:n)
        do ll=1,nbeta
        beta_tmp(1:n,ll)=upf%beta(1:n,ll)
        enddo
        do ll=1,nwave
        wave_tmp(1:n,ll)=upf%chi(1:n,ll)
        oc_tmp(ll)=upf%oc(ll)
        enddo
        do j=1,nbeta
        do i=1,nbeta
        DD(i,j)=upf%dion(i,j)
        enddo
        enddo

        lchi_max=maxval(upf%lchi(1:nwave))
        do ll=1,nwave
        map_Lw(ll)=upf%lchi(ll)
        if(upf%jchi(ll)-upf%lchi(ll).lt.0) then
        map_Jw(ll)=1      ! J=1, means, J=L-1/2
        else
        map_Jw(ll)=2      ! J=1, means, J=L-1/2, for L=0, it is J=2
        endif
        enddo

        do ll=1,nbeta
        map_L(ll)=upf%lll(ll)
        if(upf%jjj(ll)-upf%lll(ll).lt.0) then
        map_J(ll)=1     ! J=1, mean, J=L-1/2
        else
        map_J(ll)=2     ! J=1, means, J=L-1/2
        endif
        enddo

        beta=0.d0
        ind=0
        do ll=1,nbeta
        ind(map_L(ll),map_J(ll))=ind(map_L(ll),map_J(ll))+1
        m=ind(map_L(ll),map_J(ll))
        beta(1:n,map_L(ll),map_J(ll),m)=beta_tmp(1:n,ll)
        ll_map(map_L(ll),map_J(ll),m)=ll
        enddo

        do i=0,3
            do j=1,2
                m=ind(i,j)
                if(m==1) then
                    beta(1:n,i,j,2)=beta(1:n,i,j,1)
                    !ll_map(i,j,2)=0
                endif
            enddo
        enddo

        wave=0.d0
        do ll=1,nwave
        wave(1:n,map_Lw(ll),map_Jw(ll))=wave_tmp(1:n,ll)
        occ(map_Lw(ll),map_Jw(ll))=oc_tmp(ll)
        enddo

       write(6,*) "original nbeta,nwfc,lmax",nbeta,nwave,upf%lmax
       DD_0=0.d0
       DD_SO=0.d0
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
       do 1000 L=1,upf%lmax

       sum1=0.d0
       sum2=0.d0
       sumw=0.d0
       sumw1=0.d0
       sumw2=0.d0
       do i=2,n-1
       sum1=sum1+beta(i,L,1,1)*beta(i,L,2,1)*(r(i+1)-r(i-1))/2
       sum2=sum2+beta(i,L,1,2)*beta(i,L,2,2)*(r(i+1)-r(i-1))/2
       sumw=sumw+wave(i,L,1)*wave(i,L,2)*(r(i+1)-r(i-1))/2
       !sumw=1.d0
       sumw1=sumw1+wave(i,L,1)**2*(r(i+1)-r(i-1))/2
       sumw2=sumw2+wave(i,L,2)**2*(r(i+1)-r(i-1))/2
       enddo

       if(sum1.gt.0) then
       isign1=1
       else
       isign1=-1
       endif
       if(sum2.gt.0) then
       isign2=1
       else
       isign2=-1
       endif
       if(sumw.gt.0) then
       isignw=1
       else
       isignw=-1
       endif
       write(6,*) "L,isign1,2,w",L,isign1,isign2,isignw


!cccccc sumJ_m
        sum1_1=0.d0
        sum2_1=0.d0
        sum1_2=0.d0
        sum2_2=0.d0

        sum1_1T=0.d0
        sum2_1T=0.d0
        sum1_2T=0.d0
        sum2_2T=0.d0

       do i=2,n-1
       x2=dsqrt(L+1.d0)/(dsqrt(L+1.d0)+dsqrt(L*1.d0))
       x1=1-x2
       waveT1=x2*beta(i,L,2,1)+x1*isign1*beta(i,L,1,1)   ! (i,L,J,m),m=1,2
       waveT2=x2*beta(i,L,2,2)+x1*isign2*beta(i,L,1,2)
       if(L.le.lchi_max) then
           waveT=x2*wave(i,L,2)+x1*isignw*wave(i,L,1)
       else
           waveT=x2*wave(i,lchi_max,2)+x1*isignw*wave(i,lchi_max,1)
       endif
!  waveT is used for the betaSO, which will include both m=1,2
!  waveT1 is used for the beta0 for reference m=1 state
!  waveT2 is used for the beta0 for reference m=2 state
!cc   some how average the L+1/2 and L-1/2 wave function
       sum1_1=sum1_1+beta(i,L,1,1)*waveT*(r(i+1)-r(i-1))/2
       sum2_1=sum2_1+beta(i,L,2,1)*waveT*(r(i+1)-r(i-1))/2
       sum1_2=sum1_2+beta(i,L,1,2)*waveT*(r(i+1)-r(i-1))/2
       sum2_2=sum2_2+beta(i,L,2,2)*waveT*(r(i+1)-r(i-1))/2

       sum1_1T=sum1_1T+beta(i,L,1,1)*waveT1*(r(i+1)-r(i-1))/2
       sum2_1T=sum2_1T+beta(i,L,2,1)*waveT1*(r(i+1)-r(i-1))/2
       sum1_2T=sum1_2T+beta(i,L,1,2)*waveT2*(r(i+1)-r(i-1))/2
       sum2_2T=sum2_2T+beta(i,L,2,2)*waveT2*(r(i+1)-r(i-1))/2
       enddo

        write(6,*) "sumj_m=beta(j,m)*wave(aveJ);sumj_mT=beta(j,m)*beta(aveJ,m)" 
        write(6,*) "Ref: wave(aveJ) for betaSO(L); beta(aveJ,m) for beta0(L,m)"

        write(6,*) "sum1_1, sum2_1  ",sum1_1,sum2_1
        write(6,*) "sum1_1T,sum2_1T ",sum1_1T,sum2_1T
        write(6,*) "sum1_2, sum2_2  ",sum1_2,sum2_2
        write(6,*) "sum1_2T,sum2_2T ",sum1_2T,sum2_2T

       sum1_1=sum1_1*DD(ll_map(L,1,1),ll_map(L,1,1))
       sum2_1=sum2_1*DD(ll_map(L,2,1),ll_map(L,2,1))
       sum1_2=sum1_2*DD(ll_map(L,1,2),ll_map(L,1,2))
       sum2_2=sum2_2*DD(ll_map(L,2,2),ll_map(L,2,2))
       sum1_1T=sum1_1T*DD(ll_map(L,1,1),ll_map(L,1,1))
       sum2_1T=sum2_1T*DD(ll_map(L,2,1),ll_map(L,2,1))
       sum1_2T=sum1_2T*DD(ll_map(L,1,2),ll_map(L,1,2))
       sum2_2T=sum2_2T*DD(ll_map(L,2,2),ll_map(L,2,2))

!cccc  for beta0(i,L,1),m=1, use beta itself as the test psi
       beta0(1:n,L,1)=sum1_1T*L/(2*L+1.d0)*beta(1:n,L,1,1)+  &
     &   (L+1.d0)/(2*L+1.d0)*sum2_1T*beta(1:n,L,2,1)

! cccc  for beta0(i,L,2),m=2, use beta itself as the test psi
       beta0(1:n,L,2)=sum1_2T*L/(2*L+1.d0)*beta(1:n,L,1,2)+  &
     &   (L+1.d0)/(2*L+1.d0)*sum2_2T*beta(1:n,L,2,2)

       betaSO(1:n,L)=2.d0/(2*L+1.d0)*(sum2_1*beta(1:n,L,2,1)-  &
     &    sum1_1*beta(1:n,L,1,1)+sum2_2*beta(1:n,L,2,2)-   &
     &    sum1_2*beta(1:n,L,1,2))
!ccccccccccccccccccccccc
       sum1=0.d0
       sum2=0.d0
       sumSO=0.d0

       sum1t=0.d0
       sum2t=0.d0
       sumSOt=0.d0
       do i=2,n-1
       x2=dsqrt(L+1.d0)/(dsqrt(L+1.d0)+dsqrt(L*1.d0))
       x1=1-x2
       if(L.le.lchi_max) then
           waveT=x2*wave(i,L,2)+x1*isignw*wave(i,L,1)
       else
           waveT=x2*wave(i,lchi_max,2)+x1*isignw*wave(i,lchi_max,1)
       endif
       waveT1=x2*beta(i,L,2,1)+x1*isign1*beta(i,L,1,1)
       waveT2=x2*beta(i,L,2,2)+x1*isign2*beta(i,L,1,2)

       sum1=sum1+beta0(i,L,1)*waveT1*(r(i+1)-r(i-1))/2
       sum2=sum2+beta0(i,L,2)*waveT2*(r(i+1)-r(i-1))/2
       sumSO=sumSO+betaSO(i,L)*waveT*(r(i+1)-r(i-1))/2

       sum1t=sum1t+beta0(i,L,1)**2*(r(i+1)-r(i-1))/2
       sum2t=sum2t+beta0(i,L,2)**2*(r(i+1)-r(i-1))/2
       sumSOt=sumSOt+betaSO(i,L)**2*(r(i+1)-r(i-1))/2
       enddo

       if(abs(sum1t).gt.1.D-30) then
           beta0(1:n,L,1)=beta0(1:n,L,1)/dsqrt(sum1t)
       endif
       if(abs(sum2t).gt.1.D-30) then
           beta0(1:n,L,2)=beta0(1:n,L,2)/dsqrt(sum2t)
       endif

       if(abs(sumSOt).gt.1.D-30) then
       betaSO(1:n,L)=betaSO(1:n,L)/dsqrt(sumSOt)
       endif

!ccccc DD_T(L,m)
!ccccc DD_SO(L)
       if(abs(sum1).gt.1.D-30) then
           DD_0(L,1)=sum1t/sum1
       endif
       if(abs(sum2).gt.1.D-30) then
           DD_0(L,2)=sum2t/sum2
       endif
       if(abs(sumSO).gt.1.D-30) then  ! sometime wave might not been defined
       DD_SO(L)=sumSOt/sumSO
       endif

       if(L.eq.1) then
       open(11,file="plot.1.SO")
       elseif(L.eq.2) then
       open(11,file="plot.2.SO")
       elseif(L.eq.3) then
       open(11,file="plot.3.SO")
       endif
       rewind(11)
       do i=1,n
       write(11,200) r(i), dsqrt(abs(DD_0(L,1)))*beta0(i,L,1),  &
      &  dsqrt(abs(DD_0(L,2)))*beta0(i,L,2), &
      &  dsqrt(abs(DD_SO(L)))*betaSO(i,L)
       enddo
200    format(6(E14.7,1x))
       close(11)

1000   continue


       sum=0.d0
       do i=1,nbeta
       do j=1,nbeta
       if(i.ne.j) then
       sum=sum+abs(DD(i,j))
       endif
       enddo
       enddo
       if(sum.gt.1.E-20) then
       write(6,*) "DD not diagonal, stop"
       stop
       endif


       !if(nbeta.ne.6.and.nbeta.ne.10.and.nbeta.ne.14) then
       !write(6,*) "nbeta.ne.6.or.10.or.14, stop", nbeta
       !stop
       !endif
       if(nbeta.ge.3) then
        !    if(map_L(1).ne.0.or.map_L(2).ne.0.or.map_L(3).eq.0) then
        !        write(6,*) "no two L=0 beta, stop", (map_L(i),i=1,nbeta)
        !        stop
        !    endif
       else if(nbeta>0 .and. nbeta<3) then
        !    if(map_L(1).ne.0.or.map_L(2).ne.0) then
        !        write(6,*) "no two L=0 beta, stop", (map_L(i),i=1,nbeta)
        !        stop
        !    endif
       else
           write(*,*) "num_of_proj must > 1"
       endif

!ccccccccccccccccccccccccccccccccccccccccccccccc
       !if (nbeta == 6) then
       if(upf%lmax==1) then
       DD_new=0.d0
       DD_new(1,1)=DD(1,1)   ! L=0,m=1
       DD_new(2,2)=DD(2,2)   ! L=0,m=2
       DD_new(3,3)=DD_0(1,1)    ! L=1,m=1
       DD_new(4,4)=DD_0(1,2)    ! L=1,m=2
       DD_new(5,5)=DD_SO(1)  
       upf%beta(1:n,3)=beta0(1:n,1,1)    ! L=1,m=1
       upf%beta(1:n,4)=beta0(1:n,1,2)    ! L=1,m=2
       upf%lll(3)=1
       upf%lll(4)=1
       upf%beta(1:n,5)=betaSO(1:n,1)
       upf%lll(5)=10
       nbeta_new=5
       end if

       !if (nbeta == 10 .or. nbeta == 14) then

        if(upf%lmax>1) then
       DD_new=0.d0
       DD_new(1,1)=DD(1,1)   ! L=0,m=1
       DD_new(2,2)=DD(2,2)   ! L=0,m=2
       DD_new(3,3)=DD_0(1,1)    ! L=1,m=1
       DD_new(4,4)=DD_0(1,2)    ! L=1,m=2
       DD_new(5,5)=DD_0(2,1)    ! L=2,m=1
       DD_new(6,6)=DD_0(2,2)    ! L=2,m=2
        upf%beta(1:n,3)=beta0(1:n,1,1)    ! L=1,m=1
        upf%beta(1:n,4)=beta0(1:n,1,2)    ! L=1,m=2
        upf%beta(1:n,5)=beta0(1:n,2,1)    ! L=2,m=1
        upf%beta(1:n,6)=beta0(1:n,2,2)    ! L=2,m=2
        upf%lll(3)=1
        upf%lll(4)=1
        upf%lll(5)=2
        upf%lll(6)=2

        !if(nbeta.eq.10) then   ! no-f-state
        if(upf%lmax==2) then
       DD_new(7,7)=DD_SO(1)     ! L=1
       DD_new(8,8)=DD_SO(2)     ! L=2
       upf%beta(1:n,7)=betaSO(1:n,1)
       upf%beta(1:n,8)=betaSO(1:n,2)
       upf%lll(7)=10     ! indicate it is spin-orbit
       upf%lll(8)=20
       nbeta_new=8
       !elseif(nbeta.eq.14) then   ! f-state
        elseif(upf%lmax==3) then   ! f-state
       DD_new(7,7)=DD_0(3,1)
       DD_new(8,8)=DD_0(3,2)
       DD_new(9,9)=DD_SO(1)     ! L=1
       DD_new(10,10)=DD_SO(2)     ! L=2
       DD_new(11,11)=DD_SO(3)     ! L=3   do not use f-SO term, PWmat do not have
!   this term. It is usually very small anyway
        upf%beta(1:n,7)=beta0(1:n,3,1)
        upf%beta(1:n,8)=beta0(1:n,3,2)
        upf%beta(1:n,9)=betaSO(1:n,1)
        upf%beta(1:n,10)=betaSO(1:n,2)
        upf%beta(1:n,11)=betaSO(1:n,3)
       upf%lll(7)=3   
       upf%lll(8)=3
       upf%lll(9)=10   
       upf%lll(10)=20
       upf%lll(11)=30
        nbeta_new=11
        write(*,*) "upflll=",upf%lll(1:11)
       endif

       end if

       nbeta=nbeta_new
       upf%nbeta=nbeta
       deallocate(upf%dion)
       allocate(upf%dion(nbeta,nbeta))
       do j=1,nbeta
       do i=1,nbeta
       upf%dion(i,j)=DD_new(i,j)    
       enddo
       enddo
!ccccccccccccccccccccccccccccccccccccccccccccccc
       ll=0
       do jj=1,nwave
       if(map_Lw(jj).eq.0) then
       ll=ll+1
       upf%chi(1:n,ll)=wave_tmp(1:n,jj)
       upf%lchi(ll)=0
       upf%oc(ll)=oc_tmp(jj)
       endif
       enddo
       
       do 2000 L=1,upf%lmax
       sumw=0.d0
       do i=2,n-1
       sumw=sumw+wave(i,L,1)*wave(i,L,2)*(r(i+1)-r(i-1))/2
       enddo
       if(sumw.gt.0) then
       isignw=1
       else
       isignw=-1
       endif
       sumt=0.d0
       x2=dsqrt(L+1.d0)/(dsqrt(L+1.d0)+dsqrt(L*1.d0))
       x1=1-x2
       do i=2,n-1
       waveT=x2*wave(i,L,2)+x1*isignw*wave(i,L,1)
       sumt=sumt+waveT**2*(r(i+1)-r(i-1))/2
       enddo

       if(sumt.gt.1.D-30) then
       ll=ll+1
       upf%chi(1:n,ll)=(x2*wave(1:n,L,2)+x1*isignw*wave(1:n,L,1))/dsqrt(sumt)
       upf%lchi(ll)=L
       upf%oc(ll)=occ(L,1)+occ(L,2)
       endif
2000   continue
       upf%nwfc=ll

       write(6,*) "NEW:nbeta,nwfc,lmax=",upf%nbeta,upf%nwfc,upf%lmax
      
         RETURN

    end subroutine upf_to_upfSO
    !! 
    !! 
    !! 
    program upf2upfSO
	use upf_module
	implicit none
	type(atompp) :: upf
	character*200 :: vwr_name
	integer :: stat
	integer :: iovwr = 1122, ioupf = 2233

	write ( 6, * ) "input the upf file name(with SO)"
	read ( 5, * ) vwr_name
	write ( 6, * ) "opening the file: "//trim(vwr_name)

	call read_upf_v2 ( vwr_name, upf, 0)  ! iflag=0, no change rho_atc,vloc
	close ( iovwr )
	!!
        call upf_to_upfSO (upf)

	write ( 6, * ) "output PP file in UPF format : ",trim(vwr_name)//".SO"
	open ( ioupf, file = trim(vwr_name)//".SO" )


	call write_upf_v2 ( ioupf, upf )
	close ( ioupf )
	write ( 6, * ) "Pseudopotential successfully written"
	!!
    end program upf2upfSO
