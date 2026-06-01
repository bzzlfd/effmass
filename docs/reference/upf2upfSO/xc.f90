      module xc
      !use data
      implicit none
!      !real(kind=8) :: rho_kin(mr_n,2)
      logical :: potential_mix
      integer :: xctype(3)
      character(512) :: xcdescr
      integer :: xcspin,xcgrad
      logical :: ishybrid, ismgga
!      real(8) :: hybridc
      end module xc
