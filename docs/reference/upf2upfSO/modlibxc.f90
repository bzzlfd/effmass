      module modlibxc

      use xc_f90_types_m

      implicit none

      type libxc_functional
        integer :: family ! unknown=-1, none=0, lda=1, gga=2, mgga=3
        integer :: id ! identifier
        integer :: nspin ! unplarized=1,polarized=2
        integer :: kind ! ex=0, co=1, ex_co=2, kinetic=3
        character(len=256) :: name ! functional name
        character(len=256) :: ref ! functional reference
        integer :: relativistic ! non=0, rel=1
        type(xc_f90_pointer_t) :: p  ! the pointer used to call the library
        type(xc_f90_pointer_t) :: info ! information about the functional
      end type 

      type(libxc_functional) :: funcs(2)

      end module modlibxc
