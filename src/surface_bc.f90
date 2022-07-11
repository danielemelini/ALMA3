! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'surface_bc'
! computes the boundary conditions at the surface
!
! Initial version DM February 24, 2020
! Fixed DM June 18, 2021       - Added 'save' to local variables
!                                to avoid memory leaks (see FMLIB manual)
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine surface_bc(n,r,gra,bs)
use fmzm
use general_parameters
implicit none
!
type(fm) :: bs(3)
integer  :: n
type(fm) :: r
type(fm) :: gra
!
type(fm), save :: kappa
!
!
!
 bs = to_fm('0.0');
!
 kappa = ( to_fm('2') * to_fm(n) + to_fm('1') ) / ( to_fm('4') * pi * r**2 ) 
!
 if ( iload == 1 )    bs(1) = - gra * kappa
!
 bs(3) = - to_fm('4') * pi * G * kappa
!
!
!
end subroutine surface_bc
