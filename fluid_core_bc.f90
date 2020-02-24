! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'fluid_core_bd'
! computes the boundary conditions at the CMB for a fluid inviscid core
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine fluid_core_bc(n,r,rho,gra,b)
use fmzm
use general_parameters
implicit none
!
type(fm) :: b(6,3)
integer  :: n
type(fm) :: r
type(fm) :: rho
type(fm) :: mu
type(fm) :: gra
!
!
!
 b = to_fm('0.0');
!
!
 b(1,1) = - to_fm('1') / gra * (r**n)
 b(1,3) =   to_fm('1')
!
 b(2,2) =   to_fm('1')
!
 b(3,3) =   rho * gra
!
 b(5,1) =   r**n
!
 b(6,1) =   to_fm('2') * ( to_fm(n) - to_fm('1') ) * r**(n-1)
 b(6,3) =   to_fm('4') * pi * G * rho
!
!
end subroutine fluid_core_bc