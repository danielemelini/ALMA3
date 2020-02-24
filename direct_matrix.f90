! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'direct_matrix'
! computes the fundamental matrix in a layer
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine direct_matrix(n,r,rho,mu,gra,Y)
use fmzm
use general_parameters
implicit none
!
type(fm) :: Y(6,6)
integer  :: n
type(fm) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
type(fm) :: r
type(fm) :: rho
type(fm) :: mu
type(fm) :: gra
!
!
 a1 = to_fm('2')*to_fm(n) + to_fm('3')                   ! 2n+3
 a2 = to_fm(n)+to_fm('1')                                ! n+1
 a3 = to_fm(n)+to_fm('3')                                ! n+3
 a4 = to_fm(n)**2 - to_fm(n) - to_fm('3')                ! n^2-n-3
 a5 = to_fm(n)+to_fm('2')                                ! n+2
 a6 = to_fm(n) - to_fm('1')                              ! n-1
 a7 = to_fm('2') * to_fm(n) + to_fm('1')                 ! 2n+1
 a8 = to_fm('2') * to_fm(n) - to_fm('1')                 ! 2n-1
 a9 = to_fm('2') - to_fm(n)                              ! 2-n
 a10= to_fm(n)**2 + to_fm('3') * to_fm(n) - to_fm('1')   ! n^2+3n-1
 a11= to_fm(n)**2 - to_fm('1')                           ! n^2-1
!
 Y = to_fm('0.0');
!
!
!
 Y(1,1) = to_fm(n)/(to_fm('2') * a1) * r**(n+1)
 Y(2,1) = a3 / ( to_fm('2') * a1 * a2 ) * r**(n+1) 
 Y(3,1) = ( to_fm(n) * rho * gra * r + 2 * a4 * mu ) / ( to_fm('2') * a1 ) * r**n
 Y(4,1) = to_fm(n) * a5 / ( a1 * a2 ) * mu * r**n
 Y(6,1) = to_fm('2') * pi * G * rho * to_fm(n) / a1 * r**(n+1)
!
 Y(1,2) = r**(n-1)
 Y(2,2) = to_fm('1') / to_fm(n) * r**(n-1)
 Y(3,2) = ( rho * gra * r + to_fm('2') * a6 * mu ) * r**(n-2)
 Y(4,2) = to_fm('2') * a6 / to_fm(n) * mu * r**(n-2)
 Y(6,2) = to_fm('4') * pi * G * rho * r**(n-1)
!
 Y(3,3) = rho * r**n
 Y(5,3) = r**n
 Y(6,3) = a7 * r**(n-1)
!
 Y(1,4) = a2 / ( to_fm('2') * a8 ) * r**(-n)
 Y(2,4) = a9 / ( to_fm('2') * to_fm(n) * a8 ) * r**(-n)
 Y(3,4) = ( a2 * rho * gra * r - to_fm('2') * a10 * mu ) / ( to_fm('2') * a8 ) * r**(-n-1)
 Y(4,4) = a11 / ( to_fm(n) * a8 ) * mu * r**(-n-1)
 Y(6,4) = to_fm('2') * pi * G * rho * a2 / a8 * r**(-n)
!
 Y(1,5) = r**(-n-2)
 Y(2,5) = - to_fm('1') / a2 * r**(-n-2)
 Y(3,5) = ( rho * gra * r - to_fm('2') * a5 * mu ) * r**(-n-3)
 Y(4,5) = to_fm('2') * a5 / a2 * mu * r**(-n-3)
 Y(6,5) = to_fm('4') * pi * G * rho * r**(-n-2)
!
 Y(3,6) = rho * r**(-n-1)
 Y(5,6) = r**(-n-1)

!
end subroutine direct_matrix