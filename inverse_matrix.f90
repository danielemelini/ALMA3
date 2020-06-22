! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'inverse_matrix'
! computes the inverse of the fundamental matrix in a layer
!
! Initial version DM February 24, 2020
! Modified by DM June 16, 2020 - Converted to type(zm) for complex LNs
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine inverse_matrix(n,r,rho,mu,gra,Yinv)
use fmzm
use general_parameters
implicit none
!
type(zm) :: Yinv(6,6)
type(zm) :: D(6,6)
integer  :: n
type(fm) :: r
type(fm) :: rho
type(zm) :: mu
type(fm) :: gra
!
type(fm) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
!
integer  :: i
!
!
 a1 = to_fm('2') * to_fm(n) + to_fm('1')                 ! 2n+1
 a2 = to_fm('2') * to_fm(n) - to_fm('1')                 ! 2n-1
 a3 = to_fm(n) + to_fm('1')                              ! n+1
 a4 = to_fm(n) - to_fm('1')                              ! n-1
 a5 = to_fm(n) + to_fm('2')                              ! n+2
 a6 = to_fm(n) + to_fm('3')                              ! n+3
 a7 = to_fm(n)**2 + to_fm('3') * to_fm(n) - to_fm('1')   ! n^2+3n-1
 a8 = to_fm(n)**2 - to_fm(n) - to_fm('3')                ! n^2-n-3
 a9 = to_fm('2') - to_fm(n)                              ! 2-n
 a10= to_fm(n)**2 - to_fm('1')                           ! n^2-1
 a11= to_fm('2') * to_fm(n) + to_fm('3')                 ! 2n+3
!
 Yinv = to_zm('0.0');
 D    = to_fm('0.0');
! 
 D(1,1) = a3 * r**(-(n+1))
 D(2,2) = to_fm(n) * a3 / ( to_fm('2') * a2 ) * r**(-(n-1))
 D(3,3) = -r**(-(n-1))
 D(4,4) = to_fm(n) * r**n
 D(5,5) = to_fm(n) * a3 / ( to_fm('2') * a11 ) * r**(n+2)
 D(6,6) = r**(n+1)
!
 do i=1,6
    D(i,i) = D(i,i) / a1	
 end do
!
 Yinv(1,1) =   rho * gra * r / mu - to_fm('2') * a5
 Yinv(2,1) = - rho * gra * r / mu + to_fm('2') * a7 / a3
 Yinv(3,1) =   4 * pi * G * rho
 Yinv(4,1) =   rho * gra * r / mu + to_fm('2') * a4
 Yinv(5,1) = - rho * gra * r / mu - to_fm('2') * a8 / to_fm(n)
 Yinv(6,1) =   4 * pi * G * rho * r
!
 Yinv(1,2) =   to_fm('2') * to_fm(n) * a5
 Yinv(2,2) = - to_fm('2') * a10
 Yinv(4,2) =   to_fm('2') * a10
 Yinv(5,2) = - to_fm('2') * to_fm(n) * a5
!
 Yinv(1,3) = - r / mu
 Yinv(2,3) =   r / mu
 Yinv(4,3) = - r / mu
 Yinv(5,3) =   r / mu
!
 Yinv(1,4) =   to_fm(n) * r / mu
 Yinv(2,4) =   a9 * r / mu
 Yinv(4,4) = - a3 * r / mu 
 Yinv(5,4) =   a6 * r / mu
!
 Yinv(1,5) =   rho * r / mu
 Yinv(2,5) = - rho * r / mu
 Yinv(4,5) =   rho * r / mu
 Yinv(5,5) = - rho * r / mu
 Yinv(6,5) =   a1
!
 Yinv(3,6) = - to_fm('1')
 Yinv(6,6) = - r
!
 Yinv = matmul(D,Yinv)
!
end subroutine inverse_matrix