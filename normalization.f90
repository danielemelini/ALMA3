subroutine normalization
use rheological_profile
use general_parameters
implicit none
!
 integer :: i,j,k,n
!
!
!
! ------- Reference scales
!
 r0     = r(nla+2)
 rho0   = rho(1) 
 mu0    = mu(1)
 t0     = to_fm('1000') * to_fm('365.25') * to_fm('24') * to_fm('3600')
 eta0   = mu0*t0
 mass0  = rho0 * r0**3
!
! 
! ------- Normalizes model parameters
 r      = r   / r0
 rho    = rho / rho0
 mu     = mu  / mu0
 eta    = eta / eta0
!
!
! ------- Normalized Newton constant
 G      = Gnt * rho0**2 * r0**2 / mu0
!
!
! ------- Normalized mass
 mass   = mass / mass0;
 mlayer = mlayer / mass0;
!
!
! ------- Normalized gravity at internal interfaces
 gra    = gra * r0**2 * (G/Gnt) / mass0
! 
!
!
end subroutine normalization
!
!
