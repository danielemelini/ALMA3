! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'complex_rigidity'
! computes mu(s) for various rheologies
!
! Initial version DM February 24, 2020
! Modified by DM June 11, 2020 - Burgers and Andrade rheologies
! Modified by DM June 16, 2020 - Complex LNs  - converted to type(zm)
! Fixed DM June 18, 2021       - Added 'save' to local variables
!                                to avoid memory leaks (see FMLIB manual)
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine complex_rigidity(s,mu,eta,code,par,mu_s)
use fmzm
implicit none
!
type(zm) :: s
type(fm) :: mu
type(fm) :: eta
integer  :: code
type(fm) :: par(5)
type(zm) :: mu_s
!
type(fm), save :: eta2
type(fm), save :: mu2
type(fm), save :: alpha
type(fm), save :: gam
!
!
!
 if( code==1 ) then         ! Elastic
    mu_s = mu
 elseif( code==2 ) then     ! Maxwell
    mu_s = mu * s / ( s + mu/eta )
 elseif( code==3 ) then     ! Newton
    mu_s = eta * s
 elseif( code==4 ) then     ! Kelvin
    mu_s = mu + eta * s
 elseif( code==5 ) then     ! Burgers
    mu2  = par(1)*mu
    eta2 = par(2)*eta
    mu_s = mu * s * ( s + mu2/eta2 ) / &
	       ( s**2 + s*(mu/eta + (mu+mu2)/eta2) + (mu*mu2)/(eta*eta2) )   
 elseif( code==6 ) then     ! Andrade
    alpha = par(1)
	gam   = par(2)
    mu_s  = 1/mu + 1/(eta*s) + gam * (1/mu) * (s*eta/mu)**(-alpha)
	mu_s  = 1/mu_s
 else
    write(*,*) ' ERROR: Invalid rheology (code=',code,').'
    stop
 end if
!
!
end subroutine complex_rigidity
