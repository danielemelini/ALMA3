! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'complex_rigidity'
! computes mu(s) for various rheologies
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine complex_rigidity(s,mu,eta,code,mu_s)
use fmzm
implicit none
!
type(fm) :: s
type(fm) :: mu
type(fm) :: eta
integer  :: code
type(fm) :: mu_s
!
!
!
 if( code==1 ) then         ! Elastic
    mu_s = mu
 elseif( code==2 ) then     ! Maxwell
    mu_s = mu * s / ( s + mu/eta )
 elseif( code==3 ) then     ! Newton
    mu_s = eta * s
 else
    write(*,*) ' ERROR: Invalid rheology (code=',code,').'
    stop
 end if
!
!
end subroutine complex_rigidity
