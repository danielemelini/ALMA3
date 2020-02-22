!
module general_parameters
use fmzm
implicit none 
!
integer :: order                    ! Order of the Gaver sequence
integer :: nd                       ! Numerical precision
!
type(fm) :: pi                      ! Pi
type(fm) :: G                       ! Normalized Newton constant
type(fm) :: Gnt                     ! Un-normalized Newton constant    
!
integer  :: iload                   ! Load type ('0'->tidal, '1'->loading)
integer  :: itime                   ! Time scale type ('0'->linear, '1'->log)
integer  :: lmin, lmax, lstep       ! Min/max harmonic degree and step
!
type(fm), allocatable :: t(:)       ! Time steps
!
real(4) :: m1,m2                    ! Time range is 10^(m1:m2)
integer :: p                        ! Number of times minus one
!
type(fm), allocatable :: zeta(:)    ! Salzer weights
!
character(100) :: file_rheol        ! File with rheological model
!
end module
!
module rheological_profile
use fmzm
implicit none
!
integer :: nla                      ! Number of mantle layers
!
type(fm), allocatable :: r(:)       ! Normalized radii of interfaces
type(fm), allocatable :: gra(:)     ! Normalized gravity acceleration
type(fm), allocatable :: rho(:)     ! Normalized density
type(fm), allocatable :: mu(:)      ! Normalized rigidity
type(fm), allocatable :: eta(:)     ! Normalized viscosity
!
type(fm), allocatable :: mlayer(:)  ! Normalized mass of each layers
type(fm) :: mass                    ! Normalized planet mass
!
type(fm) :: r0                      ! Reference length scale
type(fm) :: rho0                    ! Reference density scale
type(fm) :: mu0                     ! Reference rigidity scale
type(fm) :: eta0                    ! Reference viscosity scale
type(fm) :: t0                      ! Reference time scale
type(fm) :: mass0                   ! Reference mass
!
integer, allocatable :: irheol(:)   ! Rheology
!
end module
!
!
!
program alma
use fmzm
use general_parameters
use rheological_profile
implicit none
!
!
 character(100) :: cfg_f              ! Configuration file from cmd-line
!
 character(50), parameter :: alma_log='alma.log'
!
 integer :: i,j,k
 integer :: it, ik, idx
 type(fm) :: s
 type(fm) :: hh, ll, kk
 integer  :: n
 integer  :: ndeg
 type(fm), allocatable :: h_love(:,:), l_love(:,:), k_love(:,:)
 type(fm) :: f
!
!
!
!
!
!
!
 write(*,*) ''
 write(*,*) ' ****'
 write(*,*) ' **** This is ALMA'
 write(*,*) ' **** (the plAnetary Love nuMbers cAlculator)'
 write(*,*) ' ****'
 write(*,*) ''
!
!
!---------------------------------------------- Open and read configuration file
 call getarg(1,cfg_f)
 open(90,file=trim(cfg_f),status='old')
!
 write(*,*) ' - Parsing configuration file: ', trim(cfg_f)
!
 call config
!
 close(90)
!
!
!---------------------------------------------- Initialize the multi-precision package
 write(*,*) ' - Initializing the multi-precision libraries'
 call FM_SET (nd)
!
 pi  = to_fm('2') * asin(to_fm('1'))
 Gnt = to_fm('6.674e-11')
!
!
!---------------------------------------------- Initialize the multi-precision package
 write(*,*) " - Opening the log file '"//trim(alma_log)//"'"
 open(99,file=trim(alma_log),status='unknown')
!
!
!---------------------------------------------- Build the rheological model
 write(*,*) ' - Building the rheological model'
!
 call build_model
 call write_log(1)
!
 call normalization
 call write_log(2)
!
!
!---------------------------------------------- Build the time steps
 write(*,*) ' - Building the time steps'
!
 call time_steps
!
!
!
!
!
!---------------------------------------------- Compute the LNs
!
 call salzer_weights
!
!
 ndeg = 0
 do n=lmin,lmax,lstep
    ndeg = ndeg+1
 end do
!
 allocate( h_love(ndeg,p+1) )
 allocate( l_love(ndeg,p+1) )
 allocate( k_love(ndeg,p+1) )
!
 h_love = to_fm('0.0')
 l_love = to_fm('0.0')
 k_love = to_fm('0.0')
!
 idx = 0
! 
 do n=lmin,lmax,lstep
!
 write(*,*) ' - Harmonic degree n = ',n 
!
 idx = idx+1
!
 do it=1,p+1
!
 f = log( to_fm('2') ) / t(it) 
! 
 do ik=1,2*order
!
 s = f * to_fm(ik)
!
 call love_numbers(n,s,hh,ll,kk)
!
 h_love(idx,it) = h_love(idx,it) + hh * zeta(ik) * f
 l_love(idx,it) = l_love(idx,it) + ll * zeta(ik) * f
 k_love(idx,it) = k_love(idx,it) + kk * zeta(ik) * f 
! 
 end do
! 
 end do
!
 end do
!
 idx=0
! 
 do n=lmin,lmax,lstep
!
 idx = idx + 1
!
 do it=1,p+1
    print *,n,it, to_dp(h_love(idx,it)), &
	              to_dp(l_love(idx,it)), &
				  to_dp(k_love(idx,it))
 end do
! 
 end do
! 
 

!
!
! n = 20
! s = to_fm('0.5')
! call love_numbers(n,s,hh,ll,kk)
! print *, to_dp(hh)
! print *, to_dp(ll)
! print *, to_dp(kk)
! 
!
!
! j=5
! call surface_bc(5,r(nla+2),gra(nla+2),b)
! do i=1,3
!     print *,to_dp(b(i))
!end do
!
! j=5
! call direct_matrix(j,r(2),rho(2),mu(2),gra(2),YYd)
! call inverse_matrix(j,r(2),rho(2),mu(2),gra(2),YYi)
! YY = matmul(YYd,YYi)
! do i=1,6
!    print *, (to_dp(YY(i,j)),j=1,6)
! end do
!
!
!
end
