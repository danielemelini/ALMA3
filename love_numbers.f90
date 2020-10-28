! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'love_numbers'
! computes the Love numbers h,l,k in the Laplace domain for
! a given value of s
!
! Initial version DM February 24, 2020
! Modified DM June 11, 2020 - Burgers and Andrade rheologies
! Modified DM June 16, 2020 - Complex LNs
!                             (the 1/s factor is now outside this module)
! Fixed DM October 27, 2020 - Wrong call to complex_rigidty for core BC
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine love_numbers(n,s,hh,ll,kk)
use fmzm
use general_parameters
use rheological_profile
implicit none
!
integer  :: n
type(zm) :: s
type(zm) :: hh, ll, kk
!
integer :: i,j,k
!
type(zm) :: mu_s
type(zm) :: Ydir(6,6), Yinv(6,6)
type(zm) :: lambda(6,6)
type(zm) :: bc(6,3)
type(fm) :: bs(3)
type(zm) :: rr(3,3)
type(zm) :: qq(3,3)
type(zm) :: x(3)
type(zm) :: prod(6,3)
!
type(zm) :: aa(3,3)
type(zm) :: bb(3)
integer  :: d, code
integer  :: indx(3)
!
!
 lambda = to_fm('0.0')
 do i=1,6
     lambda(i,i) = to_fm('1.0')
 end do 
! 
!
! ---- Build the propagator product
!
 do j=(nla+1),1,-1
!
    call complex_rigidity(s,mu(j),eta(j),irheol(j),par(j,:),mu_s)
!
    call direct_matrix (n,r(j+1),rho(j),mu_s,gra(j+1),Ydir) 
    call inverse_matrix(n,r(j  ),rho(j),mu_s,gra(j  ),Yinv)
 !
    lambda=matmul(lambda,matmul(Ydir,Yinv))
!
 end do
!
!
! ---- Compute the boundary conditions
!
 if( irheol(0)==0 ) then
    call fluid_core_bc (n,r(1),rho(0),gra(1),bc)
 else
    call complex_rigidity(s,mu(0),eta(0),irheol(0),par(0,:),mu_s)
    call direct_matrix (n,r(1),rho(0),mu_s,gra(1),Ydir) 
    bc(:,1) = Ydir(:,1)
    bc(:,2) = Ydir(:,2)
    bc(:,3) = Ydir(:,3)    
 end if
! 
 call surface_bc    (n,r(nla+2),gra(nla+2),bs)
!
!
! ---- Compute the 'R' and 'Q' arrays
!
 prod = matmul(lambda,bc)
!
 rr(1,:) = prod(3,:)
 rr(2,:) = prod(4,:)
 rr(3,:) = prod(6,:) 
!
 qq(1,:) = prod(1,:)
 qq(2,:) = prod(2,:)
 qq(3,:) = prod(5,:) 
!
! 
! ---- Compute the solution vector 
!
 aa = rr
 bb = bs 
 call ludcmp(aa,3,3,indx,d,code)
 call lubksb(aa,3,3,indx,bb)
 x = matmul(qq,bb)
!
 hh = x(1) * mass / r(nla+2) 
 ll = x(2) * mass / r(nla+2) 
 kk = ( -to_fm('1') - x(3) * mass / ( r(nla+2) * gra(nla+2) ) ) 
! 
! 
! 
end subroutine love_numbers
