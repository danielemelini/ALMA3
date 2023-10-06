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
! Fixed DM June 18, 2021    - Added 'save' to local variables
!                             to avoid memory leaks (see FMLIB manual)
! DM November 29, 2022      - Degree 1 LNs
! Modified Oct  5, 2023     - changed the numbering scheme for the layers
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
type(zm), save :: mu_s
type(zm), save :: Ydir(6,6), Yinv(6,6)
type(zm), save :: lambda(6,6)
type(zm), save :: bc(6,3)
type(fm), save :: bs(3)
type(zm), save :: rr(3,3)
type(zm), save :: qq(3,3)
type(zm), save :: x(3)
type(zm), save :: prod(6,3)
!
type(zm), save :: aa(3,3)
type(zm), save :: bb(3)
!
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
 do j=nla,2,-1
!
    call complex_rigidity(s,mu(j),eta(j),irheol(j),par(j,:),mu_s)
!
    call direct_matrix (n,r(j  ),rho(j),mu_s,gra(j  ),Ydir) 
    call inverse_matrix(n,r(j-1),rho(j),mu_s,gra(j-1),Yinv)
 !
    lambda=matmul(lambda,matmul(Ydir,Yinv))
!
 end do
!
!
! ---- Compute the boundary conditions
!
 if( irheol(1)==0 ) then
    call fluid_core_bc (n,r(1),rho(1),gra(1),bc)
 else
    call complex_rigidity(s,mu(1),eta(1),irheol(1),par(1,:),mu_s)
    call direct_matrix (n,r(1),rho(1),mu_s,gra(1),Ydir) 
    bc(:,1) = Ydir(:,1)
    bc(:,2) = Ydir(:,2)
    bc(:,3) = Ydir(:,3)    
 end if
! 
 call surface_bc    (n,r(nla),gra(nla),bs)
!
!
! ---- Compute the 'R' and 'Q' arrays
!
 prod = matmul(lambda,bc)
!
 if (n==1) then 
    rr(1,:) = prod(3,:)
    rr(2,:) = prod(4,:)
    rr(3,:) = prod(5,:) 
 else
    rr(1,:) = prod(3,:)
    rr(2,:) = prod(4,:)
    rr(3,:) = prod(6,:) 
 end if
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
 hh = x(1) * mass / r(nla) 
 ll = x(2) * mass / r(nla) 
 kk = ( -to_fm('1') - x(3) * mass / ( r(nla) * gra(nla) ) ) 
! 
! 
! 
end subroutine love_numbers
