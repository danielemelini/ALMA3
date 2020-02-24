! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'salzer_weights'
! computes the weights zeta(i) to be used in the Salzer accelerated 
! Post-Widder Laplace inversions scheme
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine salzer_weights
use fmzm
use general_parameters
implicit none
!
 integer :: m
 integer :: j1, j2, j ,k
 type(fm) :: q1, q2, q3
 type(fm) :: fattm
!
!
 m = order
 allocate( zeta(2*m) )
! 
 do k=1,2*m
   j1 = floor( (k+1.0)/2.0 )
   j2 = min( k, m )
   zeta(k) = to_fm( '0.0' )
   do j=j1,j2
     call fmfact( to_fm(m), fattm )
     call fmcomb( to_fm(m), to_fm(j), q1 )
     call fmcomb( to_fm(2*j), to_fm(j), q2 )
     call fmcomb( to_fm(j), to_fm(k-j), q3 )
     zeta(k) = zeta(k) + ( to_fm(j) ** (m+1) ) / fattm * q1 * q2 * q3
   end do
   if (mod(m+k,2) .ne. 0) zeta(k) = -zeta(k)
 end do
!
end subroutine salzer_weights