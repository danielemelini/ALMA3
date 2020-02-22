 subroutine time_steps
! 
! Generates p+1 time points between t=10^{M1} and t=10^{M2} kyrs according to 
! linearly or logarithmically spaced scales. Love numbers will be computed at
! these time steps. 
!
 use fmzm
 use general_parameters
 implicit none
!
 integer :: i
!
 type(fm) :: t1, t2, dt
 type(fm) :: a, a1, a2, da
!
!
! 
 allocate( t(p+1) )
!
 if (itime==0) then        ! linear
!
    t1=to_fm('10')**m1
    t2=to_fm('10')**m2
!
    if( p==0 ) then
       t(1) = t1
    else
       dt = (t2-t1) / to_fm(p)
       do i=1,(p+1)
          t(i) = t1 + dt * ( to_fm(i) - to_fm('1') )
       end do
    endif
!
 elseif (itime==1) then    ! log
!
        a1 = to_fm(m1)
        a2 = to_fm(m2)
!               
        if( p==0 ) then
           t(1) = to_fm('10')**a1
 	else  
           da = (a2-a1) / to_fm(p)
           do i=1,(p+1) 
              a = a1 + da * ( to_fm(i) - to_fm('1') )
              t(i) = to_fm('10')**a
           end do 
        end if
!
!
 end if
!
end subroutine time_steps
!
!
!