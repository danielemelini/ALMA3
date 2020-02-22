subroutine build_model
use rheological_profile
use general_parameters
implicit none
!
 integer :: i,j,k,n
 integer, parameter :: nh=4            ! Number of header lines in rheology file
 character(20) :: buffer(5)
 character(3)  :: code    
!
 character(20), external :: to_uppercase
!
!
! ----- Open the model file
 open(10,file=trim(file_rheol),status='old')
!
 call count_rows(10,n)
 if( n.ne.(nla+2+nh) ) then
    write(*,*) " - ERROR: Number of rows in '"//trim(file_rheol)//"' is not consistent with"
    write(*,*) "          the number of mantle layers."
    stop
 end if
! 
 rewind(10)
!
 do i=1,nh
    read(10,*)
 end do
! 
 allocate(r  (0:nla+2))
 allocate(gra(0:nla+2))
 allocate(rho(0:nla+1))
 allocate(mu (0:nla+1))
 allocate(eta(0:nla+1))
 allocate(irheol(0:nla+1))
!
 allocate(mlayer(0:nla+1))
!
 r(0)=to_fm('0.0')
!
 do i=nla+1,0,-1
 !
     read(10,*) (buffer(j), j=1,5)
 !
     r(i+1) = to_fm(buffer(1))
     rho(i) = to_fm(buffer(2))
     mu (i) = to_fm(buffer(3))
     eta(i) = to_fm(buffer(4))
!
     code=to_uppercase(buffer(5))
!
     if (code=='FLU') then
         irheol(i)=0
     elseif (code=='ELA') then
         irheol(i)=1
     elseif (code=='MAX') then
         irheol(i)=2
     else
         write(*,*) " - ERROR: Rheology '"//trim(buffer(5))//"' is unknown."
         stop
     end if
!
 end do
!
 close(10)
!
!
!
! -------- Mass of the planet
!
 mlayer(0) = rho(0) * r(1)**3
 do i=1,nla+1
    mlayer(i) = rho(i) * ( r(i+1)**3 - r(i)**3 )
 end do
 mlayer = mlayer * to_fm('4.0') / to_fm('3.0') * pi
!
 mass = sum(mlayer)
!
!
! ------- Gravity at the interface boundaries
!
 gra(0)=to_fm('0.0')
 do i=1, nla+2
     gra(i) = gnt * sum(mlayer(0:(i-1))) / r(i)**2
 end do
!
!
!
!
end subroutine build_model
!
!
!
!
    subroutine count_rows(lun,n)
    implicit none
!
! --- Returns as N the number of rows in text file opened at LUN
!     DM 28.04.2014
!
    integer lun,n
!
    rewind(lun)
!
    n=0
901 read(lun,*,end=902)
    n=n+1
    go to 901
    902 continue
!
    return
!
    end subroutine count_rows
