!
! ALMA
! the plAnetary Love nuMbers cAlculator
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2020 Daniele Melini
!
! Based on a re-writing of ALMA 2.2
! Initial version February 24, 2020
! Modified June 11, 2020 for Burgers rheology
! Modified June 16, 2020 for complex LNs
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
! +-------------------------------+
! |   MODULE GENERAL_PARAMETERS   |
! +-------------------------------+
!
! Stores a number of parameters used throughout alma routines
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
type(fm), allocatable :: t(:)       ! Time steps (or period)
!
real(4) :: m1,m2                    ! Time range is 10^(m1:m2)
integer :: p                        ! Number of times minus one
!
type(fm), allocatable :: zeta(:)    ! Salzer weights
!
character(100) :: file_rheol        ! File with rheological model
character(100) :: file_h            ! File for 'h' LN
character(100) :: file_l            ! File for 'l' LN
character(100) :: file_k            ! File for 'k' LN
character(100) :: file_log          ! Log file
!
integer :: ifmt                     ! Output file format ('1'->LN vs n, '2'->LN vs t)
integer :: itype                    ! Mode switch ('0'-> Real LNs, '1'->Complex LNs)
!
end module
!
!
!
! +--------------------------------+
! |   MODULE RHEOLOGICAL_PROFILE   |
! +--------------------------------+
!
! Stores a information about the rheological model
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
type(fm), allocatable :: par(:,:)   ! Extra rheology parameters
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
! ##################################################################
! #                                                                #
! #                          MAIN PROGRAM                          #
! #                                                                #
! ##################################################################
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
 integer :: i,j,k
 integer :: it, ik, idx
 type(zm) :: s
 type(zm) :: hh, ll, kk
 integer  :: n
 integer  :: ndeg
 type(zm), allocatable :: h_love(:,:), l_love(:,:), k_love(:,:)
 type(fm) :: f
 type(fm) :: omega
 type(zm) :: iota
!
 real(4) :: time_0, time_1, time_2 
 logical :: file_exists
 integer :: iu
 character, parameter :: ch(3)= (/ 'h', 'l', 'k' /)
 character(50) :: timestamp
 integer(8) :: idate(8)
!
!
!
! ########################## Execution starts here ##########################
!
!
 call print_logo
!
 write(*,*) ''
 write(*,*) ' ****'
 write(*,*) ' **** This is ALMA'
 write(*,*) ' **** (the plAnetary Love nuMbers cAlculator)'
 write(*,*) ' ****'
 write(*,*) ''
!
!
!---------------------------------------------- Check command line arguments
!
 if (iargc().ne.1) then
    write(*,*) ' - USAGE: alma.exe <configuration_file> '
    write(*,*) ''
    stop
 end if
!
!---------------------------------------------- Obtain initial time-mark
 call cpu_time(time_0)
!
!---------------------------------------------- Open and read configuration file
 call getarg(1,cfg_f)
!
 inquire(file=trim(cfg_f), exist=file_exists)
 if( .not. file_exists ) then
     write(*,*) " - ERROR: configuration file '"//trim(cfg_f)//"' does not exist."
     write(*,*) ""
     stop
 end if
!
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
 pi   = to_fm('2') * asin(to_fm('1'))
 iota = to_zm('(0,1)')
 Gnt  = to_fm('6.674e-11')
!
!
!---------------------------------------------- Open the log file
 write(*,*) " - Opening the log file '"//trim(file_log)//"'"
 open(99,file=trim(file_log),status='unknown')
 call write_log(0)
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
 call write_log(3)
!
!
!---------------------------------------------- Compute the LNs
!
 call salzer_weights
!
 call cpu_time(time_1)
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
!============================================== Real LNs
 if (itype==0 ) then        
!
    do n=lmin,lmax,lstep
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
            h_love(idx,it) = h_love(idx,it) + (hh / s) * zeta(ik) * f 
            l_love(idx,it) = l_love(idx,it) + (ll / s) * zeta(ik) * f 
            k_love(idx,it) = k_love(idx,it) + (kk / s) * zeta(ik) * f 
! 
          end do
! 
       end do
!
       call cpu_time(time_2)
!
       write(*,*) ' - Harmonic degree n = ',n,'(',time_2-time_1,'s)' 
       time_1 = time_2
!
    end do
!
!============================================== Complex LNs
 elseif (itype==1) then
!  
    do n=lmin,lmax,lstep
!
       idx = idx+1
!
       do it=1,p+1
!
          omega = 2 * pi / t(it)
          s     = iota * omega
! 
          call love_numbers(n,s,hh,ll,kk)
!
          h_love(idx,it) = hh
          l_love(idx,it) = ll
		  k_love(idx,it) = kk
! 
       end do
!
       call cpu_time(time_2)
!
       write(*,*) ' - Harmonic degree n = ',n,'(',time_2-time_1,'s)' 
       time_1 = time_2	   
!
    end do
!     
 end if
!
!
!
!---------------------------------------------- Write outputs
!
 write(*,*) " - Writing output files '"//trim(file_h)//"', '"// & 
           trim(file_l)//"', '"//trim(file_k)//"'"
!
 call date_and_time(values=idate)
 write(timestamp,'(i4,2(a1,i2.2),1x,i2.2,2(a1,i2.2))') &
    idate(1),'-',idate(2),'-',idate(3), &
    idate(5),':',idate(6),':',idate(7)
!
 open(71,file=trim(file_h),status='unknown')
 open(72,file=trim(file_l),status='unknown')
 open(73,file=trim(file_k),status='unknown')
!
 do iu=71,73 
    write(iu,'(a)') '# ----------------------------------------------------------' 
    if(iload==0) write(iu,'(a)') '# '//ch(iu-70)//' tidal Love number '
    if(iload==1) write(iu,'(a)') '# '//ch(iu-70)//' load Love number '
    write(iu,'(a)') '# Created by ALMA on '//trim(timestamp)
    write(iu,'(a)') '# ----------------------------------------------------------' 
    write(iu,'(a)') '# '  
 end do
!
!
 if( ifmt==1 ) then
!
   idx=0
! 
   do n=lmin,lmax,lstep
!
      idx = idx + 1
!
      if (itype==0) then
         write(71,'(i4,1x,3048(e19.8))') n,(to_dp(h_love(idx,it)),it=1,p+1) 
         write(72,'(i4,1x,3048(e19.8))') n,(to_dp(l_love(idx,it)),it=1,p+1) 
         write(73,'(i4,1x,3048(e19.8))') n,(to_dp(k_love(idx,it)),it=1,p+1) 
      elseif (itype==1) then
         write(71,'(i4,1x,3048(e19.8))') n,(to_dpz(h_love(idx,it)),it=1,p+1)
         write(72,'(i4,1x,3048(e19.8))') n,(to_dpz(l_love(idx,it)),it=1,p+1)
         write(73,'(i4,1x,3048(e19.8))') n,(to_dpz(k_love(idx,it)),it=1,p+1)
	  end if
!
   end do
!
 elseif( ifmt==2 ) then
!
   do it=1,p+1
!  
      if (itype==0) then
         write(71,'(3048(e19.8))') to_dp(t(it)),(to_dp(h_love(idx,it)),idx=1,ndeg) 
         write(72,'(3048(e19.8))') to_dp(t(it)),(to_dp(l_love(idx,it)),idx=1,ndeg) 
         write(73,'(3048(e19.8))') to_dp(t(it)),(to_dp(k_love(idx,it)),idx=1,ndeg) 
      elseif (itype==1) then
         write(71,'(3048(e19.8))') to_dp(t(it)),(to_dpz(h_love(idx,it)),idx=1,ndeg)
         write(72,'(3048(e19.8))') to_dp(t(it)),(to_dpz(l_love(idx,it)),idx=1,ndeg) 
         write(73,'(3048(e19.8))') to_dp(t(it)),(to_dpz(k_love(idx,it)),idx=1,ndeg)	   
      end if
!
   end do
!
 end if
! 
 close(71)
 close(72)
 close(73) 
!
!
!---------------------------------------------- Close the log file
!
 write(*,*) " - Closing the log file '"//trim(file_log)//"'"
 close(99)
! 
!
!
!---------------------------------------------- Release dynamic arrays
!
 deallocate(t)
 deallocate(h_love)
 deallocate(l_love)
 deallocate(k_love)
 deallocate(r)
 deallocate(rho)
 deallocate(mu)
 deallocate(eta)
 deallocate(irheol)
! 
!
!
!---------------------------------------------- All done
!
 call cpu_time(time_1)
 write(*,*) ' - ALMA job completed. Time elapsed: ',time_1-time_0,'s'
 write(*,*) ''
!
!
!
end
!