! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'config'
! reads the configuration file from unit 90
!
! Initial version DM February 24, 2020
! Modified by DM on June 16, 2020 - complex LNs
! Modified by DM on Aug 3, 2020 - config file now contains the
!                                 total number of layers
! Modified by DM on Sep 11, 2020 - external time steps
! Modified by DM on June 21, 2021 - linear time-history
! Modified by DM on Aug 11, 2021 - LN rate of change
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
subroutine config
use general_parameters
use rheological_profile
implicit none
!
 character(200) :: buffer             ! Buffer for cfg file parsing
 character(50) :: cjunk
 character(50), external :: to_uppercase
!
!
!---------------------------------------------- Parsing configuration file
 call read_data_line(90,buffer)   ;   read(buffer,*) nd
 call read_data_line(90,buffer)   ;   read(buffer,*) order
!
 call read_data_line(90,buffer)   ;   read(buffer,*) cjunk
 cjunk=to_uppercase(cjunk)
 if( trim(adjustl(cjunk)) == 'LOADING' ) then
    iload=1
 elseif( trim(adjustl(cjunk)) == 'TIDAL' ) then
    iload=0
 else
    write(*,*) " - ERROR: Unknown load type '"//trim(adjustl(cjunk))//"'"
	stop
 end if
!
 call read_data_line(90,buffer)   ;   read(buffer,*) lmin
 call read_data_line(90,buffer)   ;   read(buffer,*) lmax
 call read_data_line(90,buffer)   ;   read(buffer,*) lstep
!
 call read_data_line(90,buffer)   ;   read(buffer,*) cjunk
 cjunk=to_uppercase(cjunk)
 if( trim(adjustl(cjunk)) == 'EXT' ) then
    itime=2
 elseif( trim(adjustl(cjunk)) == 'LOG' ) then
    itime=1
 elseif( trim(adjustl(cjunk)) == 'LIN' ) then
    itime=0
 else
    write(*,*) " - ERROR: Unknown time scale '"//trim(adjustl(cjunk))//"'"
	stop
 end if
!
 call read_data_line(90,buffer)  ;   read(buffer,*) p
 call read_data_line(90,buffer)  ;   read(buffer,*) m1,m2
!
 call read_data_line(90,buffer)  ;   read(buffer,*) cjunk
 cjunk=to_uppercase(cjunk)
 if( trim(adjustl(cjunk)) == 'STEP' ) then
    ihist = 1
 elseif( trim(adjustl(cjunk)) == 'RAMP' ) then
    ihist = 2
 else 
    write(*,*) " - ERROR: Unknown load history '"//trim(adjustl(cjunk))//"'"
    stop
 end if
!
 call read_data_line(90,buffer)  ;   read(buffer,*) tau
!
 call read_data_line(90,buffer)  ;   read(buffer,*) nla 
 if (nla.le.0) then
    write(*,*) " - ERROR: The model must contain at least one layer"
	stop
 end if
!
 call read_data_line(90,buffer)  ;   file_rheol = adjustl(trim(buffer))
!
 call read_data_line(90,buffer)  ;   file_log   = adjustl(trim(buffer))
!
 call read_data_line(90,buffer)  ;   read(buffer,*) cjunk
 cjunk=to_uppercase(cjunk)
 if( trim(adjustl(cjunk)) == 'REAL' ) then
    itype=0
 elseif( trim(adjustl(cjunk)) == 'COMPLEX' ) then
    itype=1
 elseif( trim(adjustl(cjunk)) == 'RATE' ) then
    itype=2
 else
    write(*,*) " - ERROR: Unknown LN mode '"//trim(adjustl(cjunk))//"'"
    stop
 end if
!
 if( itype==2 ) then
    flag_rate=.true.
 else
    flag_rate=.false.
 end if   
!
 call read_data_line(90,buffer)  ;   read(buffer,*) cjunk
 cjunk=to_uppercase(cjunk)
 if( trim(adjustl(cjunk)) == 'LN_VS_N' ) then
    ifmt=1
 elseif( trim(adjustl(cjunk)) == 'LN_VS_T' ) then
    ifmt=2
 else
    write(*,*) " - ERROR: Unknown output format '"//trim(adjustl(cjunk))//"'"
	stop
 end if
!
!
 call read_data_line(90,buffer)  ;   file_h     = adjustl(trim(buffer))
 call read_data_line(90,buffer)  ;   file_l     = adjustl(trim(buffer))
 call read_data_line(90,buffer)  ;   file_k     = adjustl(trim(buffer))
!
!
!
!
end subroutine config
