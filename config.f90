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
!---------------------------------------------- Open log file
 !open(99,file=alma_log,status='unknown')
 !call date_and_time(values=date_vec)
 !write(99,*)
 !write(99,'(a26,2(i2.2,a1),i4,1x,i2.2,2(a1,i2.2))') &
 !           '  ALMA log file opened on ', &
 !           date_vec(3),'-',date_vec(2),'-',date_vec(1), &
 !			date_vec(5),':',date_vec(5),':',date_vec(7)
 ! write(99,*)
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
 if( trim(adjustl(cjunk)) == 'LOG' ) then
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
 call read_data_line(90,buffer)  ;   read(buffer,*) nla
 call read_data_line(90,buffer)  ;   file_rheol = adjustl(trim(buffer))
!
!
!---------------------------------------------- Write cfg info to logfile
 !write(99,*) ' Numerical precision is ',nd
 !write(99,*) ' Order of the Gaver sequence is ',order
!
!
!
!
!
!
!
end subroutine config