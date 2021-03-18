! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'read_data_line'
! Reads from file opened at LUN the next line not starting with '!' or '#'
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
    subroutine read_data_line(lun,buffer)
    implicit none
    character(*)   :: buffer
    character :: ch
    integer :: lun
    integer :: idx
!
910 read(lun,'(a)') buffer
    ch = buffer(1:1)
    if ((ch.eq.'#').or.(ch.eq.'!')) go to 910
!
    if( (index( buffer, '#' ).ne.0) .or. (index( buffer,'!' ).ne.0) ) then
! 
       if( index(buffer,'#').eq.0 ) then
           idx = index(buffer,'!')
       elseif( index(buffer,'!').eq.0 ) then
           idx = index(buffer,'#')
       else
           idx = min( index( buffer, '#' ), index( buffer, '!' ) )
       endif	   
!
       buffer = buffer( 1:(idx-1) )
!
    end if
!
    return
!
    end
!
!    
