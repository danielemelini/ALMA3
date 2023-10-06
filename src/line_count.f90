! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'line_count'
! counts the number of lines and the number of header lines 
! (i.e. those starting with '#' or '!') in a data file
!
! DM October 5, 2023
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
 subroutine line_count( iu, nh, nl )
 implicit none
!
 integer :: nh, nl
 integer :: iu
 character :: c
 logical :: header
!
!
!
 nh = 0
 nl = 0
 header = .true.
!
 rewind(iu)
!
 100 read( iu, *, end=200 ) c

     if (header) then
        if( ( c.eq.'#' ) .or. ( c.eq.'!' ) ) then
           nh=nh+1
        else
           header=.false.
        end if
     end if
        
     nl = nl+1

     go to 100

 200 continue
!
 rewind(iu)
!
!
end subroutine line_count
