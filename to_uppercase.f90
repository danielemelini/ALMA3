function to_uppercase(string) result(upper)
character(len=*), intent(in) :: string
character(len=len(string)) :: upper
integer :: j
do j = 1,len(string)
  if(string(j:j) >= "a" .and. string(j:j) <= "z") then
       upper(j:j) = achar(iachar(string(j:j)) - 32)
  else
       upper(j:j) = string(j:j)
  end if
end do
end function to_uppercase
