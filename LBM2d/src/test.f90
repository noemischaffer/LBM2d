program test

implicit none

integer, dimension(5,5) :: a
integer :: i, j,d
d = 0

do i=1,5
do j=1,5 
  a(i,j)=d
d = d+1
enddo
enddo
open(unit=12, file='test.dat', action='write', status='replace')
write(12,*) a
close(12)

write(*,*) a



end
