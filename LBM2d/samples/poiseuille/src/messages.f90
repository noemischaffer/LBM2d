! Contains messages
!**************************************************************
module Messages
  use Cdata
  implicit none
  private
  character(LEN=2*labellen) :: scaller=''
  public :: fatal_error
contains
!**************************************************************
subroutine fatal_error(location,message,force)
  ! copied from pencil-code
  character(len=*), optional :: location
  character(len=*)           :: message
  logical,          optional :: force
!
  if (present(location)) scaller=location
  if (scaller=='') then
     write (*,*) trim(message)
  else
     write (*,*) trim(scaller) // ": " // trim(message)
  endif
  STOP
endsubroutine fatal_error
!***************************************************************
endmodule Messages
