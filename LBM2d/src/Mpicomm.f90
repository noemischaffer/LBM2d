module Mpicomm
  use Sub
  implicit none
  integer :: lroot=.true.
  public
!***********************************************************************
    subroutine die_gracefully
!
!  Stop having shutdown MPI neatly.
!  With at least some MPI implementations, this only stops if all
!  processors agree to call die_gracefully.
!
! COPIED FROM pencil-code      
!  29-jun-05/tony: coded
!
!  Tell the world something went wrong -- mpirun may not propagate
!  an error status.
!
      call touch_file('ERROR')
!
      call mpifinalize
      if (lroot) then
        STOP 1                    ! Return nonzero exit status
      else
        STOP
      endif
!
    endsubroutine die_gracefully  
!****************************************************************************
    subroutine touch_file(file)
!
!  Touches a given file (used for code locking).
!       
!! COPIED FROM pencil-code      
!  25-may-03/axel: coded
!  24-mar-10/Bourdin.KIS: moved here from sub.f90
!
      character(len=*) :: file
!
      integer :: unit = 1
!
      open (unit, FILE=file)
      close (unit)
!
    endsubroutine touch_file
!***********************************************************************    
endmodule Mpicomm
