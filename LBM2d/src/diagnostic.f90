!***************************************************************
module Diagnostic
use Cdata
use Avg
use messages
implicit none
private
public :: initialize_diag,calc_diag, read_param_diag
integer :: iuxmax=0,iuymax=0,iusqr=0
namelist /diag/& 
     Nts
!***************************************************************
contains
!*********************************!
subroutine read_param_diag(unit,iostat)
  integer :: jdiag
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  ! At this moment this namelist is not being read
  ! but hardcoded. 
!  read(unit, NML=diag, IOSTAT=iostat)

endsubroutine read_param_diag
!***************************************************************
subroutine initialize_diag()
  integer :: jdiag=0
  Nts=3
  allocate(ts_name(Nts))
  allocate(ts_data(Nts))
  jdiag=1
  iuxmax=jdiag; ts_name(iuxmax) = 'uxmax'; jdiag=jdiag+1 
  iuymax=jdiag; ts_name(iuymax) = 'uymax'; jdiag=jdiag+1
  iusqr=jdiag;  ts_name(iusqr) = 'usqr'  ; jdiag=jdiag+1
!  write(*,*) 'iuxmax',iuxmax,ts_name(iuxmax)
 ! write(*,*) 'iuymax',iuymax,ts_name(iuymax)
  !write(*,*) 'iusqr',iusqr,ts_name(iusqr)
endsubroutine initialize_diag
!***************************************************************
subroutine calc_diag()
  integer :: ix
  if (iuxmax.ne.0) ts_data(iuxmax) = maxval(uu(:,:,1))
  if (iuymax.ne.0) ts_data(iuymax) = maxval(uu(:,:,2))
  if (iusqr.ne.0)  then 
     ts_data(iusqr) = sum(uu(2:Nx+1,2:Ny+1,1)*uu(2:Nx+1,2:Ny+1,1)+&
          uu(2:Nx+1,2:Ny+1,2)*uu(2:Nx+1,2:Ny+1,2))
     ts_data(iusqr) = ts_data(iusqr)/dfloat(Nx*Ny)
  endif
!  write(*,*) 'DM',uu(2:5,3,1)
!  write(*,*) 'DM',uu(2:5,2,1)
 
endsubroutine calc_diag
!***************************************************************
endmodule Diagnostic
