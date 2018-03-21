!***************************************************************
module Diagnostic
use Cdata
use Avg
use messages
implicit none
private
public :: initialize_diag,calc_diag, read_param_diag
integer :: iuxmax=0,iuymax=0,iurms=0
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
  jdiag=0
  iuxmax=jdiag+1; ts_name(iuxmax) = 'uxmax'; jdiag=jdiag+1 
  iuymax=jdiag+1; ts_name(iuymax) = 'uymax'; jdiag=jdiag+1
  iurms=jdiag+1;  ts_name(iurms) = 'urms'  ; jdiag=jdiag+1
endsubroutine initialize_diag
!***************************************************************
subroutine calc_diag()
  if (iuxmax.eq.1) ts_data(iuxmax) = maxval(uu(:,:,1))
  if (iuymax.eq.2) ts_data(iuymax) = maxval(uu(:,:,2))
  if (iurms.eq.3)  then 
    ts_data(iurms) = sqrt(sum(uu(:,:,1)*uu(:,:,1)+uu(:,:,2)*uu(:,:,2)))
    ts_data(iurms) = ts_data(iurms)/dfloat(Nx*Ny)
  endif
endsubroutine calc_diag
!***************************************************************
endmodule Diagnostic
