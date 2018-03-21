module IO

use Cdata
use Avg
use messages 
implicit none
private
public:: read_input,write_input,write_model_param,write_snap
public :: write_immersed_boundary, write_ts
character (len=labellen) :: datadir='data'
character(len=labellen) :: input_param_file='input.param'
character(len=labellen) :: model_param_file='model.param'
character(len=labellen) :: usnap_file='usnap'
character(len=labellen) :: rsnap_file='rho_snap'
character(len=labellen) :: comment_char='#'
character(len=labellen) :: is_solid_file='solid.now'
character(len=labellen) :: ts_file='ts'
namelist /io_param/& 
     datadir,input_param_file,model_param_file, &
     usnap_file,rsnap_file, comment_char, is_solid_file
contains
!*********************************!
subroutine read_input()
use InitCond
use BoundCond
use Force
integer :: unit,iostat=0
open(unit=10,file='input',status='old')
call rparam_cdata(10,iostat)
call read_param_init(10,iostat)
call read_bc(10,iostat)
call read_fpars(10,iostat)
!call read_param_diag(10,iostat)
call read_param_io(10,iostat)
close(10)
endsubroutine read_input
!*********************************!
subroutine read_param_io(unit,iostat)
!read the namelist io
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=io_param, IOSTAT=iostat)
endsubroutine read_param_io
!*********************************!
subroutine write_input()
  logical :: exists
  integer :: unit=1
!
! First check if datadir exists
  inquire(FILE=datadir,EXIST=exists)
  if (exists) then
     write(*,*) 'simulation with:'
     write(*,*) 'Nx=',Nx,'Ny=',Ny,'Lx=',Lx,'Ly=',Ly
     write(*,*) 'unit of velocity, vunit=',vunit
     write(*,*) 'collision time scale tau=',tau
     write(*,*) 'Number of time steps, iTMAX=', iTMAX
     write(*,*) 'Writing parameters to file:',trim(datadir),'/',trim(input_param_file)
     open(unit,FILE=trim(datadir)//'/'//trim(input_param_file),FORM='formatted')
     write(unit,*) Nx,Ny,Lx,Ly,vunit,tau,iTMAX
     close(unit)
     write(*,*) '.......done.     '
  else
! if not abort
     write(*,*) 'datadir=',datadir,'not found!'
     call fatal_error('write_input','ABORTING')
  endif
  endsubroutine write_input
!*********************************!
subroutine write_model_param()
  integer :: q,i,j
  integer :: unit=1
  character (len=labellen) :: model_file='model.param'
  write(*,*) 'Writing model parameters to file:',trim(datadir),'/',trim(model_param_file)
  open(unit,FILE=trim(datadir)//'/'//trim(model_param_file),FORM='formatted')
  q=1
  write(unit,*) comment_char,'----ee vectors ------'
  do j=-1,1; do i=-1,1
     write(unit,*) q,ee_int(1,q),ee_int(2,q)
     write(unit,*) q,ee(1,q),ee(1,q)
     q=q+1
  enddo; enddo
  write(unit,*) comment_char,'----weights, mirrorq --------'
  do q=1,qmom
     write(unit,*)q,weight(q),mirrorq(q)
  enddo
  close(unit)
!
endsubroutine write_model_param
!*********************************!
subroutine write_snap(itime)
  integer,intent(in) :: itime
  integer :: ix,iy
  integer :: unit=1
  write(*,*) 'writing snapshot to files:'
  write(*,*) trim(datadir)//'/'//trim(usnap_file)
  write(*,*) trim(datadir)//'/'//trim(rsnap_file)
  open(unit,FILE=trim(datadir)//'/'//trim(usnap_file),FORM='formatted')
  do iy=1,Ny
     write(unit,*)(uu(ix,iy,1),uu(ix,iy,2),ix=1,Nx)
  enddo
  close(unit)
  open(unit,FILE=trim(datadir)//'/'//trim(rsnap_file),FORM='formatted')
  do iy=1,Ny
     write(unit,*)(rho(ix,iy),ix=1,Nx)
  enddo
  close(unit)
  write(*,*) '....done.'
  write(*,*) '================================='
endsubroutine write_snap
!*********************************!
subroutine write_immersed_boundary(itime)
  integer,intent(in) :: itime
  integer :: ix,iy
  integer :: unit=1
  write(*,*) 'writing is_solid array to file:'
  write(*,*) trim(datadir)//'/'//trim(is_solid_file)
  open(unit,FILE=trim(datadir)//'/'//trim(is_solid_file),FORM='formatted')
  do iy=1,Ny
     write(unit,*)(is_solid(ix,iy),ix=1,Nx)
  enddo
  close(unit)
endsubroutine write_immersed_boundary
!*********************************!
subroutine write_ts(itime)
  integer,intent(in) :: itime
  integer :: unit=3
  integer :: its
  if (itime.eq.0) then
     write(*,*) 'it     ',(ts_name(its),its=1,Nts)
     open(unit,FILE=trim(datadir)//'/'//trim(ts_file),FORM='formatted')
     write(unit,*) 'it',(ts_name(its),its=1,Nts)
  close(unit)
  endif
  write(*,*) itime,(ts_data(its),its=1,Nts)
  open(unit,FILE=trim(datadir)//'/'//trim(ts_file),FORM='formatted',POSITION='append')
  write(unit,*) itime,(ts_data(its),its=1,Nts)
  close(unit)
endsubroutine write_ts
!*********************************!
endmodule IO
