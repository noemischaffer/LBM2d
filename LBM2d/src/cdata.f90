module Cdata
  use Sub
  implicit none
  public
!***********************!
  real,parameter :: d_pi = 3.14159265358979323846264338327950288419716939937510D+00
  integer,parameter :: labellen=25
  integer,parameter :: qmom = 9 ! we are solving the D2Q9 model
  double precision, dimension(qmom) :: weight
  integer, dimension(qmom) :: mirrorq 
  double precision, parameter :: wrest=4.0d0/9.0d0, wudlr=1.0d0/9.0d0,wcorner=1.0d0/36.0d0
  integer :: Nx=16
  integer :: Ny=16
  double precision :: Lx,Ly
  double precision, dimension(2,qmom) :: ee ! unit vectors of the lattice
  integer, dimension(2,qmom) :: ee_int ! integer unit vectors of the lattice  
  double precision,dimension(2) :: xhat=[1.0d0,0.0d0],yhat=[0.0d0,1.0d0]
  double precision :: tau=1.0d0
  double precision :: xzero=0.0d0,yzero=0.0d0
  double precision :: dx=1.0d0,dy=1.0d0
  double precision,allocatable, dimension(:,:,:) :: ff,fftemp,ffEq
  integer, allocatable, dimension(:,:) :: is_solid
  logical::lffaloc=.false.
  logical :: lstart=.true.
  logical :: lstream=.true.
  integer :: iTMAX=0
  integer, allocatable, dimension(:,:) :: surface
  double precision, allocatable, dimension(:,:,:,:) :: sigma
  double precision, allocatable, dimension(:,:) :: curl_uu
  double precision, allocatable, dimension(:) ::xx,yy
  integer :: Nts=0
  integer :: ndiag=1
  double precision, allocatable,dimension(:) :: ts_data
  character(len=labellen), allocatable,dimension(:) :: ts_name
!
  integer :: Nsurf
!
namelist /cdata_pars/ &
     Nx,Ny,tau,iTMAX,lstart,ndiag,lstream
!
contains
!***********************!
subroutine rparam_cdata(unit,iostat)
! read the namelist grid_pars
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=cdata_pars, IOSTAT=iostat)
!----------------------
endsubroutine rparam_cdata
!***********************!
subroutine allocate_cdata()
integer :: i,j
!---------------------------------------
! The array is_solid is occupied in the following
! manner. If it is  1  it is solid point
!                   0  it is a boundary point
!                   -1 it is a fluid point
!---------------------------------------
allocate(is_solid(Nx+2,Ny+2))
!
! first assume that all points are fluid
! this array is overwritten when the obstacles are set
is_solid=-1
allocate(ff(Nx+2,Ny+2,qmom))
ff=0.0d0
ff(:,:,5) = 1.
!
allocate(fftemp(Nx+2,Ny+2,qmom))
fftemp=ff
allocate(ffEq(Nx+2,Ny+2,qmom))
ffEq=ff
Lx=dx*dfloat(Nx)
Ly=dy*dfloat(Ny)
allocate(xx(Nx+2))
do i=1,Nx+1
  xx(i)=xzero+(i-2)*dx
enddo
allocate(yy(Ny))
do j=1, Ny
  yy(j)=yzero+(j-2)*dy
enddo
lffaloc=.true.
endsubroutine allocate_cdata
!***************************************************************  
subroutine free_cdata()
  if (lffaloc .eqv. .true.) then
     deallocate(ff)
     deallocate(fftemp)
     deallocate(ffEq)
  endif
endsubroutine free_cdata
!***********************!
subroutine set_model()
  integer :: q,i,j
!-------------------------------------!

!--------------------------------------------!
  q=1
  do j=-1,1; do i=-1,1
     ee_int(1,q)=i
     ee_int(2,q)=j
     ee(1,q) = dfloat(ee_int(1,q)) 
     ee(2,q) = dfloat(ee_int(2,q))
     q=q+1
  enddo; enddo
  weight(1)=wcorner; weight(3)=wcorner; weight(7)=wcorner; weight(9)=wcorner
  weight(2)=wudlr ;  weight(4)=wudlr;   weight(6)=wudlr; weight(8)=wudlr
  weight(5)=wrest
  mirrorq(1)=9;mirrorq(2)=8;mirrorq(3)=7
  mirrorq(4)=6;mirrorq(5)=5;mirrorq(6)=4
  mirrorq(7)=3;mirrorq(8)=2;mirrorq(9)=1
!---------------------------------------------!
endsubroutine set_model
!***********************!
endmodule Cdata
