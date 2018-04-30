! LBM2D
!
!  calculates evolution of ff
!
!The distribution function ff can have the following structures:
!  Array of Structures (AoS) ff(9,Nx+2,Ny+2)
!  Structure of Arrays (SoA) ff(Nx+2,Ny+2,9)
!  Bundle                    ff(3,Nx+2,Ny+2,3)
! Accordingly to Mattila et al 2007 the last one is the fastest
! in terms of cache misses ! But as the SoA is closest to the pencil
! structures we use that now. 
!***************************************************************
!
!The convention is that:
!  if a loop goes between 1-Nx+2 (1-Ny+2) the loop index is i (j)
!  if a loop goes between 2-Nx+1 (2-Ny+1) the loop index is k (l)
!
!***************************************************************
module Force
  use Messages
  use Cdata
  use Sub
  use Avg
  implicit none
  private
  public :: get_uforce, read_fpars
  character (len=labellen) :: iforce='none'
  double precision :: famp=0.0d0,gm=0.0d0,kk=0.0d0
  namelist /force_pars/ &
     iforce,famp,gm,kk
!***************************************************************
contains
!***************************************************************
subroutine read_fpars(unit,iostat)
! read the namelist force_pars
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=force_pars, IOSTAT=iostat)
endsubroutine read_fpars
  !***************************************************************
subroutine get_uforce(ii,jj,uout)
  integer, intent(in) :: ii,jj
  double precision,dimension(2), intent(out) :: uout
  double precision, dimension(2) :: gg,force
  double precision :: rhoin
  double precision, dimension(2) :: uin

  uin=uu(ii,jj,:)
  rhoin=rho(ii,jj)
  select case (iforce)
  case('none')
     uout=uin
  case('xgravity')
    gg(1)=gm
    gg(2)=0.0d0
    force=gg*rhoin
  case('kolmogorov_x')
    force(1)=famp*sin(2.0d0*d_pi*kk*yy(jj)/Ly)
    force(2)=0.0d0
  case('kolmogorov_y')
    force(1)=0.0d0
    force(2)=famp*sin(2.0d0*d_pi*kk*xx(ii)/Lx)
  case default
     call fatal_error("force","iforce does not match")
  endselect
  uout=uin+(tau*force/rhoin)
endsubroutine get_uforce
!***************************************************************
endmodule Force
