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

module Evolve
  use Cdata
  use Sub
  implicit none
  private
  public :: stream,comp_equilibrium_BGK,collision,get_feq 
  integer,allocatable,dimension(:) :: iin_pencil
!------the following are public ------------
!-------------------------------------------
contains
!***************************************************************
subroutine initialize_evolve()

endsubroutine initialize_evolve
!***************************************************************
subroutine finalize_evolve()

endsubroutine finalize_evolve
!***************************************************************
subroutine stream()
  integer :: m,n,q,ix,iy,p
  fftemp=ff
  do q=1,qmom
     do iy=2,Ny+1
        n=iy-ee_int(2,q)
        do ix=2,Nx+1
           m=ix-ee_int(1,q)
           select case(is_solid(ix,iy))
           case(1) ! do nothing
           case(0)
              p=mirrorq(q)
              ff(ix,iy,p) = fftemp(ix,iy,q)
              ! stream and bounce-back at immersed boundary
           case(-1)
              ff(ix,iy,q) = fftemp(m,n,q) ! stream
             ! write(*,*) fftemp(i,j,q), i, j, q, 'streamed from', iin, jin, q
           endselect
        enddo
     enddo
  enddo
endsubroutine stream
!***************************************************************
subroutine comp_equilibrium_BGK()
  use Avg
  use Force
  integer :: q,ix,iy
  double precision,dimension(2) :: uf
  double precision :: edotu,usqr
  do q=1,qmom
     do iy=2,Ny+1
        do ix=2,Nx+1
           if(is_solid(ix,iy).eq.-1) then
              call get_uforced(ix,iy,uf)
              call get_feq(q,uf,rho(ix,iy),ffEq(ix,iy,q))
           endif
        enddo
     enddo
  enddo
endsubroutine comp_equilibrium_BGK
!***************************************************************
subroutine get_feq(q,uin,rhoin,fEq)
  integer,intent(in) :: q
  double precision,dimension(2),intent(in) :: uin
  double precision,intent(in) :: rhoin
  double precision,intent(out) :: fEq
  double precision :: edotu,usqr,esqr
  edotu=dot2d(ee(:,q),uin)
  usqr=dot2d(uin,uin)
  fEq = weight(q)*rhoin*(1.+3.*edotu &
                   +(9./2.)*(edotu**2) &
                   -(3./2.)*usqr & ! is there is misprint here in the book ?
                   )
endsubroutine get_feq
!***************************************************************
subroutine collision()
  integer :: q,ix,iy,p
!
! If a point is solid(1)   : do nothing
!               surface(0) : bounce-back
!               fluid(-1)  : collision  
! 
  do q=1,qmom
     do iy=2,Ny+1
        do ix=2,Nx+1
           if (is_solid(ix,iy).eq.-1) then
              ff(ix,iy,q) = ff(ix,iy,q) + (ffEq(ix,iy,q)-ff(ix,iy,q))/tau
!              write(*,*) ix,iy,q,fftemp(ix,iy,q),ffEq(ix,iy,q),tau
           endif
        enddo
     enddo
  enddo
!
endsubroutine collision
!***************************************************************
endmodule Evolve
