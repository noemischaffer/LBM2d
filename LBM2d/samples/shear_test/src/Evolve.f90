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
  public :: stream,comp_equilibrium_BGK,collision 
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
  integer :: k,l,q,m,n
  do q=1,qmom
     do l=2,Ny+1
        n=l-ee_int(2,q)
        do k=2,Nx+1
           if(is_solid(k,l).ne.1) then
              m=k-ee_int(1,q)
              fftemp(k,l,q) = ff(m,n,q)
             ! write(*,*) fftemp(i,j,q), i, j, q, 'streamed from', iin, jin, q
           endif
        enddo
     enddo
  enddo
endsubroutine stream
!***************************************************************
subroutine comp_equilibrium_BGK()
  use Avg
  use Force
  integer :: q,k,l
  double precision,dimension(2) :: ueq
  double precision :: edotu,usqr
  do q=1,qmom
     do l=2,Ny+1
        do k=2,Nx+1
           if(is_solid(k,l).ne.1) then
              call get_ueq(uu(k-1,l-1,:),rho(k-1,l-1),ueq)
              edotu=dot2d(ee(:,q),ueq)
              usqr=dot2d(ueq,ueq)
              ffEq(k,l,q) = weight(q)*rho(k-1,l-1)*(1.+3.*edotu/(vunit**2) &
                                       +(9./2.)*(edotu**2)/(vunit**4) &
                                       -(3./2.)*usqr/(vunit**2) &
               )
               !write(*,*) ffEq(i,j,q), i,j,q
               !write(*,*) '*********ueq*******'
               !write(*,*) ueq
               !write(*,*) '****************'
           endif
        enddo
     enddo
  enddo
endsubroutine comp_equilibrium_BGK
!***************************************************************
subroutine collision()
  integer :: q,k,l
  do q=1,qmom
     do l=2,Ny+1
        do k=2,Nx+1
           if(is_solid(k,l).ne.1) then
              ff(k,l,q) = fftemp(k,l,q) + (ffEq(k,l,q)-fftemp(k,l,q))/tau
              !write(*,*) ff(i,j,q), i, j
           endif
        enddo
     enddo
  enddo



endsubroutine collision
!***************************************************************
endmodule Evolve
