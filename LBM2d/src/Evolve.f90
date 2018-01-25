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
  integer :: i,j,q,iin,jin
  do q=1,qmom
     do j=2,Ny+1
        jin=j-ee_int(2,q)
        do i=2,Nx+1
           if(is_solid(i,j).ne.1) then
              iin=i-ee_int(1,q)
              fftemp(i,j,q) = ff(iin,jin,q)
              !write(*,*) fftemp(i,j,q), i, j, q
              !write(*,*) '***************************'
              !write(*,*) ff(iin, jin, q), iin, jin, q
              !write(*,*) '************new step*************'
           endif
        enddo
     enddo
  enddo
endsubroutine stream
!***************************************************************
subroutine comp_equilibrium_BGK()
  use Avg
  use Force
  integer :: q,i,j
  double precision,dimension(2) :: ueq
  double precision :: edotu,usqr
  do q=1,qmom
     do j=2,Ny+1
        do i=2,Nx+1
           if(is_solid(i,j).ne.1) then
              call get_ueq(uu(i-1,j-1,:),rho(i-1,j-1),ueq)
              edotu=dot2d(ee(:,q),ueq)
              usqr=dot2d(ueq,ueq)
              ffEq(i,j,q) = weight(q)*rho(i-1,j-1)*(1.+3.*edotu/(vunit**2) &
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
  integer :: q,i,j
  do q=1,qmom
     do j=2,Ny+1
        do i=2,Nx+1
           if(is_solid(i,j).ne.1) then
              ff(i,j,q) = fftemp(i,j,q) + (ffEq(i,j,q)-fftemp(i,j,q))/tau
           endif
        enddo
     enddo
  enddo
endsubroutine collision
!***************************************************************
endmodule Evolve
