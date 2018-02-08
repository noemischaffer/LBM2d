!LBM2D

!Calculates the wall shear stress according to Equation 12 in Jäger+2017
!**************************************************************************
module ShearStress

use Cdata
use Evolve
implicit none
private
public :: allocate_shear,  wall_shear, free_shear

double precision, allocatable, dimension(:,:,:) :: ffNonEq
double precision, allocatable, dimension(:,:,:,:) :: sigma
logical :: lshear=.false.
!**************************************************************************
contains
!**************************************************************************
subroutine allocate_shear()
  allocate(ffNonEq(Nx+2,Ny+2,qmom))
  ffNonEq=1.0d0
  allocate(sigma(Nx,Ny,2,2))
  sigma=0d0
  lshear=.true.
endsubroutine allocate_shear
!*************************************************************************
subroutine wall_shear()

integer :: q,i,j,mu,nu

do q=1,qmom
   do i=1,Nx
     do j=1,Ny
     ffNonEq(i,j,q)=ff(i,j,q)-ffEq(i,j,q)
       do mu=1,2
         do nu=1,2
           sigma(i,j,mu,nu)=(1.0d0 - (1.0d0/(2*tau))) &
                            !*ffNonEq(i,j,q)*ee(nu,q)*ee(mu,q) &
                            + sigma(i,j,mu,nu)
         enddo
       enddo
     enddo
   enddo
enddo

!write(*,*) sigma(5,5,1,1)
write(*,*) ffNonEq(5,5,5)
endsubroutine wall_shear
!*************************************************************************
subroutine free_shear()
  if (lshear .eqv. .true.) then
     deallocate(ffNonEq)
  endif
endsubroutine free_shear
!*************************************************************************

endmodule ShearStress
