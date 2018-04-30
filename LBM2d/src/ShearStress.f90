!LBM2D

!Calculates the wall shear stress according to Equation 12 in Jäger+2017
!**************************************************************************
module ShearStress

use Cdata
use Evolve
implicit none
private
public :: allocate_shear,  wall_shear, free_shear

logical :: lshear=.false.
!**************************************************************************
contains
!**************************************************************************
subroutine allocate_shear()
  allocate(sigma(Nx,Ny,2,2))
  sigma=0.0d0
  lshear=.true.
endsubroutine allocate_shear
!*************************************************************************
subroutine wall_shear()


use Sub

integer :: q,mu,nu,jsurf,k,l,m,n
double precision :: x0,y0,center_x,center_y,radius,cos_theta,sin_theta
double precision, dimension(2) :: tang_vect, norm_vect
double precision, dimension(Nx,Ny):: lift
double precision, dimension(Nx,Ny):: lift_sum
double precision, dimension(Nx,Ny):: drag
double precision, dimension(Nx,Ny):: drag_sum
double precision, dimension(Nx,Ny,2):: lift_vect
double precision, dimension(Nx,Ny,2):: drag_vect


tang_vect(1) = 0.0d0
tang_vect(2) = 0.0d0
norm_vect(1) = 0.0d0
norm_vect(2) = 0.0d0

center_x = 128.0d0
center_y = 128.0d0
radius = 30.0d0
x0 = center_x
y0 = center_y
sin_theta = 0.0d0
cos_theta = 0.0d0


do q=1,qmom
  do jsurf=1, Nsurf
    k=surface(jsurf,1)
    l=surface(jsurf,2)
    sin_theta = abs(l-y0) / radius
    cos_theta = abs(k-x0) / radius
    tang_vect(1) = - sin_theta ! -sin(theta)
    tang_vect(2) =   cos_theta !  cos(theta)
    norm_vect(1) = - cos_theta ! -cos(theta)
    norm_vect(2) = - sin_theta ! -sin(theta)
    do mu=1,2
      do nu=1,2
        sigma(k,l,mu,nu)=(1.0d0 - (1.0d0/(2*tau))) &
                        *(ff(k,l,q)-ffEq(k,l,q))*ee(mu,q)*ee(nu,q) &
                        + sigma(k,l,mu,nu)

       lift(k,l) = sigma(k,l,mu,nu)*norm_vect(nu)
       lift_sum(k,l) = lift_sum(k,l)+lift(k,l)
       drag(k,l) = sigma(k,l,mu,nu)*tang_vect(nu)
       drag_sum(k,l) = drag_sum(k,l)+drag(k,l)
       enddo
       lift_vect(k,l,mu) = lift_sum(k,l)
       drag_vect(k,l,mu) = drag_sum(k,l)
     enddo
  enddo
enddo

  open(unit=13, file='drag.dat', action='write', status='replace')
  do k=1, Nx
    write(13,*) (drag_vect(k,l,1), l=1,Ny)
  enddo
  close(13)


endsubroutine wall_shear
!*************************************************************************
!subroutine drag_lift()

!do q=1,qmom
!  do jsurf=1, Nsurf
!    k=surface(jsurf,1)
!    l=surface(jsurf,2)
!    drag (k,l) = (f(k+1,l,1)-f(k+1,l,3)+f(k+1,l,5)-f(k+1,l,6)-f(k+1,l,7)+f(k+1,l,8)) &
!               + (

! enddo
!enddo



!endsubroutine drag_lift
!*************************************************************************
subroutine free_shear()
  if (lshear .eqv. .true.) then
     deallocate(sigma)
  endif
endsubroutine free_shear
!*************************************************************************

endmodule ShearStress
