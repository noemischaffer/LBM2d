!LBM2D

!Calculates the wall shear stress according to Equation 12 in Jäger+2017
!**************************************************************************
module ShearStress

use Cdata
use Evolve
use Avg
use InitCond
implicit none
private
public :: wall_shear, theo_solution
!**************************************************************************
contains
!**************************************************************************
subroutine wall_shear()

!
! The wall shear force is calculated at tau=normal vector dot sigma
!
use Sub

integer :: q,mu,nu,k,l,m,n 
!double precision, dimension(Nx+2,Ny+2,2):: force
double precision, dimension(2) :: force_tmp, unnorm
double precision :: fdrag,flift,magn, unnorm_x, unnorm_y, magn_sq
double precision, dimension(qmom) :: norm_x

forces=0.0d0

  do l=2,Ny+1
    do k=2,Nx+1
      unnorm_x = 0.0d0
      unnorm_y = 0.0d0
      if(is_solid(k,l).eq.0) then
        do q=1,qmom
          n=l-ee_int(2,q)
          m=k-ee_int(1,q)
          if((is_solid(m,n).eq.-1)) then 
            call dotmv(sigma(m,n,:,:),ee(:,q),force_tmp)
            forces(k,l,:)=forces(k,l,:)+force_tmp
            unnorm_x = unnorm_x + ee(1,q)
            unnorm_y = unnorm_y + ee(2,q)
          endif
        enddo
        magn_sq = unnorm_x**2+unnorm_y**2
        forces(k,l,:)=forces(k,l,:)/sqrt(magn_sq) 
        fdrag=fdrag+forces(k,l,1)
        flift=flift+forces(k,l,2)
      endif
    enddo
  enddo


  write(*,*) fdrag, ' drag force'
  write(*,*) flift, ' lift force'

  open(unit=9, file='sigma_xx.dat', action='write', status='replace')
  do k=1, Nx
    write(9,*) (sigma(k,l,1,1), l=1,Ny)
  enddo
  close(9)

  open(unit=8, file='sigma_yx.dat', action='write', status='replace')
  do k=1, Nx
    write(8,*) (sigma(k,l,2,1), l=1,Ny)
  enddo
  close(8)

  open(unit=7, file='sigma_xy.dat', action='write', status='replace')
  do k=1, Nx
    write(7,*) (sigma(k,l,1,2), l=1,Ny)
  enddo
  close(7)

  open(unit=10, file='sigma_yy.dat', action='write', status='replace')
  do k=1, Nx
    write(10,*) (sigma(k,l,2,2), l=1,Ny)
  enddo
  close(10)

  open(unit=12, file='drag.dat', action='write', status='replace')
  do k=1, Nx
    write(12,*) (forces(k,l,1), l=1,Ny)
  enddo
  close(12)
 
  open(unit=13, file='lift.dat', action='write', status='replace')
  do k=1, Nx
    write(13,*) (forces(k,l,2), l=1,Ny)
  enddo
  close(13)

endsubroutine wall_shear
!*************************************************************************
subroutine theo_solution()

! Theoretical stream function using Eq 2.11 from Proudman and Pearson 1956

integer :: l,k
double precision :: polar_r,kk,ll, sin_theta
double precision, dimension(Nx+2,Ny+2) :: stream_theo 
double precision, dimension(Nx+2,Ny+2) :: ux_theo, uy_theo
double precision :: Re, ratio

! vx= - d(psi)/d(y)
! vy=  d(psi)/d(x)

  Re = (v_const*radius)/((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))
  write(*,*) Re, 'this is the reynolds number'
 
  do l=2,Ny+1
    ll=l
    do k=2,Nx+1
      kk=k
      ratio = (ll-center_y)/(kk-center_x) 
      sin_theta = ratio/sqrt(1.0d0+(ratio*ratio))
      polar_r = (sqrt(abs(kk-center_x)**2 + abs(ll-center_y)**2))/radius

      stream_theo(k,l) = - (1.0d0 / (2.0d0*log10(Re))) * ((1.0d0/polar_r) - polar_r + &
                 (2.0d0*polar_r*log10(polar_r/(1.0d0+polar_r*Re))))*sin_theta
      enddo
  enddo

  do l=2,Ny
    do k=2,Nx
      ux_theo(k,l) = - (stream_theo(k,l+1)-stream_theo(k,l)) / dy
      uy_theo(k,l) =  (stream_theo(k+1,l)-stream_theo(k,l)) / dx
    enddo
  enddo

  open(unit=1, file='stream_theo.dat', action='write', status='replace')
  do k=1, Nx
    write(1,*) (stream_theo(k,l), l=1,Ny)
  enddo
  close(1)

  open(unit=2, file='ux_theo.dat', action='write', status='replace')
  do k=1, Nx
    write(2,*) (ux_theo(k,l), l=1,Ny)
  enddo
 close(2)

  open(unit=3, file='uy_theo.dat', action='write', status='replace')
  do k=1, Nx
    write(3,*) (uy_theo(k,l), l=1,Ny)
  enddo
  close(3)

endsubroutine theo_solution
!*************************************************************************
endmodule ShearStress
