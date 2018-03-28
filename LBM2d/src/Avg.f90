! LBM2D
!
!  calculates averages, i.e., velocity and density
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
module Avg
  use Cdata
  use Sub
  implicit none
  private
!------the following are public ------------
  public :: allocate_avg,calc_avg,free_avg, vorticity, rwrite_density_uu
  public :: uu,rho
!-------------------------------------------
  double precision, allocatable, dimension(:,:,:) :: uu
  double precision, allocatable, dimension(:,:) :: rho
!  double precision, allocatable, dimension(:,:) :: curl_uu
  logical :: lavg=.false.
contains
!***************************************************************
subroutine allocate_avg()
  allocate(uu(Nx+2,Ny+2,2))
  uu=0.0d0
  allocate(rho(Nx+2,Ny+2))
  rho=0.0d0
  allocate(curl_uu(Nx+2,Ny+2))
  curl_uu=0.0d0
  lavg=.true.
endsubroutine allocate_avg
!***************************************************************
subroutine calc_avg()
  integer :: q 
  double precision :: uxmax,uymax
  uu=0.0d0;rho=0.0d0

  do q=1,qmom
    rho(:,:) = rho(:,:)+ff(:,:,q)
  enddo
  uu(:,:,1) = (ff(:,:,3)+ff(:,:,6)+ff(:,:,9)&
              -ff(:,:,1)-ff(:,:,4)-ff(:,:,7))
  uu(:,:,2) = (ff(:,:,7)+ff(:,:,8)+ff(:,:,9)&
              -ff(:,:,1)-ff(:,:,2)-ff(:,:,3))
  uu(:,:,1) = uu(:,:,1)/rho(:,:)
  uu(:,:,2) = uu(:,:,2)/rho(:,:)
  !  write(*,*) maxval(ff(:,:,:))
endsubroutine calc_avg
!***************************************************************
subroutine vorticity()

integer :: i,j,q,m,n
  curl_uu=0.0d0
  do j=2,Ny-1
    do i=2, Nx-1
      if(is_solid(i,j).eq.1) then
        curl_uu(i,j)=0.0d0
        do q=1,qmom
          m=i-ee_int(1,q)
          n=j-ee_int(2,q)
          if(is_solid(m,n).eq.1) then
             curl_uu(m,n)=0.0d0
          endif
        enddo 
      else
         curl_uu(i,j) = (uu(i+1,j,2)-uu(i-1,j,2))/(2.0d0*dx) - (uu(i,j+1,1)-uu(i,j-1,1))/(2.0d0*dy)
     endif
    enddo
  enddo


!  do j=2,Ny-1
!    do i=2, Nx-1
!      curl_uu(i,j) = (uu(i+1,j,2)-uu(i-1,j,2))/(2.0d0*dx) - (uu(i,j+1,1)-uu(i,j-1,1))/(2.0d0*dy)
!    enddo
!  enddo



endsubroutine vorticity 
!***************************************************************  
subroutine rwrite_density_uu()
!  Reads and registers print parameters

  integer :: k,l 

  open(unit=2, file='dimensions.txt', action='write', status='replace')
  write(2,*) Nx
  write(2,*) Ny
  close(2)
  open(unit=1, file='vorticity.txt', action='write', status='replace')
  !write(11,*) rho(:,1)
  !write(11,*) rho(1,:)
  !write(1,*) Nx
  !write(1,*) Ny
  do l=2,Ny-1
    do k=2,Nx-1
      write(1,*) k,l,curl_uu(k,l)
    enddo
  enddo
  close(1)

  open(unit=10, file='vorticity_idl.txt', action='write', status='replace')
  !write(11,*) rho(:,1)
  !write(11,*) rho(1,:)
  write(10,*) Nx
  write(10,*) Ny
  write(10,*) curl_uu(:,:)
  close(10)

  open(unit=11, file='velocity.txt', action='write', status='replace')

  write(11,*) Nx
  write(11,*) Ny
  write(11,*) uu(:,:,1)
  write(11,*) uu(:,:,2)
  !write(*,*) maxval(uu(:,:,1)),maxval(uu(:,:,2))
  close(11)

  open(unit=12, file='velocity_math_vx.txt', action='write', status='replace')
  write(12,*) uu(:,:,1)
  close(12)
  open(unit=13, file='velocity_math_vy.txt', action='write', status='replace')
  write(13,*) uu(:,:,2)
  close(13)


endsubroutine rwrite_density_uu
!***************************************************************  
subroutine free_avg()
  if (lavg .eqv. .true.) then
     deallocate(uu); deallocate(rho)
  endif
endsubroutine free_avg
!***************************************************************  
endmodule Avg  
