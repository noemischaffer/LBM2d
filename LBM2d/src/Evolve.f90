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
  public :: stream,comp_equilibrium_BGK,get_feq,collision, allocate_shear,free_shear
  public :: allocate_force, free_force
  integer,allocatable,dimension(:) :: iin_pencil
  logical :: lshear=.false.
  logical :: lforce=.false.
!------the following are public ------------
!-------------------------------------------
contains
!***************************************************************
subroutine initialize_evolve()

endsubroutine initialize_evolve
!***************************************************************
subroutine finalize_evolve()

endsubroutine finalize_evolve
!**************************************************************************
subroutine allocate_shear()

  allocate(sigma(Nx+2,Ny+2,2,2))
  sigma=0.0d0
  lshear=.true.

endsubroutine allocate_shear
!**************************************************************************
subroutine allocate_force()

  allocate(forces(Nx+2,Ny+2,2))
  forces=0.0d0
  lforce=.true.

endsubroutine allocate_force
!*************************************************************************
subroutine stream()

  integer :: m,n,q,k,l,p

  fftemp=ff
  do q=1,qmom
     do l=2,Ny+1
        n=l-ee_int(2,q)
        do k=2,Nx+1
           m=k-ee_int(1,q)
           select case(is_solid(k,l))
           case(1) ! do nothing if solid point
           case(-1)
              ff(k,l,q) = fftemp(m,n,q) ! stream if fluid point
           case(0)
              p=mirrorq(q)
              ff(k,l,p) = fftemp(k,l,q) ! reflect back if surface point
              ! Put porosity here
           endselect
        enddo
     enddo
  enddo
endsubroutine stream
!***************************************************************
subroutine comp_equilibrium_BGK()

  use Avg
  use Force
  integer :: q,kk,ll
  double precision,dimension(2) :: uforced


  do q=1,qmom
    do ll=2,Ny+1
       do kk=2,Nx+1
          if(is_solid(kk,ll).eq.-1) then
             call get_uforce(kk,ll,uforced)
             call get_feq(q,uforced,rho(kk,ll),ffEq(kk,ll,q))
          endif
       enddo
    enddo
  enddo

endsubroutine comp_equilibrium_BGK
!***************************************************************
subroutine get_feq(q,uin,rhoin,fEq)

  double precision :: uxy, usqr
  integer, intent(in) :: q
  double precision,dimension(2),intent(in) :: uin
  double precision,intent(in) :: rhoin
  double precision,intent(out) :: fEq

  uxy=dot2d(ee(:,q),uin)
  usqr=dot2d(uin,uin)
  fEq = weight(q)*rhoin*(1.0d0+3.0d0*uxy &
              +(9.0d0/2.0d0)*(uxy**2) &
              -(3.0d0/2.0d0)*usqr &
               )
endsubroutine get_feq
!***************************************************************
subroutine collision()

  integer :: q,k,l,p,nu,mu
!
! If a point is solid(1)   : do nothing
!               surface(0) : bounce-back
!               fluid(-1)  : collision  
!
! Sigma is the shear stress tensor from Jäger et al 2017
  
  sigma=0.0d0
  do q=1,qmom
     do l=2,Ny+1
        do k=2,Nx+1
           select case(is_solid(k,l))
           case(1)
           case(0)
           case(-1)
             ff(k,l,q) = ff(k,l,q) + (ffEq(k,l,q)-ff(k,l,q))/tau
             do mu=1,2; do nu=1,2
                sigma(k,l,mu,nu)=(1.0d0 - (1.0d0/(2.0d0*tau))) &
                           *(ff(k,l,q)-ffEq(k,l,q))*ee(mu,q)*ee(nu,q) &
                           + sigma(k,l,mu,nu)
            enddo;enddo 
          endselect
        enddo
     enddo
  enddo

endsubroutine collision
!***************************************************************
subroutine free_shear()

  if (lshear .eqv. .true.) then
     deallocate(sigma)
  endif

endsubroutine free_shear
!***************************************************************
subroutine free_force()

  if (lforce .eqv. .true.) then
     deallocate(forces)
  endif

endsubroutine free_force
!***************************************************************
endmodule Evolve
