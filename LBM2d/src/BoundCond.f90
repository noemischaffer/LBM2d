!
!The convention is that:
!  if a loop goes between 1-Nx+2 (1-Ny+2) the loop index is i (j)
!  if a loop goes between 2-Nx+1 (2-Ny+1) the loop index is k (l)
!
!***************************************************************
module BoundCond
use Cdata
use Evolve
use Avg
use messages
implicit none
private
public :: read_bc,initialize_bc, boundary_condition, const_vleft,const_vright
character (len=labellen) :: bc_left='periodic',bc_right='periodic'
character (len=labellen) :: bc_bot='noslip',bc_top='noslip'
double precision :: radius
! It is left, right, bottom, top
namelist /bc_pars/ &
   bc_left,bc_right,bc_bot,bc_top, radius

!***************************************************************
contains
!***************************************************************
subroutine read_bc(unit,iostat)
! read the namelist bc_pars
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=bc_pars, IOSTAT=iostat)
endsubroutine read_bc
!***************************************************************
subroutine initialize_bc()
integer :: q
double precision :: v
!
! left boundary
!
  select case(bc_left)
     case('periodic')
        call pbcx
     case('dirichlet')
       !v = (Re*((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))) / (2*radius)
       !call const_vleft(v)
     case('const_v')
     case default
       call fatal_error('boundary condition','nothing other than pbc coded') 
  endselect
!
! bc_right should go here once we have something
! other than pbc.
!
  select case(bc_right)
     case('periodic')
     case('right_copied')
    !    call right_copy()
     case('dirichlet')
       !v = (Re*((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))) / (2*radius)
       !call const_vright(v)
     case default
       call fatal_error('boundary before','nothing other than pbc coded')
  endselect
!
! bottom boundary
!
  select case(bc_bot)
  case('periodic')
     call pbcy
  case('noslip')
     do q=1,qmom
       ff(:,1,q)=weight(q)
     enddo
     is_solid(:,1)=1
  case default
    call fatal_error('boundary before','bc not found')
  endselect
!
! top boundary
!
  select case(bc_top)
    case('periodic')
    case('noslip')
       do q=1,qmom
         ff(:,Ny+2,q)=weight(q)
       enddo
       is_solid(:,Ny+2)=1
  case default
    call fatal_error('boundary before','bc not found')
  endselect

endsubroutine initialize_bc
!***************************************************************
subroutine pbcy()
  ff(:,1,:)=ff(:,Ny+1,:)
  is_solid(:,1)=is_solid(:,Ny+1)
  ff(:,Ny+2,:)=ff(:,2,:)
  is_solid(:,Ny+2)=is_solid(:,2)
endsubroutine pbcy
!***************************************************************
subroutine pbcx()
  ff(1,:,:)=ff(Nx+1,:,:)
  is_solid(1,:)=is_solid(Nx+1,:)
  ff(Nx+2,:,:)=ff(2,:,:)
  is_solid(Nx+2,:)=is_solid(2,:)
endsubroutine pbcx
!***************************************************************
subroutine right_copy()
  ff(Nx+2,:,:)=ff(Nx+1,:,:)
  is_solid(Nx+2,:)=is_solid(Nx+1,:)
endsubroutine right_copy
!***************************************************************
subroutine boundary_condition()
!
  double precision :: v
  double precision, dimension(2) :: vel
  double precision :: rho0
!
! bc_left
!
  select case(bc_left)
  case('periodic') 
     call pbcx()
  case('dirichlet')
     v = (Re*((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))) / (2*radius)
     call const_vleft(v)
  case('const_v')
     vel(1) = (Re*((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))) / (2*radius)
     vel(2) = 0.0d0
     call const_velocity(vel)
  case default
     call fatal_error('boundary condition','nothing other than pbc coded') 
  endselect
!
! bc_right
!
  select case(bc_right)
  case('periodic')
  case('dirichlet')
    ! v = (Re*((1.0d0/3.0d0)*(tau-(1.0d0/2.0d0)))) / (2*radius)
    rho0=1.0d0
     call const_vright(rho0)
  case('right_copied')
     call right_copy()
  case default
     call fatal_error('boundary before','nothing other than pbc coded')
  endselect
!
! bottom boundary
!
  select case(bc_bot)
  case('periodic')
     call pbcy
  case('noslip')
     call stream_bot()
  case default
     call fatal_error('boundary before','nothing other than pbc coded')
  endselect
!
! top boundary
!
  select case(bc_top)
    case('periodic')
    case('noslip')
       call stream_top()
  case default
    call fatal_error('boundary before','bc not found')
 endselect
!

!call immersed_boundary()

endsubroutine boundary_condition
!***************************************************************
subroutine immersed_boundary()
  integer :: k,l,q,p
  do q=1,qmom
     do l=2,Ny+1
        do k=2,Nx+1
           if (is_solid(k,l).eq.0) then
              p=mirrorq(q)
              ff(k,l,p)=ff(k,l,q) ! this is bounce-back
           endif
        enddo
     enddo
  enddo
endsubroutine immersed_boundary
!****************************************************************
subroutine stream_bot()
  integer :: k,l,q,m,n
!
! At the bottom boundary we should stream in (except at the two corners)
! three components of momentum : 1,2,3  
!
  do q=1,3
     l=1
     n=l-ee_int(2,q)
     do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the left corner  
!
  k=1
  do q=1,2
     l=1
     n=l-ee_int(2,q)
     !do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     !enddo
  enddo
!
! For the right corner  
!
  k=Nx+2
  do q=2,3
     l=1
     n=l-ee_int(2,q)
     !do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     !enddo
  enddo
endsubroutine stream_bot
!***************************************************************
subroutine stream_top()
  integer :: k,l,q,m,n
!
! At the top boundary we should stream in (except at the two corners)
! three components of momentum : 7,8,9  
!
  do q=7,9
     l=Ny+2
     n=l-ee_int(2,q)
    do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the left corner  
!
  k=1
  do q=7,8
     l=Ny+2
     n=l-ee_int(2,q)
   !  do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
    ! enddo
  enddo
!
! For the right corner  
! 
  k=Nx+2
  do q=8,9
     l=Ny+2
     n=l-ee_int(2,q)
     !do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     !enddo
  enddo
endsubroutine stream_top
!***************************************************************
subroutine const_vleft(v)


integer :: k,l,m,n,q
double precision :: v
!
! At the left boundary we should stream in (except at the two corners)
! three components of momentum : 1, 4, 7  
!

write(*,*) 'Dirichlet constant veclocity on left'

  k=1
  do q=1,7,3
     m=k-ee_int(1,q)
    do l=2,Ny+1
        n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo

!
!For top corner
!
  l=Ny+2
  k=1
  do q=4,7,3
     m=k-ee_int(1,q)
     n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
     if (is_solid(m,n).eq.-1) then
        ff(k,l,q) = ff(m,n,q)
     endif
  enddo

!
!For bottom corner
!
  l=1
  k=1
  do q=1,4,3
     m=k-ee_int(1,q)
     n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
     if (is_solid(m,n).eq.-1) then
        ff(k,l,q) = ff(m,n,q)
     endif
  enddo

!
!Then the rest of the fs and the density on the wall is
!calculated based on Sukop and Thorne.
!Here, v is the x component of the in/outlet velocity.
!The y component is 0.
!v = (Re*nu)/L where L is the object diameter
!
   k=1
   do l=2,Ny+1
      rho(k,l) = (2.0d0*(ff(k,l,1)+ff(k,l,4)+ff(k,l,7))+ff(k,l,2)+ &
                 ff(k,l,5)+ff(k,l,8))/(1.0d0-v)
      ff(k,l,6) = ff(k,l,4)+(2.0d0/3.0d0)*v*rho(k,l)
      ff(k,l,3) = ((1.0d0/6.0d0)*v*rho(k,l))+ff(k,l,7)- &
                 (1.0d0/2.0d0)*(ff(k,l,2)-ff(k,l,8))
      ff(k,l,9) = ((1.0d0/6.0d0)*v*rho(k,l))+ff(k,l,1)- &
                 (1.0d0/2.0d0)*(ff(k,l,8)-ff(k,l,2))
   enddo

endsubroutine const_vleft
!***************************************************************
subroutine const_vright(rho0)


integer :: k,l,m,n,q
double precision :: v
double precision :: rho0
!
! At the left boundary we should stream in (except at the two corners)
! three components of momentum : 3, 6, 9  
!

!  write(*,*) 'the in/out-flow velocty is ', v

write(*,*) 'Dirichlet constant pressure on right'
  k=Nx+2
  do q=3,9,3
     m=k-ee_int(1,q)
    do l=2,Ny+1
        n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           ff(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo

!
!For top corner
!
  l=Ny+2
  k=Nx+2
  do q=6,9,3
     m=k-ee_int(1,q)
     n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
     if (is_solid(m,n).eq.-1) then
        ff(k,l,q) = ff(m,n,q)
     endif
  enddo

!
!For bottom corner
!
  l=1
  k=Nx+2
  do q=3,6,3
     m=k-ee_int(1,q)
     n=l-ee_int(2,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
     if (is_solid(m,n).eq.-1) then
        ff(k,l,q) = ff(m,n,q)
     endif
  enddo

!
!Then the rest of the fs and the density on the wall is
!calculated based on Sukop and Thorne.
!Here, v is the x component of the in/outlet velocity.
!The y component is 0.
!v = (Re*nu)/L where L is the object diameter
!
   k=Nx+2
   do l=2,Ny+1
      !rho(k,l) = (2.0d0*(ff(k,l,3)+ff(k,l,6)+ff(k,l,9))+ff(k,l,2)+ &
       !          ff(k,l,5)+ff(k,l,8))/(1.0d0+v)
      uu(k,l,1) = -1.0d0+(2.0d0*(ff(k,l,3)+ff(k,l,6)+ff(k,l,9))+ff(k,l,2)+ &
                 ff(k,l,5)+ff(k,l,8))/rho0
      ff(k,l,4) = ff(k,l,6)-((2.0d0/3.0d0)*uu(k,l,1)*rho0)
      ff(k,l,7) = -((1.0d0/6.0d0)*uu(k,l,1)*rho0)+ff(k,l,3)- &
                 (1.0d0/2.0d0)*(ff(k,l,8)-ff(k,l,2))

      ff(k,l,1) = -((1.0d0/6.0d0)*uu(k,l,1)*rho0)+ff(k,l,9)- &
                 (1.0d0/2.0d0)*(ff(k,l,2)-ff(k,l,8))
   enddo

endsubroutine const_vright
!***************************************************************
subroutine const_velocity(vel)
double precision, dimension(2) :: vel
double precision :: rho0
integer :: q,k,l

rho0 = 1.0d0
k=1
do q=1,qmom
  do l=2,Ny+1
    call get_feq(q,vel,rho0,ff(k,l,q))
  enddo
enddo
  


endsubroutine const_velocity
!***************************************************************
endmodule BoundCond
