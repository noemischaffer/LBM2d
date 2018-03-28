!
!The convention is that:
!  if a loop goes between 1-Nx+2 (1-Ny+2) the loop index is i (j)
!  if a loop goes between 2-Nx+1 (2-Ny+1) the loop index is k (l)
!
!***************************************************************
module BoundCond
use Cdata
use Avg
use messages
implicit none
private
public :: read_bc, initialize_bc, boundary_condition
character (len=labellen) :: bc_left='periodic',bc_right='periodic'
character (len=labellen) :: bc_bot='noslip',bc_top='noslip'
! It is left, right, bottom, top
namelist /bc_pars/ &
   bc_left,bc_right,bc_bot,bc_top

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
  !
! left boundary
!
  select case(bc_left)
     case('periodic')
        call pbcx
     case default
       call fatal_error('boundary condition','nothing other than pbc coded') 
  endselect
!
! bc_right should go here once we have something
! other than pbc.
!
  select case(bc_right)
     case('periodic')
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
subroutine boundary_condition()
  integer :: k,q,m,n,qnext,kq
  integer,dimension(3) :: qs
!
  select case(bc_left)
  case('periodic')
     call pbcx
  case default
     call fatal_error('boundary condition','nothing other than pbc coded') 
  endselect
!
! bc_right
!
  select case(bc_right)
  case('periodic')
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
     call bounceback('bot')
  case default
    call fatal_error('boundary before','bc not found')
  endselect
!
! top boundary
!
  select case(bc_top)
    case('periodic')
    case('noslip')
       call bounceback('top')
  case default
    call fatal_error('boundary before','bc not found')
 endselect
!
! call immersed_boundary()
!
endsubroutine boundary_condition
!***************************************************************
subroutine immersed_boundary()
  integer :: ix,iy,q,p
  do q=1,qmom
     do iy=2,Ny+1
        do ix=2,Nx+1
           if (is_solid(ix,iy).eq.0) then
              p=mirrorq(q)
              ff(ix,iy,p)=ff(ix,iy,q) ! this is bounce-back
           endif
        enddo
     enddo
  enddo
endsubroutine immersed_boundary
!***************************************************************
subroutine bounceback(boundary)
  character(len=3),intent(in) :: boundary
  integer :: q,m,n,ix,p,iq,iy
  integer,dimension(3) :: threeqs

  select case(boundary)
  case('bot')
     threeqs(1)=1
     threeqs(2)=2
     threeqs(3)=3
     iy=1
  case('top')
     threeqs(1)=7
     threeqs(2)=8
     threeqs(3)=9
     iy=Ny+2
  endselect
!
! At the bottom boundary we should stream in (except at the two corners)
! three components of momentum : 1,2,3  
!
  do iq=1,3
     q=threeqs(iq)
     p=mirrorq(q)
     do ix=2,Nx+1
        n = iy - ee_int(2,q)
        m = ix - ee_int(1,q)
        ff(ix,iy,p) = ff(m,n,q)
!           write(*,*) 'implementing bounce-back',ix,iy,m,n,q,ee_int(:,q)
     enddo
  enddo
!
! For the left corner  
!
  q=threeqs(1)
  p=mirrorq(q)
  n=iy-ee_int(2,q)
  m=1-ee_int(1,q)
  ff(1,iy,p) = ff(m,n,q)
!
! For the right corner  
!
  q=threeqs(3)
  p=mirrorq(q)
  n=iy-ee_int(2,q)
  m=Nx+2-ee_int(1,q)
  ff(Nx+2,iy,p) = ff(m,n,q)
!  
endsubroutine bounceback
!***************************************************************
endmodule BoundCond

