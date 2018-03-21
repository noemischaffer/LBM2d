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
!
! left boundary
!
  select case(bc_left)
     case('periodic')
        call pbcy
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
     call pbcx
  case('noslip')
     ff(:,1,:)=0.0d0
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
       ff(:,Ny+2,:)=0.
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
  is_solid(Nx+2,:)=is_solid(1,:)
endsubroutine pbcx
!***************************************************************
subroutine boundary_condition()
  integer :: k,q,m,n,qnext,kq
!
  select case(bc_left)
  case('periodic')
     call pbcy
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
     call pbcx
  case('noslip')
     call stream_bot()
  case default
    call fatal_error('boundary before','bc not found')
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
endsubroutine boundary_condition
!***************************************************************
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
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the left corner  
!
  do q=1,2
     l=1
     n=l-ee_int(2,q)
     do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the right corner  
!
  do q=2,3
     l=1
     n=l-ee_int(2,q)
     do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
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
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the left corner  
!
  do q=7,8
     l=Ny+2
     n=l-ee_int(2,q)
     do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
!
! For the right corner  
!
  do q=8,9
     l=Ny+2
     n=l-ee_int(2,q)
     do k=2,Nx+1
        m=k-ee_int(1,q)
        !
        !stream only if the point you are
        ! streaming from is a fluid point(-1)
        !
        if (is_solid(m,n).eq.-1) then
           fftemp(k,l,q) = ff(m,n,q)
        endif
     enddo
  enddo
endsubroutine stream_top
!***************************************************************
endmodule BoundCond

