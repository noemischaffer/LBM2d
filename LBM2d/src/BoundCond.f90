!
!The convention is that:
!  if a loop goes between 1-Nx+2 (1-Ny+2) the loop index is i (j)
!  if a loop goes between 2-Nx+1 (2-Ny+1) the loop index is k (l)
!
!***************************************************************
module BoundCond
use Cdata
use InitCond
use Avg
use messages
implicit none
private
public :: read_bc, initialize_bc, set_boundary_before,set_boundary_after
character (len=labellen) :: bc_left='periodic',bc_right='periodic'
character (len=labellen) :: bc_bot='noslip',bc_top='noslip'
character (len=labellen) :: bc_obstacle='none'

! It is left, right, bottom, top
namelist /bc_pars/ &
   bc_left,bc_right,bc_bot,bc_top,bc_obstacle

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
! Boundary types: rb, lb, tb, bb


integer :: i,j

  if (bc_left.eq.'periodic') then
    if (bc_right.eq.'periodic') then
      dx=Lx/Nx
      do i=1, Nx
        xx(i)=i*dx
      enddo
    else
    call fatal_error('initialize_bc','bc left and right dont match')
   endif 
  endif
  if (bc_bot.eq.'periodic') then
    if (bc_top.eq.'periodic') then
      dy=Ly/Ny
      do j=1, Ny
        yy(j)=j*dy
      enddo
    else
    call fatal_error('initialize_bc','bc bot and top dont match')
   endif 
  endif

endsubroutine initialize_bc
!***************************************************************
subroutine set_boundary_before()
  integer,dimension(3) :: qq 

!
! left boundary
!
  select case(bc_left)
     case('periodic')
      ff(Nx+2,:,:)=ff(2,:,:)
      ff(1,:,:)=ff(Nx+1,:,:)
     case default
       call fatal_error('boundary before','nothing other than pbc coded') 
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
!
! bottom boundary
!
  select case(bc_bot)
  case('periodic')
    ff(:,Ny+2,:)=ff(:,2,:)
    ff(:,1,:)=ff(:,Ny+1,:)
  case('noslip')
    qq=[1,2,3]
    call noslipy_before(qq,1)
  case default
    call fatal_error('boundary before','bc not found')
  endselect
!
! top boundary
!
  select case(bc_top)
    case('periodic')
    case('noslip')
     qq=[7,8,9]
     call noslipy_before(qq,Ny+2)
  case default
    call fatal_error('boundary before','bc not found')
  endselect

  select case(bc_obstacle)
    case('none')
    case('staircase')
    call noslip_obstacle_before()
  case default
     call fatal_error('boundary obstacle','bc not found')
  endselect

endsubroutine set_boundary_before
!***************************************************************
subroutine noslipy_before(qarray,jy)
  integer,dimension(3),intent(in) :: qarray
  integer, intent(in) :: jy 
  integer :: k,q,m,n,qnext,kq
  k=1
  do kq=1,2
    q=qarray(kq)
    m=k-ee_int(1,q)
    n=jy-ee_int(2,q)
    ff(k,jy,q)=ff(m,n,q)
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '****************'
    !write(*,*) ff(m,n,q), m, n, q
    !write(*,*) '*******new step********'
  enddo
  do k=2,Nx+1
    do kq=1,3
      q=qarray(kq)
      m=k-ee_int(1,q)
      n=jy-ee_int(2,q)
      ff(k,jy,q)=ff(m,n,q)
  !    write(*,*) ff(i,jy,q), i, jy, q
  !    write(*,*) '****************'
  !    write(*,*) ff(m,n,q), m, n, q
  !    write(*,*) '*******new step********'
    enddo
  enddo
  k=Nx+2
  do kq=2,3
    q=qarray(kq)
    m=k-ee_int(1,q)
    n=jy-ee_int(2,q)
    ff(k,jy,q)=ff(m,n,q)
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '****************'
    !write(*,*) ff(m,n,q), m, n, q
    !write(*,*) '*******new step********'
  enddo

endsubroutine noslipy_before
!***************************************************************
subroutine noslip_obstacle_before()
  integer :: k,l,q,m,n,jsurf

open(unit=15, file='stream_to.txt', action='write', status='replace')
open(unit=16, file='stream_from.txt', action='write', status='replace')


  do q=1, qmom
    do jsurf=1,Nsurf
      k=surface(jsurf,1)
      l=surface(jsurf,2)
      m=k-ee_int(1,q)
      n=l-ee_int(2,q)
   !   if(is_solid(m,n).ne.1) then
        ff(k,l,q)=ff(m,n,q)
        write(15,*) k, ',', l
        write(16,*) m, ',', n
         !refl_point(Nsurf,q)=1
         ! if(jsurf.eq.Nsurf-2) then
         !  write(*,*) i,j,q, 'streamed', ff(m,n,q),  ' from', m, n,q
         ! endif
   !   endif
     enddo
   enddo

close(15)
close(16)

endsubroutine noslip_obstacle_before
!***************************************************************
subroutine set_boundary_after()
  integer,dimension(3) :: qq

!
! bottom boundary
!
  select case(bc_bot)
  case('periodic')
  case('noslip')
    qq=[1,2,3]
    call noslipy_after(qq,1)
  case default
     call fatal_error('boundary after','bc not found')
  endselect
!
! top boundary
!
  select case(bc_top)
  case('periodic')
  case('noslip')
    qq=[7,8,9]
    call noslipy_after(qq,Ny+2)
  case default
    call fatal_error('boundary after','bc not found')
  endselect

  select case(bc_obstacle)
    case('none')
    case('staircase')
      call noslip_obstacle_after()
    case default
      call fatal_error('obstacle boundary after','bc not found')
  endselect
endsubroutine set_boundary_after
!***************************************************************
subroutine noslipy_after(qarray,jy)
  integer,dimension(3),intent(in) :: qarray
  integer, intent(in) :: jy
  integer :: k,q,p,kq
  k=1
  do kq=1,2
    q=qarray(kq)
    p=mirrorq(q)
    ff(k,jy,p)=ff(k,jy,q)
    !write(*,*) ff(i,jy,p), i, jy, p
    !write(*,*) '***************'
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '**********new step**********'
  enddo
  do k=2,Nx+1
    do kq=1,3
      q=qarray(kq)
      p=mirrorq(q)
      ff(k,jy,p)=ff(k,jy,q)
   ! write(*,*) ff(i,jy,p), i, jy, p
   ! write(*,*) '***************'
   ! write(*,*) ff(i,jy,q), i, jy, q
   ! write(*,*) '**********new step**********'
    enddo
  enddo
  k=Nx+2
  do kq=2,3
    q=qarray(kq)
    p=mirrorq(q)
    ff(k,jy,p)=ff(k,jy,q)
    !write(*,*) ff(i,jy,p), i, jy, p
    !write(*,*) '***************'
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '**********new step**********'
  enddo
endsubroutine noslipy_after
!***************************************************************
subroutine noslip_obstacle_after()
  integer :: k,l,q,p,jsurf

  open(unit=10, file='bounce_back.txt', action='write', status='replace')
   
  do q=1,qmom 
    do jsurf=1,Nsurf
      k=surface(jsurf,1)
      l=surface(jsurf,2)
   !   if(refl_point(jsurf,q).eq.1) then
         p=mirrorq(q)
         ff(k,l,p)=ff(k,l,q)
   !      if((q.eq.4).or.(q.eq.6)) then
         !write(*,*) i, j, p, 'reflected from',  q
   !      endif
         write(10,*) k, ',', l
   !   endif
     enddo
  enddo
  close(10)
endsubroutine noslip_obstacle_after
!***************************************************************
endmodule BoundCond

