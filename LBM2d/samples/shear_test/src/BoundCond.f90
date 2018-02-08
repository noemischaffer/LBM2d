module BoundCond
use Cdata
use Avg
use messages
implicit none
private
public :: read_bc, initialize_bc, set_boundary_before,set_boundary_after
character (len=labellen) :: bc_left='periodic',bc_right='periodic'
character (len=labellen) :: bc_bot='noslip',bc_top='noslip'
character (len=labellen) :: bc_obstacle='staircase'
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
  if (bc_left.eq.'periodic') then
    if (bc_right.eq.'periodic') then
      dx=Lx/Nx
    else
    call fatal_error('initialize_bc','bc left and right dont match')
   endif 
  endif
  if (bc_bot.eq.'periodic') then
    if (bc_top.eq.'periodic') then
      dy=Ly/Ny
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
    case('noslip')
     qq=[7,8,9]
     call noslipy_before(qq,Ny+2)
  case default
     call fatal_error('boundary before','bc not found')
  endselect

  select case(bc_obstacle)
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
  integer :: i,q,m,n,qnext,iq
  i=1
  do iq=1,2
    q=qarray(iq)
    m=i-ee_int(1,q)
    n=jy-ee_int(2,q)
    ff(i,jy,q)=ff(m,n,q)
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '****************'
    !write(*,*) ff(m,n,q), m, n, q
    !write(*,*) '*******new step********'
  enddo
  do i=2,Nx+1
    do iq=1,3
      q=qarray(iq)
      m=i-ee_int(1,q)
      n=jy-ee_int(2,q)
      ff(i,jy,q)=ff(m,n,q)
  !    write(*,*) ff(i,jy,q), i, jy, q
  !    write(*,*) '****************'
  !    write(*,*) ff(m,n,q), m, n, q
  !    write(*,*) '*******new step********'
    enddo
  enddo
  i=Nx+2
  do iq=2,3
    q=qarray(iq)
    m=i-ee_int(1,q)
    n=jy-ee_int(2,q)
    ff(i,jy,q)=ff(m,n,q)
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '****************'
    !write(*,*) ff(m,n,q), m, n, q
    !write(*,*) '*******new step********'
  enddo

endsubroutine noslipy_before
!***************************************************************
subroutine noslip_obstacle_before()
  integer :: i,j,q,m,n
  do q=1, qmom
    do j=2,Ny+1
      do i=2,Nx+1
        if(is_solid(i,j).eq.1) then
          m=i-ee_int(1,q)
          n=j-ee_int(2,q)
          if((is_solid(i,n).eq.1).or.(is_solid(m,j).eq.1).or.(is_solid(m,n).eq.1)) then
             ff(i,j,q)=ff(m,n,q)
          endif
        endif
  !    write(*,*) ff(i,jy,q), i, jy, q
  !    write(*,*) '****************'
  !    write(*,*) ff(m,n,q), m, n, q
  !    write(*,*) '*******new step********'
      enddo
    enddo
  enddo
endsubroutine noslip_obstacle_before
!***************************************************************
subroutine set_boundary_after()
  integer,dimension(3) :: qq

!
! bottom boundary
!
  select case(bc_bot)
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
    case('noslip')
     qq=[7,8,9]
     call noslipy_after(qq,Ny+2)
    case default
     call fatal_error('boundary after','bc not found')
  endselect

  select case(bc_obstacle)
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
  integer :: i,q,p,iq
  i=1
  do iq=1,2
    q=qarray(iq)
    p=mirrorq(q)
    ff(i,jy,p)=ff(i,jy,q)
    !write(*,*) ff(i,jy,p), i, jy, p
    !write(*,*) '***************'
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '**********new step**********'
  enddo
  do i=2,Nx+1
    do iq=1,3
      q=qarray(iq)
      p=mirrorq(q)
      ff(i,jy,p)=ff(i,jy,q)
   ! write(*,*) ff(i,jy,p), i, jy, p
   ! write(*,*) '***************'
   ! write(*,*) ff(i,jy,q), i, jy, q
   ! write(*,*) '**********new step**********'
    enddo
  enddo
  i=Nx+2
  do iq=2,3
    q=qarray(iq)
    p=mirrorq(q)
    ff(i,jy,p)=ff(i,jy,q)
    !write(*,*) ff(i,jy,p), i, jy, p
    !write(*,*) '***************'
    !write(*,*) ff(i,jy,q), i, jy, q
    !write(*,*) '**********new step**********'
  enddo
endsubroutine noslipy_after
!***************************************************************
subroutine noslip_obstacle_after()
  integer :: i,j,q,p
  do q=1,qmom 
    do j=2,Ny+1
      do i=2,Nx+1
        if(is_solid(i,j).eq.1) then
          p=mirrorq(q)
          ff(i,j,p)=ff(i,j,q)
   ! write(*,*) ff(i,jy,p), i, jy, p
   ! write(*,*) '***************'
   ! write(*,*) ff(i,jy,q), i, jy, q
   ! write(*,*) '**********new step**********'
        endif
      enddo
    enddo
  enddo
endsubroutine noslip_obstacle_after
!***************************************************************
endmodule BoundCond

