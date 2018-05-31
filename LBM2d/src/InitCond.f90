!***************************************************************
!
!The convention is that:
!  if a loop goes between 1 --  Nx+2 (1 -- Ny+2) the loop index is i (j)
!  if a loop goes between 2 --- Nx+1 (2 -- Ny+1) the loop index is k (l)
!
!***************************************************************
!
! A = (point_lb_x, point_lb_y)
! B = (point_lt_x, point_lt_y)
! D = (point_rb_x, point_rb_y)
!
module InitCond

use Cdata
use Avg
use messages
implicit none
private
public::read_param_init, initff,init_obstacle,update_surface,construct_surface

character (len=labellen) :: init_type='static'
character (len=labellen) :: obstacle_type='none'
!double precision :: radius=0.0d0, center_x=0.0d0, center_y=0.0d0
double precision :: point_lb_x=0.0d0,point_lb_y=0.0d0
double precision :: point_lt_x=0.0d0,point_lt_y=0.0d0
double precision :: point_rb_x=0.0d0,point_rb_y=0.0d0

namelist /init_cond/& 
    init_type, obstacle_type, radius, center_x, center_y, &
      point_lb_x,point_lb_y,point_lt_x,point_lt_y, &
      point_rb_x,point_rb_y

contains
!***************************************************************
subroutine read_param_init(unit,iostat)

!read the namelist init_con
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=init_cond, IOSTAT=iostat)

endsubroutine read_param_init
!***************************************************************
subroutine initff()

select case (init_type)
case('static')
   call static()
case default
   call fatal_error('initialize_ff',&
       'error in ff initialization! ')
endselect

endsubroutine initff
!***************************************************************
subroutine static()

  integer :: k,l,q

   do q=1,qmom
     do l=1,Ny+2
       do k=1,Nx+2
         ff(k,l,q) = weight(q)
         fftemp(k,l,q)=weight(q)
        enddo
      enddo
    enddo

endsubroutine static
!***************************************************************
subroutine init_obstacle()

  integer :: k, l, q
  double precision :: dist_square
  double precision :: AMdotAB, ABdotAB, AMdotAD, ADdotAD
!
! A = (point_lb_x, point_lb_y)
! B = (point_lt_x, point_lt_y)
! D = (point_rb_x, point_rb_y)
!

  select case (obstacle_type)
    case('none')
    case('point')
      do l=1,Ny+2
        do k=1,Nx+2
          is_solid(2,2) = 0
       enddo
     enddo
  case('circle')
    do l=2,Ny+1
      do k=2,Nx+1
        if(((xx(k-1)-center_x)**2 + (yy(l-1)-center_y)**2) .lt. radius**2) then
          is_solid(k,l)=1 ! this is a solid point
        endif
      enddo
    enddo
   case('rectangle')
     do l=2,Ny+1
       do k=2,Nx+1
         AMdotAB = (xx(k-1)-point_lb_x)*(point_lt_x-point_lb_x) + &
              (yy(l-1)-point_lb_y)*(point_lt_y-point_lb_y)
         ABdotAB = (point_lt_x-point_lb_x)*(point_lt_x-point_lb_x) + &
              (point_lt_y-point_lb_y)*(point_lt_y-point_lb_y)
         AMdotAD = (xx(k-1)-point_lb_x)*(point_rb_x-point_lb_x) + &
              (yy(l-1)-point_lb_y)*(point_rb_y-point_lb_y)
         ADdotAD = (point_rb_x-point_lb_x)*(point_rb_x-point_lb_x) + &
              (point_rb_y-point_lb_y)*(point_rb_y-point_lb_y)
         if((0.lt.AMdotAB).and.(AMdotAB.lt.ABdotAB).and.(0.lt.AMdotAD).and.(AMdotAD.lt.ADdotAD)) then
            is_solid(k,l)=1
         endif
      enddo
   enddo
   case default
     call fatal_error('initialize_obstacle',&
        'obstacle_type not coded! ')
  endselect

endsubroutine init_obstacle
!***************************************************************
subroutine construct_surface()

  allocate(surface(Nlarge,2))
  call update_surface()

endsubroutine construct_surface
!***************************************************************
subroutine update_surface()

  integer :: q,kk,m,n,k,l
  Nsurf=0

  do l=2,Ny+1
    do k=2,Nx+1
      !
      ! we check if any neighbour of a solid point (1)
      ! is a fluid point (-1). The the solid point is turned
      ! into a boundary point (0).
      !
      if(is_solid(k,l).eq.1) then
        kk=0
        do q=1,qmom
          m=k-ee_int(1,q)
          n=l-ee_int(2,q)
          if(is_solid(m,n).eq.-1) then
            if(kk.eq.0) then            
              Nsurf=Nsurf+1
              surface(Nsurf,1)=k
              surface(Nsurf,2)=l
              kk=kk+1
            endif
            is_solid(k,l) = 0
          endif
        enddo
      endif
    enddo
  enddo

  open(unit=12, file='fluid_point.txt', action='write', status='replace')
  open(unit=13, file='surface_point.txt', action='write', status='replace')
  open(unit=14, file='solid_point.txt', action='write', status='replace')

  do l=2,Ny-1
    do k=2, Nx-1
      if(is_solid(k,l).eq.-1) then
        write(12,*) k, ',', l
      endif 
      if(is_solid(k,l).eq.0) then
        write(13,*) k, ',', l
      endif 
      if(is_solid(k,l).eq.1) then
        write(14,*) k, ',', l
      endif 
    enddo
  enddo

  open(unit=13, file='surface_pointx.dat', action='write', status='replace')
  do l=2,Ny-1
    do k=2, Nx-1
      if(is_solid(k,l).eq.0) then
        write(13,*) k
      endif
    enddo
  enddo

  open(unit=12, file='surface_pointy.dat', action='write', status='replace')
    do l=2,Ny-1
      do k=2, Nx-1
        if(is_solid(k,l).eq.0) then
          write(12,*) l
        endif
      enddo
   enddo

  close(12)
  close(13)
  close(14)

endsubroutine update_surface
!***************************************************************
endmodule InitCond
