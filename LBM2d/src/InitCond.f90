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
public::mass_loss, eroded_obstacle

character (len=labellen) :: init_type='static'

namelist /init_cond/& 
    init_type, obstacle_type, radius, center_x, center_y, &
      point_lb_x,point_lb_y,point_lt_x,point_lt_y, &
      point_rb_x,point_rb_y, &
      overlap, &
      erosion_param, aa, bb
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
        is_solid(2,2)=0
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
  case('ellipse')
    do l=2,Ny+1
      do k=2,Nx+1
        if(((((xx(k-1)-center_x)**2)/aa**2) + (((yy(l-1)-center_y)**2)/bb**2)) .lt. 1) then
          is_solid(k,l)=1 ! this is a solid point
        endif
      enddo
    enddo
  case('67P')
    do l=2,Ny+1
      do k=2,Nx+1
        if(((xx(k-1)-(center_x+(radius-overlap)))**2 + (yy(l-1)-center_y)**2).lt. (radius)**2.or. &
          ((xx(k-1)-(center_x-radius))**2 + (yy(l-1)-center_y)**2).lt. (radius)**2) then
          is_solid(k,l)=1 ! this is a solid point
        endif
      enddo
    enddo  
    write(*,*) overlap, 'this is the overlap'
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
   case('restart')
     call eroded_obstacle()
   case default
     call fatal_error('initialize_obstacle',&
        'obstacle_type not coded! ')
  endselect

endsubroutine init_obstacle
!***************************************************************
subroutine construct_surface()

  allocate(surface(Nlarge,2))
  allocate(tau_magn(Nx+1,Ny+1))
  call update_surface()

endsubroutine construct_surface
!***************************************************************
subroutine update_surface()

use Evolve

  integer :: q,kk,m,n,k,l,mu,nu,j, counter
  double precision :: dm_crit
  double precision, dimension(Nx,Ny) :: theta

  dm_crit = 1.0d0 !This is the mass of one node
  counter=0


  open(unit=2, file='wall_shear_stress.txt', action='write', status='replace')
  open(unit=3, file='wall_shear_stress_theta.txt', action='write', status='replace')
  open(unit=1, file='circle_theta.txt', action='write', status='replace')
!  do l=2,Ny+1
!    do k=2,Nx+1
!      if(is_solid(k,l).eq.0) then
!        write(2,*) tau_magn(k,l)
!      endif
!    enddo
!  enddo

  do l=1,Ny
    do k=1,Nx
      if((is_solid(k,l).eq.1).and.(is_solid(k,l).eq.-1)) then
        theta(k,l)=0.0
      else
        theta(k,l) = atan2(yy(l)-center_y, xx(k)-center_x)
      endif
    enddo
  enddo

  do k=1, Nx
    write(2,*) (tau_magn(k,l),l=1,Ny)
  enddo

  do l=1, Ny
    do k=1, Nx
      write(3,*) tau_magn(k,l)
      write(1,*) theta(k,l)
    enddo
  enddo

  close(1)
  close(2)
  close(3)
  
 !
 ! If the mass loss at surface is greater than limit
 ! make boundary point into fluid point
 !
 ! do l=2,Ny+1
 !   do k=2,Nx+1
 !     if (dm(k,l).ne.0.0d0) then
 !     endif
 !     if(abs(dm(k,l)).ge.dm_crit) then
 !       is_solid(k,l)=-1
 !       dm_total=dm_total+dm_crit
 !     endif
 !   enddo
 ! enddo

!  open(unit=9, file='mass_loss_er2.txt', status='replace',action='write')
!  write(9,*) dm_total
!  close(9)

 !
 ! we check if any neighbour of a solid point (1)
 ! is a fluid point (-1). The the solid point is turned
 ! into a boundary point (0).
 !
  do l=2,Ny+1
    do k=2,Nx+1
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
            is_solid(k,l)=0
          endif
        enddo
      endif
    enddo
  enddo

 !
 ! Count number of solid and surface points to get the total mass
 !
  do l=2,Ny-1
    do k=2, Nx-1
      if(is_solid(k,l).ne.0) then
        counter=counter+1
      endif
    enddo
  enddo
  !
  !Multiply number of solid and surface points with the mass of one node
  !
  counter=counter*dm_crit

 ! open(unit=8, file='total_mass.txt', action='write', status='replace')
 ! write(8,*) counter
 ! close(8)


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
subroutine mass_loss()

use Evolve
  integer :: q,m,n,k,l,mu,nu
  double precision :: kappa
  double precision, dimension(Nx,Ny) :: tau_er
 
  kappa=1.0d0

  !
  ! Calculate the mass loss as the integral of sigma
  !
    do l=2,Ny+1
      do k=2,Nx+1
        if(is_solid(k,l).eq.0) then
          tau_magn(k,l)=tau_magn(k,l)+(sqrt(forces(k,l,1)**2+forces(k,l,2)**2))
          !tau_magn(k,l)=sqrt(forces(k,l,1)**2+forces(k,l,2)**2)
          !tau_er(k,l)=tau_magn(k,l)/erosion_param
          !dm(k,l)=dm(k,l)+(kappa*(tau_magn(k,l)-tau_er(k,l))) !force into magn and subtract tau_er
          !write(*,*) dm(k,l)
         endif
       enddo
    enddo

endsubroutine mass_loss
!***************************************************************
subroutine eroded_obstacle()

  integer :: m,n,k,l,point_type

  !
  !Start the next run with the eroded surface boundary
  !
  open (unit=99, file='updated_boundary.txt', status='old', action='read')
  do l=2,Ny-1
    do k=2,Nx-1
     read(99, *) point_type
     is_solid(k,l)=point_type
    enddo
  enddo
  close(99)


  open(unit=12, file='fluid_point_initial.txt', action='write', status='replace')
  open(unit=13, file='surface_point_initial.txt', action='write', status='replace')
  open(unit=14, file='solid_point_initial.txt', action='write', status='replace')

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


endsubroutine eroded_obstacle
!***************************************************************
endmodule InitCond


