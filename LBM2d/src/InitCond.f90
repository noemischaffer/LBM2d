!***************************************************************
!
!The convention is that:
!  if a loop goes between 1-Nx+2 (1-Ny+2) the loop index is i (j)
!  if a loop goes between 2-Nx+1 (2-Ny+1) the loop index is k (l)
!
!***************************************************************
module InitCond

use Cdata
use Avg
use messages 
implicit none
private
public:: write_boundary,read_param_init, initff, init_obstacle, construct_surface, update_surface

character (len=labellen) :: init_type='static'
character (len=labellen) :: obstacle_type='circle'
double precision :: radius=1.0d0, center_x=5.0d0, center_y=5.0d0


namelist /init_cond/& 
    init_type, obstacle_type, radius, center_x, center_y

contains
!***************************************************************
subroutine read_param_init(unit,iostat)
!read the namelist init_cond
  integer, intent(in) :: unit
  integer, intent(out) :: iostat
  read(unit, NML=init_cond, IOSTAT=iostat)
endsubroutine read_param_init
!***************************************************************
subroutine initff()
integer :: k,l,q


allocate(refl_point(Nlarge,qmom))


select case (init_type)
case('static')
  do q=1,qmom
     do l=2,Ny+1
        do k=2,Nx+1
           if(is_solid(k,l).ne.1) then
              if(q.eq.5) then
                ff(k,l,q) = 1.0d0
              else
                ff(k,l,q) = 0.0d0
              endif
           endif
        enddo
     enddo
  enddo
case default
        call fatal_error('initialize_force',&
            'error in force initialization! ')
endselect
endsubroutine initff
!***************************************************************
subroutine init_obstacle()

integer :: k, l, q
double precision :: dist_square

select case (obstacle_type)
case('circle')
     do l=2,Ny+1
        do k=2,Nx+1
           if((xx(k-1)-center_x)**2 + (yy(l-1)-center_y)**2.le.radius**2) then
             is_solid(k,l)=1
             !if(is_solid(i,j).eq.1) then
             !write(*,*) is_solid(i,j), i, j
             !endif
             ff(k,l,:) = 0.0d0
             ff(k,l,5)=1.0d0
           endif
        enddo
     enddo
case('rectangle')
case default
        call fatal_error('initialize_obstacle',&
            'error in obstacle initialization! ')
endselect
open(unit=11, file='solid_points.txt', action='write', status='replace')
do l=2, Ny+1
do k=2, Nx+1
if(is_solid(k,l).eq.1) then
write(11,*) k, ',', l
endif
enddo
enddo
close(11)

endsubroutine init_obstacle
!***************************************************************
subroutine construct_surface()



allocate(surface(Nlarge,2))

call update_surface()

refl_point(Nsurf,qmom)=0

endsubroutine construct_surface
!***************************************************************
subroutine update_surface()

integer ::q,kk,m,n, k, l, l_up, l_down, k_left, k_right
Nsurf=0

open(unit=12, file='refl_points.txt', action='write', status='replace')
!open(unit=15, file='refl_points_q.txt', action='write', status='replace')


do l=2,Ny+1
   do k=2,Nx+1
      if(is_solid(k,l).eq.1) then
        kk=0 
        do q=1,qmom
           m=k-ee_int(1,q)
           n=l-ee_int(2,q)
           if(is_solid(m,n).ne.1) then
             !write(15,*) q
             if(kk.eq.0) then            
               Nsurf=Nsurf+1
               surface(Nsurf,1)=k
               surface(Nsurf,2)=l
               kk=kk+1
             endif
             refl_point(Nsurf,q)=1
             write(12,*) k,',', l
             !write(15,*) 'compare with'
             !write(15,*) q
           endif
         enddo
      endif
   enddo
enddo
!close(15)
close(12)
call write_boundary()

endsubroutine update_surface
!***************************************************************
subroutine write_boundary()

integer::jsurf,q

open(unit=11, file='surface.txt', action='write', status='replace')
!do q=1, qmom
do jsurf=1, Nsurf
write(11,*) surface(jsurf,1),',', surface(jsurf,2)
!write(*,*) refl_point(jsurf,q)
!write(*,*) refl_point(jsurf,q)
enddo
!enddo
close(11)

endsubroutine write_boundary

endmodule InitCond
