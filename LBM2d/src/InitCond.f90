module InitCond

use Cdata
use Avg
use messages 
implicit none
private
public:: read_param_init, initff, init_obstacle, construct_interface

character (len=labellen) :: init_type='static'
character (len=labellen) :: obstacle_type='circle'
integer :: radius=4, center_x=5, center_y=5
integer, allocatable, dimension(:,:) :: interface_coord
logical::linterf=.false.

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
integer :: i,j,q

select case (init_type)
case('static')
  do q=1,qmom
     do j=2,Ny+1
        do i=2,Nx+1
           if(is_solid(i,j).ne.1) then
              if(q.eq.5) then
                ff(i,j,q) = 1.0d0
              else
                ff(i,j,q) = 0.0d0
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

integer :: dist_square, i, j, q

select case (obstacle_type)
case('circle')
  do q=1,qmom
     do j=2,Ny+1
        do i=2,Nx+1
           dist_square=(center_x-i)**2 + (center_y-j)**2
           if(dist_square.le.radius**2) then
              is_solid(i,j)=1
           endif
           if(is_solid(i,j).eq.1) then
                ff(i,j,q) = 0.0d0
           endif
        enddo
     enddo
  enddo
case default
        call fatal_error('initialize_obstacle',&
            'error in obstacle initialization! ')
endselect
endsubroutine init_obstacle
!***************************************************************
subroutine construct_interface()

integer :: Ninterface=0, i, j, j_up, j_down, i_left, i_right, kint=1
do j=1,Ny
   do i=1,Nx
      if(is_solid(i,j).eq.1) then
         j_up=j+1
         j_down=j-1
         i_left=i-1
         i_right=i+1
         if((is_solid(i_left,j_up).eq.1).and.(is_solid(i_left,j_down).eq.1) &
            .and.(is_solid(i_right,j_up).eq.1).and.(is_solid(i_right,j_down).eq.1)) then
            continue
         else
            Ninterface=Ninterface+1
         endif
      endif
   enddo
enddo

allocate(interface_coord(2,Ninterface))

do j=1,Ny
   do i=1,Nx
      if(is_solid(i,j).eq.1) then
         j_up=j+1
         j_down=j-1
         i_left=i-1
         i_right=i+1
         if((is_solid(i_left,j_up).eq.1).and.(is_solid(i_left,j_down).eq.1) &
            .and.(is_solid(i_right,j_up).eq.1).and.(is_solid(i_right,j_down).eq.1)) then
            continue
         else
            interface_coord(1,kint)=i
            interface_coord(2,kint)=j
            kint=kint+1
         endif
      endif
   enddo
enddo

!      write(*,*) 'x coordinates of interface', interface_coord(1,:)
!      write(*,*) 'y coordinates of interface', interface_coord(2,:)

if (linterf .eqv. .true.) then
     deallocate(interface_coord)
endif

endsubroutine construct_interface
!***************************************************************
endmodule InitCond
