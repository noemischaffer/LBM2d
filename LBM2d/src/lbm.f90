program lbm
!---------------------!
use Cdata
use IO  
use Avg
use Diagnostic
use Evolve
use Force
use Messages
use Sub
use BoundCond
use InitCond
use ShearStress
!---------------------!
implicit none
integer :: it,jouter,iin,k,l, dt
character (len=25) :: filename
character (len=25) :: filenama
double precision :: tflow
!---------------------!
write(*,*) '***************************************'
write(*,*) '    STARTING LBM2d                     '
write(*,*) '***************************************'
write(*,*) '========= Reading input data ...======='
call read_input()
write(*,*) ' done....'
call allocate_cdata()
call write_input()
write(*,*) '============================================'
write(*,*) '                Model D2Q9                  '
write(*,*) '========Setting model parameters...========='
call set_model()
call write_model_param()
write(*,*) '.....done '
write(*,*) '============================================'
write(*,*) '============Allocating large arrays========='
call allocate_avg()
call allocate_shear()
call allocate_force()
write(*,*) '..done.'
write(*,*) '============================================'
!
! If we start from scratch then we implement initial condition
! otherwise we read the data from a file.
!
if (lstart) then
   write(*,*) '==== Initializing .....====================='
   write(*,*) ' First we set all points as fluid points ...'
   call initff()
   write(*,*) ' Then we initialize boundary condition .... '
   call initialize_bc()
   write(*,*) ' Then we set the immersed objects .....     '
   call init_obstacle()
   call construct_surface()
   write(*,*) 'and write down the is_solid array ...        '
   call write_immersed_boundary(0)
   call boundary_condition()
   call update_surface()
  !Save the original surface points
  open(unit=13, file='surface_pointx_original_eroded.dat', action='write', status='replace')
  do l=2,Ny-1
    do k=2, Nx-1
      if(is_solid(k,l).eq.0) then
        write(13,*) k
      endif
    enddo
  enddo
  open(unit=12, file='surface_pointy_original_eroded.dat', action='write', status='replace')
    do l=2,Ny-1
      do k=2, Nx-1
        if(is_solid(k,l).eq.0) then
          write(12,*) l
        endif
      enddo
   enddo
   call initialize_diag()
   call calc_avg()
   call write_ts(0)
   call write_snap(0)
   write(*,*) '... Initial set up is done.                  '
else 
   call fatal_error('LBM', 'restarting runs not set up')
   write(*,*) '============================================'
endif
write(*,*) '=== Starting time stepping ================='
!do it=1,iTMAX
it=0
tflow=radius/v_const
do while (it .le. tflow)
  it=it+1
  call boundary_condition()
  call stream()
  call calc_avg()
  call comp_equilibrium_BGK()
  call collision()
  call wall_shear()
  !call update_surface()
  call mass_loss()
  call calc_diag()
  call write_ts(it)
  if (mod(it,1) .eq. 0) then
    WRITE (filename, '(a,I5.5,a)') 'velocitx_series',it,'.dat'
    open(unit=10, file=filename, action='write', status='replace')
     do k=1, Nx+2
       write(10,*) (uu(k,l,1), l=1,Ny+2)
     enddo
    close(10)
    WRITE (filenama, '(a,I5.5,a)') 'velocity_series',it,'.dat'
    open(unit=11, file=filenama, action='write', status='replace')
     do k=1, Nx+2
       write(11,*) (uu(k,l,2), l=1,Ny+2)
     enddo
    close(11)
  endif
enddo
call update_surface()
!Save surface boundary for next run
open(unit=19, file='updated_boundary.txt', action='write', status='replace')
  do l=2,Ny-1
    do k=2, Nx-1
      write(19,*) is_solid(k,l)
    enddo
  enddo
close(19)
write(*,*) '......done'
write(*,*) '===writing the final snapshot ============'
!call write_snap(iTMAX)
call rwrite_density_uu()
!call wall_shear()
call theo_solution()
write(*,*) '=========================================='
call free_avg()
call free_shear()
call free_force()
call free_cdata()
endprogram lbm
!----------------------------------!
