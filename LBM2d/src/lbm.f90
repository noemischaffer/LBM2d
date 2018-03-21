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
integer :: it,jouter,iin
character (len=25) :: filename
!---------------------!
write(*,*) '***************************************'
write(*,*) '    STARTING LBM2d                     '
write(*,*) '***************************************'
write(*,*) '========= Reading input data ...======='
call read_input()
write(*,*) ' done....'
call write_input()
write(*,*) '============================================'
write(*,*) '                Model D2Q9                  '
write(*,*) '========Setting model parameters...========='
call set_model()
call write_model_param()
write(*,*) '.....done '
write(*,*) '============================================'
write(*,*) '============Allocating large arrays========='
call allocate_cdata()
call allocate_avg()
write(*,*) '..done.'
write(*,*) '============================================'
!call allocate_shear()
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
   write(*,*) 'and write down the is_solid array ...        '
   call write_immersed_boundary(0)
   call boundary_condition()
   write(*,*) '... Initial set up is done.                  '
else 
   call fatal_error('LBM', 'restarting runs not set up')
   write(*,*) '============================================'
endif
write(*,*) '=== Starting time stepping ================='
call initialize_diag()
jouter=iTMAX/ndiag
do it=1,jouter
   do iin=1,ndiag
     call calc_avg()
     call stream()
     !  call vorticity()
     !  call rwrite_density_uu()
     call comp_equilibrium_BGK()
     call collision()
     !  call wall_shear()
     !  call update_surface()
     call boundary_condition()
  enddo
  !call calc_diag()
  call calc_diag()
  call write_ts(it*ndiag)
enddo
write(*,*) '......done'
write(*,*) '===writing the final snapshot ============'
call write_snap(iTMAX)
write(*,*) '=========================================='
!call free_shear()
call rwrite_density_uu()
call free_avg()
call free_cdata()
endprogram lbm
!----------------------------------!
