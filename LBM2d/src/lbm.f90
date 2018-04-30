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
integer :: it,jouter,iin,k,l
character (len=25) :: filename
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
!   call wall_shear()
   write(*,*) 'and write down the is_solid array ...        '
   call write_immersed_boundary(0)
   call boundary_condition()
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
do it=1,iTMAX
  call boundary_condition()
  call stream()
  call calc_avg()
  call comp_equilibrium_BGK()
  call collision()
  call update_surface()
!  call wall_shear()
  call calc_diag()
  call write_ts(it)
  if (mod(it,1) .eq. 0) then
    WRITE (filename, '(a,I5.5,a)') 'density_series',it,'.dat'
    open(unit=10, file=filename, action='write', status='replace')
     do k=1, Nx+2
       write(10,*) (ff(k,l,3), l=1,Ny+2)
     enddo
    close(10)
  endif
enddo
write(*,*) '......done'
write(*,*) '===writing the final snapshot ============'
call write_snap(iTMAX)
call rwrite_density_uu()
write(*,*) '=========================================='
call free_avg()
call free_cdata()
endprogram lbm
!----------------------------------!
