!=======================================================================
! PFLOTRAN v1.0 LA-CC-06-093
!=======================================================================

! The software titled "PFLOTRAN v 1.0" has been assigned LA-CC-06-093. 
! The software is unclassified and does not contain Unclassified 
! Controlled Nuclear Information (UCNI).

! The software is under review to be released as open source, which 
! requires review and approval by the appropriate DOE Program Office. 
! If DOE declines to release this software as publicly available open 
! source, it would be subject to export control under Department of 
! Commerce regulations, and classified as ECCN 
! (Export Control Classification Number) EAR99. 
 
! If released as open source, the software will be publicly available and 
! not subject to export control. Until DOE approval is obtained, the 
! software should be treated as export controlled. DOE reserves the 
! right to release or deny release of the software as open source.

! Send all bug reports/questions/comments to:
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-6, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM
!=======================================================================
  program pflow
  
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Solver_module
  use Grid_module
  use Option_module
  use Init_module
  use Output_module  
  
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petsclog.h"

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  integer :: myrank, commsize

  integer :: ierr
  PetscInt :: stage(10)
  PetscTruth :: restartflag
  character(len=MAXSTRINGLENGTH) :: restartfile
  PetscInt :: chkptfreq  
    ! The number of flow steps between checkpoints.
  PetscTruth :: chkptflag
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflowin

  
  type(simulation_type), pointer :: simulation
  type(solver_type), pointer :: solver
  type(stepper_type), pointer :: stepper
  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)

  simulation => SimulationCreate()
  realization => simulation%realization
  option => realization%option
  stepper => simulation%stepper

  option%myrank = myrank
  option%commsize = commsize

! Register stages for profiling.
! We identify three stages here: setup, output, and cleanup.
! Note that we do NOT identify a separate stage for the calculations -- 
! we include those as part of the default "Main" stage (so we need 
! to pop all other stages off when calculations are proceeding).
  call PetscLogStageRegister(stage(1), "PFLOW setup", ierr)
  call PetscLogStageRegister(stage(2), "PFLOW output", ierr)
  call PetscLogStageRegister(stage(3), "PFLOW cleanup", ierr)

! Starting PFLOW setup -- push this onto the log stack.
  call PetscLogStagePush(stage(1), ierr)
  
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflowin", &
                             pflowin, option_found, ierr)
  if(option_found /= PETSC_TRUE) pflowin = "pflow.in"
  
  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)

  call OptionCheckCommandLine(option)

  call PflowInit(simulation,pflowin)

! Done with PFLOW setup.
  call PetscLogStagePop(ierr)

! Now we do some initial output, so push this onto the log stack.
  call PetscLogStagePush(stage(2), ierr)

  if (option%restart_flag == PETSC_FALSE) then
    if(option%imode /= OWG_MODE) then
      call Output(realization)
    else
   !   call pflow_var_output(grid,kplt,iplot)
    endif
  endif
                          
  call PetscLogStagePop(ierr)

#if 0
 ! still needs to be implemented                    
 if (grid%iread_init==2 .and. option%restartflag == PETSC_FALSE)then
      call Read_init_field(grid, kplt)
      call initAccumulation(realization)
      print *, 'Restart from ASCII file: pflow_init.dat'
  endif
#endif  
           
  call StepperRun(realization,stepper,stage)
  
! Clean things up.
  call PetscLogStagePush(stage(3), ierr)
  if (option%ibrkcrv > 0) close(IUNIT4)
  call SimulationDestroy(simulation)
  call PetscLogStagePop(ierr)

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  call PetscGetTime(timex_wall(2), ierr)
  
  if (myrank == 0) then
  
!       write(*,'("steps: ",13x,i6,/, &
!   &   "total newton iters: ",i6,/, &
!   &   "total PETSc gmres:  ",i6,/, &
!   &   "total cuts: ",11x,i3)') steps,newtcum,itpetscum,icutcum
        
        write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0

        write(IUNIT2,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

        write(IUNIT2,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    &   1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0
  endif

  close(IUNIT2)

  call PetscFinalize (ierr)

  end program pflow
