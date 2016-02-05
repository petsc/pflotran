module Factory_Surface_module

  use Simulation_Surface_class

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: SurfaceInitialize, &
            SurfaceInitializePostPETSc, &
!            HijackSurfaceSimulation, &
            SurfaceJumpStart

contains

! ************************************************************************** !

subroutine SurfaceInitialize(simulation_base,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  use Input_Aux_module
  use Timestepper_Base_class
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(simulation_surface_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SurfaceSimulationCreate(option)
  call SurfaceInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine SurfaceInitialize

! ************************************************************************** !

subroutine SurfaceInitializePostPETSc(simulation, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  use Init_Common_module
  
  implicit none
  
  class(simulation_surface_type) :: simulation
  type(option_type), pointer :: option
  
  call SurfaceJumpStart(simulation)

end subroutine SurfaceInitializePostPETSc

! ************************************************************************** !

subroutine SurfaceJumpStart(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Realization_Surface_class
  use Option_module
  use Timestepper_Surface_class
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module  
  use Condition_Control_module
  use Checkpoint_Surface_module
  use Output_Surface_module, only : OutputSurface, OutputSurfaceInit

  implicit none

  type(simulation_surface_type) :: simulation
  
  class(realization_surface_type), pointer :: surf_realization
  class(timestepper_surface_type), pointer :: master_timestepper
  class(timestepper_surface_type), pointer :: surf_flow_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: surf_flow_read
  PetscBool :: failure
  PetscErrorCode :: ierr
  PetscReal :: surf_flow_prev_dt

  surf_realization => simulation%surf_realization
  
  select type(ts => simulation%surf_flow_process_model_coupler%timestepper)
    class is(timestepper_surface_type)
      surf_flow_timestepper => ts
  end select
  nullify(master_timestepper)
  
  option => surf_realization%option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  master_timestepper => surf_flow_timestepper

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  surf_flow_read = PETSC_FALSE
  failure = PETSC_FALSE
  
#if 0
  if (option%restart_flag) then
    call SurfaceRestart(surf_realization,surf_flow_prev_dt,surf_flow_read)

    if (option%time /= option%surf_flow_time) then
      option%io_buffer = 'option%time does not match option%surf_flow_time' // &
        ' while restarting simulation. Check the restart files.'
      call printErrMsg(option)
    endif

    if (surf_flow_read) then
      surf_flow_timestepper%prev_dt = surf_flow_prev_dt
      surf_flow_timestepper%target_time = option%surf_flow_time
      call TSSetTime(surf_flow_timestepper%solver%ts,option%surf_flow_time, &
                     ierr);CHKERRQ(ierr)
    endif

  endif
#endif

  ! pushed in Init()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in TimestepperFinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_timestepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputSurfaceInit(master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
!    call OutputSurface(surf_realization,plot_flag,transient_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (.not.associated(surf_flow_timestepper%cur_waypoint)) then
    option%io_buffer = &
      'Null flow waypoint list; final time likely equal to start time.'
    call printMsg(option)
    return
  else
    surf_flow_timestepper%dt_max = surf_flow_timestepper%cur_waypoint%dt_max
  endif
          
  surf_flow_timestepper%start_time_step = surf_flow_timestepper%steps + 1
  
  if (surf_realization%debug%print_couplers) then
  !  call OutputPrintSurfaceCouplers(surf_realization,ZERO_INTEGER)
  endif

end subroutine SurfaceJumpStart

! ************************************************************************** !

end module Factory_Surface_module
