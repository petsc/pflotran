#ifdef SURFACE_FLOW
module Surface_Factory_module

  use Surface_Simulation_class

  implicit none

  private

#include "definitions.h"

  public :: SurfaceInitialize, &
            SurfaceInitializePostPETSc, &
            HighjackSurfaceSimulation, &
            SurfaceJumpStart

contains

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
subroutine SurfaceInitialize(simulation_base,option)

  use Option_module
  use Input_module
  use Timestepper_Base_class
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(surface_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SurfaceSimulationCreate(option)
  call SurfaceInitializePostPETSc(simulation,option)
  
  simulation_base => simulation

end subroutine SurfaceInitialize

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
subroutine SurfaceInitializePostPETSc(simulation, option)

  use Simulation_module
  use Option_module
  use Init_module
  
  implicit none
  
  class(surface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: simulation_old
  
  !! process command line arguments specific to surface
  !!call SurfInitCommandLineSettings(option)
  
  simulation_old => SimulationCreate(option)
  call Init(simulation_old)
  call HighjackSurfaceSimulation(simulation_old,simulation)
  ! no longer need simulation
  deallocate(simulation_old)
  call SurfaceJumpStart(simulation)
  
end subroutine SurfaceInitializePostPETSc

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/27/13
! ************************************************************************** !
subroutine HighjackSurfaceSimulation(simulation_old,simulation)

  use Simulation_module
  use Surface_Realization_class
  use Option_module
  
  use PMC_Base_class
  use PMC_Surface_class
  use Simulation_Base_class
  use Process_Model_Surface_Flow_class
  use Process_Model_Base_class
  use Process_Model_module
  use Timestepper_Surface_class

  implicit none
  
  type(simulation_type) :: simulation_old
  class(surface_simulation_type) :: simulation
  
  class(pmc_surface_type), pointer :: surf_flow_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(surface_realization_type), pointer :: surf_realization
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  surf_realization => simulation_old%surf_realization
  option => surf_realization%option

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => surf_realization%output_option
  simulation%option => surf_realization%option
  nullify(cur_process_model)

  nullify(surf_flow_process_model_coupler)

  ! Create Surface-flow ProcessModel & ProcessModelCoupler
  if (option%nsurfflowdof > 0) then
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        cur_process_model => PMSurfaceFlowCreate()
      case(TH_MODE)
        !cur_process_model => PMSurfaceTHCreate()
    end select
    cur_process_model%option => surf_realization%option
    cur_process_model%output_option => surf_realization%output_option

    surf_flow_process_model_coupler => PMCSurfaceCreate()
    surf_flow_process_model_coupler%option => option
    surf_flow_process_model_coupler%pm_list => cur_process_model
    surf_flow_process_model_coupler%pm_ptr%ptr => cur_process_model
    call HijackTimestepper(simulation_old%surf_flow_stepper, &
                           surf_flow_process_model_coupler%timestepper)
    nullify(cur_process_model)
  endif

  ! Add the ProcessModelCouplers in a list
  if (associated(surf_flow_process_model_coupler)) then
    simulation%process_model_coupler_list => surf_flow_process_model_coupler%CastToBase()
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoints => surf_realization%waypoints
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pm_list
      do
        if (.not.associated(cur_process_model)) exit
        select type(cur_process_model)
          class is (pm_surface_flow_type)
            call cur_process_model%PMSurfaceFlowSetRealization(surf_realization)
            call cur_process_model_coupler%SetTimestepper( &
                    surf_flow_process_model_coupler%timestepper)
            surf_flow_process_model_coupler%timestepper%dt = option%surf_flow_dt
        end select

        call cur_process_model%Init()
        select type(cur_process_model)
          class is (pm_surface_flow_type)
            select type(ts => cur_process_model_coupler%timestepper)
              class is (timestepper_surface_type)
                call TSSetRHSFunction( &
                          ts%solver%ts, &
                          cur_process_model%residual_vec, &
                          PMRHSFunction, &
                          cur_process_model_coupler%pm_ptr, &
                          ierr)
            end select
        end select
        cur_process_model => cur_process_model%next
      enddo
      cur_process_model_coupler => cur_process_model_coupler%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo

  simulation%surf_realization => surf_realization
  simulation%surf_flow_process_model_coupler => surf_flow_process_model_coupler
  simulation%regression => simulation_old%regression
  surf_flow_process_model_coupler%surf_realization => surf_realization

  ! point the top process model coupler to Output
  !simulation%process_model_coupler_list%Output => Output

end subroutine HighjackSurfaceSimulation

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfaceJumpStart(simulation)

  use Surface_Realization_class
  use Option_module
  use Timestepper_Surface_class
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module  
  use Condition_Control_module
  use Surface_Checkpoint_module
  use Output_Surface_module, only : OutputSurface, OutputSurfaceInit

  implicit none

  type(surface_simulation_type) :: simulation
  
  class(surface_realization_type), pointer :: surf_realization
  class(timestepper_surface_type), pointer :: master_stepper
  class(timestepper_surface_type), pointer :: surf_flow_stepper
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
      surf_flow_stepper => ts
  end select
  nullify(master_stepper)
  
  option => surf_realization%option
  output_option => surf_realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported for surface-flow.'
    call printErrMsg(option)
    return
  endif
  
  master_stepper => surf_flow_stepper

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  surf_flow_read = PETSC_FALSE
  failure = PETSC_FALSE
  
  if (option%restart_flag) then
    call SurfaceRestart(surf_realization,surf_flow_prev_dt,surf_flow_read)

    if(option%time /= option%surf_flow_time) then
      option%io_buffer = 'option%time does not match option%surf_flow_time' // &
        ' while restarting simulation. Check the restart files.'
      call printErrMsg(option)
    endif

    if (surf_flow_read) then
      surf_flow_stepper%prev_dt = surf_flow_prev_dt
      surf_flow_stepper%target_time = option%surf_flow_time
      call TSSetTime(surf_flow_stepper%solver%ts,option%surf_flow_time,ierr)
    endif

  endif

  ! pushed in Init()
  call PetscLogStagePop(ierr)
  option%init_stage = PETSC_FALSE

  ! popped in TimeStepperFinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)

  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_stepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputSurfaceInit(master_stepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_stepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
!    call OutputSurface(surf_realization,plot_flag,transient_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_stepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
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

  if (.not.associated(surf_flow_stepper%cur_waypoint)) then
    option%io_buffer = &
      'Null flow waypoint list; final time likely equal to start time.'
    call printMsg(option)
    return
  else
    surf_flow_stepper%dt_max = surf_flow_stepper%cur_waypoint%dt_max
  endif
          
  surf_flow_stepper%start_time_step = surf_flow_stepper%steps + 1
  
  if (surf_realization%debug%print_couplers) then
  !  call OutputPrintSurfaceCouplers(surf_realization,ZERO_INTEGER)
  endif

end subroutine SurfaceJumpStart

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine HijackTimestepper(stepper_old,stepper_base)

  use Timestepper_Surface_class
  use Timestepper_Base_class
  use Timestepper_module

  implicit none
  
  type(stepper_type), pointer :: stepper_old
  class(stepper_base_type), pointer :: stepper_base
  
  class(timestepper_surface_type), pointer :: stepper
  
  stepper => TimeStepperSurfaceCreate()
  
  stepper%steps = stepper_old%steps
  stepper%num_constant_time_steps = stepper_old%num_constant_time_steps

  stepper%max_time_step = stepper_old%max_time_step
  stepper%max_time_step_cuts = stepper_old%max_time_step_cuts
  stepper%constant_time_step_threshold = stepper_old%constant_time_step_threshold
  stepper%cumulative_time_step_cuts = stepper_old%cumulative_time_step_cuts 
  stepper%cumulative_solver_time = stepper_old%cumulative_solver_time

  stepper%start_time = stepper_old%start_time
  stepper%start_time_step = stepper_old%start_time_step
  stepper%time_step_tolerance = stepper_old%time_step_tolerance
  stepper%target_time = stepper_old%target_time
  
  stepper%prev_dt = stepper_old%prev_dt
!  stepper%dt = stepper_old%dt
  stepper%dt_min = stepper_old%dt_min
  stepper%dt_max = stepper_old%dt_max
  stepper%cfl_limiter = stepper_old%cfl_limiter
  stepper%cfl_limiter_ts = stepper_old%cfl_limiter_ts
  
  stepper%time_step_cut_flag = stepper_old%time_step_cut_flag

  stepper%init_to_steady_state = stepper_old%init_to_steady_state
  stepper%steady_state_rel_tol = stepper_old%steady_state_rel_tol
  stepper%run_as_steady_state = stepper_old%run_as_steady_state

  stepper%solver => stepper_old%solver
  nullify(stepper_old%solver)
  stepper%cur_waypoint => stepper_old%cur_waypoint
  nullify(stepper_old%cur_waypoint)
  
  
!  stepper%prev_waypoint => stepper_old%prev_waypoint
!  nullify(stepper_old%prev_waypoint)
  
!  stepper%revert_dt = stepper_old%revert_dt
!  stepper%num_contig_revert_due_to_sync = &
!  stepper_old%num_contig_revert_due_to_sync

  stepper_base => stepper
  
end subroutine HijackTimestepper

end module Surface_Factory_module
#endif
