module Subsurface_Factory_module

  use Subsurface_Simulation_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  public :: SubsurfaceInitialize, &
            SubsurfaceInitializePostPETSc, &
            HijackSimulation, &
            SubsurfaceJumpStart

contains

! ************************************************************************** !

subroutine SubsurfaceInitialize(simulation_base,option)
  ! 
  ! Sets up PFLOTRAN subsurface simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  use Option_module
  use Input_Aux_module
  use Simulation_Base_class
  
  implicit none
  
  class(simulation_base_type), pointer :: simulation_base
  type(option_type), pointer :: option

  class(subsurface_simulation_type), pointer :: simulation

  ! NOTE: PETSc must already have been initialized here!
  simulation => SubsurfaceSimulationCreate(option)
  call SubsurfaceInitializePostPetsc(simulation,option)
  
  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE
  
  simulation_base => simulation

end subroutine SubsurfaceInitialize

! ************************************************************************** !

subroutine SubsurfaceInitializePostPetsc(simulation, option)
  ! 
  ! Sets up PFLOTRAN subsurface simulation
  ! framework after to PETSc initialization
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/07/13
  ! 

  use Simulation_module
  use Option_module
  use Init_module
  
  implicit none
  
  class(subsurface_simulation_type) :: simulation
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: simulation_old
  
  ! process command line arguments specific to subsurface
  call SubsurfInitCommandLineSettings(option)
  
  simulation_old => SimulationCreate(option)
  call Init(simulation_old)
  call HijackSimulation(simulation_old,simulation)
  call SubsurfaceJumpStart(simulation)
  
  ! no longer need simulation
  ! nullify realization and regression so that it is not destroyed
  nullify(simulation_old%realization)
  nullify(simulation_old%regression)
  call SimulationDestroy(simulation_old)
  
end subroutine SubsurfaceInitializePostPetsc

! ************************************************************************** !

subroutine SubsurfInitCommandLineSettings(option)
  ! 
  ! Initializes PFLTORAN subsurface output
  ! filenames, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif
  
end subroutine SubsurfInitCommandLineSettings

! ************************************************************************** !

subroutine HijackSimulation(simulation_old,simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Simulation_module
  use Realization_class
  use Option_module
  use Output_module, only : Output
  
  use PMC_Base_class
  use PMC_Subsurface_class  
  use PMC_Material_class
  use Simulation_Base_class
  use PM_Base_class
  use PM_General_class
  use PM_Flash2_class
  use PM_Immis_class
  use PM_Mphase_class
  use PM_Miscible_class
  use PM_Richards_class
  use PM_RT_class
  use PM_Subsurface_class
  use PM_TH_class
  use PM_Base_Pointer_module
  use Timestepper_BE_class
  use Logging_module
  use Strata_module
  
  implicit none
  
  type(simulation_type) :: simulation_old
  class(subsurface_simulation_type) :: simulation
  
  class(pmc_subsurface_type), pointer :: flow_process_model_coupler
  class(pmc_subsurface_type), pointer :: tran_process_model_coupler
  class(pmc_material_type), pointer :: material_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(realization_type), pointer :: realization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  realization => simulation_old%realization
  option => realization%option

  simulation%waypoint_list => RealizCreateSyncWaypointList(realization)

  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => realization%output_option
  simulation%option => realization%option
  nullify(cur_process_model)

  nullify(flow_process_model_coupler)
  nullify(tran_process_model_coupler)
  
  ! Create Subsurface-flow ProcessModel & ProcessModelCoupler
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(G_MODE)
        cur_process_model => PMGeneralCreate()
      case(FLASH2_MODE)
        cur_process_model => PMFlash2Create()
      case(IMS_MODE)
        cur_process_model => PMImmisCreate()
      case(MPH_MODE)
        cur_process_model => PMMphaseCreate()
      case(MIS_MODE)
        cur_process_model => PMMiscibleCreate()
      case(RICHARDS_MODE)
        cur_process_model => PMRichardsCreate()
      case(TH_MODE)
        cur_process_model => PMTHCreate()
    end select
    cur_process_model%option => realization%option
    cur_process_model%output_option => realization%output_option

    flow_process_model_coupler => PMCSubsurfaceCreate()
    flow_process_model_coupler%option => option
    flow_process_model_coupler%pms => cur_process_model
    flow_process_model_coupler%pm_ptr%ptr => cur_process_model
!    flow_process_model_coupler%timestepper => simulation_old%flow_timestepper
    flow_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%flow_timestepper, &
                           flow_process_model_coupler%timestepper)
    ! set up logging stage
    string = trim(cur_process_model%name) // 'Flow'
    call LoggingCreateStage(string,flow_process_model_coupler%stage)
    nullify(cur_process_model)
  endif

  ! Create Subsurface transport ProcessModel & ProcessModelCoupler
  if (option%ntrandof > 0) then
    cur_process_model => PMRTCreate()
    cur_process_model%output_option => realization%output_option
    cur_process_model%option => realization%option
   
    tran_process_model_coupler => PMCSubsurfaceCreate()
    tran_process_model_coupler%option => option
    tran_process_model_coupler%pms => cur_process_model
    tran_process_model_coupler%pm_ptr%ptr => cur_process_model
!    tran_process_model_coupler%timestepper => simulation_old%tran_timestepper
    tran_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%tran_timestepper, &
                           tran_process_model_coupler%timestepper)
    ! set up logging stage
    string = 'Reactive Transport'
    call LoggingCreateStage(string,tran_process_model_coupler%stage)
    nullify(cur_process_model)
  endif

  ! Add the ProcessModelCouplers in a list
  if (associated(flow_process_model_coupler)) then
    simulation%process_model_coupler_list => &
      flow_process_model_coupler%CastToBase()
    if (associated(tran_process_model_coupler)) then
      flow_process_model_coupler%below => &
        tran_process_model_coupler%CastToBase()
    endif
  else
    simulation%process_model_coupler_list => &
      tran_process_model_coupler%CastToBase()
  endif
  
  if (StrataEvolves(realization%patch%strata)) then
    material_process_model_coupler => PMCMaterialCreate()
    material_process_model_coupler%option => option
    material_process_model_coupler%realization => realization
    ! place the material process model as %next for the top pmc
    simulation%process_model_coupler_list%next => &
      material_process_model_coupler
  endif

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoint_list => realization%waypoint_list
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pms
      do
        if (.not.associated(cur_process_model)) exit
        ! set realization
        select type(cur_process_model)
          class is (pm_subsurface_type)
            call cur_process_model%PMSubsurfaceSetRealization(realization)
          class is (pm_rt_type)
            call cur_process_model%PMRTSetRealization(realization)
        end select
        ! set time stepper
        select type(cur_process_model)
          class is (pm_rt_type)
            call cur_process_model_coupler%SetTimestepper( &
                   tran_process_model_coupler%timestepper)
            tran_process_model_coupler%timestepper%dt = option%tran_dt
          class default ! otherwise flow
            call cur_process_model_coupler%SetTimestepper( &
                   flow_process_model_coupler%timestepper)
            flow_process_model_coupler%timestepper%dt = option%flow_dt
        end select

        call cur_process_model%Init()
        select type(cur_process_model)
          class default
            select type(ts => cur_process_model_coupler%timestepper)
              class is(timestepper_BE_type)
              call SNESSetFunction( &
                             ts%solver%snes, &
                             cur_process_model%residual_vec, &
                             PMResidual, &
                             cur_process_model_coupler%pm_ptr, &
                                   ierr);CHKERRQ(ierr)
              call SNESSetJacobian( &
                             ts%solver%snes, &
                             ts%solver%J, &
                             ts%solver%Jpre, &
                             PMJacobian, &
                             cur_process_model_coupler%pm_ptr, &
                                   ierr);CHKERRQ(ierr)
            end select
        end select
        cur_process_model => cur_process_model%next
      enddo
      cur_process_model_coupler => cur_process_model_coupler%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo
  
  simulation%realization => realization
  simulation%flow_process_model_coupler => flow_process_model_coupler
  simulation%rt_process_model_coupler => tran_process_model_coupler
  simulation%regression => simulation_old%regression
  
  ! point the top process model coupler to Output
  simulation%process_model_coupler_list%Output => Output

end subroutine HijackSimulation

! ************************************************************************** !

subroutine SubsurfaceJumpStart(simulation)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Realization_class
  use Option_module
  use Timestepper_Base_class
  use Timestepper_BE_class
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Condition_Control_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption  

  implicit none

  type(subsurface_simulation_type) :: simulation
  
  class(realization_type), pointer :: realization
  class(timestepper_base_type), pointer :: master_timestepper
  class(timestepper_BE_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, transient_plot_flag
  PetscBool :: activity_coefs_read
  PetscBool :: flow_read
  PetscBool :: transport_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  realization => simulation%realization
  
  if (associated(simulation%flow_process_model_coupler)) then
    select type(ts => simulation%flow_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        flow_timestepper => ts
    end select
  else
    nullify(flow_timestepper)
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    select type(ts => simulation%rt_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        tran_timestepper => ts
    end select
  else
    nullify(tran_timestepper)
  endif
  nullify(master_timestepper)
  
  option => realization%option
  output_option => realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr);CHKERRQ(ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported in refactored code.'
    call printErrMsg(option)
#if 0    
    call StepperRunSteadyState(realization,flow_timestepper,tran_timestepper)
#endif    
    ! do not want to run through time stepper
    option%status = DONE
    return 
  endif
  
  if (associated(flow_timestepper)) then
    master_timestepper => flow_timestepper
  else
    master_timestepper => tran_timestepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  failure = PETSC_FALSE

#if 0
!geh: moved to within PMInitialize routines
!geh: removed 8/11
  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
    call CondControlAssignFlowInitCond(realization)
  endif

  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
  if (associated(simulation%flow_process_model_coupler)) then
    call simulation%flow_process_model_coupler%UpdateSolution()
  endif
#endif
 
!geh: now performed in PMRTInitializeRun()
!  if (associated(simulation%rt_process_model_coupler)) then
!    call simulation%rt_process_model_coupler%UpdateSolution()
!  endif

  if (option%transport%jumpstart_kinetic_sorption .and. &
      option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. option%restart_flag) then
      option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(option)
    endif
    call RTJumpStartKineticSorption(realization)
  endif
#if 0  
!geh: removed 8/11
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
    option%status = DONE
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
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
    option%status = DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_timestepper)) then
    if (.not.associated(flow_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      flow_timestepper%dt_max = flow_timestepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_timestepper)) then
    if (.not.associated(tran_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      tran_timestepper%dt_max = tran_timestepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_timestepper)) &
    flow_timestepper%start_time_step = flow_timestepper%steps + 1
  if (associated(tran_timestepper)) &
    tran_timestepper%start_time_step = tran_timestepper%steps + 1
  
  if (realization%debug%print_couplers) then
    call OutputPrintCouplers(realization,ZERO_INTEGER)
  endif
#endif
end subroutine SubsurfaceJumpStart

! ************************************************************************** !

subroutine HijackTimestepper(timestepper_old,timestepper_base)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Timestepper_BE_class
  use Timestepper_Base_class
  use Timestepper_module

  implicit none
  
  type(timestepper_type), pointer :: timestepper_old
  class(timestepper_base_type), pointer :: timestepper_base
  
  class(timestepper_BE_type), pointer :: timestepper

  timestepper => TimestepperBECreate()
  
  timestepper%steps = timestepper_old%steps
  timestepper%num_newton_iterations = timestepper_old%num_newton_iterations
  timestepper%num_linear_iterations = timestepper_old%num_linear_iterations
  timestepper%num_constant_time_steps = timestepper_old%num_constant_time_steps

  timestepper%max_time_step = timestepper_old%max_time_step
  timestepper%max_time_step_cuts = timestepper_old%max_time_step_cuts
  timestepper%constant_time_step_threshold = timestepper_old%constant_time_step_threshold
  timestepper%iaccel = timestepper_old%iaccel

  timestepper%cumulative_newton_iterations = timestepper_old%cumulative_newton_iterations
  timestepper%cumulative_linear_iterations = timestepper_old%cumulative_linear_iterations
  timestepper%cumulative_time_step_cuts = timestepper_old%cumulative_time_step_cuts
  timestepper%cumulative_solver_time = timestepper_old%cumulative_solver_time

  timestepper%start_time = timestepper_old%start_time
  timestepper%start_time_step = timestepper_old%start_time_step
  timestepper%time_step_tolerance = timestepper_old%time_step_tolerance
  timestepper%target_time = timestepper_old%target_time
  
  timestepper%prev_dt = timestepper_old%prev_dt
!  stepper%dt = timestepper_old%dt
  timestepper%dt_min = timestepper_old%dt_min
  timestepper%dt_max = timestepper_old%dt_max
  timestepper%cfl_limiter = timestepper_old%cfl_limiter
  timestepper%cfl_limiter_ts = timestepper_old%cfl_limiter_ts
  
  timestepper%time_step_cut_flag = timestepper_old%time_step_cut_flag

  timestepper%ntfac = timestepper_old%ntfac
  timestepper%tfac => timestepper_old%tfac
  nullify(timestepper_old%tfac)
  
  timestepper%init_to_steady_state = timestepper_old%init_to_steady_state
  timestepper%steady_state_rel_tol = timestepper_old%steady_state_rel_tol
  timestepper%run_as_steady_state = timestepper_old%run_as_steady_state

  timestepper%solver => timestepper_old%solver
  nullify(timestepper_old%solver)
  timestepper%convergence_context => timestepper_old%convergence_context
  nullify(timestepper_old%convergence_context)
  timestepper%cur_waypoint => timestepper_old%cur_waypoint
  nullify(timestepper_old%cur_waypoint)
!  stepper%prev_waypoint => timestepper_old%prev_waypoint
!  nullify(timestepper_old%prev_waypoint)
  
!  stepper%revert_dt = timestepper_old%revert_dt
!  stepper%num_contig_revert_due_to_sync = &
!  timestepper_old%num_contig_revert_due_to_sync

  timestepper_base => timestepper
  
end subroutine HijackTimestepper

end module Subsurface_Factory_module
