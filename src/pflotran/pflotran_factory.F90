module PFLOTRAN_Factory_module

  use Subsurface_Simulation_class
  use PMC_Subsurface_class
  
  implicit none

  private

#include "definitions.h"

  public :: PFLOTRANInitialize, &
            PFLOTRANInitializePrePETSc, &
            PFLOTRANInitializePostPETSc, &
            PFLOTRANRun, &
            PFLOTRANFinalize

contains

! ************************************************************************** !
!
! PFLOTRANInitialize: Sets up PFLOTRAN subsurface simulation 
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PFLOTRANInitialize(simulation,option)

  use Option_module
  use Input_module
  use Timestepper_module
  use Logging_module
  
  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  type(option_type), pointer :: option

  call PFLOTRANInitializePrePETSc(option)
  call OptionInitPetsc(option)
  call LoggingCreate()
  call PFLOTRANInitializePostPETSc(simulation,option)

end subroutine PFLOTRANInitialize

! ************************************************************************** !
!
! PFLOTRANInitializePrePETSc: Sets up PFLOTRAN subsurface simulation 
!                             framework prior to PETSc initialization
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANInitializePrePETSc(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  
  ! NOTE: Cannot add anything that requires PETSc in this routins as PETSc 
  !       has not yet been initialized.
  
  call PFLOTRANInitCommandLineSettings(option)
  
end subroutine PFLOTRANInitializePrePETSc

! ************************************************************************** !
!
! PFLOTRANInitializePostPETSc: Sets up PFLOTRAN subsurface simulation 
!                              framework after to PETSc initialization
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANInitializePostPETSc(simulation, option)

  use Simulation_module
  use Timestepper_module
  use Option_module
  use Init_module
  
  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: simulation_old
  PetscInt :: init_status
  
  call OptionBeginTiming(option)
  simulation_old => SimulationCreate(option)
  simulation => SubsurfaceSimulationCreate(option)
  call Init(simulation_old)
  call HighjackSimulation(simulation_old,simulation)
  ! no longer need simulation
  deallocate(simulation_old)
  call SubsurfaceJumpStart(simulation,init_status)
  
end subroutine PFLOTRANInitializePostPETSc

! ************************************************************************** !
!
! PFLOTRANRun: Runs the PFLOTRAN simulation
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANRun(simulation, master_stepper, init_status)

  use Simulation_module
  use Timestepper_module
  
  implicit none
  
  type(simulation_type) :: simulation
  type(stepper_type), pointer :: master_stepper
  PetscInt :: init_status

end subroutine PFLOTRANRun

! ************************************************************************** !
!
! PFLOTRANFinalize: Destroys PFLOTRAN subsurface simulation framework
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANFinalize(simulation,option)

  use Regression_module
  use Option_module
  
  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  type(option_type) :: option
  
  call simulation%FinalizeRun()
  call SubsurfaceSimulationDestroy(simulation)
  call OptionEndTiming(option)
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif

end subroutine PFLOTRANFinalize

! ************************************************************************** !
!
! PFLOTRANInitCommandLineSettings: Initializes PFLOTRAN output filenames, etc.
! author: Glenn Hammond
! date: 06/06/13
!
! ************************************************************************** !
subroutine PFLOTRANInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  PetscInt :: i
  PetscErrorCode :: ierr
  
  ! check for non-default input filename
  option%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename, &
                                 pflotranin_option_found,option)
  string = '-input_prefix'
  call InputGetCommandLineString(string,option%input_prefix, &
                                 input_prefix_option_found,option)
  
  if (pflotranin_option_found .and. input_prefix_option_found) then
    option%io_buffer = 'Cannot specify both "-pflotranin" and ' // &
      '"-input_prefix" on the command lines.'
    call printErrMsg(option)
  else if (pflotranin_option_found) then
    !TODO(geh): replace this with StringSplit()
    i = index(option%input_filename,'.',PETSC_TRUE)
    if (i > 1) then
      i = i-1
    else
      ! for some reason len_trim doesn't work on MS Visual Studio in 
      ! this location
      i = len(trim(option%input_filename)) 
    endif
    option%input_prefix = option%input_filename(1:i)
  else if (input_prefix_option_found) then
    option%input_filename = trim(option%input_prefix) // '.in'
  endif
  
  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)
  if (.not.option_found) option%global_prefix = option%input_prefix  
  
  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-v'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%verbosity = 1
 
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%simulation_type = MULTISIMULATION_SIM_TYPE

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%simulation_type = STOCHASTIC_SIM_TYPE

  ! this will get overwritten later if stochastic
  string = '-realization_id'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%id = i

end subroutine PFLOTRANInitCommandLineSettings

! ************************************************************************** !
!
! HighjackSimulation: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HighjackSimulation(simulation_old,simulation)

  use Simulation_module
  use Realization_class
  use Option_module
  
  use PMC_Base_class
  use Simulation_Base_class
  use Process_Model_Richards_class
  use Process_Model_RT_class
  use Process_Model_TH_class
  use Process_Model_THC_class
  use Process_Model_Base_class
  use Process_Model_module

  implicit none
  
  type(simulation_type) :: simulation_old
  class(subsurface_simulation_type) :: simulation
  
  class(pmc_subsurface_type), pointer :: flow_process_model_coupler
  class(pmc_subsurface_type), pointer :: tran_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pm_base_type), pointer :: cur_process_model
  
  class(realization_type), pointer :: realization
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  realization => simulation_old%realization
  option => realization%option

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
      case(RICHARDS_MODE)
        cur_process_model => PMRichardsCreate()
      case(TH_MODE)
        cur_process_model => PMTHCreate()
      case(THC_MODE)
        cur_process_model => PMTHCCreate()
    end select
    cur_process_model%option => realization%option
    cur_process_model%output_option => realization%output_option

    flow_process_model_coupler => PMCSubsurfaceCreate()
    flow_process_model_coupler%option => option
    flow_process_model_coupler%pm_list => cur_process_model
    flow_process_model_coupler%pm_ptr%ptr => cur_process_model
    flow_process_model_coupler%timestepper => simulation_old%flow_stepper
    nullify(cur_process_model)
  endif

  ! Create Subsurface transport ProcessModel & ProcessModelCoupler
  if (option%ntrandof > 0) then
    cur_process_model => PMRTCreate()
    cur_process_model%output_option => realization%output_option
    cur_process_model%option => realization%option
   
    tran_process_model_coupler => PMCSubsurfaceCreate()
    tran_process_model_coupler%option => option
    tran_process_model_coupler%pm_list => cur_process_model
    tran_process_model_coupler%pm_ptr%ptr => cur_process_model
    tran_process_model_coupler%timestepper => simulation_old%tran_stepper
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

  ! For each ProcessModel, set:
  ! - realization (subsurface or surface),
  ! - stepper (flow/trans/surf_flow),
  ! - SNES functions (Residual/Jacobain), or TS function (RHSFunction)
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler_top%waypoints => realization%waypoints
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%pm_list
      do
        if (.not.associated(cur_process_model)) exit
        select type(cur_process_model)
          class is (pm_richards_type)
            call cur_process_model%PMRichardsSetRealization(realization)
            call cur_process_model_coupler%SetTimestepper(simulation_old%flow_stepper)
            simulation_old%flow_stepper%dt = option%flow_dt
          class is (pm_rt_type)
            call cur_process_model%PMRTSetRealization(realization)
            call cur_process_model_coupler%SetTimestepper(simulation_old%tran_stepper)
            simulation_old%tran_stepper%dt = option%tran_dt
          class is (pm_th_type)
            call cur_process_model%PMTHSetRealization(realization)
            call cur_process_model_coupler%SetTimestepper(simulation_old%flow_stepper)
            simulation_old%flow_stepper%dt = option%flow_dt
          class is (pm_thc_type)
            call cur_process_model%PMTHCSetRealization(realization)
            call cur_process_model_coupler%SetTimestepper(simulation_old%flow_stepper)
            simulation_old%flow_stepper%dt = option%flow_dt
        end select

        call cur_process_model%Init()
        select type(cur_process_model)
          class default
            call SNESSetFunction( &
                           cur_process_model_coupler%timestepper%solver%snes, &
                           cur_process_model%residual_vec, &
                           PMResidual, &
                           cur_process_model_coupler%pm_ptr,ierr)
            call SNESSetJacobian( &
                           cur_process_model_coupler%timestepper%solver%snes, &
                           cur_process_model_coupler%timestepper%solver%J, &
                           cur_process_model_coupler%timestepper%solver%Jpre, &
                           PMJacobian, &
                           cur_process_model_coupler%pm_ptr,ierr)
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

end subroutine HighjackSimulation

! ************************************************************************** !
!
! SubsurfaceJumpStart: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SubsurfaceJumpStart(simulation, init_status)

  use Realization_class
  use Option_module
  use Timestepper_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputPrintCouplers
  use Logging_module  
  use Condition_Control_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption  

  implicit none

  type(subsurface_simulation_type) :: simulation
  PetscInt :: init_status
  
  class(realization_type), pointer :: realization
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
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
    flow_stepper => simulation%flow_process_model_coupler%timestepper
  else
    nullify(flow_stepper)
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    tran_stepper => simulation%rt_process_model_coupler%timestepper
  else
    nullify(tran_stepper)
  endif
  nullify(master_stepper)
  
  option => realization%option
  output_option => realization%output_option

  init_status = TIMESTEPPER_INIT_PROCEED

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported in refactored code.'
    call printErrMsg(option)
#if 0    
    call StepperRunSteadyState(realization,flow_stepper,tran_stepper)
#endif    
    ! do not want to run through time stepper
    init_status = TIMESTEPPER_INIT_DONE
    return 
  endif
  
  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
  else
    master_stepper => tran_stepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  failure = PETSC_FALSE
  
  if (option%restart_flag) then
    call SubsurfaceRestart(realization,flow_stepper,tran_stepper, &
                           flow_read,transport_read,activity_coefs_read)
#if 0      
  else if (master_stepper%init_to_steady_state) then
    option%print_screen_flag = OptionPrintToScreen(option)
    option%print_file_flag = OptionPrintToFile(option)
    if (associated(flow_stepper)) then
      if (flow_stepper%init_to_steady_state) then
        option%flow_dt = master_stepper%dt_min
        call FlowStepperStepToSteadyState(realization,flow_stepper,failure)
        if (failure) then ! if flow solve fails, exit
          if (OptionPrintToScreen(option)) then
            write(*,*) ' ERROR: steady state solve failed!!!'
          endif
          if (OptionPrintToFile(option)) then
            write(option%fid_out,*) ' ERROR: steady state solve failed!!!'
          endif
          init_status = TIMESTEPPER_INIT_FAIL
          return 
        endif
        option%flow_dt = master_stepper%dt_min
      endif
      if (flow_stepper%run_as_steady_state) then
        master_stepper => tran_stepper
      endif
    endif

    if (associated(tran_stepper)) then
      if (tran_stepper%init_to_steady_state) then
    ! not yet functional
    !    step_to_steady_state = PETSC_TRUE
    !    option%tran_dt = master_stepper%dt_min
    !    call StepperStepTransportDT(realization,tran_stepper, &
    !                                tran_timestep_cut_flag, &
    !                                idum,step_to_steady_state,failure)
    !    if (failure) return ! if flow solve fails, exit
    !    option%tran_dt = master_stepper%dt_min
      endif
    endif
! #if 0
#endif
  endif

  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif

  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
  if (associated(simulation%flow_process_model_coupler)) then
    call simulation%flow_process_model_coupler%UpdateSolution()
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    call simulation%rt_process_model_coupler%UpdateSolution()
  endif

  if (option%jumpstart_kinetic_sorption .and. option%time < 1.d-40) then
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
    init_status = TIMESTEPPER_INIT_DONE
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(master_stepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_stepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
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
    init_status = TIMESTEPPER_INIT_DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      init_status = TIMESTEPPER_INIT_FAIL
      return
    else
      flow_stepper%dt_max = flow_stepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_stepper)) then
    if (.not.associated(tran_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      init_status = TIMESTEPPER_INIT_FAIL
      return
    else
      tran_stepper%dt_max = tran_stepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_stepper)) &
    flow_stepper%start_time_step = flow_stepper%steps + 1
  if (associated(tran_stepper)) &
    tran_stepper%start_time_step = tran_stepper%steps + 1
  
  if (realization%debug%print_couplers) then
    call OutputPrintCouplers(realization,ZERO_INTEGER)
  endif

end subroutine SubsurfaceJumpStart

! ************************************************************************** !
!
! SubsurfaceRestart: Calls appropriate routines to read checkpoint file and 
!                    restart
! author: Glenn Hammond
! date: 03/07/08; 06/11/13
!
! ************************************************************************** !
subroutine SubsurfaceRestart(realization,flow_stepper,tran_stepper, &
                             flow_read,transport_read,activity_coefs_read)

  use Realization_class
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Checkpoint_module
  use Option_module
  use Timestepper_module
  use Waypoint_module
  
  use Flash2_module, only: Flash2UpdateAuxVars
  use Mphase_module, only: MphaseUpdateAuxVars
  use Immis_module, only: ImmisUpdateAuxVars
  use Miscible_module, only: MiscibleUpdateAuxVars
  use Richards_module, only : RichardsUpdateAuxVars
  use TH_module, only : THUpdateAuxVars
  use THC_module, only : THCUpdateAuxVars
  use THMC_module, only : THMCUpdateAuxVars
  use General_module, only : GeneralUpdateAuxVars  

  implicit none

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  PetscBool :: activity_coefs_read
  PetscBool :: flow_read
  PetscBool :: transport_read

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_cumulative_newton_iterations, &
              flow_cumulative_time_step_cuts, flow_cumulative_linear_iterations ,&
              flow_num_constant_time_steps, flow_num_newton_iterations
  PetscInt :: tran_steps, tran_cumulative_newton_iterations,  &
              tran_cumulative_time_step_cuts, tran_cumulative_linear_iterations, &
              tran_num_constant_time_steps, tran_num_newton_iterations
  PetscReal :: flow_cum_solver_time, flow_prev_dt
  PetscReal :: tran_cum_solver_time, tran_prev_dt
  
  option => realization%option

  call Restart(realization, &
               flow_steps,flow_cumulative_newton_iterations, &
               flow_cumulative_time_step_cuts,flow_cumulative_linear_iterations, &
               flow_num_constant_time_steps,flow_num_newton_iterations, &
               flow_cum_solver_time,flow_prev_dt, &
               tran_steps,tran_cumulative_newton_iterations, &
               tran_cumulative_time_step_cuts,tran_cumulative_linear_iterations, &
               tran_num_constant_time_steps,tran_num_newton_iterations, &
               tran_cum_solver_time,tran_prev_dt, &
               flow_read,transport_read,activity_coefs_read)
  if (option%restart_time < -998.d0) then
    option%time = max(option%flow_time,option%tran_time)
    if (associated(flow_stepper) .and. flow_read) then
      flow_stepper%steps = flow_steps
      flow_stepper%cumulative_newton_iterations = flow_cumulative_newton_iterations
      flow_stepper%cumulative_time_step_cuts = flow_cumulative_time_step_cuts
      flow_stepper%cumulative_linear_iterations = flow_cumulative_linear_iterations
      flow_stepper%cumulative_solver_time = flow_cum_solver_time
      flow_stepper%prev_dt = flow_prev_dt
      if (.not.flow_stepper%run_as_steady_state) then
        flow_stepper%num_constant_time_steps = flow_num_constant_time_steps
        flow_stepper%num_newton_iterations = flow_num_newton_iterations
      endif
    endif
    if (associated(tran_stepper) .and. transport_read) then
      tran_stepper%steps = tran_steps
      tran_stepper%cumulative_newton_iterations = tran_cumulative_newton_iterations
      tran_stepper%cumulative_time_step_cuts = tran_cumulative_time_step_cuts
      tran_stepper%cumulative_linear_iterations = tran_cumulative_linear_iterations
      tran_stepper%cumulative_solver_time = tran_cum_solver_time
      tran_stepper%num_constant_time_steps = tran_num_constant_time_steps
      tran_stepper%num_newton_iterations = tran_num_newton_iterations
      tran_stepper%prev_dt = tran_prev_dt
    endif
  else
    option%time = option%restart_time
    option%flow_time = option%restart_time
    option%tran_time = option%restart_time
    if (associated(flow_stepper)) then
      option%flow_dt = flow_stepper%dt_min
      flow_stepper%steps = 0
      flow_stepper%cumulative_newton_iterations = 0
      flow_stepper%cumulative_time_step_cuts = 0
      flow_stepper%cumulative_linear_iterations = 0
      flow_stepper%cumulative_solver_time = 0.d0
      flow_stepper%num_constant_time_steps = 0
      flow_stepper%num_newton_iterations = 0
      flow_stepper%prev_dt = 0.d0
    endif
    if (associated(tran_stepper)) then
      option%tran_dt = tran_stepper%dt_min
      tran_stepper%steps = 0
      tran_stepper%cumulative_newton_iterations = 0
      tran_stepper%cumulative_time_step_cuts = 0
      tran_stepper%cumulative_linear_iterations = 0
      tran_stepper%cumulative_solver_time = 0.d0
      tran_stepper%num_constant_time_steps = 0
      tran_stepper%num_newton_iterations = 0
      tran_stepper%prev_dt = 0.d0
    endif
    option%match_waypoint = PETSC_FALSE
    realization%output_option%plot_number = 0
  endif
  
  if (associated(flow_stepper)) flow_stepper%cur_waypoint => &
    WaypointSkipToTime(realization%waypoints,option%time)
  if (associated(tran_stepper)) tran_stepper%cur_waypoint => &
    WaypointSkipToTime(realization%waypoints,option%time)

  if (flow_read) then
    flow_stepper%target_time = option%flow_time
    select case(option%iflowmode)
      case(FLASH2_MODE)
        call Flash2UpdateAuxVars(realization)
      case(IMS_MODE)
        call ImmisUpdateAuxVars(realization)
      case(MPH_MODE)
        call MphaseUpdateAuxVars(realization)
      case(MIS_MODE)
        call MiscibleUpdateAuxVars(realization)
      case(TH_MODE)
        call THUpdateAuxVars(realization)
      case(THC_MODE)
        call THCUpdateAuxVars(realization)
      case(THMC_MODE)
        call THMCUpdateAuxVars(realization)
      case(RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case(G_MODE)                            ! do not update state
        call GeneralUpdateAuxVars(realization,PETSC_FALSE)
    end select    
  endif

  if (transport_read) then
    tran_stepper%target_time = option%tran_time
    ! This is here since we need to recalculate the secondary complexes
    ! if they exist.  DO NOT update activity coefficients!!! - geh
    if (realization%reaction%use_full_geochemistry) then
                                            ! cells     bcs        act coefs.
      call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)
    endif
  endif  
    
end subroutine SubsurfaceRestart

end module PFLOTRAN_Factory_module
