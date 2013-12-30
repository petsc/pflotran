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
!
! SubsurfaceInitialize: Sets up PFLOTRAN subsurface simulation 
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine SubsurfaceInitialize(simulation_base,option)

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
!
! SubsurfaceInitializePostPetsc: Sets up PFLOTRAN subsurface simulation 
!                                framework after to PETSc initialization
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine SubsurfaceInitializePostPetsc(simulation, option)

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
  ! no longer need simulation
  deallocate(simulation_old)
  call SubsurfaceJumpStart(simulation)
  
end subroutine SubsurfaceInitializePostPetsc

! ************************************************************************** !
!
! SubsurfInitCommandLineSettings: Initializes PFLTORAN subsurface output 
!                                 filenames, etc.
! author: Glenn Hammond
! date: 06/06/13
!
! ************************************************************************** !
subroutine SubsurfInitCommandLineSettings(option)

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
!
! HijackSimulation: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HijackSimulation(simulation_old,simulation)

  use Simulation_module
  use Realization_class
  use Option_module
  use Output_module, only : Output
  
  use PMC_Base_class
  use PMC_Subsurface_class  
  use Simulation_Base_class
  use Process_Model_General_class
  use Process_Model_Flash2_class
  use Process_Model_Immis_class
  use Process_Model_Mphase_class
  use Process_Model_Miscible_class
  use Process_Model_Richards_class
  use Process_Model_RT_class
  use Process_Model_TH_class
  use Process_Model_THC_class
  use Process_Model_Base_class
  use Process_Model_module
  use Timestepper_BE_class
  
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
      case(THC_MODE)
        cur_process_model => PMTHCCreate()
    end select
    cur_process_model%option => realization%option
    cur_process_model%output_option => realization%output_option

    flow_process_model_coupler => PMCSubsurfaceCreate()
    flow_process_model_coupler%option => option
    flow_process_model_coupler%pm_list => cur_process_model
    flow_process_model_coupler%pm_ptr%ptr => cur_process_model
!    flow_process_model_coupler%timestepper => simulation_old%flow_stepper
    flow_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%flow_stepper, &
                           flow_process_model_coupler%timestepper)
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
!    tran_process_model_coupler%timestepper => simulation_old%tran_stepper
    tran_process_model_coupler%realization => realization
    call HijackTimestepper(simulation_old%tran_stepper, &
                           tran_process_model_coupler%timestepper)
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
        ! set realization
        select type(cur_process_model)
          class is (pm_general_type)
            call cur_process_model%PMGeneralSetRealization(realization)
          class is (pm_flash2_type)
            call cur_process_model%PMFlash2SetRealization(realization)
          class is (pm_immis_type)
            call cur_process_model%PMImmisSetRealization(realization)
          class is (pm_mphase_type)
            call cur_process_model%PMMphaseSetRealization(realization)
          class is (pm_miscible_type)
            call cur_process_model%PMMiscibleSetRealization(realization)
          class is (pm_richards_type)
            call cur_process_model%PMRichardsSetRealization(realization)
          class is (pm_rt_type)
            call cur_process_model%PMRTSetRealization(realization)
          class is (pm_th_type)
            call cur_process_model%PMTHSetRealization(realization)
          class is (pm_thc_type)
            call cur_process_model%PMTHCSetRealization(realization)
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
              class is(stepper_BE_type)
              call SNESSetFunction( &
                             ts%solver%snes, &
                             cur_process_model%residual_vec, &
                             PMResidual, &
                             cur_process_model_coupler%pm_ptr,ierr)
              call SNESSetJacobian( &
                             ts%solver%snes, &
                             ts%solver%J, &
                             ts%solver%Jpre, &
                             PMJacobian, &
                             cur_process_model_coupler%pm_ptr,ierr)
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
!
! SubsurfaceJumpStart: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SubsurfaceJumpStart(simulation)

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
  class(stepper_base_type), pointer :: master_stepper
  class(stepper_BE_type), pointer :: flow_stepper
  class(stepper_BE_type), pointer :: tran_stepper
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
      class is(stepper_BE_type)
        flow_stepper => ts
    end select
  else
    nullify(flow_stepper)
  endif
  if (associated(simulation%rt_process_model_coupler)) then
    select type(ts => simulation%rt_process_model_coupler%timestepper)
      class is(stepper_BE_type)
        tran_stepper => ts
    end select
  else
    nullify(tran_stepper)
  endif
  nullify(master_stepper)
  
  option => realization%option
  output_option => realization%output_option

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr)
                             
  if (option%steady_state) then
    option%io_buffer = 'Running in steady-state not yet supported in refactored code.'
    call printErrMsg(option)
#if 0    
    call StepperRunSteadyState(realization,flow_stepper,tran_stepper)
#endif    
    ! do not want to run through time stepper
    option%status = DONE
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
 
!geh: now performed in PMRTInitializeRun()
!  if (associated(simulation%rt_process_model_coupler)) then
!    call simulation%rt_process_model_coupler%UpdateSolution()
!  endif

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
    option%status = DONE
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
    option%status = DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      option%status = FAIL
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
      option%status = FAIL
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
! HijackTimestepper: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HijackTimestepper(stepper_old,stepper_base)

  use Timestepper_BE_class
  use Timestepper_Base_class
  use Timestepper_module

  implicit none
  
  type(stepper_type), pointer :: stepper_old
  class(stepper_base_type), pointer :: stepper_base
  
  class(stepper_BE_type), pointer :: stepper

  stepper => TimestepperBECreate()
  
  stepper%steps = stepper_old%steps
  stepper%num_newton_iterations = stepper_old%num_newton_iterations
  stepper%num_linear_iterations = stepper_old%num_linear_iterations
  stepper%num_constant_time_steps = stepper_old%num_constant_time_steps

  stepper%max_time_step = stepper_old%max_time_step
  stepper%max_time_step_cuts = stepper_old%max_time_step_cuts
  stepper%constant_time_step_threshold = stepper_old%constant_time_step_threshold
  stepper%iaccel = stepper_old%iaccel

  stepper%cumulative_newton_iterations = stepper_old%cumulative_newton_iterations
  stepper%cumulative_linear_iterations = stepper_old%cumulative_linear_iterations
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

  stepper%ntfac = stepper_old%ntfac
  stepper%tfac => stepper_old%tfac
  nullify(stepper_old%tfac)
  
  stepper%init_to_steady_state = stepper_old%init_to_steady_state
  stepper%steady_state_rel_tol = stepper_old%steady_state_rel_tol
  stepper%run_as_steady_state = stepper_old%run_as_steady_state

  stepper%solver => stepper_old%solver
  nullify(stepper_old%solver)
  stepper%convergence_context => stepper_old%convergence_context
  nullify(stepper_old%convergence_context)
  stepper%cur_waypoint => stepper_old%cur_waypoint
  nullify(stepper_old%cur_waypoint)
!  stepper%prev_waypoint => stepper_old%prev_waypoint
!  nullify(stepper_old%prev_waypoint)
  
!  stepper%revert_dt = stepper_old%revert_dt
!  stepper%num_contig_revert_due_to_sync = &
!  stepper_old%num_contig_revert_due_to_sync

  stepper_base => stepper
  
end subroutine HijackTimestepper

end module Subsurface_Factory_module
