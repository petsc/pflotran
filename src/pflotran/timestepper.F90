module Timestepper_module
 
  use Solver_module
  use Waypoint_module 
  use Convergence_module 
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public :: stepper_type
  
    PetscInt :: steps         ! The number of time steps taken by the code.
    PetscInt :: num_constant_time_steps   ! number of contiguous time_steps of constant size
    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step

    PetscInt :: max_time_step                ! Maximum number of time steps to be taken by the code.
    PetscInt :: max_time_step_cuts           ! Maximum number of timestep cuts within one time step.
    PetscInt :: constant_time_step_threshold        ! Steps needed after cutting to increase time step
    PetscInt :: iaccel        ! Accelerator index

    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations
    PetscInt :: cumulative_time_step_cuts       ! Total number of cuts in the timestep taken.    
    PetscReal :: cumulative_solver_time
    
    PetscReal :: dt_min
    PetscReal :: dt_max
    PetscReal :: prev_dt
    PetscReal :: cfl_limiter
    PetscReal :: cfl_limiter_ts
    
    PetscBool :: init_to_steady_state
    PetscBool :: run_as_steady_state
    PetscReal :: steady_state_rel_tol
    
    PetscBool :: time_step_cut_flag  ! flag toggled if timestep is cut
    
    PetscInt :: start_time_step ! the first time step of a given run
    PetscReal :: time_step_tolerance ! scalar used in determining time step size
    PetscReal :: target_time    ! time at end of "synchronized" time step 

    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac
            
    type(solver_type), pointer :: solver
    
    type(waypoint_type), pointer :: cur_waypoint

    type(convergence_context_type), pointer :: convergence_context

#ifdef SURFACE_FLOW
    PetscInt :: steps_surf_flow         ! The number of time steps taken by the code.
    !PetscInt :: num_constant_time_steps   ! number of contiguous time_steps of constant size
    PetscInt :: num_newton_iterations_surf_flow ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations_surf_flow ! number of linear solver iterations in a time step

    PetscInt :: cumulative_newton_iterations_surf_flow       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations_surf_flow     ! Total number of linear iterations
    PetscInt :: cumulative_time_step_cuts_surf_flow       ! Total number of cuts in the timestep taken.    
    PetscReal :: cumulative_solver_time_surf_flow
#endif
  end type stepper_type
  
  public :: TimestepperCreate, TimestepperDestroy, StepperRun, &
            TimestepperRead, TimestepperPrintInfo, TimestepperReset

contains

! ************************************************************************** !
!
! TimestepperCreate: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function TimestepperCreate()

  implicit none
  
  type(stepper_type), pointer :: TimestepperCreate
  
  type(stepper_type), pointer :: stepper
  
  allocate(stepper)
  stepper%steps = 0
  stepper%num_newton_iterations = 0
  stepper%num_linear_iterations = 0
  stepper%num_constant_time_steps = 0

  stepper%max_time_step = 999999
  stepper%max_time_step_cuts = 16
  stepper%constant_time_step_threshold = 5
  stepper%iaccel = 5

  stepper%cumulative_newton_iterations = 0
  stepper%cumulative_linear_iterations = 0
  stepper%cumulative_time_step_cuts = 0    
  stepper%cumulative_solver_time = 0.d0

  stepper%start_time_step = 0
  stepper%time_step_tolerance = 0.1d0
  stepper%target_time = 0.d0
  
  stepper%dt_min = 1.d0
  stepper%dt_max = 3.1536d6 ! One-tenth of a year.  
  stepper%prev_dt = 0.d0
  stepper%cfl_limiter = -999.d0
  stepper%cfl_limiter_ts = 1.d20
  
  stepper%time_step_cut_flag = PETSC_FALSE

  stepper%ntfac = 13
  allocate(stepper%tfac(13))
  stepper%tfac(1)  = 2.0d0; stepper%tfac(2)  = 2.0d0
  stepper%tfac(3)  = 2.0d0; stepper%tfac(4)  = 2.0d0
  stepper%tfac(5)  = 2.0d0; stepper%tfac(6)  = 1.8d0
  stepper%tfac(7)  = 1.6d0; stepper%tfac(8)  = 1.4d0
  stepper%tfac(9)  = 1.2d0; stepper%tfac(10) = 1.0d0
  stepper%tfac(11) = 1.0d0; stepper%tfac(12) = 1.0d0
  stepper%tfac(13) = 1.0d0
  
  stepper%init_to_steady_state = PETSC_FALSE
  stepper%steady_state_rel_tol = 1.d-8
  stepper%run_as_steady_state = PETSC_FALSE
  
  nullify(stepper%solver)
  nullify(stepper%convergence_context)
  nullify(stepper%cur_waypoint)
  
  stepper%solver => SolverCreate()
  
  TimeStepperCreate => stepper
  
end function TimestepperCreate 

! ************************************************************************** !
!
! TimestepperReset: Resets time stepper back to initial settings
! author: Glenn Hammond
! date: 01/27/11
!
! ************************************************************************** !
subroutine TimestepperReset(stepper,dt_min)

  implicit none

  type(stepper_type) :: stepper
  PetscReal :: dt_min

  stepper%steps = 0
  stepper%num_newton_iterations = 0
  stepper%num_linear_iterations = 0
  stepper%num_constant_time_steps = 0

  stepper%cumulative_newton_iterations = 0
  stepper%cumulative_linear_iterations = 0
  stepper%cumulative_time_step_cuts = 0
  stepper%cumulative_solver_time = 0.d0

  stepper%start_time_step = 0
  stepper%target_time = 0.d0

  stepper%dt_min = dt_min
  stepper%prev_dt = 0.d0
  stepper%cfl_limiter = -999.d0
  stepper%cfl_limiter_ts = 1.d20

  stepper%time_step_cut_flag = PETSC_FALSE

end subroutine TimestepperReset

! ************************************************************************** !
!
! TimestepperRead: Reads parameters associated with time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperRead(stepper,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none

  type(stepper_type) :: stepper
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))

      case('NUM_STEPS_AFTER_TS_CUT')
        call InputReadInt(input,option,stepper%constant_time_step_threshold)
        call InputDefaultMsg(input,option,'num_constant_time_steps_after_ts_cut')

      case('MAX_STEPS')
        call InputReadInt(input,option,stepper%max_time_step)
        call InputDefaultMsg(input,option,'max_time_step')
  
      case('TS_ACCELERATION')
        call InputReadInt(input,option,stepper%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

      case('MAX_TS_CUTS')
        call InputReadInt(input,option,stepper%max_time_step_cuts)
        call InputDefaultMsg(input,option,'max_time_step_cuts')
        
      case('CFL_LIMITER')
        call InputReadDouble(input,option,stepper%cfl_limiter)
        call InputDefaultMsg(input,option,'cfl limiter')
        
      case('DT_FACTOR')
        string='time_step_factor'
        call UtilityReadArray(stepper%tfac,NEG_ONE_INTEGER,string,input, &
            option)
        stepper%ntfac = size(stepper%tfac)

      case('INITIALIZE_TO_STEADY_STATE')
        stepper%init_to_steady_state = PETSC_TRUE
        call InputReadDouble(input,option,stepper%steady_state_rel_tol)
        call InputDefaultMsg(input,option,'steady state convergence relative tolerance')

      case('RUN_AS_STEADY_STATE')
        stepper%run_as_steady_state = PETSC_TRUE

      case('MAX_PRESSURE_CHANGE')
        call InputReadDouble(input,option,option%dpmxe)
        call InputDefaultMsg(input,option,'dpmxe')

      case('MAX_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,option%dtmpmxe)
        call InputDefaultMsg(input,option,'dtmpmxe')
  
      case('MAX_CONCENTRATION_CHANGE')
        call InputReadDouble(input,option,option%dcmxe)
        call InputDefaultMsg(input,option,'dcmxe')

      case('MAX_SATURATION_CHANGE')
        call InputReadDouble(input,option,option%dsmxe)
        call InputDefaultMsg(input,option,'dsmxe')

      case('PRESSURE_DAMPENING_FACTOR')
        call InputReadDouble(input,option,option%pressure_dampening_factor)
        call InputErrorMsg(input,option,'PRESSURE_DAMPENING_FACTOR', &
                           'TIMESTEPPER')

      case('SATURATION_CHANGE_LIMIT')
        call InputReadDouble(input,option,option%saturation_change_limit)
        call InputErrorMsg(input,option,'SATURATION_CHANGE_LIMIT', &
                           'TIMESTEPPER')
                           
      case('PRESSURE_CHANGE_LIMIT')
        call InputReadDouble(input,option,option%pressure_change_limit)
        call InputErrorMsg(input,option,'PRESSURE_CHANGE_LIMIT', &
                           'TIMESTEPPER')
                           
      case('TEMPERATURE_CHANGE_LIMIT')
        call InputReadDouble(input,option,option%temperature_change_limit)
        call InputErrorMsg(input,option,'TEMPERATURE_CHANGE_LIMIT', &
                           'TIMESTEPPER')

      case default
        option%io_buffer = 'Timestepper option: '//trim(keyword)// &
                           ' not recognized.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine TimestepperRead

! ************************************************************************** !
!
! StepperRun: Runs the time step loop
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine StepperRun(realization,surf_realization,flow_stepper,tran_stepper,surf_flow_stepper)
#else
subroutine StepperRun(realization,flow_stepper,tran_stepper)
#endif

  use Realization_module

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputVectorTecplot, &
                            OutputPrintCouplers
  use Logging_module  
  use Discretization_module
  use Condition_Control_module
#ifdef SURFACE_FLOW
  use Surface_Flow_module
  use Surface_Realization_module
#endif
  implicit none
  
#include "finclude/petscdef.h"
#include "finclude/petsclog.h"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
#ifdef SURFACE_FLOW
  type(stepper_type), pointer :: surf_flow_stepper
  type(surface_realization_type), pointer :: surf_realization
  PetscBool :: plot_flag_surf, transient_plot_flag_surf
!  PetscReal :: surf_flow_time,surf_flow_dt,surf_flow_target_time
#endif
  
  type(stepper_type), pointer :: master_stepper
  type(stepper_type), pointer :: null_stepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(waypoint_type), pointer :: prev_waypoint  

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: plot_flag, stop_flag, transient_plot_flag
  PetscBool :: activity_coefs_read
  PetscBool :: flow_read
  PetscBool :: transport_read
  PetscBool :: step_to_steady_state
  PetscBool :: run_flow_as_steady_state
  PetscBool :: failure, surf_failure
  PetscLogDouble :: start_time, end_time
  PetscReal :: tran_dt_save, flow_t0
  PetscReal :: dt_cfl_1, flow_to_tran_ts_ratio

  PetscLogDouble :: stepper_start_time, current_time, average_step_time
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option

  nullify(master_stepper,null_stepper)

  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-vecload_block_size", & 
                           failure, ierr)
                             
  if (option%steady_state) then
    call StepperRunSteadyState(realization,flow_stepper,tran_stepper)
    return 
  endif
  
  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
  else
    master_stepper => tran_stepper
  endif

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  stop_flag = PETSC_FALSE
  activity_coefs_read = PETSC_FALSE
  flow_read = PETSC_FALSE
  transport_read = PETSC_FALSE
  step_to_steady_state = PETSC_FALSE
  run_flow_as_steady_state = PETSC_FALSE
  failure = PETSC_FALSE
  
  if (option%restart_flag) then
    call StepperRestart(realization,flow_stepper,tran_stepper, &
                        flow_read,transport_read,activity_coefs_read)
    if (associated(flow_stepper)) flow_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)
    if (associated(tran_stepper)) tran_stepper%cur_waypoint => &
      WaypointSkipToTime(realization%waypoints,option%time)

    if (flow_read) then
      flow_stepper%target_time = option%flow_time
      call StepperUpdateFlowAuxVars(realization)
    endif

    if (transport_read) then
      tran_stepper%target_time = option%tran_time
      ! This is here since we need to recalculate the secondary complexes
      ! if they exist.  DO NOT update activity coefficients!!! - geh
      if (realization%reaction%use_full_geochemistry) then
        call StepperUpdateTranAuxVars(realization)
        !call StepperSandbox(realization)
      endif
    endif

  else if (master_stepper%init_to_steady_state) then
    option%print_screen_flag = OptionPrintToScreen(option)
    option%print_file_flag = OptionPrintToFile(option)
    if (associated(flow_stepper)) then
      if (flow_stepper%init_to_steady_state) then
        step_to_steady_state = PETSC_TRUE
        option%flow_dt = master_stepper%dt_min
        call StepperStepFlowDT(realization,flow_stepper,step_to_steady_state, &
                               failure)
        if (failure) then ! if flow solve fails, exit
          if (OptionPrintToScreen(option)) then
            write(*,*) ' ERROR: steady state solve failed!!!'
          endif
          if (OptionPrintToFile(option)) then
            write(option%fid_out,*) ' ERROR: steady state solve failed!!!'
          endif
          return 
        endif
        option%flow_dt = master_stepper%dt_min
        run_flow_as_steady_state = flow_stepper%run_as_steady_state
      endif
      if (run_flow_as_steady_state) then
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
  endif

  if (flow_read .and. option%overwrite_restart_flow) then
    call RealizationRevertFlowParameters(realization)
  endif

  if (transport_read .and. option%overwrite_restart_transport) then
    call CondControlAssignTranInitCond(realization)  
  endif

  ! turn on flag to tell RTUpdateSolution that the code is not timestepping
#ifdef SURFACE_FLOW
  call StepperUpdateSolution(realization,surf_realization)
#else
  call StepperUpdateSolution(realization)
#endif

  if (option%jumpstart_kinetic_sorption .and. option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. option%restart_flag) then
      option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(option)
    endif
    call StepperJumpStart(realization)
  endif
  
  call PetscLogStagePop(ierr)
  option%init_stage = PETSC_FALSE
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
  call OutputInit(realization,master_stepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_stepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(realization,plot_flag,transient_plot_flag)
#ifdef SURFACE_FLOW
    plot_flag_surf = PETSC_TRUE
    transient_plot_flag_surf = PETSC_TRUE
    call Output(surf_realization,realization,plot_flag_surf,transient_plot_flag_surf)
#endif
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 0, print out initial condition only
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

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
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
      return
    else
      tran_stepper%dt_max = tran_stepper%cur_waypoint%dt_max
    endif
  endif
           
  ! ensure that steady_state flag is off
  step_to_steady_state = PETSC_FALSE
  call PetscGetTime(stepper_start_time, ierr)
  if (associated(flow_stepper)) &
    flow_stepper%start_time_step = flow_stepper%steps + 1
  if (associated(tran_stepper)) &
    tran_stepper%start_time_step = tran_stepper%steps + 1
  
  if (realization%debug%print_couplers) then
    call OutputPrintCouplers(realization,ZERO_INTEGER)
  endif

  do

    if (OptionPrintToScreen(option) .and. &
        mod(master_stepper%steps,output_option%screen_imod) == 0) then
      option%print_screen_flag = PETSC_TRUE
    else
      option%print_screen_flag = PETSC_FALSE
    endif

    if (OptionPrintToFile(option) .and. &
        mod(master_stepper%steps,output_option%output_file_imod) == 0) then
      option%print_file_flag = PETSC_TRUE
    else
      option%print_file_flag = PETSC_FALSE
    endif

    prev_waypoint => master_stepper%cur_waypoint
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE
    
#ifdef SURFACE_FLOW

    if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then

      ! Solve surface-flow model
      if (associated(surf_flow_stepper) .and. .not.run_flow_as_steady_state) then

        ! Update model coupling time
        call SetSurfaceSubsurfaceCouplingTime(flow_stepper,tran_stepper,surf_flow_stepper, &
                            option,plot_flag,transient_plot_flag)

        ! Update subsurface pressure of top soil layer for surface flow model
        if (surf_realization%option%subsurf_surf_coupling == SEQ_COUPLED) then
          call SurfaceRealizationUpdateSurfaceBC(realization,surf_realization)
        endif

        surf_failure = PETSC_FALSE

        do ! loop on surface-flow until it reaches the target time

          ! Compute flux between surface-subsurface model
          if (surf_realization%option%subsurf_surf_coupling == SEQ_COUPLED) then
            call SurfaceRealizationComputeSurfaceSubsurfFlux(realization,surf_realization)
          endif
          
          ! Solve surface flow
          call StepperStepSurfaceFlowDT(surf_realization,surf_flow_stepper, &
                                        surf_failure)

          if(surf_failure) return
          ! update time for surface flow model
          option%surf_flow_time=option%surf_flow_time+option%surf_flow_dt

          ! if target time reached, we are done
          if((option%surf_subsurf_coupling_time-surf_flow_stepper%target_time) &
              /option%surf_flow_dt<1.d-10) exit

          ! if still stepping, update solution
          call StepperUpdateSurfaceFlowSolution(surf_realization)
          call StepperUpdateSurfaceFlowDT(surf_flow_stepper,option)

          ! Set new target time for surface model
          call StepperSetSurfaceFlowTargetTimes(surf_flow_stepper, &
                                          option,plot_flag,transient_plot_flag)
        enddo
      endif

      ! Solve subsurface model
      if (associated(surf_flow_stepper).and. &
          (surf_realization%option%subsurf_surf_coupling==SEQ_COUPLED)) then

        flow_t0 = option%flow_time
        call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)

        ! Update source/sink condition for subsurface flow model
        call SurfaceRealizationUpdateSubsurfaceBC(realization,surf_realization, &
              option%surf_subsurf_coupling_time-option%flow_time)

        do
          ! Solve subsurface flow
          call StepperStepFlowDT(realization,flow_stepper,step_to_steady_state, &
                                 failure)
          
          if (failure) return ! if flow solve fails, exit
          option%flow_time = flow_stepper%target_time

          ! If target time reached, we are done
          if((option%surf_subsurf_coupling_time-flow_stepper%target_time) &
              /option%flow_dt<1.d-10)exit
          
          ! If still stepping, update the solution and update dt
          call StepperUpdateFlowSolution(realization)
          call StepperUpdateDT(flow_stepper,tran_stepper,option)

          ! Set new target time for subsurface model
          call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                                     surf_flow_stepper, &
                                     option,plot_flag, &
                                     transient_plot_flag)
        enddo
        call PetscLogStagePop(ierr)
      
      else
        if(.not.associated(surf_flow_stepper)) then
          flow_t0 = option%flow_time
          call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                                     surf_flow_stepper, &
                                     option,plot_flag, &
                                     transient_plot_flag)
          call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
          call StepperStepFlowDT(realization,flow_stepper,step_to_steady_state, &
                                failure)
          call PetscLogStagePop(ierr)
          if (failure) return ! if flow solve fails, exit
          option%flow_time = flow_stepper%target_time
        endif
      endif ! if SEQ_COUPLED
   endif ! associated(flow_stepper)

#else

    call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                               option,plot_flag, &
                               transient_plot_flag)

    ! flow solution
    if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then

      flow_t0 = option%flow_time
      call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
      call StepperStepFlowDT(realization,flow_stepper,step_to_steady_state, &
                              failure)
      call PetscLogStagePop(ierr)
      if (failure) return ! if flow solve fails, exit
      option%flow_time = flow_stepper%target_time
    endif
#endif
    
    ! (reactive) transport solution
    if (associated(tran_stepper)) then
      call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
      tran_dt_save = -999.d0
      ! reset transpor time step if flow time step is cut
      if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
        if (flow_stepper%time_step_cut_flag) then
          tran_stepper%target_time = flow_stepper%target_time
          option%tran_dt = min(option%tran_dt,option%flow_dt)
        endif
      endif
      ! CFL Limiting
      if (tran_stepper%cfl_limiter > 0.d0) then
        call TimestepperCheckCFLLimit(tran_stepper,realization)
        call TimestepperEnforceCFLLimit(tran_stepper,option,output_option)
      endif
      ! if transport time step is less than flow step, but greater than
      ! half the flow time step, might as well cut it down to half the 
      ! flow time step size
      if (associated(flow_stepper) .and. .not.run_flow_as_steady_state) then
        flow_to_tran_ts_ratio = option%flow_dt / option%tran_dt
        if (flow_to_tran_ts_ratio > 1.d0 .and. &
            flow_to_tran_ts_ratio < 2.d0) then
          option%tran_dt = option%flow_dt * 0.5d0
        endif
      endif
      do ! loop on transport until it reaches the target time
        if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
          !global implicit
          call StepperStepTransportDT_GI(realization,tran_stepper, &
                                         run_flow_as_steady_state, &
                                         flow_t0,option%flow_time,failure)
        else
          !operator splitting
          call StepperStepTransportDT_OS(realization,tran_stepper, &
                                         run_flow_as_steady_state, &
                                         flow_t0,option%flow_time,failure)
        endif
        if (failure) then ! if transport solve fails, exit
          call PetscLogStagePop(ierr)
          return 
        endif
        
        ! update transport time
        option%tran_time = option%tran_time + option%tran_dt

        ! if target time reached, we are done
        if ((tran_stepper%target_time-option%tran_time) / &
            option%tran_time < 1.d-10) exit

        call StepperUpdateTransportSolution(realization)
        
        ! if still stepping, update the time step size based on convergence
        ! criteria
        call StepperUpdateDT(null_stepper,tran_stepper,option)
        ! check CFL limit if turned on
        if (tran_stepper%cfl_limiter > 0.d0) then
          call TimestepperEnforceCFLLimit(tran_stepper,option,output_option)
        endif         
        ! if current step exceeds (factoring in tolerance) exceeds the target
        ! time, set the step size to the difference
        if ((option%tran_time + &
             (1.d0+tran_stepper%time_step_tolerance)*option%tran_dt) > &
            tran_stepper%target_time) then
          ! tran_dt_save enables us to regain the original tran_dt if
          ! it is cut to match the target time
          tran_dt_save = option%tran_dt
          option%tran_dt = tran_stepper%target_time - option%tran_time
        endif
      enddo
      ! if substepping occured and time step size was modified to 
      ! to match the stepper target time, need to set the transport
      ! step size back to the pre-modified value
      if (tran_dt_save > -998.d0) option%tran_dt = tran_dt_save
      option%tran_time = tran_stepper%target_time
      call PetscLogStagePop(ierr)
    endif

    if (realization%debug%print_couplers) then
      call OutputPrintCouplers(realization,master_stepper%steps)
    endif  

    ! update solution variables
    
    option%time = master_stepper%target_time
#ifdef SURFACE_FLOW
    call StepperUpdateSolution(realization,surf_realization)
#else
    call StepperUpdateSolution(realization)
#endif

    ! if a time step cut has occured, need to set the below back to original values
    ! if they changed. 
    if (master_stepper%time_step_cut_flag) then
      master_stepper%cur_waypoint => prev_waypoint
      plot_flag = PETSC_FALSE
    endif
    ! however, if we are using the modulus of the output_option%imod, we may still print
    if (mod(master_stepper%steps, &
            output_option%periodic_output_ts_imod) == 0) then
      plot_flag = PETSC_TRUE
    endif
    if (plot_flag .or. mod(master_stepper%steps, &
                           output_option%periodic_tr_output_ts_imod) == 0) then
      transient_plot_flag = PETSC_TRUE
    endif

!   deprecated - geh
!    if (option%compute_mass_balance) then
!      call MassBalanceUpdate(realization,flow_stepper%solver, &
!                             tran_stepper%solver)
!    endif
#ifdef SURFACE_FLOW
    plot_flag_surf = plot_flag
    transient_plot_flag_surf = transient_plot_flag
#endif
    call Output(realization,plot_flag,transient_plot_flag)
    
    call StepperUpdateDTMax(flow_stepper,tran_stepper,option)
    call StepperUpdateDT(flow_stepper,tran_stepper,option)

#ifdef SURFACE_FLOW
    call Output(surf_realization,realization,plot_flag_surf,transient_plot_flag_surf)
    if(associated(surf_flow_stepper)) then
      call StepperUpdateSurfaceFlowDT(surf_flow_stepper,option)
    endif
#endif

    ! if a simulation wallclock duration time is set, check to see that the
    ! next time step will not exceed that value.  If it does, print the
    ! checkpoint and exit.
    if (option%wallclock_stop_flag) then
      call PetscGetTime(current_time, ierr)
      average_step_time = (current_time-stepper_start_time)/ &
                          real(master_stepper%steps-&
                               master_stepper%start_time_step+1) &
                          *2.d0  ! just to be safe, double it
      if (average_step_time + current_time > option%wallclock_stop_time) then
        call printMsg(option,"Wallclock stop time exceeded.  Exiting!!!")
        call printMsg(option,"")
        stop_flag = PETSC_TRUE
      endif
    endif

    if (option%checkpoint_flag .and. &
        mod(master_stepper%steps,option%checkpoint_frequency) == 0) then
      call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                             master_stepper%steps)  
    endif
    
    ! if at end of waypoint list (i.e. cur_waypoint = null), we are done!
    if (.not.associated(master_stepper%cur_waypoint) .or. stop_flag) exit

  enddo

  if (master_stepper%steps >= master_stepper%max_time_step) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
  endif

  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                           NEG_ONE_INTEGER)  
  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%cumulative_newton_iterations, &
            flow_stepper%cumulative_linear_iterations,flow_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(*,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%cumulative_newton_iterations, &
            tran_stepper%cumulative_linear_iterations,tran_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(*,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif
  
  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            flow_stepper%steps,flow_stepper%cumulative_newton_iterations, &
            flow_stepper%cumulative_linear_iterations,flow_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') flow_stepper%cumulative_solver_time
      write(option%fid_out,*) 'FLOW SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            tran_stepper%steps,tran_stepper%cumulative_newton_iterations, &
            tran_stepper%cumulative_linear_iterations,tran_stepper%cumulative_time_step_cuts
      write(string,'(f12.1)') tran_stepper%cumulative_solver_time
      write(option%fid_out,*) 'TRAN SNES time = ' // trim(adjustl(string)) // ' seconds'
    endif            
  endif

  call PetscLogStagePop(ierr)

end subroutine StepperRun

! ************************************************************************** !
!
! StepperUpdateDT: Updates time step
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperUpdateDT(flow_stepper,tran_stepper,option)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(option_type) :: option
  
  PetscReal :: time, dt
  PetscReal :: fac,dtt,up,utmp,uc,ut,uus,dt_tfac,dt_p
  PetscBool :: update_time_step
  PetscInt :: ifac

  ! FLOW
  update_time_step = PETSC_TRUE
  if (associated(flow_stepper)) then
    if (flow_stepper%run_as_steady_state) then
      update_time_step = PETSC_FALSE
    else
      ! if the time step was cut, set number of constant time steps to 1
      if (flow_stepper%time_step_cut_flag) then
        flow_stepper%time_step_cut_flag = PETSC_FALSE
        flow_stepper%num_constant_time_steps = 1
      ! otherwise, only increment if teh constant time step counter was
      ! initialized to 1
      else if (flow_stepper%num_constant_time_steps > 0) then
        flow_stepper%num_constant_time_steps = &
          flow_stepper%num_constant_time_steps + 1
      endif
    
      ! num_constant_time_steps = 0: normal time stepping with growing steps
      ! num_constant_time_steps > 0: restriction of constant time steps until
      !                              constant_time_step_threshold is met
      if (flow_stepper%num_constant_time_steps > &
          flow_stepper%constant_time_step_threshold) then
        flow_stepper%num_constant_time_steps = 0
      else if (flow_stepper%num_constant_time_steps > 0) then
        ! do not increase time step size
        update_time_step = PETSC_FALSE
      endif
    endif
    
    if (update_time_step) then
      
      time = option%flow_time
      dt = option%flow_dt

      if (flow_stepper%iaccel == 0) return

      select case(option%iflowmode)
        case(FLASH2_MODE)   
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uus)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(IMS_MODE)   
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uus)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(MIS_MODE)   
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            uc = option%dcmxe/(option%dcmax+1.d-6)
            ut = min(up,uc)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(MPH_MODE)   
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uc = option%dcmxe/(option%dcmax+1.d-6)
            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uc,uus)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(THC_MODE)
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uus)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(THMC_MODE)
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uus)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case(RICHARDS_MODE)
          if (flow_stepper%iaccel > 0) then
            fac = 0.5d0
            if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
              fac = 0.33d0
              ut = 0.d0
            else
              up = option%dpmxe/(option%dpmax+0.1)
              ut = up
            endif
            dtt = fac * dt * (1.d0 + ut)
          else
            ifac = max(min(flow_stepper%num_newton_iterations, &
                           flow_stepper%ntfac),1)
            dt_tfac = flow_stepper%tfac(ifac) * dt

            fac = 0.5d0
            up = option%dpmxe/(option%dpmax+0.1)
            dt_p = fac * dt * (1.d0 + up)

            dtt = min(dt_tfac,dt_p)
          endif
        case(G_MODE)   
          fac = 0.5d0
          if (flow_stepper%num_newton_iterations >= flow_stepper%iaccel) then
            fac = 0.33d0
            ut = 0.d0
          else
            up = option%dpmxe/(option%dpmax+0.1)
            utmp = option%dtmpmxe/(option%dtmpmax+1.d-5)
            uc = option%dcmxe/(option%dcmax+1.d-6)
!            uus= option%dsmxe/(option%dsmax+1.d-6)
            ut = min(up,utmp,uc)
          endif
          dtt = fac * dt * (1.d0 + ut)
        case default
          dtt = dt
          if (flow_stepper%num_newton_iterations <= flow_stepper%iaccel .and. &
              flow_stepper%num_newton_iterations <= size(flow_stepper%tfac)) then
            if (flow_stepper%num_newton_iterations == 0) then
              dtt = flow_stepper%tfac(1) * dt
            else
              dtt = flow_stepper%tfac(flow_stepper%num_newton_iterations) * dt
            endif
          endif
          
      end select

      if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
      if (dtt > flow_stepper%dt_max) dtt = flow_stepper%dt_max
      ! geh: There used to be code here that cut the time step if it is too
      !      large relative to the simulation time.  This has been removed.
      dt = dtt

      option%flow_dt = dt
    
    endif
        
  endif
  
  ! TRANSPORT
  update_time_step = PETSC_TRUE
  if (associated(tran_stepper)) then
    ! if the time step was cut, set number of constant time steps to 1
    if (tran_stepper%time_step_cut_flag) then
      tran_stepper%time_step_cut_flag = PETSC_FALSE
      tran_stepper%num_constant_time_steps = 1
    ! otherwise, only increment if teh constant time step counter was
    ! initialized to 1
    else if (tran_stepper%num_constant_time_steps > 0) then
      tran_stepper%num_constant_time_steps = &
        tran_stepper%num_constant_time_steps + 1
    endif
    
    ! num_constant_time_steps = 0: normal time stepping with growing steps
    ! num_constant_time_steps > 0: restriction of constant time steps until
    !                              constant_time_step_threshold is met
    if (tran_stepper%num_constant_time_steps > &
        tran_stepper%constant_time_step_threshold) then
      tran_stepper%num_constant_time_steps = 0
    else if (tran_stepper%num_constant_time_steps > 0) then
      ! do not increase time step size
      update_time_step = PETSC_FALSE
    endif
    
    if (update_time_step) then
    
      time = option%tran_time
      dt = option%tran_dt

      if (tran_stepper%iaccel == 0) return

      dtt = dt
      if (tran_stepper%num_newton_iterations <= tran_stepper%iaccel) then
        if (tran_stepper%num_newton_iterations <= size(tran_stepper%tfac)) then
          dtt = tran_stepper%tfac(tran_stepper%num_newton_iterations) * dt
        else
          dtt = 0.5d0 * dt
        endif
      else
!       dtt = 2.d0 * dt
        dtt = 0.5d0 * dt
      endif

      if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
      if (dtt > tran_stepper%dt_max) dtt = tran_stepper%dt_max
      ! geh: see comment above under flow stepper
      dt = dtt

      option%tran_dt = dt
        
    endif
    
  endif

  ! ensure that transport time step is not larger than flow time step
  if (associated(flow_stepper) .and. associated(tran_stepper)) then
    if (.not.flow_stepper%run_as_steady_state) then
      option%tran_dt = min(option%tran_dt,option%flow_dt)
    endif
  endif

end subroutine StepperUpdateDT

! ************************************************************************** !
!> This subroutine sets updates dt for surface flow.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/31/12
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine StepperUpdateSurfaceFlowDT(surf_flow_stepper,option)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: surf_flow_stepper
  type(option_type) :: option

  PetscReal :: time, dt
  PetscReal :: fac,dtt,up,utmp,uc,ut,uus,dt_tfac,dt_p
  PetscBool :: update_time_step
  PetscInt :: ifac

  ! FLOW
  update_time_step = PETSC_TRUE

  ! if the time step was cut, set number of constant time steps to 1
  if(surf_flow_stepper%time_step_cut_flag) then
    surf_flow_stepper%time_step_cut_flag=PETSC_FALSE
    surf_flow_stepper%num_constant_time_steps=1
  ! otherwise, only increment if the constant time step counter was
  ! initialized to 1
  else if (surf_flow_stepper%num_constant_time_steps > 0) then
    surf_flow_stepper%num_constant_time_steps = &
      surf_flow_stepper%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (surf_flow_stepper%num_constant_time_steps > &
      surf_flow_stepper%constant_time_step_threshold) then
    surf_flow_stepper%num_constant_time_steps = 0
  else if (surf_flow_stepper%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif

  if(update_time_step) then

    time=option%surf_flow_time
    dt=option%surf_flow_dt

    if(surf_flow_stepper%iaccel==0) return
    
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        fac=0.5d0
        if(surf_flow_stepper%num_newton_iterations>=surf_flow_stepper%iaccel) then
          fac=0.33d0
          ut=0.d0
        else
          !up = option%dpmxe/(option%dpmax+0.1)
          up = 2.d0
          ut = up
        endif
        dtt = fac * dt * (1.d0 + ut)
    end select

    if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
    if (dtt > surf_flow_stepper%dt_max) dtt = surf_flow_stepper%dt_max

    dt = dtt
    option%surf_flow_dt = dt

  endif

end subroutine StepperUpdateSurfaceFlowDT
#endif
! ************************************************************************** !
!
! StepperUpdateDTMax: Updates maximum time step specified by the current
!                     waypoint after the completion of a time step
! author: Glenn Hammond
! date: 05/23/12
!
! ************************************************************************** !
subroutine StepperUpdateDTMax(flow_stepper,tran_stepper,option)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(option_type) :: option
  
  PetscReal :: dt_max
  PetscBool :: flag
  type(waypoint_type), pointer :: cur_waypoint

  ! this just works through the logic of whether the flow stepper should
  ! be used to set dt_max.  If flow is not present, or is to
  ! be run as steady state, use the transport parameters
  flag = PETSC_TRUE
  if (associated(flow_stepper)) then
    if (flow_stepper%run_as_steady_state) then
      flag = PETSC_FALSE
    endif
  else
    flag = PETSC_FALSE
  endif

  if (flag) then ! flow stepper will govern the target time
    cur_waypoint => flow_stepper%cur_waypoint
  else
    cur_waypoint => tran_stepper%cur_waypoint
  endif
  
  ! update maximum time step size to current waypoint value
  if (associated(cur_waypoint)) then
    dt_max = cur_waypoint%dt_max
    
    ! target time will always be dictated by the flow solver, if present
    if (associated(flow_stepper)) then
      flow_stepper%dt_max = dt_max
    endif
    if (associated(tran_stepper)) then
      tran_stepper%dt_max = dt_max
    endif

  endif
  
end subroutine StepperUpdateDTMax

! ************************************************************************** !
!
! StepperSetTargetTimes: Sets target time for flow and transport solvers
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperSetTargetTimes(flow_stepper,tran_stepper, &
#ifdef SURFACE_FLOW
                                 surf_flow_stepper, &
#endif
                                 option,plot_flag, &
                                 transient_plot_flag)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: flow_stepper, tran_stepper
  type(option_type) :: option
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
#ifdef SURFACE_FLOW
  type(stepper_type), pointer :: surf_flow_stepper
#endif
  
  PetscReal :: target_time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: cumulative_time_steps
  PetscInt :: max_time_step
  PetscBool :: flag
  PetscReal :: tolerance
  type(waypoint_type), pointer :: cur_waypoint

  ! target time will always be dictated by the flow solver, if present
  ! this does not mean that the transport solver cannot take smaller steps
  
  ! this just works through the logic of whether the flow stepper should
  ! be used to set the target time.  If flow is not present, or is to
  ! be run as steady state, use the transport parameters
  flag = PETSC_TRUE
  if (associated(flow_stepper)) then
    if (flow_stepper%run_as_steady_state) then
      flag = PETSC_FALSE
    endif
  else
    flag = PETSC_FALSE
  endif

  ! this flag allows one to take a shorter or slightly larger step than normal
  ! in order to synchronize with the waypoint time
  if (option%match_waypoint) then
    ! if the maximum time step size decreased in the past step, need to set
    ! the time step size to the minimum of the stepper%prev_dt and stepper%dt_max
    if (associated(flow_stepper)) option%flow_dt = min(flow_stepper%prev_dt, &
                                                       flow_stepper%dt_max)
    if (associated(tran_stepper)) option%tran_dt = min(tran_stepper%prev_dt, &
                                                       tran_stepper%dt_max)
    option%match_waypoint = PETSC_FALSE
  endif

  if (flag) then ! flow stepper will govern the target time
!    time = option%flow_time + option%flow_dt
    dt = option%flow_dt
    ! dt_max must be set from current waypoint and not updated below
    cur_waypoint => flow_stepper%cur_waypoint
    dt_max = cur_waypoint%dt_max
    cumulative_time_steps = flow_stepper%steps
    max_time_step = flow_stepper%max_time_step
    tolerance = flow_stepper%time_step_tolerance
    target_time = flow_stepper%target_time + option%flow_dt
  else
!    time = option%tran_time + option%tran_dt
    dt = option%tran_dt
    cur_waypoint => tran_stepper%cur_waypoint
    ! dt_max must be set from current waypoint and not updated below
    dt_max = cur_waypoint%dt_max
    cumulative_time_steps = tran_stepper%steps
    max_time_step = tran_stepper%max_time_step
    tolerance = tran_stepper%time_step_tolerance
    target_time = tran_stepper%target_time + option%tran_dt
  endif
  
  ! For the case where the second waypoint is a printout after the first time step
  ! we must increment the waypoint beyond the first (time=0.) waypoint.  Otherwise
  ! the second time step will be zero. - geh
  if (cur_waypoint%time < 1.d-40) then
    cur_waypoint => cur_waypoint%next
  endif
  
  if (associated(flow_stepper)) flow_stepper%prev_dt = option%flow_dt
  if (associated(tran_stepper)) tran_stepper%prev_dt = option%tran_dt

#ifdef SURFACE_FLOW
  if(associated(surf_flow_stepper)) then
    if (option%surf_subsurf_coupling_flow_dt>0.d0) then

      !if next waypoint is a very close, increase the timestep to match it
      if(cur_waypoint%print_output) then
        if(cur_waypoint%time>target_time.and.(cur_waypoint%time-target_time)/dt<0.2) then
          dt=dt+(cur_waypoint%time-target_time)
          target_time=cur_waypoint%time
        endif
      endif

      ! if target time exceed model coupling time, shrink it back
      if((target_time-option%surf_subsurf_coupling_time) &
        /option%flow_dt>1.d-10) then
        dt=option%surf_subsurf_coupling_time-option%flow_time
        target_time=option%surf_subsurf_coupling_time
      endif
    endif

    if(option%subsurf_surf_coupling==DECOUPLED) then
      target_time=surf_flow_stepper%target_time
      dt = option%surf_flow_dt
    endif
  endif
#endif

! If a waypoint calls for a plot or change in src/sinks, adjust time step
! to match waypoint.
  do ! we cycle just in case the next waypoint is beyond the target_time
    if (target_time + tolerance*dt >= cur_waypoint%time .and. &
        (cur_waypoint%update_conditions .or. &
         cur_waypoint%print_output .or. &
         cur_waypoint%print_tr_output .or. &
         cur_waypoint%final)) then
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on waypoint time
      dt = cur_waypoint%time - target_time
      if (dt > dt_max .and. dabs(dt-dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
        dt = dt_max                    ! error from waypoint%time - time
        target_time = target_time + dt
      else
        target_time = cur_waypoint%time
        if (cur_waypoint%print_output) plot_flag = PETSC_TRUE
        if (cur_waypoint%print_tr_output) transient_plot_flag = PETSC_TRUE
        option%match_waypoint = PETSC_TRUE
        cur_waypoint => cur_waypoint%next
      endif
      exit
    else if (target_time > cur_waypoint%time) then
      cur_waypoint => cur_waypoint%next
    else
      exit
    endif
  enddo
  ! subtract 1 from max_time_steps since we still have to complete the current
  ! time step

  if (cumulative_time_steps >= max_time_step-1) then
    plot_flag = PETSC_TRUE
    nullify(cur_waypoint)
  endif

  ! update maximum time step size to current waypoint value
  if (associated(cur_waypoint)) then
    dt_max = cur_waypoint%dt_max
  endif
  
  ! target time will always be dictated by the flow solver, if present
  if (associated(flow_stepper)) then
    option%flow_dt = dt
    flow_stepper%target_time = target_time
    !geh: dt_max now updated in StepperUpdateDTMax() at the end of 
    !     a time step to avoid premature update if cuts or sub-stepping
    !     for transport occur during the time step.
    !flow_stepper%dt_max = dt_max
    flow_stepper%cur_waypoint => cur_waypoint
  endif
  if (associated(tran_stepper)) then
    if (flag) then ! flow stepper governs the target time
      ! if transport step is within the tolerance of the flow step, let it 
      ! try to reach the target time with a slightly larger step
      if (option%match_waypoint .and. &
          dt <= (1.d0+tolerance)*option%tran_dt) then
        option%tran_dt = dt
      endif
      !geh: Need to ensure that tran_dt <= flow_dt
      if (option%tran_dt > option%flow_dt) option%tran_dt = option%flow_dt
    else ! transport only
      option%tran_dt = dt
    endif
    tran_stepper%target_time = target_time
    !geh: see note on flow_stepper%dt_max above.
    !tran_stepper%dt_max = dt_max
    tran_stepper%cur_waypoint => cur_waypoint
  endif
  
end subroutine StepperSetTargetTimes

! ************************************************************************** !
!> This subroutine sets target time for surface flow.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/31/12
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine StepperSetSurfaceFlowTargetTimes(surf_flow_stepper, &
                                 option,plot_flag, &
                                 transient_plot_flag)

  use Option_module

  implicit none

  type(stepper_type), pointer :: surf_flow_stepper
  type(option_type) :: option
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag

  PetscReal :: target_time
  PetscReal :: dt
  PetscReal :: dt_max

  dt = option%surf_flow_dt
  target_time = surf_flow_stepper%target_time + option%surf_flow_dt
  surf_flow_stepper%target_time = target_time

  ! Cut timestep to avoid going pass the surface-subsurface coupling time
  if(option%subsurf_surf_coupling==SEQ_COUPLED) then
    if((surf_flow_stepper%target_time-option%surf_subsurf_coupling_time) &
        /option%surf_flow_dt>1.d-10) then
      option%surf_flow_dt=option%surf_subsurf_coupling_time-option%surf_flow_time
      surf_flow_stepper%target_time=option%surf_subsurf_coupling_time
    endif
  endif

end subroutine StepperSetSurfaceFlowTargetTimes
#endif

! ************************************************************************** !
!> This subroutine sets target time for model coupling.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 10/31/12
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine SetSurfaceSubsurfaceCouplingTime(flow_stepper,tran_stepper,surf_flow_stepper, &
                                        option,plot_flag,transient_plot_flag)

  use Option_module

  implicit none

  type(stepper_type), pointer :: surf_flow_stepper
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(option_type) :: option
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag

  if(option%surf_subsurf_coupling_flow_dt>0.d0) then
    ! Case-I: Coupling time is specified in the input deck, so use it
    option%surf_subsurf_coupling_time=option%surf_subsurf_coupling_time + &
                                      option%surf_subsurf_coupling_flow_dt

    ! Set new target time for surface model
    call StepperSetSurfaceFlowTargetTimes(surf_flow_stepper,option, &
              plot_flag,transient_plot_flag)

    if (option%subsurf_surf_coupling==SEQ_COUPLED) then
      ! Set new target time for subsurface model
      call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                                surf_flow_stepper, &
                                option,plot_flag, &
                                transient_plot_flag)
    else
      option%io_buffer='Coupling time specified for a DECOUPLED system. '
      call printErrMsg(option)
    endif
  else
    if(option%subsurf_surf_coupling==SEQ_COUPLED) then
      ! Case-II: Coupling time not explicitly specified in the input deck, but
      !          the simulation is for a sequentially coupled surface-subsurface
      !          model. Model coupling will be performed at each subsurface
      !          time step.

      ! Set new target time for subsurface model
      call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                                surf_flow_stepper, &
                                option,plot_flag, &
                                transient_plot_flag)

      ! Set coupling time to be subsurface target time
      option%surf_subsurf_coupling_time=flow_stepper%target_time

      ! Set new target time for surface model
      call StepperSetSurfaceFlowTargetTimes(surf_flow_stepper,option, &
                plot_flag,transient_plot_flag)

    else
      ! Case-III: Coupling time not explicitly specified in the input deck, with
      !          only surface simulation 
      ! Set new target time for surface model
      call StepperSetSurfaceFlowTargetTimes(surf_flow_stepper,option, &
                plot_flag,transient_plot_flag)

      option%surf_subsurf_coupling_time=surf_flow_stepper%target_time

      ! Set new target time for subsurface model
      call StepperSetTargetTimes(flow_stepper,tran_stepper, &
                                surf_flow_stepper, &
                                option,plot_flag, &
                                transient_plot_flag)

    endif
  endif

end subroutine SetSurfaceSubsurfaceCouplingTime
#endif

! ************************************************************************** !
!
! StepperStepFlowDT: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepFlowDT(realization,stepper,step_to_steady_state,failure)

  use Flash2_module, only : Flash2MaxChange, Flash2InitializeTimestep, &
                           Flash2TimeCut, Flash2UpdateReason
  use MPHASE_module, only : MphaseMaxChange, MphaseInitializeTimestep, &
                           MphaseTimeCut, MPhaseUpdateReason
  use Immis_module, only : ImmisMaxChange, ImmisInitializeTimestep, &
                           ImmisTimeCut, ImmisUpdateReason
  use Miscible_module, only : MiscibleMaxChange, MiscibleInitializeTimestep, &
                           MiscibleTimeCut
  use Richards_module, only : RichardsMaxChange, RichardsInitializeTimestep, &
                             RichardsTimeCut, RichardsResidual
  use THC_module, only : THCMaxChange, THCInitializeTimestep, THCTimeCut
  use THMC_module, only : THMCMaxChange, THMCInitializeTimestep, THMCTimeCut

  use General_module, only : GeneralMaxChange, GeneralInitializeTimestep, &
                             GeneralTimeCut
  use Global_module

  use Output_module, only : Output
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module, only : STRUCTURED_GRID_MIMETIC
  
  implicit none

#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscBool :: step_to_steady_state
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason, tmp_int
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm, prev_norm, dif_norm, rel_norm
  PetscReal :: tempreal, tempreal2
  Vec :: update_vec
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(discretization_type), pointer :: discretization 
  type(solver_type), pointer :: solver
  

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver


! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  icut = 0
  sum_newton_iterations = 0
  sum_linear_iterations = 0
  prev_norm = 1.d20

  do ! this loop is for steady state initial condition
  
    ! Perform some global-to-local scatters to update the ghosted vectors.
    ! We have to do this so that the routines for calculating the residual
    ! and the Jacobian will have the ghost points they need.
    ! Note that we don't do the global-to-local scatter for the pressure 
    ! vector, as that needs to be done within the residual calculation routine
    ! because that routine may get called several times during one Newton step
    ! if a method such as line search is being used.
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                    field%iphas_loc,ONEDOF)
    
    if (option%print_screen_flag) then
      if (step_to_steady_state) then
        write(*,'(/,2("=")," Initialize FLOW to Steady-State ",25("="))')
      else
        write(*,'(/,2("=")," FLOW ",52("="))')
      endif
    endif

#ifdef DASVYAT
!     write(*,*) "Flow begins"
!     read(*,*)
#endif

    if (option%ntrandof > 0) then ! store initial saturations for transport
      call GlobalUpdateAuxVars(realization,TIME_T)
    endif
    
    select case(option%iflowmode)
      case(THMC_MODE)
        call THMCInitializeTimestep(realization)
      case(THC_MODE)
        call THCInitializeTimestep(realization)
      case(RICHARDS_MODE)
        call RichardsInitializeTimestep(realization)
      case(MPH_MODE)
        call MphaseInitializeTimestep(realization)
      case(MIS_MODE)
        call MiscibleInitializeTimestep(realization)
      case(IMS_MODE)
        call ImmisInitializeTimestep(realization)
      case(FLASH2_MODE)
        call Flash2InitializeTimestep(realization)
      case(G_MODE)
        call GeneralInitializeTimestep(realization)
    end select

    
    do
      
      call PetscGetTime(log_start_time, ierr)

      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE,THMC_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)
        case(RICHARDS_MODE)
          if (discretization%itype == STRUCTURED_GRID_MIMETIC) then 
            call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx_faces, ierr)
          else 
            call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)
          end if

#if DASVYAT_DEBUG

    call PetscViewerASCIIOpen(realization%option%mycomm,'timestepp_flow_xx_after.out', &
                              viewer,ierr)
    if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
!             call VecView(field%flow_xx_faces, viewer, ierr)
     !       call VecView(field%flow_xx, viewer, ierr)
     !        call VecView(field%flow_r_faces, viewer, ierr)
              call VecNorm(field%flow_r_faces, NORM_2, tempreal, ierr)

             write(*,*) "MFD residual", tempreal
    else
            call VecView(field%flow_xx, viewer, ierr)
    end if

!    stop

    call RichardsResidual(solver%snes,field%flow_xx, field%flow_r,realization,ierr)

    call VecView(field%flow_r, viewer, ierr)
    call VecNorm(field%flow_r, NORM_2, tempreal2, ierr)


    write(*,*) "FV residual", tempreal2

    call PetscViewerDestroy(viewer,ierr)

    write(*,*) "After SNESSolve" 

!    if (tempreal2/tempreal > 1e+4) stop

!     stop
      read(*,*)   
     
#endif



      end select
      call PetscGetTime(log_end_time, ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
                                       (log_end_time - log_start_time)

  ! do we really need all this? - geh 
      call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
      call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
      call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

      sum_newton_iterations = sum_newton_iterations + num_newton_iterations
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      update_reason = 1
      
      if (snes_reason >= 0) then
        select case(option%iflowmode)
          case(IMS_MODE)
            call ImmisUpdateReason(update_reason,realization)
          case(MPH_MODE)
!           call MPhaseUpdateReason(update_reason,realization)
            if (option%use_mc) then
              option%sec_vars_update = PETSC_TRUE
            endif
          case(FLASH2_MODE)
!           call Flash2UpdateReason(update_reason,realization)
          case(THC_MODE)
            update_reason=1
            if (option%use_mc) then
              option%sec_vars_update = PETSC_TRUE
            endif
          case(THMC_MODE)
            update_reason=1
          case (MIS_MODE)
            update_reason=1
          case(RICHARDS_MODE,G_MODE)
            update_reason=1
         end select   
        !if (option%print_screen_flag) print *,'update_reason: ',update_reason
      endif
   
  !******************************************************************
      
      if (snes_reason <= 0 .or. update_reason <= 0) then
        ! The Newton solver diverged, so try reducing the time step.
        icut = icut + 1
        stepper%time_step_cut_flag = PETSC_TRUE
        option%out_of_table = PETSC_FALSE

        if (icut > stepper%max_time_step_cuts .or. option%flow_dt<1.d-20) then
          if (option%print_screen_flag) then
            print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                    stepper%max_time_step_cuts, "t= ", &
                    stepper%target_time/realization%output_option%tconv, " dt= ", &
                    option%flow_dt/realization%output_option%tconv
            print *,"Stopping execution!"
          endif
          realization%output_option%plot_name = 'flow_cut_to_failure'
          plot_flag = PETSC_TRUE
          transient_plot_flag = PETSC_FALSE
          call Output(realization,plot_flag,transient_plot_flag)
          failure = PETSC_TRUE
          return
        endif

 
       stepper%target_time = stepper%target_time - option%flow_dt

#ifdef DASVYAT_TEST_CUT
        option%flow_dt = 1d0*option%flow_dt
#else
        option%flow_dt = 0.5d0 * option%flow_dt  
#endif
      
        if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
          &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
          &   1pe12.5)')  snes_reason,icut,stepper%cumulative_time_step_cuts, &
              option%flow_time/realization%output_option%tconv, &
              option%flow_dt/realization%output_option%tconv

        stepper%target_time = stepper%target_time + option%flow_dt

        select case(option%iflowmode)
          case(THC_MODE)
            call THCTimeCut(realization)
          case(THMC_MODE)
            call THMCTimeCut(realization)
          case(RICHARDS_MODE)
            call RichardsTimeCut(realization)
          case(MPH_MODE)
            call MphaseTimeCut(realization)
          case(MIS_MODE)
            call MiscibleTimeCut(realization)
          case(IMS_MODE)
            call ImmisTimeCut(realization)
          case(FLASH2_MODE)
            call Flash2TimeCut(realization)
          case(G_MODE)
            call GeneralTimeCut(realization)
        end select
        call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr)

      else
        ! The Newton solver converged, so we can exit.
        exit
      endif
    enddo
    
    if (.not.step_to_steady_state) then
      exit
    else
    
      call SNESGetSolutionUpdate(solver%snes,update_vec,ierr)
      call VecStrideNorm(update_vec,ZERO_INTEGER,NORM_INFINITY,inorm,ierr)
      dif_norm = inorm-prev_norm
      rel_norm = dif_norm/prev_norm
      if (sum_newton_iterations > 20 .and. &
!geh: These norms don't seem to be the best, sticking with inorm
!           dabs(dif_norm) < stepper%steady_state_rel_tol) .or. &
!          dabs(rel_norm) < stepper%steady_state_rel_tol) then
          inorm < stepper%steady_state_rel_tol) then
        if (option%print_file_flag) then
          write(option%fid_out,*) 'Steady state solve converged after ', &
            sum_newton_iterations, ' iterations'
        endif
        exit
      endif
      
      prev_norm = inorm
      
      if (sum_newton_iterations > 1000) then
        option%io_buffer = 'Steady-state solve exceeded 1000 Newton iterations.'
        call printMsg(option)
        if (option%print_file_flag) then
          write(option%fid_out,*) trim(option%io_buffer)
          write(option%fid_out,*) 'inorm: ', inorm
          write(option%fid_out,*) 'prev_norm: ', prev_norm
          write(option%fid_out,*) 'rel_norm: ', rel_norm
          write(option%fid_out,*) 'dif_norm: ', dif_norm
          write(option%fid_out,*) 'relative tolerance: ', &
                                  stepper%steady_state_rel_tol
        endif
        failure = PETSC_TRUE
        exit
      endif

      ! zero out dt to prevent potential error in mass balance calc while 
      ! updating solution
      tempreal = option%flow_dt
      option%flow_dt = 0.d0  
      ! take next step at larger dt      
      call StepperUpdateFlowSolution(realization)
      option%flow_dt = tempreal
      if (num_newton_iterations < 4) then
        option%flow_dt = option%flow_dt*2.d0
      else if (num_newton_iterations < 5) then
        option%flow_dt = option%flow_dt*1.5d0
      else if (num_newton_iterations < 6) then
        option%flow_dt = option%flow_dt*1.d0
      else if (num_newton_iterations < 7) then
        option%flow_dt = option%flow_dt*0.9d0
      else if (num_newton_iterations < 8) then
        option%flow_dt = option%flow_dt*0.8d0
      else if (num_newton_iterations < 9) then
        option%flow_dt = option%flow_dt*0.7d0
      else
        option%flow_dt = option%flow_dt*0.6d0
      endif

      if (option%print_file_flag) then
        write(option%fid_out,*) 'Dt: ', option%flow_dt
        write(option%fid_out,*) 'Inf Norm: ', inorm
        write(option%fid_out,*) 'Relative Norm: ', rel_norm
        write(option%fid_out,*) 'Linear its/Newton It: ', &
          float(num_linear_iterations)/num_newton_iterations
      endif
      if (option%print_screen_flag) then
        write(*,*) 'Dt: ', option%flow_dt
        write(*,*) 'Inf Norm: ', inorm
        write(*,*) 'Relative Norm: ', rel_norm
        write(*,*) 'Linear its/Newton It: ', &
          float(num_linear_iterations)/num_newton_iterations
      endif
    endif

  enddo ! end of steady-state loop
  
  if (.not.step_to_steady_state) then
    stepper%steps = stepper%steps + 1      
  endif
  stepper%cumulative_newton_iterations = &
    stepper%cumulative_newton_iterations + sum_newton_iterations
  stepper%cumulative_linear_iterations = &
    stepper%cumulative_linear_iterations + sum_linear_iterations
  stepper%cumulative_time_step_cuts = &
    stepper%cumulative_time_step_cuts + icut

  stepper%num_newton_iterations = num_newton_iterations
  stepper%num_linear_iterations = num_linear_iterations

  if (option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(realization,TIME_TpDT)
  endif
    
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts

    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  if (option%print_file_flag) then
    write(option%fid_out, '(" FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts
  endif
  
  select case(option%iflowmode)
    case(THC_MODE)
      call THCMaxChange(realization)
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax, option%dcmax
      endif
    case(THMC_MODE)
      call THMCMaxChange(realization)
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax, option%dcmax
      endif
      if (option%print_file_flag) then 
        write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax
      endif
    case(RICHARDS_MODE)
      call RichardsMaxChange(realization)
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
      endif
      if (option%print_file_flag) then
        write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
      endif
    case(MPH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
      select case(option%iflowmode)
        case(MPH_MODE)
          call MphaseMaxChange(realization)
        case(IMS_MODE)
          call ImmisMaxChange(realization)
        case(MIS_MODE)
          call MiscibleMaxChange(realization)
        case(FLASH2_MODE)
          call FLASH2MaxChange(realization)
        case(G_MODE)
          call GeneralMaxChange(realization)
      end select
      ! note use mph will use variable switching, the x and s change is not meaningful 
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
      endif
      if (option%print_file_flag) then  
        write(option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax,option%dcmax,option%dsmax
      endif
  end select


#ifdef DASVYAT_DEBUG
    write(*,*) "End FLOW" 
    read(*,*)    
#endif

  if (option%print_screen_flag) print *, ""
  
  ! option%flow_time is updated outside this subroutine

end subroutine StepperStepFlowDT

! ************************************************************************** !
!
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine StepperStepSurfaceFlowDT(surf_realization,stepper,failure)
  
  use Surface_Realization_module
  use Surface_Flow_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Surface_Field_module
  use Grid_module
  use Output_module, only : Output
  
  implicit none
  
#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(surface_realization_type) :: surf_realization
  type(stepper_type)     :: stepper
  PetscBool              :: failure
  
  PetscInt :: icut ! Tracks the number of time step reductions applied
  PetscInt :: update_reason, tmp_int
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm, prev_norm, dif_norm, rel_norm
  SNESConvergedReason :: snes_reason 
  PetscErrorCode :: ierr
  !PetscLogDouble :: start_time, end_time
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer          :: option
  type(surface_field_type), pointer   :: surf_field 
  type(discretization_type), pointer  :: discretization 
  type(solver_type), pointer          :: solver
  character(len=MAXSTRINGLENGTH)      :: string

  option         => surf_realization%option
  discretization => surf_realization%discretization
  surf_field     => surf_realization%surf_field
  solver         => stepper%solver

  icut = 0
  sum_newton_iterations = 0
  sum_linear_iterations = 0
  prev_norm = 1.d20

  if (option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE_FLOW ",52("="))')
  endif

  select case(option%iflowmode)
    case (RICHARDS_MODE)
      call SurfaceFlowInitializeTimestep(surf_realization)
    case default
      option%io_buffer = 'ERROR: Incorrect iflowmode in SurfaceFlow'
      call printErrMsgByRank(option)
  end select

  do
    call PetscGetTime(log_start_time,ierr)
    
    select case(option%iflowmode)
      case (RICHARDS_MODE)
        call SNESSolve(solver%snes,PETSC_NULL_OBJECT,surf_field%flow_xx,ierr)
      case default
        option%io_buffer = 'ERROR: Incorrect iflowmode in SurfaceFlow'
        call printErrMsgByRank(option)
    end select

    call PetscGetTime(log_end_time,ierr)


    stepper%cumulative_solver_time_surf_flow =  &
                                  stepper%cumulative_solver_time_surf_flow + &
                                  (log_end_time - log_start_time)
    call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
    call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
    update_reason = 1

    if (snes_reason >= 0) then
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          update_reason=1
      end select
    !if (option%print_screen_flag) print *,'update_reason: ',update_reason
    endif

    if (snes_reason <= 0 .or. update_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut=icut+1
      stepper%time_step_cut_flag=PETSC_TRUE

      if (icut > stepper%max_time_step_cuts .or. option%flow_dt<1.d-20) then
        if (option%print_screen_flag) then
          print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                  stepper%max_time_step_cuts, "t= ", &
                  stepper%target_time/surf_realization%output_option%tconv, " dt= ", &
                  option%flow_dt/surf_realization%output_option%tconv
          print *,"Stopping execution!"
        endif
        surf_realization%output_option%plot_name = 'flow_cut_to_failure'
        !call Output(surf_realization,PETSC_TRUE,PETSC_FALSE)
        failure = PETSC_TRUE
        return
      endif

      ! Revert back the target time and decrease the time step
      stepper%target_time = stepper%target_time - option%flow_dt
      option%surf_flow_dt = 0.5d0 * option%surf_flow_dt

      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,stepper%cumulative_time_step_cuts, &
            option%surf_flow_time/surf_realization%output_option%tconv, &
            option%surf_flow_dt/surf_realization%output_option%tconv

      ! Set new target time
      stepper%target_time = stepper%target_time + option%surf_flow_dt

      select case(option%iflowmode)
        case(RICHARDS_MODE)
          call SurfaceFlowTimeCut(surf_realization)
        case default
          option%io_buffer = 'ERROR: TimeCut for this iflowmode in SurfaceFlow ' // &
            ' not incorporated.'
          call printErrMsg(option)
      end select
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  stepper%steps_surf_flow = stepper%steps_surf_flow + 1
  stepper%cumulative_newton_iterations_surf_flow = &
    stepper%cumulative_newton_iterations_surf_flow + sum_newton_iterations
  stepper%cumulative_linear_iterations_surf_flow = &
    stepper%cumulative_linear_iterations_surf_flow + sum_linear_iterations
  stepper%cumulative_time_step_cuts_surf_flow = &
    stepper%cumulative_time_step_cuts_surf_flow + icut

  stepper%num_newton_iterations_surf_flow = num_newton_iterations
  stepper%num_linear_iterations_surf_flow = num_linear_iterations

! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(surf_field%flow_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," SURFACE FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps_surf_flow, &
      !stepper%target_time/surf_realization%output_option%tconv, &
      option%surf_flow_time+option%surf_flow_dt, &
      option%surf_flow_dt/surf_realization%output_option%tconv, &
      surf_realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations_surf_flow,sum_linear_iterations, &
      stepper%cumulative_linear_iterations_surf_flow,icut, &
      stepper%cumulative_time_step_cuts_surf_flow

    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  select case(option%iflowmode)
    case(RICHARDS_MODE)
    call SurfaceFlowMaxChange(surf_realization)
    if (option%print_screen_flag) then
      write(*,'("  --> max chng: dpmx= ",1pe12.4)') option%dpmax
    endif
  end select

  if (option%print_screen_flag) print *, ""

end subroutine StepperStepSurfaceFlowDT
#endif

! ************************************************************************** !
!
! StepperStepTransportDT_GI: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepTransportDT_GI(realization,stepper, &
                                     steady_flow,flow_t0,flow_t1, &
                                     failure)
  
  use Reactive_Transport_module
  use Output_module, only : Output
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module
  use Level_module
  use Patch_module
  use Global_module  
  
  implicit none

#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscReal :: flow_t0, flow_t1
  PetscBool :: steady_flow
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscInt :: n, nmax_inf
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscBool :: plot_flag  
  PetscBool :: transient_plot_flag  
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)
  PetscReal :: final_tran_time
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver
  type(patch_type), pointer :: cur_patch
  type(level_type), pointer :: cur_level

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

#ifdef DASVYAT
!write(*,*) "Beginning of StepperStepTransportDT"
!read(*,*)
!stop
#endif

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                  field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                  field%tortuosity_loc,ONEDOF)

  ! interpolate flow parameters/data
  ! this must remain here as these weighted values are used by both
  ! RTInitializeTimestep and RTTimeCut (which calls RTInitializeTimestep)
  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call TimestepperSetTranWeights(option,flow_t0, flow_t1)
    ! set densities and saturations to t
    call GlobalUpdateDenAndSat(realization,option%tran_weight_t0)
  endif

  call RTInitializeTimestep(realization)
  ! note: RTUpdateTransportCoefs() is called within RTInitializeTimestep()

  ! set densities and saturations to t+dt
  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call GlobalUpdateDenAndSat(realization,option%tran_weight_t1)
  endif

  call RTUpdateTransportCoefs(realization)
  
  sum_newton_iterations = 0
  sum_linear_iterations = 0
  icut = 0

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT ",47("="))')

  do
   
    final_tran_time = option%tran_time + option%tran_dt
   
    if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
      call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
!       The below is set within RTUpdateAuxVarsPatch() when PETSC_TRUE,PETSC_TRUE,* are passed
!       patch%aux%RT%aux_vars_up_to_date = PETSC_TRUE 
    endif
    if (realization%reaction%use_log_formulation) then
      call VecCopy(field%tran_xx,field%tran_log_xx,ierr)
      call VecLog(field%tran_log_xx,ierr)

      call PetscGetTime(log_start_time, ierr)
      call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_log_xx, ierr)
      call PetscGetTime(log_end_time, ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
        (log_end_time - log_start_time)          
        
      call VecCopy(field%tran_log_xx,field%tran_xx,ierr)
      call VecExp(field%tran_xx,ierr)
    else
      call PetscGetTime(log_start_time, ierr)
      call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_xx, ierr)
      call PetscGetTime(log_end_time, ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
        (log_end_time - log_start_time)          
    endif

! do we really need all this? - geh 
    call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
    call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      
    if (snes_reason <= 0) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      stepper%time_step_cut_flag = PETSC_TRUE

      if (icut > stepper%max_time_step_cuts .or. option%tran_dt<1.d-20) then
        if (option%print_screen_flag) then
          print *,"--> max_time_step_cuts exceeded: icut/icutmax= ", &
                  icut,stepper%max_time_step_cuts, &
                  "t= ",final_tran_time/realization%output_option%tconv, " dt= ", &
                  option%tran_dt/realization%output_option%tconv
          print *,"Stopping execution!"
        endif
        realization%output_option%plot_name = 'tran_cut_to_failure'
        plot_flag = PETSC_TRUE
        transient_plot_flag = PETSC_FALSE
        call Output(realization,plot_flag,transient_plot_flag)
        failure = PETSC_TRUE
        return
      endif

      option%tran_dt = 0.5d0 * option%tran_dt
    
      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,stepper%cumulative_time_step_cuts, &
            option%tran_time/realization%output_option%tconv, &
            option%tran_dt/realization%output_option%tconv


      ! recompute weights
      if (option%nflowdof > 0 .and. .not.steady_flow) then
        call TimestepperSetTranWeights(option,flow_t0, flow_t1)
      endif
      call RTTimeCut(realization)

    else
      ! increment time step number
      stepper%steps = stepper%steps + 1      
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  ! these are the number of iterations that it takes to
  ! solve the time step
  stepper%cumulative_newton_iterations = &
    stepper%cumulative_newton_iterations + sum_newton_iterations
  stepper%cumulative_linear_iterations = &
    stepper%cumulative_linear_iterations + sum_linear_iterations
  stepper%cumulative_time_step_cuts = stepper%cumulative_time_step_cuts + icut

  ! these are the number of iterations that it takes to
  ! achieve convergence
  stepper%num_newton_iterations = num_newton_iterations
  stepper%num_linear_iterations = num_linear_iterations

  ! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
  
    if (option%nflowdof > 0 .and. .not.steady_flow) then

    write(*, '(/," TRAN ",i6," Time= ",1pe12.5," Target= ",1pe12.5, &
      & " Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      flow_t1/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts

    else

    write(*, '(/," TRAN ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts

    endif
    
    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax   
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" TRAN ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a1,"]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      realization%output_option%tunit,snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts
  endif
  
  call RTMaxChange(realization)
  if (option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4," [mol/s]")') &
      option%dcmax,option%dcmax/option%tran_dt
  endif
  if (option%print_file_flag) then  
    write(option%fid_out,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ", &
      & 1pe12.4," [mol/s]")') &
      option%dcmax,option%dcmax/option%tran_dt
  endif

  if (option%print_screen_flag) print *, ""
 
  ! option%tran_time is updated outside this subroutine

end subroutine StepperStepTransportDT_GI
! ************************************************************************** !
!
! StepperStepTransportDT_OS: Steps forward one step in time (operator split)
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperStepTransportDT_OS(realization,stepper, &
                                     steady_flow, flow_t0,flow_t1, &
                                     failure)

  use Reactive_Transport_module, only : RTUpdateRHSCoefs, RTUpdateAuxVars, &
        RTCalculateRHS_t0, RTUpdateTransportCoefs, RTCalculateRHS_t1, &
        RTCalculateTransportMatrix, RTReact, RTMaxChange, RTExplicitAdvection
  use Output_module, only : Output
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module
  use Level_module
  use Patch_module
  use Global_module  

  implicit none

#include "finclude/petsclog.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscReal :: flow_t0, flow_t1
  PetscBool :: steady_flow
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  KSPConvergedReason :: ksp_reason 
  PetscInt :: sum_linear_iterations, num_linear_iterations
  PetscInt :: idof
  PetscBool :: plot_flag  
  PetscBool :: transient_plot_flag  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: inf_norm, euclid_norm
  PetscReal :: final_tran_time
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

  call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                  field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                  field%tortuosity_loc,ONEDOF)

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT (OS) ",43("="))')

  final_tran_time = option%tran_time + option%tran_dt
   
  sum_linear_iterations = 0
  ksp_reason = 0
  
  ! do we need this...don't think so - geh
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

  call PetscGetTime(log_start_time, ierr)

  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call TimestepperSetTranWeights(option,flow_t0,flow_t1)
    ! set densities and saturations to t
    call GlobalUpdateDenAndSat(realization,option%tran_weight_t0)
  endif

  ! update time derivative on RHS
  call RTUpdateRHSCoefs(realization)
  if (option%itranmode == EXPLICIT_ADVECTION) then
    ! update the totals (must be at t0)
                                   ! cells      bcs         act. coefs.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)
    ! set densities and saturations to t+dt
    if (option%nflowdof > 0 .and. .not.steady_flow) then
      call GlobalUpdateDenAndSat(realization,option%tran_weight_t1)
    endif
    call RTExplicitAdvection(realization)
  else
  
    ! Between the next 3 subroutine calls:
    ! RTUpdateRHSCoefs() -- now calculated prior to conditional
    ! RTUpdateAuxVars()
    ! RTCalculateRHS_t0()
    ! the t0 portion of the accumulation term is calculated for the RHS vector
  
    ! calculate total component concentrations based on t0 densities
                                   ! cells      bcs         act. coefs.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)
    call RTCalculateRHS_t0(realization)
    
    ! set densities and saturations to t+dt
    if (option%nflowdof > 0 .and. .not.steady_flow) then
      call GlobalUpdateDenAndSat(realization,option%tran_weight_t1)
    endif

    ! update diffusion/dispersion coefficients
    call RTUpdateTransportCoefs(realization)
    ! RTCalculateRHS_t1() updates aux vars (cell and boundary) to k+1 and calculates 
    ! RHS fluxes and src/sinks
    call RTCalculateRHS_t1(realization)

    if (realization%debug%vecview_residual) then
      call PetscViewerASCIIOpen(realization%option%mycomm,'Trhs.out', &
                                viewer,ierr)
      call VecView(field%tran_rhs,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)
    endif

    ! RTCalculateTransportMatrix() calculates flux coefficients and the
    ! t^(k+1) coefficient in accumulation term
    call RTCalculateTransportMatrix(realization,solver%J) 
    call KSPSetOperators(solver%ksp,solver%J,solver%Jpre, &
                         SAME_NONZERO_PATTERN,ierr)

  !  call VecGetArrayF90(field%tran_xx,vec_ptr,ierr)
  !  call VecRestoreArrayF90(field%tran_xx,vec_ptr,ierr)

    ! loop over chemical component and transport
    do idof = 1, option%ntrandof

  ! for debugging
#if 0    
      call RealizationGetDataset(realization,field%work,TOTAL_MOLARITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr)
  
      call RealizationGetDataset(realization,field%work,TOTAL_MOLALITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr)
  
      call RealizationGetDataset(realization,field%work,PRIMARY_MOLALITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr)
#endif
    
      call VecStrideGather(field%tran_rhs,idof-1,field%work,INSERT_VALUES,ierr)
      option%rt_idof = idof
      call KSPSolve(solver%ksp,field%work,field%work,ierr)
      ! tran_xx will contain transported totals
      ! tran_xx_loc will still contain free-ion from previous solution
      call VecStrideScatter(field%work,idof-1,field%tran_xx,INSERT_VALUES,ierr)

  ! for debugging
#if 0
      call VecGetArrayF90(field%work,vec_ptr,ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr)
#endif      

      call KSPGetIterationNumber(solver%ksp,num_linear_iterations,ierr)
      call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
    enddo
    sum_linear_iterations = int(dble(sum_linear_iterations) / option%ntrandof)
    stepper%cumulative_linear_iterations = &
      stepper%cumulative_linear_iterations + sum_linear_iterations
  endif ! if (EXPLICIT_ADVECTION)

  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Txx.out', &
                              viewer,ierr)
    call VecView(field%tran_xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  ! activity coefficients are updated within RReact!  DO NOT updated
  ! here as doing so will cause errors in the t0 portion of the
  ! accumulation term for equilibrium sorbed species
  call RTReact(realization)
  
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RTxx.out', &
                              viewer,ierr)
    call VecView(field%tran_xx,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
  endif

  call PetscBarrier(solver%ksp,ierr)
  call PetscGetTime(log_end_time, ierr)
  stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
                                   (log_end_time - log_start_time)          

  stepper%num_newton_iterations = 1
  stepper%num_linear_iterations = sum_linear_iterations

  ! increment time step number
  stepper%steps = stepper%steps + 1     
  
  if (option%print_screen_flag) then
    write(*, '(" TRAN ",i6," Time= ",1pe12.5," Dt= ", &
          & 1pe12.5," [",a1,"]"," ksp_conv_reason: ",i4,/," linear = ",i5, &
          & " [",i10,"]")') stepper%steps, &
!geh        option%tran_time/realization%output_option%tconv, &
        final_tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,ksp_reason,sum_linear_iterations, &
        stepper%cumulative_linear_iterations
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" TRAN ",i6," Time= ",1pe12.5," Dt= ", &
          & 1pe12.5," [",a1,"]"," ksp_conv_reason = ",i4,/," linear = ",i5, &
          & " [",i10,"]")') stepper%steps, &
!geh        option%tran_time/realization%output_option%tconv, &
        final_tran_time/realization%output_option%tconv, &
        option%tran_dt/realization%output_option%tconv, &
        realization%output_option%tunit,ksp_reason,sum_linear_iterations, &
        stepper%cumulative_linear_iterations
  endif

  call RTMaxChange(realization)
  if (option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4," [mol/s]")') &
      option%dcmax,option%dcmax/option%tran_dt
  endif
  if (option%print_file_flag) then  
    write(option%fid_out,'("  --> max chng: dcmx= ",1pe12.4," dc/dt= ",1pe12.4," [mol/s]")') &
          option%dcmax,option%dcmax/option%tran_dt
  endif

  if (option%print_screen_flag) print *, ""  

end subroutine StepperStepTransportDT_OS

! ************************************************************************** !
!
! StepperRunSteadyState: Solves steady state solution for flow and transport
! author: Glenn Hammond
! date: 03/10/09
!
! ************************************************************************** !
subroutine StepperRunSteadyState(realization,flow_stepper,tran_stepper)

  use Realization_module

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit, OutputVectorTecplot
  use Logging_module
  use Discretization_module

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper

  PetscBool :: transient_plot_flag
  PetscBool :: plot_flag
  PetscBool :: failure
  PetscLogDouble :: start_time, end_time
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  failure = PETSC_FALSE

  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)

  ! print initial condition output if not a restarted sim
  if (associated(flow_stepper)) then
    call OutputInit(realization,flow_stepper%steps)
  else
    call OutputInit(realization,tran_stepper%steps)
  endif
  transient_plot_flag = PETSC_FALSE
  plot_flag = PETSC_TRUE
  if (output_option%print_initial) then
    transient_plot_flag = PETSC_FALSE
    call Output(realization,plot_flag,transient_plot_flag)
  endif

  if (OptionPrintToScreen(option)) then
    option%print_screen_flag = PETSC_TRUE
  else
    option%print_screen_flag = PETSC_FALSE
  endif

  if (OptionPrintToFile(option)) then
    option%print_file_flag = PETSC_TRUE
  else
    option%print_file_flag = PETSC_FALSE
  endif

  plot_flag = PETSC_TRUE
    
  if (associated(flow_stepper)) then
    call PetscGetTime(start_time, ierr)
    call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr)
    call StepperSolveFlowSteadyState(realization,flow_stepper,failure)
    call PetscLogStagePop(ierr)
    call PetscGetTime(end_time, ierr)
    if (OptionPrintToScreen(option)) then
      write(*, &
         &  '(/,1pe12.4," seconds to solve steady state flow problem",/)') &
        end_time-start_time
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out, &
         &  '(/,1pe12.4," seconds to solve steady state flow problem",/)') &
        end_time-start_time
    endif
    if (failure) return ! if flow solve fails, exit
  endif

  if (associated(tran_stepper)) then
    call PetscGetTime(start_time, ierr)
    call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr)
    call StepperSolveTranSteadyState(realization,tran_stepper,failure)
    call PetscLogStagePop(ierr)
    call PetscGetTime(end_time, ierr)
    if (OptionPrintToScreen(option)) then
      write(*, &
         &'(/,1pe12.4," seconds to solve steady state transport problem",/)') &
        end_time-start_time
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out, &
         &'(/,1pe12.4," seconds to solve steady state transport problem",/)') &
        end_time-start_time
    endif
    if (failure) return ! if transport solve fails, exit
  endif

  if (output_option%print_initial) then
    output_option%plot_number = 1
    call Output(realization,plot_flag,transient_plot_flag)
  endif
   
  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_stepper,tran_stepper, &
                           NEG_ONE_INTEGER)  
  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_stepper%cumulative_newton_iterations,flow_stepper%cumulative_linear_iterations
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_stepper%cumulative_newton_iterations,tran_stepper%cumulative_linear_iterations
    endif            
  endif

  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_stepper%cumulative_newton_iterations, &
            flow_stepper%cumulative_linear_iterations
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_stepper%cumulative_newton_iterations, &
            tran_stepper%cumulative_linear_iterations
    endif            
  endif

  call PetscLogStagePop(ierr)


end subroutine StepperRunSteadyState
 
! ************************************************************************** !
!
! StepperSolveFlowSteadyState: Solves the steady-state flow equation
! author: Glenn Hammond
! date: 03/10/09
!
! ************************************************************************** !
subroutine StepperSolveFlowSteadyState(realization,stepper,failure)

  use Global_module, only : GlobalUpdateAuxVars
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscBool :: failure

  PetscErrorCode :: ierr
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: snes_reason
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(solver_type), pointer :: solver
  type(discretization_type), pointer :: discretization 

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

  num_newton_iterations = 0
  num_linear_iterations = 0


  call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                  field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                  field%tortuosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                  field%icap_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                  field%ithrm_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
    
  if (option%print_screen_flag) write(*,'(/,2("=")," FLOW (STEADY STATE) ",37("="))')

  call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%flow_xx, ierr)

  call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in FLOW, reason: ', snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif
  
  stepper%cumulative_newton_iterations = num_newton_iterations
  stepper%cumulative_linear_iterations = num_linear_iterations

  if (option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(realization,TIME_T)
    call GlobalUpdateAuxVars(realization,TIME_TpDT)
  endif
    
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  if (option%print_screen_flag) print *, ""

end subroutine StepperSolveFlowSteadyState

! ************************************************************************** !
!
! StepperSolveTranSteadyState: Steps forward one step in time
! author: Glenn Hammond
! date: 02/19/08
!
! ************************************************************************** !
subroutine StepperSolveTranSteadyState(realization,stepper,failure)
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  use Patch_module
  use Level_module
  use Grid_module
    
  use Global_module, only : GlobalUpdateDenAndSat
  use Reactive_Transport_module, only : RTUpdateAuxVars  

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  type(realization_type) :: realization
  type(stepper_type) :: stepper
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  SNESConvergedReason :: snes_reason 
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver
  type(patch_type), pointer :: cur_patch
  type(level_type), pointer :: cur_level

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                  field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                  field%tortuosity_loc,ONEDOF)

  call GlobalUpdateDenAndSat(realization,1.d0)
  num_newton_iterations = 0
  num_linear_iterations = 0

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT (STEADY STATE) ",32("="))')

  if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
  endif

  if (realization%reaction%use_log_formulation) then
    call VecCopy(field%tran_xx,field%tran_log_xx,ierr)
    call VecLog(field%tran_log_xx,ierr)
    call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_log_xx, ierr)
    call VecCopy(field%tran_log_xx,field%tran_xx,ierr)
    call VecExp(field%tran_xx,ierr)
  else
    call SNESSolve(solver%snes, PETSC_NULL_OBJECT, field%tran_xx, ierr)
  endif

  ! do we really need all this? - geh 
  call SNESGetIterationNumber(solver%snes,num_newton_iterations, ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in TRAN, reason: ', snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif
      
  stepper%cumulative_newton_iterations = num_newton_iterations
  stepper%cumulative_linear_iterations = num_linear_iterations

  ! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax   
    else
       scaled_fnorm = fnorm
    endif
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif

  if (option%print_screen_flag) print *, ""

end subroutine StepperSolveTranSteadyState

! ************************************************************************** !
!
! StepperUpdateSolution: Updates the solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine StepperUpdateSolution(realization,surf_realization)
#else
subroutine StepperUpdateSolution(realization)
#endif

  use Realization_module
  use Option_module
#ifdef SURFACE_FLOW
  use Surface_Realization_module
#endif

  implicit none
  
  type(realization_type) :: realization
#ifdef SURFACE_FLOW
  type(surface_realization_type) :: surf_realization
#endif
  
  ! update solution variables
  call RealizationUpdate(realization)
  if (realization%option%nflowdof > 0) &
    call StepperUpdateFlowSolution(realization)
  if (realization%option%ntrandof > 0) &
    call StepperUpdateTransportSolution(realization)

#ifdef SURFACE_FLOW
  call SurfaceRealizationUpdate(surf_realization)
  if (surf_realization%option%nsurfflowdof > 0) &
    call StepperUpdateSurfaceFlowSolution(surf_realization)
#endif
    
end subroutine StepperUpdateSolution

! ************************************************************************** !
!
! StepperUpdateFlowSolution: Updates the flow solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
subroutine StepperUpdateFlowSolution(realization)
  
  use Flash2_module, only: Flash2UpdateSolution
  use MPHASE_module, only: MphaseUpdateSolution
  use Immis_module, only: ImmisUpdateSolution
  use Miscible_module, only: MiscibleUpdateSolution 
  use Richards_module, only : RichardsUpdateSolution
  use THC_module, only : THCUpdateSolution
  use THMC_module, only : THMCUpdateSolution
  use General_module, only : GeneralUpdateSolution

  use Realization_module
  use Option_module

  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr
  
  option => realization%option
  
  select case(option%iflowmode)
    case(MPH_MODE)
      call MphaseUpdateSolution(realization)
    case(IMS_MODE)
      call ImmisUpdateSolution(realization)
    case(MIS_MODE)
      call MiscibleUpdateSolution(realization)
    case(FLASH2_MODE)
      call Flash2UpdateSolution(realization)
    case(THC_MODE)
      call THCUpdateSolution(realization)
    case(THMC_MODE)
      call THMCUpdateSolution(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateSolution(realization)
    case(G_MODE)
      call GeneralUpdateSolution(realization)
  end select    

end subroutine StepperUpdateFlowSolution

#ifdef SURFACE_FLOW
! ************************************************************************** !
!> This subroutine updates the surface flow solution variables
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/22/12
! ************************************************************************** !
subroutine StepperUpdateSurfaceFlowSolution(surf_realization)

  use Surface_Flow_module

  use Surface_Realization_module
  use Option_module

  implicit none

  type(surface_realization_type) :: surf_realization

  type(option_type), pointer :: option

  PetscErrorCode :: ierr

  option => surf_realization%option

  select case(option%iflowmode)
    case(RICHARDS_MODE)
      call SurfaceFlowUpdateSolution(surf_realization)
  end select

end subroutine StepperUpdateSurfaceFlowSolution
#endif


! ************************************************************************** !
!
! StepperUpdateTransportSolution: Updates the transport solution variables
! author: Glenn Hammond
! date: 02/19/08 
!
! ************************************************************************** !
subroutine StepperUpdateTransportSolution(realization)

  use Realization_module
  use Reactive_Transport_module, only : RTUpdateSolution

  implicit none

  type(realization_type) :: realization

  PetscErrorCode :: ierr
  
  call RTUpdateSolution(realization)
  if (realization%reaction%update_porosity .or. &
      realization%reaction%update_tortuosity .or. &
      realization%reaction%update_permeability .or. &
      realization%reaction%update_mineral_surface_area) then
    call RealizationUpdateProperties(realization)
  endif

end subroutine StepperUpdateTransportSolution

! ************************************************************************** !
!
! StepperJumpStart: Sets kinetic sorbed concentrations
! author: Glenn Hammond
! date: 08/05/09 
!
! ************************************************************************** !
subroutine StepperJumpStart(realization)

  use Realization_module
  use Reactive_Transport_module, only : RTJumpStartKineticSorption

  implicit none

  type(realization_type) :: realization

  call RTJumpStartKineticSorption(realization)

end subroutine StepperJumpStart

! ************************************************************************** !
!
! StepperUpdateFlowAuxVars: Updates the flow auxilliary variables
! author: Glenn Hammond
! date: 10/11/08 
!
! ************************************************************************** !
subroutine StepperUpdateFlowAuxVars(realization)
  
  use Flash2_module, only: Flash2UpdateAuxVars
  use MPHASE_module, only: MphaseUpdateAuxVars
  use Immis_module, only: ImmisUpdateAuxVars
  use Miscible_module, only: MiscibleUpdateAuxVars
  use Richards_module, only : RichardsUpdateAuxVars
  use THC_module, only : THCUpdateAuxVars
  use THMC_module, only : THMCUpdateAuxVars
  use General_module, only : GeneralUpdateAuxVars

  use Realization_module
  use Option_module

  implicit none

  type(realization_type) :: realization

  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr
  
  option => realization%option
  
  select case(option%iflowmode)
    case(FLASH2_MODE)
      call Flash2UpdateAuxVars(realization)
    case(IMS_MODE)
      call ImmisUpdateAuxVars(realization)
    case(MPH_MODE)
      call MphaseUpdateAuxVars(realization)
    case(MIS_MODE)
      call MiscibleUpdateAuxVars(realization)
    case(THC_MODE)
      call THCUpdateAuxVars(realization)
    case(THMC_MODE)
      call THMCUpdateAuxVars(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateAuxVars(realization)
    case(G_MODE)
      call GeneralUpdateAuxVars(realization)
  end select    

end subroutine StepperUpdateFlowAuxVars

! ************************************************************************** !
!
! StepperUpdateTranAuxVars: Updates the flow auxilliary variables
! author: Glenn Hammond
! date: 10/11/08 
!
! ************************************************************************** !
subroutine StepperUpdateTranAuxVars(realization)
  
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Realization_module

  implicit none

  type(realization_type) :: realization

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)

end subroutine StepperUpdateTranAuxVars

! ************************************************************************** !
!
! StepperSandbox: Sandbox for temporary miscellaneous operations
! author: Glenn Hammond
! date: 06/27/11
!
! ************************************************************************** !
subroutine StepperSandbox(realization)
  
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Realization_module
  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Option_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use String_module
  use Secondary_Continuum_module

  implicit none

  type(realization_type) :: realization
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option

  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal, pointer :: volume_p(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend
  PetscInt :: species_offset
  PetscInt :: num_iterations
  PetscErrorCode :: ierr
  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim

  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction
  option => realization%option

  rt_aux_vars => patch%Aux%RT%aux_vars
  global_aux_vars => patch%Aux%Global%aux_vars
  rt_sec_transport_vars => patch%Aux%RT%sec_transport_vars

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  call GridVecGetArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecGetArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecGetArrayF90(grid,field%volume,volume_p,ierr)
  
  vol_frac_prim = 1.d0

  do species_offset = 1, reaction%naqcomp
    if (StringCompare(reaction%primary_species_names(species_offset), &
                      'UO2++',FIVE_INTEGER)) then
      ! decrement
      exit
    endif
  enddo
  species_offset = species_offset - 1

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    
    iend = local_id*reaction%naqcomp
    istart = iend-reaction%naqcomp+1
    
    if (option%use_mc) then
      vol_frac_prim = rt_sec_transport_vars(ghosted_id)%epsilon
    endif
    
    tran_xx_p(istart:iend) = rt_aux_vars(ghosted_id)%total(:,ONE_INTEGER)
    ! scale uo2++ total concentration by 0.25
    tran_xx_p(istart + species_offset) = 0.25d0 * &
                                         tran_xx_p(istart + species_offset)

    call RReact(rt_aux_vars(ghosted_id),global_aux_vars(ghosted_id), &
                tran_xx_p(istart:iend),volume_p(local_id), &
                porosity_loc_p(ghosted_id), &
                num_iterations,reaction,option,vol_frac_prim)
    tran_xx_p(istart:iend) = rt_aux_vars(ghosted_id)%pri_molal
  enddo

  call GridVecRestoreArrayF90(grid,field%tran_xx,tran_xx_p,ierr)
  call GridVecRestoreArrayF90(grid,field%porosity_loc, porosity_loc_p, ierr)  
  call GridVecRestoreArrayF90(grid,field%volume,volume_p,ierr)
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

end subroutine StepperSandbox

! ************************************************************************** !
!
! StepperCheckpoint: Calls appropriate routines to write a checkpoint file
! author: Glenn Hammond
! date: 03/07/08 
!
! ************************************************************************** !
subroutine StepperCheckpoint(realization,flow_stepper,tran_stepper,id)

  use Realization_module
  use Checkpoint_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  PetscInt :: num_const_timesteps, num_newton_iterations  
  PetscInt :: id

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_cumulative_newton_iterations, &
              flow_cumulative_time_step_cuts, flow_cumulative_linear_iterations, &
              flow_num_const_time_steps, flow_num_newton_iterations
  PetscInt :: tran_steps, tran_cumulative_newton_iterations, &
              tran_cumulative_time_step_cuts, tran_cumulative_linear_iterations, &
              tran_num_const_time_steps, tran_num_newton_iterations
  PetscReal :: flow_cumulative_solver_time, flow_prev_dt
  PetscReal :: tran_cumulative_solver_time,tran_prev_dt
  
  option => realization%option

  if (associated(flow_stepper)) then
    flow_steps = flow_stepper%steps
    flow_cumulative_newton_iterations = flow_stepper%cumulative_newton_iterations
    flow_cumulative_time_step_cuts = flow_stepper%cumulative_time_step_cuts
    flow_cumulative_linear_iterations = flow_stepper%cumulative_linear_iterations
    flow_num_const_time_steps = flow_stepper%num_constant_time_steps
    flow_num_newton_iterations = flow_stepper%num_newton_iterations
    flow_cumulative_solver_time = flow_stepper%cumulative_solver_time
    flow_prev_dt = flow_stepper%prev_dt
  endif
  if (associated(tran_stepper)) then
    tran_steps = tran_stepper%steps
    tran_cumulative_newton_iterations = tran_stepper%cumulative_newton_iterations
    tran_cumulative_time_step_cuts = tran_stepper%cumulative_time_step_cuts
    tran_cumulative_linear_iterations = tran_stepper%cumulative_linear_iterations
    tran_num_const_time_steps = tran_stepper%num_constant_time_steps
    tran_num_newton_iterations = tran_stepper%num_newton_iterations
    tran_cumulative_solver_time = tran_stepper%cumulative_solver_time
    tran_prev_dt = tran_stepper%prev_dt
  endif
  
  call Checkpoint(realization, &
                  flow_steps,flow_cumulative_newton_iterations, &
                  flow_cumulative_time_step_cuts,flow_cumulative_linear_iterations, &
                  flow_num_const_time_steps,flow_num_newton_iterations, &
                  flow_cumulative_solver_time,flow_prev_dt, &
                  tran_steps,tran_cumulative_newton_iterations, &
                  tran_cumulative_time_step_cuts,tran_cumulative_linear_iterations, &
                  tran_num_const_time_steps,tran_num_newton_iterations, &
                  tran_cumulative_solver_time,tran_prev_dt, &
                  id)
                      
end subroutine StepperCheckpoint

! ************************************************************************** !
!
! StepperRestart: Calls appropriate routines to read checkpoint file and
!                 restart
! author: Glenn Hammond
! date: 03/07/08 
!
! ************************************************************************** !
subroutine StepperRestart(realization,flow_stepper,tran_stepper, &
                          flow_read,transport_read,activity_coefs_read)

  use Realization_module
  use Checkpoint_module
  use Option_module

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
    
end subroutine StepperRestart

! ************************************************************************** !
!
! TimestepperGetTranWeight: Sets the weights at t0 or t1 for transport
! author: Glenn Hammond
! date: 01/17/11
!
! ************************************************************************** !
subroutine TimestepperSetTranWeights(option,flow_t0,flow_t1)

  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscReal :: flow_t0
  PetscReal :: flow_t1

  ! option%tran_time is the time at t0
  option%tran_weight_t0 = (option%tran_time-flow_t0)/ &
                          (flow_t1-flow_t0)
  option%tran_weight_t1 = (option%tran_time+option%tran_dt-flow_t0)/ &
                          (flow_t1-flow_t0)

end subroutine TimestepperSetTranWeights

! ************************************************************************** !
!
! TimestepperCheckCFLLimit: Checks CFL limit specified by the user
! author: Glenn Hammond
! date: 01/17/11
!
! ************************************************************************** !
subroutine TimestepperCheckCFLLimit(stepper,realization)

  use Realization_module
  
  implicit none

  type(stepper_type) :: stepper
  type(realization_type) :: realization
  
  PetscReal :: dt_cfl_1
  
  call RealizationCalculateCFL1Timestep(realization,dt_cfl_1)
  stepper%cfl_limiter_ts = dt_cfl_1*stepper%cfl_limiter
 
end subroutine TimestepperCheckCFLLimit

! ************************************************************************** !
!
! TimestepperEnforceCFLLimit: Enforces a CFL limit specified by the user
! author: Glenn Hammond
! date: 01/17/11
!
! ************************************************************************** !
subroutine TimestepperEnforceCFLLimit(stepper,option,output_option)

  use Option_module
  use Output_Aux_module

  implicit none

  type(stepper_type) :: stepper
  type(option_type) :: option
  type(output_option_type) :: output_option
  
  if (stepper%cfl_limiter_ts < option%tran_dt) then
    option%tran_dt = stepper%cfl_limiter_ts
    if (OptionPrintToScreen(option)) then
      write(*,'(" CFL Limiting: ",1pe12.4," [",a1,"]")') &
            stepper%cfl_limiter_ts/output_option%tconv,output_option%tunit
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out,'(/," CFL Limiting: ",1pe12.4," [",a1,"]",/)') &
            stepper%cfl_limiter_ts/output_option%tconv,output_option%tunit
    endif        
  endif    

end subroutine TimestepperEnforceCFLLimit

! ************************************************************************** !
!
! TimestepperPrintInfo: Prints information about time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperPrintInfo(stepper,fid,header,option)

  use Option_module
  
  implicit none
  
  type(stepper_type) :: stepper
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(string,*) stepper%max_time_step
    write(*,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) stepper%constant_time_step_threshold
    write(*,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) stepper%max_time_step_cuts
    write(*,'("max cuts:",x,a)') trim(adjustl(string))
  endif
  if (OptionPrintToFile(option)) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(string,*) stepper%max_time_step
    write(fid,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) stepper%constant_time_step_threshold
    write(fid,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) stepper%max_time_step_cuts
    write(fid,'("max cuts:",x,a)') trim(adjustl(string))
  endif    

end subroutine TimestepperPrintInfo

! ************************************************************************** !
!
! TimestepperDestroy: Deallocates a time stepper
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine TimestepperDestroy(stepper)

  implicit none
  
  type(stepper_type), pointer :: stepper
  
  if (.not.associated(stepper)) return
    
  call SolverDestroy(stepper%solver)
  call ConvergenceContextDestroy(stepper%convergence_context)
  
  if (associated(stepper%tfac)) deallocate(stepper%tfac)
  nullify(stepper%tfac)
  deallocate(stepper)
  nullify(stepper)
  
end subroutine TimestepperDestroy

end module Timestepper_module
