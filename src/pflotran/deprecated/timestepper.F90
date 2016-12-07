module Timestepper_module
 
  use Solver_module
  use Waypoint_module 
  use Convergence_module 
  use Material_module
  use Material_Aux_class
  use Variables_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: TIMESTEPPER_INIT_PROCEED = 0
  PetscInt, parameter, public :: TIMESTEPPER_INIT_DONE = 1
  PetscInt, parameter, public :: TIMESTEPPER_INIT_FAIL = 2
 
  type, public :: timestepper_type
  
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
    
    PetscReal :: dt_init
    PetscReal :: dt_min
    PetscReal :: dt_max
    PetscReal :: prev_dt
    PetscReal :: cfl_limiter
    PetscReal :: cfl_limiter_ts
    
    PetscBool :: init_to_steady_state
    PetscBool :: run_as_steady_state
    PetscReal :: steady_state_rel_tol
    
    PetscBool :: time_step_cut_flag  ! flag toggled if timestep is cut

    PetscLogDouble :: start_time    
    PetscInt :: start_time_step ! the first time step of a given run
    PetscReal :: time_step_tolerance ! scalar used in determining time step size
    PetscReal :: target_time    ! time at end of "synchronized" time step 

    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac
            
    type(solver_type), pointer :: solver
    
    type(waypoint_type), pointer :: cur_waypoint

    type(convergence_context_type), pointer :: convergence_context

  end type timestepper_type
  
  public :: TimestepperCreate, TimestepperDestroy, &
#if 0
            TimestepperExecuteRun, &
            TimestepperInitializeRun, &
            TimestepperFinalizeRun, &
#endif            
#if 0
            FlowStepperStepToSteadyState, &
            StepperCheckpoint, &
            StepperJumpStart, &
            StepperRunSteadyState, &
            StepperSetTargetTimes, &
            StepperStepFlowDT, &
            StepperStepTransportDT_GI, &
            StepperStepTransportDT_OS, &
            StepperUpdateDT, &
            StepperUpdateDTMax, &
            StepperUpdateSolution, &
            StepperUpdateTransportSolution, &
            TimestepperCheckCFLLimit, &
            TimestepperEnforceCFLLimit, &
            TimestepperRestart, &
#endif
            TimestepperRead, TimestepperPrintInfo, TimestepperReset
        

contains

! ************************************************************************** !

function TimestepperCreate()
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(timestepper_type), pointer :: TimestepperCreate
  
  type(timestepper_type), pointer :: stepper
  
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

  stepper%start_time = 0.d0  
  stepper%start_time_step = 0
  stepper%time_step_tolerance = 0.1d0
  stepper%target_time = 0.d0
  
  stepper%dt_init = 1.d0
  stepper%dt_min = 1.d-20   ! Ten zeptoseconds.
  stepper%dt_max = 3.1536d6 ! One-tenth of a year.  
  stepper%prev_dt = 0.d0
  stepper%cfl_limiter = UNINITIALIZED_DOUBLE
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
  
  TimestepperCreate => stepper
  
end function TimestepperCreate 

! ************************************************************************** !

subroutine TimestepperReset(stepper,dt_init)
  ! 
  ! Resets time stepper back to initial settings
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/27/11
  ! 

  implicit none

  type(timestepper_type) :: stepper
  PetscReal :: dt_init

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

  stepper%dt_init = dt_init
  stepper%prev_dt = 0.d0
  stepper%cfl_limiter = UNINITIALIZED_DOUBLE
  stepper%cfl_limiter_ts = 1.d20

  stepper%time_step_cut_flag = PETSC_FALSE

end subroutine TimestepperReset

! ************************************************************************** !

subroutine TimestepperRead(stepper,input,option)
  ! 
  ! Reads parameters associated with time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none

  type(timestepper_type) :: stepper
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

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

#if 0

subroutine StepperUpdateDT(flow_timestepper,tran_timestepper,option)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 

  use Option_module
  
  implicit none

  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
  type(option_type) :: option
  
  PetscReal :: time, dt
  PetscReal :: fac,dtt,up,utmp,uc,ut,uus,dt_tfac,dt_p
  PetscBool :: update_time_step
  PetscInt :: ifac

  ! FLOW
  update_time_step = PETSC_TRUE
  if (associated(flow_timestepper)) then
    if (flow_timestepper%run_as_steady_state) then
      update_time_step = PETSC_FALSE
    else
      ! if the time step was cut, set number of constant time steps to 1
      if (flow_timestepper%time_step_cut_flag) then
        flow_timestepper%time_step_cut_flag = PETSC_FALSE
        flow_timestepper%num_constant_time_steps = 1
      ! otherwise, only increment if teh constant time step counter was
      ! initialized to 1
      else if (flow_timestepper%num_constant_time_steps > 0) then
        flow_timestepper%num_constant_time_steps = &
          flow_timestepper%num_constant_time_steps + 1
      endif
    
      ! num_constant_time_steps = 0: normal time stepping with growing steps
      ! num_constant_time_steps > 0: restriction of constant time steps until
      !                              constant_time_step_threshold is met
      if (flow_timestepper%num_constant_time_steps > &
          flow_timestepper%constant_time_step_threshold) then
        flow_timestepper%num_constant_time_steps = 0
      else if (flow_timestepper%num_constant_time_steps > 0) then
        ! do not increase time step size
        update_time_step = PETSC_FALSE
      endif
    endif
    
    if (update_time_step) then
      
      time = option%flow_time
      dt = option%flow_dt

      if (flow_timestepper%iaccel == 0) return

      select case(option%iflowmode)
        case(FLASH2_MODE)   
          fac = 0.5d0
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
        case(TH_MODE)
          fac = 0.5d0
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
          if (flow_timestepper%iaccel > 0) then
            fac = 0.5d0
            if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
              fac = 0.33d0
              ut = 0.d0
            else
              up = option%dpmxe/(option%dpmax+0.1)
              ut = up
            endif
            dtt = fac * dt * (1.d0 + ut)
          else
            ifac = max(min(flow_timestepper%num_newton_iterations, &
                           flow_timestepper%ntfac),1)
            dt_tfac = flow_timestepper%tfac(ifac) * dt

            fac = 0.5d0
            up = option%dpmxe/(option%dpmax+0.1)
            dt_p = fac * dt * (1.d0 + up)

            dtt = min(dt_tfac,dt_p)
          endif
        case(G_MODE)   
          fac = 0.5d0
          if (flow_timestepper%num_newton_iterations >= flow_timestepper%iaccel) then
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
          if (flow_timestepper%num_newton_iterations <= flow_timestepper%iaccel .and. &
              flow_timestepper%num_newton_iterations <= size(flow_timestepper%tfac)) then
            if (flow_timestepper%num_newton_iterations == 0) then
              dtt = flow_timestepper%tfac(1) * dt
            else
              dtt = flow_timestepper%tfac(flow_timestepper%num_newton_iterations) * dt
            endif
          endif
          
      end select

      if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
      if (dtt > flow_timestepper%dt_max) dtt = flow_timestepper%dt_max
      ! geh: There used to be code here that cut the time step if it is too
      !      large relative to the simulation time.  This has been removed.
      dt = dtt

      option%flow_dt = dt
    
    endif
        
  endif
  
  ! TRANSPORT
  update_time_step = PETSC_TRUE
  if (associated(tran_timestepper)) then
    ! if the time step was cut, set number of constant time steps to 1
    if (tran_timestepper%time_step_cut_flag) then
      tran_timestepper%time_step_cut_flag = PETSC_FALSE
      tran_timestepper%num_constant_time_steps = 1
    ! otherwise, only increment if teh constant time step counter was
    ! initialized to 1
    else if (tran_timestepper%num_constant_time_steps > 0) then
      tran_timestepper%num_constant_time_steps = &
        tran_timestepper%num_constant_time_steps + 1
    endif
    
    ! num_constant_time_steps = 0: normal time stepping with growing steps
    ! num_constant_time_steps > 0: restriction of constant time steps until
    !                              constant_time_step_threshold is met
    if (tran_timestepper%num_constant_time_steps > &
        tran_timestepper%constant_time_step_threshold) then
      tran_timestepper%num_constant_time_steps = 0
    else if (tran_timestepper%num_constant_time_steps > 0) then
      ! do not increase time step size
      update_time_step = PETSC_FALSE
    endif
    
    if (update_time_step) then
    
      time = option%tran_time
      dt = option%tran_dt

      if (tran_timestepper%iaccel == 0) return

      dtt = dt
      if (tran_timestepper%num_newton_iterations <= tran_timestepper%iaccel) then
        if (tran_timestepper%num_newton_iterations <= size(tran_timestepper%tfac)) then
          dtt = tran_timestepper%tfac(tran_timestepper%num_newton_iterations) * dt
        else
          dtt = 0.5d0 * dt
        endif
      else
!       dtt = 2.d0 * dt
        dtt = 0.5d0 * dt
      endif

      if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
      if (dtt > tran_timestepper%dt_max) dtt = tran_timestepper%dt_max
      ! geh: see comment above under flow stepper
      dt = dtt

      option%tran_dt = dt
        
    endif
    
  endif

  ! ensure that transport time step is not larger than flow time step
  if (associated(flow_timestepper) .and. associated(tran_timestepper)) then
    if (.not.flow_timestepper%run_as_steady_state) then
      option%tran_dt = min(option%tran_dt,option%flow_dt)
    endif
  endif

end subroutine StepperUpdateDT

! ************************************************************************** !

subroutine StepperUpdateDTMax(flow_timestepper,tran_timestepper,option)
  ! 
  ! Updates maximum time step specified by the current
  ! waypoint after the completion of a time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/23/12
  ! 

  use Option_module
  
  implicit none

  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
  type(option_type) :: option
  
  PetscReal :: dt_max
  PetscBool :: flag
  type(waypoint_type), pointer :: cur_waypoint

  ! this just works through the logic of whether the flow stepper should
  ! be used to set dt_max.  If flow is not present, or is to
  ! be run as steady state, use the transport parameters
  flag = PETSC_TRUE
  if (associated(flow_timestepper)) then
    if (flow_timestepper%run_as_steady_state) then
      flag = PETSC_FALSE
    endif
  else
    flag = PETSC_FALSE
  endif

  if (flag) then ! flow stepper will govern the target time
    cur_waypoint => flow_timestepper%cur_waypoint
  else
    cur_waypoint => tran_timestepper%cur_waypoint
  endif
  
  ! update maximum time step size to current waypoint value
  if (associated(cur_waypoint)) then
    dt_max = cur_waypoint%dt_max
    
    ! target time will always be dictated by the flow solver, if present
    if (associated(flow_timestepper)) then
      flow_timestepper%dt_max = dt_max
    endif
    if (associated(tran_timestepper)) then
      tran_timestepper%dt_max = dt_max
    endif

  endif
  
end subroutine StepperUpdateDTMax

! ************************************************************************** !

subroutine StepperSetTargetTimes(flow_timestepper,tran_timestepper,option, &
                                 snapshot_plot_flag,observation_plot_flag, &
                                 observation_plot_flag,checkpoint_flag)
  !
  ! Sets target time for flow and transport solvers
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! Note: Modified by Jenn Frederick, 2/23/2016
  !

  use Option_module
  
  implicit none

  type(timestepper_type), pointer :: flow_timestepper, tran_timestepper
  type(option_type) :: option
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag                                                               
  
  PetscReal :: target_time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: cumulative_time_steps
  PetscInt :: max_time_step
  PetscBool :: flag
  PetscReal :: tolerance
  type(waypoint_type), pointer :: cur_waypoint

  checkpoint_flag = PETSC_FALSE
  ! target time will always be dictated by the flow solver, if present
  ! this does not mean that the transport solver cannot take smaller steps
  
  ! this just works through the logic of whether the flow stepper should
  ! be used to set the target time.  If flow is not present, or is to
  ! be run as steady state, use the transport parameters
  flag = PETSC_TRUE
  if (associated(flow_timestepper)) then
    if (flow_timestepper%run_as_steady_state) then
      flag = PETSC_FALSE
    endif
  else
    flag = PETSC_FALSE
  endif

  ! this flag allows one to take a shorter or slightly larger step than normal
  ! in order to synchronize with the waypoint time
  if (option%match_waypoint) then
    ! if the maximum time step size decreased in the past step, need to set
    ! the time step size to the minimum of the stepper%prev_dt and 
    ! stepper%dt_max
    if (associated(flow_timestepper)) option%flow_dt = &
      min(flow_timestepper%prev_dt,flow_timestepper%dt_max)
    if (associated(tran_timestepper)) option%tran_dt = &
      min(tran_timestepper%prev_dt,tran_timestepper%dt_max)
    option%match_waypoint = PETSC_FALSE
  endif

  if (flag) then ! flow stepper will govern the target time
!    time = option%flow_time + option%flow_dt
    dt = option%flow_dt
    ! dt_max must be set from current waypoint and not updated below
    cur_waypoint => flow_timestepper%cur_waypoint
    dt_max = cur_waypoint%dt_max
    cumulative_time_steps = flow_timestepper%steps
    max_time_step = flow_timestepper%max_time_step
    tolerance = flow_timestepper%time_step_tolerance
    target_time = flow_timestepper%target_time + option%flow_dt
  else
!    time = option%tran_time + option%tran_dt
    dt = option%tran_dt
    cur_waypoint => tran_timestepper%cur_waypoint
    ! dt_max must be set from current waypoint and not updated below
    dt_max = cur_waypoint%dt_max
    cumulative_time_steps = tran_timestepper%steps
    max_time_step = tran_timestepper%max_time_step
    tolerance = tran_timestepper%time_step_tolerance
    target_time = tran_timestepper%target_time + option%tran_dt
  endif
  
  ! For the case where the second waypoint is a printout after the first 
  ! time step we must increment the waypoint beyond the first (time=0.) 
  ! waypoint.  Otherwise the second time step will be zero. - geh
  if (cur_waypoint%time < 1.d-40) then
    cur_waypoint => cur_waypoint%next
  endif
  
  if (associated(flow_timestepper)) flow_timestepper%prev_dt = option%flow_dt
  if (associated(tran_timestepper)) tran_timestepper%prev_dt = option%tran_dt

! If a waypoint calls for a plot or change in src/sinks, adjust time step
! to match waypoint.
  do ! we cycle just in case the next waypoint is beyond the target_time
    if (target_time + tolerance*dt >= cur_waypoint%time .and. &
        (cur_waypoint%sync .or. &
         cur_waypoint%update_conditions .or. &
         cur_waypoint%print_snap_output .or. &
         cur_waypoint%print_checkpoint .or. &
         cur_waypoint%print_obs_output .or. &
         cur_waypoint%print_msbl_output .or. &
         cur_waypoint%final)) then
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on waypoint time
      dt = cur_waypoint%time - target_time
      if (dt > dt_max .and. dabs(dt-dt_max) > 1.d0) then 
        ! 1 sec tolerance to avoid cancellation
        dt = dt_max                    ! error from waypoint%time - time
        target_time = target_time + dt
      else
        target_time = cur_waypoint%time
        if (cur_waypoint%print_snap_output) snapshot_plot_flag = PETSC_TRUE
        if (cur_waypoint%print_obs_output) observation_plot_flag = PETSC_TRUE
        if (cur_waypoint%print_msbl_output) massbal_plot_flag = PETSC_TRUE
        if (cur_waypoint%print_checkpoint) checkpoint_flag = PETSC_TRUE
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
  if (associated(flow_timestepper)) then
    option%flow_dt = dt
    flow_timestepper%target_time = target_time
    !geh: dt_max now updated in StepperUpdateDTMax() at the end of 
    !     a time step to avoid premature update if cuts or sub-stepping
    !     for transport occur during the time step.
    !flow_timestepper%dt_max = dt_max
    flow_timestepper%cur_waypoint => cur_waypoint
  endif
  if (associated(tran_timestepper)) then
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
    tran_timestepper%target_time = target_time
    !geh: see note on flow_timestepper%dt_max above.
    !tran_timestepper%dt_max = dt_max
    tran_timestepper%cur_waypoint => cur_waypoint
  endif
  
end subroutine StepperSetTargetTimes

! ************************************************************************** !

subroutine StepperStepFlowDT(realization,stepper,failure)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08, 03/11/13
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use Flash2_module, only : Flash2MaxChange, Flash2InitializeTimestep, &
                           Flash2TimeCut, Flash2UpdateReason
  use Mphase_module, only : MphaseMaxChange, MphaseInitializeTimestep, &
                           MphaseTimeCut, MPhaseUpdateReason
  use Immis_module, only : ImmisMaxChange, ImmisInitializeTimestep, &
                           ImmisTimeCut, ImmisUpdateReason
  use Miscible_module, only : MiscibleMaxChange, MiscibleInitializeTimestep, &
                           MiscibleTimeCut
  use Richards_module, only : RichardsMaxChange, RichardsInitializeTimestep, &
                             RichardsTimeCut, RichardsResidual
  use TH_module, only : THMaxChange, THInitializeTimestep, THTimeCut

  use General_module, only : GeneralInitializeTimestep, GeneralTimeCut
  use Global_module

  use Output_module, only : Output
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  PetscInt :: icut ! Tracks the number of time step reductions applied
  SNESConvergedReason :: snes_reason 
  PetscInt :: update_reason
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscReal :: fnorm, scaled_fnorm, inorm
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscLogDouble :: log_start_time, log_end_time

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(discretization_type), pointer :: discretization 
  type(solver_type), pointer :: solver
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

  icut = 0
  sum_newton_iterations = 0
  sum_linear_iterations = 0

  ! Perform some global-to-local scatters to update the ghosted vectors.
  ! We have to do this so that the routines for calculating the residual
  ! and the Jacobian will have the ghost points they need.
  ! Note that we don't do the global-to-local scatter for the pressure 
  ! vector, as that needs to be done within the residual calculation routine
  ! because that routine may get called several times during one Newton step
  ! if a method such as line search is being used.
  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                  field%icap_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                  field%ithrm_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif
    
  if (option%print_screen_flag) then
    write(*,'(/,2("=")," FLOW ",72("="))')
  endif

  if (option%ntrandof > 0) then ! store initial saturations for transport
    call GlobalUpdateAuxVars(realization,TIME_T, &
                             stepper%target_time-option%flow_dt)
  endif
    
  select case(option%iflowmode)
    case(TH_MODE)
      call THInitializeTimestep(realization)
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
      
    call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

    select case(option%iflowmode)
      case(MPH_MODE,TH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
        call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                       ierr);CHKERRQ(ierr)
      case(RICHARDS_MODE)
        call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                       ierr);CHKERRQ(ierr)
    end select
    call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
    stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
                                      (log_end_time - log_start_time)

! do we really need all this? - geh 
    call SNESGetIterationNumber(solver%snes,num_newton_iterations,  &
                                ierr);CHKERRQ(ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,  &
                                      ierr);CHKERRQ(ierr)
    call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
    update_reason = 1
      
    if (snes_reason >= 0) then
      select case(option%iflowmode)
        case(IMS_MODE)
          call ImmisUpdateReason(update_reason,realization)
        case(MPH_MODE)
        case(FLASH2_MODE)
        case(TH_MODE)
          update_reason=1
        case (MIS_MODE)
          update_reason=1
        case(RICHARDS_MODE,G_MODE)
          update_reason=1
        end select   
    endif
   
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
        snapshot_plot_flag = PETSC_TRUE
        observation_plot_flag = PETSC_FALSE
        massbal_plot_flag = PETSC_FALSE
        call Output(realization,snapshot_plot_flag,observation_plot_flag, &
                    massbal_plot_flag)
        failure = PETSC_TRUE
        return
      endif
 
      stepper%target_time = stepper%target_time - option%flow_dt

      option%flow_dt = 0.5d0 * option%flow_dt  
      
      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,stepper%cumulative_time_step_cuts, &
            option%flow_time/realization%output_option%tconv, &
            option%flow_dt/realization%output_option%tconv

      stepper%target_time = stepper%target_time + option%flow_dt

      select case(option%iflowmode)
        case(TH_MODE)
          call THTimeCut(realization)
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
      call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr);CHKERRQ(ierr)

    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo
  
  stepper%steps = stepper%steps + 1      
  stepper%cumulative_newton_iterations = &
    stepper%cumulative_newton_iterations + sum_newton_iterations
  stepper%cumulative_linear_iterations = &
    stepper%cumulative_linear_iterations + sum_linear_iterations
  stepper%cumulative_time_step_cuts = &
    stepper%cumulative_time_step_cuts + icut

  stepper%num_newton_iterations = num_newton_iterations
  stepper%num_linear_iterations = num_linear_iterations

  if (option%ntrandof > 0) then ! store final saturations, etc. for transport
    call GlobalUpdateAuxVars(realization,TIME_TpDT,stepper%target_time)
  endif
    
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
    write(*, '(/," FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
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
    write(option%fid_out, '(" FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a,"]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts
  endif
  
  select case(option%iflowmode)
    case(TH_MODE)
      call THMaxChange(realization)
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax
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
!          call GeneralMaxChange(realization)
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

  if (option%print_screen_flag) print *, ""
  
end subroutine StepperStepFlowDT

! ************************************************************************** !

subroutine FlowStepperStepToSteadyState(realization,stepper,failure)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/12/13
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Global_module
  use Output_module, only : Output
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  PetscReal :: inorm, prev_norm, dif_norm, rel_norm
  PetscReal :: tempreal

  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(discretization_type), pointer :: discretization 
  type(solver_type), pointer :: solver
  
  Vec :: update_vec

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

  prev_norm = 1.d20

  if (option%print_screen_flag) then
    write(*,'(/,2("=")," Initialize FLOW to Steady-State (Begin)",38("="))')
  endif
  if (option%print_file_flag) then
    write(option%fid_out,'(/,2("="), &
      &" Initialize FLOW to Steady-State (Begin)",38("="))')
  endif

  do ! this loop is for steady state initial condition
   
    ! flow stepper increments time step. Need to decrement prior to step so that
    ! printout reflects time step 0
    stepper%steps = stepper%steps - 1      
    call StepperStepFlowDT(realization,stepper,failure)
  
    call SNESGetSolutionUpdate(solver%snes,update_vec,ierr);CHKERRQ(ierr)
    call VecStrideNorm(update_vec,ZERO_INTEGER,NORM_INFINITY,inorm, &
                       ierr);CHKERRQ(ierr)
    dif_norm = inorm-prev_norm
    rel_norm = dif_norm/prev_norm
    if (stepper%cumulative_newton_iterations > 20 .and. &
!geh: These norms don't seem to be the best, sticking with inorm
!           dabs(dif_norm) < stepper%steady_state_rel_tol) .or. &
!          dabs(rel_norm) < stepper%steady_state_rel_tol) then
        inorm < stepper%steady_state_rel_tol) then
      if (option%print_file_flag) then
        write(option%fid_out,'(/,80("="))')
        write(option%fid_out,*) 'Steady state solve converged after ', &
          stepper%cumulative_newton_iterations, ' iterations'
      endif
      exit
    endif
      
    prev_norm = inorm
      
    if (stepper%cumulative_newton_iterations > 1000) then
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
    if (stepper%num_newton_iterations < 4) then
      option%flow_dt = option%flow_dt*2.d0
    else if (stepper%num_newton_iterations < 5) then
      option%flow_dt = option%flow_dt*1.5d0
    else if (stepper%num_newton_iterations < 6) then
      option%flow_dt = option%flow_dt*1.d0
    else if (stepper%num_newton_iterations < 7) then
      option%flow_dt = option%flow_dt*0.9d0
    else if (stepper%num_newton_iterations < 8) then
      option%flow_dt = option%flow_dt*0.8d0
    else if (stepper%num_newton_iterations < 9) then
      option%flow_dt = option%flow_dt*0.7d0
    else
      option%flow_dt = option%flow_dt*0.6d0
    endif

    if (option%print_file_flag) then
      write(option%fid_out,*) 'Dt: ', option%flow_dt
      write(option%fid_out,*) 'Inf Norm: ', inorm
      write(option%fid_out,*) 'Relative Norm: ', rel_norm
      write(option%fid_out,*) 'Linear its/Newton It: ', &
        float(stepper%num_linear_iterations) / stepper%num_newton_iterations
    endif
    if (option%print_screen_flag) then
      write(*,*) 'Dt: ', option%flow_dt
      write(*,*) 'Inf Norm: ', inorm
      write(*,*) 'Relative Norm: ', rel_norm
      write(*,*) 'Linear its/Newton It: ', &
        float(stepper%num_linear_iterations) / stepper%num_newton_iterations
    endif

  enddo ! end of steady-state loop

  if (option%print_screen_flag) then
    write(*,'(/,2("=")," Initialize FLOW to Steady-State (End)",40("="),/, &
          & 80("="),/)')
  endif
  if (option%print_file_flag) then
    write(option%fid_out,'(/,2("="), &
      &" Initialize FLOW to Steady-State (End)",40("="),/,80("="),/)')
  endif

end subroutine FlowStepperStepToSteadyState

#if 0

! ************************************************************************** !

subroutine StepperStepFlowDT(realization,stepper,step_to_steady_state,failure)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Flash2_module, only : Flash2MaxChange, Flash2InitializeTimestep, &
                           Flash2TimeCut, Flash2UpdateReason
  use Mphase_module, only : MphaseMaxChange, MphaseInitializeTimestep, &
                           MphaseTimeCut, MPhaseUpdateReason
  use Immis_module, only : ImmisMaxChange, ImmisInitializeTimestep, &
                           ImmisTimeCut, ImmisUpdateReason
  use Miscible_module, only : MiscibleMaxChange, MiscibleInitializeTimestep, &
                           MiscibleTimeCut
  use Richards_module, only : RichardsMaxChange, RichardsInitializeTimestep, &
                             RichardsTimeCut, RichardsResidual
  use TH_module, only : THMaxChange, THInitializeTimestep, THTimeCut

  use General_module, only : GeneralMaxChange, GeneralInitializeTimestep, &
                             GeneralTimeCut
  use Global_module

  use Output_module, only : Output
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
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
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
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
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                    field%iphas_loc,ONEDOF)
  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif
    
    if (option%print_screen_flag) then
      if (step_to_steady_state) then
        write(*,'(/,2("=")," Initialize FLOW to Steady-State ",25("="))')
      else
        write(*,'(/,2("=")," FLOW ",52("="))')
      endif
    endif

    if (option%ntrandof > 0) then ! store initial saturations for transport
      call GlobalUpdateAuxVars(realization,TIME_T, &
                               stepper%target_time-option%flow_dt)
    endif
    
    select case(option%iflowmode)
      case(TH_MODE)
        call THInitializeTimestep(realization)
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
      
      call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

      select case(option%iflowmode)
        case(MPH_MODE,TH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
          call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                         ierr);CHKERRQ(ierr)
        case(RICHARDS_MODE)
          call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                         ierr);CHKERRQ(ierr)
      end select
      call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
                                       (log_end_time - log_start_time)

  ! do we really need all this? - geh 
      call SNESGetIterationNumber(solver%snes,num_newton_iterations,  &
                                  ierr);CHKERRQ(ierr)
      call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,  &
                                        ierr);CHKERRQ(ierr)
      call SNESGetConvergedReason(solver%snes, snes_reason,  &
                                  ierr);CHKERRQ(ierr)

      sum_newton_iterations = sum_newton_iterations + num_newton_iterations
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
      update_reason = 1
      
      if (snes_reason >= 0) then
        select case(option%iflowmode)
          case(IMS_MODE)
            call ImmisUpdateReason(update_reason,realization)
          case(MPH_MODE)
!           call MPhaseUpdateReason(update_reason,realization)
          case(FLASH2_MODE)
!           call Flash2UpdateReason(update_reason,realization)
          case(TH_MODE)
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

        if (icut > stepper%max_time_step_cuts .or. option%flow_dt<1.d-20) then
          if (option%print_screen_flag) then
            print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                    stepper%max_time_step_cuts, "t= ", &
                    stepper%target_time/realization%output_option%tconv, " dt= ", &
                    option%flow_dt/realization%output_option%tconv
            print *,"Stopping execution!"
          endif
          realization%output_option%plot_name = 'flow_cut_to_failure'
          snapshot_plot_flag = PETSC_TRUE
          observation_plot_flag = PETSC_FALSE
          massbal_plot_flag = PETSC_FALSE
          call Output(realization,snapshot_plot_flag,observation_plot_flag, &
                      massbal_plot_flag)
          failure = PETSC_TRUE
          return
        endif

 
       stepper%target_time = stepper%target_time - option%flow_dt

        option%flow_dt = 0.5d0 * option%flow_dt  
      
        if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
          &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
          &   1pe12.5)')  snes_reason,icut,stepper%cumulative_time_step_cuts, &
              option%flow_time/realization%output_option%tconv, &
              option%flow_dt/realization%output_option%tconv

        stepper%target_time = stepper%target_time + option%flow_dt

        select case(option%iflowmode)
          case(TH_MODE)
            call THTimeCut(realization)
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
        call VecCopy(field%iphas_old_loc, field%iphas_loc, ierr);CHKERRQ(ierr)

      else
        ! The Newton solver converged, so we can exit.
        exit
      endif
    enddo
    
    if (.not.step_to_steady_state) then
      exit
    else
    
      call SNESGetSolutionUpdate(solver%snes,update_vec,ierr);CHKERRQ(ierr)
      call VecStrideNorm(update_vec,ZERO_INTEGER,NORM_INFINITY,inorm, &
                         ierr);CHKERRQ(ierr)
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
    call GlobalUpdateAuxVars(realization,TIME_TpDT,stepper%target_time)
  endif
    
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
    write(*, '(/," FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
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
    write(option%fid_out, '(" FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      " [",a, "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      stepper%target_time/realization%output_option%tconv, &
      option%flow_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts
  endif
  
  select case(option%iflowmode)
    case(TH_MODE)
      call THMaxChange(realization)
      if (option%print_screen_flag) then
        write(*,'("  --> max chng: dpmx= ",1pe12.4, &
          & " dtmpmx= ",1pe12.4)') &
          option%dpmax,option%dtmpmax
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

  if (option%print_screen_flag) print *, ""
  
  ! option%flow_time is updated outside this subroutine

end subroutine StepperStepFlowDT
#endif

! ************************************************************************** !

subroutine StepperStepTransportDT_GI(realization,stepper, &
                                     steady_flow,flow_t0,flow_t1, &
                                     failure)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 
  
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Reactive_Transport_module
  use Output_module, only : Output
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module
  use Patch_module
  use Global_module  
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  
  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
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
  PetscBool :: snapshot_plot_flag  
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscReal, pointer :: r_p(:), xx_p(:), log_xx_p(:)
  PetscReal :: final_tran_time
  PetscLogDouble :: log_start_time, log_end_time

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field  
  type(solver_type), pointer :: solver
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif

  ! interpolate flow parameters/data
  ! this must remain here as these weighted values are used by both
  ! RTInitializeTimestep and RTTimeCut (which calls RTInitializeTimestep)
  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call TimestepperSetTranWeights(option,flow_t0, flow_t1)
    ! set densities and saturations to t
    call GlobalWeightAuxVars(realization,option%tran_weight_t0)
  endif

  call RTInitializeTimestep(realization)
  ! note: RTUpdateTransportCoefs() is called within RTInitializeTimestep()

  ! set densities and saturations to t+dt
  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call GlobalWeightAuxVars(realization,option%tran_weight_t1)
  endif

  call RTUpdateTransportCoefs(realization)
  
  sum_newton_iterations = 0
  sum_linear_iterations = 0
  icut = 0

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT ",67("="))')

  do
   
    final_tran_time = option%tran_time + option%tran_dt
   
    if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
      call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
!       The below is set within RTUpdateAuxVarsPatch() when PETSC_TRUE,PETSC_TRUE,* are passed
!       patch%aux%RT%auxvars_up_to_date = PETSC_TRUE 
    endif
    if (realization%reaction%use_log_formulation) then
      call VecCopy(field%tran_xx,field%tran_log_xx,ierr);CHKERRQ(ierr)
      call VecLog(field%tran_log_xx,ierr);CHKERRQ(ierr)

      call PetscTime(log_start_time, ierr);CHKERRQ(ierr)
      call SNESSolve(solver%snes, PETSC_NULL_VEC, field%tran_log_xx,  &
                     ierr);CHKERRQ(ierr)
      call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
        (log_end_time - log_start_time)          
        
      call VecCopy(field%tran_log_xx,field%tran_xx,ierr);CHKERRQ(ierr)
      call VecExp(field%tran_xx,ierr);CHKERRQ(ierr)
    else
      call PetscTime(log_start_time, ierr);CHKERRQ(ierr)
      call SNESSolve(solver%snes, PETSC_NULL_VEC, field%tran_xx,  &
                     ierr);CHKERRQ(ierr)
      call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
      stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
        (log_end_time - log_start_time)          
    endif

! do we really need all this? - geh 
    call SNESGetIterationNumber(solver%snes,num_newton_iterations,  &
                                ierr);CHKERRQ(ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,  &
                                      ierr);CHKERRQ(ierr)
    call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

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
        snapshot_plot_flag = PETSC_TRUE
        observation_plot_flag = PETSC_FALSE
        massbal_plot_flag = PETSC_FALSE
        call Output(realization,snapshot_plot_flag,observation_plot_flag, &
                    massbal_plot_flag)
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
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
  
    if (option%nflowdof > 0 .and. .not.steady_flow) then

    write(*, '(/," TRAN ",i6," Time= ",1pe12.5," Target= ",1pe12.5, &
      & " Dt= ",1pe12.5," [",a,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      flow_t1/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
      stepper%cumulative_newton_iterations,sum_linear_iterations, &
      stepper%cumulative_linear_iterations,icut, &
      stepper%cumulative_time_step_cuts

    else

    write(*, '(/," TRAN ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
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
      & " [",a,"]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i6,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      stepper%steps, &
      final_tran_time/realization%output_option%tconv, &
      option%tran_dt/realization%output_option%tconv, &
      trim(realization%output_option%tunit),snes_reason,sum_newton_iterations, &
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

subroutine StepperStepTransportDT_OS(realization,stepper, &
                                     steady_flow, flow_t0,flow_t1, &
                                     failure)
  ! 
  ! Steps forward one step in time (operator split)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Reactive_Transport_module, only : RTUpdateRHSCoefs, RTUpdateAuxVars, &
        RTCalculateRHS_t0, RTUpdateTransportCoefs, RTCalculateRHS_t1, &
        RTCalculateTransportMatrix, RTReact, RTMaxChange, RTExplicitAdvection
  use Output_module, only : Output
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Grid_module
  use Patch_module
  use Global_module  

  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
  PetscReal :: flow_t0, flow_t1
  PetscBool :: steady_flow
  PetscBool :: failure
  
  PetscErrorCode :: ierr
  KSPConvergedReason :: ksp_reason 
  PetscInt :: sum_linear_iterations, num_linear_iterations
  PetscInt :: idof
  PetscBool :: snapshot_plot_flag  
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
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

  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT (OS) ",43("="))')

  final_tran_time = option%tran_time + option%tran_dt
   
  sum_linear_iterations = 0
  ksp_reason = 0
  
  ! do we need this...don't think so - geh
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

  if (option%nflowdof > 0 .and. .not.steady_flow) then
    call TimestepperSetTranWeights(option,flow_t0,flow_t1)
    ! set densities and saturations to t
    call GlobalWeightAuxVars(realization,option%tran_weight_t0)
  endif

  ! update time derivative on RHS
  call RTUpdateRHSCoefs(realization)
  if (option%itranmode == EXPLICIT_ADVECTION) then
    ! update the totals (must be at t0)
                                   ! cells      bcs         act. coefs.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)
    ! set densities and saturations to t+dt
    if (option%nflowdof > 0 .and. .not.steady_flow) then
      call GlobalWeightAuxVars(realization,option%tran_weight_t1)
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
      call GlobalWeightAuxVars(realization,option%tran_weight_t1)
    endif

    ! update diffusion/dispersion coefficients
    call RTUpdateTransportCoefs(realization)
    ! RTCalculateRHS_t1() updates aux vars (cell and boundary) to k+1 and calculates 
    ! RHS fluxes and src/sinks
    call RTCalculateRHS_t1(realization)

    if (realization%debug%vecview_residual) then
      call PetscViewerASCIIOpen(realization%option%mycomm,'Trhs.out', &
                                viewer,ierr);CHKERRQ(ierr)
      call VecView(field%tran_rhs,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    ! RTCalculateTransportMatrix() calculates flux coefficients and the
    ! t^(k+1) coefficient in accumulation term
    call RTCalculateTransportMatrix(realization,solver%J) 
    call KSPSetOperators(solver%ksp,solver%J,solver%Jpre, &
                         SAME_NONZERO_PATTERN,ierr);CHKERRQ(ierr)

  !  call VecGetArrayF90(field%tran_xx,vec_ptr,ierr)
  !  call VecRestoreArrayF90(field%tran_xx,vec_ptr,ierr)

    ! loop over chemical component and transport
    do idof = 1, option%ntrandof

  ! for debugging
#if 0    
      call RealizationGetVariable(realization,field%work,TOTAL_MOLARITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  
      call RealizationGetVariable(realization,field%work,TOTAL_MOLALITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  
      call RealizationGetVariable(realization,field%work,PRIMARY_MOLALITY,idof)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
#endif
    
      call VecStrideGather(field%tran_rhs,idof-1,field%work,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
      option%rt_idof = idof
      call KSPSolve(solver%ksp,field%work,field%work,ierr);CHKERRQ(ierr)
      ! tran_xx will contain transported totals
      ! tran_xx_loc will still contain free-ion from previous solution
      call VecStrideScatter(field%work,idof-1,field%tran_xx,INSERT_VALUES, &
                            ierr);CHKERRQ(ierr)

  ! for debugging
#if 0
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
#endif      

      call KSPGetIterationNumber(solver%ksp,num_linear_iterations, &
                                 ierr);CHKERRQ(ierr)
      call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr);CHKERRQ(ierr)
      sum_linear_iterations = sum_linear_iterations + num_linear_iterations
    enddo
    sum_linear_iterations = int(dble(sum_linear_iterations) / option%ntrandof)
    stepper%cumulative_linear_iterations = &
      stepper%cumulative_linear_iterations + sum_linear_iterations
  endif ! if (EXPLICIT_ADVECTION)

  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'Txx.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  ! activity coefficients are updated within RReact!  DO NOT updated
  ! here as doing so will cause errors in the t0 portion of the
  ! accumulation term for equilibrium sorbed species
  call RTReact(realization)
  
  if (realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(realization%option%mycomm,'RTxx.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecView(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call PetscBarrier(solver%ksp,ierr);CHKERRQ(ierr)
  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)
  stepper%cumulative_solver_time = stepper%cumulative_solver_time + &
                                   (log_end_time - log_start_time)          

  stepper%num_newton_iterations = 1
  stepper%num_linear_iterations = sum_linear_iterations

  ! increment time step number
  stepper%steps = stepper%steps + 1     
  
  if (option%print_screen_flag) then
    write(*, '(" TRAN ",i6," Time= ",1pe12.5," Dt= ", &
          & 1pe12.5," [",a,"]"," ksp_conv_reason: ",i4,/," linear = ",i5, &
          & " [",i10,"]")') stepper%steps, &
!geh        option%tran_time/realization%output_option%tconv, &
          final_tran_time/realization%output_option%tconv, &
          option%tran_dt/realization%output_option%tconv, &
          trim(realization%output_option%tunit),ksp_reason, &
          sum_linear_iterations,stepper%cumulative_linear_iterations
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" TRAN ",i6," Time= ",1pe12.5," Dt= ", &
          & 1pe12.5," [",a,"]"," ksp_conv_reason = ",i4,/," linear = ",i5, &
          & " [",i10,"]")') stepper%steps, &
!geh        option%tran_time/realization%output_option%tconv, &
          final_tran_time/realization%output_option%tconv, &
          option%tran_dt/realization%output_option%tconv, &
          trim(realization%output_option%tunit),ksp_reason, &
          sum_linear_iterations,stepper%cumulative_linear_iterations
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

subroutine StepperRunSteadyState(realization,flow_timestepper,tran_timestepper)
  ! 
  ! Solves steady state solution for flow and transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/09
  ! 

  use Realization_class

  use Option_module
  use Output_Aux_module
  use Output_module, only : Output, OutputInit
  use Logging_module
  use Discretization_module

  class(realization_type) :: realization
  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper

  PetscBool :: observation_plot_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: failure
  PetscLogDouble :: start_time, end_time
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  PetscErrorCode :: ierr

  option => realization%option
  output_option => realization%output_option

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  failure = PETSC_FALSE

  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

  ! print initial condition output if not a restarted sim
  if (associated(flow_timestepper)) then
    call OutputInit(flow_timestepper%steps)
  else
    call OutputInit(tran_timestepper%steps)
  endif

  if (output_option%print_initial_obs) observation_plot_flag = PETSC_TRUE
  if (output_option%print_initial_snap) snapshot_plot_flag = PETSC_TRUE
  if (output_option%print_initial_massbal) massbal_plot_flag = PETSC_FALSE
  call Output(realization,snapshot_plot_flag,observation_plot_flag, &
              massbal_plot_flag)

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

  snapshot_plot_flag = PETSC_TRUE
    
  if (associated(flow_timestepper)) then
    call PetscTime(start_time, ierr);CHKERRQ(ierr)
    call PetscLogStagePush(logging%stage(FLOW_STAGE),ierr);CHKERRQ(ierr)
    call StepperSolveFlowSteadyState(realization,flow_timestepper,failure)
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
    call PetscTime(end_time, ierr);CHKERRQ(ierr)
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

  if (associated(tran_timestepper)) then
    call PetscTime(start_time, ierr);CHKERRQ(ierr)
    call PetscLogStagePush(logging%stage(TRAN_STAGE),ierr);CHKERRQ(ierr)
    call StepperSolveTranSteadyState(realization,tran_timestepper,failure)
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
    call PetscTime(end_time, ierr);CHKERRQ(ierr)
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

  if (output_option%print_initial_snap) then
    output_option%plot_number = 1
    call Output(realization,snapshot_plot_flag,observation_plot_flag, &
                massbal_plot_flag)
  endif
   
  if (option%checkpoint_flag) then
    call StepperCheckpoint(realization,flow_timestepper,tran_timestepper, &
                           NEG_ONE_INTEGER)  
  endif

  if (OptionPrintToScreen(option)) then
    if (option%nflowdof > 0) then
      write(*,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_timestepper%cumulative_newton_iterations,flow_timestepper%cumulative_linear_iterations
    endif
    if (option%ntrandof > 0) then
      write(*,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_timestepper%cumulative_newton_iterations,tran_timestepper%cumulative_linear_iterations
    endif            
  endif

  if (OptionPrintToFile(option)) then
    if (option%nflowdof > 0) then
      write(option%fid_out,'(/," FLOW newton = ",i8," linear = ",i10)') &
            flow_timestepper%cumulative_newton_iterations, &
            flow_timestepper%cumulative_linear_iterations
    endif
    if (option%ntrandof > 0) then
      write(option%fid_out,'(/," TRAN newton = ",i8," linear = ",i10)') &
            tran_timestepper%cumulative_newton_iterations, &
            tran_timestepper%cumulative_linear_iterations
    endif            
  endif

  call PetscLogStagePop(ierr);CHKERRQ(ierr)


end subroutine StepperRunSteadyState

! ************************************************************************** !

subroutine StepperSolveFlowSteadyState(realization,stepper,failure)
  ! 
  ! Solves the steady-state flow equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/10/09
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Global_module, only : GlobalUpdateAuxVars
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  use Richards_module, only : RichardsInitializeTimestep

  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
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


  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                  field%icap_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                  field%ithrm_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif
    
  if (option%print_screen_flag) write(*,'(/,2("=")," FLOW (STEADY STATE) ",37("="))')

  select case(option%iflowmode)
    case(RICHARDS_MODE)
      call RichardsInitializeTimestep(realization)
  end select

  select case(option%iflowmode)
    case(MPH_MODE,TH_MODE,IMS_MODE,MIS_MODE,FLASH2_MODE,G_MODE)
      call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                     ierr);CHKERRQ(ierr)
    case(RICHARDS_MODE)
      call SNESSolve(solver%snes, PETSC_NULL_VEC, field%flow_xx,  &
                     ierr);CHKERRQ(ierr)
  end select

  call SNESGetIterationNumber(solver%snes,num_newton_iterations,  &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,  &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

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
    call GlobalUpdateAuxVars(realization,TIME_T, &
                             stepper%target_time-option%flow_dt)
    call GlobalUpdateAuxVars(realization,TIME_TpDT,stepper%target_time)
  endif
    
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(field%flow_r,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
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

subroutine StepperSolveTranSteadyState(realization,stepper,failure)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 
  
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Realization_class
  use Discretization_module
  use Option_module
  use Solver_module
  use Field_module
  
  use Patch_module
  use Grid_module
    
  use Global_module, only : GlobalWeightAuxVars
  use Reactive_Transport_module, only : RTUpdateAuxVars  
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF

  implicit none

  class(realization_type) :: realization
  type(timestepper_type) :: stepper
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

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  solver => stepper%solver

! PetscReal, pointer :: xx_p(:), conc_p(:), press_p(:), temp_p(:)

  if (field%porosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%porosity_loc, &
                                    field%porosity_loc,ONEDOF)
  endif
  if (field%tortuosity_loc /= 0) then
    call DiscretizationLocalToLocal(discretization,field%tortuosity_loc, &
                                    field%tortuosity_loc,ONEDOF)
  endif
  if (.not.option%use_refactored_material_auxvars) then
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif

  call GlobalWeightAuxVars(realization,1.d0)
  num_newton_iterations = 0
  num_linear_iterations = 0

  if (option%print_screen_flag) write(*,'(/,2("=")" TRANSPORT (STEADY STATE) ",32("="))')

  if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
  endif

  if (realization%reaction%use_log_formulation) then
    call VecCopy(field%tran_xx,field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(field%tran_log_xx,ierr);CHKERRQ(ierr)
    call SNESSolve(solver%snes, PETSC_NULL_VEC, field%tran_log_xx,  &
                   ierr);CHKERRQ(ierr)
    call VecCopy(field%tran_log_xx,field%tran_xx,ierr);CHKERRQ(ierr)
    call VecExp(field%tran_xx,ierr);CHKERRQ(ierr)
  else
    call SNESSolve(solver%snes, PETSC_NULL_VEC, field%tran_xx,  &
                   ierr);CHKERRQ(ierr)
  endif

  ! do we really need all this? - geh 
  call SNESGetIterationNumber(solver%snes,num_newton_iterations,  &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,  &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

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
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(field%tran_r,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
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
subroutine StepperUpdateSolution(realization,update_kinetics)

!
! StepperUpdateSolution: Updates the solution variables
! Author: Glenn Hammond
! Date: 02/19/08 
!
  use Realization_class
  use Option_module

  implicit none
  
  class(realization_type) :: realization
  PetscBool :: update_kinetics
  
  ! update solution variables
  call RealizationUpdate(realization)
  if (realization%option%nflowdof > 0) &
    call StepperUpdateFlowSolution(realization)
  if (realization%option%ntrandof > 0) &
    call StepperUpdateTransportSolution(realization,update_kinetics)

end subroutine StepperUpdateSolution

! ************************************************************************** !

subroutine StepperUpdateFlowSolution(realization)
  ! 
  ! Updates the flow solution variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 
  
  use Flash2_module, only: Flash2UpdateSolution
  use Mphase_module, only: MphaseUpdateSolution
  use Immis_module, only: ImmisUpdateSolution
  use Miscible_module, only: MiscibleUpdateSolution 
  use Richards_module, only : RichardsUpdateSolution
  use TH_module, only : THUpdateSolution
  use General_module, only : GeneralUpdateSolution

  use Realization_class
  use Option_module

  implicit none

  class(realization_type) :: realization

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
    case(TH_MODE)
      call THUpdateSolution(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateSolution(realization)
    case(G_MODE)
      call GeneralUpdateSolution(realization)
  end select    

end subroutine StepperUpdateFlowSolution

! ************************************************************************** !

subroutine StepperUpdateTransportSolution(realization,update_kinetics)
  ! 
  ! Updates the transport solution variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/19/08
  ! 

  use Realization_class
  use Reactive_Transport_module, only : RTUpdateEquilibriumState, &
                                        RTUpdateKineticState, &
                                        RTUpdateMassBalance
  implicit none

  class(realization_type) :: realization
  PetscBool :: update_kinetics

  PetscErrorCode :: ierr
  
  call RTUpdateEquilibriumState(realization)
  if (update_kinetics) &
    call RTUpdateKineticState(realization)
  if (realization%reaction%update_porosity .or. &
      realization%reaction%update_tortuosity .or. &
      realization%reaction%update_permeability .or. &
      realization%reaction%update_mineral_surface_area) then
    call RealizationUpdateProperties(realization)
  endif
  
  if (realization%option%compute_mass_balance_new) then
    call RTUpdateMassBalance(realization)
  endif  

end subroutine StepperUpdateTransportSolution

! ************************************************************************** !

subroutine StepperJumpStart(realization)
  ! 
  ! Sets kinetic sorbed concentrations
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/09
  ! 

  use Realization_class
  use Reactive_Transport_module, only : RTJumpStartKineticSorption

  implicit none

  class(realization_type) :: realization

  call RTJumpStartKineticSorption(realization)

end subroutine StepperJumpStart

! ************************************************************************** !

subroutine StepperUpdateFlowAuxVars(realization)
  ! 
  ! Updates the flow auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/08
  ! 
  
  use Flash2_module, only: Flash2UpdateAuxVars
  use Mphase_module, only: MphaseUpdateAuxVars
  use Immis_module, only: ImmisUpdateAuxVars
  use Miscible_module, only: MiscibleUpdateAuxVars
  use Richards_module, only : RichardsUpdateAuxVars
  use TH_module, only : THUpdateAuxVars
  use General_module, only : GeneralUpdateAuxVars

  use Realization_class
  use Option_module

  implicit none

  class(realization_type) :: realization

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
    case(TH_MODE)
      call THUpdateAuxVars(realization)
    case(RICHARDS_MODE)
      call RichardsUpdateAuxVars(realization)
    case(G_MODE)
                                            ! do not update state
      call GeneralUpdateAuxVars(realization,PETSC_FALSE)
  end select    

end subroutine StepperUpdateFlowAuxVars

! ************************************************************************** !

subroutine StepperUpdateTranAuxVars(realization)
  ! 
  ! Updates the flow auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/08
  ! 
  
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Realization_class

  implicit none

  class(realization_type) :: realization

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)

end subroutine StepperUpdateTranAuxVars

! ************************************************************************** !

subroutine StepperSandbox(realization)
  ! 
  ! Sandbox for temporary miscellaneous operations
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/27/11
  ! 
  
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Realization_class
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
  use Secondary_Continuum_Aux_module

  implicit none

  class(realization_type) :: realization
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option

  PetscReal, pointer :: tran_xx_p(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
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

  rt_auxvars => patch%Aux%RT%auxvars
  global_auxvars => patch%Aux%Global%auxvars
  rt_sec_transport_vars => patch%Aux%SC_RT%sec_transport_vars

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

  call VecGetArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  
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
      vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
    endif
    
    tran_xx_p(istart:iend) = rt_auxvars(ghosted_id)%total(:,ONE_INTEGER)
    ! scale uo2++ total concentration by 0.25
    tran_xx_p(istart + species_offset) = 0.25d0 * &
                                         tran_xx_p(istart + species_offset)

  option%io_buffer = 'StepperSandbox() no longer used.'
  call printErrMsg(option)
!geh: this subroutine no longer used
!    call RReact(rt_auxvars(ghosted_id),global_auxvars(ghosted_id), &
!                tran_xx_p(istart:iend),volume_p(local_id), &
!                porosity_loc_p(ghosted_id), &
!                num_iterations,reaction,option,vol_frac_prim)
    tran_xx_p(istart:iend) = rt_auxvars(ghosted_id)%pri_molal
  enddo

  call VecRestoreArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)

                                   ! cells     bcs        act coefs.
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)

end subroutine StepperSandbox

! ************************************************************************** !

subroutine StepperCheckpoint(realization,flow_timestepper,tran_timestepper,id, id_stamp)
  ! 
  ! Calls appropriate routines to write a checkpoint file
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/08
  ! 

  use Realization_class
  use Checkpoint_module
  use Option_module
  use String_module, only : StringNull

  implicit none

  class(realization_type) :: realization
  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
  PetscInt :: num_const_timesteps, num_newton_iterations  
  PetscInt :: id
  character(len=MAXWORDLENGTH), optional, intent(in) :: id_stamp

  type(option_type), pointer :: option
  PetscInt :: flow_steps, flow_cumulative_newton_iterations, &
              flow_cumulative_time_step_cuts, flow_cumulative_linear_iterations, &
              flow_num_const_time_steps, flow_num_newton_iterations
  PetscInt :: tran_steps, tran_cumulative_newton_iterations, &
              tran_cumulative_time_step_cuts, tran_cumulative_linear_iterations, &
              tran_num_const_time_steps, tran_num_newton_iterations
  PetscReal :: flow_cumulative_solver_time, flow_prev_dt
  PetscReal :: tran_cumulative_solver_time,tran_prev_dt
  character(len=MAXWORDLENGTH) :: id_string
  
  option => realization%option

  if (associated(flow_timestepper)) then
    flow_steps = flow_timestepper%steps
    flow_cumulative_newton_iterations = flow_timestepper%cumulative_newton_iterations
    flow_cumulative_time_step_cuts = flow_timestepper%cumulative_time_step_cuts
    flow_cumulative_linear_iterations = flow_timestepper%cumulative_linear_iterations
    flow_num_const_time_steps = flow_timestepper%num_constant_time_steps
    flow_num_newton_iterations = flow_timestepper%num_newton_iterations
    flow_cumulative_solver_time = flow_timestepper%cumulative_solver_time
    flow_prev_dt = flow_timestepper%prev_dt
  endif
  if (associated(tran_timestepper)) then
    tran_steps = tran_timestepper%steps
    tran_cumulative_newton_iterations = tran_timestepper%cumulative_newton_iterations
    tran_cumulative_time_step_cuts = tran_timestepper%cumulative_time_step_cuts
    tran_cumulative_linear_iterations = tran_timestepper%cumulative_linear_iterations
    tran_num_const_time_steps = tran_timestepper%num_constant_time_steps
    tran_num_newton_iterations = tran_timestepper%num_newton_iterations
    tran_cumulative_solver_time = tran_timestepper%cumulative_solver_time
    tran_prev_dt = tran_timestepper%prev_dt
  endif
  
  ! default null id_string --> global_prefix-restart.chk
  id_string = ''
  if (present(id_stamp)) then
     id_string = id_stamp
  else if (id >= 0) then
     write(id_string,'(i8)') id
  end if

  call Checkpoint(realization, &
                  flow_steps,flow_cumulative_newton_iterations, &
                  flow_cumulative_time_step_cuts,flow_cumulative_linear_iterations, &
                  flow_num_const_time_steps,flow_num_newton_iterations, &
                  flow_cumulative_solver_time,flow_prev_dt, &
                  tran_steps,tran_cumulative_newton_iterations, &
                  tran_cumulative_time_step_cuts,tran_cumulative_linear_iterations, &
                  tran_num_const_time_steps,tran_num_newton_iterations, &
                  tran_cumulative_solver_time,tran_prev_dt, &
                  id_string)
                      
end subroutine StepperCheckpoint

! ************************************************************************** !

subroutine TimestepperRestart(realization,flow_timestepper,tran_timestepper, &
                              flow_read,transport_read,activity_coefs_read)
  ! 
  ! Calls appropriate routines to read checkpoint file and
  ! restart
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/08
  ! 

  use Realization_class
  use Checkpoint_module
  use Option_module

  implicit none

  class(realization_type) :: realization
  type(timestepper_type), pointer :: flow_timestepper
  type(timestepper_type), pointer :: tran_timestepper
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
  if (Uninitialized(option%restart_time)) then
    option%time = max(option%flow_time,option%tran_time)
    if (associated(flow_timestepper) .and. flow_read) then
      flow_timestepper%steps = flow_steps
      flow_timestepper%cumulative_newton_iterations = flow_cumulative_newton_iterations
      flow_timestepper%cumulative_time_step_cuts = flow_cumulative_time_step_cuts
      flow_timestepper%cumulative_linear_iterations = flow_cumulative_linear_iterations
      flow_timestepper%cumulative_solver_time = flow_cum_solver_time
      flow_timestepper%prev_dt = flow_prev_dt
      if (.not.flow_timestepper%run_as_steady_state) then
        flow_timestepper%num_constant_time_steps = flow_num_constant_time_steps
        flow_timestepper%num_newton_iterations = flow_num_newton_iterations
      endif
    endif
    if (associated(tran_timestepper) .and. transport_read) then
      tran_timestepper%steps = tran_steps
      tran_timestepper%cumulative_newton_iterations = tran_cumulative_newton_iterations
      tran_timestepper%cumulative_time_step_cuts = tran_cumulative_time_step_cuts
      tran_timestepper%cumulative_linear_iterations = tran_cumulative_linear_iterations
      tran_timestepper%cumulative_solver_time = tran_cum_solver_time
      tran_timestepper%num_constant_time_steps = tran_num_constant_time_steps
      tran_timestepper%num_newton_iterations = tran_num_newton_iterations
      tran_timestepper%prev_dt = tran_prev_dt
    endif
  else
    option%time = option%restart_time
    option%flow_time = option%restart_time
    option%tran_time = option%restart_time
    if (associated(flow_timestepper)) then
      option%flow_dt = flow_timestepper%dt_init
      flow_timestepper%steps = 0
      flow_timestepper%cumulative_newton_iterations = 0
      flow_timestepper%cumulative_time_step_cuts = 0
      flow_timestepper%cumulative_linear_iterations = 0
      flow_timestepper%cumulative_solver_time = 0.d0
      flow_timestepper%num_constant_time_steps = 0
      flow_timestepper%num_newton_iterations = 0
      flow_timestepper%prev_dt = 0.d0
    endif
    if (associated(tran_timestepper)) then
      option%tran_dt = tran_timestepper%dt_init
      tran_timestepper%steps = 0
      tran_timestepper%cumulative_newton_iterations = 0
      tran_timestepper%cumulative_time_step_cuts = 0
      tran_timestepper%cumulative_linear_iterations = 0
      tran_timestepper%cumulative_solver_time = 0.d0
      tran_timestepper%num_constant_time_steps = 0
      tran_timestepper%num_newton_iterations = 0
      tran_timestepper%prev_dt = 0.d0
    endif
    option%match_waypoint = PETSC_FALSE
    realization%output_option%plot_number = 0
  endif
  
  if (associated(flow_timestepper)) flow_timestepper%cur_waypoint => &
    WaypointReturnAtTime(realization%waypoint_list,option%time)
  if (associated(tran_timestepper)) tran_timestepper%cur_waypoint => &
    WaypointReturnAtTime(realization%waypoint_list,option%time)

  if (flow_read) then
    flow_timestepper%target_time = option%flow_time
    call StepperUpdateFlowAuxVars(realization)
  endif

  if (transport_read) then
    tran_timestepper%target_time = option%tran_time
    ! This is here since we need to recalculate the secondary complexes
    ! if they exist.  DO NOT update activity coefficients!!! - geh
    if (realization%reaction%use_full_geochemistry) then
      call StepperUpdateTranAuxVars(realization)
      !call StepperSandbox(realization)
    endif
  endif  
    
end subroutine TimestepperRestart

! ************************************************************************** !

subroutine TimestepperSetTranWeights(option,flow_t0,flow_t1)
  ! 
  ! TimestepperGetTranWeight: Sets the weights at t0 or t1 for transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/17/11
  ! 

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

subroutine TimestepperCheckCFLLimit(stepper,realization)
  ! 
  ! Checks CFL limit specified by the user
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/17/11
  ! 

  use Realization_class
  
  implicit none

  type(timestepper_type) :: stepper
  class(realization_type) :: realization
  
  PetscReal :: dt_cfl_1
  
  call RealizationCalculateCFL1Timestep(realization,dt_cfl_1)
  stepper%cfl_limiter_ts = dt_cfl_1*stepper%cfl_limiter
 
end subroutine TimestepperCheckCFLLimit

! ************************************************************************** !

subroutine TimestepperEnforceCFLLimit(stepper,option,output_option)
  ! 
  ! Enforces a CFL limit specified by the user
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/17/11
  ! 

  use Option_module
  use Output_Aux_module

  implicit none

  type(timestepper_type) :: stepper
  type(option_type) :: option
  type(output_option_type) :: output_option
  
  if (stepper%cfl_limiter_ts < option%tran_dt) then
    option%tran_dt = stepper%cfl_limiter_ts
    if (OptionPrintToScreen(option)) then
      write(*,'(" CFL Limiting: ",1pe12.4," [",a,"]")') &
            stepper%cfl_limiter_ts/output_option%tconv,trim(output_option%tunit)
    endif
    if (OptionPrintToFile(option)) then
      write(option%fid_out,'(/," CFL Limiting: ",1pe12.4," [",a,"]",/)') &
            stepper%cfl_limiter_ts/output_option%tconv,trim(output_option%tunit)
    endif        
  endif    

end subroutine TimestepperEnforceCFLLimit
#endif ! ifndef PROCESS_MODEL

! ************************************************************************** !

subroutine TimestepperPrintInfo(stepper,fid,header,option)
  ! 
  ! Prints information about time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  
  implicit none
  
  type(timestepper_type) :: stepper
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

subroutine TimestepperDestroy(stepper)
  ! 
  ! Deallocates a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(timestepper_type), pointer :: stepper
  
  if (.not.associated(stepper)) return
    
  call SolverDestroy(stepper%solver)
  call ConvergenceContextDestroy(stepper%convergence_context)
  
  if (associated(stepper%tfac)) deallocate(stepper%tfac)
  nullify(stepper%tfac)
  deallocate(stepper)
  nullify(stepper)
  
end subroutine TimestepperDestroy

end module Timestepper_module
