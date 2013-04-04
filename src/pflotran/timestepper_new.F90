module Timestepper_module
 
  use Solver_module
  use Waypoint_module 
  use Convergence_module 
 
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter, public :: TIMESTEPPER_INIT_PROCEED = 0
  PetscInt, parameter, public :: TIMESTEPPER_INIT_DONE = 1
  PetscInt, parameter, public :: TIMESTEPPER_INIT_FAIL = 2
 
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
    
    PetscReal :: dt
    PetscReal :: prev_dt
    PetscReal :: dt_min
    PetscReal :: dt_max
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
    type(waypoint_type), pointer :: prev_waypoint
    PetscBool :: match_waypoint

    type(convergence_context_type), pointer :: convergence_context

  end type stepper_type
  
  public :: TimestepperCreate, TimestepperDestroy, &
            StepperSetTargetTime, StepperStepDT, &
            StepperUpdateDT, &
            TimestepperRead, TimestepperPrintInfo

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

  stepper%start_time = 0.d0  
  stepper%start_time_step = 0
  stepper%time_step_tolerance = 0.1d0
  stepper%target_time = 0.d0
  
  stepper%prev_dt = 0.d0
  stepper%dt = 1.d0
  stepper%dt_min = 1.d0
  stepper%dt_max = 3.1536d6 ! One-tenth of a year.  
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
  nullify(stepper%prev_waypoint)
  stepper%match_waypoint = PETSC_FALSE
  
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
  stepper%match_waypoint = PETSC_FALSE

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
! StepperUpdateDT: Updates time step
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine StepperUpdateDT(timestepper,process_model)

  use Process_Model_Base_class
  
  implicit none

  type(stepper_type) :: timestepper
  class(process_model_base_type) :: process_model
  
  PetscBool :: update_time_step
  
  update_time_step = PETSC_TRUE
  if (timestepper%time_step_cut_flag) then
    timestepper%num_constant_time_steps = 1
  else if (timestepper%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    timestepper%num_constant_time_steps = &
      timestepper%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (timestepper%num_constant_time_steps > &
      timestepper%constant_time_step_threshold) then
    timestepper%num_constant_time_steps = 0
  else if (timestepper%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif
    
  if (update_time_step .and. timestepper%iaccel /= 0) then
      
    call process_model%UpdateTimestep(timestepper%dt, &
                                      timestepper%dt_max, &
                                      timestepper%iaccel, &
                                      timestepper%num_newton_iterations, &
                                      timestepper%tfac)
    
  endif

end subroutine StepperUpdateDT

! ************************************************************************** !
!
! StepperSetTargetTime: Sets target time for timestepper
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine StepperSetTargetTime(timestepper,sync_time,option,stop_flag, &
                                plot_flag,transient_plot_flag)

  use Option_module
  
  implicit none

  type(stepper_type), pointer :: timestepper
  PetscReal :: sync_time
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  
  PetscReal :: target_time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: cumulative_time_steps
  PetscInt :: max_time_step
  PetscReal :: max_time
  PetscReal :: tolerance
  PetscBool :: force_to_match_waypoint
  PetscBool :: equal_to_or_exceeds_waypoint
  PetscBool :: equal_to_or_exceeds_sync_time
  type(waypoint_type), pointer :: cur_waypoint

  option%io_buffer = 'StepperSetTargetTime()'
  call printMsg(option)
  
  if (timestepper%time_step_cut_flag) then
    timestepper%time_step_cut_flag = PETSC_FALSE
    timestepper%match_waypoint = PETSC_FALSE ! reset back to false
    timestepper%cur_waypoint => timestepper%prev_waypoint
  else
    ! if the maximum time step size decreased in the past step, need to set
    ! the time step size to the minimum of the stepper%prev_dt and stepper%dt_max
    if (timestepper%match_waypoint) then
      timestepper%dt = min(timestepper%prev_dt,timestepper%dt_max)
      timestepper%match_waypoint = PETSC_FALSE
    endif
  endif
  
  dt = timestepper%dt
  timestepper%prev_dt = dt
  cur_waypoint => timestepper%cur_waypoint
  ! need previous waypoint for reverting back on time step cut
  timestepper%prev_waypoint => timestepper%cur_waypoint
  ! dt_max must be set from current waypoint and not updated below
  dt_max = cur_waypoint%dt_max
  cumulative_time_steps = timestepper%steps
  max_time_step = timestepper%max_time_step
  tolerance = timestepper%time_step_tolerance
  target_time = timestepper%target_time + dt

  !TODO(geh): move to process model initialization stage
  ! For the case where the second waypoint is a printout after the first time step
  ! we must increment the waypoint beyond the first (time=0.) waypoint.  Otherwise
  ! the second time step will be zero. - geh
  if (cur_waypoint%time < 1.d-40) then
    cur_waypoint => cur_waypoint%next
  endif
  
  ! If a waypoint calls for a plot or change in src/sinks, adjust time step
  ! to match waypoint.
  force_to_match_waypoint = WaypointForceMatchToTime(cur_waypoint)
  equal_to_or_exceeds_waypoint = target_time + tolerance*dt >= cur_waypoint%time
  equal_to_or_exceeds_sync_time = target_time + tolerance*dt >= sync_time
  do ! we cycle just in case the next waypoint is beyond the target_time
    if (equal_to_or_exceeds_sync_time .or. &
        (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint)) then
      max_time = min(sync_time,cur_waypoint%time)
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
      if (dt > dt_max .and. dabs(dt-dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
        dt = dt_max                    ! error from waypoint%time - time
        target_time = target_time + dt
      else
        target_time = max_time
        if (equal_to_or_exceeds_waypoint) then
          if (force_to_match_waypoint) timestepper%match_waypoint = PETSC_TRUE
          if (cur_waypoint%print_output) plot_flag = PETSC_TRUE
          if (cur_waypoint%print_tr_output) transient_plot_flag = PETSC_TRUE
        endif
        if (max_time >= cur_waypoint%time) then
          cur_waypoint => cur_waypoint%next
        endif
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
    nullify(cur_waypoint)
  endif

  ! update maximum time step size to current waypoint value
  if (associated(cur_waypoint)) then
    dt_max = cur_waypoint%dt_max
  else
    stop_flag = 1 ! stop after end of time step
  endif
  
  option%refactor_dt = dt
  timestepper%dt = dt
  timestepper%dt_max = dt_max
  timestepper%target_time = target_time
  timestepper%cur_waypoint => cur_waypoint
  
 end subroutine StepperSetTargetTime

! ************************************************************************** !
!
! StepperStepDT: Steps forward one step in time
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine StepperStepDT(timestepper,process_model,stop_flag)

  use Process_Model_Base_class
  use Option_module
  use Output_module, only : Output
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  type(stepper_type) :: timestepper
  class(process_model_base_type) :: process_model
  PetscInt :: stop_flag
  
  SNESConvergedReason :: snes_reason
  PetscInt :: icut
  
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: sum_newton_iterations
  PetscInt :: sum_linear_iterations
  character(len=2) :: tunit
  PetscReal :: tconv
  PetscReal :: fnorm, inorm, scaled_fnorm
  PetscBool :: plot_flag, transient_plot_flag
  PetscErrorCode :: ierr
  
  solver => timestepper%solver
  option => process_model%option
  
  write(process_model%option%io_buffer,'(f12.2)') timestepper%dt
  process_model%option%io_buffer = 'StepperStepDT(' // &
    trim(adjustl(process_model%option%io_buffer)) // ')'
  call printMsg(process_model%option)  

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_newton_iterations = 0
  icut = 0
  
  option%dt = timestepper%dt
  option%time = timestepper%target_time-timestepper%dt

  call process_model%InitializeTimestep()
    
  do
      
    call process_model%PreSolve()
    
    call PetscTime(log_start_time, ierr)

    call SNESSolve(solver%snes,PETSC_NULL_OBJECT, &
                   process_model%solution_vec,ierr)

    call PetscTime(log_end_time, ierr)

    timestepper%cumulative_solver_time = &
      timestepper%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call SNESGetIterationNumber(solver%snes,num_newton_iterations,ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations,ierr)
    call SNESGetConvergedReason(solver%snes,snes_reason,ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
    if (snes_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      timestepper%time_step_cut_flag = PETSC_TRUE

      if (icut > timestepper%max_time_step_cuts .or. &
          timestepper%dt < 1.d-20) then
        if (option%print_screen_flag) then
          print *,"--> max_time_step_cuts exceeded: icut/icutmax= ",icut, &
                  timestepper%max_time_step_cuts, "t= ", &
                  timestepper%target_time/tconv, &
                  " dt= ", &
                  timestepper%dt/tconv
          print *,"Stopping execution!"
        endif
        process_model%output_option%plot_name = 'flow_cut_to_failure'
        plot_flag = PETSC_TRUE
        transient_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,plot_flag,transient_plot_flag)
        stop_flag = 2
        return
      endif
 
      timestepper%target_time = timestepper%target_time - timestepper%dt

      timestepper%dt = 0.5d0 * timestepper%dt  
      
      if (option%print_screen_flag) write(*,'('' -> Cut time step: snes='',i3, &
        &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
        &   1pe12.5)')  snes_reason,icut,timestepper%cumulative_time_step_cuts, &
            option%time/tconv, &
            timestepper%dt/tconv

      timestepper%target_time = timestepper%target_time + timestepper%dt
      option%dt = timestepper%dt
      call process_model%TimeCut()
  
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  timestepper%steps = timestepper%steps + 1      
  timestepper%cumulative_newton_iterations = &
    timestepper%cumulative_newton_iterations + sum_newton_iterations
  timestepper%cumulative_linear_iterations = &
    timestepper%cumulative_linear_iterations + sum_linear_iterations
  timestepper%cumulative_time_step_cuts = &
    timestepper%cumulative_time_step_cuts + icut

  timestepper%num_newton_iterations = num_newton_iterations
  timestepper%num_linear_iterations = num_linear_iterations  
  
! print screen output
  call SNESGetFunctionNorm(solver%snes,fnorm,ierr)
  call VecNorm(process_model%residual_vec,NORM_INFINITY,inorm,ierr)
  if (option%print_screen_flag) then
    write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]", &
      & " snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      timestepper%steps, &
      timestepper%target_time/tconv, &
      timestepper%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      timestepper%cumulative_newton_iterations,sum_linear_iterations, &
      timestepper%cumulative_linear_iterations,icut, &
      timestepper%cumulative_time_step_cuts

#if 0    
    if (associated(discretization%grid)) then
       scaled_fnorm = fnorm/discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif
#endif
    scaled_fnorm = fnorm

    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  if (option%print_file_flag) then
    write(option%fid_out, '(" Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a1, &
      & "]"," snes_conv_reason: ",i4,/,"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      timestepper%steps, &
      timestepper%target_time/tconv, &
      timestepper%dt/tconv, &
      tunit,snes_reason,sum_newton_iterations, &
      timestepper%cumulative_newton_iterations,sum_linear_iterations, &
      timestepper%cumulative_linear_iterations,icut, &
      timestepper%cumulative_time_step_cuts
  endif  
  
  option%time = timestepper%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""  
  
end subroutine StepperStepDT

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
