module Timestepper_Base_class
 
  use Waypoint_module 
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"
 
  type, public :: stepper_base_type
  
    PetscInt :: steps         ! The number of time steps taken by the code.
    PetscInt :: num_constant_time_steps   ! number of contiguous time_steps of constant size

    PetscInt :: max_time_step                ! Maximum number of time steps to be taken by the code.
    PetscInt :: max_time_step_cuts           ! Maximum number of timestep cuts within one time step.
    PetscInt :: constant_time_step_threshold        ! Steps needed after cutting to increase time step

    PetscInt :: cumulative_time_step_cuts       ! Total number of cuts in the timestep taken.    
    PetscReal :: cumulative_solver_time
    
    PetscReal :: dt
    PetscReal :: prev_dt
    PetscReal :: dt_min
    PetscReal :: dt_max
    PetscReal :: cfl_limiter
    PetscReal :: cfl_limiter_ts
    PetscBool :: revert_dt
    PetscInt :: num_contig_revert_due_to_sync
    
    PetscBool :: init_to_steady_state
    PetscBool :: run_as_steady_state
    PetscReal :: steady_state_rel_tol
    
    PetscBool :: time_step_cut_flag  ! flag toggled if timestep is cut

    PetscLogDouble :: start_time    
    PetscInt :: start_time_step ! the first time step of a given run
    PetscReal :: time_step_tolerance ! scalar used in determining time step size
    PetscReal :: target_time    ! time at end of "synchronized" time step 

    type(waypoint_type), pointer :: cur_waypoint
    type(waypoint_type), pointer :: prev_waypoint

  contains
    
    procedure, public :: ReadInput => TimestepperBaseRead
    procedure, public :: Init => TimestepperBaseInit
    procedure, public :: SetTargetTime => TimestepperBaseSetTargetTime
    procedure, public :: StepDT => TimestepperBaseStepDT
    procedure, public :: UpdateDT => TimestepperBaseUpdateDT
    procedure, public :: Checkpoint => TimestepperBaseCheckpoint
    procedure, public :: Restart => TimestepperBaseRestart
    procedure, public :: FinalizeRun => TimestepperBaseFinalizeRun
    procedure, public :: Strip => TimestepperBaseStrip
    procedure, public :: Destroy => TimestepperBaseDestroy
    
  end type stepper_base_type
  
  type, public :: stepper_base_header_type
    real*8 :: time
    real*8 :: dt
    real*8 :: prev_dt
    integer*8 :: num_steps
    integer*8 :: cumulative_time_step_cuts
    integer*8 :: num_constant_time_steps
    integer*8 :: num_contig_revert_due_to_sync
    integer*8 :: revert_dt
  end type stepper_base_header_type
  
  public :: TimestepperBaseCreate, TimestepperBasePrintInfo, &
            TimestepperBaseProcessKeyword, &
            TimestepperBaseStrip, &
            TimestepperBaseInit, &
            TimestepperBaseSetHeader, &
            TimestepperBaseGetHeader, &
            TimestepperBaseRegisterHeader

contains

! ************************************************************************** !
!
! TimestepperBaseCreate: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function TimestepperBaseCreate()

  implicit none
  
  class(stepper_base_type), pointer :: TimestepperBaseCreate
  
  class(stepper_base_type), pointer :: this
  
  allocate(this)
  call this%Init()
  
  TimestepperBaseCreate => this
  
end function TimestepperBaseCreate

! ************************************************************************** !
!
! TimestepperBaseInit: Allocates and initializes a new Timestepper object
! author: Glenn Hammond
! date: 07/01/13
!
! ************************************************************************** !
subroutine TimestepperBaseInit(this)

  implicit none
  
  class(stepper_base_type) :: this
  
  this%steps = 0
  this%num_constant_time_steps = 0

  this%max_time_step = 999999
  this%max_time_step_cuts = 16
  this%constant_time_step_threshold = 5

  this%cumulative_time_step_cuts = 0    
  this%cumulative_solver_time = 0.d0

  this%start_time = 0.d0  
  this%start_time_step = 0
  this%time_step_tolerance = 0.1d0
  this%target_time = 0.d0
  
  this%prev_dt = 0.d0
  this%dt = 1.d0
  this%dt_min = 1.d0
  this%dt_max = 3.1536d6 ! One-tenth of a year.  
  this%cfl_limiter = -999.d0
  this%cfl_limiter_ts = 1.d20
  
  this%time_step_cut_flag = PETSC_FALSE
  
  this%init_to_steady_state = PETSC_FALSE
  this%steady_state_rel_tol = 1.d-8
  this%run_as_steady_state = PETSC_FALSE
  
  nullify(this%cur_waypoint)
  nullify(this%prev_waypoint)
  this%revert_dt = PETSC_FALSE
  this%num_contig_revert_due_to_sync = 0
  
end subroutine TimestepperBaseInit

! ************************************************************************** !
!
! TimestepperBaseRead: Reads parameters associated with time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperBaseRead(this,input,option)

  use Option_module
  use Input_module
  
  implicit none

  class(stepper_base_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseRead not supported.  Requires extension.'
  call printErrMsg(option)

end subroutine TimestepperBaseRead

! ************************************************************************** !
!
! TimestepperBaseProcessKeyword: Updates time step
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine TimestepperBaseProcessKeyword(this,input,option,keyword)

  use Option_module
  use String_module
  use Input_module
  
  implicit none
  
  class(stepper_base_type) :: this
  character(len=MAXWORDLENGTH) :: keyword
  type(input_type) :: input
  type(option_type) :: option

  select case(trim(keyword))

    case('NUM_STEPS_AFTER_TS_CUT')
      call InputReadInt(input,option,this%constant_time_step_threshold)
      call InputDefaultMsg(input,option,'num_constant_time_steps_after_ts_cut')

    case('MAX_STEPS')
      call InputReadInt(input,option,this%max_time_step)
      call InputDefaultMsg(input,option,'max_time_step')
  
    case('MAX_TS_CUTS')
      call InputReadInt(input,option,this%max_time_step_cuts)
      call InputDefaultMsg(input,option,'max_time_step_cuts')
        
    case('CFL_LIMITER')
      call InputReadDouble(input,option,this%cfl_limiter)
      call InputDefaultMsg(input,option,'cfl limiter')

    case('INITIALIZE_TO_STEADY_STATE')
      this%init_to_steady_state = PETSC_TRUE
      call InputReadDouble(input,option,this%steady_state_rel_tol)
      call InputDefaultMsg(input,option,'steady state convergence relative tolerance')

    case('RUN_AS_STEADY_STATE')
      this%run_as_steady_state = PETSC_TRUE

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

end subroutine TimestepperBaseProcessKeyword

! ************************************************************************** !
!
! TimestepperBaseUpdateDT: Updates time step
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine TimestepperBaseUpdateDT(this,process_model)

  use Process_Model_Base_class
  use Option_module
  
  implicit none

  class(stepper_base_type) :: this
  class(pm_base_type) :: process_model
  
  process_model%option%io_buffer = 'TimestepperBaseStepDT must be extended.'
  call printErrMsg(process_model%option)

end subroutine TimestepperBaseUpdateDT

! ************************************************************************** !
!
! TimestepperBaseSetTargetTime: Sets target time for timestepper
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine TimestepperBaseSetTargetTime(this,sync_time,option, &
                                        stop_flag,plot_flag, &
                                        transient_plot_flag)

  use Option_module
  
  implicit none

  class(stepper_base_type) :: this
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
  PetscBool :: revert_due_to_waypoint
  PetscBool :: revert_due_to_sync_time
  type(waypoint_type), pointer :: cur_waypoint

!geh: for debugging
!  option%io_buffer = 'StepperSetTargetTime()'
!  call printMsg(option)
  
  if (this%time_step_cut_flag) then
    this%time_step_cut_flag = PETSC_FALSE
    !geh: pointing the cur_waypoint back may cause problems in the checkpoint
    !     file.  There is no way of knowing whether prev_waypoint is different
    !     from cur_waypoint as most of the time it will be identical.  I believe
    !     the only way around this is to check associated(cur,prev) to see if
    !     they differ and set a flag in the checkpoint file.  But even this will
    !     not work if more than one waypoint previous.
    this%cur_waypoint => this%prev_waypoint
  else
    ! If the maximum time step size decreased in the past step, need to set
    ! the time step size to the minimum of the this%prev_dt and 
    ! this%dt_max.  However, if we have to revert twice in a row, throw 
    ! away the old time step and move on.
    if (this%revert_dt .and. &
        this%num_contig_revert_due_to_sync < 2) then
      this%dt = min(this%prev_dt,this%dt_max)
    endif
  endif
  this%revert_dt = PETSC_FALSE ! reset back to false
  revert_due_to_waypoint = PETSC_FALSE
  revert_due_to_sync_time = PETSC_FALSE
  
  dt = this%dt
  this%prev_dt = dt
  cur_waypoint => this%cur_waypoint
  ! need previous waypoint for reverting back on time step cut
  this%prev_waypoint => this%cur_waypoint
  ! dt_max must be set from current waypoint and not updated below
  dt_max = cur_waypoint%dt_max
  cumulative_time_steps = this%steps
  max_time_step = this%max_time_step
  tolerance = this%time_step_tolerance
  target_time = this%target_time + dt

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
      if (dt > dt_max .and. &
          dabs(dt-dt_max) > 1.d0) then ! 1 sec tolerance to avoid cancellation
        dt = dt_max                    ! error from waypoint%time - time
        target_time = target_time + dt
      else
        target_time = max_time
        if (equal_to_or_exceeds_waypoint) then
          ! Since the time step was cut to match the waypoint, we want to set 
          ! the time step back to its prior value after the waypoint is met.
          ! %revert_dt is a flag that does so above.
          if (force_to_match_waypoint) revert_due_to_waypoint = PETSC_TRUE
          if (cur_waypoint%print_output) plot_flag = PETSC_TRUE
          if (cur_waypoint%print_tr_output) transient_plot_flag = PETSC_TRUE
        endif
        if (equal_to_or_exceeds_sync_time) then
          ! If the time step was cut to match the sync time, we want to set
          ! the time step back to its prior value.  However, if the time step
          ! is close to its full previous value, this constraint is unnecessary
          ! and limits the ability of process model couplers "below" to catch up
          ! with those above.  Thus the conditional (dt <= .5 prev_dt) below.
          !-Also note that if this timestepper is at a depth in the process 
          ! model coupler greater than 1 (not the top process model coupler)
          ! the timestepper will constantly be reverting to sync due to the
          ! tolerance applied above without the underlying conditional.
!          if (dt < 0.99d0 * this%prev_dt) then
          if (dt <= 0.5d0 * this%prev_dt) then
            revert_due_to_sync_time = PETSC_TRUE
          endif
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

  if (revert_due_to_sync_time .or. revert_due_to_waypoint) then
    this%revert_dt = PETSC_TRUE
    if (revert_due_to_sync_time) then
      this%num_contig_revert_due_to_sync = &
        this%num_contig_revert_due_to_sync + 1
    endif
  else
    this%num_contig_revert_due_to_sync = 0
  endif

  
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
  this%dt = dt
  this%dt_max = dt_max
  this%target_time = target_time
  this%cur_waypoint => cur_waypoint

 end subroutine TimestepperBaseSetTargetTime

! ************************************************************************** !
!
! TimestepperBaseStepDT: Steps forward one step in time
! author: Glenn Hammond
! date: 03/20/13
!
! ************************************************************************** !
subroutine TimestepperBaseStepDT(this,process_model,stop_flag)

  use Process_Model_Base_class
  use Option_module
  use Output_module, only : Output
  
  implicit none

  class(stepper_base_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag
  
  type(option_type), pointer :: option

  option => process_model%option
  
  option%io_buffer = 'TimestepperBaseStepDT must be extended.'
  call printErrMsg(option)
  
end subroutine TimestepperBaseStepDT

! ************************************************************************** !
!
! TimestepperBasePrintInfo: Prints information about time stepper
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine TimestepperBasePrintInfo(this,fid,header,option)

  use Option_module
  
  implicit none
  
  class(stepper_base_type) :: this
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header
  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(string,*) this%max_time_step
    write(*,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(*,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(*,'("max cuts:",x,a)') trim(adjustl(string))
  endif
  if (OptionPrintToFile(option)) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(string,*) this%max_time_step
    write(fid,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(fid,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(fid,'("max cuts:",x,a)') trim(adjustl(string))
  endif    

end subroutine TimestepperBasePrintInfo

! ************************************************************************** !
!
! TimestepperBaseCheckpoint: Checkpoints parameters/variables associated with 
!                            a time stepper.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBaseCheckpoint(this,viewer,option)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"

  class(stepper_base_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseCheckpoint must be extended.'
  call printErrMsg(option)  
    
end subroutine TimestepperBaseCheckpoint

! ************************************************************************** !
!
! TimestepperBaseRegisterHeader: Register header entries.
! author: Glenn Hammond
! date: 07/30/13
!
! ************************************************************************** !
subroutine TimestepperBaseRegisterHeader(this,bag,header)

  use Option_module

  implicit none
  
#include "finclude/petscbag.h"  

  class(stepper_base_type) :: this
  class(stepper_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  ! bagsize = 8 * 8 bytes = 64 bytes
  call PetscBagRegisterReal(bag,header%time,0,"time","",ierr)
  call PetscBagRegisterReal(bag,header%dt,0,"dt","",ierr)
  call PetscBagRegisterReal(bag,header%prev_dt,0,"prev_dt","",ierr)
  call PetscBagRegisterInt(bag,header%num_steps,0,"num_steps","",ierr)
  call PetscBagRegisterInt(bag,header%cumulative_time_step_cuts,0, &
                           "cumulative_time_step_cuts","",ierr)
  call PetscBagRegisterInt(bag,header%num_constant_time_steps,0, &
                           "num_constant_time_steps","",ierr)
  call PetscBagRegisterInt(bag,header%num_contig_revert_due_to_sync,0, &
                           "num_contig_revert_due_to_sync","",ierr)
  call PetscBagRegisterInt(bag,header%revert_dt,0, &
                           "revert_dt","",ierr)
    
end subroutine TimestepperBaseRegisterHeader

! ************************************************************************** !
!
! TimestepperBaseSetHeader: Sets values in checkpoint header.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBaseSetHeader(this,bag,header)

  use Option_module

  implicit none
  
#include "finclude/petscbag.h"  

  class(stepper_base_type) :: this
  class(stepper_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr

  header%time = this%target_time
  header%dt = this%dt
  header%prev_dt = this%prev_dt
  header%num_steps = this%steps
  header%cumulative_time_step_cuts = this%cumulative_time_step_cuts
  header%num_constant_time_steps = this%num_constant_time_steps
  header%num_contig_revert_due_to_sync = this%num_contig_revert_due_to_sync
  header%revert_dt = ZERO_INTEGER
  if (this%revert_dt) then
    header%revert_dt = ONE_INTEGER
  endif
    
end subroutine TimestepperBaseSetHeader

! ************************************************************************** !
!
! TimestepperBaseRestart: Restarts parameters/variables associated with 
!                         a time stepper.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBaseRestart(this,viewer,option)

  use Option_module

  implicit none

#include "finclude/petscviewer.h"

  class(stepper_base_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseRestart must be extended.'
  call printErrMsg(option)  
    
end subroutine TimestepperBaseRestart

! ************************************************************************** !
!
! TimestepperBaseGetHeader: Gets values in checkpoint header.
! author: Glenn Hammond
! date: 07/25/13
!
! ************************************************************************** !
subroutine TimestepperBaseGetHeader(this,header)

  use Option_module

  implicit none
  
#include "finclude/petscbag.h"  

  class(stepper_base_type) :: this
  class(stepper_base_header_type) :: header
  
  this%target_time = header%time
  this%dt = header%dt
  this%prev_dt = header%prev_dt
  this%steps = header%num_steps
  this%cumulative_time_step_cuts = header%cumulative_time_step_cuts
  this%num_constant_time_steps = header%num_constant_time_steps
  this%num_contig_revert_due_to_sync = header%num_contig_revert_due_to_sync
  this%revert_dt = (header%revert_dt == ONE_INTEGER)
    
end subroutine TimestepperBaseGetHeader

! ************************************************************************** !
!
! TimestepperBaseFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
recursive subroutine TimestepperBaseFinalizeRun(this,option)

  use Option_module
  
  implicit none
  
  class(stepper_base_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  call printMsg(option,'TSBE%FinalizeRun()')
  
  if (OptionPrintToScreen(option)) then
    write(*,'(/," TS Base steps = ",i6," cuts = ",i6)') &
            this%steps, &
            this%cumulative_time_step_cuts
    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,*) 'TS Base solver time = ' // trim(adjustl(string)) // ' seconds'
  endif
  
end subroutine TimestepperBaseFinalizeRun

! ************************************************************************** !
!
! TimestepperBaseStrip: Deallocates members of a time stepper
! author: Glenn Hammond
! date: 07/22/13
!
! ************************************************************************** !
subroutine TimestepperBaseStrip(this)

  implicit none
  
  class(stepper_base_type) :: this
  
end subroutine TimestepperBaseStrip

! ************************************************************************** !
!
! TimestepperBaseDestroy: Deallocates a time stepper
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine TimestepperBaseDestroy(this)

  implicit none
  
  class(stepper_base_type) :: this
  
  call TimestepperBaseStrip(this)
    
!  deallocate(this)
!  nullify(this)
  
end subroutine TimestepperBaseDestroy

end module Timestepper_Base_class
