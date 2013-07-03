#ifdef SURFACE_FLOW

module Timestepper_Surface_class

  use Timestepper_Base_class
  use Solver_module
  use Waypoint_module

  implicit none

#include "definitions.h"

  private

  type, public, extends(stepper_base_type) :: timestepper_surface_type
  contains
    procedure, public :: Init => TimeStepperSurfaceInit
    !procedure, public :: SetTargetTime => TimeStepperSurfaceSetTargetTime
    procedure, public :: TimeStepperSurfaceSetTargetTime
  end type timestepper_surface_type

  public TimeStepperSurfaceSetTargetTime, &
         TimeStepperSurfaceCreate

contains

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/03/13
! ************************************************************************** !
function TimeStepperSurfaceCreate()

  implicit none
  
  class(timestepper_surface_type), pointer :: TimeStepperSurfaceCreate
  
  class(timestepper_surface_type), pointer :: surf_stepper
  
  allocate(surf_stepper)
  call surf_stepper%Init()
  
  surf_stepper%solver => SolverCreate()
  
  TimeStepperSurfaceCreate => surf_stepper
  
end function TimeStepperSurfaceCreate

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/03/13
! ************************************************************************** !
subroutine TimeStepperSurfaceInit(stepper)

  implicit none
  
  class (timestepper_surface_type) :: stepper

  call TimestepperBaseInit(stepper)

end subroutine TimeStepperSurfaceInit

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/02/13
! ************************************************************************** !
subroutine TimeStepperSurfaceSetTargetTime(timestepper,sync_time,dt_max,option, &
                                        stop_flag,plot_flag, &
                                        transient_plot_flag)

  use Option_module

  implicit none

  class(timestepper_surface_type) :: timestepper
  PetscReal :: sync_time
  PetscReal :: dt_max
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag

  PetscReal :: dt
  PetscReal :: target_time
  PetscReal :: max_time
  PetscReal :: tolerance
  PetscBool :: equal_to_or_exceeds_waypoint
  PetscBool :: equal_to_or_exceeds_sync_time
  PetscBool :: force_to_match_waypoint
  type(waypoint_type), pointer :: cur_waypoint
  
  write(*,*),'HERE-1: ',timestepper%dt
  cur_waypoint => timestepper%cur_waypoint
  if (.not.associated(cur_waypoint)) write(*,*),'NOT'
  write(*,*),'HERE-2'

  dt = dt_max
  target_time = timestepper%target_time + dt
  tolerance = timestepper%time_step_tolerance

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  
  write(*,*),'WaypointForceMatchToTime: '
  force_to_match_waypoint = WaypointForceMatchToTime(cur_waypoint)
  write(*,*),'WaypointForceMatchToTime: -- done'
  equal_to_or_exceeds_waypoint = target_time + tolerance*dt >= cur_waypoint%time
  equal_to_or_exceeds_sync_time = target_time + tolerance*dt >= sync_time

  if(equal_to_or_exceeds_sync_time .and. &
      (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint)) then

      max_time = min(sync_time,cur_waypoint%time)
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time

      if(max_time == cur_waypoint%time) then
        if(cur_waypoint%print_output) plot_flag = PETSC_TRUE
      endif

  else

    if (equal_to_or_exceeds_sync_time) then
      max_time = sync_time
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
    endif

    if (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint) then
      max_time = cur_waypoint%time
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
      if(cur_waypoint%print_output) plot_flag = PETSC_TRUE
    endif
  
  endif
  
  timestepper%dt = dt
  timestepper%target_time = target_time
  timestepper%cur_waypoint => cur_waypoint
  write(*,*),'timestepper: ',dt,target_time


  option%io_buffer = 'TimeStepperSurfaceSetTargetTime()'
  call printMsg(option)

  call printErrMsg(option,'debugging in TimeStepperSurfaceSetTargetTime')

end subroutine TimeStepperSurfaceSetTargetTime

end module Timestepper_Surface_class

#endif