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
    procedure, public :: SetTargetTime2 => TimeStepperSurfaceSetTargetTime
    procedure, public :: StepDT2 => TimeStepperSurfaceStepDT
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
subroutine TimeStepperSurfaceSetTargetTime(timestepper,sync_time,dt_max_allowable, &
                                        option, &
                                        stop_flag,plot_flag, &
                                        transient_plot_flag)

  use Option_module

  implicit none

  class(timestepper_surface_type) :: timestepper
  PetscReal :: sync_time
  PetscReal :: dt_max_allowable
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
  
  cur_waypoint => timestepper%cur_waypoint

  dt = min(dt_max_allowable,timestepper%dt_max)
  target_time = timestepper%target_time + dt
  tolerance = timestepper%time_step_tolerance

  plot_flag = PETSC_FALSE
  transient_plot_flag = PETSC_FALSE
  
  force_to_match_waypoint = WaypointForceMatchToTime(cur_waypoint)
  equal_to_or_exceeds_waypoint = target_time + tolerance*dt >= cur_waypoint%time
  equal_to_or_exceeds_sync_time = target_time + tolerance*dt >= sync_time

  if(equal_to_or_exceeds_sync_time .and. &
      (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint)) then

      max_time = min(sync_time,cur_waypoint%time)
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
      target_time = target_time + dt

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
      target_time = target_time + dt

    endif

    if (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint) then
      max_time = cur_waypoint%time
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time

      target_time = target_time + dt

      if(cur_waypoint%print_output) plot_flag = PETSC_TRUE
    endif
  
  endif

  timestepper%dt = dt
  timestepper%target_time = target_time
  timestepper%cur_waypoint => cur_waypoint

end subroutine TimeStepperSurfaceSetTargetTime

! ************************************************************************** !
!> This is a dummy routine added to be extended in timestepper_surface_type
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 07/03/13
! ************************************************************************** !
subroutine TimeStepperSurfaceStepDT(timestepper,process_model,stop_flag)

  use Process_Model_Base_class
  use Process_Model_Surface_Flow_class
  use Option_module
  use Output_module, only : Output
  use Surface_Flow_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  class(timestepper_surface_type) :: timestepper
  class(pm_surface_flow_type) :: process_model
  PetscInt :: stop_flag

  PetscReal :: time
  PetscReal :: dtime
  PetscReal :: tmp
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  solver => timestepper%solver
  option => process_model%option

  call process_model%PreSolve()
  
  if(option%subsurf_surf_coupling==SEQ_COUPLED .and. &
     associated(process_model%subsurf_realization)) then
    call SurfaceFlowSurf2SubsurfFlux(process_model%subsurf_realization, &
                                     process_model%surf_realization,tmp)
   endif
  
  call TSSetTimeStep(solver%ts,option%surf_flow_dt,ierr)
  call TSSolve(solver%ts,process_model%solution_vec,ierr)
  call TSGetTime(solver%ts,time,ierr)
  call TSGetTimeStep(solver%ts,dtime,ierr)

  call process_model%PostSolve()

  timestepper%steps = timestepper%steps + 1

  if (option%print_screen_flag) then
    write(*, '(" SURFACE FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]")') &
      timestepper%steps, &
      time/process_model%surf_realization%output_option%tconv, &
      dtime/process_model%surf_realization%output_option%tconv, &
      process_model%surf_realization%output_option%tunit
  endif

end subroutine TimeStepperSurfaceStepDT

end module Timestepper_Surface_class

#endif