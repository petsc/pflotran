module Simulation_module

  use Realization_module
  use Timestepper_module
  use Solver_module

  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type

    type(realization_type), pointer :: realization
    type(stepper_type), pointer :: flow_stepper
    type(stepper_type), pointer :: tran_stepper

  end type simulation_type
  
  public :: SimulationCreate, SimulationDestroy
  
contains

! ************************************************************************** !
!
! SimulationCreate: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate()

  use Option_module

  implicit none
  
  type(simulation_type), pointer :: SimulationCreate
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
  simulation%realization => RealizationCreate()
  simulation%flow_stepper => TimestepperCreate()
  simulation%tran_stepper => TimestepperCreate()
  
  SimulationCreate => simulation
  
end function SimulationCreate  

! ************************************************************************** !
!
! SimulationDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SimulationDestroy(simulation)

  use Richards_module, only : RichardsDestroy
  use Reactive_Transport_module, only : RTDestroy

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return

  if (simulation%realization%option%nflowdof > 0) then
    select case(simulation%realization%option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsDestroy(simulation%realization)
    end select
  endif

  if (simulation%realization%option%ntrandof > 0) then
    call RTDestroy(simulation%realization)
  endif

  call RealizationDestroy(simulation%realization)
  call TimestepperDestroy(simulation%flow_stepper)
  call TimestepperDestroy(simulation%tran_stepper)
  
end subroutine SimulationDestroy
  
end module Simulation_module
