module Simulation_module

#ifndef TRANSPORT
  use Realization_module
  use Timestepper_module
  use Solver_module
#endif  
#ifdef TRANSPORT
  use Transport_Realization_module
  use Transport_Timestepper_module
  use Transport_Solver_module
#endif

  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type

#ifndef TRANSPORT
    type(realization_type), pointer :: realization
    type(stepper_type), pointer :: stepper
#endif    
#ifdef TRANSPORT
    type(tr_realization_type), pointer :: tr_realization
    type(tr_stepper_type), pointer :: tr_stepper
#endif

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

  implicit none
  
  type(simulation_type), pointer :: SimulationCreate
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
#ifndef TRANSPORT
  simulation%realization => RealizationCreate()
  simulation%stepper => TimestepperCreate()
#endif
#ifdef TRANSPORT
  simulation%tr_realization => TrRealizationCreate()
  simulation%tr_stepper => TrTimestepperCreate()
#endif
  
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

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return

#ifndef TRANSPORT
  call RealizationDestroy(simulation%realization)
  call TimestepperDestroy(simulation%stepper)
#endif
#ifdef TRANSPORT  
  call TrRealizationDestroy(simulation%tr_realization)
  call TrTimestepperDestroy(simulation%tr_stepper)
#endif  
  
end subroutine SimulationDestroy
  
end module Simulation_module
