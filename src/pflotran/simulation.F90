module Simulation_module

  use Realization_class
  use Timestepper_module
  use Solver_module
  use Regression_module

#ifdef SURFACE_FLOW
  use Surface_Realization_class
#endif
  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type

    type(realization_type), pointer :: realization
    type(stepper_type), pointer :: flow_stepper
    type(stepper_type), pointer :: tran_stepper
#ifdef SURFACE_FLOW
    type(stepper_type), pointer :: surf_flow_stepper
#endif
#ifdef SURFACE_FLOW
    type(surface_realization_type), pointer :: surf_realization
#endif
    type(regression_type), pointer :: regression
  end type simulation_type
  
  interface SimulationCreate
    module procedure SimulationCreate1
    module procedure SimulationCreate2
  end interface
  
  public :: SimulationCreate, &
            SimulationDestroy
  
contains

! ************************************************************************** !
!
! SimulationCreate1: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate1()

  use Option_module
  
  implicit none
  
  type(simulation_type), pointer :: SimulationCreate1
  
  type(simulation_type), pointer :: simulation
  type(option_type), pointer :: option
  
  nullify(option)
  SimulationCreate1 => SimulationCreate2(option)
  
end function SimulationCreate1

! ************************************************************************** !
!
! SimulationCreate2: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate2(option)

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: SimulationCreate2
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
  simulation%realization => RealizationCreate(option)
  simulation%flow_stepper => TimestepperCreate()
  simulation%tran_stepper => TimestepperCreate()
#ifdef SURFACE_FLOW
  simulation%surf_flow_stepper => TimestepperCreate()
  simulation%surf_realization => SurfRealizCreate(option)
#endif
  nullify(simulation%regression)
  
  SimulationCreate2 => simulation
  
end function SimulationCreate2

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
  use General_module, only : GeneralDestroy

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return

  if (simulation%realization%option%nflowdof > 0) then
    select case(simulation%realization%option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsDestroy(simulation%realization)
      case(G_MODE)
        call GeneralDestroy(simulation%realization)
    end select
  endif

  if (simulation%realization%option%ntrandof > 0) then
    call RTDestroy(simulation%realization)
  endif

  call RealizationDestroy(simulation%realization)
  call TimestepperDestroy(simulation%flow_stepper)
  call TimestepperDestroy(simulation%tran_stepper)
#ifdef SURFACE_FLOW
  call TimestepperDestroy(simulation%surf_flow_stepper)
#endif

#ifdef SURFACE_FLOW
  call SurfRealizDestroy(simulation%surf_realization)
#endif

  call RegressionDestroy(simulation%regression)

  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationDestroy
  
end module Simulation_module
