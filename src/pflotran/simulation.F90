module Simulation_module

  use Realization_class
  use Timestepper_module
  use Solver_module
  use Regression_module

  use Surface_Realization_class

  use Geomechanics_Realization_class

  use PFLOTRAN_Constants_module


  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public :: simulation_type

    type(realization_type), pointer :: realization
    type(timestepper_type), pointer :: flow_timestepper
    type(timestepper_type), pointer :: tran_timestepper
    type(timestepper_type), pointer :: surf_flow_timestepper
    type(surface_realization_type), pointer :: surf_realization
    type(geomech_realization_type), pointer :: geomech_realization
    type(timestepper_type), pointer :: geomech_timestepper
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

function SimulationCreate1()
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Option_module
  
  implicit none
  
  type(simulation_type), pointer :: SimulationCreate1
  
  type(simulation_type), pointer :: simulation
  type(option_type), pointer :: option
  
  nullify(option)
  SimulationCreate1 => SimulationCreate2(option)
  
end function SimulationCreate1

! ************************************************************************** !

function SimulationCreate2(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: SimulationCreate2
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
  simulation%realization => RealizationCreate(option)
  simulation%flow_timestepper => TimestepperCreate()
  simulation%tran_timestepper => TimestepperCreate()
  simulation%surf_flow_timestepper => TimestepperCreate()
  simulation%surf_realization => SurfRealizCreate(option)
  simulation%geomech_realization => GeomechRealizCreate(option)
  simulation%geomech_timestepper => TimestepperCreate()
  nullify(simulation%regression)
  
  SimulationCreate2 => simulation
  
end function SimulationCreate2

! ************************************************************************** !

subroutine SimulationDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Richards_module, only : RichardsDestroy
  use Reactive_Transport_module, only : RTDestroy
  use General_module, only : GeneralDestroy

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return

  if (associated(simulation%realization)) then
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

    ! fmy: the following appears NOT deallocating 'timestepper' allocated in
    ! 'subsurface_factory.F90' Line 233. (TODO - but don't know why)
    call RealizationDestroyLegacy(simulation%realization)
  endif
  call TimestepperDestroy(simulation%flow_timestepper)
  call TimestepperDestroy(simulation%tran_timestepper)
  call TimestepperDestroy(simulation%surf_flow_timestepper)
  call SurfRealizDestroy(simulation%surf_realization)

  call TimestepperDestroy(simulation%geomech_timestepper)
  call GeomechRealizDestroy(simulation%geomech_realization)

  call RegressionDestroy(simulation%regression)

  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationDestroy
  
end module Simulation_module
