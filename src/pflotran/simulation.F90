module Simulation_module

  use Solution_module
  use Timestepper_module
  use Solver_module

  implicit none
  
  private

  type, public :: simulation_type

    type(solution_type), pointer :: solution
    type(stepper_type), pointer :: stepper

  end type simulation_type
  
  public :: SimulationCreate, SimulationDestroy
  
contains

! ************************************************************************** !
!
! SolutionCreate: Allocates and initializes a new Solution object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate()

  implicit none
  
  type(simulation_type), pointer :: SimulationCreate
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
  simulation%solution => SolutionCreate()
  simulation%stepper => TimestepperCreate()
  
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
    
  call SolutionDestroy(simulation%solution)
  call TimestepperDestroy(simulation%stepper)
  
end subroutine SimulationDestroy
  
end module Simulation_module
