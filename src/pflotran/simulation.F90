module Simulation_module

  use Solution_module
  use Timestepper_module
  use Solver_module

  implicit none
  
  private

  type, public :: simulation_type

    type(solution_type), pointer :: solution
    type(stepper_type), pointer :: stepper
    type(solver_type), pointer :: solver

  end type simulation_type
  
  public :: createSimulation
  
contains

! ************************************************************************** !
!
! createSolution: Allocates and initializes a new Solution object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function createSimulation()

  implicit none
  
  type(simulation_type), pointer :: createSimulation
  
  allocate(createSimulation)
  createSimulation%solution => createSolution()
  nullify(createSimulation%stepper) ! nullify these since we are not sure 
  nullify(createSimulation%solver)  ! whether we need them yet
  
end function createSimulation  
  
end module Simulation_module
