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
  
  public :: createSimulation, destroySimulation
  
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
  
end function createSimulation  

! ************************************************************************** !
!
! destroySimulation: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine destroySimulation(simulation)

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return
    
  call destroySolution(simulation%solution)
  call destroyTimestepper(simulation%stepper)
  
end subroutine destroySimulation
  
end module Simulation_module
