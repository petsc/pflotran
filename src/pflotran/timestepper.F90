module Timestepper_module
 
  use Solver_module
 
  implicit none

  private
 
  type, public :: stepper_type
    type(solver_type), pointer :: solver
  end type stepper_type
  
end module Timestepper_module
