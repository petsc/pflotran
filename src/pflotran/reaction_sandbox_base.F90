module Reaction_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: reaction_sandbox_base_type
    type(reaction_sandbox_base_type), pointer :: next
  contains
    procedure, public :: Init => Base_Init
    procedure, public :: ReadInput => Base_Read
    procedure, public :: Evaluate => Base_React
    procedure, public :: Destroy => Base_Destroy
  end type reaction_sandbox_base_type

contains

! ************************************************************************** !
!
! Base_Init: Initializes reaction sandbox at beginning of simulation
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine Base_Init(reaction_sandbox)

  implicit none
  
  class(reaction_sandbox_base_type) :: reaction_sandbox
      
end subroutine Base_Init

! ************************************************************************** !
!
! Base_Read: Reads input deck for reaction sandbox parameters
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine Base_Read(reaction_sandbox,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_base_type) :: reaction_sandbox
  type(input_type) :: input
  type(option_type) :: option

end subroutine Base_Read

! ************************************************************************** !
!
! Base_React: Evaluates reaction storing residual and/or Jacobian
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine Base_React(reaction_sandbox, &
                      Residual,Jacobian,compute_derivative,rt_auxvar, &
                      global_auxvar,porosity,volume,reaction,option)

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  type(reaction_sandbox_base_type) :: reaction_sandbox  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  PetscReal :: porosity
  PetscReal :: volume
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar

end subroutine Base_React

! ************************************************************************** !
!
! Base_Destroy: Destroys allocatable or pointer objects created in this 
!                 module
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine Base_Destroy(reaction_sandbox)

  implicit none
  
  class(reaction_sandbox_base_type) :: reaction_sandbox  

end subroutine Base_Destroy

end module Reaction_Sandbox_Base_class
