module Reaction_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "definitions.h"

  type, abstract, public :: reaction_sandbox_base_type
    class(reaction_sandbox_base_type), pointer :: next
  contains
    procedure(Base_Init), public, deferred, nopass :: Init 
    procedure(Base_Read), public, deferred, nopass :: ReadInput
    procedure(Base_React), public, deferred, nopass :: Evaluate
    procedure(Base_Destroy), public, deferred, nopass :: Destroy
  end type reaction_sandbox_base_type
  
  abstract interface
  
    subroutine Base_Init()

      implicit none
  
    end subroutine Base_Init 

    subroutine Base_Read(input,option)
    
      use Option_module
      use Input_module
  
      implicit none
  
      type(input_type) :: input
      type(option_type) :: option
  
    end subroutine Base_Read 
    
    subroutine Base_React(Residual,Jacobian,compute_derivative,rt_auxvar, &
                          global_auxvar,porosity,volume,reaction,option)
      use Option_module
      use Reaction_Aux_module
      use Reactive_Transport_Aux_module
      use Global_Aux_module
  
      implicit none
  
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
      
    end subroutine
    
    subroutine Base_Destroy()

      implicit none
  
    end subroutine Base_Destroy   
    
  end interface

end module Reaction_Sandbox_Base_class
