module Reaction_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: reaction_sandbox_base_type
    type(reaction_sandbox_base_type), pointer :: next
  contains
    procedure, public :: Init => RSandboxInit
    procedure, public :: Destroy => RSandboxDestroy
  end type reaction_sandbox_base_type

contains

! ************************************************************************** !
!
! RSandboxInit: Initializes reaction sandbox at beginning of simulation
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxInit(reaction_sandbox)

  implicit none
  
  class(reaction_sandbox_base_type) :: reaction_sandbox
      
  nullify(reaction_sandbox%next)

end subroutine RSandboxInit

! ************************************************************************** !
!
! RSandboxDestroy: Destroys allocatable or pointer objects created in this 
!                  module
! author: Glenn Hammond
! date: 11/08/12
!
! ************************************************************************** !
subroutine RSandboxDestroy(reaction_sandbox)

  implicit none

  class(reaction_sandbox_base_type) :: reaction_sandbox
  
end subroutine RSandboxDestroy

end module Reaction_Sandbox_Base_class
