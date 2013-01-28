module Reaction_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: reaction_sandbox_base_type
    class(reaction_sandbox_base_type), pointer :: next
  contains
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
