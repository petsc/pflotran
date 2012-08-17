module Solid_Solution_Aux_module
  
  use Mineral_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: solid_solution_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: num_stoich_solids
    PetscInt :: num_end_members
    type(stoichiometric_solid_type), pointer :: stoich_solid
    type(solid_solution_type), pointer :: next
  end type solid_solution_type

  type, public :: stoichiometric_solid_type
    type(mineral_type) :: stoich_solid
    type(mineral_type), pointer :: end_members
  end type stoichiometric_solid_type
    
  type, public :: solid_solution_rxn_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    type(solid_solution_type), pointer :: list
    type(mineral_rxn_type), pointer :: mineral_reaction
  end type solid_solution_rxn_type

  public :: SolidSolutionReactionCreate, &
            SolidSolutionReactionDestroy
             
contains

! ************************************************************************** !
!
! SolidSolutionReactionCreate: Allocate and initialize solid solution object
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
function SolidSolutionReactionCreate()

  use Option_module

  implicit none
  
  type(solid_solution_rxn_type), pointer :: SolidSolutionReactionCreate
  
  type(solid_solution_rxn_type), pointer :: solid_solution

  allocate(solid_solution)
  
  nullify(solid_solution%list)

  SolidSolutionReactionCreate => solid_solution
  
end function SolidSolutionReactionCreate

! ************************************************************************** !
!
! SolidSolutionReactionDestroy: Deallocates a solid solution object
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
subroutine SolidSolutionReactionDestroy(solid_solution)

  implicit none

  type(solid_solution_rxn_type), pointer :: solid_solution
  
  deallocate(solid_solution)
  nullify(solid_solution)

end subroutine SolidSolutionReactionDestroy

end module Solid_Solution_Aux_module
