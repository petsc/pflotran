module Realization_Base_class

  use Patch_module

  use Discretization_module
  use Option_module
  use Input_Aux_module
  use Debug_module
  use Output_Aux_module
  use Field_module
  use Reaction_Aux_module
  use Mass_Transfer_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
  type, public :: realization_base_type

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    type(patch_list_type), pointer :: patch_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(input_type), pointer :: input
    type(field_type), pointer :: field
    type(debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    type(mass_transfer_type), pointer :: flow_mass_transfer_list
    type(mass_transfer_type), pointer :: rt_mass_transfer_list
    
    type(reaction_type), pointer :: reaction
    
  end type realization_base_type
  
  public :: RealizationBaseInit, &
            RealizationGetVariable, &
            RealizGetVariableValueAtCell, &
            RealizationSetVariable

contains

! ************************************************************************** !
!
! RealizationBaseInit: Initializes variables/objects in base realization class
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine RealizationBaseInit(realization_base,option)

  implicit none
  
  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  
  realization_base%id = 0
  if (associated(option)) then
    realization_base%option => option
  else
    realization_base%option => OptionCreate()
  endif
  nullify(realization_base%input)
  realization_base%discretization => DiscretizationCreate()
  realization_base%field => FieldCreate()
  realization_base%debug => DebugCreate()
  realization_base%output_option => OutputOptionCreate()

  realization_base%patch_list => PatchCreateList()

  nullify(realization_base%reaction)

  nullify(realization_base%patch)
  nullify(realization_base%flow_mass_transfer_list)
  nullify(realization_base%rt_mass_transfer_list)

end subroutine RealizationBaseInit

! ************************************************************************** !
!
! RealizationGetVariable: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationGetVariable(realization_base,vec,ivar,isubvar,isubvar1)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  call PatchGetVariable(realization_base%patch,realization_base%field, &
                       realization_base%reaction,realization_base%option, &
                       realization_base%output_option,vec,ivar,isubvar,isubvar1)

end subroutine RealizationGetVariable

! ************************************************************************** !
!
! RealizGetVariableValueAtCell: Extracts variables indexed by ivar and isubvar
!                              from a realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function RealizGetVariableValueAtCell(realization_base,ivar,isubvar,ghosted_id, &
                                     isubvar1)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetVariableValueAtCell
  class(realization_base_type) :: realization_base
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: ghosted_id
  
  PetscReal :: value
  
  value = PatchGetVariableValueAtCell(realization_base%patch,realization_base%field, &
                                     realization_base%reaction, &
                                     realization_base%option, &
                                     realization_base%output_option, &
                                     ivar,isubvar,ghosted_id,isubvar1)
  RealizGetVariableValueAtCell = value

end function RealizGetVariableValueAtCell

! ************************************************************************** !
!
! RealizationSetVariable: Sets variables indexed by ivar and isubvar in a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationSetVariable(realization_base,vec,vec_format,ivar,isubvar)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetVariable(realization_base%patch,realization_base%field, &
                       realization_base%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetVariable

end module Realization_Base_class
