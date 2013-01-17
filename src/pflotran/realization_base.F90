module Realization_Base_class

  use Level_module
  use Patch_module

  use Discretization_module
  use Option_module
  use Input_module
  use Debug_module
  use Output_Aux_module
  use Field_module
  use Reaction_Aux_module

  implicit none

  private

#include "definitions.h"
  type, public :: realization_base_type

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    type(level_list_type), pointer :: level_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(input_type), pointer :: input
    type(field_type), pointer :: field
    type(flow_debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    
    type(reaction_type), pointer :: reaction
    
  end type realization_base_type
  
  public :: RealizationBaseInit, &
            RealizationGetDataset, &
            RealizGetDatasetValueAtCell, &
            RealizationSetDataset

contains

! ************************************************************************** !
!
! RealizationBaseInit: Initializes variables/objects in base realization class
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine RealizationBaseInit(realization,option)

  implicit none
  
  class(realization_base_type) :: realization
  type(option_type), pointer :: option
  
  realization%id = 0
  if (associated(option)) then
    realization%option => option
  else
    realization%option => OptionCreate()
  endif
  nullify(realization%input)
  realization%discretization => DiscretizationCreate()
  realization%field => FieldCreate()
  realization%debug => DebugCreateFlow()
  realization%output_option => OutputOptionCreate()

  realization%level_list => LevelCreateList()

  nullify(realization%reaction)

  nullify(realization%patch)
  
end subroutine RealizationBaseInit

! ************************************************************************** !
!
! RealizationGetDataset: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationGetDataset(realization,vec,ivar,isubvar,isubvar1)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  call PatchGetDataset(realization%patch,realization%field, &
                       realization%reaction,realization%option, &
                       realization%output_option,vec,ivar,isubvar,isubvar1)

end subroutine RealizationGetDataset

! ************************************************************************** !
!
! RealizGetDatasetValueAtCell: Extracts variables indexed by ivar and isubvar
!                              from a realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id, &
                                     isubvar1)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetDatasetValueAtCell
  class(realization_base_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: ghosted_id
  
  PetscReal :: value
  
  value = PatchGetDatasetValueAtCell(realization%patch,realization%field, &
                                     realization%reaction, &
                                     realization%option, &
                                     realization%output_option, &
                                     ivar,isubvar,ghosted_id,isubvar1)
  RealizGetDatasetValueAtCell = value

end function RealizGetDatasetValueAtCell

! ************************************************************************** !
!
! RealizationSetDataset: Sets variables indexed by ivar and isubvar in a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationSetDataset(realization,vec,vec_format,ivar,isubvar)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetDataset(realization%patch,realization%field, &
                       realization%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetDataset

end module Realization_Base_class
