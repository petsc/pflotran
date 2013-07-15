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
  use Mass_Transfer_module

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
    type(mass_transfer_type), pointer :: flow_mass_transfer_list
    type(mass_transfer_type), pointer :: rt_mass_transfer_list
    
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
  realization_base%debug => DebugCreateFlow()
  realization_base%output_option => OutputOptionCreate()

  realization_base%level_list => LevelCreateList()

  nullify(realization_base%reaction)

  nullify(realization_base%patch)
  nullify(realization_base%flow_mass_transfer_list)
  nullify(realization_base%rt_mass_transfer_list)

end subroutine RealizationBaseInit

! ************************************************************************** !
!
! RealizationGetDataset: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationGetDataset(realization_base,vec,ivar,isubvar,isubvar1)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  call PatchGetDataset(realization_base%patch,realization_base%field, &
                       realization_base%reaction,realization_base%option, &
                       realization_base%output_option,vec,ivar,isubvar,isubvar1)

end subroutine RealizationGetDataset

! ************************************************************************** !
!
! RealizGetDatasetValueAtCell: Extracts variables indexed by ivar and isubvar
!                              from a realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function RealizGetDatasetValueAtCell(realization_base,ivar,isubvar,ghosted_id, &
                                     isubvar1)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetDatasetValueAtCell
  class(realization_base_type) :: realization_base
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: ghosted_id
  
  PetscReal :: value
  
  value = PatchGetDatasetValueAtCell(realization_base%patch,realization_base%field, &
                                     realization_base%reaction, &
                                     realization_base%option, &
                                     realization_base%output_option, &
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
subroutine RealizationSetDataset(realization_base,vec,vec_format,ivar,isubvar)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetDataset(realization_base%patch,realization_base%field, &
                       realization_base%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetDataset

end module Realization_Base_class
