module Realization_Base_class

  use Patch_module

  use Discretization_module
  use Option_module
  use Input_Aux_module
  use Debug_module
  use Output_Aux_module
  use Field_module
  use Reaction_Aux_module
  use Data_Mediator_Base_class
  use Communicator_Base_module
  use Waypoint_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
  type, public :: realization_base_type

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    class(communicator_type), pointer :: comm1
    type(patch_list_type), pointer :: patch_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(input_type), pointer :: input
    type(field_type), pointer :: field
    type(debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    type(checkpoint_option_type), pointer :: checkpoint_option
    class(data_mediator_base_type), pointer :: flow_data_mediator_list
    class(data_mediator_base_type), pointer :: tran_data_mediator_list
    
    type(reaction_type), pointer :: reaction
    type(waypoint_list_type), pointer :: waypoint_list
    
  end type realization_base_type
  
  public :: RealizationBaseInit, &
            RealizationGetVariable, &
            RealizGetVariableValueAtCell, &
            RealizationSetVariable, &
            RealizCreateTranMassTransferVec, &
            RealizCreateFlowMassTransferVec, &
            RealizCreateSyncWaypointList, &
            RealizationBaseStrip

contains

! ************************************************************************** !

subroutine RealizationBaseInit(realization_base,option)
  ! 
  ! Initializes variables/objects in base realization class
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

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
  nullify(realization_base%comm1)  
  realization_base%field => FieldCreate()
  realization_base%debug => DebugCreate()
  realization_base%output_option => OutputOptionCreate()

  realization_base%patch_list => PatchCreateList()

  nullify(realization_base%reaction)
  nullify(realization_base%waypoint_list)
  nullify(realization_base%checkpoint_option)

  nullify(realization_base%patch)
  nullify(realization_base%flow_data_mediator_list)
  nullify(realization_base%tran_data_mediator_list)

end subroutine RealizationBaseInit

! ************************************************************************** !

subroutine RealizationGetVariable(realization_base,vec,ivar,isubvar,isubvar1)
  ! 
  ! Extracts variables indexed by ivar and isubvar from a
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

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

function RealizGetVariableValueAtCell(realization_base,ivar,isubvar,ghosted_id, &
                                     isubvar1)
  ! 
  ! Extracts variables indexed by ivar and isubvar
  ! from a realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

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

subroutine RealizationSetVariable(realization_base,vec,vec_format,ivar,isubvar)
  ! 
  ! Sets variables indexed by ivar and isubvar in a
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetVariable(realization_base%patch,realization_base%field, &
                       realization_base%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetVariable

! ************************************************************************** !

subroutine RealizCreateFlowMassTransferVec(this)
  ! 
  ! Creates the Vec where mass transfer is summed prior to being added to
  ! the reactive transport residual.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/15
  ! 
  implicit none
  
  class(realization_base_type) :: this
  
  PetscInt :: ierr
  
  if (this%field%flow_mass_transfer == 0) then
    call VecDuplicate(this%field%flow_xx,this%field%flow_mass_transfer, &
                      ierr);CHKERRQ(ierr)
  endif

end subroutine RealizCreateFlowMassTransferVec

! ************************************************************************** !

subroutine RealizCreateTranMassTransferVec(this)
  ! 
  ! Creates the Vec where mass transfer is summed prior to being added to
  ! the reactive transport residual.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/15
  ! 
  implicit none
  
  class(realization_base_type) :: this
  
  PetscInt :: ierr
  
  if (this%field%tran_mass_transfer == 0) then
    call VecDuplicate(this%field%tran_xx,this%field%tran_mass_transfer, &
                      ierr);CHKERRQ(ierr)
  endif

end subroutine RealizCreateTranMassTransferVec

! ************************************************************************** !

function RealizCreateSyncWaypointList(realization)
  !
  ! Creates a list of waypoints for outer synchronization of simulation process
  ! model couplers
  !
  ! Author: Glenn Hammond
  ! Date: 10/08/14
  !

  use Option_module
  use Waypoint_module
  use Time_Storage_module

  implicit none

  class(realization_base_type) :: realization

  type(waypoint_list_type), pointer :: RealizCreateSyncWaypointList

  type(waypoint_list_type), pointer :: new_waypoint_list
  type(waypoint_type), pointer :: cur_waypoint
  type(waypoint_type), pointer :: new_waypoint

  new_waypoint_list => WaypointListCreate()

  cur_waypoint => realization%waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%sync .or. cur_waypoint%final) then
      new_waypoint => WaypointCreate(cur_waypoint)
      call WaypointInsertInList(new_waypoint,new_waypoint_list)
      if (cur_waypoint%final) exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  RealizCreateSyncWaypointList => new_waypoint_list

end function RealizCreateSyncWaypointList

! ************************************************************************** !

subroutine RealizationBaseStrip(this)
  ! 
  ! Deallocates members of base realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  ! 
  use Data_Mediator_module
  
  implicit none
  
  class(realization_base_type) :: this
  
  call FieldDestroy(this%field)

!  call OptionDestroy(realization%option) !geh it will be destroy externally
  call OutputOptionDestroy(this%output_option)
  call CheckpointOptionDestroy(this%checkpoint_option)
  
  call DiscretizationDestroy(this%discretization)
  
  if (associated(this%comm1)) then
    call this%comm1%Destroy()
    deallocate(this%comm1)
  endif
  nullify(this%comm1)
  
  call PatchDestroyList(this%patch_list)
  nullify(this%patch)

  call DebugDestroy(this%debug)
  
  call DataMediatorDestroy(this%flow_data_mediator_list)
  call DataMediatorDestroy(this%tran_data_mediator_list)

  call WaypointListDestroy(this%waypoint_list)

end subroutine RealizationBaseStrip

end module Realization_Base_class
