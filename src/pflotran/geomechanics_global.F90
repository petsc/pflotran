#ifdef GEOMECH

module Geomechanics_Global_module

  use Geomechanics_Global_Aux_module

  implicit none

  private
  
#include "definitions.h"

  public GeomechGlobalSetup
!         GeomechGlobalSetAuxVarScalar, &
!         GeomechGlobalSetAuxVarVecLoc, &
!         GeomechGlobalUpdateAuxVars

contains

! ************************************************************************** !
!
! GeomechGlobalSetup: Set up global aux vars in a realization
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetup(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  
  implicit none

  type(geomech_realization_type)           :: geomech_realization
  
  ! There is only one patch in each realization
  call GeomechGlobalSetupPatch(geomech_realization)
  
end subroutine GeomechGlobalSetup

! ************************************************************************** !
!
! GeomechGlobalSetupPatch: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetupPatch(geomech_realization)

  use Geomechanics_Realization_module
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Coupler_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
 
  implicit none
  
  type(geomech_realization_type)            :: geomech_realization

  type(option_type), pointer                :: option
  type(geomech_patch_type),pointer          :: patch
  type(geomech_grid_type), pointer          :: grid
  type(geomech_coupler_type), pointer       :: boundary_condition
  type(geomech_coupler_type), pointer       :: source_sink

  PetscInt :: ghosted_id
  type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  PetscInt                                  :: ivertex
  
  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  patch%geomech_aux%GeomechGlobal => GeomechGlobalAuxCreate()
  
  allocate(aux_vars(grid%ngmax_node))
  do ghosted_id = 1, grid%ngmax_node
    call GeomechGlobalAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%geomech_aux%GeomechGlobal%aux_vars => aux_vars
  patch%geomech_aux%GeomechGlobal%num_aux = grid%ngmax_node
   
end subroutine GeomechGlobalSetupPatch

#if 0
! ************************************************************************** !
!
! GeomechGlobalSetAuxVarScalar: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetAuxVarScalar(geomech_realization,value,ivar)

  use Geomechanics_Realization_module
  use Level_module
  use Patch_module

  implicit none

  type(geomech_realization_type) :: geomech_realization
  PetscReal :: value
  PetscInt :: ivar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => geomech_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      geomech_realization%patch => cur_patch
      call GeomechGlobalSetAuxVarScalarPatch(geomech_realization,value,ivar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeomechGlobalSetAuxVarScalar

! ************************************************************************** !
!
! GeomechGlobalSetAuxVarScalarPatch: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetAuxVarScalarPatch(geomech_realization,value,ivar)

  use Geomechanics_Realization_module
  use Option_module
  use Patch_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

  type(geomech_realization_type) :: geomech_realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => geomech_realization%patch
  option => geomech_realization%option  
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%head = value
      enddo
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux_bc
        patch%geomech_aux%GeomechGlobal%aux_vars_bc(i)%head = value
      enddo
    case(SURFACE_LIQUID_TEMPERATURE)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%temp = value
      enddo
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux_bc
        patch%geomech_aux%GeomechGlobal%aux_vars_bc(i)%temp = value
      enddo
    case(SURFACE_LIQUID_DENSITY)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%den_kg(option%liquid_phase) = value
      enddo
      do i=1, geomech_realization%patch%geomech_aux%GeomechGlobal%num_aux_bc
        patch%geomech_aux%GeomechGlobal%aux_vars_bc(i)%den_kg(option%liquid_phase) = value
      enddo
  end select
  
end subroutine GeomechGlobalSetAuxVarScalarPatch

! ************************************************************************** !
!
! GeomechGlobalSetAuxVarVecLoc: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetAuxVarVecLoc(geomech_realization,vec_loc,ivar,isubvar)

  use Geomechanics_Realization_module
  use Level_module
  use Patch_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  

  type(geomech_realization_type) :: geomech_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => geomech_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      geomech_realization%patch => cur_patch
      call GeomechGlobalSetAuxVarVecLocPatch(geomech_realization,vec_loc,ivar,isubvar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GeomechGlobalSetAuxVarVecLoc

! ************************************************************************** !
!
! GeomechGlobalSetAuxVarVecLocPatch: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalSetAuxVarVecLocPatch(geomech_realization,vec_loc,ivar,isubvar)

  use Geomechanics_Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_realization_type) :: geomech_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => geomech_realization%patch
  grid => patch%grid
  option => geomech_realization%option
  
  call GridVecGetArrayF90(grid,vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%geomech_aux%GeomechGlobal%aux_vars(ghosted_id)%head(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_TEMPERATURE)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%geomech_aux%GeomechGlobal%aux_vars(ghosted_id)%temp(1) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_DENSITY)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%geomech_aux%GeomechGlobal%aux_vars(ghosted_id)%den_kg(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
  end select

  call GridVecRestoreArrayF90(grid,vec_loc,vec_loc_p,ierr)

end subroutine GeomechGlobalSetAuxVarVecLocPatch


! ************************************************************************** !
!
! GeomechGlobalUpdateAuxVars: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/17/13
! 
! ************************************************************************** !
subroutine GeomechGlobalUpdateAuxVars(geomech_realization,time_level)

  use Geomechanics_Realization_module
  use Geomechanics_Field_module
  use Option_module
  use Discretization_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               LIQUID_DENSITY
  
  type(geomech_realization_type) :: geomech_realization
  PetscInt :: time_level
  
  type(geomech_field_type), pointer :: geomech_field
  type(option_type), pointer :: option
  
  option => geomech_realization%option
  geomech_field => geomech_realization%geomech_field
  
  ! liquid density
  call geomechRealizGetDataset(geomech_realization,geomech_field%work,LIQUID_DENSITY, &
                             ZERO_INTEGER)
  call DiscretizationGlobalToLocal(geomech_realization%discretization, &
                                   geomech_field%work,geomech_field%work_loc,ONEDOF)
  call GeomechGlobalSetAuxVarVecLoc(geomech_realization,geomech_field%work_loc, &
                                    LIQUID_DENSITY,time_level)

  select case(option%iflowmode)
    case(TH_MODE)
      ! head
      call geomechRealizGetDataset(geomech_realization,geomech_field%work, &
              SURFACE_LIQUID_HEAD,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(geomech_realization%discretization, &
                                  geomech_field%work,geomech_field%work_loc,ONEDOF)
      call GeomechGlobalSetAuxVarVecLoc(geomech_realization,geomech_field%work_loc, &
              SURFACE_LIQUID_HEAD,time_level)
 
      ! temperature
      call geomechRealizGetDataset(geomech_realization,geomech_field%work, &
              SURFACE_LIQUID_TEMPERATURE, ZERO_INTEGER)
      call DiscretizationGlobalToLocal(geomech_realization%discretization, &
                                   geomech_field%work,geomech_field%work_loc,ONEDOF)
      call GeomechGlobalSetAuxVarVecLoc(geomech_realization,geomech_field%work_loc, &
              SURFACE_LIQUID_TEMPERATURE,time_level)
      

  end select

end subroutine GeomechGlobalUpdateAuxVars
#endif

end module Geomechanics_Global_module

#endif
