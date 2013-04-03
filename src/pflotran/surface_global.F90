#ifdef SURFACE_FLOW

module Surface_Global_module

  use Surface_Global_Aux_module

  implicit none

  private
  
#include "definitions.h"

  public SurfaceGlobalSetup, &
         SurfaceGlobalSetAuxVarScalar, &
         SurfaceGlobalSetAuxVarVecLoc, &
         SurfaceGlobalUpdateAuxVars

contains

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalSetup(surf_realization)

  use Surface_Realization_class
  use Level_module
  use Patch_module
  
  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceGlobalSetupPatch(surf_realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceGlobalSetup

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalSetupPatch(surf_realization)

  use Surface_Realization_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  type(surface_realization_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection
  type(surface_global_auxvar_type), pointer :: aux_vars(:)
  type(surface_global_auxvar_type), pointer :: aux_vars_bc(:)
  type(surface_global_auxvar_type), pointer :: aux_vars_ss(:)
  
  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid

  patch%surf_aux%SurfaceGlobal => SurfaceGlobalAuxCreate()
  
  ! allocate aux_var data structures for all grid cells  
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SurfaceGlobalAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%surf_aux%SurfaceGlobal%aux_vars => aux_vars
  patch%surf_aux%SurfaceGlobal%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! aux_var data structures for them  
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(aux_vars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceGlobalAuxVarInit(aux_vars_bc(iconn),option)
    enddo
    patch%surf_aux%SurfaceGlobal%aux_vars_bc => aux_vars_bc
  endif
  patch%surf_aux%SurfaceGlobal%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! aux_var data structures for them  
  source_sink => patch%source_sinks%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    sum_connection = sum_connection + &
                     source_sink%connection_set%num_connections
    source_sink => source_sink%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(aux_vars_ss(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceGlobalAuxVarInit(aux_vars_ss(iconn),option)
    enddo
    patch%surf_aux%SurfaceGlobal%aux_vars_ss => aux_vars_ss
  endif
  patch%surf_aux%SurfaceGlobal%num_aux_ss = sum_connection

  option%iflag = 0
  
end subroutine SurfaceGlobalSetupPatch

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalSetAuxVarScalar(surf_realization,value,ivar)

  use Surface_Realization_class
  use Level_module
  use Patch_module

  implicit none

  type(surface_realization_type) :: surf_realization
  PetscReal :: value
  PetscInt :: ivar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceGlobalSetAuxVarScalarPatch(surf_realization,value,ivar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceGlobalSetAuxVarScalar

! ************************************************************************** !
! ************************************************************************** !
subroutine SurfaceGlobalSetAuxVarScalarPatch(surf_realization,value,ivar)

  use Surface_Realization_class
  use Option_module
  use Patch_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

  type(surface_realization_type) :: surf_realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => surf_realization%patch
  option => surf_realization%option  
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%aux_vars(i)%head = value
      enddo
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%aux_vars_bc(i)%head = value
      enddo
    case(SURFACE_LIQUID_TEMPERATURE)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%aux_vars(i)%temp = value
      enddo
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%aux_vars_bc(i)%temp = value
      enddo
    case(SURFACE_LIQUID_DENSITY)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%aux_vars(i)%den_kg(option%liquid_phase) = value
      enddo
      do i=1, surf_realization%patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%aux_vars_bc(i)%den_kg(option%liquid_phase) = value
      enddo
  end select
  
end subroutine SurfaceGlobalSetAuxVarScalarPatch

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalSetAuxVarVecLoc(surf_realization,vec_loc,ivar,isubvar)

  use Surface_Realization_class
  use Level_module
  use Patch_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  

  type(surface_realization_type) :: surf_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      surf_realization%patch => cur_patch
      call SurfaceGlobalSetAuxVarVecLocPatch(surf_realization,vec_loc,ivar,isubvar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfaceGlobalSetAuxVarVecLoc

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalSetAuxVarVecLocPatch(surf_realization,vec_loc,ivar,isubvar)

  use Surface_Realization_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(surface_realization_type) :: surf_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  
  call GridVecGetArrayF90(grid,vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%aux_vars(ghosted_id)%head(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_TEMPERATURE)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%aux_vars(ghosted_id)%temp(1) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_DENSITY)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%aux_vars(ghosted_id)%den_kg(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
  end select

  call GridVecRestoreArrayF90(grid,vec_loc,vec_loc_p,ierr)

end subroutine SurfaceGlobalSetAuxVarVecLocPatch


! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceGlobalUpdateAuxVars(surf_realization,time_level)

  use Surface_Realization_class
  use Surface_Field_module
  use Option_module
  use Discretization_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               LIQUID_DENSITY
  
  type(surface_realization_type) :: surf_realization
  PetscInt :: time_level
  
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  
  ! liquid density
  call SurfRealizGetDataset(surf_realization,surf_field%work,LIQUID_DENSITY, &
                             ZERO_INTEGER)
  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%work,surf_field%work_loc,ONEDOF)
  call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
                                    LIQUID_DENSITY,time_level)

  select case(option%iflowmode)
    case(TH_MODE)
      ! head
      call SurfRealizGetDataset(surf_realization,surf_field%work, &
              SURFACE_LIQUID_HEAD,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                  surf_field%work,surf_field%work_loc,ONEDOF)
      call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
              SURFACE_LIQUID_HEAD,time_level)
 
      ! temperature
      call SurfRealizGetDataset(surf_realization,surf_field%work, &
              SURFACE_LIQUID_TEMPERATURE, ZERO_INTEGER)
      call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%work,surf_field%work_loc,ONEDOF)
      call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
              SURFACE_LIQUID_TEMPERATURE,time_level)
      

  end select

end subroutine SurfaceGlobalUpdateAuxVars

end module Surface_Global_module

#endif
