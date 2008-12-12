module Global_module

  use Global_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"
  
  public GlobalSetup, GlobalSetAuxVarScalar, GlobalSetAuxVarVecLoc, &
         GlobalUpdateDenAndSat

contains

! ************************************************************************** !
!
! GlobalSetup: 
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine GlobalSetup(realization)

  use Realization_module
  use Level_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GlobalSetupPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GlobalSetup
  
! ************************************************************************** !
!
! GlobalSetupPatch: Creates arrays for auxilliary variables
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine GlobalSetupPatch(realization)

  use Realization_module
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: ghosted_id, iconn, sum_connection
  type(global_auxvar_type), pointer :: aux_vars(:), aux_vars_bc(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Global => GlobalAuxCreate()
  
  ! allocate aux_var data structures for all grid cells  
  allocate(aux_vars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call GlobalAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%aux%Global%aux_vars => aux_vars
  patch%aux%Global%num_aux = grid%ngmax
  
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
  allocate(aux_vars_bc(sum_connection))
  do iconn = 1, sum_connection
    call GlobalAuxVarInit(aux_vars_bc(iconn),option)
  enddo
  patch%aux%Global%aux_vars_bc => aux_vars_bc
  patch%aux%Global%num_aux_bc = sum_connection
  
end subroutine GlobalSetupPatch

! ************************************************************************** !
!
! GlobalSetAuxVarScalar: Sets values of auxvar data using a scalar value. 
! author: Glenn Hammond
! date: 11/19/08
!
! ************************************************************************** !
subroutine GlobalSetAuxVarScalar(realization,value,ivar)

  use Realization_module
  use Level_module
  use Patch_module

  implicit none

  type(realization_type) :: realization
  PetscReal :: value
  PetscInt :: ivar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GlobalSetAuxVarScalarPatch(realization,value,ivar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GlobalSetAuxVarScalar

! ************************************************************************** !
!
! GlobalSetAuxVarScalarPatch: Updates the auxilliary variables associated with 
!                             the Global problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine GlobalSetAuxVarScalarPatch(realization,value,ivar)

  use Realization_module
  use Option_module
  use Patch_module
  
  implicit none

  type(realization_type) :: realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => realization%patch
  option => realization%option  
  
  select case(ivar)
    case(PRESSURE)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%aux_vars(i)%pres = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%aux_vars_bc(i)%pres = value
      enddo
    case(TEMPERATURE)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%aux_vars(i)%temp = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%aux_vars_bc(i)%temp = value
      enddo
    case(LIQUID_DENSITY)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%aux_vars(i)%den_kg(option%liquid_phase) = value
        patch%aux%Global%aux_vars(i)%den(option%liquid_phase) = value/ &
                                                                FMWH2O
      enddo
      do i=1, realization%patch%aux%Global%num_aux_bc
        patch%aux%Global%aux_vars_bc(i)%den_kg(option%liquid_phase) = value
        patch%aux%Global%aux_vars_bc(i)%den(option%liquid_phase) = value/ &
                                                                   FMWH2O
      enddo
    case(LIQUID_SATURATION)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%aux_vars(i)%sat(option%liquid_phase) = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%aux_vars_bc(i)%sat(option%liquid_phase) = value
      enddo
  end select
  
end subroutine GlobalSetAuxVarScalarPatch

! ************************************************************************** !
!
! GlobalSetAuxVarVecLoc: Sets values of auxvar data using a vector. 
! author: Glenn Hammond
! date: 11/19/08
!
! ************************************************************************** !
subroutine GlobalSetAuxVarVecLoc(realization,vec_loc,ivar,isubvar)

  use Realization_module
  use Level_module
  use Patch_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"  

  type(realization_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GlobalSetAuxVarVecLocPatch(realization,vec_loc,ivar,isubvar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine GlobalSetAuxVarVecLoc

! ************************************************************************** !
!
! GlobalSetAuxVarVecPatch: Updates the auxilliary variables associated with 
!                             the Global problem
! author: Glenn Hammond
! date: 12/10/07
!
! ************************************************************************** !
subroutine GlobalSetAuxVarVecLocPatch(realization,vec_loc,ivar,isubvar)

  use Realization_module
  use Patch_module
  use Grid_module
  use Option_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  call GridVecGetArrayF90(grid,vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(PRESSURE)
      do ghosted_id=1, grid%ngmax
        patch%aux%Global%aux_vars(ghosted_id)%pres = vec_loc_p(ghosted_id)
      enddo
    case(TEMPERATURE)
      do ghosted_id=1, grid%ngmax
        patch%aux%Global%aux_vars(ghosted_id)%temp = vec_loc_p(ghosted_id)
      enddo
    case(LIQUID_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%den_kg_store(option%liquid_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%den_kg_store(option%liquid_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%den_kg(option%liquid_phase) = vec_loc_p(ghosted_id)
            patch%aux%Global%aux_vars(ghosted_id)%den(option%liquid_phase) = &
              vec_loc_p(ghosted_id)/FMWH2O
          enddo
        end select
    case(LIQUID_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%sat_store(option%liquid_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%sat_store(option%liquid_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%aux_vars(ghosted_id)%sat(option%liquid_phase) = &
              vec_loc_p(ghosted_id)
          enddo
        end select
  end select

  call GridVecRestoreArrayF90(grid,vec_loc,vec_loc_p,ierr)

end subroutine GlobalSetAuxVarVecLocPatch

! ************************************************************************** !
!
! GlobalUpdateDenAndSat: Updates the densities and saturations in auxilliary 
!                    variables associated with reactive transport
! author: Glenn Hammond
! date: 11/03/08
!
! ************************************************************************** !
subroutine GlobalUpdateDenAndSat(realization,weight)

  use Realization_module
  use Level_module
  use Patch_module

  type(realization_type) :: realization
  PetscReal :: weight
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call GlobalUpdateDenAndSatPatch(realization,weight)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine GlobalUpdateDenAndSat

! ************************************************************************** !
!
! GlobalUpdateDenAndSatPatch: Updates the densities and saturations in auxilliary 
!                         variables associated with reactive transport
! author: Glenn Hammond
! date: 11/03/08
!
! ************************************************************************** !
subroutine GlobalUpdateDenAndSatPatch(realization,weight)

  use Realization_module
  use Patch_module
  use Option_module
  
  implicit none

  type(realization_type) :: realization
  PetscReal :: weight
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id
  
  option => realization%option
  patch => realization%patch  
  
  do ghosted_id = 1, patch%aux%Global%num_aux
    ! interpolate density and saturation based on weight
    patch%aux%Global%aux_vars(ghosted_id)%den_kg(:) = &
      (weight*patch%aux%Global%aux_vars(ghosted_id)%den_kg_store(:,TIME_TpDT)+ &
       (1.d0-weight)*patch%aux%Global%aux_vars(ghosted_id)%den_kg_store(:,TIME_T))
    patch%aux%Global%aux_vars(ghosted_id)%sat(:) = &
      (weight*patch%aux%Global%aux_vars(ghosted_id)%sat_store(:,TIME_TpDT)+ &
       (1.d0-weight)*patch%aux%Global%aux_vars(ghosted_id)%sat_store(:,TIME_T))
  enddo     
  
end subroutine GlobalUpdateDenAndSatPatch

end module Global_module
