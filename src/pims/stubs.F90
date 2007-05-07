! this is a stub function so that we don't have to always
! link in SAMRAI
subroutine create_samrai_vec(p_hierarchy, dof, use_ghost, vec)
implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

PetscFortranAddr :: p_hierarchy
integer :: dof
PetscTruth :: use_ghost
Vec :: vec
end subroutine create_samrai_vec

! this stub does not make use of any SAMR information
! it simply allocates space for one pflow_localpatch_info
! object
subroutine allocate_patch_info(p_samr_hierarchy, patchlevel_info)
  use pflow_gridtype_module
  implicit none

  PetscFortranAddr :: p_samr_hierarchy
  type(PatchLevelInfoPtr), dimension(:), pointer :: patchlevel_info
  allocate(patchlevel_info(1))
  allocate(patchlevel_info(1)%patches(1))
  allocate(patchlevel_info(1)%patches(1)%patch_ptr)
  
end subroutine allocate_patch_info
