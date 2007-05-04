subroutine allocate_patch_info(p_samr_hierarchy, patchlevel_info)
  use pflow_gridtype_module
  implicit none

#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_hierarchy
  type(PatchLevelInfoPtr), dimension(:), pointer, intent(inout) :: patchlevel_info

  interface
     integer function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   integer function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
   end function level_number_patches
  end interface

  integer :: nlevels
  integer :: npatches
  integer :: ln

  nlevels =  hierarchy_number_levels(p_samr_hierarchy)
  
  allocate(patchlevel_info(nlevels))

  do ln=0,nlevels-1
     npatches = level_number_patches(p_samr_hierarchy, ln )
     allocate(patchlevel_info(ln+1)%patches(npatches))
  end do
  
end subroutine allocate_patch_info
