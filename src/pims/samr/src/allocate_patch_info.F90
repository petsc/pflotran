subroutine allocate_patch_info(p_samr_hierarchy, patchlevel_info)
  use pflow_gridtype_module
  implicit none

#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_hierarchy
!  type(PatchLevelInfoPtr), dimension(:), pointer, intent(inout) :: patchlevel_info
  type(PatchLevelInfoPtr), dimension(:), pointer :: patchlevel_info

  interface
     integer function hierarchy_number_levels(p_hierarchy)
     PetscFortranAddr, intent(inout) :: p_hierarchy
   end function hierarchy_number_levels

   integer function level_number_patches(p_hierarchy, ln)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
   end function level_number_patches

   logical function is_local_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function is_local_patch

   PetscFortranAddr function hierarchy_get_patch(p_hierarchy, ln, pn)
     PetscFortranAddr, intent(inout) :: p_hierarchy
     integer, intent(in) :: ln
     integer, intent(in) :: pn
   end function hierarchy_get_patch
   
  end interface

  integer :: nlevels
  integer :: npatches
  integer :: ln
  integer :: pn
  logical :: islocal

  nlevels =  hierarchy_number_levels(p_samr_hierarchy)
  
  allocate(patchlevel_info(nlevels))

  do ln=0,nlevels-1
     npatches = level_number_patches(p_samr_hierarchy, ln )
     allocate(patchlevel_info(ln+1)%patches(npatches))
     do pn=0,npatches-1
        islocal = is_local_patch(p_samr_hierarchy, ln, pn);
        if(islocal) then
           allocate(patchlevel_info(ln+1)%patches(pn+1)%patch_ptr)
           patchlevel_info(ln+1)%patches(pn+1)%patch_ptr%p_samr_patch = hierarchy_get_patch(p_samr_hierarchy, ln, pn)
        else              
           nullify(patchlevel_info(ln+1)%patches(pn+1)%patch_ptr)
        endif
     end do
  end do
  
end subroutine allocate_patch_info
