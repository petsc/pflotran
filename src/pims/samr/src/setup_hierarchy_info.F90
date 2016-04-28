subroutine f_setup_hierarchy_data(pgrid, pintegrator)

 use pflow_gridtype_module
 use pflow_grid_module

 implicit none
 
#include "include/finclude/petsc.h"
#include "../../definitions.h"

 external pflow_setup_index
 type(pflowGrid), pointer :: pgrid
 type(time_stepping_context), intent(inout) :: pintegrator

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
 
  nlevels =  hierarchy_number_levels(pgrid%p_samr_hierarchy)
 
!  do ln=0,nlevels-1
!     npatches = level_number_patches(pgrid%p_samr_hierarchy, ln )
!     do pn=0,npatches-1
!        if(associated(pgrid%patchlevel_info(ln+1)%patches(pn+1)%patch_ptr) .eq. .TRUE. ) then
!           call pflow_setup_index(pgrid, pgrid%patchlevel_info(ln+1)%patches(pn+1)%patch_ptr)
!        endif
!     end do
!  end do

!  allocate(pgrid%ibndtyp(MAXBCREGIONS))
!  allocate(pgrid%iface(MAXBCREGIONS))

  do ln=0,nlevels-1
     npatches = level_number_patches(pgrid%p_samr_hierarchy, ln )
     do pn=0,npatches-1
        if(associated(pgrid%patchlevel_info(ln+1)%patches(pn+1)%patch_ptr) .eq. .TRUE. ) then
           call pflowGrid_setup(pgrid, pintegrator, pgrid%patchlevel_info(ln+1)%patches(pn+1)%patch_ptr, "pflow.in")
        endif
     end do
  end do
end subroutine f_setup_hierarchy_data
