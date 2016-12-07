! this is a stub function so that we don't have to always
! link in SAMRAI




integer function samr_patch_at_bc(p_patch, axis, dim)

#include "include/finclude/petsc.h"

PetscFortranAddr :: p_patch
integer :: axis,dim

end function samr_patch_at_bc


subroutine create_samrai_vec(p_hierarchy, dof, use_ghost, vec)
#include "include/finclude/petscvec.h"
use petscvec
implicit none

PetscFortranAddr :: p_hierarchy
integer :: dof
PetscTruth :: use_ghost
Vec :: vec
end subroutine create_samrai_vec

subroutine samr_patch_get_corners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
implicit none

#include "include/finclude/petsc.h"

PetscFortranAddr :: p_patch
integer :: nxs, nys, nzs, nlx, nly, nlz

end subroutine samr_patch_get_corners

subroutine samr_patch_get_ghostcorners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
implicit none

#include "include/finclude/petsc.h"

PetscFortranAddr :: p_patch
integer :: nxs, nys, nzs, nlx, nly, nlz

end subroutine samr_patch_get_ghostcorners

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

subroutine pims_vecgetarrayf90(grid, patch, vec, f90ptr, ierr)
#include "include/finclude/petscvec.h"
  use petscvec
  use pflow_gridtype_module
 
  implicit none

  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info) :: patch
  Vec :: vec
  PetscScalar, dimension(:), pointer :: f90ptr
  integer :: ierr

  PetscScalar, dimension(:), pointer :: f90iptr

  nullify(f90iptr)

  call VecGetArrayF90(vec, f90iptr, ierr)

  f90ptr => f90iptr 

end subroutine pims_vecgetarrayf90
