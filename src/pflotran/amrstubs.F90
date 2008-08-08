


! this is a stub function so that we don't have to always
! link with SAMRAI

integer function samr_patch_at_bc(p_patch, axis, dim)

#include "include/finclude/petsc.h"

PetscFortranAddr :: p_patch
integer :: axis,dim
samr_patch_at_bc=-1
end function samr_patch_at_bc


subroutine create_samrai_vec(p_application, dof, use_ghost, vec)
implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

PetscFortranAddr :: p_application
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

integer function hierarchy_number_levels(p_application)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  hierarchy_number_levels=-1
end function hierarchy_number_levels

integer function level_number_patches(p_application, ln)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  integer, intent(in) :: ln

  level_number_patches=-1
end function level_number_patches

logical function is_local_patch(p_application, ln, pn)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  integer, intent(in) :: ln
  integer, intent(in) :: pn

  is_local_patch = .TRUE.

end function is_local_patch

PetscFortranAddr function hierarchy_get_patch(p_application, ln, pn)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  integer, intent(in) :: ln
  integer, intent(in) :: pn

  hierarchy_get_patch = 0
  
end function hierarchy_get_patch

subroutine samr_physical_dimensions(p_application, nx, ny, nz)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  integer, intent(inout) :: nx
  integer, intent(inout) :: ny
  integer, intent(inout) :: nz
end subroutine samr_physical_dimensions

subroutine samr_get_origin(p_application, x0, y0, z0)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscReal, intent(inout) :: x0
  PetscReal, intent(inout) :: y0
  PetscReal, intent(inout) :: z0
  
end subroutine samr_get_origin

subroutine assign_c_ptr(return_arg, pointer_arg)
 use cf90interface_module
implicit none
#include "include/finclude/petsc.h"
  type(f90ptrwrap), pointer :: pointer_arg
  PetscFortranAddr :: return_arg

end subroutine assign_c_ptr

subroutine samr_vecgetarrayf90(patch, petscvec, f90wrap)
implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  PetscFortranAddr, intent(inout):: patch
  Vec:: petscvec
  PetscFortranAddr :: f90wrap
end subroutine samr_vecgetarrayf90

subroutine samr_patch_get_spacing(p_samr_patch, dx, dy, dz)
implicit none
#include "include/finclude/petsc.h"
  PetscFortranAddr :: p_samr_patch
  PetscReal :: dx, dy, dz

end subroutine samr_patch_get_spacing

subroutine SAMRCreateMatrix(p_application, ndof, stencilsize, flowortransport, p_matrix)
implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
  PetscFortranAddr :: p_application
  PetscInt :: ndof
  PetscInt :: stencilsize
  PetscInt :: flowortransport
  Mat :: p_matrix

end subroutine SAMRCreateMatrix

subroutine SAMRGlobalToLocal(p_application, gvec, lvec, ierr)
implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  PetscFortranAddr :: p_application
  Vec :: lvec
  Vec :: gvec
  PetscInt :: ndof
  PetscInt :: ierr

end subroutine SAMRGlobalToLocal

subroutine SAMRLocalToLocal(p_application, gvec, lvec, ierr)
implicit none
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  PetscFortranAddr :: p_application
  Vec :: lvec
  Vec :: gvec
  PetscInt :: ndof
  PetscInt :: ierr

end subroutine SAMRLocalToLocal

subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "include/finclude/petsc.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"

Mat :: mat
PetscFortranAddr :: patch
end subroutine SAMRSetCurrentJacobianPatch
