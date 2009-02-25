subroutine create_samrai_vec(p_application, dof, centering, use_ghost, use_components, vec)
implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

PetscFortranAddr :: p_application
PetscInt :: dof
PetscInt :: centering
PetscTruth :: use_ghost
PetscTruth :: use_components
Vec :: vec
end subroutine create_samrai_vec

PetscInt function samr_patch_at_bc(p_patch, axis, dim)
implicit none

#include "finclude/petsc.h"

PetscFortranAddr :: p_patch
PetscInt :: axis,dim
samr_patch_at_bc=-1
end function samr_patch_at_bc

subroutine samr_patch_get_corners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
implicit none

#include "finclude/petsc.h"

PetscFortranAddr :: p_patch
PetscInt :: nxs, nys, nzs, nlx, nly, nlz

end subroutine samr_patch_get_corners

subroutine samr_patch_get_origin(p_patch, xs, ys, zs)
implicit none

#include "finclude/petsc.h"

PetscFortranAddr, intent(inout) :: p_patch
PetscReal, intent(inout) :: xs
PetscReal, intent(inout) :: ys
PetscReal, intent(inout) :: zs

end subroutine samr_patch_get_origin

subroutine samr_patch_get_ghostcorners(p_patch, nxs, nys, nzs, nlx, nly, nlz)
implicit none

#include "finclude/petsc.h"

PetscFortranAddr :: p_patch
PetscInt :: nxs, nys, nzs, nlx, nly, nlz

end subroutine samr_patch_get_ghostcorners

PetscInt function hierarchy_number_levels(p_application)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  hierarchy_number_levels=-1
end function hierarchy_number_levels

PetscInt function level_number_patches(p_application, ln)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscInt, intent(in) :: ln

  level_number_patches=-1
end function level_number_patches

logical function is_local_patch(p_application, ln, pn)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscInt, intent(in) :: ln
  PetscInt, intent(in) :: pn

  is_local_patch = .TRUE.

end function is_local_patch

PetscFortranAddr function hierarchy_get_patch(p_application, ln, pn)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscInt, intent(in) :: ln
  PetscInt, intent(in) :: pn

  hierarchy_get_patch = 0
  
end function hierarchy_get_patch

subroutine samr_physical_dimensions(p_application, nx, ny, nz)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscInt, intent(inout) :: nx
  PetscInt, intent(inout) :: ny
  PetscInt, intent(inout) :: nz
end subroutine samr_physical_dimensions

subroutine samr_get_origin(p_application, x0, y0, z0)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscReal, intent(inout) :: x0
  PetscReal, intent(inout) :: y0
  PetscReal, intent(inout) :: z0
  
end subroutine samr_get_origin

subroutine samr_get_upper_corner(p_application, x0, y0, z0)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr, intent(inout) :: p_application
  PetscReal, intent(inout) :: x0
  PetscReal, intent(inout) :: y0
  PetscReal, intent(inout) :: z0
  
end subroutine samr_get_upper_corner

subroutine assign_c_array_ptr(return_arg, pointer_arg)
 use cf90interface_module
implicit none

#include "finclude/petsc.h"
  type(f90ptrwrap), pointer :: pointer_arg
  PetscFortranAddr :: return_arg

end subroutine assign_c_array_ptr

subroutine samr_vecgetarraycellf90(patch, petscvec, f90wrap)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr, intent(inout):: patch
  Vec:: petscvec
  PetscFortranAddr :: f90wrap
end subroutine samr_vecgetarraycellf90

subroutine samr_vecgetarraysidef90(patch, axis, petscvec, f90wrap)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr, intent(inout):: patch
  PetscInt :: axis
  Vec:: petscvec
  PetscFortranAddr :: f90wrap
end subroutine samr_vecgetarraysidef90

subroutine samr_patch_get_spacing(p_samr_patch, dx, dy, dz)
implicit none
#include "finclude/petsc.h"
  PetscFortranAddr :: p_samr_patch
  PetscReal :: dx, dy, dz

end subroutine samr_patch_get_spacing

subroutine SAMRCreateMatrix(p_application, ndof, stencilsize, flowortransport, p_matrix)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
  PetscFortranAddr :: p_application
  PetscInt :: ndof
  PetscInt :: stencilsize
  PetscInt :: flowortransport
  Mat :: p_matrix

end subroutine SAMRCreateMatrix

subroutine SAMRGlobalToLocal(p_application, gvec, lvec, ierr)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr :: p_application
  Vec :: lvec
  Vec :: gvec
  PetscInt :: ierr

end subroutine SAMRGlobalToLocal

subroutine SAMRLocalToLocal(p_application, gvec, lvec, ierr)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr :: p_application
  Vec :: lvec
  Vec :: gvec
  PetscInt :: ierr

end subroutine SAMRLocalToLocal

subroutine SAMRCoarsenFaceFluxes(p_application, vec, ierr)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr :: p_application
  Vec :: vec
  PetscInt :: ierr

end subroutine SAMRCoarsenFaceFluxes

subroutine SAMRSetCurrentJacobianPatch(mat,patch) 
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

Mat :: mat
PetscFortranAddr :: patch
end subroutine SAMRSetCurrentJacobianPatch

subroutine SAMRSetJacobianSourceOnPatch(which_pc, index, val, application, patch) 
#include "finclude/petsc.h"
PetscInt :: which_pc
PetscInt :: index
PetscReal :: val
PetscFortranAddr :: application
PetscFortranAddr :: patch
end subroutine SAMRSetJacobianSourceOnPatch

subroutine samrpetscobjectstateincrease(vec)
implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  Vec :: vec
end subroutine samrpetscobjectstateincrease


subroutine samr_mpi_min(x,y,z)
  PetscScalar, intent(inout) :: x,y,z
end subroutine samr_mpi_min

subroutine samr_mpi_max(x,y,z)
  PetscScalar, intent(inout) :: x,y,z
end subroutine samr_mpi_max

subroutine SAMRBarrier()
end subroutine SAMRBarrier 

subroutine SAMRCopyVecToVecComponent(vec,svec, component)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  Vec :: vec, svec
  PetscInt :: component
end subroutine SAMRCopyVecToVecComponent

subroutine SAMRRegisterForViz(ptr,vec,component,dname,dnamec)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  PetscFortranAddr :: ptr
  Vec :: vec
  PetscInt :: component
  PetscInt :: dname, dnamec
end subroutine SAMRRegisterForViz

subroutine SAMRWritePlotData(ptr, time)
#include "finclude/petsc.h"
  PetscFortranAddr :: ptr
  PetscReal :: time
end subroutine SAMRWritePlotData

subroutine SAMRInitializePreconditioner(p_application, which_pc, pc)
#include "finclude/petsc.h"
#include "finclude/petscpc.h"
  PC :: pc
  PetscFortranAddr :: p_application
  PetscInt :: which_pc
end subroutine SAMRInitializePreconditioner
