subroutine pims_vecgetarrayf90(grid, patch, vec, f90ptr, ierr)
#include "include/finclude/petscvec.h"
 use petscvec
 use pflow_gridtype_module
 use cf90interface_module
 implicit none 

  type(pflowGrid), intent(inout) :: grid
  type(pflow_localpatch_info) :: patch
  Vec vec
  PetscScalar, pointer :: f90ptr(:)
  integer :: ierr

  type(f90ptrwrap), pointer :: ptr
  PetscFortranAddr :: cptr

  if(grid%Samrai_drive==PETSC_TRUE) then
     allocate(ptr)
     nullify(ptr%f90ptr)
     call pflow_loc(cptr, ptr)
     call samr_vecgetarrayf90(patch%p_samr_patch, vec, cptr)
     f90ptr => ptr%f90ptr
     deallocate(ptr)
  else
     call VecGetArrayF90(vec, f90ptr, ierr)
  endif

end subroutine pims_vecgetarrayf90

subroutine cf90bridge(array, len, f90wrap)
  use  cf90interface_module
  implicit none
#include "include/finclude/petsc.h"

  PetscInt                        :: len
  PetscScalar, intent(in), target :: array(len)
  type(f90ptrwrap), intent(inout) :: f90wrap
  
  f90wrap%f90ptr => array

end subroutine cf90bridge

