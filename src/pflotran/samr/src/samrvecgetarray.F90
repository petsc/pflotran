subroutine pflotran_vecgetarrayf90(p_samr_patch, vec, f90ptr, ierr)
 use cf90interface_module
 implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  PetscFortranAddr, intent(inout):: p_samr_patch
  Vec:: vec
  PetscScalar, pointer :: f90ptr(:)
  integer :: ierr

  type(f90ptrwrap), pointer :: ptr
  PetscFortranAddr :: cptr

  if(p_samr_patch .eq. 0) then
     call VecGetArrayF90(vec, f90ptr, ierr)
  else
     ierr=0
     allocate(ptr)
     nullify(ptr%f90ptr)
     call assign_c_ptr(cptr, ptr)
     call samr_vecgetarrayf90(p_samr_patch, vec, cptr)
     f90ptr => ptr%f90ptr
     deallocate(ptr)
  endif

end subroutine pflotran_vecgetarrayf90
  
subroutine cf90bridge(array, len, f90wrap)

  use cf90interface_module
  implicit none
#include "include/finclude/petsc.h"

  PetscInt                        :: len
  PetscScalar, intent(in), target :: array(len)
  type(f90ptrwrap), intent(inout) :: f90wrap
  
  f90wrap%f90ptr => array

end subroutine cf90bridge

