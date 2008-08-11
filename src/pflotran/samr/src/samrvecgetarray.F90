subroutine cf90bridge(array, len, f90wrap)

  use cf90interface_module
  implicit none

#include "include/finclude/petsc.h"

  PetscInt                        :: len
  PetscReal, intent(in), target   :: array(len)
  type(f90ptrwrap), intent(inout) :: f90wrap
  
  f90wrap%f90ptr => array

end subroutine cf90bridge

