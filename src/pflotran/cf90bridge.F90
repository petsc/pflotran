module cf90interface_module
implicit none
private
#include "include/finclude/petsc.h"

type, public :: f90ptrwrap
   PetscReal, pointer, dimension(:) :: f90ptr
end type f90ptrwrap

end module cf90interface_module
