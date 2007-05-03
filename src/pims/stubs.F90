! this is a stub function so that we don't have to always
! link in SAMRAI
subroutine create_samrai_vec(p_hierarchy, dof, use_ghost, vec)
implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

PetscFortranAddr :: p_hierarchy
integer :: dof
PetscTruth :: use_ghost
Vec :: vec
end subroutine create_samrai_vec
