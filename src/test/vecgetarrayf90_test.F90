program test

  implicit none
  
#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: vec
  
  PetscReal, pointer :: vec_ptr(:)
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  call VecCreate(PETSC_COMM_SELF,vec,ierr)
  if (rank == 0) then
    call VecSetSizes(vec,16,PETSC_DECIDE,ierr)  
  else
    call VecSetSizes(vec,0,PETSC_DECIDE,ierr)  
  endif
  call VecSetFromOptions(vec,ierr)

  call VecGetArrayF90(vec,vec_ptr,ierr)
  if (rank == 0) then
    print *, vec_ptr(1)
  endif
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

  call VecDestroy(vec, ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test
