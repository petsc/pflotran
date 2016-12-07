program test
#include "finclude/petscksp.h"
  use petscksp
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: x
  Vec :: b
  Mat :: A
  KSP :: ksp
  PC :: pc
  PetscReal :: value
  PetscReal :: tolerance
  PetscInt :: i, n
  PetscInt :: offset
  PetscReal, pointer :: x_ptr(:)
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  n = 2/size
  call MatCreateAIJ(PETSC_COMM_WORLD,n,n, &
                    PETSC_DETERMINE,PETSC_DETERMINE, &
                    1,PETSC_NULL_INTEGER, &
                    0,PETSC_NULL_INTEGER,A,ierr)
  call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DETERMINE,x,ierr)
  call VecDuplicate(x,b,ierr)
  offset = rank * n
  value = 1.d-16
  do i = 1, n
    call MatSetValue(A,i-1+offset,i-1+offset,value,INSERT_VALUES,ierr)
    call VecSetValue(b,i-1+offset,value,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetErrorIfNotConverged(ksp,PETSC_TRUE,ierr)
  call KSPSetOperators(ksp,A,A,ierr)
  call KSPGetPC(ksp,pc,ierr)
  call KSPSetFromOptions(ksp,ierr)
  call PCSetFromOptions(pc,ierr)
!  call KSPSetup(ksp,ierr)

  call PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS,ierr)
  tolerance = 1.d-20
  call PCFactorSetZeroPivot(pc,tolerance,ierr)

!  call KSPSetup(ksp,ierr)
  call KSPSolve(ksp,b,x,ierr)

  call VecGetArrayF90(x,x_ptr,ierr)
  print *, 'These values should be ~1: ', x_ptr 
  call VecRestoreArrayF90(x,x_ptr,ierr)

  call KSPDestroy(ksp,ierr)
  call MatDestroy(A,ierr)
  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test
