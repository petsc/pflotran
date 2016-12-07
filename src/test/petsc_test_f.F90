program test
#include "finclude/petscmat.h"
  use petscmat
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  Vec :: vec
  PetscViewer :: viewer
  
  character(len=32) :: filename
  
  PetscInt :: ndof
  PetscInt :: my_ndof
  
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  ndof = 100000
  my_ndof = int(ndof / size)
  if (mod(ndof,size) > rank) my_ndof = my_ndof + 1
  call VecCreateMPI(PETSC_COMM_WORLD,my_ndof, ndof, vec, ierr)
  call VecSet(vec, -999.d0, ierr)
  filename = 'vec.bin'
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &
                             viewer, ierr)  
!  call PetscViewerBinarySetFlowControl(viewer,2,ierr)
  if (rank == 0) print *, 'Before VecView'
  call VecView(vec, viewer, ierr)
  if (rank == 0) print *, 'After VecView'
  call PetscViewerDestroy(viewer, ierr)             
  call VecDestroy(vec, ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test
