subroutine vectompizero(vin, aout, ierr)

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec, intent(in) :: vin
  real*8 :: aout(1)
  integer, intent(out) :: ierr

  integer :: myrank, commsize
  integer :: mysize, mylow, myhigh
  integer, allocatable :: localsize(:), displacement(:)
!BEGIN_PETSCF90ARRAY
  real*8, pointer :: v_p(:)
!END_PETSCF90ARRAY

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, commsize, ierr)

  if(myrank == 0) then
    allocate(displacement(commsize))
    allocate(localsize(commsize))
  endif
  
  ! Find out the size of the local portion of the vector on each
  ! processor.
  call VecGetLocalSize(vin, mysize, ierr)
  call VecGetOwnershipRange(vin, mylow, myhigh, ierr)
  call MPI_Gather(mysize, 1, MPI_INTEGER, &
    localsize, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
  call MPI_Gather(mylow, 1, MPI_INTEGER, &
    displacement, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
  
  ! Now gather the local arrays onto processor zero.
  call VecGetArrayF90(vin, v_p, ierr)
  call MPI_Gatherv(v_p(1), mysize, MPI_DOUBLE_PRECISION, &
    aout, localsize, displacement, MPI_DOUBLE_PRECISION, &
    0, PETSC_COMM_WORLD, ierr)

  if(myrank == 0) then
    deallocate(displacement)
    deallocate(localsize)
  endif
  call VecRestoreArrayF90(vin, v_p, ierr)
end subroutine vectompizero

subroutine svectompizero(vin, aout, ierr)

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  Vec, intent(in) :: vin
  real*4 :: aout(1)
  integer, intent(out) :: ierr

  integer :: myrank, commsize
  integer :: mysize, mylow, myhigh
  integer, allocatable :: localsize(:), displacement(:)
  real*4, allocatable :: sbuff(:)
!BEGIN_PETSCF90ARRAY
  real*8, pointer :: v_p(:)
!END_PETSCF90ARRAY

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, commsize, ierr)

  if(myrank == 0) then
    allocate(displacement(commsize))
    allocate(localsize(commsize))
  endif
  
  ! Find out the size of the local portion of the vector on each
  ! processor.
  call VecGetLocalSize(vin, mysize, ierr)
  call VecGetOwnershipRange(vin, mylow, myhigh, ierr)
  call MPI_Gather(mysize, 1, MPI_INTEGER, &
    localsize, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
  call MPI_Gather(mylow, 1, MPI_INTEGER, &
    displacement, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
  
  ! Now gather the local arrays onto processor zero.
  ! Note that we need to copy what is in v_p to sbuff, since we 
  ! want to gather to a single-precision array.
  call VecGetArrayF90(vin, v_p, ierr)
  allocate(sbuff(mysize))
  do i = 1, mysize
    sbuff(i) = real(v_p(i),4)
  enddo
  call MPI_Gatherv(sbuff(1), mysize, MPI_REAL, &
    aout, localsize, displacement, MPI_REAL, &
    0, PETSC_COMM_WORLD, ierr)
  deallocate(sbuff)

  if(myrank == 0) then
    deallocate(displacement)
    deallocate(localsize)
  endif
  call VecRestoreArrayF90(vin, v_p, ierr)
end subroutine svectompizero
