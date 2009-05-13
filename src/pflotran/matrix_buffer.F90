module Matrix_Buffer_module
 
  implicit none

  private

#include "definitions.h"

  type, public :: matrix_buffer_type
    PetscInt :: itype
    PetscInt :: i_width
    PetscInt :: j_width
    PetscInt :: k_width
    PetscReal, pointer :: values(:,:)
  end type matrix_buffer_type
 
  public :: MatrixBufferCreate, MatrixBufferDestroy
 
contains

! ************************************************************************** !
!
! MatrixBufferCreate: Creates a matrix object
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
function MatrixBufferCreate()
  
  implicit none

  type(matrix_buffer_type), pointer :: MatrixBufferCreate
  
  type(matrix_buffer_type), pointer :: matrix_buffer
  
  allocate(matrix_buffer)
  matrix_buffer%itype = 0
  matrix_buffer%i_width = 0
  matrix_buffer%j_width = 0
  matrix_buffer%k_width = 0
  nullify(matrix_buffer%values)
  MatrixBufferCreate => matrix_buffer

end function MatrixBufferCreate

! ************************************************************************** !
!
! MatrixBufferInit: Initializes matrix buffer object
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
subroutine MatrixBufferInit(matrix_buffer,gnx,gny,gnz)

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  
  PetscInt :: gnx, gny, gnz

  matrix_buffer%i_width = 0
  matrix_buffer%j_width = gnx
  matrix_buffer%k_width = gnx*gny
  allocate(matrix_buffer%values(7,gnx*gny*gnz))
  matrix_buffer%values = 0.d0

end subroutine MatrixBufferInit

! ************************************************************************** !
!
! MatrixBufferZero: Zeros matrix buffer values
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
subroutine MatrixBufferZero(matrix_buffer)

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  

  matrix_buffer%values = 0.d0

end subroutine MatrixBufferZero

! ************************************************************************** !
!
! MatrixBufferAdd: Adds values to matrix buffer object
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
subroutine MatrixBufferAdd(matrix_buffer,irow,icol,value)

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  
  PetscInt :: irow, icol
  PetscReal :: value
  
  PetscInt :: index

  if (icol < 0) then
    index = 4
  else
    index = irow-icol
    if (index == 0) then
      index = 4
    else if (index == -matrix_buffer%i_width) then
      index = 3
    else if (index == matrix_buffer%i_width) then
      index = 5
    else if (index == -matrix_buffer%j_width) then
      index = 2
    else if (index == matrix_buffer%j_width) then
      index = 6
    else if (index == -matrix_buffer%k_width) then
      index = 1
    else if (index == matrix_buffer%k_width) then
      index = 7
    endif
  endif
  matrix_buffer%values(index,irow) = matrix_buffer%values(index,irow) + value

end subroutine MatrixBufferAdd

! ************************************************************************** !
!
! MatrixBufferSetValuesHypre: Sets values in PETSc Hypre matrix
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
subroutine MatrixBufferSetValuesHypre(A,matrix_buffer)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
  
  Mat :: A
  type(matrix_buffer_type), pointer :: matrix_buffer  

  PetscInt :: icol
  PetscErrorCode :: ierr

  do icol = 1, 7
    call MatSetValuesLocal(A,icol-1,0,0,0, &
                           matrix_buffer%values,INSERT_VALUES,ierr)
  enddo

end subroutine MatrixBufferSetValuesHypre

! ************************************************************************** !
!
! MatrixBufferDestroy: Destroys a matrix buffer object
! author: Glenn Hammond
! date: 05/13/09
!
! ************************************************************************** !
subroutine MatrixBufferDestroy(matrix_buffer)

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer
  
  if (.not.associated(matrix_buffer)) return
  
  if (associated(matrix_buffer%values)) deallocate(matrix_buffer%values)
  nullify(matrix_buffer%values)
  
  deallocate(matrix_buffer)
  nullify(matrix_buffer)
  
end subroutine MatrixBufferDestroy

end module Matrix_Buffer_module
