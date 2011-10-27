module Dataset_module
 
  use Dataset_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  public :: DatasetLoad

contains

! ************************************************************************** !
!
! DatasetLoad: Loads a dataset from file
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetLoad(dataset,option)

  use Option_module
  use HDF5_aux_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  call HDF5ReadDataset(dataset,option)
  call DatasetReorder(dataset,option)
  ! calculate min/max values
  if (associated(dataset%rarray)) then
    dataset%rmax = maxval(dataset%rarray)
    dataset%rmin = minval(dataset%rarray)
  endif

end subroutine DatasetLoad

! ************************************************************************** !
!
! DatasetReorder: If a dataset is loaded from an HDF5 file, and it was
!                 multidimensional in the HDF5 file, the array needs to be
!                 reordered fro Fortran -> C indexing.  This subroutine
!                 takes care of the reordering.
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetReorder(dataset,option)

  use Option_module
  
  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  PetscInt, allocatable :: dims(:)
  PetscReal, allocatable :: temp_real(:)
  PetscInt :: i, j, k
  PetscInt :: nx, ny, nz, nxXny, nyXnz
  PetscInt :: count, index, idim
  
  ! Not necessary for 1D arrays
  if (dataset%ndims < 2) return
  
  allocate(dims(dataset%ndims))
  dims = dataset%dims
  do idim = 1, dataset%ndims
    dataset%dims(idim) = dims(dataset%ndims-idim+1)
  enddo
  deallocate(dims)
  
  select case(dataset%ndims)
    case(TWO_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      count = 0
      allocate(temp_real(nx*ny))
      do j = 1, ny
        do i = 1, nz
          index = j+(i-1)*ny
          count = count+1
          temp_real(index) = dataset%rarray(count)
        enddo
      enddo  
    case(THREE_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      nz = dataset%dims(THREE_INTEGER)
      nyXnz = ny*nz
      count = 0
      allocate(temp_real(nx*ny*nz))
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            index = k+(j-1)*nz+(i-1)*nyXnz
            count = count+1
            temp_real(index) = dataset%rarray(count)
          enddo
        enddo
      enddo  
    case default
      write(option%io_buffer,*) dataset%ndims
      option%io_buffer = 'Dataset reordering not yet supported for rank ' // &
                         trim(adjustl(option%io_buffer)) // ' array.'
      call printErrMsg(option)
  end select

  dataset%rarray = temp_real
  deallocate(temp_real)
  
end subroutine DatasetReorder

end module Dataset_module
