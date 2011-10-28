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
subroutine DatasetLoad(dataset,option,cur_time)

  use Option_module
  use HDF5_aux_module

  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  PetscReal :: cur_time
  
  call HDF5ReadDataset(dataset,option,cur_time)
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
  PetscInt :: i, j, k, t
  PetscInt :: nx, ny, nz, nt, nxXny, nyXnz, nxXnyXnz
  PetscInt :: count, index, idim
  
  ! Not necessary for 1D arrays
  if (dataset%ndims < 2) return
  
  select case(dataset%ndims)
    case(TWO_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      if (associated(dataset%time_array)) then
        ny = dataset%buffer_size
      else
        ny = dataset%dims(TWO_INTEGER)
      endif
      count = 0
      allocate(temp_real(nx*ny))
      do i = 1, nx
        do j = 0, ny-1
          index = j*nx+i
          count = count+1
          temp_real(index) = dataset%rarray(count)
        enddo
      enddo  
    case(THREE_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      if (associated(dataset%time_array)) then
        nz = dataset%buffer_size
      else
        nz = dataset%dims(THREE_INTEGER)
      endif
      nyXnz = ny*nz
      count = 0
      allocate(temp_real(nx*ny*nz))
      do i = 1, nx
        do j = 0, ny-1
          do k = 0, nz-1
            index = k*nxXny+j*nx+i
            count = count+1
            temp_real(index) = dataset%rarray(count)
          enddo
        enddo
      enddo  
    case(FOUR_INTEGER)
      nx = dataset%dims(ONE_INTEGER)
      ny = dataset%dims(TWO_INTEGER)
      nz = dataset%dims(THREE_INTEGER)
      if (associated(dataset%time_array)) then
        nt = dataset%buffer_size
      else
        nt = dataset%dims(FOUR_INTEGER)
      endif
      nxXny = nx*ny
      nxXnyXnz = nxXny*nz
      count = 0
      allocate(temp_real(nx*ny*nz*nt))
      do i = 1, nx
        do j = 0, ny-1
          do k = 0, nz-1
            do t = 0, nt-1
              index = t*nxXnyXnz+k*nxXny+j*nx+i
              count = count+1
              temp_real(index) = dataset%rarray(count)
            enddo
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
