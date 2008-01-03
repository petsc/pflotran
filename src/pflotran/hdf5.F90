module HDF5_module

  implicit none

#include "include/finclude/petsc.h"

#include "definitions.h"
  
  private
  
  PetscErrorCode :: ierr
  integer :: hdf5_err

#ifdef USE_HDF5
  
  public :: HDF5MapLocalToNaturalIndices, &
            HDF5ReadIntegerArray, &
            HDF5ReadRealArray, &
            HDF5WriteStructDataSetFromVec, &
            HDF5WriteStructuredDataSet
  
contains

! ************************************************************************** !
!
! HDF5MapLocalToNaturalIndices: Set up indices array that map local cells to 
!                               entries in HDF5 grid cell vectors
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5MapLocalToNaturalIndices(grid,option,file_id,indices)

  use hdf5
  
  use Option_module
  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer :: indices(:)
  
  character(len=MAXSTRINGLENGTH) :: string 
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: local_ghosted_id, local_id, natural_id
  integer :: index_count
  integer :: cell_count
  integer(HSIZE_T) :: num_cells_in_file
  integer ::temp_int, i
  
  integer, allocatable :: cell_ids(:)
  
  integer :: read_block_size = HDF5_READ_BUFFER_SIZE

  string = "Cell Id"
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_cells_in_file,hdf5_err)
  if (num_cells_in_file /= grid%nmax) then
    if (option%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_cells_in_file, ') does not match the dimensions of the ', &
               'domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
  
  allocate(cell_ids(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  cell_count = 0
  index_count = 0
  memory_space_id = -1
  do
    if (cell_count >= num_cells_in_file) exit
    temp_int = num_cells_in_file-cell_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = cell_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                               hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == 0) then                           
#endif
      call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,cell_ids,dims,hdf5_err, &
                     memory_space_id,file_space_id,prop_id)                     
#ifdef HDF5_BROADCAST
    endif
    if (option%commsize > 1) &
      call mpi_bcast(cell_ids,dims(1),MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
#endif     
    do i=1,dims(1)
      cell_count = cell_count + 1
      natural_id = cell_ids(i)
      local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (local_ghosted_id > 0) then
        local_id = grid%nG2L(local_ghosted_id)
        if (local_id > 0) then
          index_count = index_count + 1
          indices(index_count) = cell_count
        endif
      endif
    enddo
  enddo
  
  deallocate(cell_ids)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  if (index_count /= grid%nlmax) then
    if (option%myrank == 0) &
      print *, 'ERROR: Number of indices read (', index_count, ') does not ', &
               'match the number of local grid cells (', grid%nlmax, ').'
      call PetscFinalize(ierr)
      stop
  endif

end subroutine HDF5MapLocalToNaturalIndices

! ************************************************************************** !
!
! HDF5ReadRealArray: Read in local real values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5ReadRealArray(grid,option,file_id,num_indices,indices,string, &
                             real_array)

  use hdf5
  
  use Grid_module
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer :: num_indices
  integer :: indices(:)
  character(len=MAXSTRINGLENGTH) :: string 
  real*8 :: real_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: index_count
  integer :: real_count, prev_real_count
  integer(HSIZE_T) :: num_reals_in_file
  integer :: temp_int, i, index
  
  real*8, allocatable :: real_buffer(:)
  
  integer :: read_block_size = HDF5_READ_BUFFER_SIZE

  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_file,hdf5_err)
#if 0
  if (num_reals_in_file /= grid%nmax) then
    if (option%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_reals_in_file, ') does not match the dimensions of the ', &
               'domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  allocate(real_buffer(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  real_count = 0
  prev_real_count = 0
  index_count = 0
  memory_space_id = -1
  do i=1,num_indices
    index = indices(i)
    if (index > real_count) then
      do 
        if (index <= real_count) exit
        temp_int = num_reals_in_file-real_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
          dims(1) = temp_int
          call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
        endif
        ! offset is zero-based
        offset(1) = real_count
        length(1) = dims(1)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                   length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
        if (option%myrank == 0) then                           
#endif
          call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)
#ifdef HDF5_BROADCAST
        endif
        if (option%commsize > 1) &
          call mpi_bcast(real_buffer,dims(1),MPI_DOUBLE_PRECISION,0, &
                         PETSC_COMM_WORLD,ierr)
#endif
        prev_real_count = real_count
        real_count = real_count + length(1)                  
      enddo
    endif
    real_array(i) = real_buffer(index-prev_real_count)
  enddo

#ifdef HDF5_BROADCAST
  do 
    if (real_count >= num_reals_in_file) exit
    temp_int = num_reals_in_file-real_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = real_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    if (option%myrank == 0) then                           
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
    endif
    if (option%commsize > 1) &
      call mpi_bcast(real_buffer,dims(1),MPI_DOUBLE_PRECISION,0, &
                     PETSC_COMM_WORLD,ierr)
    real_count = real_count + length(1)                  
  enddo
#endif

  deallocate(real_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

end subroutine HDF5ReadRealArray

! ************************************************************************** !
!
! HDF5ReadIntegerArray: Read in local integer values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5ReadIntegerArray(grid,option,file_id,num_indices,indices, &
                                string,integer_array)

  use hdf5
  
  use Grid_module
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer :: num_indices
  integer :: indices(:)
  character(len=MAXSTRINGLENGTH) :: string 
  integer :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: index_count
  integer :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers_in_file
  integer :: temp_int, i, index
  
  integer, allocatable :: integer_buffer(:)
  
  integer :: read_block_size = HDF5_READ_BUFFER_SIZE

  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
#if 0
  if (num_integers_in_file /= grid%nmax) then
    if (option%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_integers_in_file, ') does not match the dimensions of ', &
               'the domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  
  allocate(integer_buffer(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  integer_count = 0
  prev_integer_count = 0
  index_count = 0
  memory_space_id = -1

  do i=1,num_indices
    index = indices(i)
    if (index > integer_count) then
      do
        if (index <= integer_count) exit
        temp_int = num_integers_in_file-integer_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
          dims(1) = temp_int
          call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
        endif
        ! offset is zero-based
        offset(1) = integer_count
        length(1) = dims(1)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                   length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
        if (option%myrank == 0) then                           
#endif
          call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)   
#ifdef HDF5_BROADCAST
        endif
        if (option%commsize > 1) &
          call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,0, &
                         PETSC_COMM_WORLD,ierr)
#endif
        prev_integer_count = integer_count
        integer_count = integer_count + length(1)                  
      enddo
    endif
    integer_array(i) = integer_buffer(index-prev_integer_count)
  enddo

#ifdef HDF5_BROADCAST
  do
    if (integer_count >= num_integers_in_file) exit
    temp_int = num_integers_in_file-integer_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = integer_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    if (option%myrank == 0) then                           
      call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)   
    endif
    if (option%commsize > 1) &
      call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
    integer_count = integer_count + length(1)                  
  enddo
#endif

  deallocate(integer_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

end subroutine HDF5ReadIntegerArray

! ************************************************************************** !
!
! HDF5WriteIntegerArray: Writes an integer array to an hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5WriteIntegerArray(grid,option,file_id,num_indices,indices, &
                                 string,integer_array)

  use hdf5
  
  use Grid_module
  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer :: num_indices
  integer :: indices(:)
  character(len=MAXSTRINGLENGTH) :: string 
  integer :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer :: rank
  integer :: index_count
  integer :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers_in_file
  integer :: temp_int, i, index
  
  integer, allocatable :: integer_buffer(:)
  
  integer :: read_block_size = HDF5_READ_BUFFER_SIZE

  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
#if 0
  if (num_integers_in_file /= grid%nmax) then
    if (option%myrank == 0) then
      print *, 'ERROR: ', trim(string), ' data space dimension (', &
               num_integers_in_file, ') does not match the dimensions of ', &
               'the domain (', grid%nmax, ').'
      call PetscFinalize(ierr)
      stop
    endif
  endif
#endif
  
  allocate(integer_buffer(read_block_size))
  
  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
  dims = 0
  integer_count = 0
  prev_integer_count = 0
  index_count = 0
  memory_space_id = -1

  do i=1,num_indices
    index = indices(i)
    if (index > integer_count) then
      do
        if (index <= integer_count) exit
        temp_int = num_integers_in_file-integer_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
          dims(1) = temp_int
          call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
        endif
        ! offset is zero-based
        offset(1) = integer_count
        length(1) = dims(1)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                   length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
        if (option%myrank == 0) then                           
#endif
          call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)   
#ifdef HDF5_BROADCAST
        endif
        if (option%commsize > 1) &
          call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,0, &
                         PETSC_COMM_WORLD,ierr)
#endif
        prev_integer_count = integer_count
        integer_count = integer_count + length(1)                  
      enddo
    endif
    integer_array(i) = integer_buffer(index-prev_integer_count)
  enddo

#ifdef HDF5_BROADCAST
  do
    if (integer_count >= num_integers_in_file) exit
    temp_int = num_integers_in_file-integer_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = integer_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    if (option%myrank == 0) then                           
      call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)   
    endif
    if (option%commsize > 1) &
      call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,0, &
                     PETSC_COMM_WORLD,ierr)
    integer_count = integer_count + length(1)                  
  enddo
#endif

  deallocate(integer_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

end subroutine HDF5WriteIntegerArray

! ************************************************************************** !
!
! HDF5WriteStructDataSetFromVec: Writes data from a PetscVec to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine HDF5WriteStructDataSetFromVec(name,realization,vec,file_id,data_type)

  use hdf5
  use Realization_module
  use Grid_module
  use Option_module
  
  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  character(len=32) :: name
  type(realization_type) :: realization
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscScalar, pointer :: vec_ptr(:)
  
  grid => realization%grid
  option => realization%option
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
!GEH - Structured Grid Dependence - Begin
  call HDF5WriteStructuredDataSet(name,vec_ptr,file_id,data_type, &
                                  grid%structured_grid%nx, &
                                  grid%structured_grid%ny, &
                                  grid%structured_grid%nz, &
                                  grid%structured_grid%nlx, &
                                  grid%structured_grid%nly, &
                                  grid%structured_grid%nlz, &
                                  grid%structured_grid%nxs, &
                                  grid%structured_grid%nys, &
                                  grid%structured_grid%nzs)
!GEH - Structured Grid Dependence - End
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine HDF5WriteStructDataSetFromVec

! ************************************************************************** !
!
! HDF5WriteStructuredDataSet: Writes data from an array into HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine HDF5WriteStructuredDataSet(name,array,file_id,data_type, &
                                      nx_global,ny_global,nz_global, &
                                      nx_local,ny_local,nz_local, &
                                      istart_local,jstart_local,kstart_local)

  use hdf5
  
  implicit none
  
  character(len=32) :: name
  real*8 :: array(:)
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer :: nx_local, ny_local, nz_local
  integer :: nx_global, ny_global, nz_global
  integer :: istart_local, jstart_local, kstart_local
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer :: rank
  
  integer, pointer :: int_array(:)
  real*8, pointer :: double_array(:)
  integer :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  integer :: ny_local_X_nz_local
  integer :: num_to_write

  ny_local_X_nz_local = ny_local*nz_local
  num_to_write = nx_local*ny_local_X_nz_local
  
  ! memory space which is a 1D vector  
  rank = 1
  dims = 0
  dims(1) = num_to_write
  if (num_to_write == 0) dims(1) = 1
  call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank = 3
#define INVERT
#ifndef INVERT
    dims(1) = nx_global
    dims(2) = ny_global
    dims(3) = nz_global
#else
! have to trick hdf5 for now with inverted ordering
    dims(3) = nx_global
    dims(2) = ny_global
    dims(1) = nz_global
#endif
  call h5screate_simple_f(rank,dims,file_space_id,hdf5_err,dims)


  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  call h5dcreate_f(file_id,name,data_type,file_space_id, &
                   data_set_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)
  
  ! create the hyperslab
#ifndef INVERT
  start(1) = istart_local
  start(2) = jstart_local
  start(3) = kstart_local
  length(1) =  nx_local
  length(2) =  ny_local
  length(3) =  nz_local
#else
  start(3) = istart_local
  start(2) = jstart_local
  start(1) = kstart_local
  length(3) =  nx_local
  length(2) =  ny_local
  length(1) =  nz_local
#endif
  if (num_to_write == 0) length(1) = 1
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err) 
  else
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif
  if (num_to_write > 0) then
    if (data_type == H5T_NATIVE_INTEGER) then
      allocate(int_array(nx_local*ny_local*nz_local))
#ifdef INVERT
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array(id) = int(array(count))
          enddo
        enddo
      enddo
#else
      do i=1,grid%nlmax
        int_array(i) = int(array(i))
      enddo
#endif
      call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      deallocate(int_array)
    else
#ifdef INVERT
      allocate(double_array(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            double_array(id) = array(count)
          enddo
        enddo
      enddo
      call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      deallocate(double_array)
#else
      call h5dwrite_f(data_set_id,data_type,array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
#endif
    endif
    call h5pclose_f(prop_id,hdf5_err)
  endif
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)

end subroutine HDF5WriteStructuredDataSet
!GEH - Structured Grid Dependence - End

#endif ! USE_HDF5

end module HDF5_module