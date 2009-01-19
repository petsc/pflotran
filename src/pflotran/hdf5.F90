module HDF5_module

#ifdef USE_HDF5
  use hdf5
#endif
  use Logging_module

  implicit none

#include "definitions.h"

  private
  
  PetscErrorCode :: ierr

  PetscTruth, public :: trick_hdf5 = PETSC_FALSE

#ifdef USE_HDF5
  PetscMPIInt :: hdf5_err

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE  
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  public :: HDF5MapLocalToNaturalIndices, &
            HDF5ReadIntegerArray, &
            HDF5ReadRealArray, &
            HDF5WriteStructDataSetFromVec, &
            HDF5WriteStructuredDataSet, &
            HDF5ReadRegionFromFile, &
            HDF5ReadMaterialsFromFile

#else

  public :: HDF5ReadRegionFromFile, &
            HDF5ReadMaterialsFromFile, &
            HDF5ReadPermeabilitiesFromFile

#endif
  
contains

#ifdef USE_HDF5
! ************************************************************************** !
!
! HDF5MapLocalToNaturalIndices: Set up indices array that maps local cells to 
!                               entries in HDF5 grid cell vectors
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5MapLocalToNaturalIndices(grid,option,file_id, &
                                        dataset_name,dataset_size, &
                                        indices,num_indices)

  use hdf5
  
  use Option_module
  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscInt :: local_ghosted_id, local_id, natural_id
  PetscInt :: index_count
  PetscInt :: cell_count
  integer(HSIZE_T) :: num_cells_in_file
  PetscInt :: temp_int, i
  
  PetscMPIInt, allocatable :: cell_ids(:)
  PetscInt, allocatable :: temp(:)
  
  PetscInt :: read_block_size
  PetscInt :: indices_array_size

  call PetscLogEventBegin(logging%event_map_indices_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  read_block_size = HDF5_READ_BUFFER_SIZE
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_cells_in_file,hdf5_err)
  if (dataset_size > 0 .and. num_cells_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimension",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_cells_in_file, dataset_size
    call printErrMsg(option)    
  endif
  
  allocate(cell_ids(read_block_size))
  if (num_indices > 0) then
    allocate(indices(num_indices))
  else
    indices_array_size = read_block_size
    allocate(indices(indices_array_size))
  endif
  
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
    if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,cell_ids,dims,hdf5_err, &
                     memory_space_id,file_space_id,prop_id)                     
      call PetscLogEventEnd(logging%event_h5dread_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
#ifdef HDF5_BROADCAST
    endif
    if (option%mycommsize > 1) &
      call mpi_bcast(cell_ids,dims(1),MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
#endif     
  call PetscLogEventBegin(logging%event_hash_map, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    do i=1,dims(1)
      cell_count = cell_count + 1
      natural_id = cell_ids(i)
      local_ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (local_ghosted_id > 0) then
        local_id = grid%nG2L(local_ghosted_id)
        if (local_id > 0) then
          index_count = index_count + 1
          if (index_count > indices_array_size .and. num_indices <= 0) then
            !reallocate if array grows too large and num_indices <= 0
            allocate(temp(index_count))
            temp(1:index_count) = indices(1:index_count)
            deallocate(indices)
            indices_array_size = 2*indices_array_size
            allocate(indices(indices_array_size))
            indices = 0
            indices(1:index_count) = temp(1:index_count)
            deallocate(temp)
          endif
          indices(index_count) = cell_count
        endif
      endif
    enddo
  enddo
  
  call PetscLogEventEnd(logging%event_hash_map, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  deallocate(cell_ids)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  if (num_indices > 0 .and. index_count /= num_indices) then
    write(option%io_buffer, &
          '("Number of indices read (",i9,") does not match the number",&
           &" of indices requested (",i9,").")') &
           index_count, num_indices
    call printErrMsg(option)        
  endif
  
  if (index_count < indices_array_size .and. num_indices <= 0) then
    ! resize to index count
    allocate(temp(index_count))
    temp(1:index_count) = indices(1:index_count)
    deallocate(indices)
    allocate(indices(index_count))
    indices(1:index_count) = temp(1:index_count)
    deallocate(temp)
  endif
  
  if (num_indices <= 0) num_indices = index_count

  call PetscLogEventEnd(logging%event_map_indices_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine HDF5MapLocalToNaturalIndices

! ************************************************************************** !
!
! HDF5ReadRealArray: Read in local real values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5ReadRealArray(option,file_id,dataset_name,dataset_size, &
                             indices,num_indices,real_array)

  use hdf5
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt :: indices(:)
  PetscInt :: num_indices
  PetscReal :: real_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscInt :: index_count
  PetscInt :: real_count, prev_real_count
  integer(HSIZE_T) :: num_reals_in_file
  PetscInt :: temp_int, i, index
  
  PetscReal, allocatable :: real_buffer(:)
  
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_read_real_array_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
  read_block_size = HDF5_READ_BUFFER_SIZE
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_file,hdf5_err)
#if 0
  if (dataset_size > 0 .and. num_reals_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_reals_in_file, grid%nmax
    call printErrMsg(option)   
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
        if (option%myrank == option%io_rank) then                           
#endif
          call PetscLogEventBegin(logging%event_h5dread_f, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
          call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)
          call PetscLogEventEnd(logging%event_h5dread_f, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
#ifdef HDF5_BROADCAST
        endif
        if (option%mycommsize > 1) &
          call mpi_bcast(real_buffer,dims(1),MPI_DOUBLE_PRECISION, &
                         option%io_rank,option%mycomm,ierr)
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
    if (option%myrank == io_rank) then 
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
    endif
    if (option%mycommsize > 1) &
      call mpi_bcast(real_buffer,dims(1),MPI_DOUBLE_PRECISION, &
                     option%io_rank,option%mycomm,ierr)
    real_count = real_count + length(1)                  
  enddo
#endif

  deallocate(real_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_real_array_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
end subroutine HDF5ReadRealArray

! ************************************************************************** !
!
! HDF5ReadIntegerArray: Read in local integer values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5ReadIntegerArray(option,file_id,dataset_name,dataset_size, &
                                indices,num_indices,integer_array)

  use hdf5
  
  use Grid_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt :: indices(:)
  PetscInt :: num_indices
  PetscInt :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscInt :: index_count
  PetscInt :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers_in_file
  PetscInt :: temp_int, i, index
  
  PetscMPIInt, allocatable :: integer_buffer(:)
  
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_read_int_array_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  read_block_size = HDF5_READ_BUFFER_SIZE
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
#if 0
  if (dataset_size > 0 .and. num_integers_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_integers_in_file,dataset_size
    call printErrMsg(option)   
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
        if (option%myrank == option%io_rank) then                           
#endif
          call PetscLogEventBegin(logging%event_h5dread_f, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
          call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)   
          call PetscLogEventEnd(logging%event_h5dread_f, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
#ifdef HDF5_BROADCAST
        endif
        if (option%mycommsize > 1) &
          call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,option%io_rank, &
                         option%mycomm,ierr)
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
    if (option%myrank == option%io_rank) then 
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)   
      call PetscLogEventEnd(logging%event_h5dread_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
    endif
    if (option%mycommsize > 1) &
      call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    integer_count = integer_count + length(1)                  
  enddo
#endif

  deallocate(integer_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_int_array_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
end subroutine HDF5ReadIntegerArray

! ************************************************************************** !
!
! HDF5WriteIntegerArray: Writes an integer array to an hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
subroutine HDF5WriteIntegerArray(option,dataset_name,dataset_size,file_id, &
                                 indices,num_indices,integer_array)

  use hdf5
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  integer(HID_T) :: file_id
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: dataset_size 
  PetscInt :: indices(:)
  PetscInt :: num_indices
  PetscInt :: integer_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscInt :: index_count
  PetscInt :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers_in_file
  PetscInt :: temp_int, i, index
  
  PetscMPIInt, allocatable :: integer_buffer(:)
  
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_write_int_array_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                        
  read_block_size = HDF5_READ_BUFFER_SIZE
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_integers_in_file, &
                                      hdf5_err)
#if 0
  if (dataset_size > 0 .and. num_integers_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_integers_in_file,dataset_size
    call printErrMsg(option)   
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
        if (option%myrank == option%io_rank) then                           
#endif
          call PetscLogEventBegin(logging%event_h5dread_f, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                  PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
          call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)   
          call PetscLogEventEnd(logging%event_h5dread_f, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                                PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
#ifdef HDF5_BROADCAST
        endif
        if (option%mycommsize > 1) &
          call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,option%io_rank, &
                         option%mycomm,ierr)
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
    if (option%myrank == option%io_rank) then   
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
                            
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)   
      call PetscLogEventEnd(logging%event_h5dread_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
    endif
    if (option%mycommsize > 1) &
      call mpi_bcast(integer_buffer,dims(1),MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    integer_count = integer_count + length(1)                  
  enddo
#endif

  deallocate(integer_buffer)
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_write_int_array_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                        
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
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  character(len=32) :: name
  type(realization_type) :: realization
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization%patch
  grid => patch%grid
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
  PetscReal :: array(:)
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: istart_local, jstart_local, kstart_local
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  PetscMPIInt :: rank
  
  PetscMPIInt, pointer :: int_array(:)
  PetscReal, pointer :: double_array(:)
  PetscInt :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: ny_local_X_nz_local
  PetscInt :: num_to_write

  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                        
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
    if (data_type == HDF_NATIVE_INTEGER) then
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
      call PetscLogEventBegin(logging%event_h5dwrite_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
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
      call PetscLogEventBegin(logging%event_h5dwrite_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)         
      call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)   
      deallocate(double_array)
#else
      call PetscLogEventBegin(logging%event_h5dwrite_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)   
      call h5dwrite_f(data_set_id,data_type,array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                         
#endif
    endif
    call h5pclose_f(prop_id,hdf5_err)
  endif
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
end subroutine HDF5WriteStructuredDataSet
!GEH - Structured Grid Dependence - End

! ************************************************************************** !
!
! HDF5ReadIndices: Reads cell indices from an hdf5 dataset
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine HDF5ReadIndices(grid,option,file_id,dataset_name,dataset_size, &
                           indices)

  use hdf5
  
  use Option_module
  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  PetscMPIInt, allocatable :: indices_i4(:)
  integer(HSIZE_T) :: num_data_in_file
  
  PetscInt :: istart, iend

  call PetscLogEventBegin(logging%event_read_indices_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                        
  istart = 0  ! this will be zero-based
  iend = 0
  
  ! first determine upper and lower bound on PETSc global array
  call mpi_exscan(grid%nlmax,istart,ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
  call mpi_scan(grid%nlmax,iend,ONE_INTEGER,MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
  if (iend /= istart + grid%nlmax) then
    call printErrMsg(option,'iend /= istart+grid%nlmax')
  endif
  
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)
  if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_data_in_file,dataset_size
    call printErrMsg(option)   
  else
    dataset_size = num_data_in_file
  endif  
  
  if (istart < num_data_in_file) then
  
    allocate(indices_i4(-1:iend-istart))
    allocate(indices(-1:iend-istart))
    indices_i4(-1) = istart
    indices_i4(0) = iend
  
    rank = 1
    offset = 0
    length = 0
    stride = 1
  
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  
    dims = 0
    dims(1) = iend-istart
    memory_space_id = -1
    call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    call PetscLogEventBegin(logging%event_h5dread_f, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
                               
    call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,indices_i4(1:iend-istart), &
                   dims,hdf5_err,memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
    indices(-1:iend-istart) = indices_i4(-1:iend-istart)                
    deallocate(indices_i4)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_indices_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                           
end subroutine HDF5ReadIndices

! ************************************************************************** !
!
! HDF5ReadArray: Read an hdf5 array into a Petsc Vec
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine HDF5ReadArray(discretization,grid,option,file_id,dataset_name, &
                         dataset_size, &
                         indices,global_vec,data_type)
                         
  use hdf5
  
  use Option_module
  use Grid_Module
  use Discretization_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  Vec :: global_vec
  integer(HID_T) :: data_type
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: rank
  integer(HSIZE_T) :: num_data_in_file
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  PetscMPIInt, allocatable :: integer_buffer(:)
  PetscInt, allocatable :: indices0(:)
  
  call PetscLogEventBegin(logging%event_read_array_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
  istart = 0
  iend = 0
  
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)

  if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_data_in_file,grid%nmax
    call printErrMsg(option)   
  endif

  rank = 1
  offset = 0
  length = 0
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif

  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)
  call VecZeroEntries(natural_vec,ierr)

  ! must initialize here to avoid error below when closing memory space
  memory_space_id = -1

  if (associated(indices)) then

    istart = indices(-1)
    iend = indices(0)

    dims = 0
    dims(1) = iend-istart
    call h5screate_simple_f(rank,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    allocate(real_buffer(iend-istart))
    if (data_type == H5T_NATIVE_DOUBLE) then
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(integer_buffer(iend-istart))
      call PetscLogEventBegin(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                              PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)                              
      do i=1,iend-istart
        real_buffer(i) = real(integer_buffer(i))
      enddo
      deallocate(integer_buffer)
    endif
    ! must convert indices to zero based for VecSetValues
    allocate(indices0(iend-istart))
    indices0 = indices(1:iend-istart)-1
    call VecSetValues(natural_vec,iend-istart,indices0, &
                      real_buffer,INSERT_VALUES,ierr) 
    deallocate(indices0)
    deallocate(real_buffer)

  endif

  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call VecAssemblyBegin(natural_vec,ierr)
  call VecAssemblyEnd(natural_vec,ierr)
  call DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec, &
                                     ONEDOF)
  call VecDestroy(natural_vec,ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                            
end subroutine HDF5ReadArray

#endif ! USE_HDF5

! ************************************************************************** !
!
! HDF5ReadRegionFromFile: Reads a region from an hdf5 file
! author: Glenn Hammond
! date: 1/3/08
!
! ************************************************************************** !
subroutine HDF5ReadRegionFromFile(realization,region,filename)

#ifdef USE_HDF5
  use hdf5
#endif
  
  use Realization_module
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  type(region_type) :: region
  character(len=MAXWORDLENGTH) :: filename

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#ifdef USE_HDF5  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
#endif

  PetscInt :: num_indices, i, local_id
  PetscInt, pointer :: indices(:)
  PetscInt, pointer :: integer_array(:)
  
#ifndef USE_HDF5
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with -DUSE_HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else


  option => realization%option
  patch => realization%patch
  grid => patch%grid

  call PetscLogEventBegin(logging%event_region_read_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
  ! create hash table for fast lookup
#ifdef HASH
  call GridCreateNaturalToGhostedHash(grid,option)
#endif

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  allocate(indices(grid%nlmax))

  ! Open the Regions group
  string = 'Regions'
  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)  
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

  ! Open the Regions group
  string = trim(region%name)
  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)  
  
  call h5gopen_f(grp_id,string,grp_id2,hdf5_err)

  ! Read Cell Ids
  string = "Cell Ids"
  ! num_indices <= 0 indicates that the array size is uncertain and
  ! the size will be returned in num_indices
  num_indices = -1
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id2,string,ZERO_INTEGER,indices, &
                                    num_indices)
  allocate(integer_array(num_indices))
  integer_array = 0
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,grp_id2,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)

  ! convert cell ids from natural to local
  call PetscLogEventBegin(logging%event_hash_map, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  do i=1,num_indices
    integer_array(i) = grid%nG2L(GridGetLocalGhostedIdFromHash(grid,integer_array(i))) 
  enddo
  call PetscLogEventEnd(logging%event_hash_map, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  region%cell_ids => integer_array
                            
  allocate(integer_array(num_indices))
  integer_array = 0
  string = "Face Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)  
  call HDF5ReadIntegerArray(option,grp_id2,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)
                            
  region%faces => integer_array
  region%num_cells = num_indices
  deallocate(indices)
  nullify(indices)

  option%io_buffer = 'Closing group: ' // trim(region%name)
  call printMsg(option)  
  call h5gclose_f(grp_id2,hdf5_err)
  option%io_buffer = 'Closing group: Regions'
  call printMsg(option)   
  call h5gclose_f(grp_id,hdf5_err)
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
   
  call h5close_f(hdf5_err)

  call GridDestroyHashTable(grid)
#endif  

  call PetscLogEventEnd(logging%event_region_read_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine HDF5ReadRegionFromFile

! ************************************************************************** !
!
! HDF5ReadMaterialsFromFile: Reads material ids from an hdf5 file
! author: Glenn Hammond
! date: 1/3/08
!
! ************************************************************************** !
subroutine HDF5ReadMaterialsFromFile(realization,filename)

#ifdef USE_HDF5
  use hdf5
#endif
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#ifdef USE_HDF5  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscInt, allocatable :: integer_array(:)
  
  Vec :: global_vec
  Vec :: local_vec
  PetscReal, pointer :: vec_ptr(:)

#ifndef USE_HDF5
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with -DUSE_HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_material_read_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  ! create hash table for fast lookup
#ifdef HASH
  call PetscGetTime(tstart,ierr)
  call GridCreateNaturalToGhostedHash(grid,option)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to create hash.")') tend-tstart
  call printMsg(option) 
#endif

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option) 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                  option)

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open the Materials group
  string = 'Materials'

  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)   
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

! new approach
#if 1
  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscGetTime(tstart,ierr)
  string = "Material Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,grp_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
#else  
  allocate(indices(grid%nlmax))
  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id,string,grid%nmax, &
                                    indices,grid%nlmax)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to map local to natural indices.")') &
    tend-tstart
  call printMsg(option)  

  ! Read Material ids
  allocate(integer_array(grid%nlmax))
  string = "Material Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call PetscGetTime(tstart,ierr)
  call HDF5ReadIntegerArray(option,grp_id,string,grid%nlmax,indices, &
                            grid%nlmax,integer_array)
  call GridCopyIntegerArrayToPetscVec(integer_array,global_vec,grid%nlmax)
  deallocate(integer_array)
#endif
  
  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)
  call GridCopyPetscVecToIntegerArray(patch%imat,local_vec,grid%ngmax)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read material ids.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  option%io_buffer = 'Closing group: Materials'
  call printMsg(option)   
  call h5gclose_f(grp_id,hdf5_err)
    
  option%io_buffer = 'Closing hdf5 file: ' // filename
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
   
  call h5close_f(hdf5_err)
  
  call GridDestroyHashTable(grid)
#endif  

  call PetscLogEventEnd(logging%event_material_read_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
end subroutine HDF5ReadMaterialsFromFile

! ************************************************************************** !
!
! HDF5ReadPermeabilitiesFromFile: Reads permeabilities from an hdf5 file
! author: Glenn Hammond
! date: 01/16/09
!
! ************************************************************************** !
subroutine HDF5ReadPermeabilitiesFromFile(realization,filename)

#ifdef USE_HDF5
  use hdf5
#endif
  
  use Realization_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#ifdef USE_HDF5  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscReal, allocatable :: real_array(:)
  
  Vec :: global_vec
  Vec :: local_vec
  PetscReal, pointer :: vec_ptr(:)

#ifndef USE_HDF5
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with -DUSE_HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_permeabilities_read_hdf5, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  ! create hash table for fast lookup
#ifdef HASH
  call PetscGetTime(tstart,ierr)
  call GridCreateNaturalToGhostedHash(grid,option)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to create hash.")') tend-tstart
  call printMsg(option) 
#endif

  ! initialize fortran hdf5 interface
  call h5open_f(hdf5_err)

  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option) 
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option)
  call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                  option)

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open the Materials group
  string = 'Permeabilities'

  option%io_buffer = 'Opening group: ' // trim(string)
  call printMsg(option)   
  call h5gopen_f(file_id,string,grp_id,hdf5_err)

! new approach
#if 1
  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscGetTime(tstart,ierr)
  string = "Permeabilities"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,grp_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
#else  
  allocate(indices(grid%nlmax))
  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5MapLocalToNaturalIndices(grid,option,grp_id,string,grid%nmax, &
                                    indices,grid%nlmax)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to map local to natural indices.")') &
    tend-tstart
  call printMsg(option)  

  ! Read Material ids
  allocate(real_array(grid%nlmax))
  string = "Permeabilities"
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call PetscGetTime(tstart,ierr)
  call HDF5ReadRealArray(option,grp_id,string,grid%nlmax,indices, &
                         grid%nlmax,real_array)
  call GridCopyRealArrayToPetscVec(real_array,global_vec,grid%nlmax)
  deallocate(real_array)
#endif
  
  call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)
  call GridCopyPetscVecToIntegerArray(patch%imat,local_vec,grid%ngmax)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read material ids.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  option%io_buffer = 'Closing group: Materials'
  call printMsg(option)   
  call h5gclose_f(grp_id,hdf5_err)
    
  option%io_buffer = 'Closing hdf5 file: ' // filename
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
   
  call h5close_f(hdf5_err)
  
  call GridDestroyHashTable(grid)
#endif  

  call PetscLogEventEnd(logging%event_permeabilities_read_hdf5, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
end subroutine HDF5ReadPermeabilitiesFromFile

end module HDF5_module
