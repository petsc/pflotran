module HDF5_module

  use Logging_module

  implicit none

#include "definitions.h"

  private
  
  PetscErrorCode :: ierr

  PetscBool, public :: trick_hdf5 = PETSC_FALSE

#if defined(PARALLELIO_LIB)
  include "piof.h"  
#endif

#if defined(PETSC_HAVE_HDF5)
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: io_rank
      
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

#endif
      
#if defined(PETSC_HAVE_HDF5)
  public :: HDF5MapLocalToNaturalIndices, &
            HDF5ReadIntegerArray, &
            HDF5ReadRealArray, &
            HDF5WriteStructDataSetFromVec, &
            HDF5WriteUnstructuredDataSetFromVec, &
            HDF5WriteStructuredDataSet, &
            HDF5ReadRegionFromFile, &       
            HDF5ReadUnstructuredGridRegionFromFile, &
            HDF5ReadCellIndexedIntegerArray, & 
            HDF5ReadCellIndexedRealArray
#else
  public :: HDF5ReadRegionFromFile, &
            HDF5ReadUnstructuredGridRegionFromFile, &
            HDF5ReadCellIndexedIntegerArray, &
            HDF5ReadCellIndexedRealArray
#endif            

contains

#if defined(PETSC_HAVE_HDF5)

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

#if defined(PARALLELIO_LIB)

! Using Parallel IO library
 
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer:: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices

  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3)
  integer:: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  PetscInt :: local_ghosted_id, local_id, natural_id
  PetscInt :: index_count
  PetscInt :: cell_count
  integer:: num_cells_in_file
  PetscInt :: temp_int, i
  PetscMPIInt :: int_mpi
  
  ! Have to use 'integer' to satisfy HDF5.  PetscMPIInt works, but only if
  ! it is defined as an integer
  integer, allocatable :: cell_ids_i4(:)
  PetscInt, allocatable :: temp(:)
  
  PetscInt :: read_block_size
  PetscInt :: indices_array_size

  call PetscLogEventBegin(logging%event_map_indices_hdf5,ierr)

  read_block_size = HDF5_READ_BUFFER_SIZE

  call parallelIO_get_dataset_size( num_cells_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  !>>>> get size of dataset call h5sget_simple_extent_npoints_f(file_space_id,num_cells_in_file,hdf5_err)
  !if (dataset_size > 0 .and. num_cells_in_file /= dataset_size) then
  !  write(option%io_buffer, &
  !        '(a," data space dimension (",i9,") does not match the dimension",&
  !         &" of the domain (",i9,").")') trim(dataset_name), &
  !         num_cells_in_file, dataset_size
  !  call printErrMsg(option)    
  !endif
  !
  allocate(cell_ids_i4(read_block_size))
  if (num_indices > 0) then
    allocate(indices(num_indices))
  else
    indices_array_size = read_block_size
    allocate(indices(indices_array_size))
  endif
  
  rank_mpi = 1 ! This is in fact number of dimensions
  offset = 0
  length = 0
  stride = 1
  
  dims = 0
  cell_count = 0
  index_count = 0
  memory_space_id = -1
  do
    if (cell_count >= num_cells_in_file) exit
    temp_int = num_cells_in_file-cell_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      dims(1) = temp_int
    endif
    offset(1) = cell_count
    length(1) = dims(1)

    call PetscLogEventBegin(logging%event_h5dread_f,ierr)

    ! rank_mpi = 1 ! This is in fact number of dimensions
    call parallelIO_read_same_sub_dataset(cell_ids_i4, PIO_INTEGER, rank_mpi, dims, & 
            offset, file_id, dataset_name, option%ioread_group_id, ierr)

    !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,cell_ids_i4,dims,hdf5_err, &
                    !memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)
        
    call PetscLogEventBegin(logging%event_hash_map,ierr)

    do i=1,dims(1)
      cell_count = cell_count + 1
      natural_id = cell_ids_i4(i)
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
  
  call PetscLogEventEnd(logging%event_hash_map,ierr)

  deallocate(cell_ids_i4)
  
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

  call PetscLogEventEnd(logging%event_map_indices_hdf5,ierr)
! End of ParallelIO library

#else 

! Default & Glenn's HDF5 Broadcast Mechanism 
 
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
  PetscMPIInt :: rank_mpi
  PetscInt :: local_ghosted_id, local_id, natural_id
  PetscInt :: index_count
  PetscInt :: cell_count
  integer(HSIZE_T) :: num_cells_in_file
  PetscInt :: temp_int, i
  PetscMPIInt :: int_mpi
  
  ! Have to use 'integer' to satisfy HDF5.  PetscMPIInt works, but only if
  ! it is defined as an integer
  integer, allocatable :: cell_ids_i4(:)
  PetscInt, allocatable :: temp(:)
  
  PetscInt :: read_block_size
  PetscInt :: indices_array_size

  call PetscLogEventBegin(logging%event_map_indices_hdf5,ierr)

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
  
  allocate(cell_ids_i4(read_block_size))
  if (num_indices > 0) then
    allocate(indices(num_indices))
  else
    indices_array_size = read_block_size
    allocate(indices(indices_array_size))
  endif
  
  rank_mpi = 1
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
    temp_int = int(num_cells_in_file)-cell_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = cell_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                               hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,cell_ids_i4,dims,hdf5_err, &
                     memory_space_id,file_space_id,prop_id)                     
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)
#ifdef HDF5_BROADCAST
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
      call MPI_Bcast(cell_ids_i4,int_mpi,MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
#endif
        
    call PetscLogEventBegin(logging%event_hash_map,ierr)

    do i=1,int(dims(1))
      cell_count = cell_count + 1
      natural_id = cell_ids_i4(i)
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
  
  call PetscLogEventEnd(logging%event_hash_map,ierr)

  deallocate(cell_ids_i4)
  
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

  call PetscLogEventEnd(logging%event_map_indices_hdf5,ierr)
! End of Default & Glenn's HDF5 Broadcast Mechanism

#endif
! PARALLELIO_LIB

end subroutine HDF5MapLocalToNaturalIndices
   
! ************************************************************************** !
!
! HDF5ReadRealArray: Read in local real values from hdf5 global file
! author: Glenn Hammond
! date: 09/21/07
!
! ************************************************************************** !
!TODO(geh): this is only used by general_grid, which should be deprecated.
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
  PetscReal, pointer :: real_array(:)
  
#if defined(PARALLELIO_LIB)    
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3)
  integer:: offset(3), length(3), stride(3)
  integer:: num_reals_in_file
#else
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  integer(HSIZE_T) :: num_reals_in_file
#endif

  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  PetscInt :: real_count, prev_real_count
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  PetscReal, allocatable :: real_buffer(:)
  
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_read_real_array_hdf5,ierr)
                          
#if defined(PARALLELIO_LIB)    
  read_block_size = HDF5_READ_BUFFER_SIZE
  ! should be a rank=1 data space (i.e., one dimensional dataset)
  call parallelIO_get_dataset_size( num_reals_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
!???? get size of dataset  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_file,hdf5_err)
#if 0
  if (dataset_size > 0 .and. num_reals_in_file /= dataset_size) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") does not match the dimensions",&
           &" of the domain (",i9,").")') trim(dataset_name), &
           num_reals_in_file, grid%nmax
    call printErrMsg(option)   
  endif
#endif
  
  rank_mpi = 1
  offset = 0
  length = 0
  stride = 1
  
  dims = 0
  real_count = 0
  prev_real_count = 0
  index_count = 0
  memory_space_id = -1
  ! if a list of indices exists, use the indexed approach
  if (num_indices > 0) then
    allocate(real_buffer(read_block_size))
    do i=1,num_indices
      index = indices(i)
      if (index > real_count) then
        do 
          if (index <= real_count) exit
          temp_int = int(num_reals_in_file-real_count)
          temp_int = min(temp_int,read_block_size)
          if (dims(1) /= temp_int) then
            dims(1) = temp_int
          endif
          ! offset is zero-based
          offset(1) = real_count
          length(1) = dims(1)

        call PetscLogEventBegin(logging%event_h5dread_f,ierr)

        ! rank_mpi = 1 ! This is in fact number of dimensions
        call parallelIO_read_same_sub_dataset(real_buffer, PIO_DOUBLE, rank_mpi, dims, & 
                offset, file_id, dataset_name, option%ioread_group_id, ierr)

        !call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                       !hdf5_err,memory_space_id,file_space_id,prop_id)
        call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
          prev_real_count = real_count
          real_count = real_count + int(length(1))
        enddo
      endif
      real_array(i) = real_buffer(index-prev_real_count)
    enddo

    deallocate(real_buffer)
  ! otherwise, read the entire array
  else
    if (.not.associated(real_array)) then
      allocate(real_array(num_reals_in_file))
    endif
    real_array = 0.d0
    dims(1) = num_reals_in_file
    offset(1) = 0
    length(1) = dims(1)
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)
      call parallelIO_read_same_sub_dataset(real_buffer, PIO_DOUBLE, rank_mpi, dims, & 
              offset, file_id, dataset_name, option%ioread_group_id, ierr)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
  endif
  
#else
!PARALLELIO_LIB is not defined

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
  
  rank_mpi = 1
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
  ! if a list of indices exists, use the indexed approach
  if (num_indices > 0) then
    allocate(real_buffer(read_block_size))
    do i=1,num_indices
      index = indices(i)
      if (index > real_count) then
        do 
          if (index <= real_count) exit
          temp_int = int(num_reals_in_file)-real_count
          temp_int = min(temp_int,read_block_size)
          if (dims(1) /= temp_int) then
            if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
            dims(1) = temp_int
            call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
          endif
          ! offset is zero-based
          offset(1) = real_count
          length(1) = dims(1)
          call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                     length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
          if (option%myrank == option%io_rank) then                           
#endif
            call PetscLogEventBegin(logging%event_h5dread_f,ierr)
            call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                           hdf5_err,memory_space_id,file_space_id,prop_id)
            call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
#ifdef HDF5_BROADCAST
          endif
          if (option%mycommsize > 1) then
            int_mpi = dims(1)
            call MPI_Bcast(real_buffer,int_mpi,MPI_DOUBLE_PRECISION, &
                           option%io_rank,option%mycomm,ierr)
          endif
#endif
          prev_real_count = real_count
          real_count = real_count + int(length(1))
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
        call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
      endif
      ! offset is zero-based
      offset(1) = real_count
      length(1) = dims(1)
      call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                 length,hdf5_err,stride,stride) 
      if (option%myrank == io_rank) then 
        call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
        call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                       hdf5_err,memory_space_id,file_space_id,prop_id)
        call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
      endif
      if (option%mycommsize > 1) then
        int_mpi = dims(1)
        call MPI_Bcast(real_buffer,int_mpi,MPI_DOUBLE_PRECISION, &
                       option%io_rank,option%mycomm,ierr)
      endif
      real_count = real_count + int(length(1))
    enddo
#endif
    deallocate(real_buffer)
  ! otherwise, read the entire array
  else
    if (.not.associated(real_array)) then
      allocate(real_array(num_reals_in_file))
    endif
    real_array = 0.d0
    dims(1) = num_reals_in_file
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
    offset(1) = 0
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
    if (option%myrank == option%io_rank) then                           
#endif
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_array,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
#ifdef HDF5_BROADCAST
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
      call MPI_Bcast(real_array,int_mpi,MPI_DOUBLE_PRECISION, &
                     option%io_rank,option%mycomm,ierr)
    endif
#endif
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
#endif
!PARALLELIO_LIB

  call PetscLogEventEnd(logging%event_read_real_array_hdf5,ierr)
                          
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

#if defined(PARALLELIO_LIB)
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  PetscInt :: indices(:)
  PetscInt :: num_indices
  PetscInt :: integer_array(:)
  
  integer :: file_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3)
  integer:: offset(3), length(3), stride(3)
  integer:: num_integers


  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  PetscInt :: integer_count, prev_integer_count
  PetscInt :: num_integers_in_file
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt :: read_block_size

  call PetscLogEventBegin(logging%event_read_int_array_hdf5,ierr)

  read_block_size = HDF5_READ_BUFFER_SIZE
  allocate(integer_buffer_i4(read_block_size))
  dims = 0
  integer_count = 0
  prev_integer_count = 0
  index_count = 0
  memory_space_id = -1
  length = 0
  num_integers_in_file = 0

  call parallelIO_get_dataset_size( num_integers, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  num_integers_in_file = int(num_integers) 
#if 0  
     if (dataset_size > 0 .and. num_integers_in_file /= dataset_size) then
         write(option%io_buffer, &
              '(a," data space dimension (",i9,") does not match the dimensions",&
               &" of the domain (",i9,").")') trim(dataset_name), &
               num_integers_in_file,dataset_size
         call printErrMsg(option)   
     endif
#endif

  rank_mpi = 1
  offset = 0
  stride = 1
                  
  do i=1,num_indices
    index = indices(i)
    if (index > integer_count) then
      do
        if (index <= integer_count) exit
        temp_int = num_integers_in_file-integer_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          dims(1) = temp_int
          length(1) = dims(1)
        endif
           ! offset is zero-based
           offset(1) = integer_count
           length(1) = dims(1)
           call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
           call parallelIO_read_same_sub_dataset(integer_buffer_i4, PIO_INTEGER, rank_mpi, dims, & 
                offset, file_id, dataset_name, option%ioread_group_id, ierr)
           !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                          !hdf5_err,memory_space_id,file_space_id,prop_id)   
           call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
        prev_integer_count = integer_count
        integer_count = integer_count + length(1)                  
      enddo
    endif
    integer_array(i) = integer_buffer_i4(index-prev_integer_count)
  enddo

  do
    if (integer_count >= num_integers_in_file) exit
    temp_int = num_integers_in_file-integer_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      dims(1) = temp_int
      length(1) = dims(1)
    endif
    if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
       ! offset is zero-based
       offset(1) = integer_count
       length(1) = dims(1)
       call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
       call parallelIO_read_same_sub_dataset(integer_buffer_i4, PIO_INTEGER, rank_mpi, dims, & 
                offset, file_id, dataset_name, option%ioread_group_id, ierr)
       !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                      !hdf5_err,memory_space_id,file_space_id,prop_id)   
       call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    endif 
    integer_count = integer_count + length(1)                  
  enddo
  deallocate(integer_buffer_i4)

  call PetscLogEventEnd(logging%event_read_int_array_hdf5,ierr)

#else
! PARALLELIO_LIB is not defined

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
  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  PetscInt :: integer_count, prev_integer_count
  integer(HSIZE_T) :: num_integers_in_file
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  integer, allocatable :: integer_buffer_i4(:)
  
  PetscInt :: read_block_size

! Default & Glenn's HDF5 Broadcast Mechanism (uses HDF5 Independent I/O mode)

  call PetscLogEventBegin(logging%event_read_int_array_hdf5,ierr)

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
  
  allocate(integer_buffer_i4(read_block_size))
  
  rank_mpi = 1
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
        temp_int = int(num_integers_in_file)-integer_count
        temp_int = min(temp_int,read_block_size)
        if (dims(1) /= temp_int) then
          if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
          dims(1) = temp_int
          call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
        endif
        ! offset is zero-based
        offset(1) = integer_count
        length(1) = dims(1)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                                   length,hdf5_err,stride,stride) 
#ifdef HDF5_BROADCAST
        if (option%myrank == option%io_rank) then                           
#endif
          call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
          call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                         hdf5_err,memory_space_id,file_space_id,prop_id)   
          call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
#ifdef HDF5_BROADCAST
        endif
        if (option%mycommsize > 1) then
          int_mpi = dims(1)
          call MPI_Bcast(integer_buffer_i4,int_mpi,MPI_INTEGER,option%io_rank, &
                         option%mycomm,ierr)
        endif
#endif
        prev_integer_count = integer_count
        integer_count = integer_count + int(length(1))
      enddo
    endif
    integer_array(i) = integer_buffer_i4(index-prev_integer_count)
  enddo

#ifdef HDF5_BROADCAST
  do
    if (integer_count >= num_integers_in_file) exit
    temp_int = num_integers_in_file-integer_count
    temp_int = min(temp_int,read_block_size)
    if (dims(1) /= temp_int) then
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      dims(1) = temp_int
      call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
    endif
    ! offset is zero-based
    offset(1) = integer_count
    length(1) = dims(1)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    if (option%myrank == option%io_rank) then 
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)   
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    endif
    if (option%mycommsize > 1) then
      int_mpi = dims(1)
      call MPI_Bcast(integer_buffer_i4,int_mpi,MPI_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
    integer_count = integer_count + length(1)                  
  enddo
#endif
  deallocate(integer_buffer_i4)
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_int_array_hdf5,ierr)

! Default & Glenn's HDF5 Broadcast Mechanism (uses HDF5 Independent I/O mode)
                          
#endif
! PARALLELIO_LIB

end subroutine HDF5ReadIntegerArray

! ************************************************************************** !
!
! HDF5WriteStructuredDataSet: Writes data from an array into HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine HDF5WriteStructuredDataSet(name,array,file_id,data_type,option, &
                                      nx_global,ny_global,nz_global, &
                                      nx_local,ny_local,nz_local, &
                                      istart_local,jstart_local,kstart_local)

  use hdf5
  use Option_module
  
  implicit none
  
  type(option_type) :: option

  character(len=*) :: name
  PetscReal :: array(:)

#if defined(PARALLELIO_LIB_WRITE)    
  integer:: file_id
  integer:: data_type
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer :: dims(3),mem_dims(3)
  integer :: start(3), length(3), stride(3)
  integer :: nx_local, ny_local, nz_local, nlmax ! Sarat added nlmax 
  integer :: nx_global, ny_global, nz_global
  integer :: istart_local, jstart_local, kstart_local
  integer :: rank_mpi,file_space_rank_mpi
  integer :: i, j, k, count, id
  integer :: ny_local_X_nz_local
  integer :: num_to_write_mpi
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3),mem_dims(3)
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: nx_global, ny_global, nz_global
  PetscInt :: istart_local, jstart_local, kstart_local
  PetscMPIInt :: rank_mpi,file_space_rank_mpi  
  PetscInt :: i, j, k, count, id
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscInt :: ny_local_X_nz_local
  PetscMPIInt :: num_to_write_mpi
#endif

  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscMPIInt :: hdf5_flag

  integer, pointer :: int_array_i4(:)
  PetscReal, pointer :: double_array(:)

#if defined(PARALLELIO_LIB_WRITE)

!  write(option%io_buffer,'(" Writing dataset block: ", a, " type - ", i, ".")') trim(name), data_type
!  call printMsg(option)
!  write(option%io_buffer,'(" HDF_NATIVE_INTEGER is ", i, " and H5T_NATIVE_DOUBLE is ", i, " and H5T_NATIVE_INTEGER is ", i, ".")') HDF_NATIVE_INTEGER, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER
!  call printMsg(option)

  name = trim(name) // CHAR(0)
  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5,ierr)
  
  ny_local_X_nz_local = ny_local*nz_local
  num_to_write_mpi = nx_local*ny_local_X_nz_local
  
  ! file space which is a 3D block
  rank_mpi = 3

  dims(1) = nx_global
  dims(2) = ny_global
  dims(3) = nz_global

  !write(option%io_buffer,'(" Writing dataset block with dimensions: ", i,i,i, ".")') dims(1), dims(2), dims(3)
  !call printMsg(option)

  ! create the hyperslab
  start(1) = istart_local
  start(2) = jstart_local
  start(3) = kstart_local
  length(1) =  nx_local
  length(2) =  ny_local
  length(3) =  nz_local

  ! Sarat added this to eliminate grid%nlmax dependency in the non-INVERT case
  ! below
  nlmax = nx_local * ny_local * nz_local

  if (num_to_write_mpi == 0) length(1) = 1
  stride = 1

    ! write the data
  if (num_to_write_mpi > 0) then
    if (data_type == H5T_NATIVE_DOUBLE) then
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
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)         
       !write(option%io_buffer, &
       !   '(a," Writing double dataset1: dimensions: ",i9,i9,i9, " Data type and ndims: ",i9, i9)') & 
       !trim(name), dims(1), dims(2), dims(3), PIO_DOUBLE, rank_mpi
       !call printMsg(option)   
       call parallelIO_write_dataset_block(double_array, PIO_DOUBLE, rank_mpi, &
              dims, length, start, file_id, name, &
              option%iowrite_group_id, ierr)
      !call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      !hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)   
      deallocate(double_array)
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(int_array_i4(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array_i4(id) = int(array(count))
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)                              
       !write(option%io_buffer, &
       !   '(a," Writing integer dataset1: dimensions: ",i9,i9,i9, " Data type and ndims: ",i9, i9)') & 
       !trim(name), dims(1), dims(2), dims(3), PIO_INTEGER, rank_mpi
       !call printMsg(option)   
      call parallelIO_write_dataset_block(int_array_i4, PIO_INTEGER, rank_mpi, &
              dims, length, start, file_id, name, &
              option%iowrite_group_id, ierr)
      !!call h5dwrite_f(data_set_id,data_type,int_array_i4,dims, &
      !                hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)                              
      deallocate(int_array_i4)
    endif
  endif

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5,ierr)

#else
! PARALLELIO_LIB_WRITE is not defined  

! Default HDF5 Write

  call PetscLogEventBegin(logging%event_write_struct_dataset_hdf5,ierr)
  
  ny_local_X_nz_local = ny_local*nz_local
  num_to_write_mpi = nx_local*ny_local_X_nz_local
  
  ! memory space which is a 1D vector  
  rank_mpi = 1
  dims = 0
  dims(1) = num_to_write_mpi
  if (num_to_write_mpi == 0) dims(1) = 1
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! file space which is a 3D block
  rank_mpi = 3
  ! have to trick hdf5 for now with inverted ordering
  dims(3) = nx_global
  dims(2) = ny_global
  dims(1) = nz_global
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,name,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then 
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,name,data_type,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)
  
  ! create the hyperslab
  start(3) = istart_local
  start(2) = jstart_local
  start(1) = kstart_local
  ! recall these are inverted
  length(3) =  nx_local
  length(2) =  ny_local
  length(1) =  nz_local
  if (num_to_write_mpi == 0) length(1) = 1
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
  if (num_to_write_mpi > 0) then
    if (data_type == HDF_NATIVE_INTEGER) then
      allocate(int_array_i4(nx_local*ny_local*nz_local))
      count = 0
      do k=1,nz_local
        do j=1,ny_local
          do i=1,nx_local
            id = k+(j-1)*nz_local+(i-1)*ny_local_X_nz_local
            count = count+1
            int_array_i4(id) = int(array(count))
          enddo
        enddo
      enddo
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)                              
      call h5dwrite_f(data_set_id,data_type,int_array_i4,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)                              
      deallocate(int_array_i4)
    else
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
      call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)         
      call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                      hdf5_err,memory_space_id,file_space_id,prop_id)  
      call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)   
      deallocate(double_array)
    endif
    call h5pclose_f(prop_id,hdf5_err)
  endif
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)

  call PetscLogEventEnd(logging%event_write_struct_dataset_hdf5,ierr)
                          
#endif
! PARALLELIO_LIB_WRITE vs previous

end subroutine HDF5WriteStructuredDataSet
      
! End of Default HDF5 Write
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

#if defined(PARALLELIO_LIB)

  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
    
#if defined(PARALLELIO_LIB)    
  integer :: file_id
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3)
  integer:: offset(3), length(3), stride(3), globaldims(3)
  integer:: num_data_in_file
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3), globaldims(3)
  integer(HSIZE_T) :: num_data_in_file
#endif
 
  PetscMPIInt :: rank_mpi
  ! seeting to MPIInt to ensure i4
  integer, allocatable :: indices_i4(:)
  
  PetscInt :: istart, iend

  call PetscLogEventBegin(logging%event_read_indices_hdf5,ierr)
                        
  istart = 0  ! this will be zero-based
  iend = 0
  
  ! first determine upper and lower bound on PETSc global array
  call MPI_Scan(grid%nlmax,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                option%mycomm,ierr)
  istart = iend - grid%nlmax
  
  ! should be a rank=1 data space
  call parallelIO_get_dataset_size( num_data_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  globaldims(1) = num_data_in_file 
!???? get size of dataset  call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)
  !if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
  !  write(option%io_buffer, &
  !        '(a," data space dimension (",i9,") does not match the dimensions",&
  !         &" of the domain (",i9,").")') trim(dataset_name), &
  !         num_data_in_file,dataset_size
  !  call printErrMsg(option)   
  !else
  !  dataset_size = num_data_in_file
  !endif  
  
  dataset_size = int(num_data_in_file)

  if (istart < num_data_in_file) then
  
    allocate(indices_i4(-1:iend-istart))
    allocate(indices(-1:iend-istart))
    indices_i4(-1) = istart
    indices_i4(0) = iend
  
    rank_mpi = 1
    offset = 0
    length = 0
    stride = 1
  
    dims = 0
    dims(1) = iend-istart
    memory_space_id = -1

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
                               
    call parallelIO_read_dataset(indices_i4(1:iend-istart), PIO_INTEGER, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, NONUNIFORM_CONTIGUOUS_READ, ierr)
    !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,indices_i4(1:iend-istart), &
                   !dims,hdf5_err,memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    indices(-1:iend-istart) = indices_i4(-1:iend-istart)                
    deallocate(indices_i4)
  endif
  
  call PetscLogEventEnd(logging%event_read_indices_hdf5,ierr)

#else
! PARALLELIO_LIB is not defined

! Default HDF5 Mechanism 

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
  PetscMPIInt :: rank_mpi
  ! seeting to MPIInt to ensure i4
  integer, allocatable :: indices_i4(:)
  integer(HSIZE_T) :: num_data_in_file
  
  PetscInt :: istart, iend

  call PetscLogEventBegin(logging%event_read_indices_hdf5,ierr)
                        
  istart = 0  ! this will be zero-based
  iend = 0
  
  ! first determine upper and lower bound on PETSc global array
  call MPI_Scan(grid%nlmax,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                option%mycomm,ierr)
  istart = iend - grid%nlmax
  
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
    dataset_size = int(num_data_in_file)
  endif  
  
  if (istart < num_data_in_file) then
  
    allocate(indices_i4(-1:iend-istart))
    allocate(indices(-1:iend-istart))
    indices_i4(-1) = istart
    indices_i4(0) = iend
  
    rank_mpi = 1
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
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
                               
    call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,indices_i4(1:iend-istart), &
                   dims,hdf5_err,memory_space_id,file_space_id,prop_id)                     
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    indices(-1:iend-istart) = indices_i4(-1:iend-istart)                
    deallocate(indices_i4)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)

  call PetscLogEventEnd(logging%event_read_indices_hdf5,ierr)
! End of Default HDF5 Mechanism  
  
#endif
! PARALLELIO_LIB
  
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

#if defined(PARALLELIO_LIB)
 
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: dataset_size
  integer(HID_T) :: file_id
  PetscInt, pointer :: indices(:)
  PetscInt :: num_indices
  Vec :: global_vec
  integer:: data_type 
  
  integer:: file_space_id
  integer:: memory_space_id
  integer:: data_set_id
  integer:: prop_id
  integer:: dims(3), globaldims(3)
  integer:: offset(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi
  integer:: num_data_in_file
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt, allocatable :: indices0(:)
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr)
                          
  istart = 0
  iend = 0
  
  ! should be a rank=1 data space

  call parallelIO_get_dataset_size( num_data_in_file, file_id, dataset_name, &
          option%ioread_group_id, ierr)
  globaldims(1) = num_data_in_file
!???? get size   call h5sget_simple_extent_npoints_f(file_space_id,num_data_in_file,hdf5_err)

  !if (dataset_size > 0 .and. num_data_in_file /= dataset_size) then
  !  write(option%io_buffer, &
  !        '(a," data space dimension (",i9,") does not match the dimensions",&
  !         &" of the domain (",i9,").")') trim(dataset_name), &
  !         num_data_in_file,grid%nmax
  !  call printErrMsg(option)   
  !endif

  rank_mpi = 1
  offset = 0
  length = 0
  stride = 1
  
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

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    allocate(real_buffer(iend-istart))
    if (data_type == H5T_NATIVE_DOUBLE) then
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
    
      call parallelIO_read_dataset(real_buffer, PIO_DOUBLE, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, NONUNIFORM_CONTIGUOUS_READ, ierr)
      !call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     !hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(integer_buffer_i4(iend-istart))
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              

      call parallelIO_read_dataset(integer_buffer_i4, PIO_INTEGER, rank_mpi, globaldims, dims, & 
            file_id, dataset_name, option%ioread_group_id, NONUNIFORM_CONTIGUOUS_READ, ierr)
      !call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                     !hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
      do i=1,iend-istart
        real_buffer(i) = real(integer_buffer_i4(i))
      enddo
      deallocate(integer_buffer_i4)
    endif
    ! must convert indices to zero based for VecSetValues
    allocate(indices0(iend-istart))
    indices0 = indices(1:iend-istart)-1
    call VecSetValues(natural_vec,iend-istart,indices0, &
                      real_buffer,INSERT_VALUES,ierr) 
    deallocate(indices0)
    deallocate(real_buffer)

  endif

  call VecAssemblyBegin(natural_vec,ierr)
  call VecAssemblyEnd(natural_vec,ierr)
  call DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec, &
                                     ONEDOF)
  call VecDestroy(natural_vec,ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr)

#else
! PARALLELIO_LIB is not defined

! Default HDF5 Mechanism 
 
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
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: num_data_in_file
  Vec :: natural_vec
  PetscInt :: i, istart, iend
  PetscReal, allocatable :: real_buffer(:)
  integer, allocatable :: integer_buffer_i4(:)
  PetscInt, allocatable :: indices0(:)
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr)
                          
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

  rank_mpi = 1
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
    call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

    ! offset is zero-based
    offset(1) = istart
    length(1) = iend-istart
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                               length,hdf5_err,stride,stride) 
    allocate(real_buffer(iend-istart))
    if (data_type == H5T_NATIVE_DOUBLE) then
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
      call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_buffer,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    else if (data_type == HDF_NATIVE_INTEGER) then
      allocate(integer_buffer_i4(iend-istart))
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
      call h5dread_f(data_set_id,HDF_NATIVE_INTEGER,integer_buffer_i4,dims, &
                     hdf5_err,memory_space_id,file_space_id,prop_id)
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
      do i=1,iend-istart
        real_buffer(i) = real(integer_buffer_i4(i))
      enddo
      deallocate(integer_buffer_i4)
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
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr)
! End of Default HDF5 Mechanism

#endif
! PARALLELIO_LIB

end subroutine HDF5ReadArray

#endif
!PETSC_HAVE_HDF5

! ************************************************************************** !
!
! HDF5ReadRegionFromFile: Reads a region from an hdf5 file
! author: Glenn Hammond
! date: 1/3/08
!
! ************************************************************************** !
subroutine HDF5ReadRegionFromFile(realization,region,filename)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_class
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  type(region_type) :: region
  character(len=MAXSTRINGLENGTH) :: filename

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
#endif

  PetscInt :: num_indices, i, local_id
  PetscInt, pointer :: indices(:)
  PetscInt, pointer :: integer_array(:)
  
  option => realization%option

#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  patch => realization%patch
  grid => patch%grid

  call PetscLogEventBegin(logging%event_region_read_hdf5,ierr)
                          
  ! create hash table for fast lookup
#ifdef HASH
  call GridCreateNaturalToGhostedHash(grid,option)
#endif

#if defined(PARALLELIO_LIB)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
      option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
      call printMsg(option)
  endif

  filename = trim(filename) // CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, FILE_READONLY, &
          file_id, ierr)
  string = '/Regions/' // trim(region%name) // '/Cell Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)

 allocate(indices(grid%nlmax))
  ! Read Cell Ids  
  string = 'Regions' // '/' // trim(region%name) // "Cell Ids" // CHAR(0)
  ! num_indices <= 0 indicates that the array size is uncertain and
  ! the size will be returned in num_indices
  num_indices = -1
  call HDF5MapLocalToNaturalIndices(grid,option,file_id,string,ZERO_INTEGER,indices, &
                                    num_indices)
  allocate(integer_array(num_indices))
  integer_array = 0
  string = '/Regions/' // trim(region%name) // '/Cell Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   
  call HDF5ReadIntegerArray(option,file_id,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)

  ! convert cell ids from natural to local
  call PetscLogEventBegin(logging%event_hash_map,ierr)
  do i=1,num_indices
    integer_array(i) = grid%nG2L(GridGetLocalGhostedIdFromHash(grid,integer_array(i))) 
  enddo
  call PetscLogEventEnd(logging%event_hash_map,ierr)
  region%cell_ids => integer_array
                            
  allocate(integer_array(num_indices))
  integer_array = 0
  string = '/Regions/' // trim(region%name) // '/Face Ids' //CHAR(0)
  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)  
  call HDF5ReadIntegerArray(option,file_id,string, &
                            ZERO_INTEGER,indices,num_indices, &
                            integer_array)
                            
  region%faces => integer_array
  region%num_cells = num_indices
  deallocate(indices)
  nullify(indices)

   if (mod(option%myrank,option%hdf5_read_group_size) == 0) then
       option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
       call printMsg(option)   
   endif
   call parallelio_close_file(file_id, option%ioread_group_id, ierr)

  call GridDestroyHashTable(grid)

! PARALLELIO_LIB
#else   

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

 allocate(indices(grid%nlmax))
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
  call PetscLogEventBegin(logging%event_hash_map,ierr)
  do i=1,num_indices
    integer_array(i) = grid%nG2L(GridGetLocalGhostedIdFromHash(grid,integer_array(i))) 
  enddo
  call PetscLogEventEnd(logging%event_hash_map,ierr)
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
! if PARALLELIO_LIB is not defined

#endif
!PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_region_read_hdf5,ierr)

end subroutine HDF5ReadRegionFromFile

! ************************************************************************** !
!
! HDF5ReadUnstructuredGridRegionFromFile: Reads a region from an hdf5 file
!     for unstructured grid
! author: Gautam Bisht
! date: 5/31/11
!
! ************************************************************************** !
! gb: modified the input argument, so it works for realization (subsurface)
!     and surface_realization (surface-flow)
!subroutine HDF5ReadUnstructuredGridRegionFromFile(realization,region,filename)
subroutine HDF5ReadUnstructuredGridRegionFromFile(option,region,filename)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_class
  use Option_module
  use Grid_module
  use Region_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  !type(realization_type)         :: realization
  type(option_type), pointer :: option
  type(region_type)              :: region
  type(region_sideset_type),pointer:: sideset
  character(len=MAXSTRINGLENGTH) :: filename

  ! local
  !type(option_type), pointer :: option
  PetscMPIInt       :: hdf5_err
  PetscMPIInt       :: rank_mpi
  PetscInt          :: ndims
  PetscInt          :: remainder
  PetscInt          :: istart, iend, ii, jj
  PetscInt,pointer  :: int_buffer(:,:)
  character(len=MAXSTRINGLENGTH) :: string

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: data_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: length(2), offset(2)
#endif

  !option => realization%option

#if !defined(PETSC_HAVE_HDF5)
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted unstructured grids.")')
  call printErrMsg(option)
#else
  ! Initialize FORTRAN predefined datatypes
  call h5open_f(hdf5_err)

  ! Setup file access property with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)

#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

  ! Open the file collectively
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)
  
  ! Open dataset
  string = 'Region/'//trim(region%name)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)

  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id,data_space_id,hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id,ndims,hdf5_err)
  if ((ndims > 2).or.(ndims < 1)) then
    option%io_buffer='Dimension of '//string//' dataset in ' // trim(filename) // &
     ' is > 2 or < 1.'
  call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id,dims_h5,max_dims_h5,hdf5_err)
  
  select case(ndims)
    case(1)
    !
    !                    Read Cell IDs
    !
    region%num_cells = int(dims_h5(1)/option%mycommsize)
      remainder = int(dims_h5(1)) - region%num_cells*option%mycommsize
    if (option%myrank < remainder) region%num_cells = region%num_cells + 1
    
    ! Find istart and iend
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_cells,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_cells,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)
  
    ! Determine the length and offset of data to be read by each processor
    length(1) = iend-istart
    offset(1) = istart

    !
    rank_mpi = 1
    memory_space_id = -1
  
    ! Create data space for dataset
    call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)
  
    ! Select hyperslab
    call h5sselect_hyperslab_f(data_space_id,H5S_SELECT_SET_F,offset,length,hdf5_err)
  
    ! Initialize data buffer
    allocate(int_buffer(length(1),1))
    
    ! Create property list
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F,hdf5_err)
#endif
  
    ! Read the dataset collectively
    call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,int_buffer,&
                   dims_h5,hdf5_err,memory_space_id,data_space_id)

    ! allocate array to store vertices for each cell
    allocate(region%cell_ids(region%num_cells))
    region%cell_ids = 0

    ! It is assumed that cell ids in the HDF5 are 1-based. Converting them to
    ! 0-based
    do ii = 1,region%num_cells
      if(int_buffer(ii,1) < 1 ) then
        write(option%io_buffer,'("Cell ids in the HDF5 for region less than 1")')
        call printErrMsg(option)
      endif
      region%cell_ids(ii) = int_buffer(ii,1) - 1
    enddo
    
  case(2)
    !
    !                  (i)  Read Vertices
    !                       OR
    !                  (ii) Cell IDs with face IDs
    !

    !gb: Cells IDs with face IDs - may not work now 12/22/11
    !region%num_verts = dims_h5(2)/option%mycommsize
    !  remainder = dims_h5(2) - region%num_verts*option%mycommsize
    region%sideset => RegionCreateSideset()
    sideset => region%sideset
    sideset%nfaces = int(dims_h5(2)/option%mycommsize)
    remainder = int(dims_h5(2)) - region%num_verts*option%mycommsize

     ! Find istart and iend
     istart = 0
     iend   = 0
     call MPI_Exscan(sideset%nfaces,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
     call MPI_Scan(sideset%nfaces,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)
  
     ! Determine the length and offset of data to be read by each processor
     length(1) = dims_h5(1)
     length(2) = iend-istart
     offset(1) = 0
     offset(2) = istart
  
     !
     rank_mpi = 2
     memory_space_id = -1
  
     ! Create data space for dataset
     call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)
  
     ! Select hyperslab
     call h5sselect_hyperslab_f(data_space_id,H5S_SELECT_SET_F,offset,length,hdf5_err)
  
     ! Initialize data buffer
     allocate(int_buffer(length(1),length(2)))
  
     ! Create property list
     call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
     call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F,hdf5_err)
#endif
  
     ! Read the dataset collectively
     call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,int_buffer,&
          dims_h5,hdf5_err,memory_space_id,data_space_id)
     
     if (dims_h5(1) == 2) then
       !
       ! Input data is: Cell IDs + Face IDs
       !
       region%num_cells = region%num_verts
       allocate(region%cell_ids(region%num_cells))
       allocate(region%faces(region%num_cells))
       region%num_verts = 0
       
       do ii = 1, region%num_cells
         region%cell_ids(ii) = int_buffer(1,ii)
         region%faces(ii) = int_buffer(2,ii)
       enddo
     else
       !
       ! Input data is list of Vertices
       !
       ! allocate array to store vertices for each cell
       allocate(sideset%face_vertices(MAX_VERT_PER_FACE,sideset%nfaces))
       sideset%face_vertices = -999
  
       do ii = 1,sideset%nfaces
        do jj = 2,int_buffer(1,ii)+1
         sideset%face_vertices(jj-1,ii) = int_buffer(jj,ii)
        enddo
       enddo
     endif
  end select
    
  deallocate(dims_h5)
  deallocate(max_dims_h5)
  deallocate(int_buffer)

  call h5pclose_f(prop_id,hdf5_err)
  call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(data_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

#endif
! if defined(PETSC_HAVE_HDF5)

end subroutine HDF5ReadUnstructuredGridRegionFromFile

! ************************************************************************** !
!
! HDF5ReadCellIndexedIntegerArray: Reads an array of integer values from an 
!                                  hdf5 file
! author: Glenn Hammond
! date: 1/3/08; 02/18/09
!
! ************************************************************************** !
subroutine HDF5ReadCellIndexedIntegerArray(realization,global_vec,filename, &
                                           group_name, &
                                           dataset_name,append_realization_id)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscInt, allocatable :: integer_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_int_read_hdf5,ierr)
  
#if defined(PARALLELIO_LIB)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
     option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
     call printMsg(option) 
  end if   
  filename = trim(filename) //CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, FILE_READONLY, &
          file_id, ierr)

  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)

  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // "/Cell Ids" // CHAR(0)
  else
    string = "/Cell Ids" // CHAR(0)
  endif  
    
  call HDF5ReadIndices(grid,option,file_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscGetTime(tstart,ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif

  string = trim(dataset_name) // adjustl(trim(string))
  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // '/' // trim(string) // CHAR(0)
  else
    string = trim(string) // CHAR(0)
  endif  

  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   

  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
  
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read integer array.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)   
    call parallelio_close_file(file_id, option%ioread_group_id, ierr)
  endif

#else
! if PARALLELIO_LIB is not defined

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

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open group if necessary
  if (len_trim(group_name) > 1) then
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)   
    call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  else
    grp_id = file_id
  endif

  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)


  call PetscGetTime(tstart,ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,grp_id,string,grid%nmax, &
                     indices,global_vec,HDF_NATIVE_INTEGER)
  
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read integer array.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (file_id /= grp_id) then
    option%io_buffer = 'Closing group: ' // trim(group_name)
    call printMsg(option)   
    call h5gclose_f(grp_id,hdf5_err)
  endif
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
#endif  
! if PARALLELIO_LIB is not defined

#endif
! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_int_read_hdf5,ierr)
                          
end subroutine HDF5ReadCellIndexedIntegerArray
    
! ************************************************************************** !
!
! HDF5ReadCellIndexedRealArray: Reads an array of real values from an hdf5 file
! author: Glenn Hammond
! date: 01/16/09, 02/18/09
!
! ************************************************************************** !
subroutine HDF5ReadCellIndexedRealArray(realization,global_vec,filename, &
                                        group_name, &
                                        dataset_name,append_realization_id)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  
  use Realization_class
  use Discretization_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization
  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  

  character(len=MAXSTRINGLENGTH) :: string 

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
#endif

  PetscLogDouble :: tstart, tend
  
  PetscInt, pointer :: indices(:)
  PetscReal, allocatable :: real_array(:)
  
#if !defined(PETSC_HAVE_HDF5)
  option => realization%option
  call printMsg(option,'')
  write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                           &"read HDF5 formatted structured grids.")')
  call printErrMsg(option)
#else

  nullify(indices)

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  call PetscLogEventBegin(logging%event_cell_indx_real_read_hdf5,ierr)

#if defined(PARALLELIO_LIB)
  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
     option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
     call printMsg(option) 
  end if   
  filename = trim(filename) //CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, FILE_READONLY, &
          file_id, ierr)


! only new approach (old approach is removed)
  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)

  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // "/Cell Ids" // CHAR(0)
  else
    string = "/Cell Ids" // CHAR(0)
  endif  
    
  call HDF5ReadIndices(grid,option,file_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscGetTime(tstart,ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif

  string = trim(dataset_name) // adjustl(trim(string))
  !if group_name exists
  if (len_trim(group_name) > 1) then
    string = trim(group_name) // '/' // trim(string) // CHAR(0)
  else
    string = trim(string) // CHAR(0)
  endif  

  option%io_buffer = 'Reading dataset: ' // trim(string)
  call printMsg(option)   

  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,H5T_NATIVE_DOUBLE)
  
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read real array.")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (mod(option%myrank,option%hdf5_read_group_size) == 0) then  
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)   
    call parallelio_close_file(file_id, option%ioread_group_id, ierr)
  endif

#else
! if PARALLELIO_LIB is not defined

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

  option%io_buffer = 'Setting up grid cell indices'
  call printMsg(option) 

  ! Open group if necessary
  if (len_trim(group_name) > 1) then
    option%io_buffer = 'Opening group: ' // trim(group_name)
    call printMsg(option)   
    call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  else
    grp_id = file_id
  endif

  ! Read Cell Ids
  call PetscGetTime(tstart,ierr)
  string = "Cell Ids"
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadIndices(grid,option,grp_id,string,grid%nmax,indices)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to set up indices")') tend-tstart
  call printMsg(option)

  call PetscGetTime(tstart,ierr)
  string = ''
  if (append_realization_id) then
    write(string,'(i6)') option%id
  endif
  string = trim(dataset_name) // adjustl(trim(string))
  if (grp_id /= file_id) then
    option%io_buffer = 'Reading dataset: ' // trim(group_name) // '/' &
                       // trim(string)
  else
    option%io_buffer = 'Reading dataset: ' // trim(string)
  endif
  call printMsg(option)   
  call HDF5ReadArray(discretization,grid,option,file_id,string,grid%nmax, &
                     indices,global_vec,H5T_NATIVE_DOUBLE)
  call PetscGetTime(tend,ierr)
  write(option%io_buffer,'(f6.2," Seconds to read real array")') &
    tend-tstart
  call printMsg(option)  

  if (associated(indices)) deallocate(indices)
  nullify(indices)

  if (file_id /= grp_id) then
    option%io_buffer = 'Closing group: ' // trim(group_name)
    call printMsg(option)   
    call h5gclose_f(grp_id,hdf5_err)
  endif
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)   
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

#endif  
! if PARALLELIO_LIB is not defined

#endif
! PETSC_HAVE_HDF5

  call PetscLogEventEnd(logging%event_cell_indx_real_read_hdf5,ierr)
                          
end subroutine HDF5ReadCellIndexedRealArray

#if defined(PETSC_HAVE_HDF5)
! ************************************************************************** !
!
! HDF5WriteStructDataSetFromVec: Writes data from a PetscVec to HDF5 file
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine HDF5WriteStructDataSetFromVec(name,realization_base,vec,file_id,data_type)

  use hdf5
  use Realization_Base_class, only : realization_base_type
  use Grid_module
  use Option_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  character(len=*) :: name
  class(realization_base_type) :: realization_base
  Vec :: vec
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: vec_ptr(:)
  
  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  
  call VecGetArrayF90(vec,vec_ptr,ierr)
!GEH - Structured Grid Dependence - Begin
  call HDF5WriteStructuredDataSet(trim(name), &
                                 vec_ptr,file_id,data_type,option, &
                                  grid%structured_grid%nx, &
                                  grid%structured_grid%ny, &
                                  grid%structured_grid%nz, &
                                  grid%structured_grid%nlx, &
                                  grid%structured_grid%nly, &
                                  grid%structured_grid%nlz, &
                                  grid%structured_grid%lxs, &
                                  grid%structured_grid%lys, &
                                  grid%structured_grid%lzs)
!GEH - Structured Grid Dependence - End
  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine HDF5WriteStructDataSetFromVec

! ************************************************************************** !
!> This routine writes data from a PETSc Vec to HDF5 file for unstructured
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/31/12
! ************************************************************************** !
!subroutine HDF5WriteUnstructuredDataSetFromVec(name,realization,vec,file_id,data_type)
subroutine HDF5WriteUnstructuredDataSetFromVec(name,option,vec,file_id,data_type)

  use hdf5
  use Realization_class
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
  
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3),istart
  PetscInt :: local_size,global_size,i
  PetscInt, pointer :: int_array(:)
  PetscReal, pointer :: double_array(:)

!  patch => realization%patch
!  grid => patch%grid
!  option => realization%option
  
  call VecGetLocalSize(vec,local_size,ierr)
  call VecGetSize(vec,global_size,ierr)
  
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 3D block
  rank_mpi = 1
  dims = 0
  dims(1) = global_size
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,name,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,name,data_type,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = local_size
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


  if (data_type == H5T_NATIVE_DOUBLE) then
    allocate(double_array(local_size))  
    call VecGetArrayF90(vec,vec_ptr,ierr)
    do i=1,local_size
      double_array(i) = vec_ptr(i)
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
    call h5dwrite_f(data_set_id,data_type,double_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

    deallocate(double_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  if (data_type == H5T_NATIVE_INTEGER) then
    allocate(int_array(local_size))
    call VecGetArrayF90(vec,vec_ptr,ierr)
    do i=1,local_size
      int_array(i) = int(vec_ptr(i))
    enddo
    call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
    call PetscLogEventBegin(logging%event_h5dwrite_f,ierr)
    call h5dwrite_f(data_set_id,data_type,int_array,dims, &
                    hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dwrite_f,ierr)

    deallocate(int_array)
    call h5pclose_f(prop_id,hdf5_err)
  endif

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  
end subroutine HDF5WriteUnstructuredDataSetFromVec

#endif
      
end module HDF5_module



