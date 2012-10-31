module HDF5_Aux_module

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  use Logging_module

  implicit none

#include "definitions.h"

  private
  
  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)
  PetscMPIInt :: hdf5_err
  PetscMPIInt :: io_rank_mpi
! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE  
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  public :: HDF5ReadNDimRealArray, &
#ifdef PARALLELIO_LIB
            HDF5ReadDatasetInteger2D, &
            HDF5ReadDatasetReal2D, &
#else
            HDF5ReadDataset, &
            HDF5ReadDatasetMap, &
            HDF5GroupExists, &
#endif ! PARALLELIO_LIB
            HDF5MakeStringCompabible

contains

#endif
  

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !
!
! HDF5ReadNDimRealArray: Read in an n-dimensional array from an hdf5 file
! author: Glenn Hammond
! date: 01/13/10
!
! ************************************************************************** !
subroutine HDF5ReadNDimRealArray(option,file_id,dataset_name,ndims,dims, &
                                 real_array)

  use hdf5
  
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: dataset_name
  integer(HID_T) :: file_id
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_array(:)
  
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: ndims_hdf5
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  PetscMPIInt :: rank_mpi
  PetscInt :: index_count
  integer(HSIZE_T) :: num_reals_in_dataset
  PetscInt :: temp_int, i, index
  PetscMPIInt :: int_mpi
  
  call PetscLogEventBegin(logging%event_read_ndim_real_array_hdf5,ierr)
                          
  call h5dopen_f(file_id,dataset_name,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  ! should be a rank=1 data space
  
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  ndims = ndims_hdf5
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  allocate(dims(ndims))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  dims = dims_h5
  call h5sget_simple_extent_npoints_f(file_space_id,num_reals_in_dataset,hdf5_err)
  temp_int = dims(1)
  do i = 2, ndims
    temp_int = temp_int * dims(i)
  enddo

  rank_mpi = 1
  offset = 0
  length = num_reals_in_dataset
  stride = 1
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err,length)

  allocate(real_array(num_reals_in_dataset))
  real_array = 0.d0
#ifdef HDF5_BROADCAST
  if (option%myrank == option%io_rank) then                           
#endif
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)
    call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,real_array,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
#ifdef HDF5_BROADCAST
  endif
  if (option%mycommsize > 1) then
    int_mpi = num_reals_in_dataset
    call MPI_Bcast(real_array,int_mpi,MPI_DOUBLE_PRECISION, &
                   option%io_rank,option%mycomm,ierr)
  endif
#endif

  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  
  deallocate(dims_h5)
  deallocate(max_dims_h5) 

  call PetscLogEventEnd(logging%event_read_ndim_real_array_hdf5,ierr)
                          
end subroutine HDF5ReadNDimRealArray

! ************************************************************************** !
!
! HDF5ReadDataset: Read in a general dataset from an hdf5 file
! author: Glenn Hammond
! date: 10/25/11
!
! ************************************************************************** !
subroutine HDF5ReadDataset(dataset,option)

  use hdf5
  use Dataset_Aux_module
  use Option_module
  use Units_module
  
  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attribute_id
  integer(HID_T) :: ndims_hdf5
  integer(HID_T) :: atype_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(3)
  integer(HSIZE_T) :: offset(4), length(4), stride(4)
  integer(HSIZE_T) :: num_data_values
  integer(SIZE_T) size_t_int
  PetscInt :: i, temp_int, temp_array(4)
  PetscInt :: num_spatial_dims, time_dim, num_times
  PetscMPIInt :: array_rank_mpi
  PetscBool :: attribute_exists
  PetscBool :: read_times
  PetscReal :: units_conversion
  character(len=MAXWORDLENGTH) :: attribute_name, dataset_name, word

  !TODO(geh): add to event log
  !call PetscLogEventBegin(logging%event_read_datset_hdf5,ierr)

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(dataset%filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(dataset%filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! the dataset is actually stored in a group.  the group contains
  ! a "data" dataset and optionally a "time" dataset.
  option%io_buffer = 'Opening group: ' // trim(dataset%h5_dataset_name)
  call printMsg(option)  
  call h5gopen_f(file_id,dataset%h5_dataset_name,grp_id,hdf5_err)

  ! only want to read on first time through
  if (dataset%data_dim == DIM_NULL) then
    ! read in attributes if they exist
    attribute_name = "Dimension"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim = 1
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdf5_err)
      size_t_int = MAXWORDLENGTH
      call h5tset_size_f(atype_id,size_t_int,hdf5_err)
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,atype_id,word,attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
      ! set dimensionality of dataset
      call DatasetSetDimension(dataset,word)
    else
      option%io_buffer = &
        'Dimension attribute must be included in hdf5 dataset file.'
      call printErrMsg(option)
    endif
    attribute_name = "Discretization"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGetNDimensions(dataset)
      allocate(dataset%discretization(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,dataset%discretization, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    else
      if(dataset%data_dim/=DIM_CELL) then
        option%io_buffer = &
          'Discretization attribute must be included in hdf5 dataset file.'
        call printErrMsg(option)
      endif
    endif
    attribute_name = "Origin"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGetNDimensions(dataset)
      allocate(dataset%origin(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,dataset%origin, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    endif  
    attribute_name = "Transient"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      read_times = PETSC_TRUE
      dataset%buffer => DatasetBufferCreate()
    else
      read_times = PETSC_FALSE
    endif  
    attribute_name = "Cell Centered"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      dataset%cell_centered = PETSC_TRUE
    else
      dataset%cell_centered = PETSC_FALSE
    endif  

    if (read_times) then
      dataset%is_transient = PETSC_TRUE
      units_conversion = 1.d0
      attribute_name = "Time Units"
      call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
      if (attribute_exists) then
        attribute_dim = 1
        call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdf5_err)
        size_t_int = MAXWORDLENGTH
        call h5tset_size_f(atype_id,size_t_int,hdf5_err)
        call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
        call h5aread_f(attribute_id,atype_id,word,attribute_dim,hdf5_err)
        call h5aclose_f(attribute_id,hdf5_err)
        units_conversion = UnitsConvertToInternal(word,option)
      endif
      ! open the "time" dataset, if it exists
      dataset_name = 'time'
      call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
      call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
      call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
      allocate(dataset%buffer%time_array(num_data_values))
      dataset%buffer%time_array = 0.d0
      call PetscLogEventBegin(logging%event_h5dread_f,ierr)
      array_rank_mpi = 1
      offset(1) = 0
      length(1) = num_data_values
      stride(1) = 1
      call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
      call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,dataset%buffer%time_array, &
                     length,hdf5_err,memory_space_id,file_space_id,prop_id)
      dataset%buffer%time_array = dataset%buffer%time_array * units_conversion
      call PetscLogEventEnd(logging%event_h5dread_f,ierr)  
      call h5pclose_f(prop_id,hdf5_err)
      if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
      call h5sclose_f(file_space_id,hdf5_err)
      call h5dclose_f(dataset_id,hdf5_err)
      dataset%buffer%num_times_total = size(dataset%buffer%time_array)
      dataset%buffer%time_offset = 0
      dataset%buffer%cur_time_index = 1
      ! if maximum buffer size has not been set in the PFLOTRAN input file
      if (dataset%max_buffer_size == 0) then
        attribute_name = "Max Buffer Size"
        call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
        if (attribute_exists) then
          attribute_dim(1) = 1
          call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
          call h5aread_f(attribute_id,H5T_NATIVE_INTEGER,temp_int, &
                         attribute_dim,hdf5_err)
          call h5aclose_f(attribute_id,hdf5_err)
          dataset%buffer%num_times_in_buffer = temp_int
        endif
        if (dataset%buffer%num_times_in_buffer == 0) then
          dataset%buffer%num_times_in_buffer = dataset%buffer%num_times_total
          if (dataset%buffer%num_times_in_buffer > 20) then
            dataset%buffer%num_times_in_buffer = 20
            option%io_buffer = 'Size of dataset buffer truncated to 20.'
            call printMsg(option)
          endif
        endif
      else
        dataset%buffer%num_times_in_buffer = dataset%max_buffer_size
      endif
    endif
  endif ! dataset%data_dim == DIM_NULL
  
  num_spatial_dims = DatasetGetNDimensions(dataset)
  time_dim = -1
  num_times = 1
  if (associated(dataset%buffer)) then
    num_times = dataset%buffer%num_times_total
    time_dim = num_spatial_dims + 1
  endif
  
  
  ! open the "data" dataset
  dataset_name = 'data'
  if (dataset%realization_dependent) then
    write(word,'(i9)') option%id
    dataset_name = trim(dataset_name) // trim(adjustl(word))
  endif
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! get dataset dimensions
  if (.not.associated(dataset%dims)) then
    call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
    allocate(dims_h5(ndims_hdf5))
    allocate(max_dims_h5(ndims_hdf5))
    call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
    allocate(dataset%dims(ndims_hdf5))
    dataset%ndims = ndims_hdf5
    ! have to invert dimensions
    do i = 1, dataset%ndims
      dataset%dims(i) = dims_h5(dataset%ndims-i+1)
    enddo
    deallocate(dims_h5)
    deallocate(max_dims_h5) 
  
    call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
  
    temp_int = dataset%dims(1)
    do i = 2, num_spatial_dims
      temp_int = temp_int * dataset%dims(i)
    enddo
    if (num_data_values/num_times /= temp_int) then
      option%io_buffer = 'Number of values in dataset does not match dimensions.'
      call printErrMsg(option)
    endif
    if (read_times .and. dataset%dims(dataset%ndims) /= num_times) then
      option%io_buffer = 'Number of times does not match last dimension of data array.'
      call printErrMsg(option)
    endif
    dataset%array_size = temp_int ! non buffered array
    allocate(dataset%rarray(dataset%array_size))
    dataset%rarray = 0.d0
    if (associated(dataset%buffer)) then
      ! buffered array
      dataset%buffer%array_size = dataset%array_size* &
                                  dataset%buffer%num_times_in_buffer 
      allocate(dataset%buffer%rarray(dataset%buffer%array_size))
      dataset%buffer%rarray = 0.d0
    endif
  endif

  
  ! call PetscLogEventBegin(logging%event_h5dread_f,ierr)
 
  array_rank_mpi = 1
  length = 1
  ! length (or size) must be adjusted according to the size of the 
  ! remaining data in the file
  if (time_dim > 0) then
    length(1) = min(dataset%buffer%array_size, &
                    (num_times-dataset%buffer%time_offset)* &
                    dataset%array_size)
  else
    length(1) = dataset%array_size
  endif
  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    

  length = 1
  stride = 1
  offset = 0
  do i = 1, num_spatial_dims
    length(i) = dataset%dims(i)
  enddo
  ! cannot read beyond end of the buffer
  if (time_dim > 0) then
    dataset%buffer%num_times_in_buffer = &
      min(dataset%buffer%num_times_in_buffer, &
          num_times-dataset%buffer%time_offset)
    length(time_dim) = dataset%buffer%num_times_in_buffer
    offset(time_dim) = dataset%buffer%time_offset
  endif
  
  !geh: for some reason, we have to invert here.  Perhaps because the
  !     dataset was generated in C???
  temp_array(1:dataset%ndims) = length(1:dataset%ndims)
  do i = 1, dataset%ndims
    length(i) = temp_array(dataset%ndims-i+1)
  enddo
  temp_array(1:dataset%ndims) = offset(1:dataset%ndims)
  do i = 1, dataset%ndims
    offset(i) = temp_array(dataset%ndims-i+1)
  enddo
  ! stride is fine
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                             hdf5_err,stride,stride) 
  if (associated(dataset%buffer)) then
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,dataset%buffer%rarray,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
  else
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,dataset%rarray,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    dataset%rmax = maxval(dataset%rarray)
    dataset%rmin = minval(dataset%rarray)
  endif
  call PetscLogEventEnd(logging%event_h5dread_f,ierr)  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  

  option%io_buffer = 'Closing group: ' // trim(dataset%h5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(dataset%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  
  !TODO(geh): add to event log
  !call PetscLogEventEnd(logging%event_read_ndim_real_array_hdf5,ierr)
                          
end subroutine HDF5ReadDataset

! ************************************************************************** !
!> This routine reads mapping data associated with dataset.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 10/26/12
! ************************************************************************** !
subroutine HDF5ReadDatasetMap(dataset,option)

  use hdf5
  use Dataset_Aux_module
  use Option_module
  use Units_module
  
  implicit none
  
  type(dataset_type) :: dataset
  type(option_type) :: option

  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attribute_id
  integer(HID_T) :: ndims_hdf5
  integer(HID_T) :: atype_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(3)
  integer(HSIZE_T) :: offset(2), length(2), stride(2)
  integer(HSIZE_T) :: num_data_values
  integer(SIZE_T) size_t_int
  PetscInt :: i, temp_int, temp_array(4)
  PetscInt :: num_spatial_dims, time_dim, num_times
  PetscMPIInt :: array_rank_mpi
  PetscBool :: attribute_exists
  PetscBool :: read_times
  PetscReal :: units_conversion
  character(len=MAXWORDLENGTH) :: attribute_name, dataset_name, word

  PetscInt :: nids_local,istart,iend
  PetscInt :: remainder

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(dataset%dataset_map%filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(dataset%dataset_map%filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! Open group
  option%io_buffer = 'Opening group: ' // trim(dataset%dataset_map%h5_dataset_map_name)
  call printMsg(option)  
  call h5gopen_f(file_id,dataset%dataset_map%h5_dataset_map_name,grp_id,hdf5_err)

  ! Open the "data" dataset
  dataset_name = 'data'
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  if (ndims_hdf5 /= 2) then
    option%io_buffer='Dimension of '// trim(dataset%dataset_map%h5_dataset_map_name) // &
      '/Data dataset in ' // trim(dataset%dataset_map%filename) // ' is not equal to 2.'
    call printErrMsg(option)
  endif

  ! Get dimensions of dataset
  allocate(dims_h5(ndims_hdf5))
  allocate(max_dims_h5(ndims_hdf5))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  
  nids_local=dims_h5(2)/option%mycommsize
  remainder =dims_h5(2)-nids_local*option%mycommsize
  if(option%myrank<remainder) nids_local=nids_local+1

  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(nids_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(nids_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
  ! Determine the length and offset of data to be read by each processor
  length(1) = dims_h5(1)
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart

  ! TOdo: Only read part of data
  option%io_buffer='Gautam: Modify code for reading HDF5 to generate mapping '//&
   'of dataset'
  call printMsg(option)
  nids_local=dims_h5(2)
  length(:) = dims_h5(:)
  offset(:) = 0
  
  ! Save dimension size
  dataset%dataset_map%map_dims_global(:) = dims_h5(:)
  dataset%dataset_map%map_dims_local(:) = length(:)
  
  ! Create data space for dataset
  array_rank_mpi=2
  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err)

  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif

  ! Select hyperslab
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,offset,length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(dataset%dataset_map%map(length(1), length(2)))

  ! Read the dataset collectively
  call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, dataset%dataset_map%map, &
                 dims_h5, hdf5_err, memory_space_id, file_space_id,prop_id)

  call h5pclose_f(prop_id,hdf5_err)

  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  

  option%io_buffer = 'Closing group: ' // trim(dataset%dataset_map%h5_dataset_map_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(dataset%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)


end subroutine HDF5ReadDatasetMap

#if defined(PARALLELIO_LIB)

! ************************************************************************** !
!
! HDF5ReadDatasetInteger2D: 
! author: Gautam Bisht
! date: 05/13/2010
!
! ************************************************************************** !
subroutine HDF5ReadDatasetInteger2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)

  use hdf5
  use Option_module
  
  implicit none
  
#if defined(PARALLELIO_LIB)
  include "piof.h"  
#endif

  ! in
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  integer                        :: read_option
  type(option_type)              :: option
  
  ! out
  PetscInt,pointer               :: data(:,:)
  PetscInt                       :: data_dims(2)
  PetscInt                       :: dataset_dims(2)
  
  ! local
  PetscInt :: file_id
  PetscInt :: ndims
  PetscInt :: ii, remainder

  PetscErrorCode :: ierr
  
  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call parallelIO_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call parallelIO_read_dataset( data, PIO_INTEGER, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)
  
  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call parallelIO_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetInteger2D
#endif ! PARALLELIO_LIB

#if defined(PARALLELIO_LIB)

! ************************************************************************** !
!
! : 
! author: 
! date: 
!
! ************************************************************************** !

subroutine HDF5ReadDatasetReal2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)

  use hdf5
  use Option_module
  
  implicit none
  
#if defined(PARALLELIO_LIB)
  include "piof.h"  
#endif

  ! in
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  integer                        :: read_option
  type(option_type)              :: option
  
  ! out
  PetscReal,pointer              :: data(:,:)
  PetscInt                       :: data_dims(2)
  PetscInt                       :: dataset_dims(2)
  
  ! local
  integer :: file_id
  integer :: ndims
  PetscInt :: ii, remainder

  PetscErrorCode :: ierr
  
  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  call parallelIO_open_file(filename, option%ioread_group_id, FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call parallelIO_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call parallelIO_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call parallelIO_read_dataset( data, PIO_DOUBLE, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)

  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call parallelIO_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetReal2D
#endif ! PARALLELIO_LIB

! ************************************************************************** !
!
! HDF5GroupExists: Returns true if a group exists
! author: Glenn Hammond
! date: 03/26/2012
!
! ************************************************************************** !
function HDF5GroupExists(filename,group_name,option)

  use hdf5
  use Option_module
  
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  type(option_type) :: option

  integer(HID_T) :: file_id  
  integer(HID_T) :: grp_id  
  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0 
  PetscBool :: group_exists
  
  PetscBool :: HDF5GroupExists

  ! open the file
  call h5open_f(hdf5_err)
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  option%io_buffer = 'Testing group: ' // trim(group_name)
  call printMsg(option)
  ! I turn off error messaging since if the group does not exist, an error
  ! will be printed, but the user does not need to see this.
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,group_name,grp_id,hdf5_err)
  group_exists = .not.(hdf5_err < 0)
  call h5eset_auto_f(ON,hdf5_err)  

  if (group_exists) then
    HDF5GroupExists = PETSC_TRUE
    call h5gclose_f(grp_id,hdf5_err)  
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" found in file.'
  else
    HDF5GroupExists = PETSC_FALSE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" not found in file.  Therefore, assuming a ' // &
      'cell-indexed dataset.'
  endif
  call printMsg(option)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)  
  
end function HDF5GroupExists

#endif ! defined(PETSC_HAVE_HDF5)

! ************************************************************************** !
!
! HDF5MakeStringCompabible: Replaces '/' in string with '_' for hdf5 names
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine HDF5MakeStringCompabible(name)

  implicit none
  
  character(len=*) :: name
  
  PetscInt :: len, ichar
  
  len = len_trim(name)
  do ichar = 1, len
    if (name(ichar:ichar) == '/') then
      name(ichar:ichar) = '_'
    endif
  enddo
  
  name = trim(name)

end subroutine HDF5MakeStringCompabible

end module HDF5_Aux_module
