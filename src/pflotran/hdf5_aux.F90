module HDF5_aux_module

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

#ifdef PARALLELIO_LIB
  public :: HDF5ReadNDimRealArray, &
            HDF5ReadDatasetInteger2D, &
            HDF5ReadDatasetReal2D
#else
  public :: HDF5ReadNDimRealArray, &
            HDF5ReadDataset
#endif ! PARALLELIO_LIB

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
  PetscInt :: real_count, prev_real_count
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
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(1)
  integer(HSIZE_T) :: offset(1), length(1), stride(1)
  integer(HSIZE_T) :: num_data_values
  PetscInt :: i, temp_int
  PetscMPIInt :: rank_mpi
  PetscBool :: attribute_exists
  PetscBool :: read_times
  character(len=MAXWORDLENGTH) :: attribute_name, dataset_name

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
   
  ! read in attributes if they exist
  attribute_name = "Discretization"
  call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
  if (attribute_exists) then
    attribute_dim(1) = 3
    allocate(dataset%discretization(3))
    call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
    call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,dataset%discretization, &
                   attribute_dim,hdf5_err)
    call h5aclose_f(attribute_id,hdf5_err)
  endif
  attribute_name = "Origin"
  call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
  if (attribute_exists) then
    attribute_dim(1) = 3
    allocate(dataset%origin(3))
    call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
    call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,dataset%origin, &
                   attribute_dim,hdf5_err)
    call h5aclose_f(attribute_id,hdf5_err)
  endif  
  attribute_name = "Transient"
  call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
  if (attribute_exists) then
    read_times = PETSC_TRUE
  else
    read_times = PETSC_FALSE
  endif  

  if (read_times) then
    ! open the "time" dataset, if it exists
    dataset_name = 'time'
    call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
    call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
    call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
    allocate(dataset%time_array(num_data_values))
    dataset%time_array = 0.d0
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)
    rank_mpi = 1
    offset(1) = 0
    length(1) = num_data_values
    stride(1) = 1
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
    call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err,length)    
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,dataset%time_array,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)  
    call h5pclose_f(prop_id,hdf5_err)
    if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
    call h5sclose_f(file_space_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)    
  endif

  ! open the "data" dataset
  dataset_name = 'data'
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! get dataset dimensions
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  allocate(dims_h5(ndims_hdf5))
  allocate(max_dims_h5(ndims_hdf5))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  allocate(dataset%dims(dataset%ndims))
  dataset%ndims = ndims_hdf5
  dataset%dims = dims_h5
  deallocate(dims_h5)
  deallocate(max_dims_h5) 
  
  call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
  
  temp_int = dataset%dims(1)
  do i = 2, dataset%ndims
    temp_int = temp_int * dataset%dims(i)
  enddo
  if (num_data_values /= temp_int) then
    option%io_buffer = 'Number of values in dataset does not match dimensions.'
    call printErrMsg(option)
  endif
  if (read_times .and. dataset%dims(dataset%ndims) /= size(dataset%time_array)) then
    option%io_buffer = 'Number of times does not match last dimension of data array.'
    call printErrMsg(option)
  endif
  
  allocate(dataset%rarray(num_data_values))
  dataset%rarray = 0.d0
  call PetscLogEventBegin(logging%event_h5dread_f,ierr)
  rank_mpi = 1
  offset(1) = 0
  length(1) = num_data_values
  stride(1) = 1
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5screate_simple_f(rank_mpi,length,memory_space_id,hdf5_err,length)    
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,dataset%rarray,length, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
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


#endif ! defined(PETSC_HAVE_HDF5)

end module HDF5_aux_module
