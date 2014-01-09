module HDF5_Aux_module

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif
  use Logging_module

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  private

  PetscInt, parameter, public :: HDF5_READ_BUFFER_SIZE = 1000000
!#define HDF5_BROADCAST

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
#ifdef SCORPIO
            HDF5ReadDatasetInteger2D, &
            HDF5ReadDatasetReal2D, &
            HDF5GroupExists, &
#else
            HDF5GroupExists, &
#endif
! SCORPIO
            HDF5MakeStringCompatible

contains

#endif
  

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine HDF5ReadNDimRealArray(option,file_id,dataset_name,ndims,dims, &
                                 real_array)
  ! 
  ! Read in an n-dimensional array from an hdf5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/10
  ! 

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
  dims = int(dims_h5)
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

#if defined(SCORPIO)

! ************************************************************************** !

subroutine HDF5ReadDatasetInteger2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)
  ! 
  ! Author: Gautam Bisht
  ! Date: 05/13/2010
  ! 

  use hdf5
  use Option_module
  
  implicit none
  
#if defined(SCORPIO)
  include "scorpiof.h"  
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
  call scorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call scorpio_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call scorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call scorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call scorpio_read_dataset( data, SCORPIO_INTEGER, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)
  
  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call scorpio_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetInteger2D
#endif
! SCORPIO

#if defined(SCORPIO)

! ************************************************************************** !

subroutine HDF5ReadDatasetReal2D(filename,dataset_name,read_option,option, &
           data,data_dims,dataset_dims)
  use hdf5
  use Option_module
  
  implicit none
  
#if defined(SCORPIO)
  include "scorpiof.h"  
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
  call scorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)

  ! Get dataset dimnesions
  call scorpio_get_dataset_ndims(ndims, file_id, dataset_name, option%ioread_group_id, ierr)
  if (ndims > 2) then
    option%io_buffer='Dimension of ' // dataset_name // ' dataset in ' // filename // &
    ' is greater than to 2.'
    call printErrMsg(option)
  endif
  
  ! Get size of each dimension
  call scorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)
  
  data_dims(1) = dataset_dims(1)/option%mycommsize
  data_dims(2) = dataset_dims(2)

  remainder = dataset_dims(1) - data_dims(1)*option%mycommsize
  if (option%myrank < remainder) data_dims(1) = data_dims(1) + 1

  
  allocate(data(data_dims(2),dataset_dims(1)))
  
  call scorpio_get_dataset_dims(dataset_dims, file_id, dataset_name, option%ioread_group_id, ierr)

  ! Read the dataset collectively
  call scorpio_read_dataset( data, SCORPIO_DOUBLE, ndims, dataset_dims, data_dims, & 
            file_id, dataset_name, option%ioread_group_id, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
  
  data_dims(1) = data_dims(1) + data_dims(2)
  data_dims(2) = data_dims(1) - data_dims(2)
  data_dims(1) = data_dims(1) - data_dims(2)

  dataset_dims(1) = dataset_dims(1) + dataset_dims(2)
  dataset_dims(2) = dataset_dims(1) - dataset_dims(2)
  dataset_dims(1) = dataset_dims(1) - dataset_dims(2)

  ! Close file
  call scorpio_close_file( file_id, option%ioread_group_id, ierr)  

end subroutine HDF5ReadDatasetReal2D
#endif

! ************************************************************************** !

function HDF5GroupExists(filename,group_name,option)
  ! 
  ! SCORPIO
  ! Returns true if a group exists
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/26/2012
  ! 

  use hdf5
  use Option_module
  
  implicit none

#if defined(SCORPIO)
  include "scorpiof.h"  
#endif
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: group_name
  type(option_type) :: option

  integer(HID_T) :: file_id  
  integer(HID_T) :: grp_id  
  integer(HID_T) :: prop_id
  PetscMPIInt, parameter :: ON=1, OFF=0 
  PetscBool :: group_exists
  
  PetscBool :: HDF5GroupExists

#if defined(SCORPIO)

  ! Open file collectively
  filename = trim(filename) // CHAR(0)
  group_name = trim(group_name) // CHAR(0)
  call scorpio_open_file(filename, option%ioread_group_id, SCORPIO_FILE_READONLY, file_id, ierr)
  call scorpio_group_exists(group_name, file_id, option%ioread_group_id, ierr)
  group_exists = (ierr == 1);

  if (group_exists) then
    HDF5GroupExists = PETSC_TRUE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" found in file.'
  else
    HDF5GroupExists = PETSC_FALSE
    option%io_buffer = 'Group "' // trim(group_name) // '" in HDF5 file "' // &
      trim(filename) // '" not found in file.  Therefore, assuming a ' // &
      'cell-indexed dataset.'
  endif
  call printMsg(option)

  call scorpio_close_file( file_id, option%ioread_group_id, ierr)  
#else
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
#endif 
!SCORPIO

end function HDF5GroupExists

! ************************************************************************** !

subroutine HDF5MakeStringCompatible(name)
  ! 
  ! Replaces '/' in string with '_' for hdf5 names
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

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

end subroutine HDF5MakeStringCompatible

#endif
! defined(PETSC_HAVE_HDF5)

end module HDF5_Aux_module
