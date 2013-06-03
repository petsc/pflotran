module Dataset_Global_class
 
  use Dataset_Base_class
  
  implicit none

  private

#include "definitions.h"

  type, public, extends(dataset_base_type) :: dataset_global_type
    character(len=MAXWORDLENGTH) :: dataset_name
    character(len=MAXSTRINGLENGTH) :: filename
!  contains
!    procedure, public :: Init => DatasetGlobalInit
!    procedure, public :: Load => DatasetGlobalLoad
  end type dataset_global_type
  
  public :: DatasetGlobalCreate, &
            DatasetGlobalInit, &
            DatasetGlobalLoad, &
            DatasetGlobalDestroy
  
#if defined(PETSC_HAVE_HDF5)   
  public :: DatasetGlobalReadTimes
#endif  
  
contains

! ************************************************************************** !
!
! DatasetGlobalCreate: Creates global dataset class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
function DatasetGlobalCreate()
  
  implicit none
  
  class(dataset_global_type), pointer :: dataset

  class(dataset_global_type), pointer :: DatasetGlobalCreate
  
  allocate(dataset)
  call DatasetGlobalInit(dataset)

  DatasetGlobalCreate => dataset
    
end function DatasetGlobalCreate

! ************************************************************************** !
!
! DatasetGlobalInit: Initializes members of global dataset class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetGlobalInit(this)
  
  implicit none
  
  class(dataset_global_type) :: this
  
  call DatasetBaseInit(this)
  this%filename = ''
  this%dataset_name = ''
    
end subroutine DatasetGlobalInit

! ************************************************************************** !
!
! DatasetGlobalLoad: Load new data into dataset buffer
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetGlobalLoad(this,discretization,grid,option)
  
#if defined(PETSC_HAVE_HDF5)    
  use hdf5, only : H5T_NATIVE_DOUBLE
#endif
  use Discretization_module
  use Grid_module
  use Option_module
  use Time_Storage_module

  implicit none
  
  class(dataset_global_type) :: this
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  
  if (.not.associated(this%time_storage)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetGlobalReadTimes(this%filename,this%dataset_name, &
                                this%time_storage,option)
#endif
  endif
  
  this%time_storage%cur_time = option%time
  ! sets correct cur_time_index
  call TimeStorageUpdate(this%time_storage)
  
  if (.not.associated(this%rarray)) then
    allocate(this%rarray(grid%nlmax))
    this%rarray = 0.d0
  endif
  
  if (this%time_storage%cur_time_index > 0 .and. &
      this%time_storage%cur_time_index >= &
        ! both of the below will be zero initially
        this%buffer_slice_offset + this%buffer_nslice &
      .and. &
      (this%time_storage%cur_time_index < &
        this%time_storage%max_time_index .or. &
       ! essentially gets the data set read if only one time slice
       .not.associated(this%rbuffer))) then
    if (.not.associated(this%rbuffer)) then ! not initialized
      this%buffer_nslice = min(2,this%time_storage%max_time_index)
      allocate(this%rbuffer(grid%nlmax*this%buffer_nslice))
      this%rbuffer = 0.d0
    endif
    this%buffer_slice_offset = this%time_storage%cur_time_index - 1
#if defined(PETSC_HAVE_HDF5)    
    call DatasetGlobalReadData(discretization,grid,option,this%filename, &
                               this%dataset_name, &
                               this%buffer_slice_offset, &
                               this%rbuffer, &
                               H5T_NATIVE_DOUBLE)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
!  call this%InterpolateTime()
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetGlobalLoad

#if defined(PETSC_HAVE_HDF5)
! ************************************************************************** !
!
! DatasetGlobalReadTimes: Read dataset times into time storage
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine DatasetGlobalReadTimes(filename,dataset_name,time_storage,option)
                         
  use hdf5
  use Time_Storage_module
  use Units_module, only : UnitsConvertToInternal
  use Option_module
  use Logging_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: dataset_name
  type(time_storage_type), pointer :: time_storage
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: atype_id  
  integer(HID_T) :: attribute_id  
  integer(HSIZE_T) :: num_times  
  integer(HSIZE_T) :: length(1)
  integer(HSIZE_T) :: attribute_dim(1)
  integer(SIZE_T) :: size_t_int  
  PetscMPIInt :: array_rank_mpi
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: attribute_name, time_units
  PetscMPIInt :: int_mpi
  PetscInt :: temp_int
  PetscMPIInt :: hdf5_err
  PetscBool :: attribute_exists
  PetscErrorCode :: ierr
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr)
  
  time_storage => TimeStorageCreate()

  if (option%myrank == option%io_rank) then
    ! open the file
    call h5open_f(hdf5_err)
    option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
    call printMsg(option)
  
    ! set read file access property
    call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
    call h5pclose_f(prop_id,hdf5_err)

    option%io_buffer = 'Opening hdf5 group: ' // trim(dataset_name)
    call printMsg(option)
    call h5gopen_f(file_id,dataset_name,grp_id,hdf5_err)

    attribute_name = "Time Units"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim = 1
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,hdf5_err)
      size_t_int = MAXWORDLENGTH
      call h5tset_size_f(atype_id,size_t_int,hdf5_err)
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,atype_id,time_units,attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    else
      option%io_buffer = 'Time Units assumed to be seconds.'
      call printWrnMsg(option)
      time_units = 's'
    endif
  
    string = 'Times'
    option%io_buffer = 'Opening data set: ' // trim(string)
    call printMsg(option)  
    call h5dopen_f(grp_id,string,dataset_id,hdf5_err)
    call h5dget_space_f(dataset_id,file_space_id,hdf5_err)
    call h5sget_simple_extent_npoints_f(file_space_id,num_times,hdf5_err)
    temp_int = num_times
  endif
  
  int_mpi = 1
  call MPI_Bcast(temp_int,int_mpi,MPI_INTEGER,option%io_rank, &
                 option%mycomm,ierr)
  num_times = temp_int
  time_storage%max_time_index = num_times
  allocate(time_storage%times(num_times))
  time_storage%times = 0.d0
  
  if (option%myrank == option%io_rank) then 
    array_rank_mpi = 1
    length(1) = num_times
    call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif  
    call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,time_storage%times, &
                    length,hdf5_err,memory_space_id,file_space_id,prop_id)

      
    call h5pclose_f(prop_id,hdf5_err)
    if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
    call h5sclose_f(file_space_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)  
    option%io_buffer = 'Closing group: ' // trim(dataset_name)
    call printMsg(option)  
    call h5gclose_f(grp_id,hdf5_err)  
    option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
    call printMsg(option)  
    call h5fclose_f(file_id,hdf5_err)
    call h5close_f(hdf5_err) 
    time_storage%times = time_storage%times * &
      UnitsConvertToInternal(time_units,option)
  endif

  int_mpi = num_times
  call MPI_Bcast(time_storage%times,int_mpi,MPI_DOUBLE_PRECISION, &
                 option%io_rank,option%mycomm,ierr)

  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr)
  
end subroutine DatasetGlobalReadTimes

! ************************************************************************** !
!
! DatasetGlobalReadData: Read an hdf5 array into a Petsc Vec
! author: Glenn Hammond
! date: 01/12/08
!
! ************************************************************************** !
subroutine DatasetGlobalReadData(discretization,grid,option,filename, &
                                 dataset_name,requested_rank2_offset,buffer, &
                                 data_type)
                         
  use hdf5
  use Logging_module
  use Option_module
  use Grid_module
  use Discretization_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

! Default HDF5 Mechanism 
 
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscInt :: requested_rank2_offset
  PetscReal, pointer :: buffer(:)
  integer(HID_T) :: data_type 
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3), max_dims(3)
  integer(HSIZE_T) :: offset(3), length(3), stride(3)
  PetscMPIInt :: ndims
  PetscMPIInt :: rank_mpi
  Vec :: natural_vec
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, istart
  PetscInt :: buffer_size, buffer_rank2_size, file_rank1_size, file_rank2_size
  integer, allocatable :: integer_buffer_i4(:)
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err
  
  call PetscLogEventBegin(logging%event_read_array_hdf5,ierr)

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  string = trim(dataset_name) // '/Data'
  option%io_buffer = 'Opening data set: ' // trim(string)
  call printMsg(option)  
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  
  call h5sget_simple_extent_ndims_f(file_space_id,ndims,hdf5_err)
  call h5sget_simple_extent_dims_f(file_space_id,dims,max_dims,hdf5_err)
  
  buffer_size = size(buffer)
  buffer_rank2_size = buffer_size / grid%nlmax
  file_rank1_size = dims(1)
  if (ndims > 1) then
    file_rank2_size = dims(2)
  else
    file_rank2_size = 1
  endif

  if (mod(buffer_size,grid%nlmax) /= 0) then
    write(option%io_buffer, &
          '(a," buffer dimension (",i9,") is not a multiple of domain",&
           &" dimension (",i9,").")') trim(dataset_name), &
           size(buffer,1), grid%nlmax
    call printErrMsg(option)   
  endif

  if (mod(file_rank1_size,grid%nlmax) /= 0) then
    write(option%io_buffer, &
          '(a," data space dimension (",i9,") is not a multiple of domain",&
           &" dimension (",i9,").")') trim(dataset_name), &
           file_rank1_size, grid%nlmax
    call printErrMsg(option)   
  endif

  call MPI_Exscan(grid%nlmax,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr)
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif

  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)
  call VecZeroEntries(natural_vec,ierr)
  call VecZeroEntries(global_vec,ierr)

  ! must initialize here to avoid error below when closing memory space
  memory_space_id = -1
  
  ! Find sizes of array.  Must use identical data types
  

  ! offset is zero-based
  offset = 0
  length = 0
  stride = 1
  
  offset(1) = istart ! istart is offset in the first dimension
  length(1) = file_rank1_size
  if (ndims > 1) then
    offset(2) = requested_rank2_offset
    length(1) = file_rank1_size
    length(2) = min(buffer_rank2_size,file_rank2_size-requested_rank2_offset)
  else
    offset(1) = requested_rank2_offset*grid%nlmax
    length(1) = min(buffer_size, &
                    file_rank1_size-requested_rank2_offset*grid%nlmax)
  endif
  call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset, &
                             length,hdf5_err,stride,stride) 

  dims = 0
  rank_mpi = 1
  if (ndims > 1) then
    dims(1) = min(buffer_size, &
                  grid%nlmax*(file_rank2_size-requested_rank2_offset))
  else
    dims(1) = min(buffer_size, &
                  file_rank1_size-requested_rank2_offset*grid%nlmax)
  endif
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  ! initialize to -999 to catch errors
  buffer = -999.d0
  
  if (data_type == H5T_NATIVE_DOUBLE) then
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
    call h5dread_f(data_set_id,H5T_NATIVE_DOUBLE,buffer,dims, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
  else if (data_type == H5T_NATIVE_INTEGER) then
    allocate(integer_buffer_i4(size(buffer)))
    call PetscLogEventBegin(logging%event_h5dread_f,ierr)                              
    call h5dread_f(data_set_id,H5T_NATIVE_INTEGER,integer_buffer_i4,dims, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
    call PetscLogEventEnd(logging%event_h5dread_f,ierr)                              
    do i = 1, min(buffer_size, &
                  grid%nlmax*(file_rank2_size-requested_rank2_offset))
      buffer(i) = real(integer_buffer_i4(i))
    enddo
    deallocate(integer_buffer_i4)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  string = trim(dataset_name) // '/Data'
  option%io_buffer = 'Closing data set: ' // trim(string)
  call printMsg(option)  
  call h5dclose_f(data_set_id,hdf5_err)
  option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err) 
  
  istart = 0
  do i = 1, min(buffer_rank2_size,file_rank2_size-requested_rank2_offset)
    call VecGetArrayF90(natural_vec,vec_ptr,ierr)
    vec_ptr(:) = buffer(istart+1:istart+grid%nlmax)
    call VecRestoreArrayF90(natural_vec,vec_ptr,ierr)
    call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                       global_vec,ONEDOF)
    call VecGetArrayF90(global_vec,vec_ptr,ierr)
    buffer(istart+1:istart+grid%nlmax) = vec_ptr(:)
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr)
    istart = istart + grid%nlmax
  enddo
  call VecDestroy(natural_vec,ierr)
  call VecDestroy(global_vec,ierr)
  
  call PetscLogEventEnd(logging%event_read_array_hdf5,ierr)
! End of Default HDF5 Mechanism

end subroutine DatasetGlobalReadData

#endif

! ************************************************************************** !
!
! DatasetGlobalDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine DatasetGlobalDestroy(this)

  implicit none
  
  class(dataset_global_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetBaseStrip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetGlobalDestroy

end module Dataset_Global_class
