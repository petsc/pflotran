module Dataset_XYZ_class
 
  use Dataset_Common_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public, extends(dataset_common_hdf5_type) :: dataset_xyz_type
    PetscBool :: is_cell_centered
    PetscInt :: data_dim
    PetscReal, pointer :: origin(:)
    PetscReal, pointer :: discretization(:)
!  contains
!    procedure, public :: Init => DatasetXYZInit
!    procedure, public :: Load => DatasetXYZLoad
  end type dataset_xyz_type
  
  PetscInt, parameter, public :: DIM_NULL = 0
  PetscInt, parameter, public :: DIM_X = 1
  PetscInt, parameter, public :: DIM_Y = 2
  PetscInt, parameter, public :: DIM_Z = 3
  PetscInt, parameter, public :: DIM_XY = 4
  PetscInt, parameter, public :: DIM_XZ = 5
  PetscInt, parameter, public :: DIM_YZ = 6
  PetscInt, parameter, public :: DIM_XYZ = 7
  
  PetscInt, parameter :: MAX_NSLICE = 4

  public :: DatasetXYZCreate, &
            DatasetXYZInit, &
            DatasetXYZCast, &
            DatasetXYZLoad, &
            DatasetXYZInterpolateReal, &
            DatasetXYZStrip, &
            DatasetXYZDestroy
  
contains

! ************************************************************************** !
!
! DatasetXYZCreate: Creates global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
function DatasetXYZCreate()
  
  implicit none
  
  class(dataset_xyz_type), pointer :: dataset

  class(dataset_xyz_type), pointer :: DatasetXYZCreate
  
  allocate(dataset)
  call DatasetXYZInit(dataset)

  DatasetXYZCreate => dataset
    
end function DatasetXYZCreate

! ************************************************************************** !
!
! DatasetXYZInit: Initializes members of global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZInit(this)
  
  implicit none
  
  class(dataset_xyz_type) :: this
  
  call DatasetCommonHDF5Init(this)
  this%is_cell_centered = PETSC_FALSE
  this%data_dim = DIM_NULL
  nullify(this%origin)
  nullify(this%discretization)
    
end subroutine DatasetXYZInit

! ************************************************************************** !
!
! DatasetXYZCast: Casts a dataset_base_type to dataset_xyz_type
! author: Glenn Hammond
! date: 08/29/13
!
! ************************************************************************** !
function DatasetXYZCast(this)
  
  use Dataset_Base_class

  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_xyz_type), pointer :: DatasetXYZCast
  
  nullify(DatasetXYZCast)
  select type (this)
    class is (dataset_xyz_type)
      DatasetXYZCast => this
  end select
    
end function DatasetXYZCast

! ************************************************************************** !
!
! DatasetXYZLoad: Load new data into dataset buffer
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZLoad(this,option)
  
  use Option_module
  use Time_Storage_module
  use Dataset_Base_class  

  implicit none
  
  class(dataset_xyz_type) :: this
  type(option_type) :: option
  
  if (DatasetCommonHDF5Load(this,option)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetXYZReadData(this,option)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetXYZLoad

#if defined(PETSC_HAVE_HDF5)    
! ************************************************************************** !
!
! DatasetXYZReadData: Read an hdf5 data into arrays
! author: Glenn Hammond
! date: 10/25/11, 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZReadData(this,option)

  use hdf5
  use Option_module
  use Units_module
  use Logging_module
  
  implicit none
  
  class(dataset_xyz_type) :: this
  type(option_type) :: option
  
  integer(HID_T) :: file_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: prop_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: attribute_id
  integer(HID_T) :: ndims_h5
  integer(HID_T) :: atype_id
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: attribute_dim(3)
  integer(HSIZE_T) :: offset(4), length(4), stride(4)
  integer(HSIZE_T) :: num_data_values
  integer(SIZE_T) size_t_int
  PetscInt :: i, temp_int, temp_array(4)
  PetscInt :: num_spatial_dims, time_dim, num_times
  PetscInt :: num_dims_in_h5_file, num_times_in_h5_file
  PetscMPIInt :: array_rank_mpi
  PetscBool :: attribute_exists
  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: attribute_name, dataset_name, word

  !TODO(geh): add to event log
  !call PetscLogEventBegin(logging%event_read_datset_hdf5,ierr)

  ! open the file
  call h5open_f(hdf5_err)
  option%io_buffer = 'Opening hdf5 file: ' // trim(this%filename)
  call printMsg(option)
  
  ! set read file access property
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  call h5fopen_f(this%filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  ! the dataset is actually stored in a group.  the group contains
  ! a "data" dataset and optionally a "time" dataset.
  option%io_buffer = 'Opening group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gopen_f(file_id,this%hdf5_dataset_name,grp_id,hdf5_err)

  ! only want to read on first time through
  if (this%data_dim == DIM_NULL) then
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
      call DatasetSetDimension(this,word)
    else
      option%io_buffer = &
        'Dimension attribute must be included in hdf5 dataset file.'
      call printErrMsg(option)
    endif
    attribute_name = "Discretization"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGetNDimensions(this)
      allocate(this%discretization(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,this%discretization, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    else
      option%io_buffer = &
        '"Discretization" attribute must be included in XYZ hdf5 dataset file.'
    endif
    attribute_name = "Origin"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      attribute_dim(1) = DatasetGetNDimensions(this)
      allocate(this%origin(attribute_dim(1)))
      call h5aopen_f(grp_id,attribute_name,attribute_id,hdf5_err)
      call h5aread_f(attribute_id,H5T_NATIVE_DOUBLE,this%origin, &
                     attribute_dim,hdf5_err)
      call h5aclose_f(attribute_id,hdf5_err)
    endif
    attribute_name = "Cell Centered"
    call H5aexists_f(grp_id,attribute_name,attribute_exists,hdf5_err)
    if (attribute_exists) then
      this%is_cell_centered = PETSC_TRUE
    endif
  endif ! this%data_dim == DIM_NULL
  
  num_spatial_dims = DatasetGetNDimensions(this)
  time_dim = -1
  num_times = 1
  if (associated(this%time_storage)) then
    num_times = this%time_storage%max_time_index
    time_dim = num_spatial_dims + 1
  endif
  
  ! open the "data" dataset
  dataset_name = 'Data'
  if (this%realization_dependent) then
    write(word,'(i9)') option%id
    dataset_name = trim(dataset_name) // trim(adjustl(word))
  endif
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! get dataset dimensions
  if (.not.associated(this%dims)) then
    call h5sget_simple_extent_ndims_f(file_space_id,ndims_h5,hdf5_err)
    allocate(dims_h5(ndims_h5))
    allocate(max_dims_h5(ndims_h5))
    call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
    if (associated(this%time_storage)) then
      ! if dataset is time dependent, need to remove that time index from 
      ! dimensions and decrement rank (we don't want to include time)
      this%rank = ndims_h5-1
      ! the first dimension of dims_h5 is the time dimension
      num_times_in_h5_file = dims_h5(1)
    else
      this%rank = ndims_h5
      num_times_in_h5_file = 0
    endif
    allocate(this%dims(this%rank))
    ! have to invert dimensions, but do not count the time dimension
    do i = 1, this%rank
      this%dims(i) = int(dims_h5(ndims_h5-i+1))
    enddo
    deallocate(dims_h5)
    deallocate(max_dims_h5) 
  
    call h5sget_simple_extent_npoints_f(file_space_id,num_data_values,hdf5_err)
  
    temp_int = this%dims(1)
    do i = 2, num_spatial_dims
      temp_int = temp_int * this%dims(i)
    enddo
    if (num_data_values/num_times /= temp_int) then
      option%io_buffer = &
        'Number of values in dataset does not match dimensions.'
      call printErrMsg(option)
    endif
    if (associated(this%time_storage) .and. &
        num_times_in_h5_file /= num_times) then
      option%io_buffer = &
        'Number of times does not match last dimension of data array.'
      call printErrMsg(option)
    endif
    if (.not.associated(this%rarray)) then
      allocate(this%rarray(temp_int))
    endif
    this%rarray = 0.d0
    if (associated(this%time_storage) .and. .not.associated(this%rbuffer)) then
      ! buffered array
      allocate(this%rbuffer(size(this%rarray)*MAX_NSLICE))
      this%rbuffer = 0.d0
    endif
  endif

  call PetscLogEventBegin(logging%event_h5dread_f,ierr)

  if (associated(this%time_storage)) then
    num_dims_in_h5_file = this%rank + 1
  else
    num_dims_in_h5_file = this%rank
  endif
  
  array_rank_mpi = 1
  length = 1
  ! length (or size) must be adjusted according to the size of the 
  ! remaining data in the file
  this%buffer_nslice = min(MAX_NSLICE,(num_times-this%buffer_slice_offset))
  if (time_dim > 0) then
    length(1) = size(this%rarray) * this%buffer_nslice
  else
    length(1) = size(this%rarray)
  endif
  call h5screate_simple_f(array_rank_mpi,length,memory_space_id,hdf5_err,length)    

  length = 1
  stride = 1
  offset = 0
  do i = 1, num_spatial_dims
    length(i) = this%dims(i)
  enddo
  ! cannot read beyond end of the buffer
  if (time_dim > 0) then
    length(time_dim) = this%buffer_nslice
    offset(time_dim) = this%buffer_slice_offset
  endif
  
  
  !geh: for some reason, we have to invert here.  Perhaps because the
  !     dataset was generated in C???
  temp_array(1:num_dims_in_h5_file) = int(length(1:num_dims_in_h5_file))
  do i = 1, num_dims_in_h5_file
    length(i) = temp_array(num_dims_in_h5_file-i+1)
  enddo
  temp_array(1:num_dims_in_h5_file) = int(offset(1:num_dims_in_h5_file))
  do i = 1, num_dims_in_h5_file
    offset(i) = temp_array(num_dims_in_h5_file-i+1)
  enddo
  ! stride is fine
  
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F,hdf5_err)
#endif
  call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F,offset,length, &
                             hdf5_err,stride,stride) 
  if (associated(this%rbuffer)) then
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%rbuffer,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
  else
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,this%rarray,length, &
                   hdf5_err,memory_space_id,file_space_id,prop_id)
!    this%rmax = maxval(this%rarray)
!    this%rmin = minval(this%rarray)
  endif
  
  call h5pclose_f(prop_id,hdf5_err)
  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  

  call PetscLogEventEnd(logging%event_h5dread_f,ierr) 

  option%io_buffer = 'Closing group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)
  
  !TODO(geh): add to event log
  !call PetscLogEventEnd(logging%event_read_ndim_real_array_hdf5,ierr)
                          
end subroutine DatasetXYZReadData

#endif

! ************************************************************************** !
!
! DatasetSetDimension: Sets the dimension of the dataset
! author: Glenn Hammond
! date: 10/24/11, 05/29/13
!
! ************************************************************************** !
subroutine DatasetSetDimension(this,word)

  use String_module

  implicit none
  
  class(dataset_xyz_type) :: this
  character(len=MAXWORDLENGTH) :: word

  call StringToUpper(word)
  select case(word)
    case('X')
      this%data_dim = DIM_X
    case('Y')
      this%data_dim = DIM_Y
    case('Z')
      this%data_dim = DIM_Z
    case('XY')
      this%data_dim = DIM_XY
    case('XZ')
      this%data_dim = DIM_XZ
    case('YZ')
      this%data_dim = DIM_YZ
    case('XYZ')
      this%data_dim = DIM_XYZ
  end select
      
end subroutine DatasetSetDimension

! ************************************************************************** !
!
! DatasetGetNDimensions: Returns the number of dimensions
! author: Glenn Hammond
! date: 10/24/11, 05/29/13
!
! ************************************************************************** !
function DatasetGetNDimensions(this)

  implicit none
  
  class(dataset_xyz_type) :: this
  
  PetscInt :: DatasetGetNDimensions

  select case(this%data_dim)
    case(DIM_X,DIM_Y,DIM_Z)
      DatasetGetNDimensions = ONE_INTEGER
    case(DIM_XY,DIM_XZ,DIM_YZ)
      DatasetGetNDimensions = TWO_INTEGER
    case(DIM_XYZ)
      DatasetGetNDimensions = THREE_INTEGER
    case default
      DatasetGetNDimensions = ZERO_INTEGER
  end select
      
end function DatasetGetNDimensions

! ************************************************************************** !
!
! DatasetXYZInterpolateReal: Interpolates data from the dataset
! author: Glenn Hammond
! date: 10/26/11, 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZInterpolateReal(this,xx,yy,zz,time,real_value,option)

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  class(dataset_xyz_type) :: this
  PetscReal, intent(in) :: xx, yy, zz
  PetscReal :: time
  PetscReal :: real_value
  type(option_type) :: option
  
  PetscInt :: spatial_interpolation_method
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  PetscReal :: x1, x2, y1, y2
  PetscReal :: v1, v2, v3, v4
  PetscInt :: index
  PetscReal :: dx, dy
  PetscInt :: nx
  character(len=MAXWORDLENGTH) :: word
  
  call DatasetXYZGetIndices(this,xx,yy,zz,i,j,k,x,y,z)
  
  spatial_interpolation_method = INTERPOLATION_LINEAR
  
  ! in the below, i,j,k,xx,yy,zz to not reflect the 
  ! coordinates of the problem domain in 3D.  They
  ! are transfored to the dimensions of the dataset
  select case(spatial_interpolation_method)
    case(INTERPOLATION_STEP)
    case(INTERPOLATION_LINEAR)
      select case(this%data_dim)
        case(DIM_X,DIM_Y,DIM_Z)
          if (i < 1 .or. i+1 > this%dims(1)) then 
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_X)
                option%io_buffer = 'Out of x bounds, i = ' // trim(word)
              case(DIM_Y)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
              case(DIM_Z)
                option%io_buffer = 'Out of z bounds, k = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          dx = this%discretization(1)
          x1 = this%origin(1) + (i-1)*dx
          if (this%is_cell_centered) x1 = x1 + 0.5d0*dx
          v1 = this%rarray(i)
          v2 = this%rarray(i+1)
          real_value = v1 + (x-x1)/dx*(v2-v1)
        case(DIM_XY,DIM_XZ,DIM_YZ)
          if (i < 1 .or. i+1 > this%dims(1)) then
            write(word,*) i
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY,DIM_XZ)
                option%io_buffer = 'Out of x bounds, i = ' // trim(word)
              case(DIM_YZ)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          if (j < 1 .or. j+1 > this%dims(2)) then
            write(word,*) j
            word = adjustl(word)
            select case(this%data_dim)
              case(DIM_XY)
                option%io_buffer = 'Out of y bounds, j = ' // trim(word)
              case(DIM_YZ,DIM_XZ)
                option%io_buffer = 'Out of z bounds, k = ' // trim(word)
            end select
            call printErrMsgByRank(option)
          endif
          dx = this%discretization(1)
          dy = this%discretization(2)
          nx = this%dims(1)

          x1 = this%origin(1) + (i-1)*dx
          if (this%is_cell_centered) x1 = x1 + 0.5d0*dx
          x2 = x1 + dx
          
          index = i + (j-1)*nx
          v1 = this%rarray(index)
          v2 = this%rarray(index+1)
          
          y1 = this%origin(2) + (j-1)*dy
          if (this%is_cell_centered) y1 = y1 + 0.5d0*dy
          y2 = y1 + dy
          
           ! really (j1-1+1)
          index = i + j*nx
          v3 = this%rarray(index)
          v4 = this%rarray(index+1)
          
          real_value = InterpolateBilinear(x,y,x1,x2,y1,y2,v1,v2,v3,v4)
        case(DIM_XYZ)
          option%io_buffer = 'Trilinear interpolation not yet supported'
          call printErrMsgByRank(option)
      end select
  end select
  
end subroutine DatasetXYZInterpolateReal

! ************************************************************************** !
!
! DatasetXYZGetIndices: Returns bounding indices for point in dataset
! author: Glenn Hammond
! date: 10/26/11, 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZGetIndices(this,xx,yy,zz,i,j,k,x,y,z)

  implicit none

  class(dataset_xyz_type) :: this
  PetscReal, intent(in)  :: xx, yy, zz
  PetscInt :: i, j, k
  PetscReal :: x, y, z
  
  select case(this%data_dim)
    ! since these are 1D array, always use first dimension
    case(DIM_X)
      x = xx
    case(DIM_Y)
      x = yy
    case(DIM_XY)
      x = xx
      y = yy
    case(DIM_XYZ)
      x = xx
      y = yy
      z = zz
    case(DIM_Z)
      x = zz
    case(DIM_XZ)
      x = xx
      y = zz
    case(DIM_YZ)
      x = yy
      y = zz
  end select

  if (this%is_cell_centered) then
    i = int((x - this%origin(1))/ &
            this%discretization(1) + 0.5d0)
    i = max(1,min(i,this%dims(1)-1))
    if (this%data_dim > DIM_Z) then ! at least 2D
      j = int((y - this%origin(2))/ &
              this%discretization(2) + 0.5d0)
      j = max(1,min(j,this%dims(2)-1))
    endif
    if (this%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - this%origin(3))/ &
              this%discretization(3) + 0.5d0)
      k = max(1,min(k,this%dims(3)-1))
    endif
  else
    i = int((x - this%origin(1))/ &
            this%discretization(1) + 1.d0)
    if (this%data_dim > DIM_Z) then ! at least 2D
      j = int((y - this%origin(2))/ &
              this%discretization(2) + 1.d0)
    endif
    if (this%data_dim > DIM_YZ) then ! at least 3D
      k = int((z - this%origin(3))/ &
              this%discretization(3) + 1.d0)
    endif
  endif
  
end subroutine DatasetXYZGetIndices

! ************************************************************************** !
!
! DatasetXYZStrip: Strips allocated objects within XYZ dataset object
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetXYZStrip(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_xyz_type)  :: this
  
  call DatasetCommonHDF5Strip(this)
  
  call DeallocateArray(this%origin)
  call DeallocateArray(this%discretization)
  
end subroutine DatasetXYZStrip

! ************************************************************************** !
!
! DatasetXYZDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetXYZDestroy(this)

  implicit none
  
  class(dataset_xyz_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetXYZStrip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetXYZDestroy

end module Dataset_XYZ_class
