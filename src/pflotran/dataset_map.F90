module Dataset_Map_class
 
  use Dataset_Common_HDF5_class
  
  implicit none

  private

#include "definitions.h"

  type, public, extends(dataset_common_hdf5_type) :: dataset_map_type
    character(len=MAXSTRINGLENGTH) :: h5_dataset_map_name
    character(len=MAXSTRINGLENGTH) :: map_filename
    PetscInt,pointer :: map(:,:)
    PetscInt         :: map_dims_global(2)
    PetscInt         :: map_dims_local(2)
    PetscInt,pointer :: datatocell_ids(:)
    PetscInt,pointer :: cell_ids_local(:)
    PetscBool        :: first_time
!  contains
!    procedure, public :: Init => DatasetMapInit
!    procedure, public :: Load => DatasetMapLoad
  end type dataset_map_type
  
  public :: DatasetMapCreate, &
            DatasetMapInit, &
            DatasetMapLoad, &
            DatasetMapDestroy
  
contains

! ************************************************************************** !
!
! DatasetMapCreate: Creates global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
function DatasetMapCreate()
  
  implicit none
  
  class(dataset_map_type), pointer :: dataset

  class(dataset_map_type), pointer :: DatasetMapCreate
  
  allocate(dataset)
  call DatasetMapInit(dataset)

  DatasetMapCreate => dataset
    
end function DatasetMapCreate

! ************************************************************************** !
!
! DatasetMapInit: Initializes members of global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetMapInit(this)
  
  implicit none
  
  class(dataset_map_type) :: this
  
  call DatasetCommonHDF5Init(this)
  this%h5_dataset_map_name = ''
  this%map_filename = ''
  nullify(this%map)
  this%map_dims_global = 0
  this%map_dims_local = 0
  nullify(this%datatocell_ids)
  nullify(this%cell_ids_local)
  this%first_time = PETSC_TRUE
    
end subroutine DatasetMapInit

! ************************************************************************** !
!
! DatasetMapLoad: Load new data into dataset buffer
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetMapLoad(this,option)
  
  use Option_module
  use Time_Storage_module
  use Dataset_Base_class

  implicit none
  
  class(dataset_map_type) :: this
  type(option_type) :: option
  
  if (DatasetCommonHDF5Load(this,option)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetMapReadData(this,option)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetMapLoad

#if defined(PETSC_HAVE_HDF5)    
! ************************************************************************** !
!
! DatasetMapReadData: Read an hdf5 array 
! author: Glenn Hammond
! date: 10/25/11, 05/29/13
!
! ************************************************************************** !
subroutine DatasetMapReadData(this,option)

  use hdf5
  use Option_module
  use Units_module
  use Logging_module
  
  implicit none
  
  class(dataset_map_type) :: this
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
  
! add dataset reading routine here.  
  
  option%io_buffer = 'Closing group: ' // trim(this%hdf5_dataset_name)
  call printMsg(option)  
  call h5gclose_f(grp_id,hdf5_err)  
  option%io_buffer = 'Closing hdf5 file: ' // trim(this%filename)
  call printMsg(option)  
  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)  
  
end subroutine DatasetMapReadData
#endif

! ************************************************************************** !
!
! DatasetMapStrip: Strips allocated objects within Map dataset object
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetMapStrip(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_map_type)  :: this
  
  call DatasetCommonHDF5Strip(this)
  
  call DeallocateArray(this%map)
  call DeallocateArray(this%datatocell_ids)
  call DeallocateArray(this%cell_ids_local)
  
end subroutine DatasetMapStrip

! ************************************************************************** !
!
! DatasetMapDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetMapDestroy(this)

  implicit none
  
  class(dataset_map_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetMapStrip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetMapDestroy

end module Dataset_Map_class
