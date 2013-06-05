module Dataset_Map_class
 
  use Dataset_Common_HDF5_class
  
  implicit none

  private

#include "definitions.h"

  type, public, extends(dataset_common_hdf5_type) :: dataset_map_type
    character(len=MAXSTRINGLENGTH) :: h5_dataset_map_name
    character(len=MAXSTRINGLENGTH) :: map_filename
    PetscInt, pointer :: mapping(:,:)
    PetscInt          :: map_dims_global(2)
    PetscInt          :: map_dims_local(2)
    PetscInt, pointer :: datatocell_ids(:)
    PetscInt, pointer :: cell_ids_local(:)
    PetscBool         :: first_time
!  contains
!    procedure, public :: Init => DatasetMapInit
!    procedure, public :: Load => DatasetMapLoad
  end type dataset_map_type
  
  public :: DatasetMapCreate, &
            DatasetMapInit, &
            DatasetMapCast, &
            DatasetMapRead, &
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
  nullify(this%mapping)
  this%map_dims_global = 0
  this%map_dims_local = 0
  nullify(this%datatocell_ids)
  nullify(this%cell_ids_local)
  this%first_time = PETSC_TRUE
    
end subroutine DatasetMapInit


! ************************************************************************** !
!
! DatasetMapCast: Initializes members of common hdf5 dataset class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
function DatasetMapCast(this)

  use Dataset_Base_class
  
  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_map_type), pointer :: DatasetMapCast
  
  nullify(DatasetMapCast)
  select type (this)
    class is (dataset_map_type)
      DatasetMapCast => this
  end select
    
end function DatasetMapCast

! ************************************************************************** !
!
! DatasetMapRead: Reads in contents of a dataset card
! author: Glenn Hammond
! date: 01/12/11, 06/04/13
! 
! ************************************************************************** !
subroutine DatasetMapRead(this,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  class(dataset_map_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DATASET')
    call StringToUpper(keyword)   
      
    call DatasetCommonHDF5ReadSelectCase(this,input,keyword,found,option)
    if (.not.found) then
      select case(trim(keyword))
        case('MAP_HDF5_DATASET_NAME') 
          call InputReadWord(input,option,this%h5_dataset_map_name,PETSC_TRUE)
          call InputErrorMsg(input,option,'map name','DATASET')
        case default
          option%io_buffer = 'Keyword: ' // trim(keyword) // &
                             ' not recognized in dataset'    
          call printErrMsg(option)
      end select
    endif
  
  enddo
  
  if (len_trim(this%hdf5_dataset_name) < 1) then
    this%hdf5_dataset_name = this%name
  endif
  
end subroutine DatasetMapRead

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
  integer(HID_T) :: ndims_hdf5
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(2), length(2)
  PetscInt :: i
  PetscMPIInt :: array_rank_mpi
  PetscMPIInt :: hdf5_err
  PetscInt :: nids_local, remainder, istart, iend
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: dataset_name

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
  
  ! Open the "data" dataset
  dataset_name = 'Data'
  call h5dopen_f(grp_id,dataset_name,dataset_id,hdf5_err)
  call h5dget_space_f(dataset_id,file_space_id,hdf5_err)

  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(file_space_id,ndims_hdf5,hdf5_err)
  if (ndims_hdf5 /= 2) then
    option%io_buffer='Dimension of '// trim(this%h5_dataset_map_name) // &
      '/Data dataset in ' // trim(this%filename) // ' is not equal to 2.'
    call printErrMsg(option)
  endif

  ! Get dimensions of dataset
  allocate(dims_h5(ndims_hdf5))
  allocate(max_dims_h5(ndims_hdf5))
  call h5sget_simple_extent_dims_f(file_space_id,dims_h5,max_dims_h5,hdf5_err)
  
  nids_local=int(dims_h5(2)/option%mycommsize)
  remainder =int(dims_h5(2))-nids_local*option%mycommsize
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
!  option%io_buffer='Gautam: Modify code for reading HDF5 to generate mapping '//&
!   'of dataset'
!  call printMsg(option)
!  nids_local=dims_h5(2)
!  length(:) = dims_h5(:)
!  offset(:) = 0
  
  ! Save dimension size
  this%map_dims_global(:) = int(dims_h5(:))
  this%map_dims_local(:) = int(length(:))
  
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
  allocate(this%mapping(length(1), length(2)))

  ! Read the dataset collectively
  call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, this%mapping, &
                 dims_h5, hdf5_err, memory_space_id, file_space_id,prop_id)

  call h5pclose_f(prop_id,hdf5_err)

  if (memory_space_id > -1) call h5sclose_f(memory_space_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)
  call h5dclose_f(dataset_id,hdf5_err)  
  
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
  
  call DeallocateArray(this%mapping)
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
