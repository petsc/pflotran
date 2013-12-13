module Dataset_module
 
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Map_HDF5_class
  use Dataset_Global_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"

  public :: DatasetRead, &
            DatasetScreenForNonCellIndexed, &
            DatasetLoad, &
            DatasetVerify, &
            DatasetUpdate, &
            DatasetFindInList, &
            DatasetIsTransient, &
            DatasetGetClass, &
            DatasetPrint, &
            DatasetDestroy

contains

! *************************************************************************** !
!
! DatasetRead: Reads a dataset from the input file
! author: Glenn Hammond
! date: 03/26/12
!
! ************************************************************************** !
subroutine DatasetRead(input,dataset,option)

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none

  type(input_type) :: input
  class(dataset_base_type), pointer :: dataset
  class(dataset_map_hdf5_type), pointer :: dataset_map_hdf5
  class(dataset_global_hdf5_type), pointer :: dataset_global_hdf5
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, word2

  ! read from the buffer on the first line the type of dataset.  If it does
  ! not match any of type in the select case, assume default and use the 
  ! word as the name.
  
  call InputReadWord(input,option,word,PETSC_TRUE)
  word2 = word
  call StringToUpper(word2)
  
  select case(word2)
    case('MAPPED')
      dataset_map_hdf5 => DatasetMapHDF5Create()
      call InputReadWord(input,option,dataset_map_hdf5%name,PETSC_TRUE)
      call InputDefaultMsg(input,option,'DATASET name') 
      call DatasetMapHDF5Read(dataset_map_hdf5,input,option)
      dataset => dataset_map_hdf5
    case('GLOBAL')
      dataset_global_hdf5 => DatasetGlobalHDF5Create()
      call InputReadWord(input,option,dataset_global_hdf5%name,PETSC_TRUE)
      call InputDefaultMsg(input,option,'DATASET name') 
      call DatasetCommonHDF5Read(dataset_global_hdf5,input,option)
      dataset => dataset_global_hdf5
    case default ! CELL_INDEXED, GLOBAL, GRIDDED
      dataset => DatasetCommonHDF5Create()
      select case(word2)
        case('CELL_INDEXED', 'GRIDDED')
          call InputReadWord(input,option,dataset%name,PETSC_TRUE)
        case default
          if (.not.InputError(input)) then
            dataset%name = trim(word)
          endif
      end select
      call InputDefaultMsg(input,option,'DATASET name') 
      call DatasetCommonHDF5Read(DatasetCommonHDF5Cast(dataset),input, &
                                 option)
  end select

end subroutine DatasetRead

! *************************************************************************** !
!
! DatasetScreenForNonCellIndexed: Recasts datasets of dataset_common_hdf5_type
!                                 to dataset_gridded_hdf5_type
! author: Glenn Hammond
! date: 03/26/12
!
! ************************************************************************** !
subroutine DatasetScreenForNonCellIndexed(datasets,option)

  use Option_module
  
  implicit none
  
  class(dataset_base_type), pointer :: datasets
  type(option_type) :: option
  
  class(dataset_base_type), pointer :: cur_dataset
  class(dataset_base_type), pointer :: prev_dataset
  class(dataset_gridded_hdf5_type), pointer :: dataset_gridded_hdf5
  PetscBool :: swapped
  
  cur_dataset => datasets
  nullify(prev_dataset)
  do
    if (.not.associated(cur_dataset)) exit
    nullify(dataset_gridded_hdf5)
    swapped = PETSC_FALSE
    select type(cur_dataset)
      class is(dataset_map_hdf5_type)
        ! do nothing
      class is(dataset_global_hdf5_type)
        ! do nothing
      class is(dataset_common_hdf5_type)
        ! if class is dataset_common_hdf5_type, the dataset should really 
        ! be dataset_gridded_hdf5 unless cell indexed
        cur_dataset%is_cell_indexed = &
          DatasetCommonHDF5IsCellIndexed(cur_dataset,option)
        ! if not cell indexed, we need to change the type to the extended
        ! dataset_gridded_hdf5 type
        if (.not.cur_dataset%is_cell_indexed) then
          swapped = PETSC_TRUE
          dataset_gridded_hdf5 => DatasetGriddedHDF5Create()
          call DatasetCommonHDF5Copy(cur_dataset,dataset_gridded_hdf5)
        endif
      class default
        option%io_buffer = &
          'Unknown dataset type in DatasetScreenForNonCellIndexed.'
        call printErrMsg(option)
    end select
    ! if we changed the derived type, we need to replace the old dataset
    ! in the linked list of datasets and destroy the old dataset.
    if (swapped) then
      ! dataset_gridded_hdf5%next is already set to cur_dataset%next
      if (associated(prev_dataset)) then
        prev_dataset%next => dataset_gridded_hdf5
      else
        datasets => dataset_gridded_hdf5
      endif
      ! just to be sure
      nullify(cur_dataset%next)
      call DatasetDestroy(cur_dataset)
      cur_dataset => dataset_gridded_hdf5
    endif
    prev_dataset => cur_dataset
    cur_dataset => cur_dataset%next
  enddo
  
end subroutine DatasetScreenForNonCellIndexed

! ************************************************************************** !
!
! DatasetVerify: Verifies that a dataset is intact and useable.
! author: Glenn Hammond
! date: 10/08/13
!
! ************************************************************************** !
subroutine DatasetVerify(dataset,default_time_storage,option)

  use Time_Storage_module
  use Option_module

  implicit none

  class(dataset_base_type), pointer :: dataset
  type(time_storage_type), pointer :: default_time_storage
  type(option_type) :: option
  
  type(time_storage_type), pointer :: default

  if (.not.associated(dataset)) return

  call TimeStorageVerify(0.d0,dataset%time_storage,default_time_storage,option)
  select type(dataset_ptr => dataset)
    class is(dataset_ascii_type)
      call DatasetAsciiVerify(dataset_ptr,option)
    class is(dataset_base_type)
      call DatasetBaseVerify(dataset_ptr,option)
    class default
      option%io_buffer = 'DatasetXXXVerify needed for unknown dataset type'
      call printErrMsg(option)
  end select
  
end subroutine DatasetVerify

! ************************************************************************** !
!
! DatasetUpdate: Updates a dataset based on type
! author: Glenn Hammond
! date: 10/09/13
!
! ************************************************************************** !
recursive subroutine DatasetUpdate(dataset,time,option)

  use Option_module
  
  implicit none
  
  PetscReal :: time
  class(dataset_base_type), pointer :: dataset
  type(option_type) :: option

  class(dataset_ascii_type), pointer :: dataset_ascii

  if (.not.associated(dataset)) return

  if (associated(dataset%time_storage)) then
    dataset%time_storage%cur_time = time
  endif
  select type (selector => dataset)
    class is (dataset_ascii_type)
      dataset_ascii => selector
      call DatasetAsciiUpdate(dataset_ascii,option)
    class is (dataset_common_hdf5_type)
      if (associated(selector%time_storage) .and. &
          .not.selector%is_cell_indexed) then
        call DatasetLoad(dataset,option)
      endif
  end select
  
end subroutine DatasetUpdate

! ************************************************************************** !
!
! DatasetLoad: Loads a dataset based on type
! author: Glenn Hammond
! date: 06/03/13
!
! ************************************************************************** !
recursive subroutine DatasetLoad(dataset,option)

  use DM_Kludge_module
  use Option_module
  
  implicit none
  
  class(dataset_base_type), pointer :: dataset
  type(option_type) :: option

  class(dataset_global_hdf5_type), pointer :: dataset_global_hdf5
  class(dataset_gridded_hdf5_type), pointer :: dataset_gridded_hdf5
  class(dataset_map_hdf5_type), pointer :: dataset_map_hdf5
  class(dataset_common_hdf5_type), pointer :: dataset_common_hdf5
  class(dataset_base_type), pointer :: dataset_base

  select type (selector => dataset)
    class is (dataset_global_hdf5_type)
      dataset_global_hdf5 => selector
      call DatasetGlobalHDF5Load(dataset_global_hdf5,option)
    class is (dataset_gridded_hdf5_type)
      dataset_gridded_hdf5 => selector
      call DatasetGriddedHDF5Load(dataset_gridded_hdf5,option)
    class is (dataset_map_hdf5_type)
      dataset_map_hdf5 => selector
      call DatasetMapHDF5Load(dataset_map_hdf5,option)
    class is (dataset_common_hdf5_type)
      dataset_common_hdf5 => selector
      if (dataset_common_hdf5%is_cell_indexed) then
        option%io_buffer = 'Cell Indexed datasets opened later.'
        call printMsg(option)
      else
        option%io_buffer = 'Unrecognized dataset that extends ' // &
          'dataset_common_hdf5 in DatasetLoad.'
        call printErrMsg(option)
      endif
    class is (dataset_ascii_type)
      ! do nothing.  dataset should already be in memory at this point
    class is (dataset_base_type)
      dataset_base => selector
      option%io_buffer = 'DatasetLoad not yet supported for base dataset class.'
      call printErrMsg(option)
  end select
  
end subroutine DatasetLoad

! ************************************************************************** !
!
! DatasetFindInList: Uses a dummy dataset with name to find the actual
!                    dataset in a list of datasets
! author: Glenn Hammond
! date: 10/07/13
!
! ************************************************************************** !
subroutine DatasetFindInList(list,dataset_base,default_time_storage, &
                             error_string,option)

  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module
  use Option_module

  implicit none

  class(dataset_base_type), pointer :: list
  class(dataset_base_type), pointer :: dataset_base
  type(time_storage_type), pointer :: default_time_storage
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: dataset_name
  PetscReal, parameter :: time = 0.d0

  ! check for dataset in flow_dataset
  if (associated(dataset_base)) then
    select type(dataset => dataset_base)
      class is(dataset_ascii_type)
        ! do nothing as the correct dataset already exists
      class default
        dataset_name = dataset%name
        ! delete the dataset since it is solely a placeholder
        call DatasetDestroy(dataset_base)
        ! get dataset from list
        dataset_base => &
          DatasetBaseGetPointer(list,dataset_name,error_string,option)
        ! once a dataset is linked to the dataset list, it needs to be loaded
        ! immediately
        call DatasetLoad(dataset_base,option)
        call DatasetVerify(dataset_base,default_time_storage,option)
        ! must update after DatasetVerify since the time interpolation method
        ! may not have been properly set during the load! Force the update.
        if (associated(dataset_base%time_storage)) then
          dataset_base%time_storage%force_update = PETSC_TRUE  
        endif
        call DatasetUpdate(dataset_base,time,option)
    end select
  endif

end subroutine DatasetFindInList

! ************************************************************************** !
!
! DatasetIsTransient: Determines whether a dataset is transient
! author: Glenn Hammond
! date: 10/10/13
!
! ************************************************************************** !
function DatasetIsTransient(dataset)

  use Time_Storage_module

  implicit none

  class(dataset_base_type), pointer :: dataset
  
  PetscBool :: DatasetIsTransient

  DatasetIsTransient = PETSC_FALSE
  if (associated(dataset)) then
    if (associated(dataset%time_storage)) then
      DatasetIsTransient = PETSC_TRUE
    endif
  endif 
  
end function DatasetIsTransient


! ************************************************************************** !
!
! DatasetPrint: Prints dataset info
! author: Glenn Hammond
! date: 10/22/13
!
! ************************************************************************** !
subroutine DatasetPrint(this,option)

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option

  write(option%fid_out,'(8x,''Dataset: '',a)') trim(this%name)
  write(option%fid_out,'(10x,''Type: '',a)') trim(DatasetGetClass(this))

  call DatasetBasePrint(this,option)
  
  select type (d=>this)
    class is (dataset_ascii_type)
      call DatasetAsciiPrint(d,option)
    class is (dataset_global_hdf5_type)
      call DatasetGlobalHDF5Print(d,option)
    class is (dataset_gridded_hdf5_type)
      call DatasetGriddedHDF5Print(d,option)
    class is (dataset_map_hdf5_type)
      call DatasetMapHDF5Print(d,option)
    class is (dataset_common_hdf5_type)
      call DatasetCommonHDF5Print(d,option)
    class default
      option%io_buffer = 'Unknown dataset type for dataset "' // &
        trim(this%name) // '" in DatasetPrint()'
  end select
            
end subroutine DatasetPrint

! ************************************************************************** !
!
! DatasetGetClass: Returns a string defining the class of dataset
! author: Glenn Hammond
! date: 10/10/13
!
! ************************************************************************** !
function DatasetGetClass(dataset)

  use Option_module

  implicit none
  
  class(dataset_base_type) :: dataset

  character(len=MAXSTRINGLENGTH) :: DatasetGetClass
  
  select type (dataset)
    class is (dataset_ascii_type)
      DatasetGetClass = 'dataset_ascii_type'
    class is (dataset_global_hdf5_type)
      DatasetGetClass = 'dataset_global_hdf5_type'
    class is (dataset_gridded_hdf5_type)
      DatasetGetClass = 'dataset_gridded_hdf5_type'
    class is (dataset_map_hdf5_type)
      DatasetGetClass = 'dataset_map_hdf5_type'
    class is (dataset_common_hdf5_type)
      DatasetGetClass = 'dataset_hdf5_type'
    class is (dataset_base_type)
      DatasetGetClass = 'dataset_base_type'
    class default
      DatasetGetClass = 'dataset class unknown'
  end select
  
end function DatasetGetClass

! ************************************************************************** !
!
! DatasetDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
recursive subroutine DatasetDestroy(dataset)

  implicit none
  
  class(dataset_base_type), pointer :: dataset
  
  class(dataset_global_hdf5_type), pointer :: dataset_global_hdf5
  class(dataset_gridded_hdf5_type), pointer :: dataset_gridded_hdf5
  class(dataset_map_hdf5_type), pointer :: dataset_map_hdf5
  class(dataset_common_hdf5_type), pointer :: dataset_common_hdf5
  class(dataset_ascii_type), pointer :: dataset_ascii
  class(dataset_base_type), pointer :: dataset_base

  if (.not.associated(dataset)) return
  
  if (associated(dataset%next)) then
    call DatasetDestroy(dataset%next)
  endif
  
  select type (selector => dataset)
    class is (dataset_ascii_type)
      dataset_ascii => selector
      call DatasetAsciiDestroy(dataset_ascii)
    class is (dataset_global_hdf5_type)
      dataset_global_hdf5 => selector
      call DatasetGlobalHDF5Destroy(dataset_global_hdf5)
    class is (dataset_gridded_hdf5_type)
      dataset_gridded_hdf5 => selector
      call DatasetGriddedHDF5Destroy(dataset_gridded_hdf5)
    class is (dataset_map_hdf5_type)
      dataset_map_hdf5 => selector
      call DatasetMapHDF5Destroy(dataset_map_hdf5)
    class is (dataset_common_hdf5_type)
      dataset_common_hdf5 => selector
      call DatasetCommonHDF5Destroy(dataset_common_hdf5)
    class is (dataset_base_type)
      dataset_base => selector
      call DatasetBaseDestroy(dataset_base)
  end select
  
end subroutine DatasetDestroy

end module Dataset_module
