module Dataset_module
 
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class
  use Dataset_Map_class
  use Dataset_Global_class
  
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

  use Input_module
  use Option_module
  use String_module
  
  implicit none

  type(input_type) :: input
  class(dataset_base_type), pointer :: dataset
  class(dataset_map_type), pointer :: dataset_map
  class(dataset_global_type), pointer :: dataset_global
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
      dataset_map => DatasetMapCreate()
      call InputReadWord(input,option,dataset_map%name,PETSC_TRUE)
      call InputDefaultMsg(input,option,'DATASET name') 
      call DatasetMapRead(dataset_map,input,option)
      dataset => dataset_map
    case('GLOBAL')
      dataset_global => DatasetGlobalCreate()
      call InputReadWord(input,option,dataset_global%name,PETSC_TRUE)
      call InputDefaultMsg(input,option,'DATASET name') 
      call DatasetCommonHDF5Read(dataset_global,input,option)
      dataset => dataset_global
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
!                                 to dataset_gridded_type
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
  class(dataset_gridded_type), pointer :: dataset_gridded
  PetscBool :: swapped
  
  cur_dataset => datasets
  nullify(prev_dataset)
  do
    if (.not.associated(cur_dataset)) exit
    nullify(dataset_gridded)
    swapped = PETSC_FALSE
    select type(cur_dataset)
      class is(dataset_map_type)
        ! do nothing
      class is(dataset_global_type)
        ! do nothing
      class is(dataset_common_hdf5_type)
        ! if class is dataset_common_hdf5_type, the dataset should really 
        ! be dataset_gridded unless cell indexed
        cur_dataset%is_cell_indexed = &
          DatasetCommonHDF5IsCellIndexed(cur_dataset,option)
        ! if not cell indexed, we need to change the type to the extended
        ! dataset_gridded type
        if (.not.cur_dataset%is_cell_indexed) then
          swapped = PETSC_TRUE
          dataset_gridded => DatasetGriddedCreate()
          call DatasetCommonHDF5Copy(cur_dataset,dataset_gridded)
        endif
      class default
        option%io_buffer = &
          'Unknown dataset type in DatasetScreenForNonCellIndexed.'
        call printErrMsg(option)
    end select
    ! if we changed the derived type, we need to replace the old dataset
    ! in the linked list of datasets and destroy the old dataset.
    if (swapped) then
      ! dataset_gridded%next is already set to cur_dataset%next
      if (associated(prev_dataset)) then
        prev_dataset%next => dataset_gridded
      else
        datasets => dataset_gridded
      endif
      ! just to be sure
      nullify(cur_dataset%next)
      call DatasetDestroy(cur_dataset)
      cur_dataset => dataset_gridded
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

  class(dataset_global_type), pointer :: dataset_global
  class(dataset_gridded_type), pointer :: dataset_gridded
  class(dataset_map_type), pointer :: dataset_map
  class(dataset_base_type), pointer :: dataset_base

  select type (selector => dataset)
    class is (dataset_global_type)
      dataset_global => selector
      call DatasetGlobalLoad(dataset_global,option)
    class is (dataset_gridded_type)
      dataset_gridded => selector
      call DatasetGriddedLoad(dataset_gridded,option)
    class is (dataset_map_type)
      dataset_map => selector
      call DatasetMapLoad(dataset_map,option)
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
    class is (dataset_global_type)
      DatasetGetClass = 'dataset_global_type'
    class is (dataset_gridded_type)
      DatasetGetClass = 'dataset_gridded_type'
    class is (dataset_map_type)
      DatasetGetClass = 'dataset_map_type'
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
  
  class(dataset_global_type), pointer :: dataset_global
  class(dataset_gridded_type), pointer :: dataset_gridded
  class(dataset_map_type), pointer :: dataset_map
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
    class is (dataset_global_type)
      dataset_global => selector
      call DatasetGlobalDestroy(dataset_global)
    class is (dataset_gridded_type)
      dataset_gridded => selector
      call DatasetGriddedDestroy(dataset_gridded)
    class is (dataset_map_type)
      dataset_map => selector
      call DatasetMapDestroy(dataset_map)
    class is (dataset_common_hdf5_type)
      dataset_common_hdf5 => selector
      call DatasetCommonHDF5Destroy(dataset_common_hdf5)
    class is (dataset_base_type)
      dataset_base => selector
      call DatasetBaseDestroy(dataset_base)
  end select
  
end subroutine DatasetDestroy

end module Dataset_module
