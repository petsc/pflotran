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
            DatasetProcessDatasets, &
            DatasetLoad, &
            DatasetFindInList, &
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
! DatasetProcessDatasets: Determine whether a dataset is indexed by cell 
!                             ids
! author: Glenn Hammond
! date: 03/26/12
!
! ************************************************************************** !
subroutine DatasetProcessDatasets(datasets,option)

  use Option_module
  
  implicit none
  
  class(dataset_base_type), pointer :: datasets
  type(option_type) :: option
  
  class(dataset_base_type), pointer :: cur_dataset
  class(dataset_base_type), pointer :: prev_dataset
  class(dataset_gridded_type), pointer :: dataset_xyz
  PetscBool :: swapped
  
  cur_dataset => datasets
  nullify(prev_dataset)
  do
    if (.not.associated(cur_dataset)) exit
    nullify(dataset_xyz)
    swapped = PETSC_FALSE
    select type(cur_dataset)
      class is(dataset_map_type)
        ! do nothing
      class is(dataset_global_type)
        ! do nothing
      class is(dataset_common_hdf5_type)
        cur_dataset%is_cell_indexed = &
          DatasetCommonHDF5IsCellIndexed(cur_dataset,option)
        if (.not.cur_dataset%is_cell_indexed) then
          swapped = PETSC_TRUE
          dataset_xyz => DatasetGriddedCreate()
          call DatasetCommonHDF5Copy(cur_dataset,dataset_xyz)
        endif
    end select
    if (swapped) then
      ! dataset_xyz%next is already set to cur_dataset%next
      if (associated(prev_dataset)) then
        prev_dataset%next => dataset_xyz
      else
        datasets => dataset_xyz
      endif
      ! just to be sure
      nullify(cur_dataset%next)
      call DatasetDestroy(cur_dataset)
      cur_dataset => dataset_xyz
    endif
    prev_dataset => cur_dataset
    cur_dataset => cur_dataset%next
  enddo
  
end subroutine DatasetProcessDatasets

! ************************************************************************** !
!
! DatasetVerify: Verifies that a dataset is intact and useable.
! author: Glenn Hammond
! date: 10/08/13
!
! ************************************************************************** !
subroutine DatasetVerify(dataset,option)

  implicit none

  class(dataset_base_type), pointer :: dataset
  type(option_type) :: option

  if (.not.associated(dataset)) return

  if (associated(dataset%time_storage)) then
  endif 

end subroutine DatasetVerify

! ************************************************************************** !
!
! DatasetLoad: Loads a dataset based on type
! author: Glenn Hammond
! date: 06/03/13
!
! ************************************************************************** !
recursive subroutine DatasetLoad(dataset, dm_wrapper, option)

  use DM_Kludge_module
  use Option_module
  
  implicit none
  
  class(dataset_base_type), pointer :: dataset
  type(dm_ptr_type), pointer :: dm_wrapper
  type(option_type) :: option

  class(dataset_global_type), pointer :: dataset_global
  class(dataset_gridded_type), pointer :: dataset_xyz
  class(dataset_map_type), pointer :: dataset_map
  class(dataset_base_type), pointer :: dataset_base

  select type (selector => dataset)
    class is (dataset_global_type)
      dataset_global => selector
      call DatasetGlobalLoad(dataset_global,dm_wrapper,option)
    class is (dataset_gridded_type)
      dataset_xyz => selector
      call DatasetGriddedLoad(dataset_xyz,option)
    class is (dataset_map_type)
      dataset_map => selector
      call DatasetMapLoad(dataset_map,option)
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
subroutine DatasetFindInList(list,dataset_base,error_string,option)

  use Dataset_Base_class
  use Dataset_Ascii_class
  use Option_module

  implicit none

  class(dataset_base_type), pointer :: list
  class(dataset_base_type), pointer :: dataset_base
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: dataset_name

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
    end select
  endif

end subroutine DatasetFindInList

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
  class(dataset_gridded_type), pointer :: dataset_xyz
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
      dataset_xyz => selector
      call DatasetGriddedDestroy(dataset_xyz)
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
