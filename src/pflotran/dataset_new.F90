module Dataset_New_module
 
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_XYZ_class
  use Dataset_Map_class
  use Dataset_Global_class
  
  implicit none

  private
  
#include "definitions.h"

  public :: DatasetProcessDatasets, &
            DatasetDestroy

contains

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
  class(dataset_xyz_type), pointer :: dataset_xyz
  PetscBool :: swapped
  
  cur_dataset => datasets
  do
    if (.not.associated(cur_dataset)) exit
    nullify(dataset_xyz)
    swapped = PETSC_FALSE
    select type(cur_dataset)
      class is(dataset_common_hdf5_type)
        cur_dataset%is_cell_indexed = &
          DatasetCommonHDF5IsCellIndexed(cur_dataset,option)
        if (.not.cur_dataset%is_cell_indexed) then
          swapped = PETSC_TRUE
          dataset_xyz => DatasetXYZCreate()
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
! DatasetDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
recursive subroutine DatasetDestroy(dataset)

  implicit none
  
  class(dataset_base_type), pointer :: dataset
  
  class(dataset_global_type), pointer :: dataset_global
  class(dataset_xyz_type), pointer :: dataset_xyz
  class(dataset_map_type), pointer :: dataset_map
  class(dataset_common_hdf5_type), pointer :: dataset_common_hdf5
  class(dataset_base_type), pointer :: dataset_base

  if (.not.associated(dataset)) return
  
  if (associated(dataset%next)) then
    call DatasetDestroy(dataset%next)
  endif
  
  select type (selector => dataset)
    class is (dataset_global_type)
      dataset_global => selector
      call DatasetGlobalDestroy(dataset_global)
    class is (dataset_xyz_type)
      dataset_xyz => selector
      call DatasetXYZDestroy(dataset_xyz)
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

end module Dataset_New_module
