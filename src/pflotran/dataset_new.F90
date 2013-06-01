module Dataset_New_module
 
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_XYZ_class
  use Dataset_Map_class
  use Dataset_Global_class
  
  implicit none

  private
  
#include "definitions.h"

  public :: DatasetDestroy

contains

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
