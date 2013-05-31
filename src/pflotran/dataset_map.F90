module Dataset_CLM_class
 
  use Dataset_XYZ_class
  
  implicit none

  private

#include "definitions.h"

  type, public, extends(dataset_xyz_type) :: dataset_clm_type
    character(len=MAXSTRINGLENGTH) :: h5_dataset_map_name
    character(len=MAXSTRINGLENGTH) :: map_filename
    PetscInt,pointer :: map(:,:)
    PetscInt         :: map_dims_global(2)
    PetscInt         :: map_dims_local(2)
    PetscInt,pointer :: datatocell_ids(:)
    PetscInt,pointer :: cell_ids_local(:)
    PetscBool        :: first_time
!  contains
!    procedure, public :: Init => DatasetCLMInit
!    procedure, public :: Load => DatasetCLMLoad
  end type dataset_clm_type
  
  public :: DatasetCLMCreate, &
            DatasetCLMInit, &
            DatasetCLMLoad, &
            DatasetCLMDestroy
  
contains

! ************************************************************************** !
!
! DatasetCLMCreate: Creates global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
function DatasetCLMCreate()
  
  implicit none
  
  class(dataset_clm_type), pointer :: dataset

  class(dataset_clm_type), pointer :: DatasetCLMCreate
  
  allocate(dataset)
  call DatasetCLMInit(dataset)

  DatasetCLMCreate => dataset
    
end function DatasetCLMCreate

! ************************************************************************** !
!
! DatasetCLMInit: Initializes members of global dataset class
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetCLMInit(this)
  
  implicit none
  
  class(dataset_clm_type) :: this
  
  call DatasetXYZInit(this)
  this%h5_dataset_map_name = ''
  this%map_filename = ''
  nullify(this%map)
  this%map_dims_global = 0
  this%map_dims_local = 0
  nullify(this%datatocell_ids)
  nullify(this%cell_ids_local)
  this%first_time = PETSC_TRUE
    
end subroutine DatasetCLMInit

! ************************************************************************** !
!
! DatasetCLMLoad: Load new data into dataset buffer
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetCLMLoad(this,option)
  
  use Option_module
  use Time_Storage_module

  implicit none
  
  class(dataset_clm_type) :: this
  type(option_type) :: option
  
  if (DatasetCommonHDF5Load(this,option)) then
#if defined(PETSC_HAVE_HDF5)    
    call DatasetCLMReadData(this,option)
#endif    
!    call this%Reorder(option)
    call DatasetBaseReorder(this,option)
  endif
  call DatasetBaseInterpolateTime(this)
    
end subroutine DatasetCLMLoad

! ************************************************************************** !
!
! DatasetCLMStrip: Strips allocated objects within CLM dataset object
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetCLMStrip(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_clm_type)  :: this
  
  call DatasetXYZStrip(this)
  
  call DeallocateArray(this%map)
  call DeallocateArray(this%datatocell_ids)
  call DeallocateArray(this%cell_ids_local)
  
end subroutine DatasetCLMStrip

! ************************************************************************** !
!
! DatasetCLMDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 05/29/13
!
! ************************************************************************** !
subroutine DatasetCLMDestroy(this)

  implicit none
  
  class(dataset_clm_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetCLMStrip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetCLMDestroy

end module Dataset_CLM_class
