module Dataset_Global_class
 
  use Dataset_Base_class
  
  implicit none

  private

#include "definitions.h"

  type, public, extends(dataset_base_type) :: dataset_global_type
    character(len=MAXWORDLENGTH) :: dataset_name
    character(len=MAXSTRINGLENGTH) :: filename
  contains
    procedure, public :: Init => DatabaseGlobalInit
    procedure, public :: Load => DatabaseGlobalLoad
  end type dataset_global_type
  
  public :: DatabaseGlobalCreate, &
            DatasetGlobalDestroy
  
contains

! ************************************************************************** !
!
! DatabaseGlobalCreate: Creates global database class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
function DatabaseGlobalCreate()
  
  implicit none
  
  class(dataset_global_type), pointer :: DatabaseGlobalCreate
  
  allocate(DatabaseGlobalCreate)
  call DatabaseGlobalCreate%Init()
    
end function DatabaseGlobalCreate

! ************************************************************************** !
!
! DatabaseGlobalInit: Initializes members of global database class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatabaseGlobalInit(this)
  
  implicit none
  
  class(dataset_global_type) :: this
  
  call DatabaseBaseInit(this)
  this%filename = ''
  this%dataset_name = ''
    
end subroutine DatabaseGlobalInit

! ************************************************************************** !
!
! DatabaseGlobalLoad: Load new data into database buffer
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatabaseGlobalLoad(this,discretization,grid,option)
  
#if defined(PETSC_HAVE_HDF5)    
  use hdf5
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
  
  this%time_storage%cur_time = option%time
  if (this%time_storage%cur_time_index >= &
        this%buffer_slice_offset + this%buffer_nslice &
      .and. &
      this%time_storage%cur_time_index < &
        this%time_storage%max_time_index) then
    if (.not.associated(this%rarray)) then
      allocate(this%rarray(grid%nlmax))
      this%buffer_nslice = 10
      allocate(this%rbuffer(grid%nlmax*this%buffer_nslice))
    endif
    this%buffer_slice_offset = this%time_storage%cur_time_index - 1
#if defined(PETSC_HAVE_HDF5)    
    call HDF5CReadGlobalArray(discretization,grid,option,this%filename, &
                              this%dataset_name, &
                              this%buffer_slice_offset, &
                              this%rbuffer, &
                              H5T_NATIVE_DOUBLE)
#endif    
    call this%Reorder(option)
    call this%InterpolateTime()
  endif
    
end subroutine DatabaseGlobalLoad

! ************************************************************************** !
!
! BasePrintMe: Prints dataset info
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine BasePrintMe(dataset,option)

  use Option_module

  implicit none
  
  class(dataset_global_type) :: dataset
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  option%io_buffer = 'TODO(geh): add DatasetPrint()'
  call printMsg(option)
            
end subroutine BasePrintMe

! ************************************************************************** !
!
! BaseGetTimes: Fills an array of times based on a dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine BaseGetTimes(dataset, option, max_sim_time, time_array)

  use Option_module

  implicit none
  
  class(dataset_global_type) :: dataset
  type(option_type) :: option
  PetscReal :: max_sim_time
  PetscReal, pointer :: time_array(:)
  
  
  if (associated(dataset%time_storage)) then
    call TimeStorageGetTimes(dataset%time_storage, option, max_sim_time, &
                             time_array)
  endif
 
end subroutine BaseGetTimes

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
  
  call this%Strip()
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetGlobalDestroy

end module Dataset_Global_class
