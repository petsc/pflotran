module Dataset_Base_class
 
  use Time_Storage_module
  
  implicit none

  private

#include "definitions.h"

  type, public :: dataset_base_type
    type(time_storage_type), pointer :: time_storage ! stores transient times
    PetscInt :: rank  ! size of dims(:)
    PetscInt, pointer :: dims(:)    ! dimensions of arrays (excludes time)
    PetscInt :: data_type
    PetscInt, pointer :: iarray(:)
    PetscInt, pointer :: ibuffer(:)
    PetscReal, pointer :: rarray(:)
    PetscReal, pointer :: rbuffer(:)
    PetscInt :: buffer_slice_offset ! index of the first time slice in the buffer
    PetscInt :: buffer_nslice ! # of time slices stored in buffer
!  contains
!    procedure, public :: Init => DatasetBaseInit
!    procedure, public :: InterpolateTime => DatasetBaseInterpolateTime
!    procedure, public :: Reorder => DatasetBaseReorder
!    procedure, public :: Strip => DatasetBaseStrip
  end type dataset_base_type

  ! dataset types
  PetscInt, parameter :: DATASET_INTEGER = 1
  PetscInt, parameter :: DATASET_REAL = 2
  
  public :: DatasetBaseInit, &
            DatasetBaseInterpolateTime, &
            DatasetBaseReorder, &
            DatasetBaseStrip, &
            DatasetBaseDestroy
contains

! ************************************************************************** !
!
! DatasetBaseInit: Initializes members of base database class
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetBaseInit(this)
  
  implicit none
  
  class(dataset_base_type) :: this
  
  this%rank = 0
  this%data_type = 0
  nullify(this%time_storage)
  nullify(this%dims)
  nullify(this%iarray)
  nullify(this%ibuffer)
  nullify(this%rarray)
  nullify(this%rbuffer)
  this%buffer_slice_offset = 0
  this%buffer_nslice = 0
    
end subroutine DatasetBaseInit

! ************************************************************************** !
!
! DatasetBaseInterpolateTime: Interpolates dataset between two buffer times
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetBaseInterpolateTime(this)

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  
  PetscInt :: array_size,i
  PetscInt :: time_interpolation_method
  PetscReal :: weight2
  PetscInt :: time1_start, time1_end, time2_start, time2_end
  
  if (.not.associated(this%rbuffer)) return
  
  time_interpolation_method = INTERPOLATION_STEP
  
  if (associated(this%time_storage)) then
    ! sets correct cur_time_index
    call TimeStorageUpdate(this%time_storage)
  endif
  
  if (this%time_storage%cur_time_index >= &
      this%time_storage%max_time_index) then
    ! dataset has reached the end of the time array and is not cyclic.
    ! no interpolation needed
    if (this%time_storage%cur_time_index_changed) then
      array_size = size(this%rarray)
      time1_end = array_size*(this%time_storage%max_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      this%rarray = this%rbuffer(time1_start:time1_end)
    endif
    return
  endif
  
  array_size = size(this%rarray)
  select case(time_interpolation_method)
    case(INTERPOLATION_STEP)
      ! if time index has not changed skip
      if (.not.this%time_storage%cur_time_index_changed) return
      time1_end = array_size*(this%time_storage%cur_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      this%rarray = this%rbuffer(time1_start:time1_end)
    case(INTERPOLATION_LINEAR)
      ! if fraction has not changed skip
      if (.not.this%time_storage%cur_time_fraction_changed) return
      time1_end = array_size*(this%time_storage%cur_time_index - &
                              this%buffer_slice_offset)
      time1_start = time1_end - array_size + 1
      time2_end = time1_end + array_size
      time2_start = time1_start + array_size
      this%rarray = (1.d0-this%time_storage%cur_time_fraction) * &
                        this%rbuffer(time1_start:time1_end) + &
                        this%time_storage%cur_time_fraction * &
                        this%rbuffer(time2_start:time2_end)
  end select

end subroutine DatasetBaseInterpolateTime

! ************************************************************************** !
!
! DatasetBaseInterpolateSpace: Interpolates data from the dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetBaseInterpolateSpace(this,xx,yy,zz,time,real_value,option)

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  class(dataset_base_type) :: this
  PetscReal, intent(in) :: xx, yy, zz
  PetscReal :: time
  PetscReal :: real_value
  type(option_type) :: option
  
end subroutine DatasetBaseInterpolateSpace

! ************************************************************************** !
!
! DatasetBaseReorder: If a dataset is loaded from an HDF5 file, and it was
!              multidimensional in the HDF5 file, the array needs to be
!              reordered fro Fortran -> C indexing.  This subroutine
!              takes care of the reordering.
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetBaseReorder(this,option)

  use Option_module
  
  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  
#if 0
  PetscReal, allocatable :: temp_real(:)
  PetscInt :: i, j, k, l
  PetscInt :: dims(4), n1, n1Xn2, n1Xn2Xn3
  PetscInt :: count, index
  PetscReal, pointer :: rarray(:)
  
  if (this%data_type == DATASET_INTEGER) then
    option%io_buffer = 'Reordering of integer data sets not yet supported.'
    call printErrMsg(option)
  endif
  
  ! set each dim to 1 by default for loop below
  dims(:) = 1
  dims(1:this%rank) = this%dims(1:this%rank)
  if (associated(this%rbuffer)) then
    rarray => this%rbuffer
    dims(this%rank+1) = this%time_storage%max_time_index
  else
    rarray => this%rarray
  endif
  
  ! Not necessary for 1D arrays
  if (maxval(dims(2:)) == 1) return
  
  l = 1
  do i = 1, size(dims)
    l = l*dims(i)
  enddo
  allocate(temp_real(l))
  
  n1 = dims(1)
  n1Xn2 = n1*dims(2)
  n1Xn2Xn3 = n1Xn2*dims(3)
  do i = 1, dims(1)
    do j = 0, dims(2)-1
      do k = 0, dims(3)-1
        do l = 0, dims(4)-1
          index = l*n1Xn2Xn3+k*n1Xn2+j*n1+i
          count = count+1
          temp_real(index) = rarray(count)
        enddo
      enddo
    enddo
  enddo  

  rarray = temp_real
  deallocate(temp_real)
#endif
  
end subroutine DatasetBaseReorder

! ************************************************************************** !
!
! DatasetBasePrintMe: Prints dataset info
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetBasePrintMe(this,option)

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  option%io_buffer = 'TODO(geh): add DatasetPrint()'
  call printMsg(option)
            
end subroutine DatasetBasePrintMe

! ************************************************************************** !
!
! DatasetBaseGetTimes: Fills an array of times based on a dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetBaseGetTimes(this, option, max_sim_time, time_array)

  use Option_module

  implicit none
  
  class(dataset_base_type) :: this
  type(option_type) :: option
  PetscReal :: max_sim_time
  PetscReal, pointer :: time_array(:)
  
  
  if (associated(this%time_storage)) then
    call TimeStorageGetTimes(this%time_storage, option, max_sim_time, &
                             time_array)
  endif
 
end subroutine DatasetBaseGetTimes

! ************************************************************************** !
!
! DatasetBaseStrip: Strips allocated objects within base dataset object
! author: Glenn Hammond
! date: 05/03/13
!
! ************************************************************************** !
subroutine DatasetBaseStrip(this)

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(dataset_base_type)  :: this
  
  call DeallocateArray(this%iarray)
  call DeallocateArray(this%rarray)
  call DeallocateArray(this%ibuffer)
  call DeallocateArray(this%rbuffer)
  call DeallocateArray(this%dims)
  
end subroutine DatasetBaseStrip

! ************************************************************************** !
!
! DatasetBaseDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11, 05/03/13
!
! ************************************************************************** !
subroutine DatasetBaseDestroy(dataset)

  implicit none
  
  class(dataset_base_type), pointer :: dataset
  
  if (.not.associated(dataset)) return
  
  call DatasetBaseStrip(dataset)
  
  deallocate(dataset)
  nullify(dataset)
  
end subroutine DatasetBaseDestroy

end module Dataset_Base_class
