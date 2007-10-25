module Condition_module
 
  implicit none

  private
  
#include "definitions.h"
 
  type, public :: condition_type
    integer :: id                                 ! id from which condition can be referenced
    integer :: itype                              ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: ctype         ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    integer :: num_values                         ! number of entries in the arrays of values
    real*8, pointer :: times(:)                   ! array of times between which linear interpolation of values occurs
    real*8, pointer :: values(:,:)                ! array of condition values, size(ndof,max_time_index)
    real*8, pointer :: cur_value(:)               ! current value of condition a time t, size(ndof)
    real*8 :: datum(3)                            ! location of reference value(s) in domain
    real*8 :: gradient(3)                         ! rate at which reference value(s) change(s) over 3D space
    integer :: cur_time_index, max_time_index     ! current and maximum time index in arrays
    type(condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type condition_type
  
  integer, save :: condition_count = 0
  
  public :: destroyCondition
  
contains

! ************************************************************************** !
!
! createCondition: Creates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createCondition(num_times,num_values)

  implicit none
  
  integer :: num_times
  integer :: num_values
  
  type(condition_type), pointer :: createCondition
  
  allocate(createCondition)
  nullify(createCondition%times)
  nullify(createCondition%values)
  nullify(createCondition%cur_value)
  nullify(createCondition%next)
  createCondition%id = 0
  createCondition%num_values = 0
  createCondition%itype = 0
  createCondition%ctype = ""
  createCondition%name = ""
  createCondition%datum = 0.d0
  createCondition%gradient = 0.d0
  createCondition%cur_time_index = 0
  createCondition%max_time_index = 0
  
  condition_count = condition_count + 1
  createCondition%id = condition_count
  createCondition%num_values = num_values
  createCondition%max_time_index = num_times
  allocate(createCondition%times(num_times))
  allocate(createCondition%values(num_values,num_times))

end function createCondition

! ************************************************************************** !
!
! destroyCondition: Deallocates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine destroyCondition(condition)

  implicit none
  
  type(condition_type), pointer :: condition
  
  if (associated(condition%times)) deallocate(condition%times)
  nullify(condition%times)
  if (associated(condition%values)) deallocate(condition%values)
  nullify(condition%values)
  nullify(condition%cur_value)
  nullify(condition%next)  
  deallocate(condition)
  nullify(condition)

end subroutine destroyCondition
  
end module Condition_module
