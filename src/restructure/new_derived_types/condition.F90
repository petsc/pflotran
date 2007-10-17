module Condition_module

 implicit none
 
! condition types
#define DIRICHLET 1
#define NEUMANN 2
#define HYDROSTATIC 3

#include "definitions.h"

  private

  type, public :: condition
    integer :: id
    integer :: itype  ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: ctype ! type name
    character(len=MAXWORDLENGTH) :: name
    integer :: num_values
    real*8, pointer :: times(:)
    real*8, pointer :: values(:,:) ! (ndof,max_time_index)
    real*8, pointer :: cur_value(:) ! (ndof)
    real*8 :: datum(3) ! location of reference values
    real*8 :: gradient(3) ! rate at which values changes over 3D space
    integer :: cur_time_index, max_time_index
    type(condition), pointer :: next
  end type condition

  ! pointer data structure required for making an array of condition pointers in F90
  type :: condition_ptr
    type(condition), pointer :: ptr
  end type condition_ptr

  ! linked list of conditions
  type(condition), pointer :: condition_list
  ! pointer to last condition; facilitates appendage of conditions
  type(condition), pointer, private :: last_condition
  ! array of pointers to conditions in list
  type(condition_ptr), pointer :: condition_array(:)
  ! number of conditions in list
  integer, private :: num_conditions


  public initConditionModule, readCondition, getConditionFromName
  
contains
  
! ************************************************************************** !
!
! initConditionModule: Initializes the condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine initConditionModule()

  implicit none

  nullify(condition_list)
  nullify(last_condition)
  num_conditions = 0

end subroutine initConditionModule

! ************************************************************************** !
!
! readCondition: add condition to condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
integer function readCondition(fid)

  implicit none
  
  integer :: fid
  
  type(condition), pointer :: new_condition
  
  new_condition => createCondition()
  print *, 'Code for condition.F90:readCondition still to be written'
  stop
  call addCondition(new_condition)
  
  readCondition = new_condition%id

end function readCondition

! ************************************************************************** !
!
! getConditionFromName: Returns the pointer of the first condition in the  
!                       condition list with a matching name
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
function getConditionFromName(name)

  use fileio_module
  
  implicit none
  
  type(condition), pointer :: getConditionFromName
  character(len=MAXWORDLENGTH) ::  name
  
  type(condition), pointer :: cur_condition

  nullify(getConditionFromName)
  
  cur_condition => condition_list
  do 
    if (.not.associated(cur_condition)) exit
    if (fiStringCompare(cur_condition%name,name,len_trim(name))) then
      getConditionFromName => cur_condition 
      return
    endif
    cur_condition => cur_condition%next
  enddo

end function getConditionFromName

! ************************************************************************** !
!
! addConditionToList: add condition to condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine addConditionToList(new_condition)

  implicit none

  type(condition), pointer :: new_condition

  new_condition%cur_time_index = 1

  num_conditions = num_conditions + 1
  new_condition%id =  num_conditions
  if (.not.associated(condition_list)) condition_list => new_condition
  if (associated(last_condition)) last_condition%next => new_condition
  last_condition => new_condition

end subroutine addConditionToList

! ************************************************************************** !
!
! createCondition: creates a new condition
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function createCondition()

  implicit none

  type(condition), pointer :: createCondition

  allocate(createCondition)
  createCondition%id = 0
  createCondition%itype = 0
  createCondition%ctype = ' '
  createCondition%num_values = 0
  nullify(createCondition%times)
  nullify(createCondition%values)
  nullify(createCondition%cur_value)
  createCondition%datum = 0
  createCondition%gradient = 0
  createCondition%cur_time_index = 0
  createCondition%max_time_index = 0

end function createCondition

! ************************************************************************** !
!
! convertConditionListToArray: Converts the linked list of conditions to an 
!                              array of conditions
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine convertConditionListToArray()

  implicit none

  integer :: count
  type(condition), pointer :: cur_condition

  allocate(condition_array(num_conditions))

  cur_condition => condition_list
  do
    if (.not.associated(cur_condition)) exit
    condition_array(cur_condition%id)%ptr => cur_condition
    cur_condition => cur_condition%next
  enddo

end subroutine convertConditionListToArray

! ************************************************************************** !
!
! UpdateConditions: Updates the value of the boundary/initial condition based
!                   on the current time
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine updateConditions(time)

  implicit none

  real*8 :: time

  type(condition), pointer :: cur_condition

  integer :: cur_time_index, next_time_index
  real*8 :: time_fraction

  cur_condition => condition_list
  do

    if (.not.associated(cur_condition)) exit

    ! if more than 10 times specified, cycle the transient condition
    ! note that 10 is just an arbitrary number
    if (cur_condition%cur_time_index == cur_condition%max_time_index .and. &
        cur_condition%max_time_index > 10) then
      do cur_time_index = 1, cur_condition%max_time_index
        cur_condition%times(cur_time_index) = &     
              cur_condition%times(cur_time_index) + &     
              cur_condition%times(cur_condition%max_time_index)
      enddo
      cur_condition%cur_time_index = 1
    endif

    cur_time_index = cur_condition%cur_time_index
    next_time_index = min(cur_condition%cur_time_index+1, &
                          cur_condition%max_time_index)

    ! ensure that condition has started
    if (time >= cur_condition%times(cur_time_index) .or. &
        abs(time-cur_condition%times(cur_time_index)) < 1.d-40) then

      ! find appropriate time interval
      do
        if (time < cur_condition%times(next_time_index) .or. &
            cur_time_index == next_time_index) &
          exit
        cur_time_index = next_time_index
        ! ensure that time index does not go beyond end of array
        if (next_time_index < cur_condition%max_time_index) &
          next_time_index = next_time_index + 1
      enddo

      ! interpolate value based on time
      cur_condition%cur_time_index = cur_time_index
      if (cur_time_index < cur_condition%max_time_index) then
        time_fraction = (time-cur_condition%times(cur_time_index)) / &
                        (cur_condition%times(next_time_index) - &
                         cur_condition%times(cur_time_index))
        cur_condition%cur_value(1:cur_condition%num_values) = &
          cur_condition%values(1:cur_condition%num_values,cur_time_index) + &
          time_fraction * &
          (cur_condition%values(1:cur_condition%num_values,next_time_index) - &
           cur_condition%values(1:cur_condition%num_values,cur_time_index))
      else
        cur_condition%cur_value(1:cur_condition%num_values) =  &
                             cur_condition%values(1:cur_condition%num_values, &
                                                  cur_condition%max_time_index)
      endif 
    endif 

    cur_condition => cur_condition%next
    
  enddo

end subroutine updateConditions

! ************************************************************************** !
!
! destroyConditionList: Deallocates the module global list and array of conditions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine destroyConditionList

  implicit none
  
  type(condition), pointer :: cur_condition, prev_condition
  
  deallocate(condition_array)
  nullify(condition_array)
  
  cur_condition => condition_list
  do 
    if (.not.associated(cur_condition)) exit
    prev_condition => cur_condition
    cur_condition => cur_condition%next
    deallocate(prev_condition)
  enddo
  
  nullify(condition_list)
  nullify(last_condition)
  num_conditions = 0

end subroutine destroyConditionList

end module Condition_module
