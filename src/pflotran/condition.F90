module Condition_module

 use pflow_gridtype_module

 implicit none
 
 public
 save

#include "definitions.h"

  type condition_type
    real*8, allocatable :: times(:)
    real*8, allocatable :: values(:)
    real*8 :: cur_value
    real*8 :: datum
    integer :: id
    integer :: itype
    character(len=MAXWORDLENGTH) :: ctype
    integer :: cur_time_index, max_time_index
    type(condition_type), pointer :: next
  end type condition_type

  type condition_type_ptr
    type(condition_type), pointer :: ptr
  end type condition_type_ptr

  type(condition_type), pointer :: condition_list
  type(condition_type_ptr), allocatable :: condition_array(:)

  integer, private :: num_conditions

  type(condition_type), pointer, private :: last_condition

  public InitializeConditionList, AddCondition, UpdateConditions, &
         UpdateBoundaryConditions
  
contains
  
! ************************************************************************** !
!
! InitializeConditionList: Initializes the condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine InitializeConditionList()

  implicit none

  nullify(condition_list)
  nullify(last_condition)
  num_conditions = 0

end subroutine InitializeConditionList

! ************************************************************************** !
!
! AddCondition: Add condition to condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine AddCondition(new_condition)

  implicit none

  type(condition_type), pointer :: new_condition

  new_condition%cur_time_index = 1

  num_conditions = num_conditions + 1
  new_condition%id =  num_conditions
  if (.not.associated(condition_list)) condition_list => new_condition
  if (associated(last_condition)) last_condition%next => new_condition
  last_condition => new_condition

end subroutine AddCondition

! ************************************************************************** !
!
! ConvertListToArray: Converts the linked list of conditions to an array
!                     of conditions
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine ConvertListToArray()

  implicit none

  type(condition_type), pointer :: cur_condition

  allocate(condition_array(num_conditions))

  cur_condition => condition_list
  do
    if (.not.associated(cur_condition)) exit
    condition_array(cur_condition%id)%ptr => cur_condition
    cur_condition => cur_condition%next
  enddo

end subroutine ConvertListToArray

! ************************************************************************** !
!
! UpdateConditions: Updates the value of the boundary/initial condition based
!                   on the current time
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine UpdateConditions(time)

  implicit none

  real*8 :: time

  type(condition_type), pointer :: cur_condition

  integer :: cur_time_index, next_time_index
  real*8 :: time_fraction

  cur_condition => condition_list
  do

    if (.not.associated(cur_condition)) exit

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
      if (cur_time_index < cur_condition%max_time_index) then
        time_fraction = (time-cur_condition%times(cur_time_index)) / &
                        (cur_condition%times(next_time_index) - &
                         cur_condition%times(cur_time_index))
        cur_condition%cur_value = cur_condition%values(cur_time_index) + &
                                  time_fraction * &
                                  (cur_condition%values(next_time_index) - &
                                   cur_condition%values(cur_time_index))
      else
        cur_condition%cur_value =  &
                             cur_condition%values(cur_condition%max_time_index) 
      endif 
    endif 

    cur_condition => cur_condition%next
    
  enddo

end subroutine UpdateConditions

! ************************************************************************** !
!
! InitializeBoundaryConditions: Initializes boundary condition type
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine InitializeBoundaryConditions(grid)

  use pflow_gridtype_module
  use fileio_module

  implicit none

  type(pflowGrid) :: grid

  integer :: icond
  character(len=MAXWORDLENGTH) :: word
  type(condition_type), pointer :: cur_condition

  call ConvertListToArray()

  do icond = 1, num_conditions
    word = condition_array(icond)%ptr%ctype 
    call fiWordToLower(word)
    if (fiStringCompare("dirichlet",word,9)) then
      condition_array(icond)%ptr%itype = 1  ! integer pointer to condition
    else if (fiStringCompare("neumann",word,7)) then
      condition_array(icond)%ptr%itype = 2
    else if (fiStringCompare("hydraulic gradient",word,18)) then
      condition_array(icond)%ptr%itype = 3
    else if (fiStringCompare("seepage face",word,12)) then
      condition_array(icond)%ptr%itype = 4
    else
      if (grid%myrank == 0) print *, 'Condition type: ', word, ' not supported'
      stop
    endif

  enddo

  ! set up native pflotran integer pointer to boundary condition type
  allocate(grid%ibndtyp(num_conditions)) ! condition (boundary) type
  do icond = 1, num_conditions
    if (condition_array(icond)%ptr%itype == 2) then
      grid%ibndtyp(icond) = 2  ! neumann
    else
      grid%ibndtyp(icond) = 3  ! pressure based dirichlet
    endif
  enddo

  ! updated to initial boundary condition values
  call UpdateBoundaryConditions(grid)

end subroutine InitializeBoundaryConditions

! ************************************************************************** !
!
! UpdateBoundaryConditions: Updates the value of the boundary condition based
!                           on the current time
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine UpdateBoundaryConditions(grid)

  use pflow_gridtype_module

  implicit none

  type(pflowGrid) :: grid

  integer :: iconnbc, icond, itype
  real*8 :: value, datum, delz, delp
  real*8 :: patm = 101325.d0

  call UpdateConditions(grid%t)

  do iconnbc = 1, grid%nconnbc
    icond = grid%ibconn(iconnbc)
    itype = condition_array(icond)%ptr%itype
    datum = condition_array(icond)%ptr%datum
    value = condition_array(icond)%ptr%cur_value
    if (itype == 3 .or. itype == 4) then ! correct for pressure gradient
      delz = grid%z(grid%nL2A(grid%mblkbc(iconnbc))+1)-datum
      delp = delz*9.81d0*998.32
      value = value - delp
    endif
    if (itype == 2) then ! neumann
      grid%velocitybc(1:grid%nphase,iconnbc) = value
      grid%xxbc(1:grid%nphase,iconnbc) = patm
    else
!      grid%pressurebc(1:grid%nphase,iconnbc) = value
      grid%xxbc(1:grid%nphase,iconnbc) = value
    endif
  enddo

end subroutine UpdateBoundaryConditions

end module Condition_module
