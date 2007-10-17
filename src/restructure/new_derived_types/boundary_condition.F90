module Boundary_Condition_module

 use Condition_module
 use Region_module

 implicit none
 
#include "definitions.h"

  private

  type, public :: boundary_condition
    integer :: id
    character(len=MAXWORDLENGTH) :: condition_name
    character(len=MAXWORDLENGTH) :: region_name
    integer :: icondition
    integer :: iregion
    type(condition), pointer :: condition_ptr
    type(region), pointer :: region_ptr
    type(boundary_condition), pointer :: next
  end type boundary_condition

  ! pointer data structure required for making an array of boundary condition pointers in F90
  type :: boundary_condition_ptr
    type(boundary_condition), pointer :: ptr
  end type boundary_condition_ptr

  ! linked list of boundary conditions
  type(boundary_condition), pointer :: boundary_condition_list
  ! pointer to last boundary condition; facilitates appendage of boundary conditions
  type(boundary_condition), pointer, private :: last_boundary_condition
  ! array of pointers to boundary conditions in list
  type(boundary_condition_ptr), pointer :: boundary_condition_array(:)
  ! number of regions in list
  integer, private :: num_boundary_conditions


  public initBoundaryConditionModule, readBoundaryCondition
  
contains
  
! ************************************************************************** !
!
! initBoundaryConditionModule: Initializes the boundary condition list
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine initBoundaryConditionModule()

  implicit none

  nullify(boundary_condition_list)
  nullify(last_boundary_condition)
  num_boundary_conditions = 0

end subroutine initBoundaryConditionModule

! ************************************************************************** !
!
! readBoundaryCondition: Reads a boundary condition from input
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine readBoundaryCondition(string)

  use fileio_module
  
  implicit none
  
  character(len=*) :: string
  type(boundary_condition), pointer :: new_boundary_condition
  integer :: ierr
  
  new_boundary_condition = createBoundaryCondition()
  
  call fiReadWord(string,new_boundary_condition%region_name,.true.,ierr)
  call fiErrorMsg('region','BCON',ierr)

  call fiReadWord(string,new_boundary_condition%condition_name,.true.,ierr)
  call fiErrorMsg('condition','BCON',ierr)

end subroutine readBoundaryCondition

! ************************************************************************** !
!
! mapBoundaryConditions: Maps boundary condition regions and conditions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine mapBoundaryConditions()

  type(boundary_condition), pointer :: cur_boundary_condition

  cur_boundary_condition => boundary_condition_list
  do
    if (.not.associated(cur_boundary_condition)) exit
    ! map condition
    cur_boundary_condition%condition_ptr => &
                 getConditionFromName(cur_boundary_condition%condition_name)
    cur_boundary_condition%icondition = cur_boundary_condition%condition_ptr%id
    ! map region
    cur_boundary_condition%region_ptr => &
                 getRegionFromName(cur_boundary_condition%region_name)
    cur_boundary_condition%iregion = cur_boundary_condition%region_ptr%id
    cur_boundary_condition => cur_boundary_condition%next
  enddo
end subroutine mapBoundaryConditions

! ************************************************************************** !
!
! addBoundaryConditionToList: add boundary condition to boundary condition list
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine addBoundaryConditionToList(new_boundary_condition)

  implicit none

  type(boundary_condition), pointer :: new_boundary_condition

  num_boundary_conditions = num_boundary_conditions + 1
  new_boundary_condition%id =  num_boundary_conditions
  if (.not.associated(boundary_condition_list)) &
    boundary_condition_list => new_boundary_condition
  if (associated(last_boundary_condition)) &
    last_boundary_condition%next => new_boundary_condition
  last_boundary_condition => new_boundary_condition

end subroutine addBoundaryConditionToList

! ************************************************************************** !
!
! createBoundaryCondition: creates a new boundary condition
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function createBoundaryCondition()

  implicit none

  type(boundary_condition), pointer :: createBoundaryCondition

  allocate(createBoundaryCondition)
  createBoundaryCondition%icondition = 0
  createBoundaryCondition%iregion = 0
  nullify(createBoundaryCondition%condition_ptr)
  nullify(createBoundaryCondition%region_ptr)

end function createBoundaryCondition

! ************************************************************************** !
!
! convertBoundaryConditionListToArray: Converts the linked list of conditions to an 
!                              array of conditions
! author: Glenn Hammond
! date: 03/15/07
!
! ************************************************************************** !
subroutine convertBoundaryConditionListToArray()

  implicit none

  type(boundary_condition), pointer :: cur_boundary_condition

  allocate(boundary_condition_array(num_boundary_conditions))

  cur_boundary_condition => boundary_condition_list
  do
    if (.not.associated(cur_boundary_condition)) exit
    boundary_condition_array(cur_boundary_condition%id)%ptr => &
                                                cur_boundary_condition
    cur_boundary_condition => cur_boundary_condition%next
  enddo

end subroutine convertBoundaryConditionListToArray

! ************************************************************************** !
!
! destroyBoundaryConditionList: Deallocates the module global list and 
!                               array of conditions
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine destroyBoundaryConditionList

  implicit none
  
  type(boundary_condition), pointer :: cur_boundary_condition, &
                                       prev_boundary_condition
  
  deallocate(boundary_condition_array)
  nullify(boundary_condition_array)
  
  cur_boundary_condition => boundary_condition_list
  do 
    if (.not.associated(cur_boundary_condition)) exit
    prev_boundary_condition => cur_boundary_condition
    cur_boundary_condition => cur_boundary_condition%next
    deallocate(prev_boundary_condition)
  enddo
  
  nullify(boundary_condition_list)
  nullify(last_boundary_condition)
  num_boundary_conditions = 0

end subroutine destroyBoundaryConditionList

end module Boundary_Condition_module
