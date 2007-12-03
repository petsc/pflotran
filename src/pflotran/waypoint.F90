module Waypoint_module
 
  use Option_module
  
  implicit none

  private

  ! linked-list for waypoints in the simulation
  type, public :: waypoint_type
    real*8 :: time
    logical :: print_output
    type(output_option_type), pointer :: output_option
    logical :: update_bcs
    logical :: update_srcs
    real*8 :: dt_max
    logical :: final  ! any waypoint after this will be deleted
    type(waypoint_type), pointer :: prev
    type(waypoint_type), pointer :: next
  end type waypoint_type
  
  type, public :: waypoint_list_type
    integer :: num_waypoints
    type(waypoint_type), pointer :: first
    type(waypoint_type), pointer :: last
    type(waypoint_type), pointer :: array(:)    
  end type waypoint_list_type
  
  
  public :: WaypointCreate, &
            WaypointListCreate, &
            WaypointInsertInList, &
            WaypointListFillIn, &
            WaypointListRemoveExtraWaypnts, &
            WaypointConvertTimes

contains

! ************************************************************************** !
!
! WaypointCreate: Creates a simulation waypoint
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function WaypointCreate()

  implicit none
  
  type(waypoint_type), pointer :: WaypointCreate
  
  type(waypoint_type), pointer :: waypoint
  
  allocate(waypoint)
  waypoint%time = 0.d0
  waypoint%print_output = .false.
  waypoint%final = .false.
!  waypoint%output_option => OutputOptionCreate()
  nullify(waypoint%output_option)
  waypoint%update_bcs = .false.
  waypoint%update_srcs = .false.
  waypoint%dt_max = 0.d0
  nullify(waypoint%next)
  nullify(waypoint%prev)
    
  WaypointCreate => waypoint
  
end function WaypointCreate 

! ************************************************************************** !
!
! WaypointListCreate: Creates a simulation waypoint list
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function WaypointListCreate()

  implicit none
  
  type(waypoint_list_type), pointer :: WaypointListCreate
  
  type(waypoint_list_type), pointer :: waypoint_list
  
  allocate(waypoint_list)
  nullify(waypoint_list%first)
  nullify(waypoint_list%last)
  nullify(waypoint_list%array)
  waypoint_list%num_waypoints = 0

  WaypointListCreate => waypoint_list
  
end function WaypointListCreate 

! ************************************************************************** !
!
! WaypointInsertInList: Correctly inserts a waypoing in a list
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointInsertInList(new_waypoint,waypoint_list)

  type(waypoint_type), pointer :: new_waypoint
  type(waypoint_list_type) :: waypoint_list

  type(waypoint_type), pointer :: waypoint

!    real*8 :: time
!    logical :: print_output
!    type(output_option_type), pointer :: output_option
!    logical :: update_bcs
!    logical :: update_srcs
!    real*8 :: dt_max
    
    ! place new waypoint in proper location within list
  waypoint => waypoint_list%first
  if (associated(waypoint)) then ! list exists
    ! if waypoint time matches another waypoint time, merge them
    if (new_waypoint%time > 0.999999d0*waypoint%time .and. &
        new_waypoint%time < 1.000001d0*waypoint%time) then ! same
      call WaypointMerge(waypoint,new_waypoint)
      return
    else
      ! if waypoint time is less than any previous, insert at beginning of list
      if (new_waypoint%time < waypoint%time) then 
        waypoint_list%first => new_waypoint
        new_waypoint%next => waypoint
        new_waypoint%next%prev => new_waypoint
      else
        ! find its location in the list
        do
          if (associated(waypoint)) then 
            if (new_waypoint%time > 0.999999d0*waypoint%time .and. &
                new_waypoint%time < 1.000001d0*waypoint%time) then ! same
              call WaypointMerge(waypoint,new_waypoint)
              return
            elseif (associated(waypoint%next)) then 
              if (new_waypoint%time > waypoint%time .and. & ! within list
                  new_waypoint%time < waypoint%next%time) then 
                new_waypoint%next => waypoint%next
                new_waypoint%next%prev => new_waypoint
                waypoint%next => new_waypoint
                new_waypoint%prev => waypoint
                waypoint_list%num_waypoints = waypoint_list%num_waypoints+1
                return
              else
                waypoint => waypoint%next
                cycle
              endif
            else ! at end of list
              waypoint%next => new_waypoint
              new_waypoint%prev => waypoint
              waypoint_list%last => new_waypoint
              exit
            endif
          endif
        enddo
      endif
    endif
  else
    waypoint_list%first => new_waypoint
    waypoint_list%last => new_waypoint 
  endif
  waypoint_list%num_waypoints = waypoint_list%num_waypoints + 1

end subroutine WaypointInsertInList

! ************************************************************************** !
!
! WaypointListFillIn: Fills in missing values (e.g. dt_max) in waypoint list
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointListFillIn(option,waypoint_list)
  
  implicit none
  
  type(option_type) :: option
  type(waypoint_list_type) :: waypoint_list
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  real*8 :: dt_max = -999.d0
  
  ! find first value of dt_max > 0.d0 in list
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    prev_waypoint => waypoint
    waypoint => waypoint%next
    if (associated(waypoint) .and. waypoint%dt_max > 1.d-40) then
      dt_max = waypoint%dt_max
      exit
    endif
  enddo

  if (dt_max <= 1.d-40) then
    call printErrMsg(option,'All values of dt_max in input file uninitialized')
  endif
  
  ! assign that value to the first waypoint, if waypoint%dt_max not already > 1.d-40
  waypoint => waypoint_list%first
  if (waypoint%dt_max < 1.d-40) waypoint%dt_max = dt_max
  
  ! fill in the rest
  do
    prev_waypoint => waypoint
    waypoint => waypoint%next
    if (.not.associated(waypoint)) exit 
    if (waypoint%dt_max < 1.d-40) then
      waypoint%dt_max = prev_waypoint%dt_max
    endif
  enddo
  
end subroutine WaypointListFillIn 

! ************************************************************************** !
!
! WaypointConvertTimes: Converts time units to seconds
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointConvertTimes(waypoint_list,time_conversion)

  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  real*8 :: time_conversion
  
  type(waypoint_type), pointer :: waypoint
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    waypoint%time = waypoint%time * time_conversion
    waypoint%dt_max = waypoint%dt_max * time_conversion
    waypoint => waypoint%next
  enddo
  
end subroutine WaypointConvertTimes 

! ************************************************************************** !
!
! WaypointListRemoveExtraWaypnts: 
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointListRemoveExtraWaypnts(option,waypoint_list)

  implicit none
  
#include "definitions.h"

  type(option_type) :: option
  type(waypoint_list_type) :: waypoint_list
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  character(len=MAXSTRINGLENGTH) :: string
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint) .or. waypoint%final) exit
    waypoint => waypoint%next
  enddo
  
  if (associated(waypoint)) then
    prev_waypoint => waypoint
    waypoint => waypoint%next
    nullify(prev_waypoint%next)
  endif
  
  do
    if (.not.associated(waypoint)) exit
    prev_waypoint => waypoint
    waypoint => waypoint%next
    write(string,'("Waypoint at time:", 1pe12.4, &
  &       " is beyond the end of simulation")') &
          prev_waypoint%time
    call printWrnMsg(option,trim(string))
    call WaypointDestroy(prev_waypoint)   
  enddo

end subroutine WaypointListRemoveExtraWaypnts 

! ************************************************************************** !
!
! WaypointMerge: Merges 2 waypoints performing an OR operation on logicals
! author: Glenn Hammond
! date: 10/28/03
!
! ************************************************************************** !
subroutine WaypointMerge(old_waypoint,new_waypoint)

  implicit none

  type(waypoint_type), pointer :: old_waypoint, new_waypoint

  new_waypoint%time = 0.d0

!    real*8 :: time
!    logical :: print_output
!    type(output_option_type), pointer :: output_option
!    logical :: update_bcs
!    logical :: update_srcs
!    real*8 :: dt_max

  if (old_waypoint%print_output .or. new_waypoint%print_output) then
    old_waypoint%print_output = .true.
  else
    old_waypoint%print_output = .false.
  endif

  if (old_waypoint%update_bcs .or. new_waypoint%update_bcs) then
    old_waypoint%update_bcs = .true.
  else
    old_waypoint%update_bcs = .false.
  endif

  if (old_waypoint%update_srcs .or. new_waypoint%update_srcs) then
    old_waypoint%update_srcs = .true.
  else
    old_waypoint%update_srcs = .false.
  endif

  if (new_waypoint%dt_max > 0.d0) then
    old_waypoint%dt_max = new_waypoint%dt_max
  endif
  
  ! deallocate new waypoint
  deallocate(new_waypoint)
  ! point new_waypoint to old
  new_waypoint => old_waypoint

end subroutine WaypointMerge

! ************************************************************************** !
!
! WaypointListDestroy: Deallocates a waypoint list
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointListDestroy(waypoint_list)

  implicit none
  
  type(waypoint_list_type), pointer :: waypoint_list
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  
  if (.not.associated(waypoint_list)) return

  if (associated(waypoint_list%array)) deallocate(waypoint_list%array)
  nullify(waypoint_list%array)
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    prev_waypoint => waypoint
    waypoint => waypoint%next
    call WaypointDestroy(prev_waypoint)
  enddo
  
  nullify(waypoint_list%first)
  nullify(waypoint_list%last)
  deallocate(waypoint_list)
  nullify(waypoint_list)
  
end subroutine WaypointListDestroy

! ************************************************************************** !
!
! WaypointDestroy: Deallocates a waypoint
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine WaypointDestroy(waypoint)

  implicit none
  
  type(waypoint_type), pointer :: waypoint
  
  if (.not.associated(waypoint)) return

  nullify(waypoint%prev)
  nullify(waypoint%next)
  deallocate(waypoint)
  nullify(waypoint)
  
end subroutine WaypointDestroy

end module Waypoint_module
