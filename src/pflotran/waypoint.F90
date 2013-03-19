module Waypoint_module
 
  use Option_module
  use Output_Aux_module
  
  implicit none
  
  private

#include "definitions.h"

  ! linked-list for waypoints in the simulation
  type, public :: waypoint_type
    PetscReal :: time
    PetscBool :: print_output
    PetscBool :: print_tr_output
!    type(output_option_type), pointer :: output_option
    PetscBool :: update_conditions
    PetscReal :: dt_max
    PetscBool :: final  ! any waypoint after this will be deleted
    type(waypoint_type), pointer :: prev
    type(waypoint_type), pointer :: next
  end type waypoint_type
  
  type, public :: waypoint_list_type
    PetscInt :: num_waypoints
    type(waypoint_type), pointer :: first
    type(waypoint_type), pointer :: last
    type(waypoint_type), pointer :: array(:)    
  end type waypoint_list_type
  
  interface WaypointCreate
    module procedure WaypointCreate1
    module procedure WaypointCreate2
  end interface  
  
  public :: WaypointCreate, &
            WaypointListCreate, &
            WaypointInsertInList, &
            WaypointDeleteFromList, &
            WaypointListFillIn, &
            WaypointListRemoveExtraWaypnts, &
            WaypointConvertTimes, &
            WaypointSkipToTime, &
            WaypointListPrint

contains

! ************************************************************************** !
!
! WaypointCreate1: Creates a simulation waypoint
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function WaypointCreate1()

  implicit none
  
  type(waypoint_type), pointer :: WaypointCreate1
  
  type(waypoint_type), pointer :: waypoint
  
  allocate(waypoint)
  waypoint%time = 0.d0
  waypoint%print_output = PETSC_FALSE
  waypoint%print_tr_output = PETSC_FALSE
  waypoint%final = PETSC_FALSE
  waypoint%update_conditions = PETSC_FALSE
  waypoint%dt_max = 0.d0
  nullify(waypoint%next)
  nullify(waypoint%prev)
    
  WaypointCreate1 => waypoint
  
end function WaypointCreate1

! ************************************************************************** !
!
! WaypointCreate2: Creates a simulation waypoint
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function WaypointCreate2(original_waypoint)

  implicit none
  
  type(waypoint_type), pointer :: original_waypoint
  
  type(waypoint_type), pointer :: WaypointCreate2
  
  type(waypoint_type), pointer :: waypoint
  
  allocate(waypoint)
  waypoint%time = original_waypoint%time
  waypoint%print_output = original_waypoint%print_output
  waypoint%print_tr_output = original_waypoint%print_tr_output
  waypoint%final = original_waypoint%final
  waypoint%update_conditions = original_waypoint%update_conditions
  waypoint%dt_max = original_waypoint%dt_max
    
  WaypointCreate2 => waypoint
  
end function WaypointCreate2

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

  use Utility_module

  type(waypoint_type), pointer :: new_waypoint
  type(waypoint_list_type) :: waypoint_list

  type(waypoint_type), pointer :: waypoint

!    PetscReal :: time
!    PetscBool :: print_output
!    type(output_option_type), pointer :: output_option
!    PetscBool :: update_bcs
!    PetscBool :: update_srcs
!    PetscReal :: dt_max
    
    ! place new waypoint in proper location within list
  waypoint => waypoint_list%first
  if (associated(waypoint)) then ! list exists
    ! if waypoint time matches another waypoint time, merge them
!geh    if ((new_waypoint%time > 0.999999d0*waypoint%time .and. &
!geh         new_waypoint%time < 1.000001d0*waypoint%time) .or. &
         ! need to account for waypoint%time = 0.d0
    if (Equal(new_waypoint%time,waypoint%time) .or. &
        (new_waypoint%time < 1.d-40 .and. &
         waypoint%time < 1.d-40)) then ! same
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
            if (Equal(new_waypoint%time,waypoint%time)) then
!geh            if (new_waypoint%time > 0.999999d0*waypoint%time .and. &
!geh                new_waypoint%time < 1.000001d0*waypoint%time) then ! same
              call WaypointMerge(waypoint,new_waypoint)
              return
            else if (associated(waypoint%next)) then 
              if (new_waypoint%time-waypoint%time > 1.d-10 .and. & ! within list
                  new_waypoint%time-waypoint%next%time < -1.d-10) then 
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
! WaypointDeleteFromList: Deletes a waypoing in a list
! author: Gautam Bisht
! date: 01/20/11
!
! ************************************************************************** !
subroutine WaypointDeleteFromList(obsolete_waypoint,waypoint_list)

  implicit none

  type(waypoint_type), pointer :: obsolete_waypoint
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  type(waypoint_list_type)     :: waypoint_list

  waypoint => waypoint_list%first

  if (associated(waypoint)) then ! list exists

    ! Is the waypoint to be deleted is the first waypoint?
    if (waypoint%time == obsolete_waypoint%time) then
      waypoint_list%first => waypoint%next
      call WaypointDestroy(waypoint)
      return
    else

      prev_waypoint => waypoint
      waypoint => waypoint%next
      do
        if (associated(waypoint)) then
          if (dabs(waypoint%time-obsolete_waypoint%time) < 1.d-10) then
            prev_waypoint%next => waypoint%next
            call WaypointDestroy(waypoint)
            return
          endif
          prev_waypoint => waypoint
          waypoint => waypoint%next
          cycle
        else
         ! at the end of the list, didn't find obsolete waypoint
          return
        endif
      enddo
    endif
  else
    ! list does not exists
    return
  endif
  
end subroutine WaypointDeleteFromList

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
  PetscReal :: dt_max = -999.d0
  
  ! find first value of dt_max > 0.d0 in list
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    if (waypoint%dt_max > 1.d-40) then
      dt_max = waypoint%dt_max
      exit
    endif
    waypoint => waypoint%next
  enddo

  if (dt_max <= 1.d-40) then
    option%io_buffer = 'All values of dt_max in input file uninitialized'
    call printErrMsg(option)
  endif
  
  ! assign that value to the first waypoint, if waypoint%dt_max not already > 1.d-40
  waypoint => waypoint_list%first
  if (waypoint%dt_max < 1.d-40) waypoint%dt_max = dt_max
  
  ! fill in missing values
  do
    prev_waypoint => waypoint
    waypoint => waypoint%next
    if (.not.associated(waypoint)) exit 
    if (waypoint%dt_max < 1.d-40) then
      waypoint%dt_max = prev_waypoint%dt_max
    endif
  enddo
  
  ! IMPORTANT NOTE:  The dt_max must be assigned to the "next" waypoint.  The
  ! "current" waypoint in the stepper is always the next waypoint .  Therefore
  ! we must shift all the dt_max entries. 
  waypoint => waypoint_list%last
  ! work backwards
  do
    prev_waypoint => waypoint%prev
    if (.not.associated(prev_waypoint)) exit 
    waypoint%dt_max = prev_waypoint%dt_max
    waypoint => prev_waypoint
  enddo
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit 
    waypoint => waypoint%next
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
  PetscReal :: time_conversion
  
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
  
  type(option_type) :: option
  type(waypoint_list_type) :: waypoint_list
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  
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
    write(option%io_buffer,'("Waypoint at time:", 1pe12.4, &
  &       " is beyond the end of simulation")') &
          prev_waypoint%time
    call printWrnMsg(option)
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

!    PetscReal :: time
!    PetscBool :: print_output
!    type(output_option_type), pointer :: output_option
!    PetscBool :: update_bcs
!    PetscBool :: update_srcs
!    PetscReal :: dt_max
!    PetscBool :: final  ! any waypoint after this will be deleted
    
  if (old_waypoint%print_output .or. new_waypoint%print_output) then
    old_waypoint%print_output = PETSC_TRUE
  else
    old_waypoint%print_output = PETSC_FALSE
  endif

  if (old_waypoint%print_tr_output .or. new_waypoint%print_tr_output) then
    old_waypoint%print_tr_output = PETSC_TRUE
  else
    old_waypoint%print_tr_output = PETSC_FALSE
  endif

  if (old_waypoint%update_conditions .or. new_waypoint%update_conditions) then
    old_waypoint%update_conditions = PETSC_TRUE
  else
    old_waypoint%update_conditions = PETSC_FALSE
  endif

  if (new_waypoint%dt_max > 0.d0) then
    old_waypoint%dt_max = new_waypoint%dt_max
  endif
  
  if (old_waypoint%final .or. new_waypoint%final) then
    old_waypoint%final = PETSC_TRUE
  else
    old_waypoint%final = PETSC_FALSE
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
! WaypointSkipToTime: Returns a pointer to the first waypoint after time
! author: Glenn Hammond
! date: 1/03/08
!
! ************************************************************************** !
function WaypointSkipToTime(list,time)

  implicit none

  type(waypoint_list_type), pointer :: list
  PetscReal :: time

  type(waypoint_type), pointer :: WaypointSkipToTime
  type(waypoint_type), pointer :: waypoint
  
  waypoint => list%first
  do 
    if (.not.associated(waypoint)) exit
    if (waypoint%time > time) exit
    waypoint => waypoint%next
  enddo

  if (associated(waypoint)) then
    WaypointSkipToTime => waypoint
  else
    nullify(WaypointSkipToTime)
  endif

end function WaypointSkipToTime

! ************************************************************************** !
!
! WaypointListPrint: Prints a waypoint
! author: Glenn Hammond
! date: 05/20/11
!
! ************************************************************************** !
subroutine WaypointListPrint(list,option,output_option)

  use Option_module
  use Output_Aux_module

  implicit none
  
  type(waypoint_list_type), pointer :: list
  type(option_type) :: option
  type(output_option_type) :: output_option

  type(waypoint_type), pointer :: cur_waypoint
  PetscInt :: icount

  100 format(/)
  110 format(a)
  20 format('  ',a20,':',10i6)

  if (OptionPrintToScreen(option)) then
    write(*,100)
    write(*,110) 'List of Waypoints:'
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,100)
    write(option%fid_out,110) 'List of Waypoints:'
    write(option%fid_out,100)
  endif

  icount = 0
  cur_waypoint => list%first
  do 
    if (.not.associated(cur_waypoint)) exit
    call WaypointPrint(cur_waypoint,option,output_option)
    icount = icount + 1
    cur_waypoint => cur_waypoint%next
  enddo

  if (OptionPrintToScreen(option)) then
    write(option%fid_out,20) 'Total Waypoints:', icount
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,20) 'Total Waypoints:', icount
    write(option%fid_out,100)
  endif

end subroutine WaypointListPrint

! ************************************************************************** !
!
! WaypointListCopy: Copies a waypoint list
! author: Glenn Hammond
! date: 03/19/13
!
! ************************************************************************** !
function WaypointListCopy(list,option,output_option)

  use Option_module
  use Output_Aux_module

  implicit none
  
  type(waypoint_list_type), pointer :: new_list
  
  type(waypoint_list_type), pointer :: list
  type(waypoint_type), pointer :: new_waypoint
  type(waypoint_type), pointer :: prev_new_waypoint
  
  new_list => WaypointListCreate()
  
  nullify(prev_new_waypoint)
  
  cur_waypoint => list%first
  do 
    if (.not.associated(cur_waypoint)) exit
    new_waypoint => WaypointCreate(cur_waypoint)
    if (associated(prev_new_waypoint)) then
      prev_new_waypoint%next => new_waypoint
    else
      new_list%first => new_waypoint
    endif
    new_list%num_waypoints = new_list%num_waypoints + 1
    prev_new_waypoint => new_waypoint
    nullify(new_waypoint)
    cur_waypoint => cur_waypoint%next
  enddo

end function WaypointListCopy

! ************************************************************************** !
!
! WaypointPrint: Prints a waypoint
! author: Glenn Hammond
! date: 05/20/11
!
! ************************************************************************** !
subroutine WaypointPrint(waypoint,option,output_option)

  use Option_module

  implicit none
  
  type(waypoint_type), pointer :: waypoint
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXSTRINGLENGTH) :: string

  10 format('  ',a20,':',10es13.5)
  20 format('  ',a20,':',10i6)
  30 format('  ',a20,':',10l)
  40 format('  ',a20,':',a20)
  100 format(/)
  110 format(a)

  if (OptionPrintToScreen(option)) then
    write(*,110) 'Waypoint:'
    write(string,*) 'Time [' // trim(adjustl(output_option%tunit)) // ']'
    write(*,10) trim(string), waypoint%time/output_option%tconv
    write(*,30) 'Print Output', waypoint%print_output
    write(*,30) 'Print Tr. Output', waypoint%print_tr_output
    write(*,30) 'Update Conditions', waypoint%update_conditions
    write(*,30) 'Print Output', waypoint%print_output
    write(string,*) 'Max DT [' // trim(adjustl(output_option%tunit)) // ']'
    write(*,10) trim(string), waypoint%dt_max/output_option%tconv
    write(*,30) 'Final', waypoint%final
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,110) 'Waypoint:'
    write(string,*) 'Time [' // trim(adjustl(output_option%tunit)) // ']'
    write(option%fid_out,10) trim(string), waypoint%time/output_option%tconv
    write(option%fid_out,30) 'Print Output', waypoint%print_output
    write(option%fid_out,30) 'Print Tr. Output', waypoint%print_tr_output
    write(option%fid_out,30) 'Update Conditions', waypoint%update_conditions
    write(option%fid_out,30) 'Print Output', waypoint%print_output
    write(string,*) 'Max DT [' // trim(adjustl(output_option%tunit)) // ']'
    write(option%fid_out,10) trim(string), waypoint%dt_max/output_option%tconv
    write(option%fid_out,30) 'Final', waypoint%final
    write(option%fid_out,100)
  endif
 
end subroutine WaypointPrint

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
