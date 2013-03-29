module Synchronizer_module

  use Process_Model_Coupler_module
  use Timestepper_module
  use Option_module
  use Output_Aux_module
  use Waypoint_module

  implicit none

#include "definitions.h"
  
  private

  type, public :: synchronizer_type
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    type(waypoint_list_type), pointer :: waypoints
    type(process_model_coupler_type), pointer :: process_model_coupler_list
  contains
    procedure, public :: ExecuteRun => SynchronizerExecuteRun
  end type synchronizer_type
  
  public :: SynchronizerCreate, &
            SynchronizerDestroy
            
  interface SynchronizerRunToTime
    module procedure RunToTime1
    module procedure RunToTime2
  end interface
  
contains

! ************************************************************************** !
!
! SynchronizerCreate: Allocates and initializes a new synchronizer object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function SynchronizerCreate()

  implicit none
  
  type(synchronizer_type), pointer :: SynchronizerCreate
  
  type(synchronizer_type), pointer :: synchronizer
  
  print *, 'SynchronizerCreate'
  
  allocate(synchronizer)
  nullify(synchronizer%option)
  nullify(synchronizer%output_option)
  nullify(synchronizer%waypoints)
  nullify(synchronizer%process_model_coupler_list)
  
  SynchronizerCreate => synchronizer  
  
end function SynchronizerCreate

! ************************************************************************** !
!
! SynchronizerExecuteRun: Executes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine SynchronizerExecuteRun(this,stop_flag)

  implicit none
  
  class(synchronizer_type) :: this
  PetscInt :: stop_flag
  
  type(waypoint_type), pointer :: cur_waypoint
  
  call printMsg(this%option,'Synchronizer%ExecuteRun()')
  
  cur_waypoint => this%waypoints%first
  do
    if (.not.associated(cur_waypoint)) exit
    call SynchronizerRunToTime(this,cur_waypoint%time,stop_flag)
    ! run to waypoint
 !   call this%RunTo(cur_waypoint%time)
    ! output plot files
 !   if (associated(cur_waypoint%print_output)) then
 !     call this%Output()
 !   endif
    cur_waypoint => cur_waypoint%next

    if (this%option%wallclock_stop_flag) then
    !TODO(geh): reformulate 
#if 0    
      call PetscGetTime(current_time, ierr)
      average_step_time = (current_time-master_stepper%start_time)/ &
                          real(master_stepper%steps-&
                               master_stepper%start_time_step+1) &
                          *2.d0  ! just to be safe, double it
      if (average_step_time + current_time > option%wallclock_stop_time) then
        call printMsg(option,"Wallclock stop time exceeded.  Exiting!!!")
        call printMsg(option,"")
        stop_flag = 1
        return
      endif
#endif    
    endif
    
  enddo
  
end subroutine SynchronizerExecuteRun

! ************************************************************************** !
!
! RunToTime1: Executes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine RunToTime1(this,target_time,stop_flag)

  implicit none
  
  type(synchronizer_type) :: this
  PetscReal :: target_time
  PetscInt :: stop_flag
  
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  
  call printMsg(this%option,'Synchronizer%RunToTime1()')
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%RunTo(target_time,stop_flag)
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

end subroutine RunToTime1

! ************************************************************************** !
!
! RunToTime2: Executes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine RunToTime2(this,process_model_coupler_list,target_time,stop_flag)

  implicit none
  
  type(synchronizer_type) :: this
  type(process_model_coupler_type), pointer :: process_model_coupler_list
  PetscReal :: target_time
  PetscInt :: stop_flag
  
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  
  call printMsg(this%option,'Synchronizer%RunToTime2()')
  
  cur_process_model_coupler => process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%RunTo(target_time,stop_flag)
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

end subroutine RunToTime2

! ************************************************************************** !
!
! Output: Executes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine Output(this)

  implicit none
  
  type(synchronizer_type) :: this
  
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  
  call printMsg(this%option,'Synchronizer%Output()')
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%Output()!PETSC_TRUE,PETSC_FALSE)
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

end subroutine Output

! ************************************************************************** !
!
! SynchronizerDestroy: Deallocates an synchronizer object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine SynchronizerDestroy(synchronizer)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  type(synchronizer_type), pointer :: synchronizer
  
  call printMsg(synchronizer%option,'SynchronizerDestroy()')
  
  if (.not.associated(synchronizer)) return

  deallocate(synchronizer)
  nullify(synchronizer)
  
end subroutine SynchronizerDestroy
  
end module Synchronizer_module
