module Process_Model_Coupler_module

  use Process_Model_Base_class
  use Timestepper_module
  use Option_module
  use Waypoint_module

  implicit none

#include "definitions.h"
  
  private

  type, public :: process_model_coupler_type
    type(option_type), pointer :: option
    type(stepper_type), pointer :: timestepper
    class(process_model_base_type), pointer :: process_model_list
    type(waypoint_list_type), pointer :: waypoints
    type(process_model_coupler_type), pointer :: below
    type(process_model_coupler_type), pointer :: next
  contains
    procedure, public :: InitializeRun
    procedure, public :: RunTo => RunToTime
    procedure, public :: FinalizeRun
    procedure, public :: Output
    procedure, public :: UpdateSolution => PMCUpdateSolution
    procedure, public :: Destroy
  end type process_model_coupler_type
  
  public :: ProcessModelCouplerCreate
  
contains

! ************************************************************************** !
!
! ProcessModelCouplerCreate: Allocates and initializes a new 
!                            process_model_coupler object.
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function ProcessModelCouplerCreate()

  implicit none
  
  type(process_model_coupler_type), pointer :: ProcessModelCouplerCreate
  
  type(process_model_coupler_type), pointer :: process_model_coupler

  print *, 'ProcessModelCoupler%Create()'
  
  allocate(process_model_coupler)
  nullify(process_model_coupler%option)
  nullify(process_model_coupler%timestepper)
  nullify(process_model_coupler%process_model_list)
  nullify(process_model_coupler%waypoints)
  nullify(process_model_coupler%below)
  nullify(process_model_coupler%next)
  
  ProcessModelCouplerCreate => process_model_coupler  
  
end function ProcessModelCouplerCreate

! ************************************************************************** !
!
! ProcModelCouplerSetTimestepper: 
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine ProcModelCouplerSetTimestepper(this,timestepper)

  use Timestepper_module 
  
  implicit none
  
  class(process_model_coupler_type) :: this
  type(stepper_type), pointer :: timestepper

  call printMsg(this%option,'ProcessModelCoupler%SetTimestepper()')
  
  this%timestepper => timestepper
  
end subroutine ProcModelCouplerSetTimestepper

! ************************************************************************** !
!
! InitializeRun: Initializes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine InitializeRun(this)

  implicit none
  
  class(process_model_coupler_type) :: this
  
  class(process_model_base_type), pointer :: cur_process_model
  
  call printMsg(this%option,'ProcessModelCoupler%InitializeRun()')
  
  cur_process_model => this%process_model_list
  do
    if (.not.associated(cur_process_model)) exit
    call cur_process_model%InitializeRun()
    cur_process_model => cur_process_model%next
  enddo
  
  if (associated(this%below)) then
    call this%below%InitializeRun()
  endif
  
end subroutine InitializeRun

! ************************************************************************** !
!
! RunToTime: Runs the actual simulation.
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine RunToTime(this,sync_time)

  use Timestepper_module
  
  implicit none
  
  class(process_model_coupler_type) :: this
  PetscReal :: sync_time
  
  PetscBool :: failure
  class(process_model_base_type), pointer :: cur_process_model
  
  write(this%option%io_buffer,'(f12.2)') sync_time
  this%option%io_buffer = 'ProcessModelCoupler%RunToTime(' // &
    trim(adjustl(this%option%io_buffer)) // ')'
  call printMsg(this%option)
  
  do
    if (this%timestepper%target_time >= sync_time) exit
    call StepperSetTargetTime(this%timestepper,sync_time,this%option)
    call StepperStepDT(this%timestepper,this%process_model_list,failure)
    cur_process_model => this%process_model_list
    do
      if (.not.associated(cur_process_model)) exit
      call StepperUpdateDT(this%timestepper,cur_process_model)
      cur_process_model => cur_process_model%next
    enddo
    if (associated(this%below)) then
      call this%below%RunTo(this%timestepper%target_time)
    endif
  enddo
  
end subroutine RunToTime

! ************************************************************************** !
!
! PMCUpdateSolution: 
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCUpdateSolution(this)

  implicit none
  
  class(process_model_coupler_type) :: this

  call printMsg(this%option,'ProcessModelCoupler%UpdateSolution()')
  
end subroutine PMCUpdateSolution

! ************************************************************************** !
!
! FinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine FinalizeRun(this)

  implicit none
  
  class(process_model_coupler_type) :: this
  
  call printMsg(this%option,'ProcessModelCoupler%FinalizeRun()')
  
end subroutine FinalizeRun

! ************************************************************************** !
!
! Output: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine Output(this)

  implicit none
  
  class(process_model_coupler_type) :: this
  
  class(process_model_base_type), pointer :: cur_process_model
  
  call printMsg(this%option,'ProcessModelCoupler%Output()')
  
  cur_process_model => this%process_model_list
  do
    if (.not.associated(cur_process_model)) exit
!    call Output(cur_process_model%realization,plot_flag,transient_plot_flag)
    cur_process_model => cur_process_model%next
  enddo
    
end subroutine Output

! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
recursive subroutine Destroy(this)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(process_model_coupler_type) :: this
  
  call printMsg(this%option,'ProcessModelCoupler%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
  if (associated(this%process_model_list)) then
    call this%process_model_list%Destroy()
  endif

!  deallocate(process_model_coupler)
!  nullify(process_model_coupler)
  
end subroutine Destroy
  
end module Process_Model_Coupler_module
