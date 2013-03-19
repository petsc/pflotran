module Process_Model_Coupler_module

  use Process_Model_Base_class
  use Timestepper_module

  implicit none

#include "definitions.h"
  
  private

  type, public :: process_model_coupler_type
    type(stepper_type), pointer :: timestepper
    class(process_model_base_type), pointer :: process_model_list
    type(process_model_coupler_type), pointer :: next
  contains
    procedure, public :: InitializeRun
    procedure, public :: RunTo
    procedure, public :: FinalizeRun
    procedure, public :: Output
    procedure, public :: UpdateSolution => PMCUpdateSolution
    procedure, public :: Destroy
  end type process_model_coupler_type
  
  public :: ProcessModelCouplerCreate, &
            ProcessModelCouplerDestroy
  
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
  
  allocate(process_model_coupler)
  nullify(process_model_coupler%timestepper)
  nullify(process_model_coupler%process_model_list)
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
  
  class(process_model_base_type) :: this
  type(stepper_type), pointer :: timestepper

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
  
  type(process_model_coupler_type), pointer :: this
  
  type(process_model_base_type), pointer :: cur_process_model
  
  cur_process_model => this%process_model_list
  do
    if (.not.associated(cur_process_model)) exit
    call cur_process_model%InitializeRun()
    cur_process_model => cur_process_model%next
  enddo  
  
end subroutine InitializeRun

! ************************************************************************** !
!
! RunToTime: Runs the actual simulation.
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine RunToTime(this,sync_time)

  implicit none
  
  type(process_model_coupler_type), pointer :: this
  PetscReal :: target_time
  
  do
    call StepperSetTargetTime(this%timestepper,sync_time,option)

  
    
  enddo
  
end subroutine RunToTime

! ************************************************************************** !
!
! FinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine FinalizeRun(this)

  implicit none
  
  type(process_model_coupler_type), pointer :: this
  
end subroutine FinalizeRun


! ************************************************************************** !
!
! Output: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine Output(this,plot_flag,transient_plot_flag)

  implicit none
  
  type(process_model_coupler_type), pointer :: this
  
  call Output(this%realization,plot_flag,transient_plot_flag)
  
end subroutine Output

! ************************************************************************** !
!
! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine Destroy(this)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  type(process_model_coupler_type), pointer :: this
  
  if (.not.associated(process_model_coupler)) return

  if (associated(this%below)) then
    call this%below%Destroy()
  endif  
  if (associated(this%next) then
    call this%next%Destroy()
  endif 
  
  if (associated(cur_process_model)) then
    call this%process_model_list%Destroy()
  endif

  deallocate(process_model_coupler)
  nullify(process_model_coupler)
  
end subroutine ProcessModelCouplerDestroy
  
end module Process_Model_Coupler_module
