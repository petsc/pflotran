module Process_Model_Coupler_module

  use Process_Model_Base_class
  use Timestepper_module
  use Option_module
  use Waypoint_module
  use Process_Model_module

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
    type(process_model_pointer_type), pointer :: pm_ptr
    PetscInt :: depth
  contains
    procedure, public :: InitializeRun
    procedure, public :: SetTimestepper => ProcModelCouplerSetTimestepper
    procedure, public :: RunTo => RunToTime
    procedure, public :: FinalizeRun
    procedure, public :: OutputLocal
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
  process_model_coupler%depth = 0
  
  allocate(process_model_coupler%pm_ptr)
  nullify(process_model_coupler%pm_ptr%ptr)
  
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
recursive subroutine RunToTime(this,sync_time,stop_flag)

  use Timestepper_module
  use Output_module, only : Output
  
  implicit none
  
  class(process_model_coupler_type) :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  class(process_model_base_type), pointer :: cur_process_model
  
  write(this%option%io_buffer,'(es12.5)') sync_time
  this%option%io_buffer = 'ProcessModelCoupler%RunToTime(' // &
    trim(adjustl(this%option%io_buffer)) // ')'
  call printMsg(this%option)
  
  local_stop_flag = 0
  do
    if (local_stop_flag > 0) exit ! end simulation
    if (this%timestepper%target_time >= sync_time) exit
    
    call SetOutputFlags(this)
    plot_flag = PETSC_FALSE
    transient_plot_flag = PETSC_FALSE
    
    call StepperSetTargetTime(this%timestepper,sync_time,this%option, &
                              local_stop_flag,plot_flag,transient_plot_flag)
    call StepperStepDT(this%timestepper,this%process_model_list, &
                       local_stop_flag)

    if (local_stop_flag > 1) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_process_model => this%process_model_list
    do
      if (.not.associated(cur_process_model)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_process_model%UpdateSolution()
      call StepperUpdateDT(this%timestepper,cur_process_model)
      cur_process_model => cur_process_model%next
    enddo
    ! Run underlying process model couplers
    if (associated(this%below)) then
      call this%below%RunTo(this%timestepper%target_time,local_stop_flag)
    endif
    
    ! only print output for process models of depth 0
    if (this%depth == 0) then
      if (this%timestepper%time_step_cut_flag) then
        plot_flag = PETSC_FALSE
      endif
      ! however, if we are using the modulus of the output_option%imod, we may
      ! still print
      if (mod(this%timestepper%steps, &
              this%process_model_list% &
                output_option%periodic_output_ts_imod) == 0) then
        plot_flag = PETSC_TRUE
      endif
      if (plot_flag .or. mod(this%timestepper%steps, &
                             this%process_model_list%output_option% &
                               periodic_tr_output_ts_imod) == 0) then
        transient_plot_flag = PETSC_TRUE
      endif
    endif
    call Output(this%process_model_list%realization_base, &
                plot_flag,transient_plot_flag)
    
  enddo
  
  stop_flag = max(stop_flag,local_stop_flag)
  
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
! SetOutputFlags: Toggles flags that determine whether output is printed
!                 to the screen and output file during a time step.
! author: Glenn Hammond
! date: 03/29/13
!
! ************************************************************************** !
subroutine SetOutputFlags(this)

  use Option_module
  use Output_Aux_module
  
  implicit none
  
  class(process_model_coupler_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  output_option => this%process_model_list%output_option

  if (OptionPrintToScreen(this%option) .and. &
      mod(this%timestepper%steps,output_option%screen_imod) == 0) then
    this%option%print_screen_flag = PETSC_TRUE
  else
    this%option%print_screen_flag = PETSC_FALSE
  endif

  if (OptionPrintToFile(this%option) .and. &
      mod(this%timestepper%steps,output_option%output_file_imod) == 0) then
    this%option%print_file_flag = PETSC_TRUE
  else
    this%option%print_file_flag = PETSC_FALSE
      
  endif
  
end subroutine SetOutputFlags

! ************************************************************************** !
!
! OutputLocal: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine OutputLocal(this)

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
    
end subroutine OutputLocal

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
