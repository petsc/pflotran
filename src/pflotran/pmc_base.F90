! Process Model Coupler Base class
module PMC_Base_class

  use Process_Model_Base_class
  use Timestepper_module
  use Option_module
  use Waypoint_module
  use Process_Model_module

  implicit none

#include "definitions.h"
  
  private
  
  ! process model coupler type
  type, public :: pmc_base_type
    type(option_type), pointer :: option
    type(stepper_type), pointer :: timestepper
    class(pm_base_type), pointer :: pm_list
    type(waypoint_list_type), pointer :: waypoints
    class(pmc_base_type), pointer :: below
    class(pmc_base_type), pointer :: next
    type(pm_pointer_type), pointer :: pm_ptr
    PetscInt :: depth
  contains
    procedure, public :: Init => PMCInit
    procedure, public :: InitializeRun
    procedure, public :: CastToBase => PMCCastToBase
    procedure, public :: SetTimestepper => ProcModelCouplerSetTimestepper
    procedure, public :: RunTo => RunToTime
    procedure, public :: FinalizeRun
    procedure, public :: OutputLocal
    procedure, public :: UpdateSolution => PMCUpdateSolution
    procedure, public :: Destroy
  end type pmc_base_type
  
  public :: PMCCreate, &
            PMCInit
  
contains

! ************************************************************************** !
!
! PMCCreate: Allocates and initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
function PMCCreate()

  implicit none
  
  class(pmc_base_type), pointer :: PMCCreate
  
  class(pmc_base_type), pointer :: pmc

  print *, 'PMC%Create()'
  
  allocate(pmc)
  call pmc%Init()

  PMCCreate => pmc  
  
end function PMCCreate

! ************************************************************************** !
!
! PMCInit: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PMCInit(this)

  implicit none
  
  class(pmc_base_type) :: this
  
  print *, 'PMC%Init()'
  
  nullify(this%option)
  nullify(this%timestepper)
  nullify(this%pm_list)
  nullify(this%waypoints)
  nullify(this%below)
  nullify(this%next)
  this%depth = 0
  
  allocate(this%pm_ptr)
  nullify(this%pm_ptr%ptr)
  
end subroutine PMCInit

! ************************************************************************** !
!
! PMCCastToBase: Initializes a new process model coupler object.
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
function PMCCastToBase(this)

  implicit none
  
  class(pmc_base_type), target :: this
  
  class(pmc_base_type), pointer :: PMCCastToBase

  PMCCastToBase => this
  
end function PMCCastToBase

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
  
  class(pmc_base_type) :: this
  type(stepper_type), pointer :: timestepper

  call printMsg(this%option,'PMC%SetTimestepper()')
  
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
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
  call printMsg(this%option,'PMC%InitializeRun()')
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%InitializeRun()
    cur_pm => cur_pm%next
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
  
  class(pmc_base_type) :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: plot_flag
  PetscBool :: transient_plot_flag
  class(pm_base_type), pointer :: cur_pm
  
  write(this%option%io_buffer,'(es12.5)') sync_time
  this%option%io_buffer = 'PMC%RunToTime(' // &
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
    call StepperStepDT(this%timestepper,this%pm_list, &
                       local_stop_flag)

    if (local_stop_flag > 1) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_pm => this%pm_list
    do
      if (.not.associated(cur_pm)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_pm%UpdateSolution()
      call StepperUpdateDT(this%timestepper,cur_pm)
      cur_pm => cur_pm%next
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
              this%pm_list% &
                output_option%periodic_output_ts_imod) == 0) then
        plot_flag = PETSC_TRUE
      endif
      if (plot_flag .or. mod(this%timestepper%steps, &
                             this%pm_list%output_option% &
                               periodic_tr_output_ts_imod) == 0) then
        transient_plot_flag = PETSC_TRUE
      endif
      call Output(this%pm_list%realization_base, &
                  plot_flag,transient_plot_flag)
    endif
    
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
  
  class(pmc_base_type) :: this

  class(pm_base_type), pointer :: cur_pm
  
  call printMsg(this%option,'PMC%UpdateSolution()')
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    cur_pm => cur_pm%next
  enddo  
  
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
  
  class(pmc_base_type) :: this
  
  character(len=MAXSTRINGLENGTH) :: string
  
  call printMsg(this%option,'PMC%FinalizeRun()')
  
  if (OptionPrintToScreen(this%option)) then
    write(*,'(/," PMC steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            this%timestepper%steps, &
            this%timestepper%cumulative_newton_iterations, &
            this%timestepper%cumulative_linear_iterations, &
            this%timestepper%cumulative_time_step_cuts
    write(string,'(f12.1)') this%timestepper%cumulative_solver_time
    write(*,*) 'PMC SNES time = ' // trim(adjustl(string)) // ' seconds'
  endif

  if (associated(this%below)) then
    call this%below%FinalizeRun()
  endif
  
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
  
  class(pmc_base_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  output_option => this%pm_list%output_option

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
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
  call printMsg(this%option,'PMC%Output()')
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
!    call Output(cur_pm%realization,plot_flag,transient_plot_flag)
    cur_pm => cur_pm%next
  enddo
    
end subroutine OutputLocal

! ************************************************************************** !
!
! PMCDestroy: Deallocates a pmc object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
recursive subroutine Destroy(this)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_base_type) :: this
  
  call printMsg(this%option,'PMC%Destroy()')
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif 
  
  if (associated(this%pm_list)) then
    call this%pm_list%Destroy()
  endif

!  deallocate(pmc)
!  nullify(pmc)
  
end subroutine Destroy
  
end module PMC_Base_class
