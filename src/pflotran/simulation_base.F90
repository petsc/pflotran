module Simulation_Base_class

  use PMC_Base_class
  use Option_module
  use Output_Aux_module
  use Output_module
  use Simulation_Aux_module
  use Waypoint_module
  
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public :: simulation_base_type
    type(option_type), pointer :: option
    type(waypoint_list_type), pointer :: waypoint_list ! for outer sync loop
    type(output_option_type), pointer :: output_option
    PetscInt :: stop_flag
    class(pmc_base_type), pointer :: process_model_coupler_list
    type(simulation_aux_type), pointer :: sim_aux
  contains
    procedure, public :: Init => SimulationBaseInit
    procedure, public :: InitializeRun => SimulationBaseInitializeRun
    procedure, public :: JumpStart => SimulationBaseJumpStart
    procedure, public :: ExecuteRun
    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SimulationBaseFinalizeRun
    procedure, public :: Strip => SimulationBaseStrip
  end type simulation_base_type
  
  public :: SimulationBaseCreate, &
            SimulationBaseInit, &
            SimulationBaseInitializeRun, &
            SimulationGetFinalWaypointTime, &
            SimulationBaseFinalizeRun, &
            SimulationBaseStrip, &
            SimulationBaseDestroy
  
contains

! ************************************************************************** !

function SimulationBaseCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Option_module

  implicit none
  
  class(simulation_base_type), pointer :: SimulationBaseCreate

  type(option_type), pointer :: option
  
  allocate(SimulationBaseCreate)
  call SimulationBaseCreate%Init(option)

end function SimulationBaseCreate

! ************************************************************************** !

subroutine SimulationBaseInit(this,option)
  ! 
  ! Initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 
  use Timestepper_Base_class, only : TS_CONTINUE
  use Option_module

  implicit none
  
  class(simulation_base_type) :: this
  type(option_type), pointer :: option

  this%option => option
  nullify(this%waypoint_list)
  nullify(this%output_option)
  nullify(this%process_model_coupler_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = TS_CONTINUE

end subroutine SimulationBaseInit

! ************************************************************************** !

subroutine SimulationBaseInitializeRun(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Logging_module

  implicit none
  
#include "finclude/petscviewer.h"  

  class(simulation_base_type) :: this

  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseInitializeRun()')
#endif
  
  if (this%option%restart_flag) then
    call this%process_model_coupler_list%Restart(viewer)
  endif
  
  ! initialize performs overwrite of restart, if applicable
  call this%process_model_coupler_list%InitializeRun()  
  call this%JumpStart()
  
  ! pushed in Init()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in FinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)
  
end subroutine SimulationBaseInitializeRun

! ************************************************************************** !

subroutine SimulationBaseJumpStart(this)
  ! 
  ! Gets the time stepping, etc. up and running
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  ! 
  use Option_module
  
  implicit none
  
  class(simulation_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseJumpStart()')
#endif

  this%option%io_buffer = 'SimulationBaseJumpStart must be extended for ' // &
    'each simulation mode.'
  call printErrMsg(this%option)
  
end subroutine SimulationBaseJumpStart

! ************************************************************************** !

subroutine ExecuteRun(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 
  use Waypoint_module

  implicit none
  
  class(simulation_base_type) :: this
  
  PetscReal :: final_time
  PetscReal :: sync_time
  type(waypoint_type), pointer :: cur_waypoint
  PetscViewer :: viewer
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseExecuteRun()')
#endif

  final_time = SimulationGetFinalWaypointTime(this)
  cur_waypoint => this%waypoint_list%first
  call WaypointSkipToTime(cur_waypoint,this%option%time)
  do
    if (.not.associated(cur_waypoint)) exit
    call this%RunToTime(min(final_time,cur_waypoint%time))
    cur_waypoint => cur_waypoint%next
  enddo
  if (this%option%checkpoint_flag) then
    call this%process_model_coupler_list%Checkpoint(viewer,-1)
  endif
  
end subroutine ExecuteRun

! ************************************************************************** !

subroutine RunToTime(this,target_time)
  ! 
  ! Executes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Option_module
  use Simulation_Aux_module

  implicit none
  
#include "finclude/petscviewer.h" 

  class(simulation_base_type) :: this
  PetscReal :: target_time
  
  class(pmc_base_type), pointer :: cur_process_model_coupler
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseRunToTime()')
#endif
  
  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine RunToTime

! ************************************************************************** !

subroutine SimulationBaseFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Logging_module
  use Timestepper_Base_class, only : TS_STOP_WALLCLOCK_EXCEEDED
  
  implicit none
  
  class(simulation_base_type) :: this
  
  PetscErrorCode :: ierr
  
  class(pmc_base_type), pointer :: cur_process_model_coupler

#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseFinalizeRun()')
#endif
  
  if (this%stop_flag == TS_STOP_WALLCLOCK_EXCEEDED) then
    call printMsg(this%option,"Wallclock stop time exceeded.  Exiting!!!")
    call printMsg(this%option,"")
  endif
  
  call this%process_model_coupler_list%FinalizeRun()
  
  ! pushed in InitializeRun()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)
  ! popped in OptionFinalize()
  call PetscLogStagePush(logging%stage(FINAL_STAGE),ierr);CHKERRQ(ierr)
  
end subroutine SimulationBaseFinalizeRun

! ************************************************************************** !

function SimulationGetFinalWaypointTime(this)
  ! 
  ! Returns the earliest final waypoint time
  ! from the top layer of process model
  ! couplers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/12/13
  ! 

  use Waypoint_module

  implicit none
  
  class(simulation_base_type) :: this
  
  PetscReal :: SimulationGetFinalWaypointTime

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscReal :: final_time
  
  SimulationGetFinalWaypointTime = 0.d0
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    final_time = WaypointListGetFinalTime(cur_process_model_coupler%waypoint_list)
    if (SimulationGetFinalWaypointTime < 1.d-40 .or. &
        final_time < SimulationGetFinalWaypointTime) then
      SimulationGetFinalWaypointTime = final_time
    endif
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

end function SimulationGetFinalWaypointTime

! ************************************************************************** !

subroutine SimulationBaseStrip(this)
  ! 
  ! Deallocates members of simulation base
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 
  use Input_Aux_module
  use Waypoint_module
  
  implicit none
  
  class(simulation_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseStrip()')
#endif
  call WaypointListDestroy(this%waypoint_list)
  call SimAuxDestroy(this%sim_aux)
  if (associated(this%process_model_coupler_list)) then
    call this%process_model_coupler_list%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%process_model_coupler_list)
    nullify(this%process_model_coupler_list)
  endif
  call InputDbaseDestroy()
  
end subroutine SimulationBaseStrip

! ************************************************************************** !

subroutine SimulationBaseDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  implicit none
  
  class(simulation_base_type), pointer :: simulation
  
#ifdef DEBUG
  call printMsg(simulation%option,'SimulationDestroy()')
#endif
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationBaseDestroy
  
end module Simulation_Base_class
