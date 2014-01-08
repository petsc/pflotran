module Simulation_Base_class

  use PMC_Base_class
  use Option_module
  use Output_Aux_module
  use Output_module
  use Simulation_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public :: simulation_base_type
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    PetscInt :: stop_flag
    class(pmc_base_type), pointer :: process_model_coupler_list
    type(simulation_aux_type), pointer :: sim_aux
  contains
    procedure, public :: Init => SimulationBaseInit
    procedure, public :: InitializeRun
    procedure, public :: ExecuteRun
    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SimulationBaseFinalizeRun
    procedure, public :: Strip => SimulationBaseStrip
  end type simulation_base_type
  
  public :: SimulationBaseCreate, &
            SimulationBaseInit, &
            SimulationGetFinalWaypointTime, &
            SimulationBaseFinalizeRun, &
            SimulationBaseStrip, &
            SimulationBaseDestroy
  
contains

! ************************************************************************** !
!
! SimulationBaseCreate: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
function SimulationBaseCreate(option)

  use Option_module

  implicit none
  
  class(simulation_base_type), pointer :: SimulationBaseCreate

  type(option_type), pointer :: option
  
  allocate(SimulationBaseCreate)
  call SimulationBaseCreate%Init(option)

end function SimulationBaseCreate

! ************************************************************************** !
!
! SimulationBaseInit: Initializes a new simulation object
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SimulationBaseInit(this,option)

  use Option_module

  implicit none
  
  class(simulation_base_type) :: this
  type(option_type), pointer :: option

  this%option => option
  nullify(this%output_option)
  nullify(this%process_model_coupler_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = 0 

end subroutine SimulationBaseInit

! ************************************************************************** !
!
! InitializeRun: Initializes simulation
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine InitializeRun(this)

  use Logging_module

  implicit none
  
#include "finclude/petscviewer.h"  

  class(simulation_base_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseInitializeRun()')
#endif
  
  call this%process_model_coupler_list%InitializeRun()  

  !TODO(geh): place logic here to stop if only initial state desired (e.g.
  !           solution composition, etc.).
  
  ! update solutions?
  
  !TODO(geh): Place I/O reoutines here
  
  !TODO(geh): place logic here to stop if only initial condition desired
  
  if (this%option%restart_flag) then
    call this%process_model_coupler_list%Restart(viewer)
  endif
  
  ! pushed in Init()
  call PetscLogStagePop(ierr)

  ! popped in FinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)
  
end subroutine InitializeRun

! ************************************************************************** !
!
! ExecuteRun: Initializes simulation
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine ExecuteRun(this)

  implicit none
  
  class(simulation_base_type) :: this
  
  PetscReal :: final_time
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseExecuteRun()')
#endif

  final_time = SimulationGetFinalWaypointTime(this)
  call this%RunToTime(final_time)
  
#if 0  
  type(waypoint_type), pointer :: cur_waypoint
  
  cur_waypoint => this%waypoints%first
  do
    if (.not.associated(cur_waypoint)) exit
    call this%RunToTime(cur_waypoint%time)
    cur_waypoint => cur_waypoint%next

    if (this%option%wallclock_stop_flag) then
    !TODO(geh): reformulate 
#if 0    
      call PetscTime(current_time, ierr)
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
#endif  
  
end subroutine ExecuteRun

! ************************************************************************** !
!
! RunToTime: Executes simulation
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine RunToTime(this,target_time)

  use Option_module
  use Simulation_Aux_module

  implicit none
  
#include "finclude/petscviewer.h" 

  class(simulation_base_type) :: this
  PetscReal :: target_time
  
  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscViewer :: viewer
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseRunToTime()')
#endif
  
  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)
  if (this%option%checkpoint_flag) then
    call this%process_model_coupler_list%Checkpoint(viewer,-1)
  endif

end subroutine RunToTime

! ************************************************************************** !
!
! SimulationBaseFinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SimulationBaseFinalizeRun(this)

  use Logging_module
  
  implicit none
  
  class(simulation_base_type) :: this
  
  PetscErrorCode :: ierr
  
  class(pmc_base_type), pointer :: cur_process_model_coupler

#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseFinalizeRun()')
#endif
  
  call this%process_model_coupler_list%FinalizeRun()
  
  ! pushed in InitializeRun()
  call PetscLogStagePop(ierr)
  ! popped in OptionFinalize()
  call PetscLogStagePush(logging%stage(FINAL_STAGE),ierr)
  
end subroutine SimulationBaseFinalizeRun

! ************************************************************************** !
!
! SimulationGetFinalWaypointTime: Returns the earliest final waypoint time
!                                 from the top layer of process model 
!                                 couplers.
! author: Glenn Hammond
! date: 06/12/13
!
! ************************************************************************** !
function SimulationGetFinalWaypointTime(this)

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
    final_time = WaypointListGetFinalTime(cur_process_model_coupler%waypoints)
    if (SimulationGetFinalWaypointTime < 1.d-40 .or. &
        final_time < SimulationGetFinalWaypointTime) then
      SimulationGetFinalWaypointTime = final_time
    endif
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

end function SimulationGetFinalWaypointTime

! ************************************************************************** !
!
! SimulationBaseStrip: Deallocates members of simulation base
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SimulationBaseStrip(this)

  implicit none
  
  class(simulation_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseStrip()')
#endif
  call SimAuxDestroy(this%sim_aux)
  
end subroutine SimulationBaseStrip

! ************************************************************************** !
!
! SimulationBaseDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SimulationBaseDestroy(simulation)

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
