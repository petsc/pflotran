module Hydrogeophysics_Simulation_class
  
  use Option_module
  use Subsurface_Simulation_class
!  use PMC_Hydrogeophysics_class

  implicit none

#include "definitions.h"
  
  private

  type, public, extends(subsurface_simulation_type) :: &
    hydrogeophysics_simulation_type
  contains
    procedure, public :: Init => HydrogeophysicsInit
    procedure, public :: InitializeRun => HydrogeophysicsInitializeRun
    procedure, public :: ExecuteRun => HydrogeophysicsExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => HydrogeophysicsFinalizeRun
    procedure, public :: Strip => HydrogeophysicsStrip
  end type hydrogeophysics_simulation_type
  
  public :: HydrogeophysicsCreate, &
            HydrogeophysicsDestroy
  
contains

! ************************************************************************** !
!
! HydrogeophysicsCreate: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
function HydrogeophysicsCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(hydrogeophysics_simulation_type), pointer :: HydrogeophysicsCreate
  
  print *, 'SimulationCreate'
  
  allocate(HydrogeophysicsCreate)
  call HydrogeophysicsCreate%Init(option)
  
end function HydrogeophysicsCreate


! ************************************************************************** !
!
! HydrogeophysicsInit: Initializes simulation values
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInit(this,option)

  use Option_module
  
  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SubsurfaceSimulationInit(this,option)
  
end subroutine HydrogeophysicsInit

! ************************************************************************** !
!
! HydrogeophysicsInitializeRun: Initializes simulation
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsInitializeRun(this)

  use Output_module
  use PMC_Base_class

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Simulation%InitializeRun()')

  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    depth = 0
    call cur_process_model_coupler%InitializeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

  ! set depth in tree
  cur_process_model_coupler_top => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    depth = 0
    cur_process_model_coupler_top%depth = depth
    cur_process_model_coupler_below => cur_process_model_coupler_top%below
    do
      if (.not.associated(cur_process_model_coupler_below)) exit
      depth = depth + 1
      cur_process_model_coupler_below%depth = depth
      cur_process_model_coupler_below => cur_process_model_coupler_below%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo

end subroutine HydrogeophysicsInitializeRun

! ************************************************************************** !
!
! HydrogeophysicsExecuteRun: 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HydrogeophysicsExecuteRun(this)

  use Simulation_Base_class

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscReal :: final_time
  PetscReal :: dt
  PetscReal :: current_time
  
  call printMsg(this%option,'HydrogeophysicsExecuteRun()')

  final_time = SimulationGetFinalWaypointTime(this)
  ! take hourly steps until final time
  current_time = 0.d0
  dt = 365.d0*24.d0*3600.d0
  do
    current_time = min(current_time + dt,final_time)
    call this%RunToTime(current_time)
    if (this%stop_flag > 0) exit
  enddo
  
end subroutine HydrogeophysicsExecuteRun
  
! ************************************************************************** !
!
! HydrogeophysicsFinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsFinalizeRun(this)

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'HydrogeophysicsFinalizeRun()')
  
  call SubsurfaceFinalizeRun(this)
  
end subroutine HydrogeophysicsFinalizeRun

! ************************************************************************** !
!
! HydrogeophysicsStrip: Deallocates members of hydrogeophysics simulation 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HydrogeophysicsStrip(this)

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  call printMsg(this%option,'HydrogeophysicsStrip()')
  
  call SubsurfaceSimulationStrip(this)
  
end subroutine HydrogeophysicsStrip

! ************************************************************************** !
!
! HydrogeophysicsDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
subroutine HydrogeophysicsDestroy(simulation)

  implicit none
  
  class(hydrogeophysics_simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine HydrogeophysicsDestroy
  
end module Hydrogeophysics_Simulation_class
