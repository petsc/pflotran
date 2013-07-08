module Hydrogeophysics_Simulation_class
  
  use Option_module
  use Subsurface_Simulation_class
  use PMC_Hydrogeophysics_class

  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(subsurface_simulation_type) :: &
    hydrogeophysics_simulation_type
    ! pointer to hydrogeophysics coupler
    class(pmc_hydrogeophysics_type), pointer :: hydrogeophysics_coupler
    PetscMPIInt :: pf_e4d_comm
    PetscMPIInt :: pf_e4d_grp
    PetscMPIInt :: pf_e4d_size
    PetscMPIInt :: pf_e4d_rank
    PetscBool :: subsurface_process
    Vec :: sigma
  contains
    procedure, public :: Init => HydrogeophysicsInit
    procedure, public :: InitializeRun => HydrogeophysicsInitializeRun
!    procedure, public :: ExecuteRun => HydrogeophysicsExecuteRun
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
  
  call printMsg(option,'HydrogeophysicsCreate()')
  
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
  nullify(this%hydrogeophysics_coupler)
  this%sigma = 0
  ! -999 denotes uninitialized
  this%pf_e4d_comm = -999
  this%pf_e4d_grp = -999
  this%pf_e4d_size = -999
  this%pf_e4d_rank = -999
  this%subsurface_process = PETSC_FALSE
   
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
  
  call printMsg(this%option,'Hydrogeophysics%InitializeRun()')

  call this%process_model_coupler_list%InitializeRun()

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
  use e4d_setup, only : setup_e4d
  use e4d_run, only: run_e4D

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscReal :: final_time
  PetscReal :: dt
  PetscReal :: current_time
  
  call printMsg(this%option,'Hydrogeophysics%ExecuteRun()')

  if (this%subsurface_process) then
    final_time = SimulationGetFinalWaypointTime(this)
    ! take hourly steps until final time
    current_time = 0.d0
    dt = 365.d0*24.d0*3600.d0
    do
      current_time = min(current_time + dt,final_time)
      call this%RunToTime(current_time)
      if (this%stop_flag > 0) exit
    enddo
  else
    call setup_e4d
    call run_e4d    
  endif
  
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
  
  call printMsg(this%option,'Hydrogeophysics%FinalizeRun()')
  
  if (this%subsurface_process) then
    call SubsurfaceFinalizeRun(this)
  endif
  
end subroutine HydrogeophysicsFinalizeRun

! ************************************************************************** !
!
! HydrogeophysicsStrip: Deallocates members of hydrogeophysics simulation 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine HydrogeophysicsStrip(this)

  use Hydrogeophysics_Wrapper_module

  implicit none
  
  class(hydrogeophysics_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Hydrogeophysics%Strip()')
  
  if (this%subsurface_process) then
    call SubsurfaceSimulationStrip(this)
  else
    call HydrogeophysicsWrapperDestroy(this%option)
  endif
  if (this%sigma /= 0) &
    call VecDestroy(this%sigma ,ierr)
  
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
  
  call printMsg(simulation%option,'Hydrogeophysics%Destroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine HydrogeophysicsDestroy
  
end module Hydrogeophysics_Simulation_class
