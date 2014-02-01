#ifdef SURFACE_FLOW
module Surf_Subsurf_Simulation_class

  use Surface_Simulation_class
  use Subsurface_Simulation_class
  use Regression_module
  use Option_module
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Surface_class
  use Realization_class
  use Surface_Realization_class

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public, extends(subsurface_simulation_type) :: surfsubsurface_simulation_type
    class(pmc_surface_type), pointer    :: surf_flow_process_model_coupler
    class(surface_realization_type), pointer :: surf_realization
  contains
    procedure, public :: Init => SurfSubsurfaceSimulationInit
    procedure, public :: InitializeRun => SurfSubsurfaceInitializeRun
    procedure, public :: FinalizeRun => SurfSubsurfaceFinalizeRun
    procedure, public :: Strip => SurfSubsurfaceSimulationStrip
    procedure, public :: ExecuteRun => SurfSubsurfaceExecuteRun
    procedure, public :: RunToTime => SurfSubsurfaceSimulationRunToTime
  end type surfsubsurface_simulation_type

  public :: SurfSubsurfaceSimulationCreate, &
            SurfSubsurfaceSimulationInit, &
            SurfSubsurfaceFinalizeRun, &
            SurfSubsurfaceSimulationStrip, &
            SurfSubsurfaceSimulationDestroy

contains

! ************************************************************************** !

function SurfSubsurfaceSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(surfsubsurface_simulation_type), pointer :: SurfSubsurfaceSimulationCreate
  
  print *, 'SurfSubsurfaceSimulationCreate'
  
  allocate(SurfSubsurfaceSimulationCreate)
  call SurfSubsurfaceSimulationCreate%Init(option)
  
end function SurfSubsurfaceSimulationCreate

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationInit(this,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Option_module
  
  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SubsurfaceSimulationInit(this,option)
  nullify(this%surf_realization)
  
end subroutine SurfSubsurfaceSimulationInit

! ************************************************************************** !

subroutine SurfSubsurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Logging_module
  use Output_module
  use PMC_Surface_class

  implicit none
  
#include "finclude/petscviewer.h"

  class(surfsubsurface_simulation_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  call printMsg(this%option,'Simulation%InitializeRun()')

  call this%process_model_coupler_list%InitializeRun()

  if (this%option%restart_flag) then
    call this%process_model_coupler_list%Restart(viewer)
    cur_process_model_coupler => this%process_model_coupler_list
    select type(pmc => cur_process_model_coupler)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (RICHARDS_MODE)
            call pmc%PMCSurfaceGetAuxDataAfterRestart()
          case (TH_MODE)
            call pmc%PMCSurfaceGetAuxDataAfterRestart()
          case default
            call printErrMsg(this%option,'SurfSubsurfaceInitializeRun ' // &
                  'not supported in current flow mode.')
        end select
    end select

  endif

end subroutine SurfSubsurfaceInitializeRun

! ************************************************************************** !

subroutine SurfSubsurfaceExecuteRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Simulation_Base_class

  implicit none
  
#include "finclude/petscviewer.h"

  class(surfsubsurface_simulation_type) :: this

  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt
  PetscViewer :: viewer

  time = 0.d0
  time = this%option%time

  final_time = SimulationGetFinalWaypointTime(this)

  call printMsg(this%option,'SurfSubsurfaceExecuteRun()')

  if(.not.associated(this%surf_realization)) then
    call this%RunToTime(final_time)

  else

    ! If simulation is decoupled surface-subsurface simulation, set
    ! dt_coupling to be dt_max
    if (this%surf_realization%dt_coupling == 0.d0) &
      this%surf_realization%dt_coupling = this%surf_realization%dt_max

    do
      if (time + this%surf_realization%dt_coupling > final_time) then
        dt = final_time-time
      else
        dt = this%surf_realization%dt_coupling
      endif

      time = time + dt
      call this%RunToTime(time)
      
      if (time >= final_time) exit
    enddo

  endif
  if (this%option%checkpoint_flag) then
    call this%process_model_coupler_list%Checkpoint(viewer,-1)
  endif

end subroutine SurfSubsurfaceExecuteRun

! ************************************************************************** !

subroutine SurfSubsurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Simulation_Base_class
  use Timestepper_Base_class

  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'SurfSubsurfaceFinalizeRun()')
  
  call SubsurfaceFinalizeRun(this)
  !call SurfaceFinalizeRun(this)
  
end subroutine SurfSubsurfaceFinalizeRun

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Simulation_Base_class

  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  
  call printMsg(this%option,'SurfSubsurfaceSimulationStrip()')
  
  call SubsurfaceSimulationStrip(this)
  !call SubsurfaceSimulationStrip(this)
  call RegressionDestroy(this%regression)
  
end subroutine SurfSubsurfaceSimulationStrip

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationRunToTime(this,target_time)
  ! 
  ! This routine executes surface-subsurface simualation
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  use Simulation_Aux_module

  implicit none

#include "finclude/petscviewer.h"

  class(surfsubsurface_simulation_type) :: this
  PetscReal :: target_time

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscViewer :: viewer

  call printMsg(this%option,'RunToTime()')
  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine SurfSubsurfaceSimulationRunToTime

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  implicit none
  
  class(surfsubsurface_simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SurfSubsurfaceSimulationDestroy

end module Surf_Subsurf_Simulation_class
#endif
