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

  implicit none

  private

#include "definitions.h"

  type, public, extends(subsurface_simulation_type) :: surfsubsurface_simulation_type
    class(pmc_surface_type), pointer    :: surf_flow_process_model_coupler
    class(surface_realization_type), pointer :: surf_realization
  contains
    procedure, public :: Init => SurfSubsurfaceSimulationInit
    procedure, public :: InitializeRun => SurfSubsurfaceInitializeRun
    procedure, public :: FinalizeRun => SurfSubsurfaceFinalizeRun
    procedure, public :: Strip => SurfSubsurfaceSimulationStrip
    procedure, public :: ExecuteRun => SurfSubsurfaceExecuteRun
  end type surfsubsurface_simulation_type

  public :: SurfSubsurfaceSimulationCreate, &
            SurfSubsurfaceSimulationInit, &
            SurfSubsurfaceFinalizeRun, &
            SurfSubsurfaceSimulationStrip, &
            SurfSubsurfaceSimulationDestroy

contains

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
function SurfSubsurfaceSimulationCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(surfsubsurface_simulation_type), pointer :: SurfSubsurfaceSimulationCreate
  
  print *, 'SurfSubsurfaceSimulationCreate'
  
  allocate(SurfSubsurfaceSimulationCreate)
  call SurfSubsurfaceSimulationCreate%Init(option)
  
end function SurfSubsurfaceSimulationCreate

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceSimulationInit(this,option)

  use Option_module
  
  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SubsurfaceSimulationInit(this,option)
  nullify(this%surf_realization)
  
end subroutine SurfSubsurfaceSimulationInit

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceInitializeRun(this)

  use Logging_module
  use Output_module

  implicit none
  
  class(surfsubsurface_simulation_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Simulation%InitializeRun()')

  call this%process_model_coupler_list%InitializeRun()

end subroutine SurfSubsurfaceInitializeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceExecuteRun(this)

  use Simulation_Base_class

  implicit none
  
  class(surfsubsurface_simulation_type) :: this

  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt

  time = 0.d0
  final_time = SimulationGetFinalWaypointTime(this)

  call printMsg(this%option,'SurfSubsurfaceExecuteRun()')

  if(.not.associated(this%surf_realization)) then
    call this%RunToTime(final_time)

  else

    do
      if(time+this%surf_realization%dt_coupling > final_time) then
        dt = final_time-time
      else
        dt = this%surf_realization%dt_coupling
      endif

      time = time + dt
      call this%RunToTime(time)
      
      if (time >= final_time) exit
    enddo

  endif

end subroutine SurfSubsurfaceExecuteRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceFinalizeRun(this)

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
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceSimulationStrip(this)

  use Simulation_Base_class

  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  
  call printMsg(this%option,'SurfSubsurfaceSimulationStrip()')
  
  call SubsurfaceSimulationStrip(this)
  !call SubsurfaceSimulationStrip(this)
  call RegressionDestroy(this%regression)
  
end subroutine SurfSubsurfaceSimulationStrip

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceSimulationDestroy(simulation)

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
