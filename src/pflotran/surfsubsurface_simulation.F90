#ifdef SURFACE_FLOW
module SurfSubsurface_Simulation_class

  use Simulation_Base_class
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

  type, public, extends(simulation_base_type) :: surfsubsurface_simulation_type
    class(pmc_subsurface_type), pointer :: flow_process_model_coupler
    class(pmc_subsurface_type), pointer :: rt_process_model_coupler
    class(pmc_surface_type), pointer    :: surf_flow_process_model_coupler
    class(realization_type), pointer :: realization
    class(surface_realization_type), pointer :: surf_realization
    type(regression_type), pointer :: regression
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
  
  call SimulationBaseInit(this,option)
  nullify(this%regression)
  
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
    cur_process_model_coupler_below => cur_process_model_coupler_top%below
    do
      if (.not.associated(cur_process_model_coupler_below)) exit
      depth = depth + 1
      cur_process_model_coupler_below => cur_process_model_coupler_below%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo
 
end subroutine SurfSubsurfaceInitializeRun

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceFinalizeRun(this)

  use Timestepper_Base_class

  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  class(stepper_base_type), pointer :: flow_stepper
  class(stepper_base_type), pointer :: tran_stepper
  class(stepper_base_type), pointer :: surf_flow_stepper

  call printMsg(this%option,'SubsurfaceFinalizeRun()')
  
  call SimulationBaseFinalizeRun(this)
  
  nullify(flow_stepper)
  nullify(tran_stepper)
  nullify(surf_flow_stepper)

  if (associated(this%flow_process_model_coupler)) &
    flow_stepper => this%flow_process_model_coupler%timestepper
  if (associated(this%rt_process_model_coupler)) &
    tran_stepper => this%rt_process_model_coupler%timestepper
  surf_flow_stepper => this%surf_flow_process_model_coupler%timestepper

  call RegressionOutput(this%regression,this%realization, &
                        flow_stepper,tran_stepper)  
  
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

  implicit none
  
  class(surfsubsurface_simulation_type) :: this
  
  call printMsg(this%option,'SurfSubsurfaceSimulationStrip()')
  
  call SimulationBaseStrip(this)
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

! ************************************************************************** !
!> This routine
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 06/28/13
! ************************************************************************** !
subroutine SurfSubsurfaceExecuteRun(this)

  implicit none
  
  class(surfsubsurface_simulation_type) :: this

  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt

  time = 0.d0
  final_time = SimulationGetFinalWaypointTime(this)

  call printMsg(this%option,'SurfSubsurfaceExecuteRun()')

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

end subroutine SurfSubsurfaceExecuteRun


end module SurfSubsurface_Simulation_class
#endif
