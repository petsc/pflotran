#ifdef SURFACE_FLOW
module Surface_Simulation_class

  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Surface_class
  use PMC_Base_class
  use Surface_Realization_class

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  private

  type, public, extends(simulation_base_type) :: surface_simulation_type
    class(pmc_surface_type), pointer :: surf_flow_process_model_coupler
    class(surface_realization_type), pointer :: surf_realization
    type(regression_type), pointer :: regression
  contains
    procedure, public :: Init => SurfaceSimulationInit
    procedure, public :: InitializeRun => SurfaceInitializeRun
    procedure, public :: FinalizeRun => SurfaceFinalizeRun
    procedure, public :: Strip => SurfaceSimulationStrip
  end type surface_simulation_type

  public :: SurfaceSimulationCreate, &
            SurfaceSimulationInit, &
            SurfaceFinalizeRun, &
            SurfaceSimulationStrip, &
            SurfaceSimulationDestroy

contains

! ************************************************************************** !

function SurfaceSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(surface_simulation_type), pointer :: SurfaceSimulationCreate
  
  print *, 'SurfaceSimulationCreate'
  
  allocate(SurfaceSimulationCreate)
  call SurfaceSimulationCreate%Init(option)
  
end function SurfaceSimulationCreate

! ************************************************************************** !

subroutine SurfaceSimulationInit(this,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  
  implicit none
  
  class(surface_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%regression)
  
end subroutine SurfaceSimulationInit

! ************************************************************************** !

subroutine SurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Logging_module
  use Output_module

  implicit none
  
  class(surface_simulation_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'SurfaceInitializeRun: Simulation%InitializeRun()')

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

end subroutine SurfaceInitializeRun

! ************************************************************************** !

subroutine SurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Timestepper_Base_class

  implicit none
  
  class(surface_simulation_type) :: this
  
  PetscErrorCode :: ierr

  class(stepper_base_type), pointer :: surf_flow_stepper

  call printMsg(this%option,'SurfaceFinalizeRun()')
  
  call SimulationBaseFinalizeRun(this)
  
  nullify(surf_flow_stepper)
  surf_flow_stepper => this%surf_flow_process_model_coupler%timestepper

end subroutine SurfaceFinalizeRun

! ************************************************************************** !

subroutine SurfaceSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(surface_simulation_type) :: this
  
  call printMsg(this%option,'SurfaceSimulationStrip()')
  
  call SimulationBaseStrip(this)
  call RegressionDestroy(this%regression)
  
end subroutine SurfaceSimulationStrip

! ************************************************************************** !

subroutine SurfaceSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(surface_simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'SurfaceSimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SurfaceSimulationDestroy

end module Surface_Simulation_class
#endif