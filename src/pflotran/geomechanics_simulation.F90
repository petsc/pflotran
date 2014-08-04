module Geomechanics_Simulation_class

  use Option_module
  use Subsurface_Simulation_class
  use Regression_module
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Geomechanics_class
  use Realization_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(subsurface_simulation_type) :: geomechanics_simulation_type
    ! pointer to geomechanics coupler
    class(pmc_geomechanics_type), pointer :: geomech_process_model_coupler
    class(geomech_realization_type), pointer :: geomech_realization
  contains
    procedure, public :: Init => GeomechanicsSimulationInit
    procedure, public :: InitializeRun => GeomechanicsSimulationInitializeRun
    procedure, public :: ExecuteRun => GeomechanicsSimulationExecuteRun
    procedure, public :: FinalizeRun => GeomechanicsSimulationFinalizeRun
    procedure, public :: Strip => GeomechanicsSimulationStrip
    !procedure, public :: RunToTime => GeomechanicsSimulationRunToTime
  end type Geomechanics_simulation_type
  
  public :: GeomechanicsSimulationCreate, &
            GeomechanicsSimulationDestroy
  
contains

! ************************************************************************** !

function GeomechanicsSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module

  implicit none

  type(option_type), pointer :: option

  class(geomechanics_simulation_type), pointer :: GeomechanicsSimulationCreate

  print *,'GeomechanicsSimulationCreate'

  allocate(GeomechanicsSimulationCreate)
  call GeomechanicsSimulationCreate%Init(option)

end function GeomechanicsSimulationCreate

! ************************************************************************** !

subroutine GeomechanicsSimulationInit(this, option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module

  implicit none

  class(geomechanics_simulation_type) :: this
  type(option_type), pointer :: option

  call SubsurfaceSimulationInit(this, option)
  nullify(this%geomech_realization)

end subroutine GeomechanicsSimulationInit

! ************************************************************************** !

subroutine GeomechanicsSimulationInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Output_module
  use PMC_Geomechanics_class

  implicit none

  class(geomechanics_simulation_type) :: this

  call printMsg(this%option,'Simulation%InitializeRun()')
  call this%process_model_coupler_list%InitializeRun()

  if (this%option%restart_flag) then
    call printErrMsg(this%option,'add code for restart of GeomechanicsSimulation')
  endif

end subroutine GeomechanicsSimulationInitializeRun

! ************************************************************************** !

subroutine GeomechanicsSimulationExecuteRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Simulation_Base_class

  implicit none

#include "finclude/petscviewer.h"

  class(geomechanics_simulation_type) :: this

  PetscReal :: final_time

  final_time = SimulationGetFinalWaypointTime(this)
  call this%RunToTime(final_time)

end subroutine GeomechanicsSimulationExecuteRun

! ************************************************************************** !

subroutine GeomechanicsSimulationFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Simulation_Base_class
  use Timestepper_Base_class

  implicit none

  class(geomechanics_simulation_type) :: this

  call printMsg(this%option,'GeomechanicsSimulationFinalizeRun')

  call SubsurfaceFinalizeRun(this)
  !call GeomechanicsFinalizeRun(this)

end subroutine GeomechanicsSimulationFinalizeRun

! ************************************************************************** !

subroutine GeomechanicsSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Simulation_Base_class

  implicit none
  
  class(geomechanics_simulation_type) :: this
  
  call printMsg(this%option,'GeomechanicsSimulationStrip()')
  
  call SubsurfaceSimulationStrip(this)
  call RegressionDestroy(this%regression)
  
end subroutine GeomechanicsSimulationStrip

! ************************************************************************** !

subroutine GeomechanicsSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(geomechanics_simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'GeomehanicsSimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine GeomechanicsSimulationDestroy

end module Geomechanics_Simulation_class
