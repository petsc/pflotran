module Subsurface_Simulation_class
  
  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Subsurface_class
  use PMC_Base_class
  use Realization_class

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public, extends(simulation_base_type) :: subsurface_simulation_type
    ! pointer to flow process model coupler
    class(pmc_subsurface_type), pointer :: flow_process_model_coupler
    ! pointer to reactive transport process model coupler
    class(pmc_subsurface_type), pointer :: rt_process_model_coupler
    ! pointer to realization object shared by flow and reactive transport
    class(realization_type), pointer :: realization 
    ! regression object
    type(regression_type), pointer :: regression
  contains
    procedure, public :: Init => SubsurfaceSimulationInit
!    procedure, public :: InitializeRun => SubsurfaceInitializeRun
!    procedure, public :: ExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SubsurfaceFinalizeRun
    procedure, public :: Strip => SubsurfaceSimulationStrip
  end type subsurface_simulation_type
  
  public :: SubsurfaceSimulationCreate, &
            SubsurfaceSimulationInit, &
            SubsurfaceFinalizeRun, &
            SubsurfaceSimulationStrip, &
            SubsurfaceSimulationDestroy
  
contains

! ************************************************************************** !
!
! SubsurfaceSimulationCreate: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SubsurfaceSimulationCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(subsurface_simulation_type), pointer :: SubsurfaceSimulationCreate
  
  print *, 'SimulationCreate'
  
  allocate(SubsurfaceSimulationCreate)
  call SubsurfaceSimulationCreate%Init(option)
  
end function SubsurfaceSimulationCreate

! ************************************************************************** !
!
! SubsurfaceSimulationInit: Initializes simulation values
! author: Glenn Hammond
! date: 04/22/13
!
! ************************************************************************** !
subroutine SubsurfaceSimulationInit(this,option)

  use Option_module
  
  implicit none
  
  class(subsurface_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%regression)
  
end subroutine SubsurfaceSimulationInit

! ************************************************************************** !
!
! SubsurfaceInitializeRun: Initializes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine SubsurfaceInitializeRun(this)

  use Logging_module
  use Output_module

  implicit none
  
  class(subsurface_simulation_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Simulation%InitializeRun()')

  call this%process_model_coupler_list%InitializeRun()

end subroutine SubsurfaceInitializeRun

! ************************************************************************** !
!
! SubsurfaceFinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine SubsurfaceFinalizeRun(this)

  use Timestepper_BE_class

  implicit none
  
  class(subsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  class(stepper_BE_type), pointer :: flow_stepper
  class(stepper_BE_type), pointer :: tran_stepper

  call printMsg(this%option,'SubsurfaceFinalizeRun()')
  
  call SimulationBaseFinalizeRun(this)
  
  nullify(flow_stepper)
  nullify(tran_stepper)
  if (associated(this%flow_process_model_coupler)) then
    select type(ts => this%flow_process_model_coupler%timestepper)
      class is(stepper_BE_type)
        flow_stepper => ts
    end select
  endif
  if (associated(this%rt_process_model_coupler)) then
    select type(ts => this%rt_process_model_coupler%timestepper)
      class is(stepper_BE_type)
        tran_stepper => ts
    end select
  endif
  
  call RegressionOutput(this%regression,this%realization, &
                        flow_stepper,tran_stepper)  
  
end subroutine SubsurfaceFinalizeRun

! ************************************************************************** !
!
! SubsurfaceSimulationStrip: Deallocates members of subsurface simulation 
! author: Glenn Hammond
! date: 06/11/13
!
! ************************************************************************** !
subroutine SubsurfaceSimulationStrip(this)

  implicit none
  
  class(subsurface_simulation_type) :: this
  
  call printMsg(this%option,'SubsurfaceSimulationStrip()')
  
  call SimulationBaseStrip(this)
  call RegressionDestroy(this%regression)
  
end subroutine SubsurfaceSimulationStrip

! ************************************************************************** !
!
! SubsurfaceSimulationDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SubsurfaceSimulationDestroy(simulation)

  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SubsurfaceSimulationDestroy
  
end module Subsurface_Simulation_class
