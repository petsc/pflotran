module Subsurface_Simulation_class
  
  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Subsurface_class
  use PMC_Base_class
  use Realization_class

  implicit none

#include "definitions.h"
  
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
    procedure, public :: InitializeRun => SubsurfaceInitializeRun
!    procedure, public :: ExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SubsurfaceFinalizeRun
  end type subsurface_simulation_type
  
  public :: SubsurfaceSimulationCreate, &
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

  type(pmc_base_type), pointer :: cur_process_model_coupler
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'Simulation%InitializeRun()')
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%InitializeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

  !TODO(geh): place logic here to stop if only initial state desired (e.g.
  !           solution composition, etc.).
  
  !TODO(geh): replace integer arguments with logical
#ifndef SIMPLIFY  
  if (this%option%restart_flag) then
    call OutputInit(1) ! number greater than 0
  else
    call OutputInit(0)
  endif
#endif
  
  if (this%output_option%print_initial) then
    cur_process_model_coupler => this%process_model_coupler_list
    do
      if (.not.associated(cur_process_model_coupler)) exit
      ! arg1 = plot_flag
      ! arg2 = transient_plot_flag
!      call cur_process_model_coupler%Output(PETSC_TRUE,PETSC_TRUE)
      cur_process_model_coupler => cur_process_model_coupler%next
    enddo
  endif
  
  !TODO(geh): place logic here to stop if only initial condition desired
  
  ! pushed in Init()
  call PetscLogStagePop(ierr)
  this%option%init_stage = PETSC_FALSE

  ! popped in FinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr)
  
end subroutine SubsurfaceInitializeRun

! ************************************************************************** !
!
! SubsurfaceFinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine SubsurfaceFinalizeRun(this)

  use Timestepper_module

  implicit none
  
  class(subsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  type(pmc_base_type), pointer :: cur_process_model_coupler
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper

  call printMsg(this%option,'Simulation%FinalizeRun()')
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%FinalizeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo
  
  nullify(flow_stepper)
  nullify(tran_stepper)
  
  if (associated(this%flow_process_model_coupler)) &
    flow_stepper => this%flow_process_model_coupler%timestepper
  if (associated(this%rt_process_model_coupler)) &
    tran_stepper => this%rt_process_model_coupler%timestepper
  
  call RegressionOutput(this%regression,this%realization, &
                        flow_stepper,tran_stepper)  
  
  ! pushed in InitializeRun()
  call PetscLogStagePop(ierr)
  
end subroutine SubsurfaceFinalizeRun

! ************************************************************************** !
!
! SubsurfaceSimulationDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SubsurfaceSimulationDestroy(simulation)

#ifndef SIMPLIFY
  use Richards_module, only : RichardsDestroy
  use Reactive_Transport_module, only : RTDestroy
  use General_module, only : GeneralDestroy
#endif

  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  class(simulation_base_type), pointer :: simulation_base
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call RegressionDestroy(simulation%regression)
  
  simulation_base => simulation
  call SimulationBaseDestroy(simulation_base)
  
end subroutine SubsurfaceSimulationDestroy
  
end module Subsurface_Simulation_class
