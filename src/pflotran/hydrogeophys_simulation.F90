module Hydrogeophys_Simulation_module
  
  use Simulation_Base_module
  use Subsurface_Simulation_module
  use Option_module
  use Process_Model_Coupler_module

  implicit none

#include "definitions.h"
  
  private

  type, public, extends(simulation_base_type) :: subsurface_simulation_type
  contains
    procedure, public :: Init => SubsurfaceSimulationInit
    procedure, public :: InitializeRun => SubsurfaceInitializeRun
!    procedure, public :: ExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SubsurfaceFinalizeRun
  end type subsurface_simulation_type
  
  public :: SubsurfaceSimulationCreate, &
            SimulationDestroy
  
  interface SimulationDestroy
    module procedure SimulationBaseDestroy
    module procedure SubsurfaceSimulationDestroy
  end interface
  
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
  
  class(subsurface_simulation_type), pointer :: simulation

  print *, 'SimulationCreate'
  
  allocate(simulation)
  call simulation%Init(option)
  
  SubsurfaceSimulationCreate => simulation
  
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

  type(process_model_coupler_type), pointer :: cur_process_model_coupler
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

  use Realization_class
  use Timestepper_module
  use Process_Model_Base_class
  use Process_Model_Richards_class
  use Process_Model_TH_class
  use Process_Model_THC_class
  use Process_Model_RT_class

  implicit none
  
  class(subsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  type(process_model_coupler_type), pointer :: cur_process_model_coupler_top
  class(process_model_base_type), pointer :: cur_process_model
  class(realization_type), pointer :: realization
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper

  call printMsg(this%option,'Simulation%FinalizeRun()')
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%FinalizeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo
  
  ! temporary pointers for regression
  nullify(flow_stepper)
  nullify(tran_stepper)
  cur_process_model_coupler_top => this%process_model_coupler_list
  do ! loop over next process model coupler
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler => cur_process_model_coupler_top
    do ! loop over process model couplers below
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%process_model_list
      do
        if (.not.associated(cur_process_model)) exit
        select type(cur_process_model)
          class is (process_model_richards_type)
            realization => cur_process_model%realization
            flow_stepper => cur_process_model_coupler%timestepper
          class is (process_model_th_type)
            realization => cur_process_model%realization
            flow_stepper => cur_process_model_coupler%timestepper
          class is (process_model_thc_type)
            realization => cur_process_model%realization
            flow_stepper => cur_process_model_coupler%timestepper
          class is (process_model_rt_type)
            realization => cur_process_model%realization
            tran_stepper => cur_process_model_coupler%timestepper
        end select
        cur_process_model => cur_process_model%next
      enddo
      cur_process_model_coupler => cur_process_model_coupler%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo
  
  call RegressionOutput(this%regression,realization, &
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
  
end module Subsurface_Simulation_module
