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
    procedure, public :: JumpStart => SubsurfaceSimulationJumpStart
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

function SubsurfaceSimulationCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(subsurface_simulation_type), pointer :: SubsurfaceSimulationCreate
  
#ifdef DEBUG
  print *, 'SimulationCreate'
#endif
  
  allocate(SubsurfaceSimulationCreate)
  call SubsurfaceSimulationCreate%Init(option)
  
end function SubsurfaceSimulationCreate

! ************************************************************************** !

subroutine SubsurfaceSimulationInit(this,option)
  ! 
  ! Initializes simulation values
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/22/13
  ! 

  use Option_module
  
  implicit none
  
  class(subsurface_simulation_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%regression)
  
end subroutine SubsurfaceSimulationInit

! ************************************************************************** !

subroutine SubsurfaceSimulationJumpStart(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  ! 
  use Logging_module
  use Output_module
  use Option_module
  use Output_Aux_module
  use Timestepper_Base_class

  implicit none
  
  class(subsurface_simulation_type) :: this

  class(stepper_base_type), pointer :: master_stepper
  class(stepper_base_type), pointer :: flow_stepper
  class(stepper_base_type), pointer :: tran_stepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  PetscBool :: plot_flag, transient_plot_flag
  
#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceSimulationJumpStart()')
#endif

  nullify(master_stepper)
  nullify(flow_stepper)
  nullify(tran_stepper)

  option => this%option
  output_option => this%output_option

  ! first time stepper is master
  master_stepper => this%process_model_coupler_list%timestepper
  if (associated(this%flow_process_model_coupler)) then
    flow_stepper => this%flow_process_model_coupler%timestepper
  endif
  if (associated(this%rt_process_model_coupler)) then
    tran_stepper => this%rt_process_model_coupler%timestepper
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_stepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    option%status = DONE
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(master_stepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_stepper%max_time_step >= 0 .and. &
      output_option%print_initial) then
    plot_flag = PETSC_TRUE
    transient_plot_flag = PETSC_TRUE
    call Output(this%realization,plot_flag,transient_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_stepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_stepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    option%status = DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_stepper)) then
    if (.not.associated(flow_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      flow_stepper%dt_max = flow_stepper%cur_waypoint%dt_max
    endif
  endif  
  if (associated(tran_stepper)) then
    if (.not.associated(tran_stepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      tran_stepper%dt_max = tran_stepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_stepper)) &
    flow_stepper%start_time_step = flow_stepper%steps + 1
  if (associated(tran_stepper)) &
    tran_stepper%start_time_step = tran_stepper%steps + 1
  
  if (this%realization%debug%print_couplers) then
    call OutputPrintCouplers(this%realization,ZERO_INTEGER)
  endif  

end subroutine SubsurfaceSimulationJumpStart

! ************************************************************************** !

subroutine SubsurfaceFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_BE_class

  implicit none
  
  class(subsurface_simulation_type) :: this
  
  PetscErrorCode :: ierr
  
  class(stepper_BE_type), pointer :: flow_stepper
  class(stepper_BE_type), pointer :: tran_stepper

#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceFinalizeRun()')
#endif
  
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

subroutine SubsurfaceSimulationStrip(this)
  ! 
  ! Deallocates members of subsurface simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  implicit none
  
  class(subsurface_simulation_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceSimulationStrip()')
#endif
  
  call SimulationBaseStrip(this)
  call RealizationStrip(this%realization)
  deallocate(this%realization)
  nullify(this%realization)
  call RegressionDestroy(this%regression)
  
end subroutine SubsurfaceSimulationStrip

! ************************************************************************** !

subroutine SubsurfaceSimulationDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  class(subsurface_simulation_type), pointer :: simulation
  
#ifdef DEBUG
  call printMsg(simulation%option,'SimulationDestroy()')
#endif
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SubsurfaceSimulationDestroy
  
end module Subsurface_Simulation_class
