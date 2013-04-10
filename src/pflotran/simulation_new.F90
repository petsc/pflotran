module Simulation_module

  use Process_Model_Coupler_module
  use Option_module
  use Output_Aux_module
  use Output_module
  use Synchronizer_module
  use Process_Model_Coupler_module
  use Regression_module
  
  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    type(synchronizer_type), pointer :: synchronizer
    type(regression_type), pointer :: regression    
    type(process_model_coupler_type), pointer :: process_model_coupler_list
  contains
    procedure, public :: Initialize
    procedure, public :: InitializeRun
    procedure, public :: ExecuteRun
    procedure, public :: FinalizeRun
  end type simulation_type
  
  public :: SimulationCreate, &
            SimulationDestroy, &
            SimulationCreateProcessorGroups
  
contains

! ************************************************************************** !
!
! SimulationCreate1: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate(option)

  use Option_module
  
  implicit none
  
  type(simulation_type), pointer :: SimulationCreate
  
  type(simulation_type), pointer :: simulation
  type(option_type), pointer :: option

  print *, 'SimulationCreate'
  
  allocate(simulation)
  simulation%option => option
  nullify(simulation%output_option)
  nullify(simulation%synchronizer)
  nullify(simulation%process_model_coupler_list)
  nullify(simulation%regression)  
  
  SimulationCreate => simulation
  
end function SimulationCreate

! ************************************************************************** !
!
! SimulationCreateProcessorGroups: Splits MPI_COMM_WORLD into N separate
!                                  processor groups
! author: Glenn Hammond
! date: 08/11/09
!
! ************************************************************************** !
subroutine SimulationCreateProcessorGroups(option,num_groups)

  use Option_module

  type(option_type) :: option
  PetscInt :: num_groups

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscErrorCode :: ierr

  local_commsize = option%global_commsize / num_groups
  remainder = option%global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (option%global_rank >= offset .and. &
        option%global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  option%mygroup_id = igroup
  mykey_mpi = option%global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor_mpi,mykey_mpi,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)

  PETSC_COMM_WORLD = option%mycomm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)

end subroutine SimulationCreateProcessorGroups

subroutine Initialize(this)

  use Logging_module
  use Process_Model_Coupler_module
  use Process_Model_Richards_class
  use Process_Model_RT_class
  use Process_Model_Base_class
  
  use Waypoint_module
  use Timestepper_module

  implicit none
  
  class(simulation_type) :: this

  type(process_model_coupler_type), pointer :: richards_process_model_coupler
  type(process_model_coupler_type), pointer :: rt_process_model_coupler
  class(process_model_base_type), pointer :: process_model
  type(synchronizer_type), pointer :: synchronizer
  type(output_option_type), pointer :: output_option
  type(waypoint_type), pointer :: waypoint
  
  type(waypoint_list_type), pointer :: waypoint_list
  type(stepper_type), pointer :: timestepper
  PetscErrorCode :: ierr
  
  ! Richards Flow
  
  ! popped in InitializeRun()
  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr)

  output_option => OutputOptionCreate()

  process_model => PMRichardsCreate()
  process_model%option => this%option

  waypoint_list => WaypointListCreate()
  waypoint => WaypointCreate()
  waypoint%time = 0.d0
  waypoint%dt_max = 1.d0
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 3.9d0
  waypoint%dt_max = 2.1d0
  waypoint%update_conditions = PETSC_TRUE
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 11.d0
  waypoint%dt_max = 3.d0
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 20.d0
  waypoint%dt_max = 999.d0
  waypoint%final = PETSC_TRUE
  call WaypointInsertInList(waypoint,waypoint_list)
  call WaypointListFillIn(this%option,waypoint_list)
  
  timestepper => TimestepperCreate()
  timestepper%cur_waypoint => waypoint_list%first

  richards_process_model_coupler => ProcessModelCouplerCreate()
  richards_process_model_coupler%option => this%option
  richards_process_model_coupler%process_model_list => process_model
  richards_process_model_coupler%waypoints => waypoint_list
  richards_process_model_coupler%timestepper => timestepper
  
  nullify(waypoint)
  nullify(waypoint_list)

  ! Reactive Transport
  process_model => PMRTCreate()
  process_model%option => this%option

  waypoint_list => WaypointListCreate()
  waypoint => WaypointCreate()
  waypoint%time = 0.d0
  waypoint%dt_max = 0.75d0
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 2.5d0
  waypoint%dt_max = 1.3d0
  waypoint%update_conditions = PETSC_TRUE
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 11.1d0
  waypoint%dt_max = 3.2d0
  call WaypointInsertInList(waypoint,waypoint_list)
  waypoint => WaypointCreate()
  waypoint%time = 20.d0
  waypoint%dt_max = 999.d0
  waypoint%final = PETSC_TRUE
  call WaypointInsertInList(waypoint,waypoint_list)
  call WaypointListFillIn(this%option,waypoint_list)
  
  timestepper => TimestepperCreate()
  timestepper%cur_waypoint => waypoint_list%first

  rt_process_model_coupler => ProcessModelCouplerCreate()
  rt_process_model_coupler%option => this%option
  rt_process_model_coupler%process_model_list => process_model
  rt_process_model_coupler%waypoints => waypoint_list
  rt_process_model_coupler%timestepper => timestepper  
  
  nullify(waypoint)
  nullify(waypoint_list)

  synchronizer => SynchronizerCreate()
  synchronizer%waypoints => WaypointListCreate()
  waypoint => WaypointCreate()
  waypoint%time = 0.d0
  call WaypointInsertInList(waypoint,synchronizer%waypoints)
  waypoint => WaypointCreate()
  waypoint%time = 10.d0
  call WaypointInsertInList(waypoint,synchronizer%waypoints)
  waypoint => WaypointCreate()
  waypoint%final = PETSC_TRUE
  waypoint%time = 20.d0
  call WaypointInsertInList(waypoint,synchronizer%waypoints)
  call WaypointListFillIn(this%option,synchronizer%waypoints)

  synchronizer%option => this%option
  synchronizer%output_option => output_option

  this%output_option => output_option
  this%synchronizer => synchronizer
  this%synchronizer%process_model_coupler_list => richards_process_model_coupler
  richards_process_model_coupler%below => rt_process_model_coupler
  this%process_model_coupler_list => richards_process_model_coupler
  
end subroutine Initialize

! ************************************************************************** !
!
! InitializeRun: Initializes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine InitializeRun(this)

  use Logging_module

  implicit none
  
  class(simulation_type) :: this

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
  
end subroutine InitializeRun

! ************************************************************************** !
!
! ExecuteRun: Initializes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine ExecuteRun(this)

  implicit none
  
  class(simulation_type) :: this
  PetscInt :: stop_flag
  
  call printMsg(this%option,'Simulation%ExecuteRun()')

  call this%synchronizer%ExecuteRun(stop_flag)
  
end subroutine ExecuteRun

! ************************************************************************** !
!
! FinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine FinalizeRun(this)

  use Realization_class
  use Timestepper_module
  use Process_Model_Base_class
  use Process_Model_Richards_class
  use Process_Model_TH_class
  use Process_Model_RT_class

  implicit none
  
  class(simulation_type) :: this
  
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
  
end subroutine FinalizeRun

! ************************************************************************** !
!
! SimulationDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SimulationDestroy(simulation)

#ifndef SIMPLIFY
  use Richards_module, only : RichardsDestroy
  use Reactive_Transport_module, only : RTDestroy
  use General_module, only : GeneralDestroy
#endif

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  
  call RegressionDestroy(simulation%regression)  
  
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationDestroy
  
end module Simulation_module
