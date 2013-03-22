module Simulation_module

  use Process_Model_Coupler_module
  use Option_module
  use Output_Aux_module
  use Synchronizer_module
  use Process_Model_Coupler_module
  
  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    type(synchronizer_type), pointer :: synchronizer
    type(process_model_coupler_type), pointer :: process_model_coupler_list
  contains
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

  allocate(simulation)
  simulation%option => option
  nullify(simulation%output_option)
  nullify(simulation%synchronizer)
  nullify(simulation%process_model_coupler_list)
  
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
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%InitializeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo

  !TODO(geh): place logic here to stop if only initial state desired (e.g.
  !           solution composition, etc.).
  
  !TODO(geh): replace integer arguments with logical
  if (this%option%restart_flag) then
    call OutputInit(1) ! number greater than 0
  else
    call OutputInit(0)
  endif
  
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

  ! popped in TimeStepperFinalizeRun()
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
  
  call SynchronizerRunToTime(this%synchronizer,0.d0)
  
end subroutine ExecuteRun

! ************************************************************************** !
!
! FinalizeRun: Finalizes simulation
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine FinalizeRun(this)

  implicit none
  
  class(simulation_type) :: this
  
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    call cur_process_model_coupler%FinalizeRun()
    cur_process_model_coupler => cur_process_model_coupler%next
  enddo
  
end subroutine FinalizeRun

! ************************************************************************** !
!
! SimulationDestroy: Deallocates a simulation
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SimulationDestroy(simulation)

  use Richards_module, only : RichardsDestroy
  use Reactive_Transport_module, only : RTDestroy
  use General_module, only : GeneralDestroy

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return
  
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationDestroy
  
end module Simulation_module
