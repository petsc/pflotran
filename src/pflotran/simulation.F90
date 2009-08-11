module Simulation_module

  use Realization_module
  use Timestepper_module
  use Solver_module

  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type

    type(realization_type), pointer :: realization
    type(stepper_type), pointer :: flow_stepper
    type(stepper_type), pointer :: tran_stepper

  end type simulation_type
  
  interface SimulationCreate
    module procedure SimulationCreate1
    module procedure SimulationCreate2
  end interface
  
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
function SimulationCreate1()

  use Option_module
  
  implicit none
  
  type(simulation_type), pointer :: SimulationCreate1
  
  type(simulation_type), pointer :: simulation
  type(option_type), pointer :: option
  
  nullify(option)
  SimulationCreate1 => SimulationCreate2(option)
  
end function SimulationCreate1

! ************************************************************************** !
!
! SimulationCreate2: Allocates and initializes a new simulation object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SimulationCreate2(option)

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  
  type(simulation_type), pointer :: SimulationCreate2
  
  type(simulation_type), pointer :: simulation
  
  allocate(simulation)
  simulation%realization => RealizationCreate(option)
  simulation%flow_stepper => TimestepperCreate()
  simulation%tran_stepper => TimestepperCreate()
  
  SimulationCreate2 => simulation
  
end function SimulationCreate2

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
  PetscMPIInt :: mycolor, mykey  
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
  mycolor = igroup
  option%mygroup_id = igroup
  mykey = option%global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor,mykey,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)

  PETSC_COMM_WORLD = option%mycomm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)

end subroutine SimulationCreateProcessorGroups

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

  implicit none
  
  type(simulation_type), pointer :: simulation
  
  if (.not.associated(simulation)) return

  if (simulation%realization%option%nflowdof > 0) then
    select case(simulation%realization%option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsDestroy(simulation%realization)
    end select
  endif

  if (simulation%realization%option%ntrandof > 0) then
    call RTDestroy(simulation%realization)
  endif

  call RealizationDestroy(simulation%realization)
  call TimestepperDestroy(simulation%flow_stepper)
  call TimestepperDestroy(simulation%tran_stepper)
  
end subroutine SimulationDestroy
  
end module Simulation_module
