module Simulation_module

  use Realization_class
  use Timestepper_module
  use Solver_module
  use Regression_module

#ifdef SURFACE_FLOW
  use Surface_Realization_class
#endif
#ifdef GEOMECH
  use Geomechanics_Realization_module
#endif
  implicit none

#include "definitions.h"
  
  private

  type, public :: simulation_type

    type(realization_type), pointer :: realization
    type(stepper_type), pointer :: flow_stepper
    type(stepper_type), pointer :: tran_stepper
#ifdef SURFACE_FLOW
    type(stepper_type), pointer :: surf_flow_stepper
#endif
#ifdef SURFACE_FLOW
    type(surface_realization_type), pointer :: surf_realization
#endif
#ifdef GEOMECH
    type(geomech_realization_type), pointer :: geomech_realization
#endif
    type(regression_type), pointer :: regression
  end type simulation_type
  
  interface SimulationCreate
    module procedure SimulationCreate1
    module procedure SimulationCreate2
  end interface
  
  public :: SimulationCreate, &
            SimulationDestroy, &
            SimulationResetTimeSteppers, &
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
#ifdef SURFACE_FLOW
  simulation%surf_flow_stepper => TimestepperCreate()
  simulation%surf_realization => SurfRealizCreate(option)
#endif
#ifdef GEOMECH
  simulation%geomech_realization => GeomechRealizCreate(option)
#endif
  nullify(simulation%regression)
  
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
! SimulationResetTimeSteppers: Sets time steppers back to initial settings
! author: Glenn Hammond
! date: 01/27/11
!
! ************************************************************************** !
subroutine SimulationResetTimeSteppers(simulation)

  use Timestepper_module

  implicit none

  type(simulation_type) :: simulation

  PetscReal :: dt_min
  PetscReal :: flow_dt_min = 0.d0
  PetscReal :: tran_dt_min = 0.d0
#ifdef SURFACE_FLOW
  PetscReal :: surf_flow_dt_min = 0.d0
#endif

  if (associated(simulation%flow_stepper)) &
    flow_dt_min = simulation%flow_stepper%dt_min
  if (associated(simulation%tran_stepper)) &
    tran_dt_min = simulation%tran_stepper%dt_min
#ifdef SURFACE_FLOW
  if (associated(simulation%surf_flow_stepper)) &
    surf_flow_dt_min = simulation%surf_flow_stepper%dt_min
#endif

  dt_min = max(flow_dt_min,tran_dt_min)
#ifdef SURFACE_FLOW
  dt_min = max(flow_dt_min,tran_dt_min,surf_flow_dt_min)
#endif

  simulation%realization%option%flow_time = 0.d0
  simulation%realization%option%flow_dt = dt_min
  simulation%realization%option%tran_time = 0.d0
  simulation%realization%option%tran_dt = dt_min
  simulation%realization%option%match_waypoint = PETSC_FALSE

  simulation%realization%output_option%plot_number = 0

  if (associated(simulation%flow_stepper)) then
    simulation%flow_stepper%cur_waypoint => &
      simulation%realization%waypoints%first
    call TimestepperReset(simulation%flow_stepper,dt_min)
  endif
  if (associated(simulation%tran_stepper)) then
    simulation%tran_stepper%cur_waypoint => &
      simulation%realization%waypoints%first
    call TimestepperReset(simulation%tran_stepper,dt_min)
  endif
#ifdef SURFACE_FLOW
  if (associated(simulation%surf_flow_stepper)) then
    simulation%surf_flow_stepper%cur_waypoint => &
      simulation%realization%waypoints%first
    call TimestepperReset(simulation%surf_flow_stepper,dt_min)

    simulation%surf_realization%option%flow_time = 0.d0
    simulation%surf_realization%option%flow_dt = dt_min
    simulation%surf_realization%option%tran_time = 0.d0
    simulation%surf_realization%option%tran_dt = dt_min
    simulation%surf_realization%option%match_waypoint = PETSC_FALSE

    simulation%surf_realization%output_option%plot_number = 0

    simulation%surf_flow_stepper%cur_waypoint => &
      simulation%surf_realization%waypoints%first
    call TimestepperReset(simulation%surf_flow_stepper,dt_min)
  endif

#endif

end subroutine SimulationResetTimeSteppers

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

  if (simulation%realization%option%nflowdof > 0) then
    select case(simulation%realization%option%iflowmode)
      case(RICHARDS_MODE)
        call RichardsDestroy(simulation%realization)
      case(G_MODE)
        call GeneralDestroy(simulation%realization)
    end select
  endif

  if (simulation%realization%option%ntrandof > 0) then
    call RTDestroy(simulation%realization)
  endif

  call RealizationDestroy(simulation%realization)
  call TimestepperDestroy(simulation%flow_stepper)
  call TimestepperDestroy(simulation%tran_stepper)
#ifdef SURFACE_FLOW
  call TimestepperDestroy(simulation%surf_flow_stepper)
#endif

#ifdef SURFACE_FLOW
  call SurfRealizDestroy(simulation%surf_realization)
#endif

#ifdef GEOMECH
  call GeomechRealizDestroy(simulation%geomech_realization)
#endif

  call RegressionDestroy(simulation%regression)

  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationDestroy
  
end module Simulation_module
