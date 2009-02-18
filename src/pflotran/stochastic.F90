module Stochastic_module

  use Stochastic_Aux_module

  implicit none
  
  private
  
#include "definitions.h"
  
  public :: StochasticInit, &
            StochasticRun
  
contains

! ************************************************************************** !
!
! StochasticInit: Initializes a stochastic simulation
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticInit(stochastic,option)

  use Option_module
  use Input_module
  
  implicit none

  type(stochastic_type) :: stochastic
  type(option_type) :: option

  PetscMPIInt :: mycolor, mykey

  PetscInt :: i
  PetscInt :: igroup, irealization
  PetscInt :: num_groups
  PetscInt :: local_commsize, offset, delta, remainder

  PetscInt :: realization_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: option_found
  PetscErrorCode :: ierr

#if 0
  ! set up global and local communicator groups, processor ranks, and group sizes
  call MPI_Init(ierr)
  global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,global_group,ierr)
  PETSC_COMM_WORLD = global_comm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

#if 1
  ! query user for number of communicator groups and realizations
  string = '-num_groups'
  call InputGetCommandLineInt(string,stochastic%num_groups,option_found,option)

  string = '-num_realizations'
  call InputGetCommandLineInt(string,stochastic%num_realizations,option_found,option)

  ! error checking
  if (stochastic%num_groups == 0) then
    option%io_buffer = 'Number of stochastic processor groups not ' // &
                       'initialized. Setting to 1.'
    call printWrnMsg(option)
    stochastic%num_groups = 1
  endif
  if (stochastic%num_realizations == 0) then
    option%io_buffer = 'Number of stochastic realizations not ' // &
                       'initialized. Setting to 1.'
    call printWrnMsg(option)
    stochastic%num_realizations = 1
  endif
#endif
  
  local_commsize = option%global_commsize / stochastic%num_groups
  remainder = option%global_commsize - stochastic%num_groups * local_commsize
  offset = 0
  do igroup = 1, stochastic%num_groups
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

  ! divvy up the realizations
  stochastic%num_local_realizations = stochastic%num_realizations / &
                                      stochastic%num_groups
  remainder = stochastic%num_realizations - stochastic%num_groups * &
                                            stochastic%num_local_realizations
  offset = 0
  do i = 1, igroup-1
    delta = stochastic%num_local_realizations
    if (i < remainder) delta = delta + 1
    offset = offset + delta
  enddo
  
  if (igroup < remainder) stochastic%num_local_realizations = &
                          stochastic%num_local_realizations + 1
  allocate(stochastic%realization_ids(stochastic%num_local_realizations))
  stochastic%realization_ids = 0
  do i = 1, stochastic%num_local_realizations
    stochastic%realization_ids(i) = offset + i
  enddo

end subroutine StochasticInit

! ************************************************************************** !
!
! StochasticRun: Runs a stochastic simulation
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticRun(stochastic,option)

  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module

  implicit none

#include "finclude/petsclog.h"
  
  type(stochastic_type), pointer :: stochastic
  type(option_type), pointer :: option

  PetscLogDouble :: timex(4), timex_wall(4)
  PetscInt :: irealization
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  call OptionCheckCommandLine(option)

  do irealization = 1, stochastic%num_local_realizations

    call LoggingCreate()

    call OptionInitRealization(option)
    simulation => SimulationCreate(option)
    realization => simulation%realization

    option%id = stochastic%realization_ids(irealization)
    write(string,'(i6)') option%id
    option%group_prefix = 'R' // trim(adjustl(string))

    call PetscGetCPUTime(timex(1), ierr)
    call PetscGetTime(timex_wall(1), ierr)
    option%start_time = timex_wall(1)

    call Init(simulation)

    call StepperRun(simulation%realization,simulation%flow_stepper, &
                    simulation%tran_stepper)

    call SimulationDestroy(simulation)

  ! Final Time
    call PetscGetCPUTime(timex(2), ierr)
    call PetscGetTime(timex_wall(2), ierr)
    
    if (option%myrank == option%io_rank) then

      if (option%print_to_screen) then
        write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
          (timex(2)-timex(1))/3600.d0

        write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
          (timex_wall(2)-timex_wall(1))/3600.d0
      endif
      if (option%print_to_file) then
        write(option%fid_out,'(/," CPU Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
          (timex(2)-timex(1))/3600.d0

        write(option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
          (timex_wall(2)-timex_wall(1))/3600.d0
      endif
    endif

    if (option%myrank == option%io_rank .and. option%print_to_file) &
      close(option%fid_out)
    if (option%myrank == option%io_rank .and. mod(irealization,10) == 0) then
      write(string,'(i6)') option%id
      print *, 'Finished with ' // trim(adjustl(string)), irealization, &
               ' of ', stochastic%num_local_realizations
    endif


    call LoggingDestroy()
    
  enddo
  
  call MPI_Barrier(option%global_comm,ierr)

end subroutine StochasticRun

end module Stochastic_module
