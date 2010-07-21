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

  use Simulation_module
  use Option_module
  use Input_module
  
  implicit none

  type(stochastic_type) :: stochastic
  type(option_type) :: option

  PetscInt :: i
  PetscInt :: offset, delta, remainder

  PetscInt :: realization_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: option_found
  PetscErrorCode :: ierr

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
  
  call SimulationCreateProcessorGroups(option,stochastic%num_groups)
  
  ! divvy up the realizations
  stochastic%num_local_realizations = stochastic%num_realizations / &
                                      stochastic%num_groups
  remainder = stochastic%num_realizations - stochastic%num_groups * &
                                            stochastic%num_local_realizations
  offset = 0
  do i = 1, option%mygroup_id-1
    delta = stochastic%num_local_realizations
    if (i < remainder) delta = delta + 1
    offset = offset + delta
  enddo
  
  if (option%mygroup_id < remainder) stochastic%num_local_realizations = &
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
  PetscInt :: status

  call OptionCheckCommandLine(option)

  ! moved outside due to errors when allocating/deallocating  over and over
  call LoggingCreate()

  do irealization = 1, stochastic%num_local_realizations

    call OptionInitRealization(option)
    simulation => SimulationCreate(option)
    realization => simulation%realization

    option%id = stochastic%realization_ids(irealization)
    write(string,'(i6)') option%id
    option%group_prefix = 'R' // trim(adjustl(string))

#if 0
    string = 'restart' // trim(adjustl(option%group_prefix)) // '.chk.info'
    open(unit=86,file=string,status="old",iostat=status)
    ! if file found, cycle
    if (status == 0) then
      close(86)
      cycle
    endif
#endif

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

  enddo
  
  ! moved outside due to errors when allocating/deallocating  over and over
  call LoggingDestroy()
    
  call MPI_Barrier(option%global_comm,ierr)

end subroutine StochasticRun

end module Stochastic_module
