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
  use Init_module
  
  implicit none

  type(stochastic_type) :: stochastic
  type(option_type) :: option

  PetscInt :: i
  PetscInt :: offset, delta, remainder

  PetscInt :: realization_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscInt, pointer :: realization_ids_from_file(:)
  character(len=MAXSTRINGLENGTH) :: filename
  type(input_type), pointer :: input
  PetscErrorCode :: ierr
  
  call StochasticReadCardFromInput(stochastic,option)

  ! query user for number of communicator groups and realizations
  string = '-num_groups'
  call InputGetCommandLineInt(string,stochastic%num_groups, &
                              option_found,option)

  string = '-num_realizations'
  call InputGetCommandLineInt(string,stochastic%num_realizations, &
                              option_found,option)

  ! read realization ids from a file - contributed by Xingyuan
  string = '-realization_ids_file'
  call InputGetCommandLineString(string,filename,option_found,option)
  if (option_found) then
    input => InputCreate(IUNIT_TEMP,filename,option)
    allocate(realization_ids_from_file(stochastic%num_realizations))
    realization_ids_from_file = 0
    string = &
      '# of realization ids read from file may be too few in StochasticInit()'
    do i = 1, stochastic%num_realizations
      call InputReadFlotranString(input,option)
      call InputReadStringErrorMsg(input,option,string)
      call InputReadInt(input,option,realization_ids_from_file(i))
      call InputErrorMsg(input,option,'realization id', &
                         'StochasticInit')
    enddo
    call InputDestroy(input)
  else
    nullify(realization_ids_from_file)
  endif
    
  ! Realization offset contributed by Xingyuan.  This allows one to specify the
  ! smallest/lowest realization id (other than zero) in a stochastic simulation
  string = '-realization_offset'
  call InputGetCommandLineInt(string,offset,option_found,option)
  if (.not.option_found) then
    offset = 0
  endif

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
  
  call InitCreateProcessorGroups(option,stochastic%num_groups)
  
  ! divvy up the realizations
  stochastic%num_local_realizations = stochastic%num_realizations / &
                                      stochastic%num_groups
  remainder = stochastic%num_realizations - stochastic%num_groups * &
                                            stochastic%num_local_realizations
  
  ! offset is initialized above after check for '-realization_offset'
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
  
  ! map ids from file - contributed by Xingyuan
  if (associated(realization_ids_from_file)) then
    do i = 1, stochastic%num_local_realizations
      stochastic%realization_ids(i) = &
        realization_ids_from_file(stochastic%realization_ids(i))
    enddo
  endif

end subroutine StochasticInit

! ************************************************************************** !
!
! StochasticReadCardFromInput: Reads stochastic card from input file
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticReadCardFromInput(stochastic,option)

  use Option_module
  use Input_module
  use Stochastic_Aux_module

  implicit none
  
  type(stochastic_type) :: stochastic
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  type(input_type), pointer :: input
  PetscBool :: print_warning
  
  input => InputCreate(IN_UNIT,option%input_filename,option)

  ! MODE information
  string = "STOCHASTIC"
  print_warning = PETSC_FALSE
  call InputFindStringInFile(input,option,string,print_warning)

  if (.not.InputError(input)) then
    call StochasticRead(stochastic,input,option)
  endif
  
  call InputDestroy(input)

end subroutine StochasticReadCardFromInput

! ************************************************************************** !
!
! StochasticRun: Runs a stochastic simulation
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticRun(stochastic,option)

  use Simulation_module
  use Realization_class
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module
  use Regression_module

  implicit none

#include "finclude/petsclog.h"
  
  type(stochastic_type), pointer :: stochastic
  type(option_type), pointer :: option

  PetscLogDouble :: timex_wall(4)
  PetscInt :: irealization
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(stepper_type), pointer :: master_stepper
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: init_status
  PetscErrorCode :: ierr
  PetscInt :: status
  
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

    call PetscTime(timex_wall(1), ierr)
    option%start_time = timex_wall(1)

    call Init(simulation)

#ifdef SURFACE_FLOW
    !call StepperRun(simulation%realization,simulation%flow_stepper, &
    !                simulation%tran_stepper,simulation%surf_flow_stepper)
    option%io_buffer = 'Stochastic mode not tested for surface-flow'
    call printErrMsgByRank(option)
#else
    call TimestepperInitializeRun(simulation%realization, &
                                  master_stepper, &
                                  simulation%flow_stepper, &
                                  simulation%tran_stepper, &
                                  init_status)
    select case(init_status)
      case(TIMESTEPPER_INIT_PROCEED)
        call  TimestepperExecuteRun(simulation%realization, &
                                    master_stepper, &
                                    simulation%flow_stepper, &
                                    simulation%tran_stepper)
        call  TimestepperFinalizeRun(simulation%realization, &
                                     master_stepper, &
                                     simulation%flow_stepper, &
                                     simulation%tran_stepper)
      case(TIMESTEPPER_INIT_FAIL)
      case(TIMESTEPPER_INIT_DONE)
    end select
#endif

    call RegressionOutput(simulation%regression,simulation%realization, &
                          simulation%flow_stepper,simulation%tran_stepper)

    call SimulationDestroy(simulation)

  ! Final Time
    call PetscTime(timex_wall(2), ierr)
    
    if (option%myrank == option%io_rank) then

      if (option%print_to_screen) then
        write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
        & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
          timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
          (timex_wall(2)-timex_wall(1))/3600.d0
      endif
      if (option%print_to_file) then
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
  
  call MPI_Barrier(option%global_comm,ierr)

end subroutine StochasticRun

end module Stochastic_module
