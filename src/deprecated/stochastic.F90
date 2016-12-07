module Stochastic_module


#include "finclude/petscsys.h"
  use petscsys
  use Stochastic_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  

  
  public :: StochasticInit, &
            StochasticRun
  
contains

! ************************************************************************** !

subroutine StochasticInit(stochastic,option)
  ! 
  ! Initializes a stochastic simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 

  use Simulation_module
  use Option_module
  use Input_Aux_module
  
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
      call InputReadPflotranString(input,option)
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
  
  call OptionCreateProcessorGroups(option,stochastic%num_groups)
  
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

subroutine StochasticReadCardFromInput(stochastic,option)
  ! 
  ! Reads stochastic card from input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 

  use Option_module
  use Input_Aux_module
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

subroutine StochasticRun(stochastic,option)
  ! 
  ! Runs a stochastic simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 
#include "finclude/petsclog.h"
  use petscsys
  use Simulation_module
  use Realization_class
  use Timestepper_module
  use Option_module
  use Init_module
  use PFLOTRAN_Factory_module
  use Logging_module
  use Geomechanics_Logging_module

  implicit none

  
  type(stochastic_type), pointer :: stochastic
  type(option_type), pointer :: option

  type(timestepper_type), pointer :: master_timestepper
  PetscInt :: irealization
  type(simulation_type), pointer :: simulation
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: init_status
  PetscErrorCode :: ierr
  PetscInt :: status
  
  call OptionInitPetsc(option)
  call LoggingCreate()
  call GeomechLoggingCreate()

  do irealization = 1, stochastic%num_local_realizations

    call OptionInitRealization(option)

    ! Set group prefix based on id
    option%id = stochastic%realization_ids(irealization)
    write(string,'(i6)') option%id
    option%group_prefix = 'R' // trim(adjustl(string))

#if 0
    ! code for restarting stochastic runs; may no longer need this.
    string = 'restart' // trim(adjustl(option%group_prefix)) // '.chk.info'
    open(unit=86,file=string,status="old",iostat=status)
    ! if file found, cycle
    if (status == 0) then
      close(86)
      cycle
    endif
#endif

    option%io_buffer = 'Stochastic mode not tested for surface-flow'
    call printErrMsgByRank(option)

    call PFLOTRANInitializePostPETSc(simulation,master_timestepper,option, &
                                     init_status)
    call PFLOTRANRun(simulation,master_timestepper,init_status)
    call PFLOTRANFinalize(simulation,option)

    call SimulationDestroy(simulation)

  enddo
  
  call LoggingDestroy()
  call GeomechLoggingDestroy()
  call MPI_Barrier(option%global_comm,ierr)

end subroutine StochasticRun

end module Stochastic_module
