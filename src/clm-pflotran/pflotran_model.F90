module pflotran_model_module
#include "finclude/petsclog.h"
  use petscsys
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Input_module
  use Init_module
  use Logging_module
  use Stochastic_module
  use Stochastic_Aux_module
  use Waypoint_module
  use Units_module

  implicit none

#include "definitions.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  PetscTruth :: truth
  PetscTruth :: option_found
  PetscTruth :: single_inputfile
  PetscInt :: i
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  type, public :: pflotran_model_type
     type(stochastic_type),  pointer :: stochastic
     type(simulation_type),  pointer :: simulation
     type(realization_type), pointer :: realization
     type(option_type), pointer :: option
  end type pflotran_model_type


  public::pflotranModelCreate,               &
       pflotranModelStepperRunInit,          &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelUpdateTopBCHomogeneous,  &
       pflotranModelStepperRunFinalize,      &
       pflotranModelInsertWaypoint,          &
       pflotranModelDestroy

contains


  ! ************************************************************************** !
  !
  ! pflotranModelCreate: Allocates and initializes the pflotranModel object.
  !             It performs the same sequence of commands as done in pflotran.F90
  !     before model integration is performed by the call to StepperRun()
  !             routine
  !
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  function pflotranModelCreate()

    use Simulation_module
    use Realization_module
    use Timestepper_module
    use Option_module
    use Input_module
    use Init_module
    use Logging_module
    use Stochastic_module
    use Stochastic_Aux_module

    implicit none


    type(pflotran_model_type), pointer :: pflotranModelCreate


    PetscLogDouble :: timex(4), timex_wall(4)

    PetscTruth :: truth
    PetscTruth :: option_found
    PetscTruth :: single_inputfile
    PetscInt   :: i
    PetscInt   :: temp_int
    PetscErrorCode :: ierr
    character(len=MAXSTRINGLENGTH)          :: string
    character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

    type(pflotran_model_type),      pointer :: pflotran_model


    allocate(pflotran_model)
    allocate(pflotran_model%stochastic)
    allocate(pflotran_model%simulation)
    allocate(pflotran_model%realization)
    allocate(pflotran_model%option)

    nullify(pflotran_model%stochastic)
    pflotran_model%option => OptionCreate()
    pflotran_model%option%fid_out = IUNIT2
    single_inputfile = PETSC_TRUE


    call MPI_Init(ierr)
    pflotran_model%option%global_comm = MPI_COMM_WORLD
    call MPI_Comm_rank(MPI_COMM_WORLD,pflotran_model%option%global_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,pflotran_model%option%global_commsize,ierr)
    call MPI_Comm_group(MPI_COMM_WORLD,pflotran_model%option%global_group,ierr)
    pflotran_model%option%mycomm = pflotran_model%option%global_comm
    pflotran_model%option%myrank = pflotran_model%option%global_rank
    pflotran_model%option%mycommsize = pflotran_model%option%global_commsize
    pflotran_model%option%mygroup = pflotran_model%option%global_group

    ! check for non-default input filename
    pflotran_model%option%input_filename = "pflotran.in"
    string = '-pflotranin'
    call InputGetCommandLineString(string,pflotran_model%option%input_filename,option_found,pflotran_model%option)

    string = '-screen_output'
    call InputGetCommandLineTruth(string,pflotran_model%option%print_to_screen,option_found,pflotran_model%option)

    string = '-file_output'
    call InputGetCommandLineTruth(string,pflotran_model%option%print_to_file,option_found,pflotran_model%option)

    string = '-v'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) pflotran_model%option%verbosity = 1

    string = '-multisimulation'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) then
       single_inputfile = PETSC_FALSE
    endif

    string = '-stochastic'
    call InputGetCommandLineTruth(string,truth,option_found,pflotran_model%option)
    if (option_found) pflotran_model%stochastic => StochasticCreate()

    call InitReadStochasticCardFromInput(pflotran_model%stochastic,pflotran_model%option)

    if (associated(pflotran_model%stochastic)) then
       !  call StochasticInit(stochastic,option)
       !  call StochasticRun(stochastic,option)
    endif

    if (single_inputfile) then
       print *,'single_inputfile'
       PETSC_COMM_WORLD = MPI_COMM_WORLD
       call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    else
       call InitReadInputFilenames(pflotran_model%option,filenames)
       temp_int = size(filenames)
       call SimulationCreateProcessorGroups(pflotran_model%option,temp_int)
       pflotran_model%option%input_filename = filenames(pflotran_model%option%mygroup_id)
       i = index(pflotran_model%option%input_filename,'.',PETSC_TRUE)
       if (i > 1) then
          i = i-1
       else
          ! for some reason len_trim doesn't work on MS Visual Studio in
          ! this location
          i = len(trim(pflotran_model%option%input_filename))
       endif
       pflotran_model%option%global_prefix = pflotran_model%option%input_filename(1:i)
       write(string,*) pflotran_model%option%mygroup_id
       pflotran_model%option%group_prefix = 'G' // trim(adjustl(string))
    endif

    if (pflotran_model%option%verbosity > 0) then
       call PetscLogBegin(ierr)
       string = '-log_summary'
       call PetscOptionsInsertString(string, ierr)
    endif
    call LoggingCreate()

    pflotran_model%simulation => SimulationCreate(pflotran_model%option)
    pflotran_model%realization => pflotran_model%simulation%realization

    call OptionCheckCommandLine(pflotran_model%option)

    call PetscGetCPUTime(timex(1), ierr)
    call PetscGetTime(timex_wall(1), ierr)
    pflotran_model%option%start_time = timex_wall(1)

    call Init(pflotran_model%simulation)

    pflotranModelCreate => pflotran_model

  end function pflotranModelCreate


  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunInit: It performs the same execution of commands
  !             that are carried out in StepperRun() before the model integration
  !             begins over the entire simulation time interval
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunInit(pflotran_model)

    type(pflotran_model_type), pointer :: pflotran_model

    call StepperRunInit(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)

  end subroutine pflotranModelStepperRunInit

  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunTillPauseTime: It performs the model integration
  !             till the specified pause_time.
  !
  ! NOTE: It is assumed 'pause_time' is in seconds
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunTillPauseTime(pflotran_model, pause_time)

    type(pflotran_model_type), pointer :: pflotran_model
    PetscReal                          :: pause_time

    call pflotranModelInsertWaypoint(pflotran_model, pause_time)
    call pflotranModelInsertWaypoint(pflotran_model, pause_time + 1.0d0)

    call StepperRunOneDT(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper, pause_time)

  end subroutine pflotranModelStepperRunTillPauseTime


  ! ************************************************************************** !
  !
  !
  ! author: Gautam Bisht
  ! date: 9/22/2010
  ! ************************************************************************** !
  subroutine pflotranModelUpdateTopBCHomogeneous(pflotran_model, flux_value)

    use Condition_module


    type(pflotran_model_type), pointer        :: pflotran_model
    type(realization_type),pointer            :: realization
    type(condition_list_type),pointer         :: condition_list
    type(flow_condition_type), pointer        :: condition
    type(flow_sub_condition_type), pointer    :: sub_condition
    type(flow_condition_dataset_type),pointer :: dataset
    type(option_type),pointer                 :: option


    PetscReal                                 :: flux_value
    PetscInt                                  :: isub_condition

    realization => pflotran_model%realization
    condition_list => realization%flow_conditions
    option => realization%option

    condition => condition_list%first
    do
       if (.not.associated(condition)) exit

       if (trim(condition%name) == 'top') then

          do isub_condition = 1, condition%num_sub_conditions

             sub_condition => condition%sub_condition_ptr(isub_condition)%ptr

             if (associated(sub_condition)) then
                dataset => sub_condition%dataset
                dataset%values(1,1) = flux_value
                print *,'In UpdateTopBCHomogeneous: ', flux_value
                !call FlowSubConditionUpdateDataset(option,1.0d0,sub_condition%dataset)
             endif

          enddo
       endif

       condition => condition%next

    enddo

    call StepperUpdateSolution(realization)

  end subroutine pflotranModelUpdateTopBCHomogeneous

  ! ************************************************************************** !
  !
  ! pflotranModelStepperRunFinalize: It performs the same execution of commands
  !             that are carried out in StepperRun() once the model integration is
  !             finished
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelStepperRunFinalize(pflotran_model)

    type(pflotran_model_type), pointer :: pflotran_model

    call StepperRunFinalize(pflotran_model%simulation%realization, &
         pflotran_model%simulation%flow_stepper, &
         pflotran_model%simulation%tran_stepper)

  end subroutine pflotranModelStepperRunFinalize



  ! ************************************************************************** !
  !
  ! pflotranModelInsertWaypoint: Inserts a waypoint within the waypoint list
  !             so that the model integration can be paused when that waypoint is
  !             reached
  !
  ! NOTE: It is assumed the 'waypoint_time' is in seconds
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelInsertWaypoint(pflotran_model, waypoint_time)

    type(pflotran_model_type),pointer :: pflotran_model
    type(waypoint_type)      ,pointer :: waypoint
    type(option_type)        ,pointer :: option
    PetscReal                         :: waypoint_time
    character(len=MAXWORDLENGTH)      :: word

    option => pflotran_model%realization%option
    word = 's'
    waypoint => WaypointCreate()
    waypoint%time        = waypoint_time * UnitsConvertToInternal(word,option)
    waypoint%update_srcs = PETSC_TRUE
    waypoint%dt_max      = 3153600

    call WaypointInsertInList(waypoint,pflotran_model%realization%waypoints)


  end subroutine pflotranModelInsertWaypoint

  ! ************************************************************************** !
  !
  ! pflotranModelDestroy: Deallocates the pflotranModel object
  !
  ! author: Gautam Bisht
  ! date: 9/10/2010
  ! ************************************************************************** !
  subroutine pflotranModelDestroy(pflotran_model)

    type(pflotran_model_type), pointer :: pflotran_model

    ! Clean things up.
    call SimulationDestroy(pflotran_model%simulation)

    ! Final Time
    call PetscGetCPUTime(timex(2), ierr)
    call PetscGetTime(timex_wall(2), ierr)

    if (pflotran_model%option%myrank == pflotran_model%option%io_rank) then

       if (pflotran_model%option%print_to_screen) then
          write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
               (timex(2)-timex(1))/3600.d0

          write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
               (timex_wall(2)-timex_wall(1))/3600.d0
       endif
       if (pflotran_model%option%print_to_file) then
          write(pflotran_model%option%fid_out,'(/," CPU Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
               (timex(2)-timex(1))/3600.d0

          write(pflotran_model%option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
               & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
               timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
               (timex_wall(2)-timex_wall(1))/3600.d0
       endif
    endif

    if (pflotran_model%option%myrank == pflotran_model%option%io_rank .and. pflotran_model%option%print_to_file) &
         close(pflotran_model%option%fid_out)

    call LoggingDestroy()

    call PetscOptionsSetValue('-options_left','no',ierr);

    call OptionDestroy(pflotran_model%option)
    call PetscFinalize(ierr)
    call MPI_Finalize(ierr)

  end subroutine pflotranModelDestroy

end module pflotran_model_module
