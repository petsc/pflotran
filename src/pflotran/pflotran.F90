!=======================================================================
! PFLOTRAN v1.0 LA-CC-06-093
!=======================================================================

! The software titled "PFLOTRAN v 1.0" has been assigned LA-CC-06-093. 
! The software is unclassified and does not contain Unclassified 
! Controlled Nuclear Information (UCNI).

! The software is under review to be released as open source, which 
! requires review and approval by the appropriate DOE Program Office. 
! If DOE declines to release this software as publicly available open 
! source, it would be subject to export control under Department of 
! Commerce regulations, and classified as ECCN 
! (Export Control Classification Number) EAR99. 
 
! If released as open source, the software will be publicly available and 
! not subject to export control. Until DOE approval is obtained, the 
! software should be treated as export controlled. DOE reserves the 
! right to release or deny release of the software as open source.

! Send all bug reports/questions/comments to:
! Peter C. Lichtner
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-6, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM
!=======================================================================
  program pflotran
  
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module
  use Stochastic_module
  use Stochastic_Aux_module  
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  PetscTruth :: truth
  PetscTruth :: option_found  
  PetscInt :: ierr
  character(len=MAXSTRINGLENGTH) :: string

  type(stochastic_type), pointer :: stochastic
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  option => OptionCreate()
  option%fid_out = IUNIT2

  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,option%global_rank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,option%global_commsize,ierr)
  call MPI_Comm_group(PETSC_COMM_WORLD,option%global_group,ierr)

  ! check for non-default input filename
  option_found = PETSC_FALSE
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflotranin", &
                             string, option_found, ierr)
  if(option_found) then
    option%input_filename = trim(string)
  else
    option%input_filename = "pflotran.in"
  endif
  
  ! check for screen output
  option_found = PETSC_FALSE
  call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-screen_output', &
                            truth,option_found, ierr)
  if (option_found) option%print_to_screen = truth

  ! check for file output
  option_found = PETSC_FALSE
  call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-file_output', &
                            truth,option_found, ierr)
  if (option_found) option%print_to_file = truth

  option_found = PETSC_FALSE
  call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-stochastic', &
                            truth,option_found, ierr)
  if (option_found) stochastic => StochasticCreate()

  call InitReadStochasticCardFromInput(stochastic,option)

  if (associated(stochastic)) then
    call StochasticInit(stochastic,option)
    call StochasticRun(stochastic,option)
  else

    option%mycomm = option%global_comm
    option%myrank = option%global_rank
    option%mycommsize = option%global_commsize
    option%mygroup = option%global_group
    
    call LoggingCreate()

    simulation => SimulationCreate(option)
    realization => simulation%realization

    call OptionCheckCommandLine(option)

    call PetscGetCPUTime(timex(1), ierr)
    call PetscGetTime(timex_wall(1), ierr)
    option%start_time = timex_wall(1)

    call Init(simulation)

    call StepperRun(simulation%realization,simulation%flow_stepper, &
                    simulation%tran_stepper)

  ! Clean things up.
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

    call LoggingDestroy()
    
  endif
  
  call OptionDestroy(option)
  call PetscFinalize (ierr)

  end program pflotran
