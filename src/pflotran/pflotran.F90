!=======================================================================
! PFLOTRAN v1.0 LA-CC-06-093
!=======================================================================

! The software titled "PFLOTRAN v 1.0" has been assigned LA-CC-06-093. 
! The software is unclassified and does not contain Unclassified 
! Controlled Nuclear Information (UCNI).

! The software titled "PFLOTRAN v 2.0" has been assigned LA-CC-09-047. 
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
  use Input_module
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
  PetscTruth :: single_inputfile
  PetscInt :: i
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)
  type(stochastic_type), pointer :: stochastic
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  nullify(stochastic)
  option => OptionCreate()
  option%fid_out = IUNIT2
  single_inputfile = PETSC_TRUE

  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,option%global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,option%global_group,ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group
#ifdef VAMSI_HDF5
  option%broadcast_size = HDF5_BROADCAST_SIZE
  option%color = option%global_rank / option%broadcast_size
  option%key = option%global_rank
  call MPI_Comm_split(option%global_comm,option%color,option%key,option%iogroup,ierr)
  call MPI_Comm_size(option%iogroup,option%localsize,ierr)
  call MPI_Comm_rank(option%iogroup,option%localrank,ierr)
  if (mod(option%global_rank,option%broadcast_size) == 0) then
	option%reader_color = 1
	option%reader_key = option%global_rank
  else
    option%reader_color = 0
    option%reader_key = option%global_rank
  endif
  call MPI_Comm_split(option%global_comm,option%reader_color,option%reader_key,option%readers,ierr)
  call MPI_Comm_size(option%readers,option%reader_size,ierr)
  call MPI_Comm_rank(option%readers,option%reader_rank,ierr)
#endif
  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename,option_found,option)

  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-v'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) option%verbosity = 1

  string = '-multisimulation'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) then
    single_inputfile = PETSC_FALSE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,truth,option_found,option)
  if (option_found) stochastic => StochasticCreate()

  call InitReadStochasticCardFromInput(stochastic,option)

  if (associated(stochastic)) then
    call StochasticInit(stochastic,option)
    call StochasticRun(stochastic,option)
  else

    if (single_inputfile) then
      PETSC_COMM_WORLD = MPI_COMM_WORLD
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    else
      call InitReadInputFilenames(option,filenames)
      call SimulationCreateProcessorGroups(option,size(filenames))
      option%input_filename = filenames(option%mygroup_id)
      i = index(option%input_filename,'.',PETSC_TRUE)
      if (i > 1) then
        i = i-1
      else
        ! for some reason len_trim doesn't work on MS Visual Studio in 
        ! this location
        i = len(trim(option%input_filename)) 
      endif
      option%global_prefix = option%input_filename(1:i)
      write(string,*) option%mygroup_id
      option%group_prefix = 'G' // trim(adjustl(string))
    endif
 
    if (option%verbosity > 0) then 
      call PetscLogBegin(ierr)
      string = '-log_summary'
      call PetscOptionsInsertString(string, ierr)
    endif
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
  
  call PetscOptionsSetValue('-options_left','no',ierr);

  call OptionDestroy(option)
  call PetscFinalize (ierr)
  call MPI_Finalize (ierr)

  end program pflotran
