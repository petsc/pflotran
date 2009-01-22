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
  
  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  PetscLogDouble :: timex(4), timex_wall(4)

  PetscMPIInt :: global_rank, global_commsize, global_comm, global_group
  PetscMPIInt :: myrank, mycommsize, mycomm, mygroup
  PetscMPIInt :: mycolor, mykey
  PetscMPIInt :: io_rank

  PetscInt :: i
  PetscInt :: out_unit
  PetscInt :: igroup, irealization
  PetscInt :: num_groups
  PetscInt :: local_commsize, offset, delta, remainder

  PetscInt :: num_realizations
  PetscInt :: num_local_realizations
  PetscInt, allocatable :: realization_ids(:)

  PetscInt :: ierr
  PetscInt :: stage(10)
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflotranin
  character(len=MAXWORDLENGTH) :: string

  
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  num_groups = 1
  num_realizations = 1
#ifdef GLENN
  ! set up global and local communicator groups, processor ranks, and group sizes
  call MPI_Init(ierr)
  global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,global_group,ierr)
  local_commsize = global_commsize / num_groups
  remainder = global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (global_rank >= offset .and. global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor = igroup
  mykey = global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor,mykey,mycomm,ierr)
  call MPI_Comm_group(mycomm,mygroup,ierr)
  PETSC_COMM_WORLD = mycomm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(mycomm,myrank, ierr)
  call MPI_Comm_size(mycomm,mycommsize,ierr)

  ! divvy up the realizations
  num_local_realizations = num_realizations / num_groups
  remainder = num_realizations - num_groups * num_local_realizations
  offset = 0
  do i = 1, igroup-1
    delta = num_local_realizations
    if (i < remainder) delta = delta + 1
    offset = offset + delta
  enddo
  
  if (igroup < remainder) num_local_realizations = num_local_realizations + 1
  allocate(realization_ids(num_local_realizations))
  realization_ids = 0
  do i = 1, num_local_realizations
    realization_ids(i) = offset + i
  enddo

#else  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  global_comm = PETSC_COMM_WORLD
  call MPI_Comm_rank(PETSC_COMM_WORLD,global_rank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,global_commsize,ierr)
  call MPI_Comm_group(PETSC_COMM_WORLD,global_group,ierr)
  mycomm = global_comm
  myrank = global_rank
  mycommsize = global_commsize
  mygroup = global_group
  num_local_realizations = num_realizations
#endif  
  
  do irealization = 1, num_local_realizations

    call LoggingCreate()

    simulation => SimulationCreate()
    realization => simulation%realization
    option => realization%option

#ifdef GLENN    
    option%id = realization_ids(irealization)
#endif

    option%fid_out = IUNIT2
    out_unit = option%fid_out

    option%global_comm = global_comm
    option%global_rank = global_rank
    option%global_commsize = global_commsize
    option%global_group = global_group

    option%mycomm = mycomm
    option%myrank = myrank
    option%mycommsize = mycommsize
    option%mygroup = mygroup

#ifdef GLENN
    if (num_realizations > 1) then
      write(string,'(i6)') realization_ids(irealization)
      option%group_prefix = 'R' // trim(adjustl(string)) // '_'
    endif
#endif

    call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflotranin", &
                               pflotranin, option_found, ierr)
    if(.not.option_found) pflotranin = "pflotran.in"
    
    call OptionCheckCommandLine(option)

    call PetscGetCPUTime(timex(1), ierr)
    call PetscGetTime(timex_wall(1), ierr)
    option%start_time = timex_wall(1)

    call Init(simulation,pflotranin)

    call StepperRun(simulation%realization,simulation%flow_stepper, &
                    simulation%tran_stepper)

  ! Clean things up.
    io_rank = option%io_rank
    call SimulationDestroy(simulation)

  ! Final Time
    call PetscGetCPUTime(timex(2), ierr)
    call PetscGetTime(timex_wall(2), ierr)
    
    if (myrank == io_rank) then

      write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

      write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0

      write(out_unit,'(/," CPU Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
        (timex(2)-timex(1))/3600.d0

      write(out_unit,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
        (timex_wall(2)-timex_wall(1))/3600.d0
    endif

    if (myrank == io_rank) close(out_unit)

    call LoggingDestroy()
    
  enddo
  
  call PetscFinalize (ierr)
  deallocate(realization_ids)

  end program pflotran
