subroutine f_initialize_simulation(simulation)

  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module

  implicit none

#include "definitions.h"
#include "include/finclude/petsclog.h"

  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  PetscLogDouble :: timex(4), timex_wall(4)

  PetscMPIInt :: myrank, commsize
  PetscInt :: ierr
  PetscInt :: stage(10)
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflotranin

  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)
  
  call LoggingCreate()
  realization => simulation%realization
  option => realization%option

  option%comm    = PETSC_COMM_WORLD
  option%myrank = myrank
  option%commsize = commsize

  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflotranin", &
                             pflotranin, option_found, ierr)
  if(.not.option_found) pflotranin = "pflotran.in"
  
  call PetscGetCPUTime(timex(1), ierr)
  call PetscGetTime(timex_wall(1), ierr)
  option%start_time = timex_wall(1)

  call OptionCheckCommandLine(option)

  call Init(simulation,pflotranin)

end subroutine f_initialize_simulation
