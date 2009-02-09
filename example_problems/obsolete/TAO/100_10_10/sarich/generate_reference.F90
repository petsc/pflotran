  program generate_reference
      
  use Simulation_module
  use Realization_module
  use Timestepper_module
  use Option_module
  use Init_module
  use Logging_module

  
  implicit none
#include "definitions.h"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/tao_solver.h"

  PetscScalar :: f
  PetscInt :: dummy
  PetscErrorCode :: ierr

  PetscScalar, pointer :: x_val(:)
  integer :: MAXLINELENGTH
  parameter (MAXLINELENGTH=128)
  character*(MAXLINELENGTH) :: buffer


  PetscLogDouble :: timex(4), timex_wall(4)

  PetscInt :: stage(10)
  PetscTruth :: option_found  ! For testing presence of a command-line option.
  character(len=MAXSTRINGLENGTH) :: pflotranin

  PetscMPIInt :: myrank, commsize
  
  type(simulation_type), pointer :: simulation
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option

  integer :: i,j,fd1,fd2,cur_size
  character(len=MAXWORDLENGTH) :: inputfilename
  
  double precision :: calculate_norm

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CHKERRQ(ierr)


  call LoggingCreate()

  call MPI_Comm_rank(PETSC_COMM_WORLD,myrank, ierr)
  CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,commsize,ierr)
  CHKERRQ(ierr)

  simulation => SimulationCreate()
  realization => simulation%realization
  option => realization%option

  option%myrank = myrank
  option%commsize = commsize

  

  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-pflotranin", &
                             pflotranin, option_found, ierr)

  if(.not.option_found) pflotranin = "pflotran_tao_reference.in"

  fd1=25
  fd2=26


  call PetscGetCPUTime(timex(1), ierr)
  CHKERRQ(ierr)
  call PetscGetTime(timex_wall(1), ierr)
  CHKERRQ(ierr)
  option%start_time = timex_wall(1)

  call OptionCheckCommandLine(option)

  call Init(simulation,pflotranin)

  call StepperRun(simulation%realization,simulation%flow_stepper,simulation%tran_stepper)

! Clean things up.
  print *,'Calling SimulationDestroy'
  call SimulationDestroy(simulation)
  print *,'Done Calling SimulationDestroy'

! Final Time
  call PetscGetCPUTime(timex(2), ierr)
  CHKERRQ(ierr)
  call PetscGetTime(timex_wall(2), ierr)
  CHKERRQ(ierr)
  
  if (myrank == 0) then
  
    write(*,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
      (timex(2)-timex(1))/3600.d0

    write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
      (timex_wall(2)-timex_wall(1))/3600.d0

    write(IUNIT2,'(/," CPU Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex(2)-timex(1), (timex(2)-timex(1))/60.d0, &
      (timex(2)-timex(1))/3600.d0

    write(IUNIT2,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
    & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
      timex_wall(2)-timex_wall(1), (timex_wall(2)-timex_wall(1))/60.d0, &
      (timex_wall(2)-timex_wall(1))/3600.d0
  endif

  close(IUNIT2)
  call PetscFinalize (ierr)
  CHKERRQ(ierr)

  end program generate_reference
