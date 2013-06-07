module PFLOTRAN_Factory_module

  implicit none

  private

#include "definitions.h"

  public :: InitializePFLOTRAN, &
            FinalizePFLOTRAN

contains

! ************************************************************************** !
!
! InitializePFLOTRAN: Sets up PFLOTRAN subsurface simulation framework
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine InitializePFLOTRAN(option)

  use Option_module
  use Logging_module
  use Init_module
  
  implicit none
  
  type(option_type), pointer :: option
  
  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr
  
  option => OptionCreate()
  
  call InitMPI(option)
  call InitPFLOTRANCommandLineSettings(option)
  call LoggingCreate()
  call PetscTime(timex_wall, ierr)
  option%start_time = timex_wall
  
end subroutine InitializePFLOTRAN

! ************************************************************************** !
!
! FinalizePFLOTRAN: Sets up PFLOTRAN subsurface simulation framework
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine FinalizePFLOTRAN(simulation)

  use Simulation_module
  use Option_module
  use Regression_module
  use Logging_module
  
  implicit none
  
  type(simulation_type), pointer :: simulation

  type(option_type), pointer :: option
  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr
  
  option => simulation%realization%option
  
  
  call RegressionOutput(simulation%regression,simulation%realization, &
                        simulation%flow_stepper,simulation%tran_stepper)

! Clean things up.
  call SimulationDestroy(simulation)

! Final Time
  call PetscTime(timex_wall, ierr)
    
  if (option%myrank == option%io_rank) then

    if (option%print_to_screen) then
      write(*,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall-option%start_time, &
        (timex_wall-option%start_time)/60.d0, &
        (timex_wall-option%start_time)/3600.d0
    endif
    if (option%print_to_file) then
      write(option%fid_out,'(/," Wall Clock Time:", 1pe12.4, " [sec] ", &
      & 1pe12.4, " [min] ", 1pe12.4, " [hr]")') &
        timex_wall-option%start_time, &
        (timex_wall-option%start_time)/60.d0, &
        (timex_wall-option%start_time)/3600.d0
    endif
  endif

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif

  call LoggingDestroy()

  call PetscOptionsSetValue('-options_left','no',ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue('-objects_left','yes',ierr)

  call OptionDestroy(option)
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)
  
end subroutine FinalizePFLOTRAN

end module PFLOTRAN_Factory_module
