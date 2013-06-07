module PFLOTRAN_Factory_module

  implicit none

  private

#include "definitions.h"

  public :: PFLOTRANInitialize, &
            PFLOTRANRun, &
            PFLOTRANFinalize

contains

! ************************************************************************** !
!
! PFLOTRANInitialize: Sets up PFLOTRAN subsurface simulation framework prior
!                     to PETSc initialization
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANInitialize(option)

  use Option_module
  use Logging_module
  
  implicit none
  
  type(option_type), pointer :: option
  
  option => OptionCreate()
  
  ! NOTE: Cannot add anything that requires PETSc in this routins as PETSc 
  !       has not yet been initialized.
  
  call OptionInitMPI(option)
  call PFLOTRANInitCommandLineSettings(option)
  
end subroutine PFLOTRANInitialize

! ************************************************************************** !
!
! PFLOTRANRun: Runs the PFLOTRAN simulation
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANRun(option)

  use Simulation_module
  use Realization_class
  use Timestepper_module
  use Option_module
  use Simulation_module
  use Regression_module
  use Init_module
  
  implicit none
  
  type(option_type), pointer :: option

  type(simulation_type), pointer :: simulation
  type(stepper_type), pointer :: master_stepper

  PetscInt :: init_status
  
  simulation => SimulationCreate(option)
  call Init(simulation)

#ifdef SURFACE_FLOW
  call TimestepperInitializeRun(simulation%realization, &
                                simulation%surf_realization, &
                                master_stepper, &
                                simulation%flow_stepper, &
                                simulation%tran_stepper, &
                                simulation%surf_flow_stepper, &
                                init_status)
  select case(init_status)
    case(TIMESTEPPER_INIT_PROCEED)
      call  TimestepperExecuteRun(simulation%realization, &
                                  simulation%surf_realization, &
                                  master_stepper, &
                                  simulation%flow_stepper, &
                                  simulation%tran_stepper, &
                                  simulation%surf_flow_stepper)
      call  TimestepperFinalizeRun(simulation%realization, &
                                    simulation%surf_realization, &
                                    master_stepper, &
                                    simulation%flow_stepper, &
                                    simulation%tran_stepper, &
                                    simulation%surf_flow_stepper)
    case(TIMESTEPPER_INIT_FAIL)
    case(TIMESTEPPER_INIT_DONE)
  end select
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

! Clean things up.
  call SimulationDestroy(simulation)
  
end subroutine PFLOTRANRun

! ************************************************************************** !
!
! PFLOTRANFinalize: Destroys PFLOTRAN subsurface simulation framework
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANFinalize(option)

  use Option_module
  use Logging_module
  
  implicit none
  
  type(option_type), pointer :: option
  
  call OptionEndLogging(option)
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif
  call OptionFinalize(option)

end subroutine PFLOTRANFinalize

! ************************************************************************** !
!
! PFLOTRANInitCommandLineSettings: Initializes PFLOTRAN output filenames, etc.
! author: Glenn Hammond
! date: 06/06/13
!
! ************************************************************************** !
subroutine PFLOTRANInitCommandLineSettings(option)

  use Option_module
  use Input_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  PetscInt :: i
  PetscErrorCode :: ierr
  
  ! check for non-default input filename
  option%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename, &
                                 pflotranin_option_found,option)
  string = '-input_prefix'
  call InputGetCommandLineString(string,option%input_prefix, &
                                 input_prefix_option_found,option)
  
  if (pflotranin_option_found .and. input_prefix_option_found) then
    option%io_buffer = 'Cannot specify both "-pflotranin" and ' // &
      '"-input_prefix" on the command lines.'
    call printErrMsg(option)
  else if (pflotranin_option_found) then
    !TODO(geh): replace this with StringSplit()
    i = index(option%input_filename,'.',PETSC_TRUE)
    if (i > 1) then
      i = i-1
    else
      ! for some reason len_trim doesn't work on MS Visual Studio in 
      ! this location
      i = len(trim(option%input_filename)) 
    endif
    option%input_prefix = option%input_filename(1:i)
  else if (input_prefix_option_found) then
    option%input_filename = trim(option%input_prefix) // '.in'
  endif
  
  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)
  if (.not.option_found) option%global_prefix = option%input_prefix  
  
  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-v'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%verbosity = 1
 
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%simulation_type = MULTISIMULATION_SIM_TYPE

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) option%simulation_type = STOCHASTIC_SIM_TYPE

  ! this will get overwritten later if stochastic
  string = '-realization_id'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%id = i

end subroutine PFLOTRANInitCommandLineSettings

end module PFLOTRAN_Factory_module
