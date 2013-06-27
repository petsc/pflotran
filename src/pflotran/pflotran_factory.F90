module PFLOTRAN_Factory_module

  implicit none

  private

#include "definitions.h"

  public :: PFLOTRANInitialize, &
            PFLOTRANInitializePrePETSc, &
#ifndef PROCESS_MODEL
            PFLOTRANInitializePostPETSc, &
            PFLOTRANRun, &
#endif
            PFLOTRANFinalize

contains

! ************************************************************************** !
!
! PFLOTRANInitialize: Sets up PFLOTRAN subsurface simulation 
! author: Glenn Hammond
! date: 06/10/13
!
! ************************************************************************** !
subroutine PFLOTRANInitialize(option)

  use Option_module
  use Input_module
  use Timestepper_module
  
  implicit none
  
  type(option_type), pointer :: option

  call PFLOTRANInitializePrePETSc(option)
  call OptionInitPetsc(option)
#ifdef PROCESS_MODEL
  call PFLOTRANInitializePostPETSc(option)
#endif

end subroutine PFLOTRANInitialize

! ************************************************************************** !
!
! PFLOTRANInitializePrePETSc: Sets up PFLOTRAN subsurface simulation 
!                             framework prior to PETSc initialization
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANInitializePrePETSc(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  
  ! NOTE: Cannot add anything that requires PETSc in this routins as PETSc 
  !       has not yet been initialized.
  
  call PFLOTRANInitCommandLineSettings(option)
  
end subroutine PFLOTRANInitializePrePETSc

! ************************************************************************** !
!
! PFLOTRANInitializePostPETSc: Sets up PFLOTRAN subsurface simulation 
!                              framework after PETSc initialization
! author: Glenn Hammond
! date: 06/17/13
!
! ************************************************************************** !
#ifdef PROCESS_MODEL
subroutine PFLOTRANInitializePostPETSc(option)

  use Option_module
  use Logging_module
  
  implicit none
  
  type(option_type) :: option

  call LoggingCreate()
  call OptionBeginTiming(option)
  
#else
subroutine PFLOTRANInitializePostPETSc(simulation, master_stepper, option, &
                                       init_status)

  use Option_module
  use Init_module

  implicit none

  type(simulation_type), pointer :: simulation
  type(stepper_type), pointer :: master_stepper
  type(option_type), pointer :: option
  PetscInt :: init_status

  call OptionBeginTiming(option)
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
#else
  call TimestepperInitializeRun(simulation%realization, &
                                master_stepper, &
                                simulation%flow_stepper, &
                                simulation%tran_stepper, &
                                init_status)
#endif 
#endif 
end subroutine PFLOTRANInitializePostPETSc

#ifndef PROCESS_MODEL
! ************************************************************************** !
!
! PFLOTRANRun: Runs the PFLOTRAN simulation
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine PFLOTRANRun(simulation, master_stepper, init_status)

  use Simulation_module
  use Timestepper_module

  implicit none

  type(simulation_type) :: simulation
  type(stepper_type), pointer :: master_stepper
  PetscInt :: init_status

  select case(init_status)
    case(TIMESTEPPER_INIT_PROCEED)
#ifdef SURFACE_FLOW
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
#else
      call  TimestepperExecuteRun(simulation%realization, &
                                  master_stepper, &
                                  simulation%flow_stepper, &
                                  simulation%tran_stepper)
      call  TimestepperFinalizeRun(simulation%realization, &
                                    master_stepper, &
                                    simulation%flow_stepper, &
                                    simulation%tran_stepper)
#endif 
    case(TIMESTEPPER_INIT_FAIL)
    case(TIMESTEPPER_INIT_DONE)
  end select

end subroutine PFLOTRANRun
#endif

! ************************************************************************** !
!
! PFLOTRANFinalize: Destroys PFLOTRAN subsurface simulation framework
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
#ifdef PROCESS_MODEL
subroutine PFLOTRANFinalize(option)

  use Option_module
  use Logging_module
  
  implicit none
  
  type(option_type) :: option
  
  call OptionEndTiming(option)
  call LoggingDestroy()
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif
#else
subroutine PFLOTRANFinalize(simulation,option)

  use Option_module
  implicit none

  type(simulation_type), pointer :: simulation
  type(option_type) :: option

  call RegressionOutput(simulation%regression,simulation%realization, &
                        simulation%flow_stepper,simulation%tran_stepper)

! Clean things up.
  call SimulationDestroy(simulation)
  call OptionEndTiming(option)
  call LoggingDestroy()
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    close(option%fid_out)
  endif
#endif

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
  use String_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string, string2
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
 
  ! this will get overwritten later if stochastic
  string = '-realization_id'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%id = i
  
  ! this will get overwritten later if stochastic
  string = '-simulation_mode'
  call InputGetCommandLineString(string,string2, &
                                 option_found,option)
  if (option_found) then
    call StringToUpper(string2)
    option%simulation_mode = string2
  endif

end subroutine PFLOTRANInitCommandLineSettings

end module PFLOTRAN_Factory_module
