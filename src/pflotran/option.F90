module Option_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"


  type, public :: option_type 
  
    PetscInt :: id                         ! id of realization
  
    PetscMPIInt :: global_comm             ! MPI_COMM_WORLD
    PetscMPIInt :: global_rank             ! rank in MPI_COMM_WORLD
    PetscMPIInt :: global_commsize         ! size of MPI_COMM_WORLD
    PetscMPIInt :: global_group            ! id of group for MPI_COMM_WORLD
  
    PetscMPIInt :: mycomm                  ! PETSC_COMM_WORLD
    PetscMPIInt :: myrank                  ! rank in PETSC_COMM_WORLD
    PetscMPIInt :: mycommsize              ! size of PETSC_COMM_WORLD
    PetscMPIInt :: mygroup                 ! id of group for PETSC_COMM_WORLD
    PetscMPIInt :: mygroup_id

! don't place a character string near here.  It causes the Windows Intel compiler
! to crash.  Don't know why....
        
    PetscMPIInt :: io_rank
    PetscTruth :: broadcast_read

#ifdef VAMSI_HDF5
    MPI_Comm :: iogroup, readers 
    PetscMPIInt :: localsize, localrank, reader_rank, reader_size
    PetscInt :: color, key, reader_color, reader_key, broadcast_size
#endif	

    character(len=MAXSTRINGLENGTH) :: io_buffer
  
    PetscInt :: fid_out
    
    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXWORDLENGTH) :: flowmode
    PetscInt :: iflowmode
    character(len=MAXWORDLENGTH) :: tranmode
    PetscInt :: itranmode

    ! vector centering, used by SAMR
    ! 0 - CELL CENTERED
    ! 1 - FACE CENTERED
    PetscInt :: ivar_centering
    PetscTruth :: use_samr

    PetscInt :: nphase
    PetscInt :: liquid_phase
    PetscInt :: gas_phase
    PetscInt :: nflowdof
    PetscInt :: nflowspec

    PetscInt :: ntrandof
  
    PetscInt :: iflag
    PetscTruth :: init_stage
    PetscTruth :: print_screen_flag
    PetscTruth :: print_file_flag
    PetscTruth :: print_to_screen
    PetscTruth :: print_to_file
    PetscInt :: verbosity  ! Values >0 indicate additional console output.
    
    PetscInt, pointer :: garbage ! for some reason, Intel will not compile without this

    PetscReal :: uniform_velocity(3)
    PetscTruth :: store_solute_fluxes

    ! Program options
    PetscTruth :: use_matrix_free  ! If true, do not form the Jacobian.
    
    PetscTruth :: use_isothermal
    
    character(len=MAXWORDLENGTH) :: generalized_grid
    PetscTruth :: use_generalized_grid
      
    PetscReal :: flow_time, tran_time, time  ! The time elapsed in the simulation.
    PetscReal :: tran_weight_t0, tran_weight_t1
    PetscReal :: flow_dt, tran_dt, dt ! The size of the time step.
    PetscTruth :: match_waypoint
    PetscReal :: prev_dt
  
      ! Basically our target number of newton iterations per time step.
    PetscReal :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    PetscReal :: dpmax,dtmpmax,dsmax,dcmax

    PetscReal :: gravity(3)
    
    PetscReal :: scale
 !   PetscReal, pointer :: dencpr(:),ckwet(:)
 !   PetscReal, pointer :: sir(:,:)

    PetscReal :: m_nacl
!    PetscReal :: difaq, delhaq, eqkair, ret=1.d0, fc=1.d0
!    PetscReal :: difsc
!    PetscReal :: difgs
!    PetscReal :: disp
    
    PetscInt :: ideriv
    PetscInt :: idt_switch
    PetscReal :: reference_temperature
    PetscReal :: reference_pressure
    PetscReal :: reference_water_density
    PetscReal :: reference_porosity
    PetscReal :: reference_saturation
    
    PetscReal :: minimum_hydrostatic_pressure
    
    PetscTruth :: calculate_porosity
    PetscTruth :: initialize_with_molality
    PetscTruth :: jumpstart_kinetic_sorption
        
!   table lookup
    PetscInt :: itable
    PetscInt :: co2eos

    PetscTruth :: restart_flag
    PetscReal :: restart_time
    character(len=MAXSTRINGLENGTH) :: restart_filename
    character(len=MAXSTRINGLENGTH) :: input_filename
    PetscTruth :: checkpoint_flag
    PetscInt :: checkpoint_frequency
    
    PetscLogDouble :: start_time
    PetscTruth :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time
    
    PetscInt :: log_stage(10)
    
    PetscTruth :: numerical_derivatives
    PetscTruth :: compute_statistics
    PetscTruth :: compute_mass_balance_new
    PetscTruth :: use_touch_options
    PetscTruth :: overwrite_restart_transport
    PetscTruth :: overwrite_restart_flow
    PetscInt :: io_handshake_buffer_size
    
    character(len=MAXSTRINGLENGTH) :: permx_filename
    character(len=MAXSTRINGLENGTH) :: permy_filename
    character(len=MAXSTRINGLENGTH) :: permz_filename
    
    character(len=MAXWORDLENGTH) :: global_prefix
    character(len=MAXWORDLENGTH) :: group_prefix
    
    PetscTruth :: steady_state
    
    PetscTruth :: use_matrix_buffer

  end type option_type
  
  type, public :: output_option_type

    character(len=2) :: tunit
    PetscReal :: tconv

    PetscTruth :: print_initial
    PetscTruth :: print_final
  
    PetscTruth :: print_hdf5
    PetscTruth :: print_hdf5_velocities
    PetscTruth :: print_hdf5_flux_velocities

    PetscTruth :: print_tecplot 
    PetscInt :: tecplot_format
    PetscTruth :: print_tecplot_velocities
    PetscTruth :: print_tecplot_flux_velocities
    
    PetscTruth :: print_vtk 
    PetscTruth :: print_vtk_velocities

    PetscTruth :: print_mad 

    PetscInt :: screen_imod
    
    PetscInt :: periodic_output_ts_imod
    PetscInt :: periodic_tr_output_ts_imod
    
    PetscReal :: periodic_output_time_incr
    PetscReal :: periodic_tr_output_time_incr
    
    PetscTruth :: print_permeability
    PetscTruth :: print_porosity
    
    PetscInt :: plot_number
    character(len=MAXWORDLENGTH) :: plot_name

  end type output_option_type

  interface printMsg
    module procedure printMsg1
    module procedure printMsg2
  end interface

  interface printErrMsg
    module procedure printErrMsg1
    module procedure printErrMsg2
  end interface
  
  interface printWrnMsg
    module procedure printWrnMsg1
    module procedure printWrnMsg2
  end interface
    
  public :: OptionCreate, &
            OutputOptionCreate, &
            OptionCheckCommandLine, &
            printErrMsg, &
            printWrnMsg, &
            printMsg, &
            OptionDestroy, &
            OptionCheckTouch, &
            OptionPrintToScreen, &
            OptionPrintToFile, &
            OutputOptionDestroy, &
            OptionInitRealization

contains

! ************************************************************************** !
!
! OptionCreate: Allocates and initializes a new Option object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function OptionCreate()

  implicit none
  
  type(option_type), pointer :: OptionCreate
  
  type(option_type), pointer :: option
  
  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide 
  ! whether the member needs initialization once for all stochastic 
  ! simulations or initialization for every realization (e.g. within multiple 
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionInitAll(option)
  OptionCreate => option
  
end function OptionCreate

! ************************************************************************** !
!
! OptionInitAll: Initializes all option variables 
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OptionInitAll(option)

  implicit none
  
  type(option_type) :: option
  
  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)
  
  option%id = 0

  option%global_comm = 0
  option%global_rank = 0
  option%global_commsize = 0
  option%global_group = 0
  
  option%mycomm = 0
  option%myrank = 0
  option%mycommsize = 0
  option%mygroup = 0
  option%mygroup_id = 0
  
  option%global_prefix = 'pflotran'
  option%group_prefix = ''
    
  option%broadcast_read = PETSC_FALSE
  option%io_rank = 0

#ifdef VAMSI_HDF5
  option%iogroup = 0
  option%readers = 0
  option%localsize = 0
  option%localrank = 0
  option%reader_rank = 0
  option%reader_size = 0
  option%color = 0
  option%key = 0
  option%reader_color = 0
  option%reader_key = 0
  option%broadcast_size = 0
#endif
  
  option%print_screen_flag = PETSC_FALSE
  option%print_file_flag = PETSC_FALSE
  option%print_to_screen = PETSC_TRUE
  option%print_to_file = PETSC_TRUE
  option%verbosity = 0

  option%input_filename = ''

  call OptionInitRealization(option)

end subroutine OptionInitAll

! ************************************************************************** !
!
! OptionInitRealization: Initializes option variables specific to a single 
!                        realization
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine OptionInitRealization(option)

  implicit none
  
  type(option_type) :: option
  
  ! These variables should be initialized once at the beginning of every 
  ! PFLOTRAN realization or simulation of a single realization
    
  option%fid_out = 0

  option%iflag = 0
  option%io_buffer = ''
  
  option%use_isothermal = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  
  option%flowmode = ""
  option%iflowmode = NULL_MODE
  option%nflowdof = 0

  option%tranmode = ""
  option%itranmode = NULL_MODE
  option%ntrandof = 0
  
  option%ivar_centering = CELL_CENTERED
  option%use_samr = PETSC_FALSE

  option%nphase = 0
  option%liquid_phase = 0
  option%gas_phase = 0
  
  option%uniform_velocity = 0.d0
  option%store_solute_fluxes = PETSC_FALSE
  
!-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  option%reference_pressure = 101325.d0
  option%reference_temperature = 25.d0
  option%reference_water_density = 0.d0
  option%reference_porosity = 0.25d0
  option%reference_saturation = 1.d0
  option%calculate_porosity = PETSC_FALSE
  option%initialize_with_molality = PETSC_FALSE
  option%jumpstart_kinetic_sorption = PETSC_FALSE
  
  option%minimum_hydrostatic_pressure = -1.d20

  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  option%ideriv = 1

  option%gravity(:) = 0.d0
  option%gravity(3) = -9.8068d0    ! m/s^2

  option%dpmxe = 5.d4
  option%dtmpmxe = 5.d0
  option%dsmxe = 0.5d0
  option%dcmxe = 1.d0

  !physical constants and defult variables
!  option%difaq = 1.d-9 ! m^2/s read from input file
!  option%difaq = 0.d0
!  option%delhaq = 12.6d0 ! kJ/mol read from input file
!  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  option%m_nacl = 0.d0
  
!  option%disp = 0.d0
  
  option%generalized_grid = ""
  option%use_generalized_grid = PETSC_FALSE

  option%restart_flag = PETSC_FALSE
  option%restart_filename = ""
  option%restart_time = -999.d0
  option%checkpoint_flag = PETSC_FALSE
  option%checkpoint_frequency = huge(option%checkpoint_frequency)
  
  option%start_time = 0.d0
  option%wallclock_stop_flag = PETSC_FALSE
  option%wallclock_stop_time = 0.d0
  
  option%log_stage = 0
  
  option%numerical_derivatives = PETSC_FALSE
  option%compute_statistics = PETSC_FALSE
  option%compute_mass_balance_new = PETSC_FALSE

  option%use_touch_options = PETSC_FALSE
  option%overwrite_restart_transport = PETSC_FALSE
  option%overwrite_restart_flow = PETSC_FALSE

  option%flow_time = 0.d0
  option%tran_time = 0.d0
  option%time = 0.d0
  option%tran_weight_t0 = 0.d0
  option%tran_weight_t1 = 0.d0
  option%flow_dt = 0.d0
  option%tran_dt = 0.d0
  option%dt = 0.d0
  option%match_waypoint = PETSC_FALSE
  option%prev_dt = 0.d0

  option%io_handshake_buffer_size = 0
  
  option%permx_filename = ""
  option%permy_filename = ""
  option%permz_filename = ""
  
  option%steady_state = PETSC_FALSE
  
  option%itable=0
  option%co2eos=EOS_SPAN_WAGNER
  option%idt_switch = -1

  option%use_matrix_buffer = PETSC_FALSE
  option%init_stage = PETSC_FALSE

end subroutine OptionInitRealization

! ************************************************************************** !
!
! OutputOptionCreate: Creates output options object
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
function OutputOptionCreate()

  implicit none
  
  type(output_option_type), pointer :: OutputOptionCreate

  type(output_option_type), pointer :: output_option
  
  allocate(output_option)
  output_option%print_hdf5 = PETSC_FALSE
  output_option%print_hdf5_velocities = PETSC_FALSE
  output_option%print_hdf5_flux_velocities = PETSC_FALSE
  output_option%print_tecplot = PETSC_FALSE
  output_option%tecplot_format = 0
  output_option%print_tecplot_velocities = PETSC_FALSE
  output_option%print_tecplot_flux_velocities = PETSC_FALSE
  output_option%print_vtk = PETSC_FALSE
  output_option%print_vtk_velocities = PETSC_FALSE
  output_option%print_mad = PETSC_FALSE
  output_option%print_initial = PETSC_TRUE
  output_option%print_final = PETSC_TRUE
  output_option%plot_number = 0
  output_option%screen_imod = 1
  output_option%periodic_output_ts_imod  = 100000000
  output_option%periodic_output_time_incr = 0.d0
  output_option%periodic_tr_output_ts_imod = 100000000
  output_option%periodic_tr_output_time_incr = 0.d0
  output_option%plot_name = ""
  output_option%print_permeability = PETSC_FALSE
  output_option%print_porosity = PETSC_FALSE
  
  output_option%tconv = 1.d0
  output_option%tunit = 's'

  OutputOptionCreate => output_option
  
end function OutputOptionCreate

! ************************************************************************** !
!
! OptionCheckCommandLine: Checks all PETSc options on input
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine OptionCheckCommandLine(option)
  
  implicit none
  
  type(option_type) :: option
  
  PetscTruth :: option_found 
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-buffer_matrix", & 
                           option%use_matrix_buffer, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
                           option%use_matrix_free, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isothermal", &
                           option%use_isothermal, ierr)
                           
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', &
                             option%restart_filename, &
                             option%restart_flag, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-chkptfreq', &
                          option%checkpoint_frequency, &
                          option%checkpoint_flag, ierr)                           
  ! check on possible modes                                                     
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richards", &
                           option_found, ierr)
  if (option_found) option%flowmode = "richards"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thc", &
                           option_found, ierr)
  if (option_found) option%flowmode = "thc"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
                           option_found, ierr)
  if (option_found) option%flowmode = "mph"                           
 
 
  option_found = PETSC_FALSE
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-realization_id', &
                          temp_int,option_found, ierr)
  if (option_found) option%id = temp_int

end subroutine OptionCheckCommandLine

! ************************************************************************** !
!
! printErrMsg1: Prints the error message from p0 and stops
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printErrMsg1(option)

  implicit none
  
  type(option_type) :: option
  
  PetscTruth :: petsc_initialized
  PetscErrorCode :: ierr
  
  if (OptionPrintToScreen(option)) then
    print *
    print *, 'ERROR: ' // trim(option%io_buffer)
    print *, 'Stopping!'
  endif    
  call PetscInitialized(petsc_initialized, ierr)
  if (petsc_initialized) call PetscFinalize(ierr)
  stop
  
end subroutine printErrMsg1

! ************************************************************************** !
!
! printErrMsg2: Prints the error message from p0 and stops
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printErrMsg2(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  PetscTruth :: petsc_initialized
  PetscErrorCode :: ierr
  
  if (OptionPrintToScreen(option)) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *, 'Stopping!'
  endif    
  call PetscInitialized(petsc_initialized, ierr)
  if (petsc_initialized) call PetscFinalize(ierr)
  stop
  
end subroutine printErrMsg2
 
! ************************************************************************** !
!
! printWrnMsg1: Prints the warning message from p0
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printWrnMsg1(option)

  implicit none
  
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) print *, 'WARNING: ' // trim(option%io_buffer)
  
end subroutine printWrnMsg1

! ************************************************************************** !
!
! printWrnMsg2: Prints the warning message from p0
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printWrnMsg2(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  if (OptionPrintToScreen(option)) print *, 'WARNING: ' // trim(string)
  
end subroutine printWrnMsg2

! ************************************************************************** !
!
! printMsg1: Prints the message from p0
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine printMsg1(option)

  implicit none
  
  type(option_type) :: option
  
  if (OptionPrintToScreen(option)) print *, trim(option%io_buffer)
  
end subroutine printMsg1

! ************************************************************************** !
!
! printMsg2: Prints the message from p0
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine printMsg2(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  if (OptionPrintToScreen(option)) print *, trim(string)
  
end subroutine printMsg2

! ************************************************************************** !
!
! OptionCheckTouch: Users can steer the code by touching files.
! author: Glenn Hammond
! date: 03/04/08
!
! ************************************************************************** !
function OptionCheckTouch(option,filename)

  implicit none

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: ios
  PetscInt :: fid = 86
  PetscTruth :: OptionCheckTouch
  PetscErrorCode :: ierr
  
  OptionCheckTouch = PETSC_FALSE

  if (option%myrank == option%io_rank) &
    open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,1,MPI_INTEGER,option%io_rank,option%mycomm,ierr)

  if (ios == 0) then
    if (option%myrank == option%io_rank) close(fid,status='delete')
    OptionCheckTouch = PETSC_TRUE
  endif

end function OptionCheckTouch

! ************************************************************************** !
!
! OptionPrintToScreen: Determines whether printing should occur
! author: Glenn Hammond
! date: 12/09/08
!
! ************************************************************************** !
function OptionPrintToScreen(option)

  implicit none

  type(option_type) :: option
  
  PetscTruth :: OptionPrintToScreen
  
  if (option%myrank == option%io_rank .and. option%print_to_screen) then
    OptionPrintToScreen = PETSC_TRUE
  else
    OptionPrintToScreen = PETSC_FALSE
  endif

end function OptionPrintToScreen

! ************************************************************************** !
!
! OptionPrintToFile: Determines whether printing to file should occur
! author: Glenn Hammond
! date: 01/29/09
!
! ************************************************************************** !
function OptionPrintToFile(option)

  implicit none

  type(option_type) :: option
  
  PetscTruth :: OptionPrintToFile
  
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    OptionPrintToFile = PETSC_TRUE
  else
    OptionPrintToFile = PETSC_FALSE
  endif

end function OptionPrintToFile

! ************************************************************************** !
!
! OutputOptionDestroy: Deallocates an output option
! author: Glenn Hammond
! date: 11/07/07
!
! ************************************************************************** !
subroutine OutputOptionDestroy(output_option)

  implicit none
  
  type(output_option_type), pointer :: output_option
  
  deallocate(output_option)
  nullify(output_option)
  
end subroutine OutputOptionDestroy

! ************************************************************************** !
!
! OptionDestroy: Deallocates an option
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine OptionDestroy(option)

  implicit none
  
  type(option_type), pointer :: option
  
  ! all kinds of stuff needs to be added here.

  ! all the below should be placed somewhere other than option.F90
#if 0
  if (associated(option%dencpr)) deallocate(option%dencpr)
  nullify(option%dencpr)
  if (associated(option%ckwet)) deallocate(option%ckwet)
  nullify(option%ckwet)
  if (associated(option%sir)) deallocate(option%sir)
  nullify(option%sir)
  if (associated(option%tfac)) deallocate(option%tfac)
  nullify(option%tfac)
#endif  
  
  deallocate(option)
  nullify(option)
  
end subroutine OptionDestroy

end module Option_module
