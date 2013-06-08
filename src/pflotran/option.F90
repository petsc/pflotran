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
    PetscMPIInt :: hdf5_read_group_size, hdf5_write_group_size
    PetscBool :: broadcast_read
    
    PetscInt :: reactive_transport_coupling

#if defined(SCORPIO)
    PetscMPIInt :: ioread_group_id, iowrite_group_id
#endif

    character(len=MAXSTRINGLENGTH) :: io_buffer
  
    PetscInt :: fid_out
    
    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXWORDLENGTH) :: flowmode
    PetscInt :: iflowmode
    character(len=MAXWORDLENGTH) :: tranmode
    PetscInt :: itranmode
    PetscInt :: tvd_flux_limiter

    PetscInt :: nphase
    PetscInt :: liquid_phase
    PetscInt :: gas_phase
    PetscInt :: nflowdof
    PetscInt :: nflowspec
    PetscInt :: rt_idof
    PetscInt :: nmechdof
    PetscInt :: nsec_cells
#ifdef SURFACE_FLOW
    PetscInt :: nsurfflowdof
    PetscInt :: subsurf_surf_coupling
    PetscInt :: surface_flow_formulation
    PetscReal :: surf_flow_time, surf_flow_dt
    PetscReal :: surf_subsurf_coupling_time
    PetscReal :: surf_subsurf_coupling_flow_dt
    PetscBool :: surf_flow_explicit
    character(len=MAXSTRINGLENGTH) :: surf_initialize_flow_filename
#endif
    PetscBool :: sec_vars_update
    PetscInt :: air_pressure_id
    PetscInt :: capillary_pressure_id
    PetscInt :: vapor_pressure_id 
    PetscInt :: water_id  ! index of water component dof
    PetscInt :: air_id  ! index of air component dof
    PetscInt :: energy_id  ! index of energy dof

    PetscInt :: ntrandof
  
    PetscBool :: variables_swapped
    
    PetscInt :: iflag
    PetscBool :: init_stage
    PetscBool :: print_screen_flag
    PetscBool :: print_file_flag
    PetscBool :: print_to_screen
    PetscBool :: print_to_file
    PetscInt :: verbosity  ! Values >0 indicate additional console output.
    
    PetscInt, pointer :: garbage ! for some reason, Intel will not compile without this

    PetscReal :: uniform_velocity(3)
    PetscBool :: store_solute_fluxes
    PetscBool :: store_flowrate

    ! Program options
    PetscBool :: use_matrix_free  ! If true, do not form the Jacobian.
    
    PetscBool :: use_isothermal
    PetscBool :: use_mc           ! If true, multiple continuum formulation is used.
    PetscBool :: set_secondary_init_temp  ! If true, then secondary init temp is different from prim. init temp
    PetscBool :: set_secondary_init_conc
    
    PetscBool :: update_flow_perm ! If true, permeability changes due to pressure    
    
    character(len=MAXWORDLENGTH) :: generalized_grid
    PetscBool :: use_generalized_grid
      
    PetscReal :: flow_time, tran_time, time  ! The time elapsed in the simulation.
    PetscReal :: tran_weight_t0, tran_weight_t1
    PetscReal :: flow_dt, tran_dt ! The size of the time step.
    PetscBool :: match_waypoint
  
      ! Basically our target number of newton iterations per time step.
    PetscReal :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    PetscReal :: dpmax,dtmpmax,dsmax,dcmax

    PetscReal :: gravity(3)
    
    PetscReal :: scale

    PetscReal :: m_nacl
    
    PetscInt :: ideriv
    PetscInt :: idt_switch
    PetscReal :: reference_temperature
    PetscReal :: reference_pressure
    PetscReal :: reference_water_density
    PetscReal :: reference_porosity
    PetscReal :: reference_saturation
    
    PetscReal :: pressure_dampening_factor
    PetscReal :: saturation_change_limit
    PetscReal :: pressure_change_limit
    PetscReal :: temperature_change_limit
    PetscReal :: stomp_norm
    PetscBool :: check_stomp_norm
    
    PetscReal :: infnorm_res_sec  ! inf. norm of secondary continuum rt residual
    
    PetscReal :: minimum_hydrostatic_pressure
    
    PetscBool :: jumpstart_kinetic_sorption
    PetscBool :: no_checkpoint_kinetic_sorption
    PetscBool :: no_restart_kinetic_sorption
    PetscBool :: no_restart_mineral_vol_frac
        
!   table lookup
    PetscInt :: itable
    PetscInt :: co2eos
    character(len=MAXSTRINGLENGTH) :: co2_database_filename

    PetscBool :: restart_flag
    PetscReal :: restart_time
    character(len=MAXSTRINGLENGTH) :: restart_filename
    character(len=MAXSTRINGLENGTH) :: input_filename
    PetscBool :: checkpoint_flag
    PetscInt :: checkpoint_frequency
    
    PetscLogDouble :: start_time
    PetscBool :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time
    
    PetscInt :: log_stage(10)
    
    PetscBool :: numerical_derivatives_flow
    PetscBool :: numerical_derivatives_rxn
    PetscBool :: numerical_derivatives_multi_coupling
    PetscBool :: compute_statistics
    PetscBool :: compute_mass_balance_new
    PetscBool :: use_touch_options
    PetscBool :: overwrite_restart_transport
    PetscBool :: overwrite_restart_flow
    PetscInt :: io_handshake_buffer_size

    character(len=MAXSTRINGLENGTH) :: initialize_flow_filename
    character(len=MAXSTRINGLENGTH) :: initialize_transport_filename
        
    character(len=MAXSTRINGLENGTH) :: input_prefix
    character(len=MAXSTRINGLENGTH) :: global_prefix
    character(len=MAXWORDLENGTH) :: group_prefix
    
    PetscBool :: steady_state
    PetscBool :: use_matrix_buffer
    PetscBool :: force_newton_iteration
    PetscBool :: mimetic
    PetscBool :: ani_relative_permeability
    PetscBool :: use_upwinding
    PetscBool :: out_of_table
    PetscBool :: print_explicit_primal_grid    ! prints primal grid if true
    PetscBool :: print_explicit_dual_grid      ! prints voronoi (dual) grid if true
    PetscInt :: secondary_continuum_solver     ! Specify secondary continuum solver
    
    PetscInt :: simulation_type

  end type option_type
  
  PetscInt, parameter, public :: SUBSURFACE_SIM_TYPE = 1
  PetscInt, parameter, public :: MULTISIMULATION_SIM_TYPE = 2
  PetscInt, parameter, public :: STOCHASTIC_SIM_TYPE = 3
  
  interface printMsg
    module procedure printMsg1
    module procedure printMsg2
  end interface

  interface printMsgAnyRank
    module procedure printMsgAnyRank1
    module procedure printMsgAnyRank2
  end interface

  interface printMsgByRank
    module procedure printMsgByRank1
    module procedure printMsgByRank2
  end interface

  interface printErrMsgByRank
    module procedure printErrMsgByRank1
    module procedure printErrMsgByRank2
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
            OptionCheckCommandLine, &
            printErrMsg, &
            printErrMsgByRank, &
            printWrnMsg, &
            printMsg, &
            printMsgAnyRank, &
            printMsgByRank, &
            OptionCheckTouch, &
            OptionPrintToScreen, &
            OptionPrintToFile, &
            OptionInitRealization, &
            OptionMeanVariance, &
            OptionMaxMinMeanVariance, &
            OptionInitMPI, &
            OptionInitPETSc, &
            OptionDivvyUpSimulations, &
            OptionCreateProcessorGroups, &
            OptionBeginTiming, &
            OptionEndTiming, &
            OptionFinalize, &
            OptionDestroy

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
  
  option%input_prefix = 'pflotran'
  option%group_prefix = ''
  option%global_prefix = ''
    
  option%broadcast_read = PETSC_FALSE
  option%io_rank = 0
  option%hdf5_read_group_size = 0
  option%hdf5_write_group_size = 0

  option%print_screen_flag = PETSC_FALSE
  option%print_file_flag = PETSC_FALSE
  option%print_to_screen = PETSC_TRUE
  option%print_to_file = PETSC_TRUE
  option%verbosity = 0

  option%input_filename = ''

  option%mimetic = PETSC_FALSE
  option%ani_relative_permeability = PETSC_FALSE

  option%use_upwinding = PETSC_TRUE

  option%out_of_table = PETSC_FALSE
  
  option%simulation_type = SUBSURFACE_SIM_TYPE
 
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
    
  option%fid_out = OUT_UNIT

  option%iflag = 0
  option%io_buffer = ''
  
  option%use_isothermal = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  option%use_mc = PETSC_FALSE
  option%set_secondary_init_temp = PETSC_FALSE
  option%set_secondary_init_conc = PETSC_FALSE
  
  option%update_flow_perm = PETSC_FALSE
  
  option%flowmode = ""
  option%iflowmode = NULL_MODE
  option%nflowdof = 0
  option%nmechdof = 0
  option%nsec_cells = 0
#ifdef SURFACE_FLOW
   option%nsurfflowdof = 0
   option%subsurf_surf_coupling = DECOUPLED
   option%surface_flow_formulation = KINEMATIC_WAVE
   option%surf_flow_dt = 0.d0
   option%surf_flow_time =0.d0
   option%surf_subsurf_coupling_time = 0.d0
   option%surf_subsurf_coupling_flow_dt = 0.d0
   option%surf_flow_explicit = PETSC_TRUE
   option%surf_initialize_flow_filename =""
#endif

  option%tranmode = ""
  option%itranmode = NULL_MODE
  option%ntrandof = 0
  option%tvd_flux_limiter = 1
  option%rt_idof = -999
  
  option%reactive_transport_coupling = GLOBAL_IMPLICIT

  option%nphase = 0
  option%liquid_phase = 0
  option%gas_phase = 0
  
  option%air_pressure_id = 0
  option%capillary_pressure_id = 0
  option%vapor_pressure_id = 0

  option%water_id = 0
  option%air_id = 0
  option%energy_id = 0
  
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
  
  option%pressure_dampening_factor = 0.d0
  option%saturation_change_limit = 0.d0
  option%pressure_change_limit = 0.d0
  option%temperature_change_limit = 0.d0
  option%stomp_norm = 0.d0
  option%check_stomp_norm = PETSC_FALSE
  
  option%infnorm_res_sec = 0.d0
  
  option%jumpstart_kinetic_sorption = PETSC_FALSE
  option%no_checkpoint_kinetic_sorption = PETSC_FALSE
  option%no_restart_kinetic_sorption = PETSC_FALSE
  option%no_restart_mineral_vol_frac = PETSC_FALSE
  
  option%minimum_hydrostatic_pressure = -1.d20

  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  option%ideriv = 1

  option%gravity(:) = 0.d0
  option%gravity(3) = -9.8068d0    ! m/s^2

  option%dpmxe = 5.d5
  option%dtmpmxe = 5.d0
  option%dsmxe = 0.5d0
  option%dcmxe = 1.d0

  option%dpmax = 0.d0
  option%dtmpmax = 0.d0
  option%dsmax = 0.d0
  option%dcmax = 0.d0

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
  
  option%numerical_derivatives_flow = PETSC_FALSE
  option%numerical_derivatives_rxn = PETSC_FALSE
  option%numerical_derivatives_multi_coupling = PETSC_FALSE
  option%compute_statistics = PETSC_FALSE
  option%compute_mass_balance_new = PETSC_FALSE
  option%store_flowrate = PETSC_FALSE
#ifdef STORE_FLOWRATES
  option%store_flowrate = PETSC_TRUE
#endif

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
  option%match_waypoint = PETSC_FALSE

  option%io_handshake_buffer_size = 0

  option%initialize_flow_filename = ''
  option%initialize_transport_filename = ''
  
  option%steady_state = PETSC_FALSE
  
  option%itable = 0
  option%co2eos = EOS_SPAN_WAGNER
  option%co2_database_filename = ''

! option%idt_switch = 1
  option%idt_switch = -1

  option%use_matrix_buffer = PETSC_FALSE
  option%init_stage = PETSC_FALSE 
  option%force_newton_iteration = PETSC_FALSE
  option%mimetic = PETSC_FALSE
  option%variables_swapped = PETSC_FALSE
  option%print_explicit_primal_grid = PETSC_FALSE
  option%print_explicit_dual_grid = PETSC_FALSE  
  option%secondary_continuum_solver = 1
  
end subroutine OptionInitRealization

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
  
  PetscBool :: option_found 
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-buffer_matrix", & 
                           option%use_matrix_buffer, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
                           option%use_matrix_free, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isothermal", &
                           option%use_isothermal, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mc", &
                           option%use_mc, ierr)
                           
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
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thmc", &
                           option_found, ierr)
  if (option_found) option%flowmode = "thmc"     
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
                           option_found, ierr)
  if (option_found) option%flowmode = "mph"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_flash2", &
                           option_found, ierr)
  if (option_found) option%flowmode = "flash2"                           
 
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
  
  call printErrMsg2(option,option%io_buffer)
  
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
  
  PetscBool :: petsc_initialized
  PetscErrorCode :: ierr
  
  if (OptionPrintToScreen(option)) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *, 'Stopping!'
  endif    
  call MPI_Barrier(option%mycomm,ierr)
  call PetscInitialized(petsc_initialized, ierr)
  if (petsc_initialized) call PetscFinalize(ierr)
  stop
  
end subroutine printErrMsg2
 
! ************************************************************************** !
!
! printErrMsgByRank1: Prints the error message from processor with error along
!               with rank
! author: Glenn Hammond
! date: 11/04/11
!
! ************************************************************************** !
subroutine printErrMsgByRank1(option)

  implicit none
  
  type(option_type) :: option
  
  call printErrMsgByRank2(option,option%io_buffer)
  
end subroutine printErrMsgByRank1

! ************************************************************************** !
!
! printErrMsgByRank2: Prints the error message from processor with error along
!               with rank
! author: Glenn Hammond
! date: 11/04/11
!
! ************************************************************************** !
subroutine printErrMsgByRank2(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  character(len=MAXWORDLENGTH) :: word
  
  write(word,*) option%myrank
  print *
  print *, 'ERROR(' // trim(adjustl(word)) // '): ' // trim(option%io_buffer)
  print *, 'Stopping!'
  stop
  
end subroutine printErrMsgByRank2
 
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
  
  call printWrnMsg2(option,option%io_buffer)
  
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
  
  call printMsg2(option,option%io_buffer)
  
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
! printMsgAnyRank1: Prints the message from any processor core
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine printMsgAnyRank1(option)

  implicit none
  
  type(option_type) :: option
  
  call printMsgAnyRank2(option%io_buffer)
  
end subroutine printMsgAnyRank1

! ************************************************************************** !
!
! printMsgAnyRank2: Prints the message from any processor core
! author: Glenn Hammond
! date: 01/12/12
!
! ************************************************************************** !
subroutine printMsgAnyRank2(string)

  implicit none
  
  character(len=*) :: string
  
  print *, trim(string)
  
end subroutine printMsgAnyRank2

! ************************************************************************** !
!
! printMsgByRank1: Prints a message from processor along with rank
! author: Glenn Hammond
! date: 03/27/12
!
! ************************************************************************** !
subroutine printMsgByRank1(option)

  implicit none
  
  type(option_type) :: option
  
  call printMsgByRank2(option,option%io_buffer)
  
end subroutine printMsgByRank1

! ************************************************************************** !
!
! printMsgByRank2: Prints a message from processor along with rank
! author: Glenn Hammond
! date: 03/27/12
!
! ************************************************************************** !
subroutine printMsgByRank2(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  character(len=MAXWORDLENGTH) :: word
  
  write(word,*) option%myrank
  print *, '(' // trim(adjustl(word)) // '): ' // trim(option%io_buffer)
  
end subroutine printMsgByRank2
 
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
  PetscBool :: OptionCheckTouch
  PetscErrorCode :: ierr
  
  OptionCheckTouch = PETSC_FALSE

  if (option%myrank == option%io_rank) &
    open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                 option%mycomm,ierr)

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
  
  PetscBool :: OptionPrintToScreen
  
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
  
  PetscBool :: OptionPrintToFile
  
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    OptionPrintToFile = PETSC_TRUE
  else
    OptionPrintToFile = PETSC_FALSE
  endif

end function OptionPrintToFile

! ************************************************************************** !
!
! OptionMaxMinMeanVariance: Calculates the maximum, minumum, mean and 
!                           optionally variance of a number across processor 
!                           cores
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine OptionMaxMinMeanVariance(value,max,min,mean,variance, &
                                    calculate_variance,option)

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: max
  PetscReal :: min
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real_in(2), temp_real_out(2)
  PetscErrorCode :: ierr
  
  temp_real_in(1) = value
  temp_real_in(2) = -1.d0*value
  call MPI_Allreduce(temp_real_in,temp_real_out,TWO_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MAX,option%mycomm,ierr)
  max = temp_real_out(1)
  min = -1.d0*temp_real_out(2)
  
  call OptionMeanVariance(value,mean,variance,calculate_variance,option)
  
end subroutine OptionMaxMinMeanVariance

! ************************************************************************** !
!
! OptionMeanVariance: Calculates the mean and optionally variance of a number
!                     across processor cores
! author: Glenn Hammond
! date: 05/29/10
!
! ************************************************************************** !
subroutine OptionMeanVariance(value,mean,variance,calculate_variance,option)

  implicit none

  type(option_type) :: option
  PetscReal :: value
  PetscReal :: mean
  PetscReal :: variance
  PetscBool :: calculate_variance

  PetscReal :: temp_real
  PetscErrorCode :: ierr
  
  call MPI_Allreduce(value,temp_real,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr)
  mean = temp_real / dble(option%mycommsize)
  
  if (calculate_variance) then
    temp_real = value-mean
    temp_real = temp_real*temp_real
    call MPI_Allreduce(temp_real,variance,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_SUM,option%mycomm,ierr)
    variance = variance / dble(option%mycommsize)
  endif
  
end subroutine OptionMeanVariance

! ************************************************************************** !
!
! OptionInitMPI: Initializes base MPI communicator
! author: Glenn Hammond
! date: 06/06/13
!
! ************************************************************************** !
subroutine OptionInitMPI(option)

  implicit none
  
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,option%global_commsize,ierr)
  call MPI_Comm_group(MPI_COMM_WORLD,option%global_group,ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

end subroutine OptionInitMPI

! ************************************************************************** !
!
! OptionInitPETSc: Initialization of PETSc.
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine OptionInitPETSc(option)

  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  PETSC_COMM_WORLD = option%mycomm
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  
  if (option%verbosity > 0) then 
    call PetscLogBegin(ierr)
    string = '-log_summary'
    call PetscOptionsInsertString(string, ierr)
  endif 
  
end subroutine OptionInitPETSc

! ************************************************************************** !
!
! OptionBeginTiming: Start outer timing.
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine OptionBeginTiming(option)

  use Logging_module
  
  implicit none
  
#include "finclude/petsclog.h"
  
  type(option_type) :: option
  
  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr
  
  call PetscTime(timex_wall, ierr)
  option%start_time = timex_wall
  
end subroutine OptionBeginTiming

! ************************************************************************** !
!
! OptionEndTiming: End timing.
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine OptionEndTiming(option)

  use Logging_module
  
  implicit none
  
#include "finclude/petsclog.h"
  
  type(option_type) :: option
  
  PetscLogDouble :: timex_wall
  PetscErrorCode :: ierr
  
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

end subroutine OptionEndTiming

! ************************************************************************** !
!
! OptionDivvyUpSimulations: Divides simulation in to multple simulations with
!                           multiple input decks
! author: Glenn Hammond
! date: 06/06/13
!
! ************************************************************************** !
subroutine OptionDivvyUpSimulations(option,filenames)

  implicit none
  
  type(option_type) :: option
  
  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)
  
  i = size(filenames) 
  call OptionCreateProcessorGroups(option,i)
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
  
end subroutine OptionDivvyUpSimulations

! ************************************************************************** !
!
! OptionCreateProcessorGroups: Splits MPI_COMM_WORLD into N separate
!                              processor groups
! author: Glenn Hammond
! date: 08/11/09
!
! ************************************************************************** !
subroutine OptionCreateProcessorGroups(option,num_groups)

  implicit none

  type(option_type) :: option
  PetscInt :: num_groups

  PetscInt :: local_commsize
  PetscInt :: offset, delta, remainder
  PetscInt :: igroup
  PetscMPIInt :: mycolor_mpi, mykey_mpi
  PetscErrorCode :: ierr

  local_commsize = option%global_commsize / num_groups
  remainder = option%global_commsize - num_groups * local_commsize
  offset = 0
  do igroup = 1, num_groups
    delta = local_commsize
    if (igroup < remainder) delta = delta + 1
    if (option%global_rank >= offset .and. &
        option%global_rank < offset + delta) exit
    offset = offset + delta
  enddo
  mycolor_mpi = igroup
  option%mygroup_id = igroup
  mykey_mpi = option%global_rank - offset
  call MPI_Comm_split(MPI_COMM_WORLD,mycolor_mpi,mykey_mpi,option%mycomm,ierr)
  call MPI_Comm_group(option%mycomm,option%mygroup,ierr)

  call MPI_Comm_rank(option%mycomm,option%myrank, ierr)
  call MPI_Comm_size(option%mycomm,option%mycommsize,ierr)

end subroutine OptionCreateProcessorGroups

! ************************************************************************** !
!
! OptionFinalize: End the simulation.
! author: Glenn Hammond
! date: 06/07/13
!
! ************************************************************************** !
subroutine OptionFinalize(option)

  use Logging_module

  implicit none
  
  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr
  
  call PetscOptionsSetValue('-options_left','no',ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue('-objects_left','yes',ierr)
  call OptionDestroy(option)
  call PetscFinalize(ierr)
  call MPI_Finalize(ierr)
  call exit(86)
  
end subroutine OptionFinalize
  
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
  
  deallocate(option)
  nullify(option)
  
end subroutine OptionDestroy

end module Option_module
