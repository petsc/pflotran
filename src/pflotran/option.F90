module Option_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"


  type, public :: option_type 
  
    PetscMPIInt :: comm                      ! alternative to PETSC_COMM_WORLD
    PetscMPIInt :: myrank                    ! rank in PETSC_COMM_WORLD
    PetscMPIInt :: commsize                  ! size of PETSC_COMM_WORLD
    PetscTruth :: broadcast_read
  
    PetscInt :: fid_in
    PetscInt :: fid_out
    
    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXWORDLENGTH) :: flowmode
    PetscInt :: iflowmode
    character(len=MAXWORDLENGTH) :: tranmode
    PetscInt :: itranmode

    PetscInt :: nphase
    PetscInt :: nflowdof
    PetscInt :: nflowspec

    PetscInt :: ntrandof
  
    PetscInt :: iflag
    
    PetscInt, pointer :: garbage ! for some reason, Intel will not compile without this

    PetscReal :: uniform_velocity(3)

    ! Program options
    PetscTruth :: use_matrix_free  ! If true, do not form the Jacobian.
    
    PetscInt :: imod
    
    PetscTruth :: use_isoth
    
    character(len=MAXWORDLENGTH) :: generalized_grid
    logical :: use_generalized_grid
      
    PetscReal :: flow_time, tran_time, time  ! The time elapsed in the simulation.
    PetscReal :: tran_weight_t0, tran_weight_t1
    PetscReal :: flow_dt, tran_dt, dt ! The size of the time step.
    PetscTruth :: match_waypoint
    PetscReal :: prev_dt
  
!    PetscReal, pointer :: tplot(:)
    PetscReal, pointer :: tfac(:)
      ! An array of multiplicative factors that specify how to increase time step.
        
      ! Basically our target number of newton iterations per time step.
    PetscReal :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    PetscReal :: dpmax,dtmpmax,dsmax,dcmax
    
    PetscReal :: scale
    PetscReal, pointer :: dencpr(:),ckwet(:)
    PetscReal, pointer :: sir(:,:)


    PetscReal :: m_nacl
    PetscReal :: difaq, delhaq, gravity(3), fmwh2o= 18.0153D0, fmwa=28.96D0, &
              fmwco2=44.0098D0, eqkair, ret=1.d0, fc=1.d0
    
    PetscInt :: ideriv
    PetscInt :: idt_switch = -1
    PetscReal :: tref
    PetscReal :: pref
    PetscReal :: den_ref
        
    PetscReal :: disp
    
!   table lookup
    PetscInt :: itable=0
    PetscInt :: co2eos=EOS_SPAN_WAGNER

    PetscTruth :: restart_flag
    PetscReal :: restart_time
    character(len=MAXWORDLENGTH) :: restart_file
    PetscTruth :: checkpoint_flag
    PetscInt :: checkpoint_frequency
    
    PetscLogDouble :: start_time
    PetscTruth :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time
    
    PetscInt :: log_stage(10)
    
    logical :: numerical_derivatives
    logical :: compute_statistics
    logical :: compute_mass_balance
    logical :: use_touch_options
    logical :: overwrite_restart_transport
    logical :: overwrite_restart_flow_params
    PetscInt :: io_handshake_buffer_size
    
    character(len=MAXWORDLENGTH) :: permx_filename
    character(len=MAXWORDLENGTH) :: permy_filename
    character(len=MAXWORDLENGTH) :: permz_filename
    
  end type option_type
  
  type, public :: output_option_type

    character(len=2) :: tunit
    PetscReal :: tconv
  
    logical :: print_hdf5
    logical :: print_hdf5_velocities
    logical :: print_hdf5_flux_velocities

    logical :: print_tecplot 
    PetscInt :: tecplot_format
    logical :: print_tecplot_velocities
    logical :: print_tecplot_flux_velocities
    
    logical :: print_vtk 
    logical :: print_vtk_velocities

    PetscInt :: plot_number
    character(len=MAXWORDLENGTH) :: plot_name

  end type output_option_type

  interface OptionDotProduct
    module procedure OptionDotProduct1
    module procedure OptionDotProduct2
    module procedure OptionDotProduct3
  end interface
  
  public :: OptionCreate, &
            OutputOptionCreate, &
            OptionCheckCommandLine, &
            printErrMsg, &
            printWrnMsg, &
            printMsg, &
            OptionDotProduct, &
            OptionDestroy, &
            OptionCheckTouch, &
            OutputOptionDestroy

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

  option%comm = 0
  option%myrank = 0
  option%commsize = 0
  option%broadcast_read = PETSC_FALSE
  
  option%fid_in = 0
  option%fid_out = 0

  option%iflag = 0
  
  option%iflag = 0
  option%use_isoth = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  
  option%flowmode = ""
  option%iflowmode = NULL_MODE
  option%nflowdof = 0

  option%tranmode = ""
  option%itranmode = NULL_MODE
  option%ntrandof = 0
  
  option%uniform_velocity = 0.d0
  option%imod = 1
   
!-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  option%pref = 101325.d0
  option%tref = 25.d0
  option%den_ref = 0.d0
  
  allocate(option%tfac(13))
      
  option%tfac(1)  = 2.0d0; option%tfac(2)  = 2.0d0
  option%tfac(3)  = 2.0d0; option%tfac(4)  = 2.0d0
  option%tfac(5)  = 2.0d0; option%tfac(6)  = 1.8d0
  option%tfac(7)  = 1.6d0; option%tfac(8)  = 1.4d0
  option%tfac(9)  = 1.2d0; option%tfac(10) = 1.0d0
  option%tfac(11) = 1.0d0; option%tfac(12) = 1.0d0
  option%tfac(13) = 1.0d0    
  
  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  option%use_matrix_free = PETSC_FALSE
  option%ideriv = 1

  option%dpmxe = 5.d4
  option%dtmpmxe = 5.d0
  option%dsmxe = 0.5d0
  option%dcmxe = 1.d0

  !physical constants and defult variables
!  option%difaq = 1.d-9 ! m^2/s read from input file
  option%difaq = 0.d0
  option%delhaq = 12.6d0 ! kJ/mol read from input file
  option%gravity(:) = 0.d0
  option%gravity(3) = -9.8068d0    ! m/s^2
  option%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  option%fmwco2 = 44.0098d0
  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  option%m_nacl = 0.d0
  
  option%disp = 0.d0
  
  option%generalized_grid = ""
  option%use_generalized_grid = PETSC_FALSE

  option%restart_flag = PETSC_FALSE
  option%restart_file = ""
  option%restart_time = -999.d0
  option%checkpoint_flag = PETSC_FALSE
  option%checkpoint_frequency = huge(option%checkpoint_frequency)
  
  option%start_time = 0.d0
  option%wallclock_stop_flag = PETSC_FALSE
  option%wallclock_stop_time = 0.d0
  
  option%log_stage = 0
  
  option%numerical_derivatives = PETSC_FALSE
  option%compute_statistics = PETSC_FALSE
  option%compute_mass_balance = PETSC_FALSE

  option%use_touch_options = PETSC_FALSE
  option%overwrite_restart_transport = PETSC_FALSE
  option%overwrite_restart_flow_params = PETSC_FALSE

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

  OptionCreate => option
  
end function OptionCreate

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
  output_option%plot_number = 0
  output_option%plot_name = ""

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
  PetscErrorCode :: ierr
  
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-snes_mf", & 
                           option%use_matrix_free, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isoth", &
                           option%use_isoth, ierr)
                           
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', option%restart_file, &
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
 
end subroutine OptionCheckCommandLine

! ************************************************************************** !
!
! printErrMsg: Prints the error message from p0 and stops
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printErrMsg(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  PetscErrorCode :: ierr
  
  if (option%myrank == 0) then
    print *
    print *, 'ERROR: ' // trim(string)
    print *, 'Stopping!'
  endif    
  call PetscFinalize(ierr)
  stop
  
end subroutine printErrMsg
 
! ************************************************************************** !
!
! printWrnMsg: Prints the warning message from p0
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine printWrnMsg(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  if (option%myrank == 0) print *, 'WARNING: ' // trim(string)
  
end subroutine printWrnMsg

! ************************************************************************** !
!
! printMsg: Prints the message from p0
! author: Glenn Hammond
! date: 11/14/07
!
! ************************************************************************** !
subroutine printMsg(option,string)

  implicit none
  
  type(option_type) :: option
  character(len=*) :: string
  
  if (option%myrank == 0) print *, trim(string)
  
end subroutine printMsg

! ************************************************************************** !
!
! OptionDotProduct1: Computes the dot product between two 3d vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function OptionDotProduct1(v1,v2)

  implicit none
  
  PetscReal :: v1(3), v2(3)
  
  PetscReal :: OptionDotProduct1
  
  OptionDotProduct1 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

end function OptionDotProduct1

! ************************************************************************** !
!
! OptionDotProduct2: Computes the dot product between two 3d vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function OptionDotProduct2(v1,v2x,v2y,v2z)

  implicit none
  
  PetscReal :: v1(3), v2x, v2y, v2z
  
  PetscReal :: OptionDotProduct2
  
  OptionDotProduct2 = v1(1)*v2x+v1(2)*v2y+v1(3)*v2z

end function OptionDotProduct2

! ************************************************************************** !
!
! OptionDotProduct3: Computes the dot product between components of two 3d 
!                    vectors
! author: Glenn Hammond
! date: 11/28/07
!
! ************************************************************************** !
function OptionDotProduct3(v1x,v1y,v1z,v2x,v2y,v2z)

  implicit none
  
  PetscReal :: v1x, v1y, v1z, v2x, v2y, v2z
  
  PetscReal :: OptionDotProduct3
  
  OptionDotProduct3 = v1x*v2x+v1y*v2y+v1z*v2z

end function OptionDotProduct3

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
  character(len=MAXWORDLENGTH) :: filename
  
  PetscInt :: ios
  PetscInt :: fid = 86
  logical :: OptionCheckTouch
  PetscErrorCode :: ierr
  
  OptionCheckTouch = PETSC_FALSE

  if (option%myrank == 0) &
    open(unit=fid,file=trim(filename),status='old',iostat=ios)
  call MPI_Bcast(ios,1,MPI_INTEGER,ZERO_INTEGER,option%comm,ierr)

  if (ios == 0) then
    if (option%myrank == 0) close(fid,status='delete')
    OptionCheckTouch = PETSC_TRUE
  endif

end function OptionCheckTouch

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
  if (associated(option%dencpr)) deallocate(option%dencpr)
  nullify(option%dencpr)
  if (associated(option%ckwet)) deallocate(option%ckwet)
  nullify(option%ckwet)
  if (associated(option%sir)) deallocate(option%sir)
  nullify(option%sir)
  if (associated(option%tfac)) deallocate(option%tfac)
  nullify(option%tfac)
  
  deallocate(option)
  nullify(option)
  
end subroutine OptionDestroy

end module Option_module
