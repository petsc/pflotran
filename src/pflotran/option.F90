module Option_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

  type, public :: option_type 
  
    PetscMPIInt :: myrank                    ! rank in PETSC_COMM_WORLD
    PetscMPIInt :: commsize                  ! size of PETSC_COMM_WORLD
  
    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXNAMELENGTH) :: mode
    PetscInt :: imode
  
    PetscInt :: nphase, nvar, ndof  ! Number of phases we are dealing with.
    PetscInt :: jh2o, jgas, jco2 ! specific phase indices
    PetscInt :: nspec, npricomp

    ! Program options
    PetscTruth :: use_numerical  ! If true, use numerical Jacobian.
    PetscTruth :: use_matrix_free  ! If true, do not form the Jacobian.
      ! Note that if 'use_numerical' is true and 'use_matrix_free' is false,
      ! the Jacobian will be computed numerically and stored.
    PetscTruth :: print_hhistory
      ! If true, and if use_matrix_free is true, then store the differencing
      ! values h and print them out at the end of the simulation.

    PetscReal, pointer :: hhistory(:)
    PetscTruth :: monitor_h
    
    PetscInt :: idt_switch
    
    PetscTruth :: use_ksp
    PetscTruth :: use_isoth
    PetscTruth :: use_debug
    PetscTruth :: print_bcinfo
    
    PetscTruth :: print_convergence
    PetscTruth :: print_detailed_convergence
    PetscTruth :: check_infinity_norm
    PetscTruth :: force_at_least_1_iteration  
    
    character(len=MAXWORDLENGTH) :: generalized_grid
    logical :: use_generalized_grid
      
    ! If run_coupled == PETSC_TRUE, then some parts of ptran_init 
    ! will not be executed, since they are made redundant by 
    ! pflowGrid_new() and pflowGrid_setup().
    PetscTruth :: run_coupled

    PetscReal :: time  ! The time elapsed in the simulation.
    PetscReal :: dt ! The size of the time step.
  
!    PetscReal, pointer :: tplot(:)
    PetscReal, pointer :: tfac(:)
      ! An array of multiplicative factors that specify how to increase time step.
!    PetscInt :: kplot      ! Printout steps.
    PetscInt :: write_init = 0 ! Flag to printout initial conditions.
    PetscInt :: imod = 1   ! screen printout modulus
    PetscInt :: itecplot = 0 ! tecplot print format (1-interchange x and z)
    PetscInt :: iblkfmt = 0 ! blocked format
    PetscInt :: isync = 0  ! Synchronize pflow and ptran time steps (1)
  
    PetscInt :: iphch
    PetscInt :: iread_init = 0 ! flag for reading initial conditions.
      ! Basically our target number of newton iterations per time step.
    PetscReal :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    PetscReal :: dpmax,dtmpmax,dsmax,dcmax
    
    PetscInt :: nldof  ! nlmax times the number of phases.
    PetscInt :: ngdof  ! ngmax times the number of phases.

    PetscInt :: iran_por=0, iread_perm=0, iread_geom =1
    PetscReal :: ran_fac=-1.d0

!   solid reaction rate
    PetscInt :: ityprxn
    PetscReal :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    PetscReal ::qu_kin, yh2o_in_co2=0.D0
    
!   breakthrough curves
    PetscInt :: ibrkcrv = 0
    PetscInt, pointer :: ibrktyp(:),ibrkface(:)
    
!   dual continuum
    PetscInt :: idcdm = 0, idcmblk = 0
    PetscReal, pointer :: fracture_aperture(:), matrix_block(:)
    
    PetscInt, pointer :: icap_reg(:),ithrm_reg(:)
    PetscReal :: scale
    PetscReal, pointer :: rock_density(:),cpr(:),dencpr(:),ckdry(:),ckwet(:), &
                       tau(:),cdiff(:),cexp(:)
    PetscReal, pointer :: swir(:),lambda(:),alpha(:),pckrm(:),pcwmax(:),pcbetac(:), &
                       pwrprm(:),sir(:,:)
    PetscInt, pointer:: icaptype(:)

    PetscReal :: m_nacl
    PetscReal :: difaq, delhaq, gravity(3), fmwh2o= 18.0153D0, fmwa=28.96D0, &
              fmwco2=44.0098D0, eqkair, ret=1.d0, fc=1.d0
    
    PetscInt :: ihydrostatic = 0,ideriv = 1
    PetscReal :: dTdz,beta,tref,pref,conc0  ! these need to go
    PetscReal :: hydro_ref_xyz(3)
    
!   table lookup
    PetscInt :: itable=0

    PetscReal, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)
    
    PetscTruth :: restart_flag
    character(len=MAXWORDLENGTH) :: restart_file
    PetscTruth :: checkpoint_flag
    PetscInt :: checkpoint_frequency
    
    PetscLogDouble :: start_time
    PetscTruth :: wallclock_stop_flag
    PetscLogDouble :: wallclock_stop_time
    
    logical :: numerical_derivatives
    logical :: inexact_newton
    
  end type 
  
  type, public :: output_option_type

    character(len=2) :: tunit
    PetscReal :: tconv
  
    logical :: print_hdf5
    logical :: print_hdf5_velocities
    logical :: print_hdf5_flux_velocities

    logical :: print_tecplot 
    logical :: print_tecplot_velocities
    logical :: print_tecplot_flux_velocities

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

#include "definitions.h"
  
  type(option_type), pointer :: OptionCreate
  
  type(option_type), pointer :: option
  
  allocate(option)
  
  option%use_ksp = PETSC_FALSE
  option%use_isoth = PETSC_FALSE
  option%print_bcinfo = PETSC_FALSE

  option%use_numerical = PETSC_FALSE
  option%use_matrix_free = PETSC_FALSE
  
  option%print_convergence = PETSC_TRUE
  option%print_detailed_convergence = PETSC_FALSE
  option%check_infinity_norm = PETSC_TRUE
  option%force_at_least_1_iteration = PETSC_TRUE
  
  option%mode = ""
  option%imode = NULL_MODE
  option%idt_switch = 0
   
  option%run_coupled = PETSC_FALSE
  
!-----------------------------------------------------------------------
      ! Initialize some parameters to sensible values.  These are parameters
      ! which should be set via the command line or the input file, but it
      ! seems good practice to set them to sensible values when a pflowGrid
      ! is created.
!-----------------------------------------------------------------------
  allocate(option%tfac(13))
      
  option%tfac(1)  = 2.0d0; option%tfac(2)  = 2.0d0
  option%tfac(3)  = 2.0d0; option%tfac(4)  = 2.0d0
  option%tfac(5)  = 2.0d0; option%tfac(6)  = 1.8d0
  option%tfac(7)  = 1.6d0; option%tfac(8)  = 1.4d0
  option%tfac(9)  = 1.2d0; option%tfac(10) = 1.0d0
  option%tfac(11) = 1.0d0; option%tfac(12) = 1.0d0
  option%tfac(13) = 1.0d0      
  option%use_matrix_free = 1

  option%ihydrostatic = 0
  option%conc0 = 1.d-6
  
  option%dpmxe = 5.d4
  option%dtmpmxe = 2.d0
  option%dsmxe = 5.d0
  option%dcmxe = 5.d0

  !physical constants and defult variables
  option%difaq = 1.d-9 ! m^2/s read from input file
  option%delhaq = 12.6d0 ! kJ/mol read from input file
  option%gravity(:) = 0.d0
  option%gravity(3) = -9.8068d0    ! m/s^2
  option%tref   = 50.D0
  option%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  option%fmwco2 = 44.0098d0
  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  option%m_nacl = 0.d0
  
  option%generalized_grid = ""
  option%use_generalized_grid = .false.

  option%restart_flag = PETSC_FALSE
  option%restart_file = ""
  option%checkpoint_flag = PETSC_FALSE
  option%checkpoint_frequency = int(1d20)
  
  option%start_time = 0.d0
  option%wallclock_stop_flag = PETSC_FALSE
  option%wallclock_stop_time = 0.d0
  
  option%numerical_derivatives = .false.
  option%inexact_newton = .false.
    
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
  output_option%print_hdf5 = .false.
  output_option%print_hdf5_velocities = .false.
  output_option%print_hdf5_flux_velocities = .false.
  output_option%print_tecplot = .false.
  output_option%print_tecplot_velocities = .false.
  output_option%print_tecplot_flux_velocities = .false.
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
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_numerical", &
                           option%use_numerical, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
                           option%print_hhistory, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
                           option%monitor_h, ierr) 
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_debug", &
                           option%use_debug, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
                           option%use_ksp, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isoth", &
                           option%use_isoth, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_bcinfo", &
                           option%print_bcinfo, ierr)
                           
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', option%restart_file, &
                             option%restart_flag, ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-chkptfreq', &
                          option%checkpoint_frequency, &
                          option%checkpoint_flag, ierr)                           
                             
  ! check on possible modes                                                     
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richards", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "richards"                           
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richards_lite", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "richards_lite"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_liquid", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "liquid"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_cond", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "cond"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_th", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "th"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_thc", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "thc"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_2ph", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "2ph"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_mph", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "mph"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_flash", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "flash"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_owg", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "owg"                           
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_vadose", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "vadose"                                                     
 
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
    print *, 'ERROR: ' // string
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
  
  if (option%myrank == 0) print *, 'WARNING: ' // string
  
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
  
  if (option%myrank == 0) print *, string
  
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
  
  deallocate(option)
  nullify(option)
  
end subroutine OptionDestroy

end module Option_module
