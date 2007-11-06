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
  
    integer :: myrank                    ! rank in PETSC_COMM_WORLD
    integer :: commsize                  ! size of PETSC_COMM_WORLD
  
    ! defines the mode (e.g. mph, richards, vadose, etc.
    character(len=MAXNAMELENGTH) :: mode
    integer :: imode
    
    PetscTruth :: restartflag
    character(len=MAXNAMELENGTH) :: restartfile
  
    integer :: nphase, nvar, ndof  ! Number of phases we are dealing with.
    integer :: jh2o, jgas, jco2 ! specific phase indices
    integer :: nspec, npricomp

    ! Program options
    PetscTruth :: use_analytical  ! If true, use analytical Jacobian.
    PetscTruth :: use_matrix_free  ! If true, do not form the Jacobian.
      ! Note that if 'use_analytical' and 'use_matrix_free' are both false,
      ! the Jacobian will be computed numerically and stored.
    PetscTruth :: print_hhistory
      ! If true, and if use_matrix_free is true, then store the differencing
      ! values h and print them out at the end of the simulation.

    PetscScalar, pointer :: hhistory(:)
    PetscTruth :: monitor_h
    
    PetscTruth :: use_ksp
    PetscTruth :: use_isoth
    PetscTruth :: use_debug
    PetscTruth :: print_bcinfo
    
#if 0
      ! If true, print the value of h at the end of each SNES iteration.
    PetscTruth :: use_liquid, use_cond, use_th, use_thc, use_2ph, &
    use_mph, use_ksp, use_owg, use_vadose, use_flash
    PetscTruth :: use_isoth, use_debug, use_richards 
#endif     
    ! If run_coupled == PETSC_TRUE, then some parts of ptran_init 
    ! will not be executed, since they are made redundant by 
    ! pflowGrid_new() and pflowGrid_setup().
    PetscTruth :: run_coupled

    real*8 :: t  ! The time elapsed in the simulation.
    real*8 :: dt ! The size of the time step.
    real*8 :: dt_min  ! Maximum size of the time step.
    real*8 :: dt_max  ! Maximum size of the time step.
    real*8 :: tconv ! Input time conversion factor
    character*2 :: tunit ! Input time units
    real*8, pointer :: tplot(:), tstep(:), dtstep(:)
    real*8, pointer :: tfac(:)
      ! An array of multiplicative factors that specify how to increase time step.
    integer :: flowsteps  ! The number of time-steps taken by the flow code.
    integer :: stepmax    ! The maximum number of time-steps taken by the flow code.
    integer :: nstpmax    ! The maximum number of time-step increments.
    integer :: kplot      ! Printout steps.
    integer :: write_init = 0 ! Flag to printout initial conditions.
    integer :: iprint = 0 ! Print level (-1-none, 0-fields, >=1-vel, 2-perm/por, 3-pflow.bc)
    logical :: print_hdf5 = .false. ! toggle for printing hdf5
    logical :: print_hdf5_velocities = .false.
    logical :: print_hdf5_flux_velocities = .false.
    logical :: print_tecplot = .false. ! toggle for printing tecplot
    logical :: print_tecplot_velocities = .false.
    logical :: print_tecplot_flux_velocities = .false.
    integer :: imod = 1   ! screen printout modulus
    integer :: itecplot = 0 ! tecplot print format (1-interchange x and z)
    integer :: iblkfmt = 0 ! blocked format
    integer :: isync = 0  ! Synchronize pflow and ptran time steps (1)
    integer :: ndtcmx = 5 ! Steps needed after cutting to increase time step
    integer :: newtcum    ! Total number of Newton steps taken.
    integer :: icutcum    ! Total number of cuts in the timestep taken.
    integer :: newton_max ! Max number of Newton steps for one time step.
    integer :: icut_max   ! Max number of dt cuts for one time step.
    integer :: iaccel,iphch
    integer :: iread_init = 0 ! flag for reading initial conditions.
      ! Basically our target number of newton iterations per time step.
    real*8 :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    real*8 :: dpmax,dtmpmax,dsmax,dcmax
    
    integer*4 :: nldof  ! nlmax times the number of phases.
    integer*4 :: ngdof  ! ngmax times the number of phases.

    integer*4, pointer :: iperm1(:), iperm2(:), ipermbc(:)

    real*8, pointer :: density_bc(:),d_p_bc(:),d_t_bc(:), d_s_bc(:),d_c_bc(:),&
                       avgmw_bc(:),avgmw_c_bc(:),&
                       hh_bc(:),h_p_bc(:),h_t_bc(:),h_s_bc(:), h_c_bc(:), &
                       viscosity_bc(:),v_p_bc(:),v_t_bc(:),&
                       uu_bc(:),u_p_bc(:),u_t_bc(:),u_s_bc(:), u_c_bc(:),&    
                       df_bc(:),df_p_bc(:),df_t_bc(:),df_s_bc(:), df_c_bc(:), &
                       hen_bc(:),hen_p_bc(:),hen_t_bc(:),hen_s_bc(:),hen_c_bc(:), &
                       pc_bc(:),pc_p_bc(:),pc_t_bc(:),pc_s_bc(:),pc_c_bc(:), &
                       kvr_bc(:),kvr_p_bc(:),kvr_t_bc(:),kvr_s_bc(:),kvr_c_bc(:)
    real*8, pointer :: xphi_co2(:),xxphi_co2(:),den_co2(:), dden_co2(:)
    
    ! Boundary conditions (BC's)
    integer*4 :: nblkbc
      ! The number of "blocks" of boundary conditions that are defined.
      ! Such a block is a specification of a set of boundary conditions.
      ! This set of boundary conditions can apply to any number of regions,
      ! so nblkbc does NOT equal the number of boundary condition regions.
!    integer*4 :: nconnbc  ! The number of interfaces along boundaries.
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1bc(:), i2bc(:), j1bc(:), j2bc(:), k1bc(:), k2bc(:)
!GEH - Structured Grid Dependence - End
!geh - now in grid
!    integer*4, pointer :: ibconn(:)
      ! ibconn(nc) specifies the id of the of boundary condition block that
      ! applies at boundary interface nc.  
    integer*4, pointer :: ibndtyp(:)
      ! ibndtyp(ibc) specifies the type of boundary condition that applies
      ! for boundary condition block ibc.
    integer*4, pointer :: iface(:)
      ! iface(ibc) specifies the face (left, right, top, bottom, etc.) on
      ! which BC block ibc lies.
    integer*4, pointer :: mblkbc(:)
      ! mblkbc(nc) gives the local, non-ghosted index of the cell that has
      ! boundary connection nc.
    integer*4, pointer :: iregbc1(:), iregbc2(:)
      ! iregbc1(ibc) and iregbc2(ibc) give the id of the first region and 
      ! last region, respectively, that utilizes the boundary conditions in 
      ! boundary condition block ibc.
    real*8, pointer :: pressurebc(:,:)
      ! For a Dirichlet BC, pressurebc(j,ibc) gives the partial pressure 
      ! for phase j along the BC block ibc.
    real*8, pointer :: velocitybc(:,:)
      ! For a Neumann BC, velocitybc(j,ibc) gives the velocity q for phase
      ! j along BC block ibc.
    real*8, pointer :: tempbc(:),concbc(:),sgbc(:),xphi_co2_bc(:),xxphi_co2_bc(:)
    real*8, pointer :: xxbc(:,:), varbc(:)
    integer, pointer:: iphasebc(:)

    !block BC values read from input
    real*8, pointer :: pressurebc0(:,:)
    real*8, pointer :: velocitybc0(:,:)
    real*8, pointer :: tempbc0(:),concbc0(:),sgbc0(:)
    real*8, pointer :: xxbc0(:,:)
    integer, pointer:: iphasebc0(:)  

    integer :: iran_por=0, iread_perm=0, iread_geom =1
    real*8 :: ran_fac=-1.d0
#if 0
!   phik
    integer :: iregperm
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1reg(:),i2reg(:),j1reg(:),j2reg(:),k1reg(:),k2reg(:)
!GEH - Structured Grid Dependence - End
    real*8, pointer :: por_reg(:),tor_reg(:),perm_reg(:,:)
#endif
#if 0
!   initial conditions
    integer :: iregini
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1ini(:),i2ini(:),j1ini(:),j2ini(:),k1ini(:),k2ini(:)
!GEH - Structured Grid Dependence - End
    real*8, pointer :: pres_ini(:),temp_ini(:),conc_ini(:),sat_ini(:), &
                       xmol_ini(:)
    real*8, pointer :: xx_ini(:,:)
    integer, pointer:: iphas_ini(:)
#endif

!   source term
    integer :: nblksrc = 0, ntimsrc = 0, isrc1 = 2
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1src(:), i2src(:), j1src(:), j2src(:), k1src(:), k2src(:)
!GEH - Structured Grid Dependence - End
    real*8, pointer :: timesrc(:,:), tempsrc(:,:), qsrc(:,:), csrc(:,:), hsrc(:,:)
    
!   solid reaction rate
    integer*4 :: ityprxn
    real*8 :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    real*8 ::qu_kin, yh2o_in_co2=0.D0
    
!   breakthrough curves
    integer :: ibrkcrv = 0
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1brk(:),i2brk(:),j1brk(:),j2brk(:),k1brk(:),k2brk(:)
!GEH - Structured Grid Dependence - End
    integer*4, pointer :: ibrktyp(:),ibrkface(:)
    
!   dual continuum
    integer :: idcdm = 0, idcmblk = 0
!GEH - Structured Grid Dependence - Begin
    integer*4, pointer :: i1dcm(:),i2dcm(:),j1dcm(:),j2dcm(:),k1dcm(:),k2dcm(:)
!GEH - Structured Grid Dependence - End
    real*8, pointer :: fracture_aperture(:), matrix_block(:)
    
    integer*4, pointer :: icap_reg(:),ithrm_reg(:)
    real*8 :: scale
    real*8, pointer :: rock_density(:),cpr(:),dencpr(:),ckdry(:),ckwet(:), &
                       tau(:),cdiff(:),cexp(:)
    real*8, pointer :: swir(:),lambda(:),alpha(:),pckrm(:),pcwmax(:),pcbetac(:), &
                       pwrprm(:),sir(:,:)
    integer, pointer:: icaptype(:)
!geh material id
    integer, pointer :: imat(:)
    real*8 :: m_nacl
    real*8 :: difaq, delhaq, gravity, fmwh2o= 18.0153D0, fmwa=28.96D0, &
              fmwco2=44.0098D0, eqkair, ret=1.d0, fc=1.d0
    
    integer :: ihydrostatic = 0,ideriv = 1
    real*8 :: dTdz,beta,tref,pref,conc0
    real*8 :: hydro_ref_xyz(3)
    
!   table lookup
    integer :: itable=0

    !-------------------------------------------------------------------
    ! Quantities defined at each grid point.
    ! NOTE: I adopt the convention that _loc indicates the local portion
    ! of any global vector.
    !-------------------------------------------------------------------

    ! One degree of freedom: Physical coordinates.
    Vec :: conc
    Vec :: porosity, porosity0, porosity_loc, tor, tor_loc
!GEH - Structured Grid Dependence - Begin
!transferred to structured_grid.F90    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings
!GEH - Structured Grid Dependence - End
    Vec :: volume  ! Volume of a cell in the grid
    Vec :: ithrm, ithrm_loc, icap, icap_loc, iphas, iphas_loc, iphas_old
    Vec :: ttemp, ttemp_loc, temp ! 1 dof
    Vec :: phis

    ! Three degrees of freedom:
!   Vec :: perm, perm_loc
    Vec :: perm_xx, perm_xx_loc, perm_yy, perm_yy_loc, perm_zz, perm_zz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    ! Multiple degrees of freedom (equal to number of phases present):
    Vec :: var,var_loc
    Vec :: ppressure, ppressure_loc, pressure, dp
    Vec :: ssat, ssat_loc, sat    ! saturation
    Vec :: xxmol, xxmol_loc, xmol ! mole fraction
    Vec :: density       ! Density at time k
    Vec :: ddensity, ddensity_loc  ! Density at time k+1
    Vec :: d_p, d_p_loc  ! dD/dp at time k+1
    Vec :: d_t, d_t_loc  ! dD/dT at time k+1
    Vec :: d_c, d_c_loc  ! dD/dT at time k+1
    Vec :: d_s, d_s_loc  ! dD/dT at time k+1
    Vec :: avgmw,avgmw_loc  ! Density at time k+1molecular weight at time k+1
    Vec :: avgmw_c,avgmw_c_loc
    Vec :: h             ! H     at time k
    Vec :: hh, hh_loc    ! H     at time k+1
    Vec :: h_p, h_p_loc  ! dH/dp at time k+1
    Vec :: h_t, h_t_loc  ! dH/dT at time k+1
    Vec :: h_c, h_c_loc  ! dD/dT at time k+1
    Vec :: h_s, h_s_loc  ! dD/dT at time k+1
    Vec :: u            ! H     at time k
    Vec :: uu, uu_loc    ! H     at time k+1
    Vec :: u_p, u_p_loc  ! dH/dp at time k+1
    Vec :: u_t, u_t_loc  ! dH/dT at time k+1
    Vec :: u_c, u_c_loc  ! dD/dT at time k+1
    Vec :: u_s, u_s_loc  ! dD/dT at time k+1
    Vec :: hen, hen_loc    ! H     at time k+1
    Vec :: hen_p, hen_p_loc  ! dH/dp at time k+1
    Vec :: hen_t, hen_t_loc  ! dH/dT at time k+1
    Vec :: hen_c, hen_c_loc  ! dD/dT at time k+1
    Vec :: hen_s, hen_s_loc  ! dD/dT at time k+1
    Vec :: df, df_loc    ! H     at time k+1
    Vec :: df_p, df_p_loc  ! dH/dp at time k+1
    Vec :: df_t, df_t_loc  ! dH/dT at time k+1
    Vec :: df_c, df_c_loc  ! dD/dT at time k+1
    Vec :: df_s, df_s_loc  ! dD/dT at time k+1
    Vec :: viscosity, viscosity_loc  !kept for early routine
   
    Vec :: v_p, v_p_loc  ! dv/dp at time k+1
    Vec :: v_t, v_t_loc  ! dv/dT at time k+1
    Vec :: pcw, pcw_loc    ! H     at time k+1
    Vec :: pc_p, pc_p_loc  ! dH/dp at time k+1
    Vec :: pc_t, pc_t_loc  ! dH/dT at time k+1
    Vec :: pc_c, pc_c_loc  ! dD/dT at time k+1
    Vec :: pc_s, pc_s_loc  ! dD/dT at time k+1
    Vec :: kvr, kvr_loc    ! H     at time k+1
    Vec :: kvr_p, kvr_p_loc  ! d/dp at time k+1
    Vec :: kvr_t, kvr_t_loc  ! dm/dT at time k+1
    Vec :: kvr_c, kvr_c_loc  ! d/d at time k+1
    Vec :: kvr_s, kvr_s_loc  ! dD/dT at time k+1
    Vec :: r             ! The residual.  (NOT the negative of the residual.)

    Vec :: vl, vvl, vg, vvg ! phase (liquid and gas) velocities stored at interfaces


    real*8, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
    real*8, pointer :: vvlbc(:), vvgbc(:)
    real*8, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)

    ! Solution vectors
    Vec :: xx, xx_loc, dxx, yy, accum
        ! Jacobian matrix
    Mat :: J
    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType  :: pc_type
    KSP   ::  ksp
    PC    ::  pc
   

    integer var_plot_num
    character*16, pointer :: var_plot_nm(:)
    integer, pointer :: var_plot_ind(:)  
   
  end type option_type
  
  public :: OptionCreate, &
            OptionCheckCommandLine, &
            printErrMsg, &
            printWrnMsg, &
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

#include "definitions.h"
  
  type(option_type), pointer :: OptionCreate
  
  type(option_type), pointer :: option
  
  allocate(option)
  
  option%use_ksp = PETSC_FALSE
  option%use_isoth = PETSC_FALSE
  option%print_bcinfo = PETSC_FALSE
  
  option%mode = ""
  option%imode = NULL_MODE
   
  option%run_coupled = PETSC_FALSE
  
!-----------------------------------------------------------------------
      ! Initialize some counter variables.
!-----------------------------------------------------------------------
  option%t = 0.d0
  option%flowsteps = 0
  option%newtcum = 0
  option%icutcum = 0

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
  option%newton_max = 16
  option%icut_max = 16
  option%iaccel = 1
  option%dt = 1.d0
  option%dt_min = 1.d0
  option%dt_max = 3.1536d6 ! One-tenth of a year.

  option%ihydrostatic = 0
  option%conc0 = 1.d-6
  
  option%dpmxe = 5.d4
  option%dtmpmxe = 2.d0
  option%dsmxe = 5.d0
  option%dcmxe = 5.d0

  !physical constants and defult variables
  option%difaq = 1.d-9 ! m^2/s read from input file
  option%delhaq = 12.6d0 ! kJ/mol read from input file
  option%gravity = 9.8068d0    ! m/s^2
  option%tref   = 50.D0
  option%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  option%fmwco2 = 44.0098d0
  option%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  OptionCreate => option
  
end function OptionCreate

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
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_analytical", &
                           option%use_analytical, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_hhistory", &
                           option%print_hhistory, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-monitor_h", &
                           option%monitor_h, ierr) 
  call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-restart', &
                             option%restartfile, &
                             option%restartflag, ierr)     
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_debug", &
                           option%use_debug, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_ksp", &
                           option%use_ksp, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_isoth", &
                           option%use_isoth, ierr)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-print_bcinfo", &
                           option%print_bcinfo, ierr)
                             
  ! check on possible modes                                                     
  option_found = PETSC_FALSE
  call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-use_richards", &
                           option_found, ierr)
  if (option_found == PETSC_TRUE) option%mode = "richards"                           
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
  
  if (option%myrank == 0) then
    print *, 'ERROR: ', string
    print *, 'Stopping!'
  endif    
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
  
  if (option%myrank == 0) print *, 'WARNING: ', string
  
end subroutine printWrnMsg

! ************************************************************************** !
!
! OptionDestroy: Deallocates an option
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine OptionDestroy(option)

  implicit none
  
  type(option_type) :: option
  
  ! all kinds of stuff needs to be added here.
  
end subroutine OptionDestroy

end module Option_module
