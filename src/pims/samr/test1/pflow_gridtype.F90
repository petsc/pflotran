module pflow_gridtype_module
#include "include/finclude/petscsnes.h"
  use petscsnes
#include "definitions.h"
private

 type, public :: time_stepping_context
  
     real*8, pointer :: tfac(:)
     real*8 :: dt_min  ! Maximum size of the time step.
     real*8 :: dt_max  ! Maximum size of the time step.
         character*2 :: tunit ! Input time units
     integer  iaccel, icut_max, nstpmax, kplot  
     real*8, pointer :: tplot(:), tstep(:), dtstep(:)
     real*8 :: dpmxe,dsmxe !maximum allowed changes in field vars.
   

 
 end   type time_stepping_context


 type, public :: pflow_localpatch_info

  
    ! Local quantities
  
   
    integer*4 :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    integer*4 :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    integer*4 :: nxs, nys, nzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    integer*4 :: ngxs, ngys, ngzs
      ! Global indices of ghosted starting corner of local domain.
    integer*4 :: nxe, nye, nze, ngxe, ngye, ngze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    integer*4 :: nlxy, nlxz, nlyz
    integer*4 :: ngxy, ngxz, ngyz
    integer*4 :: nlmax  ! Total number of non-ghosted nodes in local domain.
    integer*4 :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    integer*4 :: nldof  ! nlmax times the number of phases.
    integer*4 :: ngdof  ! ngmax times the number of phases.
    integer*4 :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the local x-index of the non-ghosted starting (lower left)
      ! corner. iend gives the local x-index of the non-ghosted ending 
      ! corner. jstart, jend correspond to y-index, kstart, kend to z-index.
   

    ! Grid connections
    integer*4 :: nconn, nconnx, nconny
    integer*4, pointer :: nd1(:), nd2(:)
      ! Nodes upstream and downstream of a connection (assuming flow in 
      ! positive direction.  These are local, ghosted indices.

    integer*4, pointer :: iperm1(:), iperm2(:), ipermbc(:)

    real*8, pointer :: dist1(:),dist2(:),distbc(:),area(:),areabc(:), &
                       delzbc(:), vlbc(:), vvlbc(:),vgbc(:),vvgbc(:)
    real*8, pointer ::  delz(:) 


    integer*4, pointer :: nL2G(:), nG2L(:), nL2A(:),nG2N(:)

    ! Boundary conditions (BC's)
    integer*4 :: nblkbc
      ! The number of "blocks" of boundary conditions that are defined.
      ! Such a block is a specification of a set of boundary conditions.
      ! This set of boundary conditions can apply to any number of regions,
      ! so nblkbc does NOT equal the number of boundary condition regions.
    integer*4 :: nconnbc  ! The number of interfaces along boundaries.
    integer*4, pointer :: ibconn(:)
      ! ibconn(nc) specifies the index of the boundary condition block that
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
          ! For a Dirichlet BC, pressurebc(j,ibc) gives the partial pressure 
      ! for phase j along the BC block ibc.
    real*8, pointer :: velocitybc(:,:)
      ! For a Neumann BC, velocitybc(j,ibc) gives the velocity q for phase
      ! j along BC block ibc.
    real*8, pointer :: xxbc(:,:), varbc(:)
	
 !   real*8, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
	real*8, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)


!    real*8, allocatable :: Resold_AR(:,:), Resold_FL(:,:)
    real*8, pointer :: var(:) 
  PetscScalar, pointer ::accum_p(:)

  PetscScalar, pointer :: r_p(:), xx_loc_p(:), xx_p(:), yy_p(:),&
                 porosity_loc_p(:), volume_p(:), &
                 phis_p(:), tor_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:)
                          
               
  PetscScalar, pointer :: pc_p(:), pc_loc_p(:),kvr_p(:), kvr_loc_p(:)

  PetscScalar, pointer :: icap_p(:),&
                          icap_loc_p(:), ithrm_loc_p(:),ithrm_p(:)


 end type pflow_localpatch_info




  type, public:: pflowGrid
! Note that preprocessor directives MUST start in the first column!
!#ifndef DEBUG
!   private
!#endif

    integer :: myrank, commsize  ! Rank in PETSC_COMM_WORLD.
     integer*4 :: npx, npy, npz ! Processor partition in each direction.
    integer*4 :: nxy, nmax     ! nx * ny, nx * ny * nz
    integer :: nphase, nvar, ndof  ! Number of phases we are dealing with.
    integer :: size_var_use, size_var_node
	integer :: jh2o, jgas, joil ! specific phase indices


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
      ! If true, print the value of h at the end of each SNES iteration.
    PetscTruth :: use_ksp, Samrai_drive
    PetscTruth :: use_isoth, use_debug	
    ! If using_pflowGrid == PETSC_TRUE, then some parts of ptran_init 
    ! will not be executed, since they are made redundant by 
    ! pflowGrid_new() and pflowGrid_setup().

    real*8 :: t  ! The time elapsed in the simulation.
    real*8 :: dt ! The size of the time step.
    real*8 :: tconv ! Input time conversion factor
   
   
      ! An array of multiplicative factors that specify how to increase time step.
    integer :: flowsteps  ! The number of time-steps taken by the flow code.
    integer :: stepmax    ! The maximum number of time-steps taken by the flow code.
    integer :: nstpmax    ! The maximum number of time-step increments.
   ! integer :: kplot      ! Printout steps.
    integer :: write_init = 0 ! Flag to printout initial conditions.
    integer :: iprint = 0 ! Print level (-1-none, 0-fields, >=1-vel, 2-perm/por, 3-pflow.bc)
    integer :: imod = 1   ! screen printout  modulus
    integer :: itecplot = 0 ! tecplot print format (1-interchange x and z)
    integer :: iblkfmt = 0 ! blocked format
    integer :: isync = 0  ! Synchronize pflow and ptran time steps (1)
    integer :: ndtcmx = 5 ! Steps needed after cutting to increase time step
    integer :: newtcum    ! Total number of Newton steps taken.
    integer :: icutcum    ! Total number of cuts in the timestep taken.
    integer :: newton_max ! Max number of Newton steps for one time step.
    integer :: icut_max   ! Max number of dt cuts for one time step.
    integer :: iphch
    integer :: iread_init = 0 ! flag for reading initial conditions.
      ! Basically our target number of newton iterations per time step.
      real*8 :: dpmax,dsmax 
        
              	
    ! Grid topology
    integer :: igeom
	integer*4 :: nx, ny, nz    ! Global domain dimensions of the grid.

      ! Arrays for indexing between local ghosted and non-ghosted, local to natural arrays.
    DA :: da_1_dof, da_3np_dof, da_ndof
	  ! DA's for 1, 3, and multiple (number of phases) degrees of freedom.
      ! da_ndof = total degrees of freedom per node

    integer*4, pointer :: i1bc(:), i2bc(:), j1bc(:), j2bc(:), k1bc(:), k2bc(:)

    integer*4, pointer :: iregbc1(:), iregbc2(:)
      ! iregbc1(ibc) and iregbc2(ibc) give the id of the first region and 
      ! last region, respectively, that utilizes the boundary conditions in 
      ! boundary condition block ibc.


    !block BC values read from input
     real*8, pointer :: velocitybc0(:,:)
     real*8, pointer :: xxbc0(:,:)
     real*8 :: radius_0
!   phik
    integer :: iregperm, iran_por=0, iread_perm=0
    real*8 :: ran_fac=-1.d0
    integer*4, pointer :: i1reg(:),i2reg(:),j1reg(:),j2reg(:),k1reg(:),k2reg(:)
    real*8, pointer :: por_reg(:),tor_reg(:),perm_reg(:,:)

!   initial conditions
    integer :: iregini
    integer*4, pointer :: i1ini(:),i2ini(:),j1ini(:),j2ini(:),k1ini(:),k2ini(:)
    real*8, pointer :: xx_ini(:,:)
	
					                      
!   source term
    integer :: nblksrc = 0, ntimsrc = 0, isrc1 = 2
    integer*4, pointer :: i1src(:), i2src(:), j1src(:), j2src(:), k1src(:), k2src(:)
    real*8, pointer :: timesrc(:,:), tempsrc(:,:), qsrc(:,:,:)

!   solid reaction rate
    integer*4 :: ityprxn
    real*8 :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    real*8 ::qu_kin, yh2o_in_co2=0.D0

!   breakthrough curves
    integer :: ibrkcrv = 0
    integer*4, pointer :: i1brk(:),i2brk(:),j1brk(:),j2brk(:),k1brk(:),k2brk(:)
    integer*4, pointer :: ibrktyp(:),ibrkface(:)

!   dual continuum
    integer :: idcdm = 0, idcmblk = 0
    integer*4, pointer :: i1dcm(:),i2dcm(:),j1dcm(:),j2dcm(:),k1dcm(:),k2dcm(:)
    real*8, pointer :: fracture_aperture(:), matrix_block(:)

    integer*4, pointer :: icap_reg(:),ithrm_reg(:)
    real*8 :: scale
    real*8, pointer :: rock_density(:),cpr(:),dencpr(:),ckdry(:),ckwet(:), &
                       tau(:),cdiff(:),cexp(:)
    real*8, pointer :: swir(:),lambda(:),alpha(:),pckrm(:),pcwmax(:),pcbetac(:), &
                       pwrprm(:),sir(:,:)
    integer, pointer:: icaptype(:)
	

    integer :: ihydrostatic = 0,ideriv = 1
    real*8 :: dTdz,beta,tref,pref, gravity 

!   table lookup
    integer :: itable=0

    real*8, pointer :: dx0(:), dy0(:), dz0(:), rd(:)
    real*8, pointer :: x(:), y(:), z(:)

    !-------------------------------------------------------------------
    ! Quantities defined at each grid point.
    ! NOTE: I adopt the convention that _loc indicates the local portion
    ! of any global vector.
    !-------------------------------------------------------------------

    ! One degree of freedom: Physical coordinates.
    Vec :: conc
    Vec :: porosity, porosity0, porosity_loc, tor, tor_loc
    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings
    Vec :: volume  ! Volume of a cell in the grid
    Vec :: ithrm, ithrm_loc, icap, icap_loc
	Vec :: ttemp, ttemp_loc, temp ! 1 dof

    ! Three degrees of freedom:
!   Vec :: perm, perm_loc
    Vec :: perm_xx, perm_xx_loc, perm_yy, perm_yy_loc, perm_zz, perm_zz_loc
	Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    ! Multiple degrees of freedom (equal to number of phases present):
     Vec :: r             ! The residual.  (NOT the negative of the residual.)

    Vec :: vl! , vvl, vg, vvg ! phase (liquid and gas) velocities stored at interfaces
	
    ! Solution vectors
    Vec :: xx, xx_loc, dxx, yy, accum
        ! Jacobian matrix
    Mat :: J
    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

 	
    real*8 :: atol, rtol, stol, dtol
      ! Absolute, relative, and "change in norm of solution" tolerances.
    integer :: maxit, maxf
      ! The maximum number of iterations and function evaluations, respectively

 !  Vec :: p_nat, t_nat, c_nat, phis_nat, por_nat, vl_nat, s_nat, x_nat !, perm_nat
    ! Holds contents in natural ordering.
!   Vec :: p_all, t_all, c_all, phis_all, por_all, vl_all, s_all !, perm_all 
    ! Used to hold all values on processor 0.

  type(pflow_localpatch_info), pointer :: locpat(:)
	
  end type pflowGrid
  
  
type, public :: pflow_solver_context   
     ! Could you move SNES stuff out
    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType  :: pc_type
    KSP   ::  ksp
    PC    ::  pc

end type pflow_solver_context


 end module pflow_gridtype_module
