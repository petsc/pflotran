module pflow_gridtype_module


private
#include "definitions.h"
! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
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

  type, public:: pflowGrid

    PetscMPIInt :: myrank, commsize  ! Rank in PETSC_COMM_WORLD.

    PetscInt :: nphase, nvar, ndof  ! Number of phases we are dealing with.
    PetscInt :: jh2o, jgas, jco2 ! specific phase indices
    PetscInt :: nspec, npricomp

    ! Program options
    PetscTruth :: use_analytical  ! If true, use analytical Jacobian.
    PetscTruth :: use_matrix_free  ! If true, do not form the Jacobian.
      ! Note that if 'use_analytical' and 'use_matrix_free' are both false,
      ! the Jacobian will be computed numerically and stored.
    PetscTruth :: print_hhistory
      ! If true, and if use_matrix_free is true, then store the differencing
      ! values h and print them out at the end of the simulation.

    PetscReal, pointer :: hhistory(:)
    PetscTruth :: monitor_h
      ! If true, print the value of h at the end of each SNES iteration.
    PetscTruth :: use_liquid, use_cond, use_th, use_thc, use_2ph, &
    use_mph, use_ksp, use_owg, use_vadose, use_flash
    PetscTruth :: use_isoth, use_debug, use_richards  
    ! If using_pflowGrid == PETSC_TRUE, then some parts of ptran_init 
    ! will not be executed, since they are made redundant by 
    ! pflowGrid_new() and pflowGrid_setup().
    PetscTruth :: using_pflowGrid = PETSC_FALSE

    PetscReal :: t  ! The time elapsed in the simulation.
    PetscReal :: dt ! The size of the time step.
    PetscReal :: dt_min  ! Maximum size of the time step.
    PetscReal :: dt_max  ! Maximum size of the time step.
    PetscReal :: tconv ! Input time conversion factor
    character*2 :: tunit ! Input time units
    PetscReal, pointer :: tplot(:), tstep(:), dtstep(:)
    PetscReal, pointer :: tfac(:)
      ! An array of multiplicative factors that specify how to increase time step.
    PetscInt :: flowsteps  ! The number of time-steps taken by the flow code.
    PetscInt :: stepmax    ! The maximum number of time-steps taken by the flow code.
    PetscInt :: nstpmax    ! The maximum number of time-step increments.
    PetscInt :: kplot      ! Printout steps.
    PetscInt :: write_init = 0 ! Flag to printout initial conditions.
    PetscInt :: iprint = 0 ! Print level (-1-none, 0-fields, >=1-vel, 2-perm/por, 3-pflow.bc)
    logical :: print_hdf5 = .false. ! toggle for printing hdf5
    logical :: print_hdf5_velocities = .false.
    logical :: print_hdf5_flux_velocities = .false.
    logical :: print_tecplot = .false. ! toggle for printing tecplot
    logical :: print_tecplot_velocities = .false.
    logical :: print_tecplot_flux_velocities = .false.
    PetscInt :: imod = 1   ! screen printout modulus
    PetscInt :: itecplot = 0 ! tecplot print format (1-interchange x and z)
    PetscInt :: iblkfmt = 0 ! blocked format
    PetscInt :: isync = 0  ! Synchronize pflow and ptran time steps (1)
    PetscInt :: ndtcmx = 5 ! Steps needed after cutting to increase time step
    PetscInt :: newtcum    ! Total number of Newton steps taken.
    PetscInt :: icutcum    ! Total number of cuts in the timestep taken.
    PetscInt :: newton_max ! Max number of Newton steps for one time step.
    PetscInt :: icut_max   ! Max number of dt cuts for one time step.
    PetscInt :: iaccel,iphch
    PetscInt :: iread_init = 0 ! flag for reading initial conditions.
      ! Basically our target number of newton iterations per time step.
    PetscReal :: dpmxe,dtmpmxe,dsmxe,dcmxe !maximum allowed changes in field vars.
    PetscReal :: dpmax,dtmpmax,dsmax,dcmax

    ! Grid topology
    PetscInt :: igeom
    
!GEH - Structured Grid Dependence - Begin
    PetscInt :: nx, ny, nz    ! Global domain dimensions of the grid.
    PetscInt :: nxy, nmax     ! nx * ny, nx * ny * nz
    PetscInt :: npx, npy, npz ! Processor partition in each direction.
    PetscInt :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    PetscInt :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    PetscInt :: nxs, nys, nzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    PetscInt :: ngxs, ngys, ngzs
      ! Global indices of ghosted starting corner of local domain.
    PetscInt :: nxe, nye, nze, ngxe, ngye, ngze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    PetscInt :: nlxy, nlxz, nlyz
    PetscInt :: ngxy, ngxz, ngyz
!GEH - Structured Grid Dependence - End

    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt :: nldof  ! nlmax times the number of phases.
    PetscInt :: ngdof  ! ngmax times the number of phases.

!GEH - Structured Grid Dependence - Begin
    PetscInt :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the local x-index of the non-ghosted starting (lower left)
      ! corner. iend gives the local x-index of the non-ghosted ending 
      ! corner. jstart, jend correspond to y-index, kstart, kend to z-index.
!GEH - Structured Grid Dependence - End

    PetscReal :: radius_0

    ! Grid connections
!GEH - Structured Grid Dependence - Begin
    PetscInt :: nconn, nconnx, nconny
!GEH - Structured Grid Dependence - End

    PetscInt, pointer :: nd1(:), nd2(:)
      ! Nodes upstream and downstream of a connection (assuming flow in 
      ! positive direction.  These are local, ghosted indices.
      
    PetscInt, pointer :: iperm1(:), iperm2(:), ipermbc(:)
    
    PetscReal, pointer :: dist1(:),dist2(:),distbc(:),area(:),areabc(:), grav_ang(:), &
                       delzbc(:), vlbc(:), vvlbc(:),vgbc(:),vvgbc(:)

    PetscReal, pointer :: density_bc(:),d_p_bc(:),d_t_bc(:), d_s_bc(:),d_c_bc(:),&
                       avgmw_bc(:),avgmw_c_bc(:),&
                       hh_bc(:),h_p_bc(:),h_t_bc(:),h_s_bc(:), h_c_bc(:), &
                       viscosity_bc(:),v_p_bc(:),v_t_bc(:),&
                       uu_bc(:),u_p_bc(:),u_t_bc(:),u_s_bc(:), u_c_bc(:),&    
                       df_bc(:),df_p_bc(:),df_t_bc(:),df_s_bc(:), df_c_bc(:), &
                       hen_bc(:),hen_p_bc(:),hen_t_bc(:),hen_s_bc(:),hen_c_bc(:), &
                       pc_bc(:),pc_p_bc(:),pc_t_bc(:),pc_s_bc(:),pc_c_bc(:), &
                       kvr_bc(:),kvr_p_bc(:),kvr_t_bc(:),kvr_s_bc(:),kvr_c_bc(:)
    PetscReal, pointer :: xphi_co2(:),xxphi_co2(:),den_co2(:), dden_co2(:)
    !nL2G :  not collective, local processor: local  =>  ghosted local  
    !nG2L :  not collective, local processor:  ghosted local => local  
    !nG2N :  collective,  ghosted local => global index , used for   
    !                     matsetvaluesblocked ( not matsetvaluesblockedlocal)  
    !nL2A :   collective, local => natural index, used for initialization   
    !                              and source/sink setup  
    PetscInt, pointer :: nL2G(:), nG2L(:), nL2A(:),nG2N(:)
    PetscInt, pointer :: nG2A(:)
      ! Arrays for indexing between local ghosted and non-ghosted, local to natural arrays.
    DA :: da_1_dof, da_nphase_dof, da_3np_dof, da_ndof
    DA :: da_NphaNcomp_dof,da_NphaNspec_dof,da_NphaNspecNcomp_dof
    DA :: da_var_dof
      ! DA's for 1, 3, and multiple (number of phases) degrees of freedom.
      ! da_ndof = total degrees of freedom per node

    ! Boundary conditions (BC's)
    PetscInt :: nblkbc
      ! The number of "blocks" of boundary conditions that are defined.
      ! Such a block is a specification of a set of boundary conditions.
      ! This set of boundary conditions can apply to any number of regions,
      ! so nblkbc does NOT equal the number of boundary condition regions.
    PetscInt :: nconnbc  ! The number of interfaces along boundaries.
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1bc(:), i2bc(:), j1bc(:), j2bc(:), k1bc(:), k2bc(:)
!GEH - Structured Grid Dependence - End
    PetscInt, pointer :: ibconn(:)
      ! ibconn(nc) specifies the id of the of boundary condition block that
      ! applies at boundary interface nc.  
    PetscInt, pointer :: ibndtyp(:)
      ! ibndtyp(ibc) specifies the type of boundary condition that applies
      ! for boundary condition block ibc.
    PetscInt, pointer :: iface(:)
      ! iface(ibc) specifies the face (left, right, top, bottom, etc.) on
      ! which BC block ibc lies.
    PetscInt, pointer :: mblkbc(:)
      ! mblkbc(nc) gives the local, non-ghosted index of the cell that has
      ! boundary connection nc.
    PetscInt, pointer :: iregbc1(:), iregbc2(:)
      ! iregbc1(ibc) and iregbc2(ibc) give the id of the first region and 
      ! last region, respectively, that utilizes the boundary conditions in 
      ! boundary condition block ibc.
    PetscReal, pointer :: pressurebc(:,:)
      ! For a Dirichlet BC, pressurebc(j,ibc) gives the partial pressure 
      ! for phase j along the BC block ibc.
    PetscReal, pointer :: velocitybc(:,:)
      ! For a Neumann BC, velocitybc(j,ibc) gives the velocity q for phase
      ! j along BC block ibc.
    PetscReal, pointer :: tempbc(:),concbc(:),sgbc(:),xphi_co2_bc(:),xxphi_co2_bc(:)
    PetscReal, pointer :: xxbc(:,:), varbc(:)
    PetscInt, pointer:: iphasebc(:)

    !block BC values read from input
    PetscReal, pointer :: pressurebc0(:,:)
    PetscReal, pointer :: velocitybc0(:,:)
    PetscReal, pointer :: tempbc0(:),concbc0(:),sgbc0(:)
    PetscReal, pointer :: xxbc0(:,:)
    PetscInt, pointer:: iphasebc0(:)  

!   phik
    PetscInt :: iregperm, iran_por=0, iread_perm=0, iread_geom =1
    PetscReal :: ran_fac=-1.d0
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1reg(:),i2reg(:),j1reg(:),j2reg(:),k1reg(:),k2reg(:)
!GEH - Structured Grid Dependence - End
    PetscReal, pointer :: por_reg(:),tor_reg(:),perm_reg(:,:)

!   initial conditions
    PetscInt :: iregini
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1ini(:),i2ini(:),j1ini(:),j2ini(:),k1ini(:),k2ini(:)
!GEH - Structured Grid Dependence - End
    PetscReal, pointer :: pres_ini(:),temp_ini(:),conc_ini(:),sat_ini(:), &
                       xmol_ini(:)
    PetscReal, pointer :: xx_ini(:,:)
    PetscInt, pointer:: iphas_ini(:)

!   source term
    PetscInt :: nblksrc = 0, ntimsrc = 0, isrc1 = 2
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1src(:), i2src(:), j1src(:), j2src(:), k1src(:), k2src(:)
!GEH - Structured Grid Dependence - End
    PetscReal, pointer :: timesrc(:,:), tempsrc(:,:), qsrc(:,:), csrc(:,:), hsrc(:,:)
    
!   solid reaction rate
    PetscInt :: ityprxn
    PetscReal :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    PetscReal ::qu_kin, yh2o_in_co2=0.D0
    
!   breakthrough curves
    PetscInt :: ibrkcrv = 0
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1brk(:),i2brk(:),j1brk(:),j2brk(:),k1brk(:),k2brk(:)
!GEH - Structured Grid Dependence - End
    PetscInt, pointer :: ibrktyp(:),ibrkface(:)
    
!   dual continuum
    PetscInt :: idcdm = 0, idcmblk = 0
!GEH - Structured Grid Dependence - Begin
    PetscInt, pointer :: i1dcm(:),i2dcm(:),j1dcm(:),j2dcm(:),k1dcm(:),k2dcm(:)
!GEH - Structured Grid Dependence - End
    PetscReal, pointer :: fracture_aperture(:), matrix_block(:)
    
    PetscInt, pointer :: icap_reg(:),ithrm_reg(:)
    PetscReal :: scale
    PetscReal, pointer :: rock_density(:),cpr(:),dencpr(:),ckdry(:),ckwet(:), &
                       tau(:),cdiff(:),cexp(:)
    PetscReal, pointer :: swir(:),lambda(:),alpha(:),pckrm(:),pcwmax(:),pcbetac(:), &
                       pwrprm(:),sir(:,:)
    PetscInt, pointer:: icaptype(:)
!geh material id
    PetscInt, pointer :: imat(:)
    PetscReal :: m_nacl
    PetscReal :: difaq, delhaq, gravity, fmwh2o= 18.0153D0, fmwa=28.96D0, &
              fmwco2=44.0098D0, eqkair, ret=1.d0, fc=1.d0
    
    PetscInt :: ihydrostatic = 0,ideriv = 1
    PetscReal :: dTdz,beta,tref,pref,conc0
    PetscReal :: hydro_ref_xyz(3)
    
!   table lookup
    PetscInt :: itable=0

!GEH - Structured Grid Dependence - Begin    
    PetscReal, pointer :: dx0(:), dy0(:), dz0(:), rd(:)
!GEH - Structured Grid Dependence - End
    PetscReal, pointer :: x(:), y(:), z(:), delz(:) 
    PetscReal :: x_max, x_min, y_max, y_min, z_max, z_min
    !-------------------------------------------------------------------
    ! Quantities defined at each grid point.
    ! NOTE: I adopt the convention that _loc indicates the local portion
    ! of any global vector.
    !-------------------------------------------------------------------

    ! One degree of freedom: Physical coordinates.
    Vec :: conc
    Vec :: porosity, porosity0, porosity_loc, tor, tor_loc
!GEH - Structured Grid Dependence - Begin
    Vec :: dx, dy, dz, dx_loc, dy_loc, dz_loc  ! Grid spacings
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


    PetscReal, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
    PetscReal, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)

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
   

    PetscReal :: atol, rtol, stol, dtol, inf_tol
    PetscReal it_norm, step_norm
    PetscReal, pointer :: steady_eps(:)
    PetscInt :: idt_switch
    PetscInt :: var_plot_num
    character*16, pointer :: var_plot_nm(:)
    PetscInt, pointer :: var_plot_ind(:)  
      ! Absolute, relative, and "change in norm of solution" tolerances.
    PetscInt :: maxit, maxf
      ! The maximum number of iterations and function evaluations, respectively

 !  Vec :: p_nat, t_nat, c_nat, phis_nat, por_nat, vl_nat, s_nat, x_nat !, perm_nat
    ! Holds contents in natural ordering.
!   Vec :: p_all, t_all, c_all, phis_all, por_all, vl_all, s_all !, perm_all 
    ! Used to hold all values on processor 0.

  end type pflowGrid
end module pflow_gridtype_module
