module Init_module

  implicit none

  private

#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

  public :: PflowInit

contains

! ************************************************************************** !
!
! PflowInit: Initializes a pflow grid object
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine PflowInit(simulation,filename)

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Realization_module
  use Material_module
  use Timestepper_module
  use Field_module
  use Connection_module
  use Coupler_module

  use span_wagner_module
  use MPHASE_module
  use Richards_module
  use pflow_convergence_module
  use pflow_solv_module 
  use pflow_vector_ops_module
  use Utility_module
    
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  type(solver_type), pointer :: solver
  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  ISColoring :: iscoloring

  integer :: mcomp, mphas
  integer :: temp_int
  PetscTruth :: iflag

  real*8, pointer :: phis_p(:)
                       
  ! needed for SNESLineSearchGetParams()/SNESLineSearchSetParams()
  real*8 :: alpha, maxstep, steptol
  
  PetscErrorCode :: ierr
  
  ! set pointers to objects
  solver => simulation%stepper%solver
  realization => simulation%realization
  option => realization%option
  field => realization%field
  
  ! read MODE,GRID,PROC,COMP,PHAS cards
  call readSelectCardsFromInput(realization,filename,mcomp,mphas)
  grid => realization%grid

  ! set the operational mode (e.g. RICHARDS_MODE, MPH_MODE, etc)
  call setMode(option,mcomp,mphas)
  ! process command line options
  call OptionCheckCommandLine(option)
                             
! hardwire to uncoupled for now
!  if (icouple == 0) then
  option%run_coupled = PETSC_FALSE
!  else
!    option%run_coupled = PETSC_TRUE
!  endif

  !set specific phase indices
  option%jh2o = 1; option%jgas =1
  select case(option%nphase)
    case(2)
      if (option%imode == TWOPH_MODE) then
        option%jgas=  2
        option%jco2 = 3
      else
        option%jco2 = 2
        option%jgas =3 
      endif
    case(3)
      option%jco2 = 2
      option%jgas =3 
  end select 

  call GridCreateDMs(grid,option)
  
 !-----------------------------------------------------------------------
 ! Create the vectors with parallel layout corresponding to the DM's,
 ! and, for vectors that need to be ghosted, create the corresponding
 ! ghosted vectors.
 !-----------------------------------------------------------------------

  ! 1 degree of freedom
  call GridCreateVector(grid,ONEDOF,field%porosity0,GLOBAL)
  call VecDuplicate(field%porosity0, grid%volume, ierr)
  call VecDuplicate(field%porosity0, field%phis, ierr)
  call VecDuplicate(field%porosity0, field%perm0_xx, ierr)
  call VecDuplicate(field%porosity0, field%perm0_yy, ierr)
  call VecDuplicate(field%porosity0, field%perm0_zz, ierr)
  call VecDuplicate(field%porosity0, field%perm_pow, ierr)
      
  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)
      call VecDuplicate(field%porosity0, field%conc, ierr)
      call VecDuplicate(field%porosity0, field%temp, ierr)
      call VecDuplicate(field%porosity0, field%ttemp, ierr)
  end select    
      
  call GridCreateVector(grid,ONEDOF,field%porosity_loc,LOCAL)
  call VecDuplicate(field%porosity_loc, field%tor_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%ithrm_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%icap_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%iphas_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%iphas_old_loc, ierr)
  
  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)  
      call VecDuplicate(field%porosity_loc, field%ttemp_loc, ierr)
  end select

  call VecDuplicate(field%porosity_loc, field%perm_xx_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%perm_yy_loc, ierr)
  call VecDuplicate(field%porosity_loc, field%perm_zz_loc, ierr)

  if (associated(grid%structured_grid)) then
    call VecDuplicate(field%porosity0, grid%structured_grid%dx, ierr)
    call VecDuplicate(field%porosity0, grid%structured_grid%dy, ierr)
    call VecDuplicate(field%porosity0, grid%structured_grid%dz, ierr)

    call VecDuplicate(field%porosity_loc, grid%structured_grid%dx_loc, ierr)
    call VecDuplicate(field%porosity_loc, grid%structured_grid%dy_loc, ierr)
    call VecDuplicate(field%porosity_loc, grid%structured_grid%dz_loc, ierr)
  endif

  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)
      ! nphase degrees of freedom
      call GridCreateVector(grid,NPHASEDOF,field%pressure,GLOBAL)
      call VecDuplicate(field%pressure, field%sat, ierr)
      call VecDuplicate(field%pressure, field%xmol, ierr)
      call VecDuplicate(field%pressure, field%ppressure, ierr)
      call VecDuplicate(field%pressure, field%ssat, ierr)
      call VecDuplicate(field%pressure, field%dp, ierr)
      call VecDuplicate(field%pressure, field%density, ierr)
      call VecDuplicate(field%pressure, field%ddensity, ierr)
      call VecDuplicate(field%pressure, field%avgmw, ierr)
      call VecDuplicate(field%pressure, field%d_p, ierr)
      call VecDuplicate(field%pressure, field%d_t, ierr)
      call VecDuplicate(field%pressure, field%h, ierr)
      call VecDuplicate(field%pressure, field%hh, ierr)
      call VecDuplicate(field%pressure, field%h_p, ierr)
      call VecDuplicate(field%pressure, field%h_t, ierr)
      call VecDuplicate(field%pressure, field%viscosity, ierr)
      call VecDuplicate(field%pressure, field%v_p, ierr)
      call VecDuplicate(field%pressure, field%v_t, ierr)
     ! xmol may not be nphase DOF, need change later 
      call VecDuplicate(field%pressure, field%xxmol, ierr)
  end select

  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)
      call GridCreateVector(grid,NPHASEDOF, field%ppressure_loc, LOCAL)
      call VecDuplicate(field%ppressure_loc, field%ssat_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%xxmol_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%ddensity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%avgmw_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%hh_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%viscosity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_t_loc, ierr)
  end select
  
  ! 3 * nphase degrees of freedom (velocity vector)
  call GridCreateVector(grid,THREENPDOF, field%vl, GLOBAL)
      
  if (option%run_coupled == PETSC_TRUE) &
    call VecDuplicate(field%vl, field%vvl, ierr)


  if (option%imode == TWOPH_MODE) then
    print *,'2ph add var'
    call VecDuplicate(field%sat, field%d_s, ierr)
    call VecDuplicate(field%sat, field%h_s, ierr)
    call VecDuplicate(field%sat, field%u, ierr)
    call VecDuplicate(field%sat, field%uu, ierr)
    call VecDuplicate(field%sat, field%u_s, ierr)
    call VecDuplicate(field%sat, field%u_p, ierr)
    call VecDuplicate(field%sat, field%u_t, ierr)

    call VecDuplicate(field%ssat_loc, field%d_s_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%h_s_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%u_s_loc, ierr)

    call VecDuplicate(field%sat, field%pcw, ierr)
    call VecDuplicate(field%sat, field%pc_p, ierr)
    call VecDuplicate(field%sat, field%pc_t, ierr)
    call VecDuplicate(field%sat, field%pc_s, ierr)
    call VecDuplicate(field%ssat_loc, field%pcw_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%pc_p_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%pc_t_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%pc_s_loc, ierr)


    call VecDuplicate(field%sat, field%kvr, ierr)
    call VecDuplicate(field%sat, field%kvr_p, ierr)
    call VecDuplicate(field%sat, field%kvr_t, ierr)
    call VecDuplicate(field%sat, field%kvr_s, ierr)
    call VecDuplicate(field%ssat_loc, field%kvr_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%kvr_p_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%kvr_t_loc, ierr)
    call VecDuplicate(field%ssat_loc, field%kvr_s_loc, ierr)

    call GridCreateVector(grid,NPHANCOMPDOF, field%h_c,GLOBAL)
    call VecDuplicate(field%h_c, field%u_c, ierr)
    call VecDuplicate(field%h_c, field%avgmw_c, ierr)
    call VecDuplicate(field%h_c, field%d_c, ierr)
    call VecDuplicate(field%h_c, field%pc_c, ierr)
    call VecDuplicate(field%h_c, field%kvr_c, ierr)
        
    call GridCreateVector(grid,NPHANCOMPDOF, field%h_c_loc,LOCAL)
    call VecDuplicate(field%h_c_loc, field%avgmw_c_loc, ierr)
    call VecDuplicate(field%h_c_loc, field%d_c_loc, ierr)
    call VecDuplicate(field%h_c_loc, field%pc_c_loc, ierr)
    call VecDuplicate(field%h_c_loc, field%kvr_c_loc, ierr)


    call GridCreateVector(grid,NPHANSPECDOF, field%hen,GLOBAL)
    call GridCreateVector(grid,NPHANSPECDOF, field%hen_loc,LOCAL)
    call VecDuplicate(field%hen, field%hen_p, ierr)
    call VecDuplicate(field%hen_loc, field%hen_p_loc, ierr)
    call VecDuplicate(field%hen, field%hen_t, ierr)
    call VecDuplicate(field%hen_loc, field%hen_t_loc, ierr)
    call VecDuplicate(field%hen, field%hen_s, ierr)
    call VecDuplicate(field%hen_loc, field%hen_s_loc, ierr)
    call VecDuplicate(field%hen, field%df, ierr)
    call VecDuplicate(field%hen_loc, field%df_loc, ierr)
    call VecDuplicate(field%df, field%df_p, ierr)
    call VecDuplicate(field%df_loc, field%df_p_loc, ierr)
    call VecDuplicate(field%df, field%df_t, ierr)
    call VecDuplicate(field%df_loc, field%df_t_loc, ierr)
    call VecDuplicate(field%df, field%df_s, ierr)
    call VecDuplicate(field%df_loc, field%df_s_loc, ierr)
  
    call GridCreateVector(grid,NPHANSPECNCOMPDOF,field%hen_c, GLOBAL)
    call GridCreateVector(grid,NPHANSPECNCOMPDOF,field%hen_c_loc,LOCAL)
      
    call VecDuplicate(field%hen_c, field%df_c, ierr)
    call VecDuplicate(field%hen_c_loc, field%df_c_loc, ierr)
  end if  

  select case(option%imode)
    case(MPH_MODE,RICHARDS_MODE,FLASH_MODE,OWG_MODE,VADOSE_MODE)
      call GridCreateVector(grid,VARDOF, field%var_loc,LOCAL)
  end select

      ! ndof degrees of freedom
  call GridCreateVector(grid,NDOF, field%xx, GLOBAL)
  call VecDuplicate(field%xx, field%yy, ierr)
  call VecDuplicate(field%xx, field%dxx, ierr)
  call VecDuplicate(field%xx, field%r, ierr)
  call VecDuplicate(field%xx, field%accum, ierr)
     
  call VecSetBlocksize(field%dxx, option%ndof, ierr)

  call GridCreateVector(grid,NDOF, field%xx_loc, LOCAL)
  
!-----------------------------------------------------------------------
! Set up PETSc nonlinear solver context.
!-----------------------------------------------------------------------
  call SNESCreate(PETSC_COMM_WORLD, solver%snes, ierr)
  CHKERRQ(ierr)
!-----------------------------------------------------------------------
  ! Set up indexing of grid ids (local to global, global to local, etc
!-----------------------------------------------------------------------
  ! set up nG2L, NL2G, etc.
  call GridMapIndices(grid)

  option%nldof = grid%nlmax * option%nphase
  option%ngdof = grid%ngmax * option%nphase
  
!-----------------------------------------------------------------------
      ! Allocate memory for allocatable arrays.
!-----------------------------------------------------------------------
  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)
      allocate(field%xphi_co2(grid%nlmax))
      allocate(field%xxphi_co2(grid%nlmax))
      allocate(field%den_co2(grid%nlmax))
      allocate(field%dden_co2(grid%nlmax))
      field%xphi_co2 = 1.d0
      field%xxphi_co2 = 1.d0
      field%den_co2 = 1.d0
      field%dden_co2 = 1.d0
  end select
  
  !set scale factor for heat equation, i.e. use units of MJ for energy
  option%scale = 1.d-6

  if (option%run_coupled == PETSC_TRUE) then
    ! necessary for water balance due to generation/consumption of H20 
    ! by chemical reactions during coupled pflow/ptran runs
    allocate(option%rtot(grid%nlmax,2))
    option%rtot=0.D0
  endif

  ! read in the remainder of the input file
  call readInput(simulation,filename)

  select case(option%imode)
    case(MPH_MODE,OWG_MODE,FLASH_MODE)
      call initialize_span_wagner(option%itable,option%myrank)  
  end select
  
  call GridComputeSpacing(grid)
  call GridComputeCoordinates(grid,option)
  call GridComputeVolumes(grid,option)

  ! set up internal connectivity, distance, etc.
  call GridComputeInternalConnect(grid,option)
  call RealizationProcessCouplers(realization)
  
  ! clip regions and set up boundary connectivity, distance  
  call GridLocalizeRegions(realization%regions,realization%grid,realization%option)

  ! connectivity between initial conditions, boundary conditions, srcs/sinks, etc and grid
  call GridComputeCouplerConnections(grid,option,realization%initial_conditions%first)
  call GridComputeCouplerConnections(grid,option,realization%boundary_conditions%first)
  call GridComputeCouplerConnections(grid,option,realization%source_sinks%first)
                                
  call assignMaterialPropToRegions(realization)
  call RealizationInitCouplerAuxVars(realization,realization%initial_conditions)
  call RealizationInitCouplerAuxVars(realization,realization%boundary_conditions)
  call assignInitialConditions(realization)

  allocate(realization%field%internal_velocities(option%nphase, &
             ConnectionGetNumberInList(realization%grid%internal_connection_list)))
  temp_int = CouplerGetNumConnectionsInList(realization%boundary_conditions)
  allocate(realization%field%boundary_velocities(option%nphase,temp_int))           

  allocate(field%xphi_co2_bc(temp_int))
  allocate(field%xxphi_co2_bc(temp_int))

  select case(option%imode)
    ! everything but RICHARDS_MODE for now
    case(MPH_MODE,COND_MODE,TWOPH_MODE,VADOSE_MODE,LIQUID_MODE,OWG_MODE, &
         FLASH_MODE,TH_MODE,THC_MODE)
      temp_int = ConnectionGetNumberInList(grid%internal_connection_list)
      allocate(field%vl_loc(temp_int))
      allocate(field%vvl_loc(temp_int))
      allocate(field%vg_loc(temp_int))
      allocate(field%vvg_loc(temp_int))
      field%vl_loc = 0.D0
      field%vvl_loc = 0.D0
      field%vg_loc = 0.D0
      field%vvg_loc = 0.D0
  end select
    
! check number of dofs and phases
  iflag = PETSC_FALSE
  select case(option%imode)
    case(COND_MODE)
      if (option%ndof .ne. 1 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
    case(TH_MODE)
      if (option%ndof .ne. 2 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
    case(THC_MODE)
      if (option%ndof .ne. 3 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
    case(TWOPH_MODE)
      if (option%ndof .ne. 4 .or. option%nphase .ne. 2) iflag = PETSC_TRUE
    case(MPH_MODE,RICHARDS_MODE,FLASH_MODE,VADOSE_MODE)
      if (option%ndof .ne. (option%nspec+1)) iflag = PETSC_TRUE
    case(OWG_MODE)
      if (option%ndof .ne. (option%nspec)) iflag = PETSC_TRUE
    case default
      if (option%ndof .ne. 1 .or. option%nphase .ne. 1) iflag = PETSC_TRUE
  end select
  
  if (iflag == PETSC_TRUE) then
    write(*,*) 'Specified number of dofs or phases not correct-stop: ', &
               trim(option%mode), 'ndof= ',option%ndof,' nph= ', &
               option%nphase
    stop
  endif


  if (option%myrank == 0) then
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    if (grid%igrid == STRUCTURED) then
      write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
        option%commsize,grid%structured_grid%npx,grid%structured_grid%npy, &
        grid%structured_grid%npz
    endif
    write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
      option%ndof,option%nphase
    select case(option%imode)
      case(COND_MODE)
        write(*,'(" mode = Conduction: T")')
      case(TH_MODE)
        write(*,'(" mode = TH: p, T")')
      case(THC_MODE)
        write(*,'(" mode = THC: p, T, C")')
      case(TWOPH_MODE)
        write(*,'(" mode = 2-PH: p, T, s, C")')
      case(MPH_MODE)
        write(*,'(" mode = MPH: p, T, s/C")')
      case(FLASH_MODE)
        write(*,'(" mode = flash: p, T, z")')
      case(VADOSE_MODE)
        write(*,'(" mode = VAD: p, T, s/C")')
      case(RICHARDS_MODE)
        write(*,'(" mode = Richards: p, T, s/C")')
      case(OWG_MODE)
        write(*,'(" mode = O+W+G: p, T, s/C")')
      case default
        write(*,'(" mode = Single Liquid Phase: p")')
    end select
  endif

  !-----------------------------------------------------------------------
  ! Set up the Jacobian matrix.  We do this here instead of in 
  ! pflowgrid_new() because we may have to parse the input file to 
  ! determine how we want to do the Jacobian (matrix vs. matrix-free, for
  ! example).
  !-----------------------------------------------------------------------
  if (option%use_analytical == PETSC_TRUE) then
  
    option%ideriv = 1
  
    call GridCreateJacobian(grid,solver,option)
  
!   if (myrank == 0) write(*,'(" analytical jacobian as ")'); &
!                    print *, grid%iblkfmt

    select case(option%imode)
#if 0    
      ! still needs to be implemented
      case(COND_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, CondJacobian, &
                           grid, ierr); CHKERRQ(ierr)
      case(TH_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, THJacobian, &
                           grid, ierr); CHKERRQ(ierr)
      case(THC_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, THCJacobian, &
                           grid, ierr); CHKERRQ(ierr)
      case(TWOPH_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, TTPHASEJacobian, &
                           grid, ierr); CHKERRQ(ierr)
      case(FLASH_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, FlashJacobian, &
                           grid, ierr); CHKERRQ(ierr)
        if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
      case(VADOSE_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, VADOSEJacobian, &
                           grid, ierr); CHKERRQ(ierr)
        if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
#endif
      case(RICHARDS_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, RichardsJacobian, &
                             realization, ierr); CHKERRQ(ierr)
        if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(realization,solver)
      case(MPH_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, MPHASEJacobian, &
                             realization, ierr); CHKERRQ(ierr)
        if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(realization,solver)
#if 0      
      ! still needs to be implemented
      case(OWG_MODE)
        call SNESSetJacobian(solver%snes, solver%J, solver%J, OWGJacobian, &
                           grid, ierr); CHKERRQ(ierr)
        if (option%use_ksp == PETSC_TRUE) call pflow_kspsolver_init(grid)
      case default
        call SNESSetJacobian(solver%snes, solver%J, solver%J, LiquidJacobian, &
                           grid, ierr); CHKERRQ(ierr)
#endif
    end select                         

  else if (option%use_matrix_free == PETSC_TRUE) then
  
    option%ideriv = 0
  
    if (option%myrank == 0) write(*,'(" Using matrix-free Newton-Krylov")')

    select case(option%imode)
      case(COND_MODE)
        call MatCreateMFFD(solver%snes,field%ttemp,solver%J,ierr)
      case(TH_MODE,THC_MODE,TWOPH_MODE,MPH_MODE,FLASH_MODE,RICHARDS_MODE, &
           VADOSE_MODE,OWG_MODE)
        call MatCreateMFFD(solver%snes,field%xx,solver%J,ierr)
      case default
        call MatCreateMFFD(solver%snes,field%ppressure,solver%J,ierr)
    end select
    
    ! It seems that I ought to call SNESSetJacobian here now, but I don't know
    ! what function I am supposed to pass to it to get it to use one of the 
    ! finite-difference routines for computing the Jacobian.  It might not
    ! actually matter if -snes_mf has been specified.
    ! Pernice thinks that perhaps the I need to provide a function which 
    ! simply calls MatAssemblyBegin/End.
    call SNESSetJacobian(solver%snes, solver%J, solver%J, &
                         SolverComputeMFJacobian, PETSC_NULL_OBJECT, ierr)

    ! Use "Walker-Pernice" differencing.
    call MatMFFDSetType(solver%J, MATMFFD_WP, ierr)

    if (option%print_hhistory == PETSC_TRUE) then
      allocate(option%hhistory(HHISTORY_LENGTH))
      call MatMFFDSetHHistory(solver%J, option%hhistory, HHISTORY_LENGTH, ierr)
    endif

    if (option%monitor_h == PETSC_TRUE) then
      call SNESMonitorSet(solver%snes, SolverMonitorH, grid, &
                          PETSC_NULL_OBJECT, ierr)
    endif
    
  else
  
    option%ideriv = 0
  
    if (option%myrank == 0) write(*,'(" numerical jacobian")')
    ! We will compute the Jacobian via finite differences and store it.
    
    ! Create matrix with correct parallel layout and nonzero structure to 
    ! hold the Jacobian.
      ! MatFDColoringCreate() currently does not support MATMPIBAIJ.
      
    temp_int = option%iblkfmt
    option%iblkfmt = 0   ! to turn off MATMPIBAIJ
    call GridCreateJacobian(grid,solver,option)
    option%iblkfmt = temp_int

    call GridCreateColoring(grid,option,iscoloring)
        
    call MatFDColoringCreate(solver%J, iscoloring,solver%matfdcoloring, ierr)
    
    call ISColoringDestroy(iscoloring, ierr)

    select case(option%imode)
#if 0
    ! need to be implemented
      case(COND_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          CondResidual,option,ierr)
      case(TH_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          THResidual,option, ierr)
      case(THC_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          THCResidual,option, ierr)
      case(TWOPH_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          TTPHASEResidual,option, ierr)
      case(FLASH_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          FlashResidual,option, ierr)
      case(VADOSE_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          VADOSEResidual,option, ierr)
#endif                                          
      case(RICHARDS_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          RichardsResidual,option, ierr)
      case(MPH_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          MPHASEResidual,option, ierr)
#if 0                                          
    ! need to be implemented
      case(OWG_MODE)
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          OWGResidual,option, ierr)
      case default
        call MatFDColoringSetFunctionSNES(solver%matfdcoloring, &
                                          LiquidResidual,option, ierr)
#endif                                          
    end select
        
    call MatFDColoringSetFromOptions(solver%matfdcoloring, ierr)
    
    call SNESSetJacobian(solver%snes, solver%J, solver%J, &
                         SNESDefaultComputeJacobianColor,  &
                         solver%matfdcoloring, ierr)
  endif

  if (option%myrank == 0) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')

  select case(option%imode)
#if 0  
    ! need to be implemented
    case(COND_MODE)
      call SNESSetFunction(solver%snes,field%r,CondResidual,realization,ierr)
    case(TH_MODE)
      call SNESSetFunction(solver%snes,field%r,THResidual,realization,ierr)
    case(THC_MODE)
      call SNESSetFunction(solver%snes,field%r,THCResidual,realization,ierr)
    case(TWOPH_MODE)
      call SNESSetFunction(solver%snes,field%r,TTPHASEResidual,realization,ierr)
    case(FLASH_MODE)
      call SNESSetFunction(solver%snes,field%r,FlashResidual,realization,ierr)
    case(VADOSE_MODE)
      call SNESSetFunction(solver%snes,field%r,VADOSEResidual,realization,ierr)
#endif      
    case(RICHARDS_MODE)
      call SNESSetFunction(solver%snes,field%r,RichardsResidual,realization,ierr)
    case(MPH_MODE)
      call SNESSetFunction(solver%snes,field%r,MPHASEResidual,realization,ierr)
#if 0 
    ! need to be implemented     
    case(OWG_MODE)
      call SNESSetFunction(solver%snes,field%r,OWGResidual,realization,ierr)
    case default
      call SNESSetFunction(solver%snes,field%r,LiquidResidual,realization,ierr)
#endif      
  end select

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(solver%snes, solver%atol, solver%rtol, solver%stol, & 
                         solver%maxit, solver%maxf, ierr)
  call SNESSetFromOptions(solver%snes, ierr)
  
  ! shell for custom convergence test.  The default SNES convergence test 
  ! is call within this function.
  call SNESSetConvergenceTest(solver%snes,PFLOWConvergenceTest,simulation,ierr)
                           
  call SNESLineSearchGetParams(solver%snes, alpha, maxstep, steptol, ierr) 
  call SNESLineSearchSetParams(solver%snes, alpha, maxstep, solver%stol, ierr) 
  call SNESGetKSP(solver%snes, solver%ksp, ierr)
  call KSPSetTolerances(solver%ksp,solver%rtol,solver%atol,solver%dtol, &
                        10000,ierr)

  if (option%myrank == 0) write(*,'("  Finished setting up of SNES ")')
 

  select case(option%imode)
    case(MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE)
!      allocate(field%varbc(1:(option%ndof+1)*(2+7*option%nphase + 2 *  &
!                                       option%nphase*option%nspec)))
    case default  
      allocate(field%density_bc(option%nphase))
      allocate(field%d_p_bc(option%nphase))
      allocate(field%d_t_bc(option%nphase))
      allocate(field%d_c_bc(option%nphase))
      allocate(field%d_s_bc(option%nphase))
      allocate(field%avgmw_bc(option%nphase))
      allocate(field%avgmw_c_bc(option%nphase*option%npricomp))
      allocate(field%hh_bc(option%nphase))
      allocate(field%h_p_bc(option%nphase))
      allocate(field%h_t_bc(option%nphase))
      allocate(field%h_c_bc(option%nphase*option%npricomp))
      allocate(field%h_s_bc(option%nphase))
      allocate(field%uu_bc(option%nphase))
      allocate(field%u_p_bc(option%nphase))
      allocate(field%u_t_bc(option%nphase))
      allocate(field%u_c_bc(option%nphase*option%npricomp))
      allocate(field%u_s_bc(option%nphase))
      allocate(field%df_bc(option%nphase*option%nspec))
      allocate(field%df_p_bc(option%nphase*option%nspec))
      allocate(field%df_t_bc(option%nphase*option%nspec))
      allocate(field%df_c_bc(option%nphase*option%nspec*option%npricomp))
      allocate(field%df_s_bc(option%nphase*option%nspec))
      allocate(field%hen_bc(option%nphase*option%nspec))
      allocate(field%hen_p_bc(option%nphase*option%nspec))
      allocate(field%hen_t_bc(option%nphase*option%nspec))
      allocate(field%hen_c_bc(option%nphase*option%nspec*option%npricomp))
      allocate(field%hen_s_bc(option%nphase*option%nspec))
      allocate(field%viscosity_bc(option%nphase))
      allocate(field%v_p_bc(option%nphase))
      allocate(field%v_t_bc(option%nphase))
      allocate(field%pc_bc(option%nphase))
      allocate(field%pc_p_bc(option%nphase))
      allocate(field%pc_t_bc(option%nphase))
      allocate(field%pc_c_bc(option%nphase*option%npricomp))
      allocate(field%pc_s_bc(option%nphase))
      allocate(field%kvr_bc(option%nphase))
      allocate(field%kvr_p_bc(option%nphase))
      allocate(field%kvr_t_bc(option%nphase))
      allocate(field%kvr_c_bc(option%nphase*option%npricomp))
      allocate(field%kvr_s_bc(option%nphase))
  end select
   

  if(option%print_bcinfo == PETSC_TRUE) then
    if (option%myrank == 0) print *, 'Someone will have to rewrite print_bcinfo'
  endif

#if 0
  ! should we still support this
  if (grid%iread_perm == 1) then
    call Read_perm_field(grid)
  endif
#endif
!geh
  nullify(field%imat)

#if 0
  ! should we still support this
  if (grid%iread_geom == 10) then 
    if (myrank == 0) print *, 'Reading structured grid from hdf5' 
    allocate(grid%imat(grid%ngmax))  ! allocate material id array
    call ReadStructuredGridHDF5(grid)
    
  else if (grid%iread_geom == -1) then 
    if (myrank == 0) print *, 'Reading unstructured grid' 
    allocate(grid%pressurebc(grid%nphase,grid%nconnbc)) 

    allocate(grid%imat(grid%ngmax))  ! allocate material id array
    grid%imat = 0      

    call ReadUnstructuredGrid(grid) 
! this call is to set up an array for zeroing inactive and isothermal cells
    call createRichardsZeroArray(grid)
    
    ! dump material ids to file in natural ordering
    call DACreateGlobalVector(grid%da_1_dof,temp_vec,ierr)
    call VecGetArrayF90(temp_vec,temp_p,ierr)
    do i=1, grid%nlmax
      temp_p(i) = grid%imat(grid%nL2G(i))*1.d0      
    enddo
    call VecRestoreArrayF90(temp_vec,temp_p,ierr)
    call VecDestroy(temp_vec,ierr)
  endif 
#endif   

  if (option%ihydrostatic == 3) then
    if (option%imode == MPH_MODE .or. &
        option%imode == VADOSE_MODE .or. &
        option%imode == FLASH_MODE .or. &
        option%imode == RICHARDS_MODE) then
        ! print *,'in nhydro'
      print *, 'fix nhydrostatic3'
      stop
#if 0
      ! needs to be reimplemented       
      call nhydrostatic3(grid)
      ! print *,'out nhydro'
#endif        
    endif
  endif
 
! Note: VecAssemblyBegin/End needed to run on the Mac - pcl (11/21/03)!
  if (field%conc /= 0) then
    call VecAssemblyBegin(field%conc,ierr)
    call VecAssemblyEnd(field%conc,ierr)
  endif
  if (field%xmol /= 0) then
    call VecAssemblyBegin(field%xmol,ierr)
    call VecAssemblyEnd(field%xmol,ierr)
  endif

  if (option%myrank == 0) write(*,'("  Finished setting up of INIT ")')

  !-----------------------------------------------------------------------
  ! Initialize field variables
  !-----------------------------------------------------------------------  
  if (option%imode /= MPH_MODE .and. option%imode /= OWG_MODE .and. &
      option%imode /= VADOSE_MODE .and. option%imode /= FLASH_MODE .and. &
      option%imode /= RICHARDS_MODE) then   

    select case(option%ndof)
      case(1)
        call VecCopy(field%pressure, field%ppressure, ierr)
        call VecCopy(field%temp, field%ttemp, ierr)
      case(2)
        call pflow_pack_xx2(field%yy, field%pressure, option%nphase, field%temp, 1, &
                            ierr)
        call VecCopy(field%yy, field%xx, ierr)     
      case(3)
        call pflow_pack_xx3(field%yy, field%pressure, option%nphase, field%temp, 1, &
                            field%conc, 1, ierr)
        call VecCopy(field%yy, field%xx, ierr)      
    end select
  endif
      
         
  select case(option%imode)
#if 0
    ! still need implementation
    case(TWOPH_MODE)
      call pflow_2phase_initadj(grid)
      call VecCopy(grid%iphas, grid%iphas_old,ierr)
      call pflow_pack_xx4(grid%yy, grid%pressure, grid%nphase, grid%temp, 1, &
           grid%xmol,grid%nphase , grid%sat,grid%nphase , ierr)
      call VecCopy(grid%yy, grid%xx, ierr)
      call pflow_update_2phase(grid)
#endif
    case(RICHARDS_MODE)
      call pflow_richards_initadj(realization)
      call VecCopy(field%iphas_loc, field%iphas_old_loc,ierr)
      call VecCopy(field%xx, field%yy, ierr)
      call pflow_update_richards(realization)
    case(MPH_MODE)
      call pflow_mphase_initadj(realization)
      call VecCopy(field%iphas_loc, field%iphas_old_loc,ierr)
      call VecCopy(field%xx, field%yy, ierr)
      call pflow_update_mphase(realization)
#if 0
    ! still need implementation
    case(FLASH_MODE)
      call pflow_flash_initadj(grid)
      call VecCopy(grid%iphas, grid%iphas_old,ierr)
      call VecCopy(grid%xx, grid%yy, ierr)
      call pflow_update_flash(grid)
    case(OWG_MODE)
      call pflow_owg_initadj(grid)
      call VecCopy(grid%iphas, grid%iphas_old,ierr)
      call VecCopy(grid%xx, grid%yy, ierr)
      call pflow_update_owg(grid)
    case(VADOSE_MODE)
      call pflow_vadose_initadj(grid)
      call VecCopy(grid%iphas, grid%iphas_old,ierr)
      call VecCopy(grid%xx, grid%yy, ierr)
      call pflow_update_vadose(grid)
#endif
  end select  
  
! zero initial velocity
  call VecSet(field%vl,0.d0,ierr)
  if (option%run_coupled == PETSC_TRUE) call VecSet(field%vvl,0.d0,ierr)
 
  if (option%myrank == 0) &
    write(*,'("  Finished setting up of INIT2 ")')
   
  call initializeSolidReaction(realization)

  if (option%myrank == 0) write(*,'("  Finished setting up ")')
  
end subroutine PflowInit

! ************************************************************************** !
!
! readSelectCardsFromInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine readSelectCardsFromInput(realization,filename,mcomp,mphas)

  use Option_module
  use Grid_module
  use Fileio_module
  use Realization_module
  
  implicit none

  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename
  integer :: mcomp, mphas

  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  integer :: igeom
  
  option => realization%option
  
  open(IUNIT1, file=filename, action="read", status="old") 
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  mphas=0
  mcomp = 0

! Read in select required cards
!.........................................................................

  ! MODE information
  string = "MODE"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! read in keyword 
  call fiReadWord(string,option%mode,.true.,ierr)
  call fiErrorMsg('mode','mode',ierr)

!.........................................................................

  ! GRID information
  string = "GRID"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! key off igeom for structured vs unstructured 
  call fiReadInt(string,igeom,ierr)
  call fiDefaultMsg('igeom',ierr)
  
  realization%grid => GridCreate(igeom) 
  grid => realization%grid

  if (grid%igrid == STRUCTURED) then ! structured
    call fiReadInt(string,grid%structured_grid%nx,ierr)
    call fiDefaultMsg('nx',ierr)
    
    call fiReadInt(string,grid%structured_grid%ny,ierr)
    call fiDefaultMsg('ny',ierr)
    
    call fiReadInt(string,grid%structured_grid%nz,ierr)
    call fiDefaultMsg('nz',ierr)
    
    grid%structured_grid%nxy = grid%structured_grid%nx*grid%structured_grid%ny
    grid%structured_grid%nmax = grid%structured_grid%nxy * &
                                grid%structured_grid%nz
    grid%nmax = grid%structured_grid%nmax
  else ! unstructured
  endif
      
  call fiReadInt(string,option%nphase,ierr)
  call fiDefaultMsg('nphase',ierr)

  call fiReadInt(string,option%nspec,ierr)
  call fiDefaultMsg('nspec',ierr)

  call fiReadInt(string,option%npricomp,ierr)
  call fiDefaultMsg('npricomp',ierr)

  call fiReadInt(string,option%ndof,ierr)
  call fiDefaultMsg('ndof',ierr)
      
  call fiReadInt(string,option%idcdm,ierr)
  call fiDefaultMsg('idcdm',ierr)

  call fiReadInt(string,option%itable,ierr)
  call fiDefaultMsg('itable',ierr)

  if (option%myrank==0) then
    if (grid%igrid == STRUCTURED) then
      write(IUNIT2,'(/," *GRID ",/, &
        &"  igeom   = ",i4,/, &
        &"  nx      = ",i4,/, &
        &"  ny      = ",i4,/, &
        &"  nz      = ",i4,/, &
        &"  nphase  = ",i4,/, &
        &"  ndof    = ",i4,/, &
        &"  idcdm   = ",i4,/, &
        &"  itable  = ",i4    &
        &   )') grid%igeom, grid%structured_grid%nx, &
                grid%structured_grid%ny, grid%structured_grid%nz, &
                option%nphase, option%ndof, option%idcdm, option%itable
    else
    endif
  endif

!.........................................................................

  if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
    
    ! PROC information
    string = "PROC"
    call fiFindStringInFile(IUNIT1,string,ierr)

    if (ierr == 0) then

      ! strip card from front of string
      call fiReadWord(string,word,.false.,ierr)
      call fiReadInt(string,grid%structured_grid%npx,ierr)
      call fiDefaultMsg('npx',ierr)
      call fiReadInt(string,grid%structured_grid%npy,ierr)
      call fiDefaultMsg('npy',ierr)
      call fiReadInt(string,grid%structured_grid%npz,ierr)
      call fiDefaultMsg('npz',ierr)
 
      if (option%myrank == 0) &
        write(IUNIT2,'(/," *PROC",/, &
          & "  npx   = ",3x,i4,/, &
          & "  npy   = ",3x,i4,/, &
          & "  npz   = ",3x,i4)') grid%structured_grid%npx, &
            grid%structured_grid%npy, grid%structured_grid%npz
  
      if (option%commsize /= grid%structured_grid%npx * &
                             grid%structured_grid%npy * &
                             grid%structured_grid%npz) then
        if (option%myrank==0) &
          write(*,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%npx*grid%structured_grid%npy* &
                       grid%structured_grid%npz,' commsize = ',option%commsize
        stop
      endif
    endif
  endif
  
!.........................................................................

  ! COMP information
  string = "COMP"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr == 0) then

    mcomp = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('COMP',ierr)
      
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','GAS',ierr)
        
      call fiWordToUpper(name) 
      select case(name(1:len_trim(name)))
        case('H2O')
            mcomp = mcomp +1
        case('TRACER')     
            mcomp = mcomp +4
        case('CO2')
            mcomp = mcomp +2   
        case('C10H22')
            mcomp = mcomp +8
        case('AIR')              
            mcomp = mcomp +16
        case('ENG')
            mcomp = mcomp +32
        case default               
           print *,' Wrong comp name::  Stop:' 
           stop
      end select 
    enddo
  endif          

!....................................................................

  ! COMP information
  string = "PHAS"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr == 0) then

    mphas = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('phase',ierr)
    
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','phase',ierr)
        
      call fiWordToUpper(name) 
         
      select case(name(1:len_trim(name)))
        case('ROCK')
          mphas = mphas +1
        case('H2O')     
          mphas = mphas +2
        case('CO2')
          mphas = mphas +4   
        case('GAS')
          mphas = mphas +8
        case('NAPL1')              
          mphas = mphas +16
        case('NAPL2')
          mphas = mphas +32
        case default               
          print *,' Wrong phase name::  Stop:' 
          stop
      end select 
    enddo
  endif
    
end subroutine readSelectCardsFromInput

! ************************************************************************** !
!
! readInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readInput(simulation,filename)

  use Simulation_module
  use Option_module
  use Field_module
  use Grid_module
  use Structured_Grid_module
  use Solver_module
  use Material_module
  use Fileio_module
  use Realization_module
  use Timestepper_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module
  use Waypoint_module

  use pckr_module 
  
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
    
  real*8, parameter:: fmwnacl = 58.44277D0, fmwh2o  = 18.01534d0
  integer :: i, i1, i2, idum, ireg, isrc, j
  integer :: ibc, ibrk, ir,np  

  logical :: continuation_flag
  character(len=1) :: backslash
  real*8 :: temp_real
  integer :: temp_int
  
  integer :: count, id
  
! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
!           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR

  type(region_type), pointer :: region
  type(condition_type), pointer :: condition
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_type), pointer :: material
  type(thermal_property_type), pointer :: thermal_property
  type(saturation_function_type), pointer :: saturation_function

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(solver_type), pointer :: solver
  type(stepper_type), pointer :: stepper
  
  realization => simulation%realization
  grid => realization%grid
  option => realization%option
  field => realization%field
  stepper => simulation%stepper
  solver => stepper%solver

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                              
  rewind(IUNIT1)  
    
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
!    call fiReadCard(word,card,ierr)
    card = trim(word)

    if (option%myrank == 0) print *, 'pflow_read:: ',card

    select case(trim(card))

!....................
      case ('MODE')

!....................
      case ('GRID')

!....................
      case ('PROC')
      
!....................
      case ('REGION','REGN')
        region => RegionCreate()
        call fiReadWord(string,region%name,.true.,ierr)
        call fiErrorMsg('regn','name',ierr) 
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('REGN',ierr)
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg('type','REGN', ierr)
        if (fiStringCompare(word,"BLOCK",5)) then ! block region
          if (grid%igrid /= STRUCTURED) then
            call printErrMsg(option,"BLOCK region not supported for &
                             &unstructured grid")
          endif
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('REGN',ierr)
          call fiReadInt(string,region%i1,ierr) 
          call fiErrorMsg('i1','REGN', ierr)
          call fiReadInt(string,region%i2,ierr)
          call fiErrorMsg('i2','REGN', ierr)
          call fiReadInt(string,region%j1,ierr)
          call fiErrorMsg('j1','REGN', ierr)
          call fiReadInt(string,region%j2,ierr)
          call fiErrorMsg('j2','REGN', ierr)
          call fiReadInt(string,region%k1,ierr)
          call fiErrorMsg('k1','REGN', ierr)
          call fiReadInt(string,region%k2,ierr)
          call fiErrorMsg('k2','REGN', ierr)
        else if (fiStringCompare(word,"LIST",4)) then
          word = ""
          call fiReadWord(string,word,.true.,ierr)
          if (fiStringCompare(word,"FILE",4)) then
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('REGN',ierr)
            call fiReadWord(string,word,.true.,ierr)
            call fiErrorMsg('filename','REGN', ierr)
            call RegionReadFromFile(region,word)
          else
            call RegionReadFromFile(region,IUNIT1)
          endif            
        else
          call printErrMsg(option,"REGION type not recognized")
        endif
        call RegionAddToList(region,realization%regions)      

!....................
      case ('CONDITION','COND')
        condition => ConditionCreate(option)
        call fiReadWord(string,condition%name,.true.,ierr)
        call fiErrorMsg('cond','name',ierr) 
        call ConditionRead(condition,option,IUNIT1)
        call ConditionAddToList(condition,realization%conditions)
      
!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1)
        call CouplerAddToList(coupler,realization%boundary_conditions)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1)
        call CouplerAddToList(coupler,realization%initial_conditions)
      
!....................
      case ('STRATIGRAPHY')
        strata => StrataCreate()
        call StrataRead(strata,IUNIT1)
        call StrataAddToList(strata,realization%strata)
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call CouplerRead(coupler,IUNIT1)
        call CouplerAddToList(coupler,realization%source_sinks)
      
!.....................
      case ('STRATA')
        strata => StrataCreate()
        call StrataRead(strata,IUNIT1)
        call StrataAddToList(strata,realization%strata)
      
!.....................
      case ('COMP') 
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo
!.....................
      case ('PHAS')
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PHASE',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo 
      
!....................

      case ('COUP')

        call fiReadStringErrorMsg('COUP',ierr)

        call fiReadInt(string,option%isync,ierr)
        call fiDefaultMsg('isync',ierr)

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') option%isync

!....................

      case ('GRAV','GRAVITY')

        call fiReadStringErrorMsg('GRAV',ierr)

        call fiReadDouble(string,temp_real,ierr)
        if (ierr /= 0) then
          call fiDefaultMsg('gravity',ierr)
        else
          call fiReadDouble(string,option%gravity(2),ierr)
          if (ierr /= 0) then
            option%gravity(:) = 0.d0
            option%gravity(3) = temp_real
          else
            option%gravity(1) = temp_real
            call fiReadDouble(string,option%gravity(3),ierr)
          endif
        endif

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,3pe12.4 &
            & )') option%gravity(1:3)

!....................

      case ('HDF5')
        realization%output_option%print_hdf5 = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_hdf5_velocities = .true.
            case('FLUX')
              realization%output_option%print_hdf5_flux_velocities = .true.
            case default
          end select
            
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *HDF5",10x,i1,/)') realization%output_option%print_hdf5

!.....................
      case ('INVERT_Z','INVERTZ')
        if (associated(grid%structured_grid)) then
          grid%structured_grid%invert_z_axis = .true.
          option%gravity(3) = -1.d0*option%gravity(3)
        endif
      
!....................

      case ('TECP')
        realization%output_option%print_tecplot = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_tecplot_velocities = .true.
            case('FLUX')
              realization%output_option%print_tecplot_flux_velocities = .true.
            case default
          end select
          
        enddo

        if (option%myrank == 0) &
          write(IUNIT2,'(/," *TECP",10x,i1,/)') realization%output_option%print_tecplot

!....................


      case ('OPTS')

        call fiReadStringErrorMsg('OPTS',ierr)

        call fiReadInt(string,option%write_init,ierr)
        call fiDefaultMsg('write_init',ierr)

        call fiReadInt(string,id,ierr)
        call fiDefaultMsg('iprint',ierr)

        call fiReadInt(string,option%imod,ierr)
        call fiDefaultMsg('mod',ierr)

        call fiReadInt(string,option%itecplot,ierr)
        call fiDefaultMsg('itecplot',ierr)

        call fiReadInt(string,option%iblkfmt,ierr)
        call fiDefaultMsg('iblkfmt',ierr)

        call fiReadInt(string,stepper%ndtcmx,ierr)
        call fiDefaultMsg('ndtcmx',ierr)

        call fiReadInt(string,option%iran_por,ierr)
        call fiDefaultMsg('iran_por',ierr)
  
        call fiReadDouble(string,option%ran_fac,ierr)
        call fiDefaultMsg('ran_fac',ierr)
    
        call fiReadInt(string,option%iread_perm,ierr)
        call fiDefaultMsg('iread_perm',ierr)
    
        call fiReadInt(string,option%iread_geom,ierr)
        call fiDefaultMsg('iread_geom',ierr)


        if (option%myrank == 0) &
          write(IUNIT2,'(/," *OPTS",/, &
            & "  write_init = ",3x,i2,/ &
            & "  imod       = ",3x,i2,/, &
            & "  itecplot   = ",3x,i2,/, &
            & "  iblkfmt    = ",3x,i2,/, &
            & "  ndtcmx     = ",3x,i2,/, &
            & "  iran_por   = ",3x,i2,/, &
            & "  ran_fac    = ",3x,1pe12.4,/, &
            & "  iread_perm = ",3x,i2,/, &
            & "  iread_geom = ",3x,i2 &
            & )') option%write_init,option%imod,option%itecplot, &
            option%iblkfmt,stepper%ndtcmx,option%iran_por,option%ran_fac, &
            option%iread_perm,option%iread_geom

!....................

      case ('TOLR')

        call fiReadStringErrorMsg('TOLR',ierr)

        call fiReadInt(string,stepper%stepmax,ierr)
        call fiDefaultMsg('stepmax',ierr)
  
        call fiReadInt(string,stepper%iaccel,ierr)
        call fiDefaultMsg('iaccel',ierr)

        call fiReadInt(string,stepper%newton_max,ierr)
        call fiDefaultMsg('newton_max',ierr)

        call fiReadInt(string,stepper%icut_max,ierr)
        call fiDefaultMsg('icut_max',ierr)

        call fiReadDouble(string,option%dpmxe,ierr)
        call fiDefaultMsg('dpmxe',ierr)

        call fiReadDouble(string,option%dtmpmxe,ierr)
        call fiDefaultMsg('dtmpmxe',ierr)
  
        call fiReadDouble(string,option%dcmxe,ierr)
        call fiDefaultMsg('dcmxe',ierr)

        call fiReadDouble(string,option%dsmxe,ierr)
        call fiDefaultMsg('dsmxe',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
          &"  flowsteps  = ",i6,/,      &
          &"  iaccel     = ",i3,/,      &
          &"  newtmx     = ",i3,/,      &
          &"  icutmx     = ",i3,/,      &
          &"  dpmxe      = ",1pe12.4,/, &
          &"  dtmpmxe    = ",1pe12.4,/, &
          &     "  dcmxe      = ",1pe12.4,/, &
          &"  dsmxe      = ",1pe12.4)') &
! For commented-out lines to work with the Sun f95 compiler, we have to 
! terminate the string in the line above; otherwise, the compiler tries to
! include the commented-out line as part of the continued string.
          stepper%stepmax,stepper%iaccel,stepper%newton_max,stepper%icut_max, &
          option%dpmxe,option%dtmpmxe,option%dcmxe, option%dsmxe

!....................

      case ('DXYZ')
      
        if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
          call StructuredGridReadDXYZ(grid%structured_grid,option)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "DXYZ" not supported for unstructured grid'
            stop
        endif

!....................

      case('ORIG','ORIGIN')
        call fiReadDouble(string,grid%origin(X_DIRECTION),ierr)
        call fiErrorMsg('X direction','Origin',ierr)
        call fiReadDouble(string,grid%origin(Y_DIRECTION),ierr)
        call fiErrorMsg('Y direction','Origin',ierr)
        call fiReadDouble(string,grid%origin(Z_DIRECTION),ierr)
        call fiErrorMsg('Z direction','Origin',ierr)
        
!....................

      case('RAD0')
    
        if (grid%igrid == STRUCTURED) then  ! look for processor decomposition
          call fiReadDouble(string,grid%structured_grid%Radius_0,ierr)
          call fiDefaultMsg('R_0',ierr)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "RAD0" not supported for unstructured grid'
            stop
        endif


      case ('DIFF')

        call fiReadStringErrorMsg('DIFF',ierr)

        call fiReadDouble(string,option%difaq,ierr)
        call fiDefaultMsg('difaq',ierr)

        call fiReadDouble(string,option%delhaq,ierr)
        call fiDefaultMsg('delhaq',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          option%difaq,option%delhaq

!....................

      case ('RCTR')

        call fiReadStringErrorMsg('RCTR',ierr)

        call fiReadInt(string,option%ityprxn,ierr)
        call fiDefaultMsg('ityprxn',ierr)

        call fiReadDouble(string,option%rk,ierr)
        call fiDefaultMsg('rk',ierr)

        call fiReadDouble(string,option%phis0,ierr)
        call fiDefaultMsg('phis0',ierr)

        call fiReadDouble(string,option%areas0,ierr)
        call fiDefaultMsg('areas0',ierr)

        call fiReadDouble(string,option%pwrsrf,ierr)
        call fiDefaultMsg('pwrsrf',ierr)

        call fiReadDouble(string,option%vbars,ierr)
        call fiDefaultMsg('vbars',ierr)

        call fiReadDouble(string,option%ceq,ierr)
        call fiDefaultMsg('ceq',ierr)

        call fiReadDouble(string,option%delHs,ierr)
        call fiDefaultMsg('delHs',ierr)

        call fiReadDouble(string,option%delEs,ierr)
        call fiDefaultMsg('delEs',ierr)

        call fiReadDouble(string,option%wfmts,ierr)
        call fiDefaultMsg('wfmts',ierr)

        if (option%myrank == 0) &
        write(IUNIT2,'(/," *RCTR",/, &
          & "  ityp   = ",3x,i3,/, &
          & "  rk     = ",3x,1pe12.4," [mol/cm^2/s]",/, &
          & "  phis0  = ",3x,1pe12.4," [-]",/, &
          & "  areas0 = ",3x,1pe12.4," [1/cm]",/, &
          & "  pwrsrf = ",3x,1pe12.4," [-]",/, &
          & "  vbars  = ",3x,1pe12.4," [cm^3/mol]",/, &
          & "  ceq    = ",3x,1pe12.4," [mol/L]",/, &
          & "  delHs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  delEs  = ",3x,1pe12.4," [J/kg]",/, &
          & "  wfmts  = ",3x,1pe12.4," [g/mol]" &
          & )') option%ityprxn,option%rk,option%phis0,option%areas0,option%pwrsrf, &
          option%vbars,option%ceq,option%delHs,option%delEs,option%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
        option%rk = option%rk * option%areas0 * 1.d3
        option%vbars = option%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
        option%delHs = option%delHs * option%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!        option%delHs = option%delHs * option%scale ! convert J/kmol -> MJ/kmol

!....................

      case ('RADN')

        call fiReadStringErrorMsg('RADN',ierr)

        call fiReadDouble(string,option%ret,ierr)
        call fiDefaultMsg('ret',ierr)

        call fiReadDouble(string,option%fc,ierr)
        call fiDefaultMsg('fc',ierr)

        if (option%myrank==0) write(IUNIT2,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          option%ret,option%fc

!....................


      case ('PHAR')

        call fiReadStringErrorMsg('PHAR',ierr)

        call fiReadDouble(string,option%qu_kin,ierr)
        call fiDefaultMsg('TransReaction',ierr)
        if (option%myrank==0) write(IUNIT2,'(/," *PHAR ",1pe12.4)')option%qu_kin
        option%yh2o_in_co2 = 0.d0
        if (option%qu_kin > 0.d0) option%yh2o_in_co2 = 1.d-2 ! check this number!
     
!......................

      case('RICH')
        call fiReadStringErrorMsg('RICH',ierr)
        call fiReadDouble(string,option%pref,ierr)
        call fiDefaultMsg('Ref. Pressure',ierr) 

!......................

      case('BRIN','BRINE')
        call fiReadStringErrorMsg('BRIN',ierr)
        call fiReadDouble(string,option%m_nacl,ierr)
        call fiDefaultMsg('NaCl Concentration',ierr) 

        call fiReadWord(string,word,.false.,ierr)
        call fiWordToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            option%m_nacl = option%m_nacl /fmwnacl/(1.D0-option%m_nacl)
          case('MOLE')    
            option%m_nacl = option%m_nacl /fmwh2o/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (option%myrank == 0) print *, option%m_nacl
!......................

      case ('NO_PRINT_CONVERGENCE')
        option%print_convergence = PETSC_FALSE

      case ('NO_INF_NORM','NO_INFINITY_NORM')
        option%check_infinity_norm = PETSC_FALSE

      case ('NO_FORCE_ITERATION')
        option%force_at_least_1_iteration = PETSC_FALSE

      case ('PRINT_DETAILED_CONVERGENCE')
        option%print_detailed_convergence = PETSC_TRUE

!....................

      case ('SOLV')
    
        call fiReadStringErrorMsg('SOLV',ierr)

!       call fiReadDouble(string,eps,ierr)
!       call fiDefaultMsg('eps',ierr)

        call fiReadDouble(string,solver%atol,ierr)
        call fiDefaultMsg('atol_petsc',ierr)

        call fiReadDouble(string,solver%rtol,ierr)
        call fiDefaultMsg('rtol_petsc',ierr)

        call fiReadDouble(string,solver%stol,ierr)
        call fiDefaultMsg('stol_petsc',ierr)
      
        solver%dtol=1.D5
!       if (option%use_ksp == 1) then
        call fiReadDouble(string,solver%dtol,ierr)
        call fiDefaultMsg('dtol_petsc',ierr)
!       endif
   
        call fiReadInt(string,solver%maxit,ierr)
        call fiDefaultMsg('maxit',ierr)
      
        call fiReadInt(string,solver%maxf,ierr)
        call fiDefaultMsg('maxf',ierr)
       
        call fiReadInt(string,option%idt_switch,ierr)
        call fiDefaultMsg('idt',ierr)
        
        solver%inf_tol = solver%atol
        call fiReadDouble(string,solver%inf_tol,ierr)
        call fiDefaultMsg('inf_tol_pflow',ierr)
 
        if (option%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
          &"  atol_petsc   = ",1pe12.4,/, &
          &"  rtol_petsc   = ",1pe12.4,/, &
          &"  stol_petsc   = ",1pe12.4,/, &
          &"  dtol_petsc   = ",1pe12.4,/, &
          &"  maxit        = ",8x,i5,/, &
          &"  maxf        = ",8x,i5,/, &
          &"  idt         = ",8x,i5 &
          &    )') &
           solver%atol,solver%rtol,solver%stol,solver%dtol,solver%maxit, &
           solver%maxf,option%idt_switch

! The line below is a commented-out portion of the format string above.
! We have to put it here because of the stupid Sun compiler.
!    &"  eps          = ",1pe12.4,/, &

!....................

      case ('THRM','THERMAL_PROPERTY','THERMAL_PROPERTIES')

        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('THRM',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',3)) exit
       
          count = count + 1
          thermal_property => ThermalPropertyCreate()
      
          call fiReadInt(string,thermal_property%id,ierr)
          call fiErrorMsg('id','THRM', ierr)

          call fiReadDouble(string,thermal_property%rock_density,ierr)
          call fiErrorMsg('rock density','THRM', ierr)

          call fiReadDouble(string,thermal_property%spec_heat,ierr)
          call fiErrorMsg('cpr','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_dry,ierr)
          call fiErrorMsg('ckdry','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_wet,ierr)
          call fiErrorMsg('ckwet','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%tort_bin_diff,ierr)
          call fiErrorMsg('tau','THRM', ierr)

          call fiReadDouble(string,thermal_property%vap_air_diff_coef,ierr)
          call fiErrorMsg('cdiff','THRM', ierr)

          call fiReadDouble(string,thermal_property%exp_binary_diff,ierr)
          call fiErrorMsg('cexp','THRM', ierr)

        !scale thermal properties
          thermal_property%spec_heat = option%scale * &
                                       thermal_property%spec_heat
          thermal_property%therm_cond_dry = option%scale * &
                                            thermal_property%therm_cond_dry
          thermal_property%therm_cond_wet = option%scale * &
                                            thermal_property%therm_cond_wet
          
          call ThermalAddPropertyToList(thermal_property, &
                                        realization%thermal_properties)
        enddo
        
        ! allocate dynamic arrays holding saturation function information
        allocate(option%rock_density(count))
        allocate(option%cpr(count))
        allocate(option%dencpr(count))
        allocate(option%ckdry(count))
        allocate(option%ckwet(count))
        allocate(option%tau(count))
        allocate(option%cdiff(count))
        allocate(option%cexp(count))
        
        ! fill arrays with values from linked list
        thermal_property => realization%thermal_properties
        do
        
          if (.not.associated(thermal_property)) exit
          
          id = thermal_property%id
          
          if (id > count) then
            call printErrMsg(option,'Thermal property id greater than &
                                    &number of thermal properties')
          endif
                    
          option%rock_density(id) = thermal_property%rock_density
          option%cpr(id) = thermal_property%spec_heat
          option%dencpr(id) = thermal_property%rock_density * &
                              thermal_property%spec_heat
          option%ckdry(id) = thermal_property%therm_cond_dry
          option%ckwet(id) = thermal_property%therm_cond_wet
          option%tau(id) = thermal_property%tort_bin_diff
          option%cdiff(id) = thermal_property%vap_air_diff_coef
          option%cexp(id) = thermal_property%exp_binary_diff
          
          thermal_property => thermal_property%next
          
        enddo
        
        do i=1,count
          if (option%rock_density(i) < 1.d-40) then
            call printErrMsg(option,'Thermal property ids must be numbered &
                             &consecutively from 1 to N')
          endif
        enddo
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *THRM: ",i3)') count
          write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, count
            write(IUNIT2,'(i4,1p7e11.4)') i,option%rock_density(i), &
            option%cpr(i),option%ckdry(i),option%ckwet(i), &
            option%tau(i),option%cdiff(i),option%cexp(i)
          enddo
        endif

!....................

      case ('PCKR','SATURATION_FUNCTION','SATURATION_FUNCTIONS')
      
        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PCKR',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',3)) exit
       
          count = count + 1
          saturation_function => SaturationFunctionCreate(option)
          
          call fiReadInt(string,saturation_function%id,ierr)
          call fiErrorMsg('id','PCKR', ierr)
          
          call fiReadInt(string,saturation_function%itype,ierr)
          call fiErrorMsg('icaptype','PCKR', ierr)
      
          select case(option%imode)
            case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
              do np=1, option%nphase
                call fiReadDouble(string,saturation_function%Sr(np),ierr)
                call fiErrorMsg('Sr','PCKR', ierr)
              enddo 
            case default
              call fiReadDouble(string,saturation_function%Sr(1),ierr)
              call fiErrorMsg('Sr','PCKR', ierr)
          end select
        
          call fiReadDouble(string,saturation_function%m,ierr)
          call fiErrorMsg('pckrm','PCKR', ierr)
          saturation_function%lambda = saturation_function%m / &
                                      (1.d0-saturation_function%m)

          call fiReadDouble(string,saturation_function%alpha,ierr)
          call fiErrorMsg('alpha','PCKR', ierr)

          call fiReadDouble(string,saturation_function%pcwmax,ierr)
          call fiErrorMsg('pcwmax','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%betac,ierr)
          call fiErrorMsg('pbetac','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%power,ierr)
          call fiErrorMsg('pwrprm','PCKR', ierr)
          
          call SaturationFunctionAddToList(saturation_function, &
                                           realization%saturation_functions)

        enddo
        
        ! allocate dynamic arrays holding saturation function information
        allocate(option%icaptype(count))
        option%icaptype = 0
  
        select case(option%imode)
          case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
            allocate(option%sir(1:option%nphase,count))
          case default
            allocate(option%swir(count))
        end select
  
        allocate(option%lambda(count))
        allocate(option%alpha(count))
        allocate(option%pckrm(count))
        allocate(option%pcwmax(count))
        allocate(option%pcbetac(count))
        allocate(option%pwrprm(count))

        ! fill arrays with values from linked list
        saturation_function => realization%saturation_functions
        do 
        
          if (.not.associated(saturation_function)) exit
          
          id = saturation_function%id
          
          if (id > count) then
            call printErrMsg(option,'Saturation function id greater than &
                                    &number of saturation functions')
          endif
          
          option%icaptype(id) = saturation_function%itype
          select case(option%imode)
            case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
              do i=1,option%nphase
                option%sir(i,id) = saturation_function%Sr(i)
              enddo
            case default
              option%swir(id) = saturation_function%Sr(1)
          end select
          option%lambda(id) = saturation_function%lambda
          option%alpha(id) = saturation_function%alpha
          option%pckrm(id) = saturation_function%m
          option%pcwmax(id) = saturation_function%pcwmax
          option%pcbetac(id) = saturation_function%betac
          option%pwrprm(id) = saturation_function%power
          
          saturation_function => saturation_function%next
          
        enddo
        
        ! check to ensure that all saturation functions were set based on id
        do id = 1,count
          if (option%icaptype(id) == 0) then
            call printErrMsg(option,'Saturation function ids must be numbered &
                               &consecutively from 1 to N')
          endif
        enddo

        if (option%imode == MPH_MODE .or. &
            option%imode == VADOSE_MODE .or. &
            option%imode == FLASH_MODE .or. &
            option%imode == RICHARDS_MODE) then
          call pckr_init(option%nphase,count,grid%nlmax, &
                         option%icaptype,option%sir, option%pckrm, &
                         option%lambda,option%alpha,option%pcwmax, &
                         option%pcbetac,option%pwrprm)
        endif 

      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *PCKR: ",i3)') ireg
          write(IUNIT2,'("  icp swir    lambda         alpha")')
          do j = 1, count
            i=option%icaptype(j)
            if (option%imode == MPH_MODE .or. &
                option%imode == OWG_MODE .or. &
                option%imode == VADOSE_MODE .or. &
                option%imode == FLASH_MODE .or. &
                option%imode == RICHARDS_MODE) then
              write(IUNIT2,'(i4,1p8e12.4)') i,(option%sir(np,i),np=1, &
                option%nphase),option%lambda(i),option%alpha(i), &
                option%pcwmax(i),option%pcbetac(i),option%pwrprm(i)
            else
              write(IUNIT2,'(i4,1p7e12.4)') i,option%swir(i), &
                option%lambda(i),option%alpha(i),option%pcwmax(i), &
                option%pcbetac(i),option%pwrprm(i)
            endif
          enddo
        end if

        if (option%imode == MPH_MODE .or. &
            option%imode == VADOSE_MODE .or. &
            option%imode == FLASH_MODE .or. &
            option%imode == RICHARDS_MODE) then
          deallocate(option%icaptype, option%pckrm, option%lambda, &
                     option%alpha,option%pcwmax, option%pcbetac, &
                     option%pwrprm)
        endif 
 

!....................
      
      case ('PHIK','MATERIAL','MATERIALS')

        count = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PHIK',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',3)) exit
       
          count = count + 1
          material => MaterialCreate()

          call fiReadWord(string,material%name,.true.,ierr)
          call fiErrorMsg('name','PHIK', ierr)
                
          call fiReadInt(string,material%id,ierr)
          call fiErrorMsg('id','PHIK', ierr)
                
          call fiReadInt(string,material%icap,ierr)
          call fiErrorMsg('icap','PHIK', ierr)
  
          call fiReadInt(string,material%ithrm,ierr)
          call fiErrorMsg('ithrm','PHIK', ierr)
  
          call fiReadDouble(string,material%porosity,ierr)
          call fiErrorMsg('por','PHIK', ierr)
          
          call fiReadDouble(string,material%tortuosity,ierr)
          call fiErrorMsg('tor','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(1,1),ierr)
          call fiErrorMsg('permx','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(2,2),ierr)
          call fiErrorMsg('permy','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(3,3),ierr)
          call fiErrorMsg('permz','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability_pwr,ierr)
          call fiErrorMsg('permpwr','PHIK', ierr)
          
          material%permeability(1:3,1:3) = material%permeability(1:3,1:3)
          
          call MaterialAddToList(material,realization%materials)
          
        enddo          

!....................
      
      case ('INIT')
    
#if 0
! INIT is deprecated by condition/region coupling
        call fiReadInt(string,option%iread_init,ierr) 
        call fiDefaultMsg('iread_init',ierr)
      
        if (option%myrank==0) then
          write(IUNIT2,'(/," *INIT: iread = ",i2)') option%iread_init
        endif
      
        if (option%iread_init == 0 .or. option%iread_init == 2) then
      
          ireg = 0
          do
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('INIT',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ireg = ireg + 1
            
!GEH - Structured Grid Dependence - Begin
            call fiReadInt(string,option%i1ini(ireg),ierr) 
            call fiDefaultMsg('i1',ierr)
            call fiReadInt(string,option%i2ini(ireg),ierr)
            call fiDefaultMsg('i2',ierr)
            call fiReadInt(string,option%j1ini(ireg),ierr)
            call fiDefaultMsg('j1',ierr)
            call fiReadInt(string,option%j2ini(ireg),ierr)
            call fiDefaultMsg('j2',ierr)
            call fiReadInt(string,option%k1ini(ireg),ierr)
            call fiDefaultMsg('k1',ierr)
            call fiReadInt(string,option%k2ini(ireg),ierr)
            call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

            select case(option%imode)
              case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)                
                call fiReadInt(string,option%iphas_ini(ireg),ierr)
                call fiDefaultMsg('iphase',ierr)
       
                do j=1,option%ndof
                  call fiReadDouble(string,field%xx_ini(j,ireg),ierr)
                  call fiDefaultMsg('xxini',ierr)
                enddo
              case default
                call fiReadDouble(string,option%pres_ini(ireg),ierr)
                call fiDefaultMsg('pres',ierr)
  
                call fiReadDouble(string,field%temp_ini(ireg),ierr)
                call fiDefaultMsg('temp',ierr)
  
                call fiReadDouble(string,field%sat_ini(ireg),ierr)
                call fiDefaultMsg('sat',ierr)
!                field%sat_ini(ireg)=1.D0 - field%sat_ini(ireg)
  
                call fiReadDouble(string,field%conc_ini(ireg),ierr)
                call fiDefaultMsg('conc',ierr)
            end select
          enddo
      
          option%iregini = ireg
      
          if (option%myrank==0) then
            write(IUNIT2,'("  ireg = ",i4)') option%iregini
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]   &
              &   ",    "sl [-]      c [mol/L]")')
            do ireg = 1, option%iregini
!GEH - Structured Grid Dependence - Begin
              select case(option%imode)
                case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)                
                  write(IUNIT2,'(7i4,1p10e12.4)') &
                    option%i1ini(ireg),option%i2ini(ireg), &
                    option%j1ini(ireg),option%j2ini(ireg), &
                    option%k1ini(ireg),option%k2ini(ireg), &
                    option%iphas_ini(ireg),(field%xx_ini(np,ireg),np =1, &
                                            option%ndof)
                case default
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    option%i1ini(ireg),option%i2ini(ireg), &
                    option%j1ini(ireg),option%j2ini(ireg), &
                    option%k1ini(ireg),option%k2ini(ireg), &
                    option%pres_ini(ireg),field%temp_ini(ireg), &
                    field%sat_ini(ireg),field%conc_ini(ireg)
              end select
!GEH - Structured Grid Dependence - End
            enddo
          endif

        else if (option%iread_init == 1) then
    
!     read in initial conditions from file: pflow_init.dat
          if (option%myrank == 0) then
            write(*,*) '--> read in initial conditions from file: &
                        &pflow_init.dat'
  
            open(IUNIT3, file='pflow_init.dat', action="read", status="old")

            ireg = 0
            do
              call fiReadFlotranString(IUNIT3,string,ierr)
!             call fiReadStringErrorMsg('INIT',ierr)

              if (string(1:1) == '.' .or. string(1:1) == '/') exit
              ireg = ireg + 1

!GEH - Structured Grid Dependence - Begin
              call fiReadInt(string,option%i1ini(ireg),ierr) 
              call fiDefaultMsg('i1',ierr)
              call fiReadInt(string,option%i2ini(ireg),ierr)
              call fiDefaultMsg('i2',ierr)
              call fiReadInt(string,option%j1ini(ireg),ierr)
              call fiDefaultMsg('j1',ierr)
              call fiReadInt(string,option%j2ini(ireg),ierr)
              call fiDefaultMsg('j2',ierr)
              call fiReadInt(string,option%k1ini(ireg),ierr)
              call fiDefaultMsg('k1',ierr)
              call fiReadInt(string,option%k2ini(ireg),ierr)
              call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

              select case(option%imode)
                case(MPH_MODE,OWG_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
                  call fiReadInt(string,option%iphas_ini(ireg),ierr)
                  call fiDefaultMsg('iphase_ini',ierr)
            
                  do j=1,option%ndof
                    call fiReadDouble(string,field%xx_ini(j,ireg),ierr)
                    call fiDefaultMsg('xx_ini',ierr)
                  enddo
                case default
  
                  call fiReadDouble(string,option%pres_ini(ireg),ierr)
                  call fiDefaultMsg('pres',ierr)
  
                  call fiReadDouble(string,field%temp_ini(ireg),ierr)
                  call fiDefaultMsg('temp',ierr)
  
                  call fiReadDouble(string,field%sat_ini(ireg),ierr)
                  call fiDefaultMsg('sat',ierr)

                  call fiReadDouble(string,field%conc_ini(ireg),ierr)
                  call fiDefaultMsg('conc',ierr)
              end select
       
            enddo
            option%iregini = ireg
            close(IUNIT3)
          endif
        endif
#endif
!....................

      case ('TIME')

        call fiReadStringErrorMsg('TIME',ierr)
      
        call fiReadWord(string,word,.false.,ierr)
      
        realization%output_option%tunit = trim(word)

        if (realization%output_option%tunit == 's') then
          realization%output_option%tconv = 1.d0
        else if (realization%output_option%tunit == 'm') then
          realization%output_option%tconv = 60.d0
        else if (realization%output_option%tunit == 'h') then
          realization%output_option%tconv = 60.d0 * 60.d0
        else if (realization%output_option%tunit == 'd') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0
        else if (realization%output_option%tunit == 'mo') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
        else if (realization%output_option%tunit == 'y') then
          realization%output_option%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
        else
          if (option%myrank == 0) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') realization%output_option%tunit
          endif
          stop
        endif

        continuation_flag = .true.
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = .false.
          if (index(string,backslash) > 0) continuation_flag = .true.
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              waypoint%print_output = .true.              
              call WaypointInsertInList(waypoint,stepper%waypoints)
            endif
          enddo
        enddo
        
        ! make last waypoint final
        waypoint%final = .true.

!....................

      case ('DTST')

        call fiReadStringErrorMsg('DTST',ierr)

        call fiReadDouble(string,stepper%dt_min,ierr)
        call fiDefaultMsg('dt_min',ierr)
            
        continuation_flag = .true.
        temp_int = 0       
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = .false.
          if (index(string,backslash) > 0) continuation_flag = .true.
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              call fiReadDouble(string,waypoint%dt_max,ierr)
              call fiErrorMsg('dt_max','dtst',ierr)
              if (temp_int == 0) stepper%dt_max = waypoint%dt_max
              call WaypointInsertInList(waypoint,stepper%waypoints)
              temp_int = temp_int + 1
            endif
          enddo
        enddo
        
        option%dt = stepper%dt_min
      
        option%dt = realization%output_option%tconv * option%dt
        stepper%dt_min = realization%output_option%tconv * stepper%dt_min
        stepper%dt_max = realization%output_option%tconv * stepper%dt_max

!....................
   
      case ('BRK')
        print *, 'BRK (breakthrough) needs to be implemented'
        stop
#if 0
! Needs implementation
        ibrk = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BRK',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ibrk = ibrk + 1

!GEH - Structured Grid Dependence - Begin
          call fiReadInt(string,option%i1brk(ibrk),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,option%i2brk(ibrk),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,option%j1brk(ibrk),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,option%j2brk(ibrk),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,option%k1brk(ibrk),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,option%k2brk(ibrk),ierr)
          call fiDefaultMsg('k2',ierr)
!GEH - Structured Grid Dependence - End

          call fiReadInt(string,option%ibrktyp(ibrk),ierr)
          call fiDefaultMsg('ibrktyp',ierr)

          call fiReadInt(string,option%ibrkface(ibrk),ierr)
          call fiDefaultMsg('ibrkface',ierr)
        enddo
        option%ibrkcrv = ibrk
            
        if (option%myrank==0) then
          write(IUNIT2,'(/," *BRK: ibrk = ",i4)') option%ibrkcrv
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2  ibrktyp  ibrkface  ")')
          do ibrk = 1, option%ibrkcrv
            write(IUNIT2,'(6i4,4x,i2,7x,i2)') &
              option%i1brk(ibrk),option%i2brk(ibrk), &
              option%j1brk(ibrk),option%j2brk(ibrk), &
              option%k1brk(ibrk),option%k2brk(ibrk),option%ibrktyp(ibrk), &
              option%ibrkface(ibrk)
          enddo
        endif

        if (option%ndof == 1) option%ibrkcrv = 0
#endif
!....................
      case('SDST')
        print *, 'SDST needs to be implemented'
        stop
#if 0
! Needs implementation         
        allocate(stepper%steady_eps(option%ndof))
        do j=1,option%ndof
          call fiReadDouble(string,stepper%steady_eps(j),ierr)
          call fiDefaultMsg('steady tol',ierr)
        enddo
        if (option%myrank==0) write(IUNIT2,'(/," *SDST ",/, &
          &"  dpdt        = ",1pe12.4,/, &
          &"  dtmpdt        = ",1pe12.4,/, &
          &"  dcdt        = ",1pe12.4)') &
          stepper%steady_eps
#endif
!....................
      case default
    
        if (option%myrank == 0) then
          print *, "Error reading input file: keyword (", trim(word), &
                   ") not found. Terminating."
        endif
        call PetscFinalize(ierr)
        stop

    end select

  enddo

  close(IUNIT1)
  
end subroutine readInput

! ************************************************************************** !
!
! initAccumulation: Initializes accumulation term?
! author: 
! date:
!
! ************************************************************************** !
subroutine initAccumulation(realization)
!  use water_eos_module
!  use TTPHASE_module
!  use Flash_module
!  use MPHASE_module
!  use OWG_module
!  use Vadose_module
 use Richards_module , only: pflow_richards_initaccum  ! for some reason intel compiler fails without "only" clause

  use Realization_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  real*8, pointer :: den_p(:), pressure_p(:), temp_p(:), h_p(:)
  real*8 :: dw_kg,dl,hl
  integer :: m, ierr
  
  option => realization%option
  
#if 0  
  ! needs to be implemented
  if ( grid%use_owg == PETSC_TRUE) then
    call pflow_owg_initaccum(grid)
  else if (grid%use_mph == PETSC_TRUE) then
    call pflow_mphase_initaccum(grid)
  else if (grid%use_richards == PETSC_TRUE) then
#endif    
  if (option%imode == RICHARDS_MODE) then
    call pflow_richards_initaccum(realization)
#if 0    
  ! needs to be implemented
  else if (grid%use_flash == PETSC_TRUE) then
    call pflow_flash_initaccum(grid)
  else if (grid%use_vadose == PETSC_TRUE) then
    call pflow_vadose_initaccum(grid)
  else if (grid%use_2ph == PETSC_TRUE) then
    call pflow_2phase_initaccum(grid)
  else if (grid%ndof > 1) then
 !   call VecSet(grid%iphas,1.d0,ierr)
    call VecGetArrayF90(grid%pressure, pressure_p, ierr)
    call VecGetArrayF90(grid%temp, temp_p, ierr)
    call VecGetArrayF90(grid%density, den_p, ierr)
    call VecGetArrayF90(grid%h, h_p, ierr)
    do m = 1, grid%nlmax
      call wateos_noderiv(temp_p(m),pressure_p(m),dw_kg,dl,hl,grid%scale,ierr)
      den_p(m) = dl
      h_p(m) = hl
    enddo
    call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
    call VecRestoreArrayF90(grid%temp, temp_p, ierr)
    call VecRestoreArrayF90(grid%density, den_p, ierr)
    call VecRestoreArrayF90(grid%h, h_p, ierr)

  else

    call VecGetArrayF90(grid%pressure, pressure_p, ierr)
    call VecGetArrayF90(grid%temp, temp_p, ierr)
    call VecGetArrayF90(grid%density, den_p, ierr)
    do m = 1, grid%nlmax
      call wateos_noderiv(temp_p(m),pressure_p(m),dw_kg,dl,hl,grid%scale,ierr)
      den_p(m) = dl
    enddo
    call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
    call VecRestoreArrayF90(grid%temp, temp_p, ierr)
    call VecRestoreArrayF90(grid%density, den_p, ierr)
#endif
  endif

end subroutine initAccumulation

! ************************************************************************** !
!
! setMode: Sets the flow mode (richards, vadose, mph, etc.)
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine setMode(option,mcomp,mphas)

  use Option_module
  use Fileio_module

  implicit none
  
#include "definitions.h"  

  type(option_type) :: option
  integer :: mcomp, mphas
  
  call fiCharsToLower(option%mode,len_trim(option%mode))
  if (fiStringCompare(option%mode,"richards",8)) then
    option%imode = RICHARDS_MODE
#if 0  
  ! needs to be implemented
  else if (fiStringCompare(option%mode,"MPH",3)) then
  else if (fiStringCompare(option%mode,"",#)) then
#endif  
  endif 
  
  if (option%imode /= LIQUID_MODE .and. &
      option%imode /= COND_MODE .and. &
      option%imode /= TH_MODE .and. &
      option%imode /= THC_MODE .and. &
      option%imode /= TWOPH_MODE .and. &
      option%imode /= MPH_MODE .and. &
      option%imode /= FLASH_MODE .and. &
      option%imode /= OWG_MODE .and. &
      option%imode /= VADOSE_MODE .and. &
      option%imode /= RICHARDS_MODE) then 

    if (mcomp >0 .and. mphas>0)then
      if (option%imode /= COND_MODE .and. mcomp ==1)then
        option%imode = COND_MODE
        option%mode = 'cond'
        option%nphase = 1; option%ndof =1
      endif
      if (option%imode /= PETSC_FALSE .and. mcomp == 32)then
        option%imode = COND_MODE
        option%mode = 'cond'
        option%nphase = 1; option%ndof =1
      endif
      if (option%imode /= TH_MODE .and. mcomp == 33 .and. mphas == 3)then
        option%imode = TH_MODE
        option%mode = 'th'
        option%nphase = 1; option%ndof =2
      endif
      if (option%imode /= THC_MODE .and. mcomp == 37)then
        option%imode = THC_MODE
        option%mode = 'thc'
        option%nphase = 1; option%ndof =3
      endif
      if (option%imode /= MPH_MODE .and. mcomp == 35)then
        option%imode = MPH_MODE
        option%mode = 'mph'
        option%nphase = 2; option%ndof =3; option%nspec =2 
      endif
      if (option%imode /= VADOSE_MODE .and. mcomp == 49)then
        option%imode = VADOSE_MODE
        option%mode = 'vadose'
        option%nphase = 2; option%ndof =3; option%nspec =2 
      endif
      if (option%imode /= RICHARDS_MODE .and. mcomp == 33 .and. mphas == 11) then
        option%imode = RICHARDS_MODE
        option%mode = 'richards'
        option%nphase = 1; option%ndof = 2
        if (option%nspec > 1) then
          option%ndof = option%nspec +1
        endif
      endif
    endif
    if (option%imode /= OWG_MODE .and. mcomp == 11)then
      option%imode = OWG_MODE
      option%mode = 'owg'
      option%nphase = 3; option%ndof =3; option%nspec =3 
    endif
  endif
  if (option%imode == NULL_MODE) then
    call printErrMsg(option,"No mode specified")
  endif       

end subroutine setMode

! ************************************************************************** !
!
! assignMaterialPropToRegions: Assigns material properties to 
!                                    associated regions in the model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine assignMaterialPropToRegions(realization)

  use Realization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none
  
  type(realization_type) :: realization
  
  PetscScalar, pointer :: icap_loc_p(:)
  PetscScalar, pointer :: ithrm_loc_p(:)
  PetscScalar, pointer :: por0_p(:)
  PetscScalar, pointer :: perm_xx_p(:)
  PetscScalar, pointer :: perm_yy_p(:)
  PetscScalar, pointer :: perm_zz_p(:)
  PetscScalar, pointer :: perm_pow_p(:)
  PetscScalar, pointer :: tor_loc_p(:)
  
  integer :: icell, local_id, ghosted_id
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  
  type(material_type), pointer :: material
  type(region_type), pointer :: region
  
  option => realization%option
  grid => realization%grid
  field => realization%field
  
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecGetArrayF90(field%porosity0,por0_p,ierr)
  call VecGetArrayF90(field%tor_loc,tor_loc_p,ierr)
  call VecGetArrayF90(field%perm0_xx,perm_xx_p,ierr)
  call VecGetArrayF90(field%perm0_yy,perm_yy_p,ierr)
  call VecGetArrayF90(field%perm0_zz,perm_zz_p,ierr)
  call VecGetArrayF90(field%perm_pow,perm_pow_p,ierr)

  strata => realization%strata%first
  do
    if (.not.associated(strata)) exit
    
    region => strata%region
    material => strata%material
    do icell=1,region%num_cells
      local_id = region%cell_ids(icell)
      ghosted_id = grid%nL2G(local_id)
      icap_loc_p(ghosted_id) = material%icap
      ithrm_loc_p(ghosted_id) = material%ithrm
      por0_p(local_id) = material%porosity
      tor_loc_p(ghosted_id) = material%tortuosity
      perm_xx_p(local_id) = material%permeability(1,1)
      perm_yy_p(local_id) = material%permeability(2,2)
      perm_zz_p(local_id) = material%permeability(3,3)
      perm_pow_p(local_id) = material%permeability_pwr
    enddo
    strata => strata%next
  enddo

  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr)
  call VecRestoreArrayF90(field%porosity0,por0_p,ierr)
  call VecRestoreArrayF90(field%perm0_xx,perm_xx_p,ierr)
  call VecRestoreArrayF90(field%perm0_yy,perm_yy_p,ierr)
  call VecRestoreArrayF90(field%perm0_zz,perm_zz_p,ierr)
  call VecRestoreArrayF90(field%perm_pow,perm_pow_p,ierr)
  call VecRestoreArrayF90(field%tor_loc,tor_loc_p,ierr)

  call GridGlobalToLocal(realization%grid,field%porosity0, &
                         field%porosity_loc,ONEDOF)
  call GridGlobalToLocal(realization%grid,field%perm0_xx, &
                         field%perm_xx_loc,ONEDOF)  
  call GridGlobalToLocal(realization%grid,field%perm0_yy, &
                         field%perm_yy_loc,ONEDOF)  
  call GridGlobalToLocal(realization%grid,field%perm0_zz, &
                         field%perm_zz_loc,ONEDOF)   
  
end subroutine assignMaterialPropToRegions

! ************************************************************************** !
!
! assignInitialConditions: Assigns initial conditions to model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine assignInitialConditions(realization)

  use Realization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  
  use MPHASE_module, only : pflow_mphase_setupini
  use Richards_module, only : pflow_Richards_setupini
  use hydrostat_module

  implicit none
  
  type(realization_type) :: realization
  
  PetscScalar, pointer :: pressure_p(:)
  PetscScalar, pointer :: temp_p(:)
  PetscScalar, pointer :: sat_p(:)
  PetscScalar, pointer :: conc_p(:)
  PetscScalar, pointer :: xmol_p(:)
  
  integer :: icell, iconn, count, jn1, jn2
  integer :: local_id, ghosted_id, iend, ibegin
  PetscScalar, pointer :: xx_p(:), iphase_loc_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(coupler_type), pointer :: initial_condition
    
  option => realization%option
  field => realization%field
    
  select case(option%imode)
  ! needs to be implemented
#if 0
    case(FLASH_MODE)
      call pflow_flash_setupini(realization)
#endif  
    case(RICHARDS_MODE)
      call pflow_richards_setupini(realization)
    case(MPH_MODE)
      call pflow_mphase_setupini(realization)
#if 0
  ! needs to be implemented
    case(OWG_MODE)
      call pflow_owg_setupini(realization)
    case(VADOSE_MODE)
      call pflow_vadose_setupini(realization)
#endif      
    case default
      call VecGetArrayF90(field%pressure,pressure_p,ierr)
      call VecGetArrayF90(field%temp,temp_p,ierr)
      call VecGetArrayF90(field%sat,sat_p,ierr)
      if (option%ndof == 3) call VecGetArrayF90(field%conc,conc_p,ierr)
      if (option%ndof == 4) call VecGetArrayF90(field%xmol,xmol_p,ierr)
    
      initial_condition => realization%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit
        
        do icell=1,initial_condition%region%num_cells
          local_id = initial_condition%region%cell_ids(icell)
          jn1 = 1+local_id*option%nphase-1
          jn2 = jn1+1
              
          count = 1
          pressure_p(jn1) = initial_condition%condition%cur_value(1)
          count = count + 1
          
          if (option%nphase>1) then
            pressure_p(jn2) = initial_condition%condition%cur_value(count)
            count = count + 1
          endif
          
          temp_p(local_id) = initial_condition%condition%cur_value(count)
          count = count + 1
          
          sat_p(jn1) = initial_condition%condition%cur_value(count)
          count = count + 1          
          
          if (option%nphase>1) then
            sat_p(jn2) = 1.d0 - sat_p(jn1)
            count = count + 1
          endif
          
          if (option%ndof == 3) then
            conc_p(local_id) = initial_condition%condition%cur_value(count)
            count = count + 1
          endif

          if (option%ndof == 4) then
            xmol_p(jn2) = initial_condition%condition%cur_value(count)
            count = count + 1
          endif
               
        enddo
  
        initial_condition => initial_condition%next
  
      enddo   
       
      call VecRestoreArrayF90(field%pressure,pressure_p,ierr)
      call VecRestoreArrayF90(field%temp,temp_p,ierr)
      call VecRestoreArrayF90(field%sat,sat_p,ierr)
      if (option%ndof == 3) call VecRestoreArrayF90(field%conc,conc_p,ierr)
      if (option%ndof == 4) call VecRestoreArrayF90(field%xmol,xmol_p,ierr)
  end select 

  ! assign initial conditions values to domain
  call VecGetArrayF90(field%xx,xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
  
  xx_p = -1.d20
  
  initial_condition => realization%initial_conditions%first
  do
  
    if (.not.associated(initial_condition)) exit

    if (.not.associated(initial_condition%connection)) then
      do icell=1,initial_condition%region%num_cells
        local_id = initial_condition%region%cell_ids(icell)
        ghosted_id = realization%grid%nL2G(local_id)
        iend = local_id*option%ndof
        ibegin = iend-option%ndof+1
        xx_p(ibegin:iend) = &
          initial_condition%condition%cur_value(1:option%ndof)
        iphase_loc_p(ghosted_id)=initial_condition%condition%iphase
      enddo
    else
      do iconn=1,initial_condition%connection%num_connections
        local_id = initial_condition%connection%id_dn(iconn)
        ghosted_id = realization%grid%nL2G(local_id)
        iend = local_id*option%ndof
        ibegin = iend-option%ndof+1
        xx_p(ibegin:iend) = &
          initial_condition%aux_real_var(1:option%ndof,iconn)
        iphase_loc_p(ghosted_id)=initial_condition%aux_int_var(1,iconn)
      enddo
    endif
    initial_condition => initial_condition%next
  enddo
  
  call VecRestoreArrayF90(field%xx,xx_p, ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)



#if 1 
  ! needs to be implemented
  ! set hydrostatic properties for initial and boundary conditions with depth
  if (option%ihydrostatic == 1) then
    select case(option%imode)
      case(MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE)
        call mhydrostatic(realization)
      case(OWG_MODE)
!        call owghydrostatic(realization)
      case default
!        call hydrostatic(realization)
    end select
  endif
#endif  
  !if (grid%iread_init==2) call Read_init_field(grid)
  
  if (option%imode == TWOPH_MODE) then
    field%pressurebc0(2,:) = field%pressurebc0(1,:)
    field%velocitybc0(2,:) = field%velocitybc0(1,:)
  endif

#if 0
  select case(option%imode)
    case(MPH_MODE,VADOSE_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE)
      deallocate(field%xx_ini)
      deallocate(option%iphas_ini)
    case default
      deallocate(option%pres_ini)
      deallocate(field%temp_ini)
      deallocate(field%sat_ini)
      deallocate(field%xmol_ini)
      deallocate(field%conc_ini)
  end select
#endif
end subroutine assignInitialConditions

! ************************************************************************** !
!
! initializeSolidReaction: Allocates and initializes arrays associated with
!                          mineral reactions
! author: Glenn Hammond
! date: 11/15/07
!
! ************************************************************************** !
subroutine initializeSolidReaction(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  integer :: icell
  PetscScalar, pointer :: phis_p(:)
  PetscErrorCode :: ierr
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  if (option%rk > 0.d0) then
    allocate(option%area_var(grid%nlmax))
    allocate(option%rate(grid%nlmax))
    call VecGetArrayF90(field%phis,phis_p,ierr)
    do icell = 1, grid%nlmax
      phis_p(icell) = option%phis0
      option%area_var(icell) = 1.d0
    enddo
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif
  
end subroutine initializeSolidReaction

end module Init_module
