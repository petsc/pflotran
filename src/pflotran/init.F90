module Init_module

  implicit none

  private

#include "definitions.h"

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscsnes.h"


  public :: Init

contains

! ************************************************************************** !
!
! Init: Initializes a pflow grid object
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine Init(simulation,filename)

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Discretization_module
  use Realization_module
  use Material_module
  use Timestepper_module
  use Field_module
  use Connection_module
  use Coupler_module
  use General_Grid_module
  use Debug_module
  use Convergence_module
  use Waypoint_module
  use Patch_module
  use Mass_Balance_module
  use Logging_module  
  use Database_module
  
  use MPHASE_module
  use Richards_module
  use THC_module
  
  use Reactive_Transport_module

  use Utility_module
    
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  type(realization_type), pointer :: realization
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(pflow_debug_type), pointer :: debug
  type(waypoint_list_type), pointer :: waypoint_list
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: temp_int
  PetscErrorCode :: ierr

  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr)
  call PetscLogEventBegin(logging%event_init,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr)
  
  ! set pointers to objects
  flow_stepper => simulation%flow_stepper
  tran_stepper => simulation%tran_stepper
  realization => simulation%realization
  discretization => realization%discretization
  option => realization%option
  field => realization%field
  debug => realization%debug
  
  nullify(flow_solver)
  nullify(tran_solver)
  
  ! read MODE,GRID,PROC,COMP,PHAS cards
  call readRequiredCardsFromInput(realization,filename)

  patch => realization%patch

  if(associated(patch%grid)) then
    grid => patch%grid
  endif

  ! process command line options
  call OptionCheckCommandLine(option)

  waypoint_list => WaypointListCreate()
  realization%waypoints => waypoint_list
  
  ! initialize flow mode
  if (len_trim(option%flowmode) > 0) then
    ! set the operational mode (e.g. THC_MODE, MPH_MODE, etc)
    call setFlowMode(option)
    flow_solver => flow_stepper%solver
  else
    call TimestepperDestroy(simulation%flow_stepper)
    nullify(flow_stepper)
  endif
  
  ! initialize transport mode
  if (option%ntrandof > 0) then
    tran_solver => tran_stepper%solver
  else
    call TimestepperDestroy(simulation%tran_stepper)
    nullify(tran_stepper)
  endif

  ! read in the remainder of the input file
  call readInput(simulation,filename)
  
  ! read reaction database
  if (associated(realization%reaction)) then
    if (option%ncmplx > 0 .or. &
        option%nmnrl > 0) then
      call DatabaseRead(realization%reaction,option)
      call BasisInit(realization%reaction,option)    
    endif
  endif

  ! create grid and allocate vectors
  call RealizationCreateDiscretization(realization)
  if (option%compute_mass_balance) then
    call MassBalanceCreate(realization)
  endif  
  
  if (option%myrank == 0) then
    ! general print statements for both flow and transport modes
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    if (realization%discretization%itype == STRUCTURED_GRID) then
      write(*,'(" number of processors = ",i5,", npx,y,z= ",3i5)') &
        option%commsize,grid%structured_grid%npx,grid%structured_grid%npy, &
        grid%structured_grid%npz
    endif
  endif

  ! update flow mode based on optional input
  if (option%nflowdof > 0) then
  
    if (flow_solver%J_mat_type == MATAIJ) then
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE)
          call printErrMsg(option,&
                           'AIJ matrix not supported for current mode: '// &
                           option%flowmode)
      end select
    endif

    if (option%myrank == 0) then
      write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
        option%nflowdof,option%nphase
      select case(option%iflowmode)
        case(MPH_MODE)
          write(*,'(" mode = MPH: p, T, s/C")')
        case(THC_MODE)
          write(*,'(" mode = Richards: p, T, s/C")')
        case(RICHARDS_MODE)
          write(*,'(" mode = Richards: p")')      
      end select
    endif

    call printMsg(option,"  Beginning set up of FLOW SNES ")

    call SolverCreateSNES(flow_solver,option%comm)  
    call SNESSetOptionsPrefix(flow_solver%snes, "flow_", ierr)
    call SolverCheckCommandLine(flow_solver)

    if (flow_solver%Jpre_mat_type == '') then
      if (flow_solver%J_mat_type /= MATMFFD) then
        flow_solver%Jpre_mat_type = flow_solver%J_mat_type
      else
        flow_solver%Jpre_mat_type = MATBAIJ
      endif
    endif

    call DiscretizationCreateJacobian(discretization,NFLOWDOF, &
                                      flow_solver%Jpre_mat_type, &
                                      flow_solver%Jpre, &
                                      option)
    call MatSetOptionsPrefix(flow_solver%Jpre, "flow_", ierr)

    if (flow_solver%J_mat_type /= MATMFFD) then
      flow_solver%J = flow_solver%Jpre
    endif

    if (flow_solver%use_galerkin_mg) then
      call DiscretizationCreateInterpolation(discretization,NFLOWDOF, &
                                             flow_solver%interpolation, &
                                             flow_solver%galerkin_mg_levels_x, &
                                             flow_solver%galerkin_mg_levels_y, &
                                             flow_solver%galerkin_mg_levels_z, &
                                             option)
    endif
    
    select case(option%iflowmode)
      case(THC_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,THCResidual, &
                             realization,ierr)
      case(RICHARDS_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,RichardsResidual, &
                             realization,ierr)
      case(MPH_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,MPHASEResidual, &
                             realization,ierr)
    end select
    
    if (flow_solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(flow_solver%snes,flow_solver%J,ierr)
    endif

    select case(option%iflowmode)
      case(THC_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             THCJacobian,realization,ierr)
      case(RICHARDS_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             RichardsJacobian,realization,ierr)
      case(MPH_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             MPHASEJacobian,realization, ierr)
    end select
    
    call SolverSetSNESOptions(flow_solver)

    string = 'Solver: ' // trim(flow_solver%ksp_type)
    call printMsg(option,string)
    string = 'Preconditioner: ' // trim(flow_solver%pc_type)
    call printMsg(option,string)

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    flow_stepper%convergence_context => ConvergenceContextCreate(flow_solver,option)
    call SNESSetConvergenceTest(flow_solver%snes,ConvergenceTest, &
                                flow_stepper%convergence_context, &
                                PETSC_NULL_FUNCTION,ierr) 
    call printMsg(option,"  Finished setting up FLOW SNES ")

  endif

  
  ! update transport mode based on optional input
  if (option%ntrandof > 0) then

    call printMsg(option,"  Beginning set up of TRAN SNES ")
    
    call SolverCreateSNES(tran_solver,option%comm)  
    call SNESSetOptionsPrefix(tran_solver%snes, "tran_", ierr)
    call SolverCheckCommandLine(tran_solver)

    if (tran_solver%Jpre_mat_type == '') then
      if (tran_solver%J_mat_type /= MATMFFD) then
        tran_solver%Jpre_mat_type = tran_solver%J_mat_type
      else
        tran_solver%Jpre_mat_type = MATBAIJ
      endif
    endif

    call DiscretizationCreateJacobian(discretization,NTRANDOF, &
                                      tran_solver%Jpre_mat_type, &
                                      tran_solver%Jpre,option)

    if (tran_solver%J_mat_type /= MATMFFD) then
      tran_solver%J = tran_solver%Jpre
    endif
    
    call MatSetOptionsPrefix(tran_solver%Jpre, "tran_", ierr)
    
    if (tran_solver%use_galerkin_mg) then
      call DiscretizationCreateInterpolation(discretization,NTRANDOF, &
                                             tran_solver%interpolation, &
                                             tran_solver%galerkin_mg_levels_x, &
                                             tran_solver%galerkin_mg_levels_y, &
                                             tran_solver%galerkin_mg_levels_z, &
                                             option)
    endif

    call SNESSetFunction(tran_solver%snes,field%tran_r,RTResidual,realization,ierr)

    if (tran_solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(tran_solver%snes,tran_solver%J,ierr)
    endif
    
    call SNESSetJacobian(tran_solver%snes,tran_solver%J,tran_solver%Jpre, &
                         RTJacobian,realization, ierr)

    call SNESLineSearchSet(tran_solver%snes,SNESLineSearchNo,PETSC_NULL_OBJECT,ierr)

    call SolverSetSNESOptions(tran_solver)

    string = 'Solver: ' // trim(tran_solver%ksp_type)
    call printMsg(option,string)
    string = 'Preconditioner: ' // trim(tran_solver%pc_type)
    call printMsg(option,string)

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    tran_stepper%convergence_context => ConvergenceContextCreate(tran_solver,option)
    call SNESSetConvergenceTest(tran_solver%snes,ConvergenceTest, &
                                tran_stepper%convergence_context, &
                                PETSC_NULL_FUNCTION,ierr) 

    call printMsg(option,"  Finished setting up TRAN SNES ")
  
  endif

  if (option%myrank == 0) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')


  call PetscLogEventBegin(logging%event_setup,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,ierr)
  ! read any regions provided in external files
  call readRegionFiles(realization)
  ! clip regions and set up boundary connectivity, distance  
  call RealizationLocalizeRegions(realization)
  ! link conditions with regions through couplers and generate connectivity
  call RealizationProcessCouplers(realization)
  call RealizationProcessConditions(realization)
  call assignMaterialPropToRegions(realization)
  call RealizationInitAllCouplerAuxVars(realization)

  ! should we still support this
  if (option%use_generalized_grid) then 
    if (option%myrank == 0) print *, 'Reading structured grid from hdf5' 
    if (.not.associated(patch%imat)) &
      allocate(patch%imat(grid%ngmax))  ! allocate material id array
    call ReadStructuredGridHDF5(realization)
  endif
  call PetscLogEventEnd(logging%event_setup,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,ierr)

  ! add waypoints associated with boundary conditions, source/sinks etc. to list
  call RealizationAddWaypointsToList(realization)
  ! fill in holes in waypoint data
  call WaypointListFillIn(option,realization%waypoints)
  call WaypointListRemoveExtraWaypnts(option,realization%waypoints)
  ! convert times from in put time to seconds
  call WaypointConvertTimes(realization%waypoints,realization%output_option%tconv)

  if (associated(flow_stepper)) flow_stepper%cur_waypoint => realization%waypoints%first
  if (associated(tran_stepper)) tran_stepper%cur_waypoint => realization%waypoints%first
  
  ! initialize FLOW
  ! set up auxillary variable arrays
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(THC_MODE)
        call THCSetup(realization)
      case(RICHARDS_MODE)
        call RichardsSetup(realization)
      case(MPH_MODE)
        call MphaseSetup(realization)
    end select
  
    ! assign initial conditions
    call RealizAssignFlowInitCond(realization)
  
    select case(option%iflowmode)
      case(THC_MODE)
        call THCUpdateAuxVars(realization)
      case(RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case(MPH_MODE)
        call MphaseUpdateAuxVars(realization)
    end select
  endif

  if (option%ntrandof > 0) then
    call RTSetup(realization)

    if (dabs(option%uniform_velocity(1)) + dabs(option%uniform_velocity(2)) + &
        dabs(option%uniform_velocity(3)) >  0.d0) then
      call RealizAssignUniformVelocity(realization)
    endif
    
    ! initialize densities and saturations
    if (option%nflowdof > 0) then
      call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
      call RealizationGetDataset(realization,global_vec,LIQUID_SATURATION,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(realization%discretization, &
                                       global_vec,field%saturation_loc,ONEDOF)   
      call RealizationGetDataset(realization,global_vec,LIQUID_DENSITY,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(realization%discretization, &
                                       global_vec,field%density_loc,ONEDOF)   
      call VecDestroy(global_vec,ierr)
    else
      call VecSet(field%saturation_loc,1.d0,ierr)
      call VecCopy(field%saturation_loc,field%saturation0_loc,ierr)
      call VecSet(field%density_loc,997.160290931658d0,ierr)
      call VecCopy(field%density_loc,field%density0_loc,ierr)
    endif

    ! map densities and saturations to reactive transport aux vars
    call RealizBridgeFlowAndTransport(realization) 
    ! initial concentrations must be assigned after densities are set !!!
    call RealizAssignTransportInitCond(realization)
    call RTUpdateAuxVars(realization)

  endif
  
  ! print info
  if (associated(flow_stepper)) then
    string = 'Flow Stepper:'
    call TimestepperPrintInfo(flow_stepper,option%fid_out,string,option)
  endif    
  if (associated(tran_stepper)) then
    string = 'Transport Stepper:'
    call TimestepperPrintInfo(tran_stepper,option%fid_out,string,option)
  endif    
  if (associated(flow_solver)) then
    string = 'Flow Newton Solver:'
    call SolverPrintNewtonInfo(flow_solver,option%fid_out,string,option%myrank)
  endif    
  if (associated(tran_solver)) then
    string = 'Transport Newton Solver:'
    call SolverPrintNewtonInfo(tran_solver,option%fid_out,string,option%myrank)
  endif    
  if (associated(flow_solver)) then
    string = 'Flow Linear Solver:'
    call SolverPrintLinearInfo(flow_solver,option%fid_out,string,option%myrank)
  endif    
  if (associated(tran_solver)) then
    string = 'Transport Linear Solver'
    call SolverPrintLinearInfo(tran_solver,option%fid_out,string,option%myrank)
  endif    

  if (debug%print_couplers) then
    call verifyAllCouplers(realization)
  endif
  
  call printMsg(option," ")
  call printMsg(option,"  Finished Initialization")

  call PetscLogEventEnd(logging%event_init,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,ierr)

end subroutine Init

! ************************************************************************** !
!
! readRequiredCardsFromInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readRequiredCardsFromInput(realization,filename)

  use Option_module
  use Discretization_module
  use Grid_module
  use Fileio_module
  use Patch_module
  use Level_module
  use Realization_module
  use AMR_Grid_module

  use Reaction_module  
  use Reaction_Aux_module  

  implicit none

  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename
  PetscInt ::  idum, i

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
  
  type(patch_type), pointer :: patch 
  type(level_type), pointer :: level
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization
  
  option%fid_in = IUNIT1
  option%fid_out = IUNIT2
  open(option%fid_in, file=filename, action="read", status="old") 
  open(option%fid_out, file='pflotran.out', action="write", status="unknown")

! we initialize the word to blanks to avoid error reported by valgrind
  do i=1,MAXWORDLENGTH
    word(i:i) = ' '
  enddo


! Read in select required cards
!.........................................................................

  ! MODE information
  string = "MODE"
  call fiFindStringInFile(option%fid_in,string,ierr)

  if (ierr == 0) then  
    ! strip card from front of string
    call fiReadWord(string,word,PETSC_FALSE,ierr)
 
    ! read in keyword 
    call fiReadWord(string,option%flowmode,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'flowmode','mode',ierr)
  endif

!.........................................................................

  ! GRID information
  string = "GRID"
  call fiFindStringInFile(option%fid_in,string,ierr)
  call fiFindStringErrorMsg(option%myrank,string,ierr)

  call DiscretizationRead(discretization,option%fid_in,PETSC_TRUE, option)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(realization%level_list)) then
        realization%level_list => LevelCreateList()
      endif
      level => LevelCreate()
      call LevelAddToList(level,realization%level_list)
      level%patch_list => PatchCreateList()
      call PatchAddToList(patch,level%patch_list)
      realization%patch => patch
    case(AMR_GRID)
      realization%level_list => AMRGridCreateLevelPatchLists(discretization%amrgrid)
      realization%patch => realization%level_list%first%patch_list%first
  end select
!.........................................................................

  if (realization%discretization%itype == STRUCTURED_GRID) then  ! look for processor decomposition
    
    ! PROC information
    string = "PROC"
    call fiFindStringInFile(option%fid_in,string,ierr)

    if (ierr == 0) then

      ! strip card from front of string
      call fiReadWord(string,word,PETSC_FALSE,ierr)
      call fiReadInt(string,grid%structured_grid%npx,ierr)
      call fiDefaultMsg(option%myrank,'npx',ierr)
      call fiReadInt(string,grid%structured_grid%npy,ierr)
      call fiDefaultMsg(option%myrank,'npy',ierr)
      call fiReadInt(string,grid%structured_grid%npz,ierr)
      call fiDefaultMsg(option%myrank,'npz',ierr)
 
      if (option%myrank == 0) &
        write(option%fid_out,'(/," *PROC",/, &
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
  string = "CHEMISTRY"
  call fiFindStringInFile(option%fid_in,string,ierr)

  if (ierr == 0) then
    realization%reaction => ReactionCreate()
    call ReactionRead(realization%reaction,option%fid_in,option)
    option%ntrandof = GetPrimarySpeciesCount(realization%reaction)
    option%comp_names => GetPrimarySpeciesNames(realization%reaction)
    option%ncomp = option%ntrandof
    option%ncmplx = GetSecondarySpeciesCount(realization%reaction)
    option%ngas = GetGasCount(realization%reaction)
!    option%gas_names => GetGasNames(realization%reaction)
    option%nmnrl = GetMineralCount(realization%reaction)
    option%mnrl_names => GetMineralNames(realization%reaction)
    option%nsorb = realization%reaction%neqsurfcmplx + &
                   realization%reaction%neqionx
    realization%reaction%ncomp = option%ntrandof
  endif

!.........................................................................

  ! COMP information
  string = "COMP"
  call fiFindStringInFile(option%fid_in,string,ierr)

  if (ierr == 0) then
    ! enter src here
  endif          

!....................................................................

  ! COMP information
  string = "PHAS"
  call fiFindStringInFile(option%fid_in,string,ierr)

  if (ierr == 0) then
    ! enter src here
  endif

!....................................................................

  ! TRAN information
  string = "TRAN"
  call fiFindStringInFile(option%fid_in,string,ierr)

  if (ierr == 0) then

    ! strip card from front of string
    call fiReadWord(string,word,PETSC_FALSE,ierr)

    call fiReadInt(string,option%ntrandof,ierr)
    call fiDefaultMsg(option%myrank,'ntrandof',ierr)
    option%ncomp = option%ntrandof

  endif          

    
end subroutine readRequiredCardsFromInput

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
  use AMR_Grid_module
  use Solver_module
  use Material_module
  use Fileio_module
  use Realization_module
  use Timestepper_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module
  use Breakthrough_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Discretization_module
 
  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXWORDLENGTH) :: filename

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXWORDLENGTH) :: card
    
  PetscReal, parameter:: fmwnacl = 58.44277D0, fmwh2o  = 18.01534d0
  PetscInt :: i, i1, i2, idum, ireg, isrc, j
  PetscInt :: ibc, ibrk, ir,np
  PetscReal :: rdum

  logical :: continuation_flag
  logical :: periodic_output_flag = PETSC_FALSE
  PetscReal :: periodic_rate = 0.d0
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscInt :: temp_int
  PetscInt :: length 
  PetscInt :: count, id
  
! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
!           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR

  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(tran_constraint_type), pointer :: tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(breakthrough_type), pointer :: breakthrough
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_type), pointer :: material
  type(thermal_property_type), pointer :: thermal_property
  type(saturation_function_type), pointer :: saturation_function

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch   
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  type(solver_type), pointer :: master_solver
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: master_stepper
  
  
  nullify(flow_stepper)
  nullify(tran_stepper)
  nullify(flow_solver)
  nullify(tran_solver)
  
  realization => simulation%realization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  tran_stepper => simulation%tran_stepper
  if (associated(tran_stepper)) tran_solver => tran_stepper%solver
  flow_stepper => simulation%flow_stepper
  if (associated(flow_stepper)) flow_solver => flow_stepper%solver

  if (associated(flow_stepper)) then
    master_stepper => flow_stepper
    master_solver => flow_solver
  else
    master_stepper => tran_stepper
    master_solver => tran_solver
  endif

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
                              
  rewind(option%fid_in)  
    
  do
    call fiReadFlotranString(option%fid_in, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,PETSC_FALSE,ierr)
    length = len_trim(word)
    call fiCharsToUpper(word,length)
!    call fiReadCard(word,card,ierr)
    card = trim(word)

    call printMsg(option,'pflotran card:: '//trim(card))

    select case(trim(card))

!....................
      case ('MODE')

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,option%fid_in,PETSC_FALSE,option)

!....................
      case ('CHEMISTRY')
        do
          call fiReadFlotranString(option%fid_in,string,ierr)
          call fiReadStringErrorMsg(option%myrank,card,ierr)
          if (fiCheckExit(string)) exit
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'word','CHEMISTRY',ierr) 
          select case(trim(word))
            case('PRIMARY_SPECIES','SECONDARY_SPECIES','GAS_SPECIES', &
                 'MINERALS')
              call fiSkipToEND(option%fid_in,option%myrank,card)
            case('MINERAL_KINETICS')
              call ReactionReadMineralKinetics(realization%reaction,option%fid_in,option)
            case('SORPTION')
              do
                call fiReadFlotranString(option%fid_in,string,ierr)
                call fiReadStringErrorMsg(option%myrank,card,ierr)
                if (fiCheckExit(string)) exit
                call fiReadWord(string,word,PETSC_TRUE,ierr)
                call fiErrorMsg(option%myrank,'word','CHEMISTRY,SORPTION',ierr) 
                select case(trim(word))
                  case('SURFACE_COMPLEXATION_RXN')
                    do
                      call fiReadFlotranString(option%fid_in,string,ierr)
                      call fiReadStringErrorMsg(option%myrank,card,ierr)
                      if (fiCheckExit(string)) exit
                      call fiReadWord(string,word,PETSC_TRUE,ierr)
                      call fiErrorMsg(option%myrank,'word','CHEMISTRY,SORPTION',ierr)
                      select case(trim(word))
                        case('COMPLEXES')
                          call fiSkipToEND(option%fid_in,option%myrank,card)
                      end select 
                    enddo
                  case('ION_EXCHANGE_RXN')
                  case('DISTRIBUTION_COEF')
                end select
              enddo
          end select
        enddo

!....................
      case ('TRAN')

!....................
      case ('UNIFORM_VELOCITY')
        call fiReadDouble(string,option%uniform_velocity(1),ierr)
        call fiErrorMsg(option%myrank,'velx','UNIFORM_VELOCITY', ierr)
        call fiReadDouble(string,option%uniform_velocity(2),ierr)
        call fiErrorMsg(option%myrank,'vely','UNIFORM_VELOCITY', ierr)
        call fiReadDouble(string,option%uniform_velocity(3),ierr)
        call fiErrorMsg(option%myrank,'velz','UNIFORM_VELOCITY', ierr)
      
!....................
      case ('DEBUG','PFLOW_DEBUG')
        call DebugRead(realization%debug,option%fid_in,option%myrank)
        
!....................
      case ('GENERALIZED_GRID')
        option%use_generalized_grid = PETSC_TRUE
        call fiReadWord(string,option%generalized_grid,PETSC_TRUE,ierr)

!....................
      case ('PROC')
      
!....................
      case ('REGION')
        region => RegionCreate()
        call fiReadWord(string,region%name,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'name','REGION',ierr) 
        call printMsg(option,region%name)
        call RegionRead(region,option%fid_in,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%regions)   
        nullify(region)   

!....................
      case ('FLOW_CONDITION')
        flow_condition => ConditionCreate(option)
        call fiReadWord(string,flow_condition%name,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'FLOW_CONDITION','name',ierr) 
        call printMsg(option,flow_condition%name)
        call ConditionRead(flow_condition,option,option%fid_in)
        call ConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)
        
!....................
      case ('TRANSPORT_CONDITION')
        tran_condition => TranConditionCreate(option)
        call fiReadWord(string,tran_condition%name,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'TRANSPORT_CONDITION','name',ierr) 
        call printMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition,realization%transport_constraints, &
                               option,option%fid_in)
        call TranConditionAddToList(tran_condition,realization%transport_conditions)
        nullify(tran_condition)

!....................
      case('CONSTRAINT')
        tran_constraint => TranConstraintCreate(option)
        call fiReadWord(string,tran_constraint%name,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'constraint','name',ierr) 
        call printMsg(option,tran_constraint%name)
        call TranConstraintRead(tran_constraint,option,option%fid_in)
        call TranConstraintAddToList(tran_constraint,realization%transport_constraints)
        nullify(tran_constraint)

!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call CouplerRead(coupler,option%fid_in,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call CouplerRead(coupler,option%fid_in,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call CouplerRead(coupler,option%fid_in,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,option%fid_in,option)
        call RealizationAddStrata(realization,strata)
        nullify(strata)
      
!.....................
      case ('DATASET') 
        call fiReadWord(string,word,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'dataset','name',ierr) 
        call printMsg(option,word)
        length = len_trim(word)
        call fiCharsToLower(word,length)        
        select case(word)
          case('permx')
            call fiReadWord(string,option%permx_filename,PETSC_TRUE,ierr)
            call fiErrorMsg(option%myrank,'dataset','permx_filename',ierr) 
          case('permy')
            call fiReadWord(string,option%permy_filename,PETSC_TRUE,ierr)
            call fiErrorMsg(option%myrank,'dataset','permy_filename',ierr) 
          case('permz')
            call fiReadWord(string,option%permz_filename,PETSC_TRUE,ierr)
            call fiErrorMsg(option%myrank,'dataset','permz_filename',ierr) 
        end select          
        
!.....................
      case ('COMP') 
        call fiSkipToEND(option%fid_in,option%myrank,card)
        
!.....................
      case ('PHAS')
        call fiSkipToEND(option%fid_in,option%myrank,card)
      
!....................

      case ('COUP')

        call printWrnMsg(option,"COUP not currently supported")
        call fiReadStringErrorMsg(option%myrank,'COUP',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'isync',ierr)

        if (option%myrank == 0) &
          write(option%fid_out,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') idum

!....................

      case ('GRAV','GRAVITY')

        call fiReadStringErrorMsg(option%myrank,'GRAV',ierr)

        call fiReadDouble(string,temp_real,ierr)
        if (ierr /= 0) then
          call fiDefaultMsg(option%myrank,'gravity',ierr)
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
          write(option%fid_out,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,3pe12.4 &
            & )') option%gravity(1:3)

!....................

      case ('HDF5')
        realization%output_option%print_hdf5 = PETSC_TRUE
        do
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          if (ierr /= 0) exit
          length = len_trim(word)
          call fiCharsToUpper(word,length)
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_hdf5_velocities = PETSC_TRUE
            case('FLUX')
              realization%output_option%print_hdf5_flux_velocities = PETSC_TRUE
            case default
          end select
            
        enddo

        if (option%myrank == 0) &
          write(option%fid_out,'(/," *HDF5",10x,l1,/)') realization%output_option%print_hdf5

!.....................
      case ('INVERT_Z','INVERTZ')
        if (associated(grid%structured_grid)) then
          grid%structured_grid%invert_z_axis = PETSC_TRUE
          option%gravity(3) = -1.d0*option%gravity(3)
        endif
      
!....................

      case ('TECP')
        realization%output_option%print_tecplot = PETSC_TRUE
        do
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          if (ierr /= 0) exit
          length = len_trim(word)
          call fiCharsToUpper(word,length)
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              realization%output_option%print_tecplot_velocities = PETSC_TRUE
            case('FLUX')
              realization%output_option%print_tecplot_flux_velocities = PETSC_TRUE
            case default
          end select
          
        enddo

        if (option%myrank == 0) &
          write(option%fid_out,'(/," *TECP",10x,l1,/)') realization%output_option%print_tecplot

!....................

      case ('IMOD')
        call fiReadInt(string,option%imod,ierr)
        call fiDefaultMsg(option%myrank,'mod',ierr)

!....................

      case ('TOLR')

        call fiReadStringErrorMsg(option%myrank,'TOLR',ierr)

        call fiReadInt(string,master_stepper%nstepmax,ierr)
        call fiDefaultMsg(option%myrank,'nstepmax',ierr)
  
        call fiReadInt(string,master_stepper%iaccel,ierr)
        call fiDefaultMsg(option%myrank,'iaccel',ierr)

        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg(option%myrank,'newton_max',ierr)

        call fiReadInt(string,master_stepper%icut_max,ierr)
        call fiDefaultMsg(option%myrank,'icut_max',ierr)

        call fiReadDouble(string,option%dpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dpmxe',ierr)

        call fiReadDouble(string,option%dtmpmxe,ierr)
        call fiDefaultMsg(option%myrank,'dtmpmxe',ierr)
  
        call fiReadDouble(string,option%dcmxe,ierr)
        call fiDefaultMsg(option%myrank,'dcmxe',ierr)

        call fiReadDouble(string,option%dsmxe,ierr)
        call fiDefaultMsg(option%myrank,'dsmxe',ierr)
        
        if (associated(tran_stepper)) then
          tran_stepper%icut_max = master_stepper%icut_max
        endif

        idum = 0
        if (option%myrank==0) write(option%fid_out,'(/," *TOLR ",/, &
          &"  steps  = ",i6,/,      &
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
          master_stepper%nstepmax,master_stepper%iaccel, &
          idum,master_stepper%icut_max, &
          option%dpmxe,option%dtmpmxe,option%dcmxe, option%dsmxe

!....................

#if 0
      case ('DXYZ')
      
        if (realization%discretization%itype == STRUCTURED_GRID) then  ! look for processor decomposition
          call StructuredGridReadDXYZ(grid%structured_grid,option)
        else if(realization%discretization%itype == AMR_GRID) then
          call AMRGridReadDXYZ(realization%discretization%amrgrid,option)
        else
          if (option%myrank == 0) &
            print *, 'ERROR: Keyword "DXYZ" not supported for unstructured grid'
            stop
        endif
#endif

!....................

#if 0
      case('ORIG','ORIGIN')
        call fiReadDouble(string,realization%discretization%origin(X_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'X direction','Origin',ierr)
        call fiReadDouble(string,realization%discretization%origin(Y_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'Y direction','Origin',ierr)
        call fiReadDouble(string,realization%discretization%origin(Z_DIRECTION),ierr)
        call fiErrorMsg(option%myrank,'Z direction','Origin',ierr)
#endif
        
!....................

      case ('DIFF')

        call fiReadStringErrorMsg(option%myrank,'DIFF',ierr)

        call fiReadDouble(string,option%difaq,ierr)
        call fiDefaultMsg(option%myrank,'difaq',ierr)

        call fiReadDouble(string,option%delhaq,ierr)
        call fiDefaultMsg(option%myrank,'delhaq',ierr)

        if (option%myrank==0) write(option%fid_out,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          option%difaq,option%delhaq

!....................

      case ('RCTR')

        call printErrMsg(option,"RCTR currently out of date.  Needs to be reimplemented")
#if 0
        call fiReadStringErrorMsg(option%myrank,'RCTR',ierr)

        call fiReadInt(string,option%ityprxn,ierr)
        call fiDefaultMsg(option%myrank,'ityprxn',ierr)

        call fiReadDouble(string,option%rk,ierr)
        call fiDefaultMsg(option%myrank,'rk',ierr)

        call fiReadDouble(string,option%phis0,ierr)
        call fiDefaultMsg(option%myrank,'phis0',ierr)

        call fiReadDouble(string,option%areas0,ierr)
        call fiDefaultMsg(option%myrank,'areas0',ierr)

        call fiReadDouble(string,option%pwrsrf,ierr)
        call fiDefaultMsg(option%myrank,'pwrsrf',ierr)

        call fiReadDouble(string,option%vbars,ierr)
        call fiDefaultMsg(option%myrank,'vbars',ierr)

        call fiReadDouble(string,option%ceq,ierr)
        call fiDefaultMsg(option%myrank,'ceq',ierr)

        call fiReadDouble(string,option%delHs,ierr)
        call fiDefaultMsg(option%myrank,'delHs',ierr)

        call fiReadDouble(string,option%delEs,ierr)
        call fiDefaultMsg(option%myrank,'delEs',ierr)

        call fiReadDouble(string,option%wfmts,ierr)
        call fiDefaultMsg(option%myrank,'wfmts',ierr)

        if (option%myrank == 0) &
        write(option%fid_out,'(/," *RCTR",/, &
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
#endif
!....................

      case ('RADN')

        call fiReadStringErrorMsg(option%myrank,'RADN',ierr)

        call fiReadDouble(string,option%ret,ierr)
        call fiDefaultMsg(option%myrank,'ret',ierr)

        call fiReadDouble(string,option%fc,ierr)
        call fiDefaultMsg(option%myrank,'fc',ierr)

        if (option%myrank==0) write(option%fid_out,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          option%ret,option%fc

!....................


      case ('PHAR')
        call printWrnMsg(option,"PHAR currently out of date.  Needs to be reimplemented")
#if 0
        call fiReadStringErrorMsg(option%myrank,'PHAR',ierr)

        call fiReadDouble(string,option%qu_kin,ierr)
        call fiDefaultMsg(option%myrank,'TransReaction',ierr)
        if (option%myrank==0) write(option%fid_out,'(/," *PHAR ",1pe12.4)')option%qu_kin
        option%yh2o_in_co2 = 0.d0
        if (option%qu_kin > 0.d0) option%yh2o_in_co2 = 1.d-2 ! check this number!
#endif     
!......................

      case('REFERENCE_PRESSURE')
        call fiReadStringErrorMsg(option%myrank,'REFERENCE_PRESSURE',ierr)
        call fiReadDouble(string,option%pref,ierr)
        call fiDefaultMsg(option%myrank,'Reference Pressure',ierr) 

!......................

      case('BRIN','BRINE')
        call fiReadStringErrorMsg(option%myrank,'BRIN',ierr)
        call fiReadDouble(string,option%m_nacl,ierr)
        call fiDefaultMsg(option%myrank,'NaCl Concentration',ierr) 

        call fiReadWord(string,word,PETSC_FALSE,ierr)
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

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call fiReadWord(string,option%restart_file,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'RESTART','Restart file name',ierr) 
        call fiReadDouble(string,option%restart_time,ierr)
        call fiDefaultMsg(option%myrank,'Restart time',ierr) 

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call fiReadInt(string,option%checkpoint_frequency,ierr)
        call fiErrorMsg(option%myrank,'CHECKPOINT','Checkpoint frequency',ierr) 

!......................

      case ('NUMERICAL_JACOBIAN')
        option%numerical_derivatives = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS','STATISTICS')
        option%compute_statistics = PETSC_TRUE

      case ('COMPUTE_MASS_BALANCE','MASS_BALANCE')
        option%compute_mass_balance = PETSC_TRUE

!....................

      case ('TIMESTEPPER')
        call fiReadWord(string,word,PETSC_FALSE,ierr)
        length = len_trim(word)
        call fiCharsToUpper(word,length)
        select case(word)
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call TimestepperRead(tran_stepper,option%fid_in,option)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
          case default
            if (associated(flow_solver)) then
              call TimestepperRead(flow_stepper,option%fid_in,option)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
        end select

!....................

      case ('LINEAR_SOLVER')
        call fiReadWord(string,word,PETSC_FALSE,ierr)
        length = len_trim(word)
        call fiCharsToUpper(word,length)
        select case(word)
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call SolverReadLinear(tran_solver,option%fid_in,option%myrank)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
          case default
            if (associated(flow_solver)) then
              call SolverReadLinear(flow_solver,option%fid_in,option%myrank)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
        end select

!....................

      case ('NEWTON_SOLVER')
        call fiReadWord(string,word,PETSC_FALSE,ierr)
        length = len_trim(word)
        call fiCharsToUpper(word,length)
        select case(word)
          case('TRAN','TRANSPORT')
            if (associated(tran_stepper)) then
              call SolverReadNewton(tran_solver,option%fid_in,option%myrank)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
          case default
            if (associated(flow_stepper)) then
              call SolverReadNewton(flow_solver,option%fid_in,option%myrank)
            else
              call fiSkipToEND(option%fid_in,option%myrank,card)
            endif
        end select

!....................

      case ('FLUID_PROPERTY','FLUID_PROPERTIES')

        realization%fluid_properties => FluidPropertyCreate(option%nphase)
        
        count = 0
        do
          call fiReadFlotranString(option%fid_in,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'FLUID_PROPERTIES',ierr)
          
          if (fiCheckExit(string)) exit
         
          count = count + 1 
          if (count > option%nphase) exit              
                        
          call fiReadDouble(string,realization%fluid_properties%diff_base(count),ierr)
          call fiErrorMsg(option%myrank,'diff_base','FLUID_PROPERTIES', ierr)          
        
          call fiReadDouble(string,realization%fluid_properties%diff_exp(count),ierr)
          call fiErrorMsg(option%myrank,'diff_exp','FLUID_PROPERTIES', ierr)          

        enddo
        
!....................

      case ('THRM','THERMAL_PROPERTY','THERMAL_PROPERTIES')

        count = 0
        do
          call fiReadFlotranString(option%fid_in,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'THRM',ierr)

          if (fiCheckExit(string)) exit
       
          count = count + 1
          thermal_property => ThermalPropertyCreate()
      
          call fiReadInt(string,thermal_property%id,ierr)
          call fiErrorMsg(option%myrank,'id','THRM', ierr)

          call fiReadDouble(string,thermal_property%rock_density,ierr)
          call fiErrorMsg(option%myrank,'rock density','THRM', ierr)

          call fiReadDouble(string,thermal_property%spec_heat,ierr)
          call fiErrorMsg(option%myrank,'cpr','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_dry,ierr)
          call fiErrorMsg(option%myrank,'ckdry','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%therm_cond_wet,ierr)
          call fiErrorMsg(option%myrank,'ckwet','THRM', ierr)
        
          call fiReadDouble(string,thermal_property%tort_bin_diff,ierr)
          call fiErrorMsg(option%myrank,'tau','THRM', ierr)

          call fiReadDouble(string,thermal_property%vap_air_diff_coef,ierr)
          call fiErrorMsg(option%myrank,'cdiff','THRM', ierr)

          call fiReadDouble(string,thermal_property%exp_binary_diff,ierr)
          call fiErrorMsg(option%myrank,'cexp','THRM', ierr)

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
          write(option%fid_out,'(/," *THRM: ",i3)') count
          write(option%fid_out,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(option%fid_out,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, count
            write(option%fid_out,'(i4,1p7e11.4)') i,option%rock_density(i), &
            option%cpr(i),option%ckdry(i),option%ckwet(i), &
            option%tau(i),option%cdiff(i),option%cexp(i)
          enddo
        endif

!....................

      case ('PCKR','SATURATION_FUNCTION','SATURATION_FUNCTIONS')
      
        count = 0
        do
          call fiReadFlotranString(option%fid_in,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'PCKR',ierr)

          if (fiCheckExit(string)) exit
       
          count = count + 1
          saturation_function => SaturationFunctionCreate(option)
          
          call fiReadInt(string,saturation_function%id,ierr)
          call fiErrorMsg(option%myrank,'id','PCKR', ierr)
          
          call fiReadInt(string,saturation_function%saturation_function_itype,ierr)
          call fiErrorMsg(option%myrank,'icaptype','PCKR', ierr)
      
          select case(option%iflowmode)
            case(MPH_MODE,THC_MODE,RICHARDS_MODE)
              do np=1, option%nphase
                call fiReadDouble(string,saturation_function%Sr(np),ierr)
                call fiErrorMsg(option%myrank,'Sr','PCKR', ierr)
              enddo 
            case default
              call fiReadDouble(string,saturation_function%Sr(1),ierr)
              call fiErrorMsg(option%myrank,'Sr','PCKR', ierr)
          end select
        
          call fiReadDouble(string,saturation_function%m,ierr)
          call fiErrorMsg(option%myrank,'pckrm','PCKR', ierr)
          saturation_function%lambda = saturation_function%m

          call fiReadDouble(string,saturation_function%alpha,ierr)
          call fiErrorMsg(option%myrank,'alpha','PCKR', ierr)

          call fiReadDouble(string,saturation_function%pcwmax,ierr)
          call fiErrorMsg(option%myrank,'pcwmax','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%betac,ierr)
          call fiErrorMsg(option%myrank,'pbetac','PCKR', ierr)
      
          call fiReadDouble(string,saturation_function%power,ierr)
          call fiErrorMsg(option%myrank,'pwrprm','PCKR', ierr)

          call SaturationFunctionComputeSpline(option,saturation_function)
          
          call SaturationFunctionAddToList(saturation_function, &
                                           realization%saturation_functions)

        enddo
        
        ! allocate dynamic arrays holding saturation function information
        allocate(option%icaptype(count))
        option%icaptype = 0
  
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,RICHARDS_MODE)
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
          
          option%icaptype(id) = saturation_function%saturation_function_itype
          select case(option%iflowmode)
            case(MPH_MODE,THC_MODE,RICHARDS_MODE)
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

      !clu removed on 05/21/08
#if 0
        if (option%iflowmode == MPH_MODE .or. &
            option%iflowmode == THC_MODE .or. &
            option%iflowmode == RICHARDS_MODE) then
          call pckr_init(option%nphase,count,grid%nlmax, &
                         option%icaptype,option%sir, option%pckrm, &
                         option%lambda,option%alpha,option%pcwmax, &
                         option%pcbetac,option%pwrprm)
        endif 
#endif
      
        if (option%myrank==0) then
          write(option%fid_out,'(/," *PCKR: ",i3)') count
          write(option%fid_out,'("  icp swir    lambda         alpha")')
          do j = 1, count
            if (option%iflowmode == MPH_MODE .or. &
                option%iflowmode == THC_MODE .or. &
                option%iflowmode == RICHARDS_MODE) then
              write(option%fid_out,'(i4,1p8e12.4)') option%icaptype(j),(option%sir(np,j),np=1, &
                option%nphase),option%lambda(j),option%alpha(j), &
                option%pcwmax(j),option%pcbetac(j),option%pwrprm(j)
            else
              write(option%fid_out,'(i4,1p7e12.4)') option%icaptype(j),option%swir(j), &
                option%lambda(j),option%alpha(j),option%pcwmax(j), &
                option%pcbetac(j),option%pwrprm(j)
            endif
          enddo
        end if

        if (option%iflowmode == MPH_MODE .or. &
            option%iflowmode == THC_MODE .or. &
            option%iflowmode == RICHARDS_MODE) then
          deallocate(option%icaptype, option%pckrm, option%lambda, &
                     option%alpha,option%pcwmax, option%pcbetac, &
                     option%pwrprm)
        endif 
 
        call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                           realization%saturation_function_array)
        
!....................
      
      case ('PHIK','MATERIAL','MATERIALS')

        count = 0
        do
          call fiReadFlotranString(option%fid_in,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'PHIK',ierr)

          if (fiCheckExit(string)) exit
       
          count = count + 1
          material => MaterialCreate()

          call fiReadWord(string,material%name,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'name','PHIK', ierr)
                
          call fiReadInt(string,material%id,ierr)
          call fiErrorMsg(option%myrank,'id','PHIK', ierr)
                
          call fiReadInt(string,material%icap,ierr)
          call fiErrorMsg(option%myrank,'icap','PHIK', ierr)
  
          call fiReadInt(string,material%ithrm,ierr)
          call fiErrorMsg(option%myrank,'ithrm','PHIK', ierr)
  
          call fiReadDouble(string,material%porosity,ierr)
          call fiErrorMsg(option%myrank,'por','PHIK', ierr)
          
          call fiReadDouble(string,material%tortuosity,ierr)
          call fiErrorMsg(option%myrank,'tor','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(1,1),ierr)
          call fiErrorMsg(option%myrank,'permx','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(2,2),ierr)
          call fiErrorMsg(option%myrank,'permy','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability(3,3),ierr)
          call fiErrorMsg(option%myrank,'permz','PHIK', ierr)
  
          call fiReadDouble(string,material%permeability_pwr,ierr)
          call fiErrorMsg(option%myrank,'permpwr','PHIK', ierr)
          
          material%permeability(1:3,1:3) = material%permeability(1:3,1:3)
          
          call MaterialAddToList(material,realization%materials)
          
        enddo          

        call MaterialConvertListToArray(realization%materials, &
                                        realization%material_array)
                                        
!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
        string = '-viewer_binary_mpiio'
!        call PetscOptionsInsertString(string,ierr)

      case ('HANDSHAKE_IO')
        call fiReadInt(string,option%io_handshake_buffer_size,ierr)
        call fiErrorMsg(option%myrank,'io_handshake_buffer_size','HANDSHAKE_IO', ierr)

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%overwrite_restart_transport = PETSC_TRUE

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%overwrite_restart_flow_params = PETSC_TRUE

      case ('TIME')

        call fiReadStringErrorMsg(option%myrank,'TIME',ierr)
      
        call fiReadWord(string,word,PETSC_FALSE,ierr)
      
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


        call fiReadWord(string,word,PETSC_FALSE,ierr)
        if (ierr == 0) then
          call fiWordToUpper(word)
          if (fiStringCompare(word,'EVERY',FIVE_INTEGER)) then
            periodic_output_flag = PETSC_TRUE
            call fiReadDouble(string,periodic_rate,ierr)
          endif
        endif

        continuation_flag = PETSC_TRUE
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(option%fid_in,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = PETSC_FALSE
          if (index(string,backslash) > 0) continuation_flag = PETSC_TRUE
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              waypoint%print_output = PETSC_TRUE              
              call WaypointInsertInList(waypoint,realization%waypoints)
            endif
          enddo
        enddo
        
        ! make last waypoint final
        waypoint%final = PETSC_TRUE
        
        if (periodic_output_flag) then
          temp_real2 = waypoint%time   ! final simulation time
          temp_real = periodic_rate
          do
            if (temp_real > temp_real2) exit
            
            waypoint => WaypointCreate()
            waypoint%time = temp_real
            waypoint%print_output = PETSC_TRUE              
            call WaypointInsertInList(waypoint,realization%waypoints)
            
            temp_real = temp_real + periodic_rate
          enddo
        endif

!....................

      case ('DTST')

        call fiReadStringErrorMsg(option%myrank,'DTST',ierr)

        call fiReadDouble(string,master_stepper%dt_min,ierr)
        call fiDefaultMsg(option%myrank,'dt_min',ierr)
            
        continuation_flag = PETSC_TRUE
        temp_int = 0       
        do
          if (.not.continuation_flag) exit
          call fiReadFlotranString(option%fid_in,string,ierr)
          if (ierr /= 0) exit
          continuation_flag = PETSC_FALSE
          if (index(string,backslash) > 0) continuation_flag = PETSC_TRUE
          ierr = 0
          do
            if (ierr /= 0) exit
            call fiReadDouble(string,temp_real,ierr)
            if (ierr == 0) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              call fiReadDouble(string,waypoint%dt_max,ierr)
              call fiErrorMsg(option%myrank,'dt_max','dtst',ierr)
              if (temp_int == 0) master_stepper%dt_max = waypoint%dt_max
              call WaypointInsertInList(waypoint,realization%waypoints)
              temp_int = temp_int + 1
            endif
          enddo
        enddo
        
        master_stepper%dt_min = realization%output_option%tconv * master_stepper%dt_min
        master_stepper%dt_max = realization%output_option%tconv * master_stepper%dt_max

        option%flow_dt = master_stepper%dt_min
        option%tran_dt = master_stepper%dt_min
      
!....................
      case ('BRK','BREAKTHROUGH')
        breakthrough => BreakthroughCreate()
        call BreakthroughRead(breakthrough,option%fid_in,option)
        call RealizationAddBreakthrough(realization,breakthrough)        
      
!....................
      case('SDST')
        print *, 'SDST needs to be implemented'
        stop
#if 0
! Needs implementation         
        allocate(master_stepper%steady_eps(option%nflowdof))
        do j=1,option%nflowdof
          call fiReadDouble(string,master_stepper%steady_eps(j),ierr)
          call fiDefaultMsg(option%myrank,'steady tol',ierr)
        enddo
        if (option%myrank==0) write(option%fid_out,'(/," *SDST ",/, &
          &"  dpdt        = ",1pe12.4,/, &
          &"  dtmpdt        = ",1pe12.4,/, &
          &"  dcdt        = ",1pe12.4)') &
          master_stepper%steady_eps
#endif

!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call fiReadDouble(string,option%wallclock_stop_time,ierr)
        call fiErrorMsg(option%myrank,'stop time','WALLCLOCK_STOP', ierr) 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*3600.d0
      
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

  close(option%fid_in)
  
end subroutine readInput

! ************************************************************************** !
!
! setFlowMode: Sets the flow mode (richards, vadose, mph, etc.)
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine setFlowMode(option)

  use Option_module
  use Fileio_module

  implicit none 

  type(option_type) :: option
  PetscInt :: length
  
  length = len_trim(option%flowmode)
  call fiCharsToUpper(option%flowmode,length)
  select case(option%flowmode)
    case('THC')
      option%iflowmode = THC_MODE
      option%nphase = 1
      option%nflowdof = 3
      option%nflowspec = 2
    case('RICHARDS')
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%nflowdof = 1
      option%nflowspec = 1
    case('MPH','MPHASE')
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
    case default
      call printErrMsg(option,'Mode: '//trim(option%flowmode)//' not recognized.')
  end select
  
end subroutine setFlowMode

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
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Fileio_module
  use Patch_module
  use Level_module

  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: por0_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_pow_p(:)
  PetscReal, pointer :: tor_loc_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, material_id
  PetscInt :: istart, iend
  PetscInt :: fid = 86, status
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: patch  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  type(material_type), pointer :: material
  type(region_type), pointer :: region
  
  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  ! loop over all strata to determine if any are inactive or
  ! have associated cell by cell material ids
  cur_level => realization%level_list%first
  do 
     if (.not.associated(cur_level)) exit
     cur_patch => cur_level%patch_list%first
     do
        if (.not.associated(cur_patch)) exit
        strata => cur_patch%strata%first
        do
           if (.not.associated(strata)) exit
           if (.not.strata%active .or. .not.associated(strata%region)) then
              if (.not.associated(cur_patch%imat)) then
                 allocate(cur_patch%imat(cur_patch%grid%ngmax))
                 cur_patch%imat = -999
              endif
              exit
           endif
           strata => strata%next
        enddo
        cur_patch => cur_patch%next
     enddo
     cur_level => cur_level%next
  enddo

  ! Read in cell by cell material ids if they exist
  cur_level => realization%level_list%first
  do 
     if (.not.associated(cur_level)) exit
     cur_patch => cur_level%patch_list%first
     do
        if (.not.associated(cur_patch)) exit
        strata => cur_patch%strata%first
        do
           if (.not.associated(strata)) exit
           if (.not.associated(strata%region) .and. strata%active) then
              call readMaterialsFromFile(realization,strata%material_name)
           endif
           strata => strata%next
        end do
        cur_patch => cur_patch%next
     enddo
     cur_level => cur_level%next
  enddo
    
  cur_level => realization%level_list%first
  do 
     if (.not.associated(cur_level)) exit
     cur_patch => cur_level%patch_list%first
     do
        if (.not.associated(cur_patch)) exit

        grid => cur_patch%grid

        if (option%nflowdof > 0) then
           call GridVecGetArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
           call GridVecGetArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)
           call GridVecGetArrayF90(grid,field%perm0_xx,perm_xx_p,ierr)
           call GridVecGetArrayF90(grid,field%perm0_yy,perm_yy_p,ierr)
           call GridVecGetArrayF90(grid,field%perm0_zz,perm_zz_p,ierr)
           call GridVecGetArrayF90(grid,field%perm_pow,perm_pow_p,ierr)
        endif
        call GridVecGetArrayF90(grid,field%porosity0,por0_p,ierr)
        call GridVecGetArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
        
        strata => cur_patch%strata%first
        do
           if (.not.associated(strata)) exit
           
           if (strata%active) then
              region => strata%region
              material => strata%material
              if (associated(region)) then
                 istart = 1
                 iend = region%num_cells
              else
                 istart = 1
                 iend = grid%nlmax
              endif
              do icell=istart, iend
                 if (associated(region)) then
                    local_id = region%cell_ids(icell)
                 else
                    local_id = icell
                 endif

                 ghosted_id = grid%nL2G(local_id)
                 if (associated(cur_patch%imat)) then
                    ! if patch%imat is allocated and the id > 0, the material id 
                    ! supercedes the material pointer for the strata
                    material_id = cur_patch%imat(ghosted_id)
                    if (material_id > 0 .and. &
                         material_id <= size(realization%material_array)) then
                       material => realization%material_array(material_id)%ptr
                    endif
                    ! otherwide set the imat value to the stratas material
                    if (material_id < -998) & ! prevent overwrite of cell already set to inactive
                         cur_patch%imat(ghosted_id) = material%id
                 endif
                 if (associated(material)) then
                    if (option%nflowdof > 0) then
                       icap_loc_p(ghosted_id) = material%icap
                       ithrm_loc_p(ghosted_id) = material%ithrm
                       perm_xx_p(local_id) = material%permeability(1,1)
                       perm_yy_p(local_id) = material%permeability(2,2)
                       perm_zz_p(local_id) = material%permeability(3,3)
                       perm_pow_p(local_id) = material%permeability_pwr
                    endif
                    por0_p(local_id) = material%porosity
                    tor_loc_p(ghosted_id) = material%tortuosity
                 endif
              enddo
           endif
           strata => strata%next
        enddo

        if (option%nflowdof > 0) then
           call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
           call GridVecRestoreArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)
           call GridVecRestoreArrayF90(grid,field%perm0_xx,perm_xx_p,ierr)
           call GridVecRestoreArrayF90(grid,field%perm0_yy,perm_yy_p,ierr)
           call GridVecRestoreArrayF90(grid,field%perm0_zz,perm_zz_p,ierr)
           call GridVecRestoreArrayF90(grid,field%perm_pow,perm_pow_p,ierr)
        endif
        call GridVecRestoreArrayF90(grid,field%porosity0,por0_p,ierr)
        call GridVecRestoreArrayF90(grid,field%tor_loc,tor_loc_p,ierr)
        
        ! read in any cell by cell data 
        if (len_trim(option%permx_filename) > 1) then
           call readVectorFromFile(realization,field%perm0_xx, &
                option%permx_filename,GLOBAL)  
        endif
        if (len_trim(option%permy_filename) > 1) then
           call readVectorFromFile(realization,field%perm0_yy, &
                option%permy_filename,GLOBAL)  
        endif
        if (len_trim(option%permz_filename) > 1) then
           call readVectorFromFile(realization,field%perm0_zz, &
                option%permz_filename,GLOBAL)  
        endif
        
        cur_patch => cur_patch%next
     enddo
     cur_level => cur_level%next
  enddo

  if (option%nflowdof > 0) then
     call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
          field%perm_xx_loc,ONEDOF)  
     call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
          field%perm_yy_loc,ONEDOF)  
     call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
          field%perm_zz_loc,ONEDOF)   
     
     call DiscretizationLocalToLocal(discretization,field%icap_loc, &
          field%icap_loc,ONEDOF)   
     call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
          field%ithrm_loc,ONEDOF)
  endif
  
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
       field%porosity_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%tor_loc, &
       field%tor_loc,ONEDOF)   

end subroutine assignMaterialPropToRegions

#if 0
! ************************************************************************** !
!
! assignInitialConditions: Assigns initial conditions to regions of realization
! author: Glenn Hammond
! date: 09/15/08
!
! ************************************************************************** !
subroutine assignInitialConditions(realization)

  use Realization_module
  use Level_module
  use Patch_module
  use Discretization_module
  use Grid_module
  use Field_module
  use Coupler_module
  use Option_module

  implicit none

  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof
  PetscInt :: local_id, ghosted_id
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  
  ! first call initialization of primary dofs
  call RealizAssignInitialConditions(realization)
  
  ! now, let's fill in the secondary variables (mineral vol fracs, etc.)
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      
      if (option%nflowdof > 0) then
      endif

      if (option%ntrandof > 0) then
      
        initial_condition => cur_patch%initial_conditions%first
        do
        
          if (.not.associated(initial_condition)) exit

          if (.not.associated(initial_condition%connection_set)) then
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              if (associated(cur_patch%imat)) then
                if (cur_patch%imat(ghosted_id) <= 0) then
                  cycle
                endif
              endif
              ! minerals              
              if (associated(initial_condition%tran_condition%cur_constraint_coupler%minerals)) then
                do idof = 1, option%nmnrl
                  cur_patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(idof) = &
                    initial_condition%tran_condition%cur_constraint_coupler% &
                      minerals%basis_mol_frac(idof)
                enddo
              endif
            enddo
          else
            do iconn=1,initial_condition%connection_set%num_connections
              local_id = initial_condition%connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (associated(cur_patch%imat)) then
                if (cur_patch%imat(ghosted_id) <= 0) then
                  cycle
                endif
              endif
              ! minerals 
              if (associated(initial_condition%tran_condition%cur_constraint_coupler%minerals)) then
                do idof = 1, option%nmnrl
                  cur_patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(idof) = &
                    initial_condition%tran_condition%cur_constraint_coupler% &
                      minerals%basis_mol_frac(idof)
                enddo
              endif
            enddo
          endif
          initial_condition => initial_condition%next
        enddo
        
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine assignInitialConditions  
#endif

! ************************************************************************** !
!
! assignUniformVelocity: Assigns uniform velocity in connection list
!                        darcy velocities
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine assignUniformVelocity(realization)

  use Realization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Grid_module
  use Patch_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch   
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: vdarcy
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
    
  ! Internal Flux Terms -----------------------------------
  cur_connection_set => grid%internal_connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = OptionDotProduct(option%uniform_velocity, &
                                cur_connection_set%dist(1:3,iconn))
      patch%internal_velocities(1,sum_connection) = vdarcy
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = OptionDotProduct(option%uniform_velocity, &
                                cur_connection_set%dist(1:3,iconn))
      patch%boundary_velocities(1,sum_connection) = vdarcy
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine assignUniformVelocity

! ************************************************************************** !
!
! verifyAllCouplers: Verifies the connectivity of a coupler
! author: Glenn Hammond
! date: 1/8/08
!
! ************************************************************************** !
subroutine verifyAllCouplers(realization)

  use Realization_module
  use Level_module
  use Patch_module
  use Coupler_module

  implicit none

  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

        call verifyCoupler(realization,cur_patch,cur_patch%initial_conditions)
        call verifyCoupler(realization,cur_patch,cur_patch%boundary_conditions)
        call verifyCoupler(realization,cur_patch,cur_patch%source_sinks)

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine verifyAllCouplers

! ************************************************************************** !
!
! verifyCoupler: Verifies the connectivity of a coupler
! author: Glenn Hammond
! date: 1/8/08
!
! ************************************************************************** !
subroutine verifyCoupler(realization,patch,coupler_list)

  use Realization_module
  use Discretization_module
  use Option_module 
  use Coupler_module
  use Condition_module
  use Grid_module
  use Output_module
  use Patch_module

  implicit none

  type(realization_type) :: realization
  type(coupler_list_type), pointer :: coupler_list

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(coupler_type), pointer :: coupler
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: iconn, icell, local_id
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid

  if (.not.associated(coupler_list)) return

  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    call VecZeroEntries(global_vec,ierr)
    call VecGetArrayF90(global_vec,vec_ptr,ierr) 
    if (associated(coupler%connection_set)) then
      do iconn = 1, coupler%connection_set%num_connections
        local_id = coupler%connection_set%id_dn(iconn)
        vec_ptr(local_id) = coupler%id
      enddo
    else
      if (associated(coupler%region)) then
        do icell = 1, coupler%region%num_cells
          local_id = coupler%region%cell_ids(icell)
          vec_ptr(local_id) = coupler%id
        enddo
      endif
    endif
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr) 
    if (len_trim(coupler%flow_condition_name) > 0) then
      dataset_name = coupler%flow_condition_name
    elseif (len_trim(coupler%tran_condition_name) > 0) then
      dataset_name = coupler%tran_condition_name
    endif
    write(word,*) patch%id
    dataset_name = trim(dataset_name) // '_' // &
                   trim(coupler%region%name) // '_' // &
                   trim(adjustl(word))
    dataset_name = dataset_name(1:28)
    filename = trim(dataset_name) // '.tec'
    call OutputVectorTecplot(filename,dataset_name,realization,global_vec)

    coupler => coupler%next
  enddo

  call VecDestroy(global_vec,ierr)

end subroutine verifyCoupler

! ************************************************************************** !
!
! readRegionFiles: Reads in grid cell ids stored in files
! author: Glenn Hammond
! date: 1/03/08
!
! ************************************************************************** !
subroutine readRegionFiles(realization)

  use Realization_module
  use Region_module
  use HDF5_module

  implicit none

  type(realization_type) :: realization
  
  type(region_type), pointer :: region
  
  region => realization%regions%first
  do 
    if (.not.associated(region)) exit
    if (len_trim(region%filename) > 1) then
      if (index(region%filename,'.h5') > 0) then
        call HDF5ReadRegionFromFile(realization,region,region%filename)
      else
        call RegionReadFromFile(region,region%filename)
      endif
    endif
    region => region%next
  enddo

end subroutine readRegionFiles

! ************************************************************************** !
!
! readMaterialsFromFile: Reads in grid cell materials
! author: Glenn Hammond
! date: 1/03/08
!
! ************************************************************************** !
subroutine readMaterialsFromFile(realization,filename)

  use Realization_module
  use Field_module
  use Grid_module
  use Option_module
  use Fileio_module
  use Patch_module
  use Logging_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (index(filename,'.h5') > 0) then
    call HDF5ReadMaterialsFromFile(realization,filename)
  else
    call GridCreateNaturalToGhostedHash(grid,option)
    status = 0
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      string = "File: " // filename // " not found."
      call printErrMsg(option,string)
    endif
    call PetscLogEventBegin(logging%event_hash_map, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
    do
      call fiReadFlotranString(fid,string,ierr)
      if (ierr /= 0) exit
      call fiReadInt(string,natural_id,ierr)
      call fiErrorMsg(option%myrank,'natural id','STRATA', ierr)
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call fiReadInt(string,material_id,ierr)
        call fiErrorMsg(option%myrank,'material id','STRATA', ierr)
        patch%imat(ghosted_id) = material_id
      endif
    enddo
    call PetscLogEventEnd(logging%event_hash_map, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
    call GridDestroyHashTable(grid)
  endif
  
end subroutine readMaterialsFromFile

! ************************************************************************** !
!
! readVectorFromFile: Reads data from a file into an associated vector
! author: Glenn Hammond
! date: 03/18/08
!
! ************************************************************************** !
subroutine readVectorFromFile(realization,vector,filename,vector_type)

  use Realization_module
  use Discretization_module
  use Field_module
  use Grid_module
  use Option_module
  use Fileio_module
  use Patch_module
  use Logging_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  Vec :: vector
  character(len=MAXWORDLENGTH) :: filename
  PetscInt :: vector_type
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr
  PetscInt :: count, read_count, i
  PetscInt, pointer :: indices(:)
  PetscReal, pointer :: values(:)
  PetscInt, parameter :: block_size = 10000
  Vec :: natural_vec, global_vec
  PetscMPIInt :: source = 0

  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (index(filename,'.h5') > 0) then
    ! to be taken care of later
  else
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      string = "File: " // filename // " not found."
      call printErrMsg(option,string)
    endif
    allocate(values(block_size))
    allocate(indices(block_size))
    call DiscretizationCreateVector(discretization,ONEDOF,natural_vec, &
                                    NATURAL,option)
    count = 0
    do
      if (count >= grid%nmax) exit
      read_count = min(block_size,grid%nmax-count)
      do i=1,read_count
        indices(i) = count+i-1 ! zero-based indexing
      enddo
      ierr = 0
      if (option%myrank == 0) read(fid,*,iostat=ierr) values(1:read_count)
      call mpi_bcast(ierr,ONE_INTEGER,MPI_INTEGER,source,option%comm,ierr)      
      if (ierr /= 0) then
        string = 'Insufficent data in file: ' // filename
        call printErrMsg(option,string)
      endif
      if (option%myrank == 0) then
        call VecSetValues(natural_vec,read_count,indices,values,INSERT_VALUES, &
                          ierr)
      endif
      count = count + read_count
    enddo
    call mpi_bcast(count,ONE_INTEGER,MPI_INTEGER,source,option%comm,ierr)      
    if (count /= grid%nmax) then
      write(string,'(a,i8,a,i8,a)') 'Number of data in file (', count, &
                                    ') does not match size of vector (', &
                                    grid%nlmax, ')'
      call printErrMsg(option,string)
    endif
    close(fid)
    deallocate(values)
    nullify(values)
    deallocate(indices)
    nullify(indices)
    call VecAssemblyBegin(natural_vec,ierr)
    call VecAssemblyEnd(natural_vec,ierr)
    select case(vector_type)
      case(LOCAL)
        call DiscretizationCreateVector(discretization,ONEDOF,global_vec, &
                                        GLOBAL,option)        
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                           global_vec,ONEDOF)  
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         vector,ONEDOF)
        call VecDestroy(global_vec,ierr)  
      case(GLOBAL)
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                          vector,ONEDOF) 
    end select 
    call VecDestroy(natural_vec,ierr)
  endif
  
end subroutine readVectorFromFile

end module Init_module
