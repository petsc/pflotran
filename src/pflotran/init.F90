module Init_module

  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"

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
  use Input_module
  
  use MPHASE_module
  use Richards_module
  use THC_module
  
  use Reactive_Transport_module
  
  use Global_module

  use water_eos_module
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
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: temp_int
  PetscErrorCode :: ierr
  PCSide:: pcside
  PetscReal :: r1, r2, r3, r4, r5, r6

  interface

     subroutine SAMRInitializePreconditioner(p_application, which_pc, pc)
#include "finclude/petsc.h"
#include "finclude/petscpc.h"
       PC :: pc
       PetscFortranAddr :: p_application
       PetscInt :: which_pc
     end subroutine SAMRInitializePreconditioner

  end interface

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
  input => realization%input
  
  nullify(flow_solver)
  nullify(tran_solver)
  

  realization%input => InputCreate(IUNIT1,filename)
  open(option%fid_out, file='pflotran.out', action="write", status="unknown")

  ! read required cards
  call readRequiredCardsFromInput(realization)

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
    option%nphase = 1
    option%liquid_phase = 1
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
  call readInput(simulation)
  call InputDestroy(realization%input)

  ! initialize reference density
  call wateos(option%reference_temperature,option%reference_pressure, &
              option%reference_density,r1,r2,r3,r4,r5,r6, &
              option%scale,ierr)
  
  ! read reaction database
  if (associated(realization%reaction)) then
    if (realization%reaction%use_full_geochemistry) then
      call DatabaseRead(realization%reaction,option)
      call BasisInit(realization%reaction,option)    
    endif
  endif

  ! create grid and allocate vectors
  call RealizationCreateDiscretization(realization)
  if (option%compute_mass_balance) then
    call MassBalanceCreate(realization)
  endif  
  
  if (OptionPrint(option)) then
    ! general print statements for both flow and transport modes
    write(*,'(/,"++++++++++++++++++++++++++++++++++++++++++++++++++++&
      &++++++++")')
    if (realization%discretization%itype == STRUCTURED_GRID) then
      write(*,'(" Requested processors and decomposition = ",i5,", npx,y,z= ",3i5)') &
        option%commsize,grid%structured_grid%npx,grid%structured_grid%npy, &
        grid%structured_grid%npz
      write(*,'(" Actual decomposition: npx,y,z= ",3i5)') &
        grid%structured_grid%npx_final,grid%structured_grid%npy_final, &
        grid%structured_grid%npz_final
    endif
  endif

  ! update flow mode based on optional input
  if (option%nflowdof > 0) then
  
    if (flow_solver%J_mat_type == MATAIJ) then
      select case(option%iflowmode)
        case(MPH_MODE,THC_MODE)
          option%io_buffer = 'AIJ matrix not supported for current mode: '// &
                             option%flowmode
          call printErrMsg(option)
      end select
    endif

    if (OptionPrint(option)) then
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

    call printMsg(option,"  Beginning setup of FLOW SNES ")

    call SolverCreateSNES(flow_solver,option%comm)  
    call SNESSetOptionsPrefix(flow_solver%snes, "flow_",ierr)
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
    call MatSetOptionsPrefix(flow_solver%Jpre,"flow_",ierr)

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
                             MPHASEJacobian,realization,ierr)
    end select
    
    call SolverSetSNESOptions(flow_solver)

    ! If we are using a structured grid, set the corresponding flow DA 
    ! as the DA for the PCEXOTIC preconditioner, in case we choose to use it.
    ! The PCExoticSetDA() call is ignored if the PCEXOTIC preconditioner is 
    ! no used.  We need to put this call after SolverCreateSNES() so that 
    ! KSPSetFromOptions() will already have been called.
    ! I also note that this preconditioner is intended only for the flow, 
    ! solver.  --RTM
    if (realization%discretization%itype == STRUCTURED_GRID) then
      call PCExoticSetDA(flow_solver%pc, &
                         realization%discretization%dm_nflowdof,ierr);
    endif

    ! setup a shell preconditioner and initialize in the case of AMR
    if(associated(discretization%amrgrid)) then
!       flow_solver%pc_type = PCSHELL
       pcside = PC_RIGHT
       if(flow_solver%pc_type==PCSHELL) then
          call KSPSetPreconditionerSide(flow_solver%ksp, pcside,ierr)
          call SAMRInitializePreconditioner(discretization%amrgrid%p_application, 0, flow_solver%pc)
       endif
    endif

    option%io_buffer = 'Solver: ' // trim(flow_solver%ksp_type)
    call printMsg(option)
    option%io_buffer = 'Preconditioner: ' // trim(flow_solver%pc_type)
    call printMsg(option)

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

    call printMsg(option,"  Beginning setup of TRAN SNES ")
    
    call SolverCreateSNES(tran_solver,option%comm)  
    call SNESSetOptionsPrefix(tran_solver%snes, "tran_",ierr)
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
    
    call MatSetOptionsPrefix(tran_solver%Jpre,"tran_",ierr)
    
    if (tran_solver%use_galerkin_mg) then
      call DiscretizationCreateInterpolation(discretization,NTRANDOF, &
                                             tran_solver%interpolation, &
                                             tran_solver%galerkin_mg_levels_x, &
                                             tran_solver%galerkin_mg_levels_y, &
                                             tran_solver%galerkin_mg_levels_z, &
                                             option)
    endif

    call SNESSetFunction(tran_solver%snes,field%tran_r,RTResidual,&
                         realization,ierr)

    if (tran_solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(tran_solver%snes,tran_solver%J,ierr)
    endif
    
    call SNESSetJacobian(tran_solver%snes,tran_solver%J,tran_solver%Jpre, &
                         RTJacobian,realization,ierr)

    ! this could be changed in the future if there is a way to ensure that the linesearch
    ! update does not perturb concentrations negative.
    call SNESLineSearchSet(tran_solver%snes,SNESLineSearchNo, &
                           PETSC_NULL_OBJECT,ierr)

    call SolverSetSNESOptions(tran_solver)

    ! setup a shell preconditioner and initialize in the case of AMR
    if(associated(discretization%amrgrid)) then
!       flow_solver%pc_type = PCSHELL
       pcside = PC_RIGHT
       if(tran_solver%pc_type==PCSHELL) then
          call KSPSetPreconditionerSide(tran_solver%ksp, pcside,ierr)
          call SAMRInitializePreconditioner(discretization%amrgrid%p_application, 1, tran_solver%pc)
       endif
    endif
    option%io_buffer = 'Solver: ' // trim(tran_solver%ksp_type)
    call printMsg(option)
    option%io_buffer = 'Preconditioner: ' // trim(tran_solver%pc_type)
    call printMsg(option)

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    tran_stepper%convergence_context => ConvergenceContextCreate(tran_solver,option)
    call SNESSetConvergenceTest(tran_solver%snes,ConvergenceTest, &
                                tran_stepper%convergence_context, &
                                PETSC_NULL_FUNCTION,ierr) 

    ! this update check must be in place, otherwise reactive transport is likely
    ! to fail
    call SNESLineSearchSetPreCheck(tran_solver%snes, &
                                   ConvergenceRTUpdateCheck, &
                                   realization,ierr)

    call printMsg(option,"  Finished setting up TRAN SNES ")
  
  endif

  if (OptionPrint(option)) write(*,'("++++++++++++++++++++++++++++++++&
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
  if (option%ntrandof > 0) then
    call RealizationInitConstraints(realization)
  endif
  call RealizationPrintCouplers(realization)

  ! should we still support this
  if (option%use_generalized_grid) then 
    call printMsg(option,'Reading structured grid from hdf5')
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
  
  ! initialize global auxilliary variable object
  call GlobalSetup(realization)
  
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
  
    ! assign initial conditionsRealizAssignFlowInitCond
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

    if (option%nflowdof == 0) then
      call GlobalSetAuxVarScalar(realization,option%reference_pressure, &
                                 PRESSURE)
      call GlobalSetAuxVarScalar(realization,option%reference_temperature, &
                                 TEMPERATURE)
      call GlobalSetAuxVarScalar(realization,option%reference_saturation, &
                                 LIQUID_SATURATION)
      call GlobalSetAuxVarScalar(realization,option%reference_density, &
                                 LIQUID_DENSITY)
    endif

    ! initial concentrations must be assigned after densities are set !!!
    call RealizAssignTransportInitCond(realization)
    call RTUpdateAuxVars(realization)
    ! at this point the auxvars have been computed with activity coef = 1.d0
    ! to use intitial condition with activity coefs /= 1.d0, must update
    ! activity coefs and recompute auxvars
    if (realization%reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
      call RTUpdateSolution(realization)    
      call RTUpdateAuxVars(realization)
    endif
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
  if (OptionPrint(option) .and. associated(flow_solver)) then
    string = 'Flow Newton Solver:'
    call SolverPrintNewtonInfo(flow_solver,option%fid_out,string)
  endif    
  if (OptionPrint(option) .and. associated(tran_solver)) then
    string = 'Transport Newton Solver:'
    call SolverPrintNewtonInfo(tran_solver,option%fid_out,string)
  endif    
  if (OptionPrint(option) .and. associated(flow_solver)) then
    string = 'Flow Linear Solver:'
    call SolverPrintLinearInfo(flow_solver,option%fid_out,string)
  endif    
  if (OptionPrint(option) .and. associated(tran_solver)) then
    string = 'Transport Linear Solver'
    call SolverPrintLinearInfo(tran_solver,option%fid_out,string)
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
subroutine readRequiredCardsFromInput(realization)

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Level_module
  use Realization_module
  use AMR_Grid_module
  use Input_module

  use Reaction_module  
  use Reaction_Aux_module  

  implicit none

  type(realization_type) :: realization

  character(len=MAXSTRINGLENGTH) :: string
  
  type(patch_type), pointer :: patch 
  type(level_type), pointer :: level
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  
  patch => realization%patch
  option => realization%option
  discretization => realization%discretization
  
  input => realization%input
  
! Read in select required cards
!.........................................................................

  ! MODE information
  string = "MODE"
  call InputFindStringInFile(input,option,string)

  if (.not.InputError(input)) then  
    ! read in keyword 
    call InputReadWord(input,option,option%flowmode,PETSC_TRUE)
    call InputErrorMsg(input,option,'flowmode','mode')
  endif

!.........................................................................

  ! GRID information
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call DiscretizationRead(discretization,input,PETSC_TRUE,option)
  
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
    call InputFindStringInFile(input,option,string)

    if (.not.InputError(input)) then

      grid => realization%patch%grid
      ! strip card from front of string
      call InputReadInt(input,option,grid%structured_grid%npx)
      call InputDefaultMsg(input,option,'npx')
      call InputReadInt(input,option,grid%structured_grid%npy)
      call InputDefaultMsg(input,option,'npy')
      call InputReadInt(input,option,grid%structured_grid%npz)
      call InputDefaultMsg(input,option,'npz')
 
      if (option%myrank == option%io_rank) then
        option%io_buffer = ' Processor Decomposition:'
        call printMsg(option)
        write(option%io_buffer,'("  npx   = ",3x,i4)') grid%structured_grid%npx
        call printMsg(option)
        write(option%io_buffer,'("  npy   = ",3x,i4)') grid%structured_grid%npy
        call printMsg(option)
        write(option%io_buffer,'("  npz   = ",3x,i4)') grid%structured_grid%npz
        call printMsg(option)
      endif
  
      if (option%commsize /= grid%structured_grid%npx * &
                             grid%structured_grid%npy * &
                             grid%structured_grid%npz) then
        write(option%io_buffer,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%npx*grid%structured_grid%npy* &
                       grid%structured_grid%npz,' commsize = ',option%commsize
        call printErrMsg(option)
      endif
    endif
  endif
  
!.........................................................................

  ! CHEMISTRY information
  string = "CHEMISTRY"
  call InputFindStringInFile(input,option,string)

  if (.not.InputError(input)) then
    reaction => ReactionCreate()
    realization%reaction => reaction
    call ReactionRead(reaction,input,option)
    reaction%primary_species_names => GetPrimarySpeciesNames(reaction)
    option%ntrandof = GetPrimarySpeciesCount(reaction)
  endif
    
end subroutine readRequiredCardsFromInput

! ************************************************************************** !
!
! readInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine readInput(simulation)

  use Simulation_module
  use Option_module
  use Field_module
  use Grid_module
  use Structured_Grid_module
  use AMR_Grid_module
  use Solver_module
  use Material_module
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
  use Reaction_Aux_module
  use Discretization_module
  use Input_module
  use String_module
 
  implicit none
  
  type(simulation_type) :: simulation

  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
    
  PetscInt :: i, i1, i2, idum, ireg, isrc, j
  PetscInt :: ibc, ibrk, ir,np
  PetscReal :: rdum

  PetscTruth :: continuation_flag
  PetscTruth :: periodic_output_flag = PETSC_FALSE
  PetscReal :: periodic_rate = 0.d0
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscInt :: temp_int
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
  type(reaction_type), pointer :: reaction
  type(input_type), pointer :: input
  
  
  nullify(flow_stepper)
  nullify(tran_stepper)
  nullify(flow_solver)
  nullify(tran_solver)
  
  realization => simulation%realization
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  input => realization%input

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
                              
  rewind(input%fid)  
    
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))

!....................
      case ('MODE')

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input, &
                                PETSC_FALSE,option)

!....................
      case ('CHEMISTRY')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','CHEMISTRY') 
          select case(trim(word))
            case('PRIMARY_SPECIES','SECONDARY_SPECIES','GAS_SPECIES', &
                 'MINERALS')
              call InputSkipToEND(input,option,card)
            case('MINERAL_KINETICS')
              call ReactionReadMineralKinetics(reaction,input,option)
            case('SORPTION')
              do
                call InputReadFlotranString(input,option)
                call InputReadStringErrorMsg(input,option,card)
                if (InputCheckExit(input,option)) exit
                call InputReadWord(input,option,word,PETSC_TRUE)
                call InputErrorMsg(input,option,'SORPTION','CHEMISTRY') 
                select case(trim(word))
                  case('SURFACE_COMPLEXATION_RXN','ION_EXCHANGE_RXN')
                    do
                      call InputReadFlotranString(input,option)
                      call InputReadStringErrorMsg(input,option,card)
                      if (InputCheckExit(input,option)) exit
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'SORPTION','CHEMISTRY')
                      select case(trim(word))
                        case('COMPLEXES','CATIONS')
                          call InputSkipToEND(input,option,word)
                      end select 
                    enddo
                  case('DISTRIBUTION_COEF')
                end select
              enddo
          end select
        enddo

!....................
      case ('TRAN')

!....................
      case ('UNIFORM_VELOCITY')
        call InputReadDouble(input,option,option%uniform_velocity(1))
        call InputErrorMsg(input,option,'velx','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,option%uniform_velocity(2))
        call InputErrorMsg(input,option,'vely','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,option%uniform_velocity(3))
        call InputErrorMsg(input,option,'velz','UNIFORM_VELOCITY')
      
!....................
      case ('DEBUG','PFLOW_DEBUG')
        call DebugRead(realization%debug,input,option)
        
!....................
      case ('GENERALIZED_GRID')
        option%use_generalized_grid = PETSC_TRUE
        call InputReadWord(input,option,option%generalized_grid,PETSC_TRUE)

!....................
      case ('PROC')
      
!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION') 
        call printMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%regions)   
        nullify(region)   

!....................
      case ('FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'FLOW_CONDITION','name') 
        call printMsg(option,flow_condition%name)
        call FlowConditionRead(flow_condition,input,option)
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)
        
!....................
      case ('TRANSPORT_CONDITION')
        tran_condition => TranConditionCreate(option)
        call InputReadWord(input,option,tran_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'TRANSPORT_CONDITION','name') 
        call printMsg(option,tran_condition%name)
        call TranConditionRead(tran_condition,realization%transport_constraints, &
                               reaction,input,option)
        call TranConditionAddToList(tran_condition,realization%transport_conditions)
        nullify(tran_condition)

!....................
      case('CONSTRAINT')
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name') 
        call printMsg(option,tran_constraint%name)
        call TranConstraintRead(tran_constraint,reaction,input,option)
        call TranConstraintAddToList(tran_constraint,realization%transport_constraints)
        nullify(tran_constraint)

!....................
      case ('BOUNDARY_CONDITION')
        coupler => CouplerCreate(BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Boundary Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)
      
!....................
      case ('INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('SOURCE_SINK')
        coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name') 
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization,coupler)
        nullify(coupler)        
      
!....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization,strata)
        nullify(strata)
      
!.....................
      case ('DATASET') 
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'dataset','name') 
        call printMsg(option,word)
        call StringToLower(word)        
        select case(word)
          case('permx')
            call InputReadWord(input,option,option%permx_filename,PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset','permx_filename') 
          case('permy')
            call InputReadWord(input,option,option%permy_filename,PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset','permy_filename') 
          case('permz')
            call InputReadWord(input,option,option%permz_filename,PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset','permz_filename') 
        end select          
        
!.....................
      case ('COMP') 
        call InputSkipToEnd(input,option,card)
        
!.....................
      case ('PHAS')
        call InputSkipToEnd(input,option,card)
      
!....................

      case ('COUP')

        call printWrnMsg(option,"COUP not currently supported")
        call InputReadStringErrorMsg(input,option,'COUP')

        call InputReadInt(input,option,idum)
        call InputDefaultMsg(input,option,'isync')

        if (option%myrank == option%io_rank) &
          write(option%fid_out,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') idum

!....................

      case ('GRAV','GRAVITY')

        call InputReadStringErrorMsg(input,option,'GRAV')

        call InputReadDouble(input,option,temp_real)
        if (InputError(input)) then
          call InputDefaultMsg(input,option,'gravity')
        else
          call InputReadDouble(input,option,option%gravity(2))
          if (InputError(input)) then
            option%gravity(:) = 0.d0
            option%gravity(3) = temp_real
          else
            option%gravity(1) = temp_real
            call InputReadDouble(input,option,option%gravity(3))
          endif
        endif

        if (option%myrank == option%io_rank) &
          write(option%fid_out,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,3pe12.4 &
            & )') option%gravity(1:3)

!....................

      case ('HDF5')
        realization%output_option%print_hdf5 = PETSC_TRUE
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          call StringToUpper(word)

          select case(word)
            case('VELO')
              realization%output_option%print_hdf5_velocities = PETSC_TRUE
            case('FLUX')
              realization%output_option%print_hdf5_flux_velocities = PETSC_TRUE
            case default
          end select
            
        enddo

        if (option%myrank == option%io_rank) &
          write(option%fid_out,'(/," *HDF5",10x,l1,/)') &
            realization%output_option%print_hdf5

!.....................
      case ('INVERT_Z','INVERTZ')
        if (associated(grid%structured_grid)) then
          grid%structured_grid%invert_z_axis = PETSC_TRUE
          option%gravity(3) = -option%gravity(3)
        endif
      
!....................

      case ('TECP','TECPLOT')
        realization%output_option%print_tecplot = PETSC_TRUE

        if (StringCompare(card,'TECPLOT',SEVEN_INTEGER)) then
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'TECPLOT','type') 
          call StringToUpper(word)
          select case(trim(word))
            case('POINT')
              realization%output_option%tecplot_format = TECPLOT_POINT_FORMAT
            case('BLOCK')
              realization%output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
            case default
              option%io_buffer = 'TECPLOT format (' // trim(word) // &
                                 ') not recongnized.'
              call printErrMsg(option)
          end select
        else
          realization%output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
        endif
        
        ! safety catch: only block supported in parallel
        if (realization%output_option%tecplot_format == TECPLOT_POINT_FORMAT .and. &
            option%commsize > 1) then
          realization%output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
        endif
          
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          call StringToUpper(word)

          select case(word)
            case('VELO')
              realization%output_option%print_tecplot_velocities = PETSC_TRUE
            case('FLUX')
              if (realization%output_option%tecplot_format == TECPLOT_POINT_FORMAT) then
                option%io_buffer = &
                    'Printing of fluxes not supported in TECPLOT POINT format.'
                call printErrMsg(option)
              endif
              realization%output_option%print_tecplot_flux_velocities = PETSC_TRUE
            case default
          end select
          
        enddo

        if (option%myrank == option%io_rank) &
          write(option%fid_out,'(/," *TECP",10x,l1,/)') &
            realization%output_option%print_tecplot

      case ('VTK')
        realization%output_option%print_vtk = PETSC_TRUE

        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          call StringToUpper(word)

          select case(word)
            case('VELO')
              realization%output_option%print_vtk_velocities = PETSC_TRUE
            case('FLUX')
              if (realization%output_option%tecplot_format == TECPLOT_POINT_FORMAT) then
                call printErrMsg(option,'Printing of fluxes not supported &
                                         &in TECPLOT POINT format.')
              endif
              realization%output_option%print_tecplot_flux_velocities = PETSC_TRUE
            case default
          end select
          
        enddo

        if (option%myrank == option%io_rank) &
          write(option%fid_out,'(/," *VTK",10x,l1,/)') &
            realization%output_option%print_vtk

!....................

      case ('IMOD')
        call InputReadInt(input,option,option%imod)
        call InputDefaultMsg(input,option,'mod')

!....................

      case ('PRINT_ACT_COEFS')
        realization%output_option%print_act_coefs = PETSC_TRUE

!....................

      case ('TOLR')

        call InputReadStringErrorMsg(input,option,'TOLR')

        call InputReadInt(input,option,master_stepper%nstepmax)
        call InputDefaultMsg(input,option,'nstepmax')
  
        call InputReadInt(input,option,master_stepper%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

        call InputReadInt(input,option,idum)
        call InputDefaultMsg(input,option,'newton_max')

        call InputReadInt(input,option,master_stepper%icut_max)
        call InputDefaultMsg(input,option,'icut_max')

        call InputReadDouble(input,option,option%dpmxe)
        call InputDefaultMsg(input,option,'dpmxe')

        call InputReadDouble(input,option,option%dtmpmxe)
        call InputDefaultMsg(input,option,'dtmpmxe')
  
        call InputReadDouble(input,option,option%dcmxe)
        call InputDefaultMsg(input,option,'dcmxe')

        call InputReadDouble(input,option,option%dsmxe)
        call InputDefaultMsg(input,option,'dsmxe')
        
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
          option%io_buffer = 'Keyword "DXYZ" not supported for unstructured grid'
          call printErrMsg(option)
        endif
#endif

!....................

#if 0
      case('ORIG','ORIGIN')
        call InputReadDouble(input,option,realization%discretization%origin(X_DIRECTION))
        call InputErrorMsg(input,option,'X direction','Origin')
        call InputReadDouble(input,option,realization%discretization%origin(Y_DIRECTION))
        call InputErrorMsg(input,option,'Y direction','Origin')
        call InputReadDouble(input,option,realization%discretization%origin(Z_DIRECTION))
        call InputErrorMsg(input,option,'Z direction','Origin')
#endif
        
!....................

      case ('DIFF')

        call InputReadStringErrorMsg(input,option,'DIFF')

        call InputReadDouble(input,option,option%difaq)
        call InputDefaultMsg(input,option,'difaq')

        call InputReadDouble(input,option,option%delhaq)
        call InputDefaultMsg(input,option,'delhaq')

        if (option%myrank==0) write(option%fid_out,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          option%difaq,option%delhaq

!....................

      case ('RCTR')

        call printErrMsg(option,"RCTR currently out of date.  Needs to be reimplemented")
#if 0
        call InputReadStringErrorMsg(input,option,'RCTR')

        call InputReadInt(input,option,option%ityprxn)
        call InputDefaultMsg(input,option,'ityprxn')

        call InputReadDouble(input,option,option%rk)
        call InputDefaultMsg(input,option,'rk')

        call InputReadDouble(input,option,option%phis0)
        call InputDefaultMsg(input,option,'phis0')

        call InputReadDouble(input,option,option%areas0)
        call InputDefaultMsg(input,option,'areas0')

        call InputReadDouble(input,option,option%pwrsrf)
        call InputDefaultMsg(input,option,'pwrsrf')

        call InputReadDouble(input,option,option%vbars)
        call InputDefaultMsg(input,option,'vbars')

        call InputReadDouble(input,option,option%ceq)
        call InputDefaultMsg(input,option,'ceq')

        call InputReadDouble(input,option,option%delHs)
        call InputDefaultMsg(input,option,'delHs')

        call InputReadDouble(input,option,option%delEs)
        call InputDefaultMsg(input,option,'delEs')

        call InputReadDouble(input,option,option%wfmts)
        call InputDefaultMsg(input,option,'wfmts')

        if (OptionPrint(option)) &
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

        call InputReadStringErrorMsg(input,option,'RADN')

        call InputReadDouble(input,option,option%ret)
        call InputDefaultMsg(input,option,'ret')

        call InputReadDouble(input,option,option%fc)
        call InputDefaultMsg(input,option,'fc')

        if (option%myrank==0) write(option%fid_out,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          option%ret,option%fc

!....................


      case ('PHAR')
        call printWrnMsg(option,"PHAR currently out of date.  Needs to be reimplemented")
#if 0
        call InputReadStringErrorMsg(input,option,'PHAR')

        call InputReadDouble(input,option,option%qu_kin)
        call InputDefaultMsg(input,option,'TransReaction')
        if (option%myrank==0) write(option%fid_out,'(/," *PHAR ",1pe12.4)')option%qu_kin
        option%yh2o_in_co2 = 0.d0
        if (option%qu_kin > 0.d0) option%yh2o_in_co2 = 1.d-2 ! check this number!
#endif     
!......................

      case('MOLAL','MOLALITY')
        option%initialize_with_molality = PETSC_TRUE

!......................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_pressure)
        call InputDefaultMsg(input,option,'Reference Pressure') 

!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_temperature)
        call InputDefaultMsg(input,option,'Reference Temperature') 

!......................

      case('REFERENCE_POROSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_porosity)
        call InputDefaultMsg(input,option,'Reference Porosity') 

!......................

      case('REFERENCE_SATURATION')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_saturation)
        call InputDefaultMsg(input,option,'Reference Saturation') 

!......................

      case('BRIN','BRINE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%m_nacl)
        call InputDefaultMsg(input,option,'NaCl Concentration') 

        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            option%m_nacl = option%m_nacl /FMWNACL/(1.D0-option%m_nacl)
          case('MOLE')    
            option%m_nacl = option%m_nacl /FMWH2O/(1.D0-option%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (OptionPrint(option)) print *, option%m_nacl
!......................

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call InputReadWord(input,option,option%restart_file,PETSC_TRUE)
        call InputErrorMsg(input,option,'RESTART','Restart file name') 
        call InputReadDouble(input,option,option%restart_time)
        call InputDefaultMsg(input,option,'Restart time') 

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call InputReadInt(input,option,option%checkpoint_frequency)
        call InputErrorMsg(input,option,'CHECKPOINT','Checkpoint frequency') 

!......................

      case ('NUMERICAL_JACOBIAN')
        option%numerical_derivatives = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS','STATISTICS')
        option%compute_statistics = PETSC_TRUE

!      case ('COMPUTE_MASS_BALANCE','MASS_BALANCE')
!        option%compute_mass_balance = PETSC_TRUE

      case ('COMPUTE_MASS_BALANCE','MASS_BALANCE','COMPUTE_MASS_BALANCE_NEW','MASS_BALANCE_NEW')
        option%compute_mass_balance_new = PETSC_TRUE

!....................

      case ('TIMESTEPPER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            if (associated(flow_solver)) then
              call TimestepperRead(flow_stepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call TimestepperRead(tran_stepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case default
            if (associated(master_stepper)) then
              call TimestepperRead(master_stepper,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
        end select

!....................

      case ('LINEAR_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            if (associated(flow_solver)) then
              call SolverReadLinear(flow_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call SolverReadLinear(tran_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case default
            if (associated(master_solver)) then
              call SolverReadLinear(master_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
        end select

!....................

      case ('NEWTON_SOLVER')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('FLOW')
            if (associated(flow_solver)) then
              call SolverReadNewton(flow_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case('TRAN','TRANSPORT')
            if (associated(tran_solver)) then
              call SolverReadNewton(tran_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
          case default
            if (associated(master_solver)) then
              call SolverReadNewton(master_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
        end select

!....................

      case ('FLUID_PROPERTY','FLUID_PROPERTIES')

        realization%fluid_properties => FluidPropertyCreate(option%nphase)
        
        count = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'FLUID_PROPERTIES')
          
          if (InputCheckExit(input,option)) exit
         
          count = count + 1 
          if (count > option%nphase) exit              
                        
          call InputReadDouble(input,option,realization%fluid_properties%diff_base(count))
          call InputErrorMsg(input,option,'diff_base','FLUID_PROPERTIES')          
        
          call InputReadDouble(input,option,realization%fluid_properties%diff_exp(count))
          call InputErrorMsg(input,option,'diff_exp','FLUID_PROPERTIES')
! hardware the diffusion coefficient of SC and gas phases for now
          option%difsc=2.13D-5           
          option%difgs=2.13D-5
        enddo
        
!....................

      case ('THRM','THERMAL_PROPERTY','THERMAL_PROPERTIES')

        count = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'THRM')

          if (InputCheckExit(input,option)) exit
       
          count = count + 1
          thermal_property => ThermalPropertyCreate()
      
          call InputReadInt(input,option,thermal_property%id)
          call InputErrorMsg(input,option,'id','THRM')

          call InputReadDouble(input,option,thermal_property%rock_density)
          call InputErrorMsg(input,option,'rock density','THRM')

          call InputReadDouble(input,option,thermal_property%spec_heat)
          call InputErrorMsg(input,option,'cpr','THRM')
        
          call InputReadDouble(input,option,thermal_property%therm_cond_dry)
          call InputErrorMsg(input,option,'ckdry','THRM')
        
          call InputReadDouble(input,option,thermal_property%therm_cond_wet)
          call InputErrorMsg(input,option,'ckwet','THRM')
        
          call InputReadDouble(input,option,thermal_property%tort_bin_diff)
          call InputErrorMsg(input,option,'tau','THRM')

          call InputReadDouble(input,option,thermal_property%vap_air_diff_coef)
          call InputErrorMsg(input,option,'cdiff','THRM')

          call InputReadDouble(input,option,thermal_property%exp_binary_diff)
          call InputErrorMsg(input,option,'cexp','THRM')

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
        allocate(option%dencpr(count))
        allocate(option%ckwet(count))
        
        ! fill arrays with values from linked list
        thermal_property => realization%thermal_properties
        do
        
          if (.not.associated(thermal_property)) exit
          
          id = thermal_property%id
          
          if (id > count) then
            call printErrMsg(option,'Thermal property id greater than &
                                    &number of thermal properties')
          endif
                    
          option%dencpr(id) = thermal_property%rock_density * &
                              thermal_property%spec_heat
          option%ckwet(id) = thermal_property%therm_cond_wet
          
          thermal_property => thermal_property%next
          
        enddo
        
!....................

      case ('PCKR','SATURATION_FUNCTION','SATURATION_FUNCTIONS')
      
        count = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'PCKR')

          if (InputCheckExit(input,option)) exit
       
          count = count + 1
          saturation_function => SaturationFunctionCreate(option)
          
          call InputReadInt(input,option,saturation_function%id)
          call InputErrorMsg(input,option,'id','PCKR')
          
          call InputReadInt(input,option,saturation_function%saturation_function_itype)
          call InputErrorMsg(input,option,'icaptype','PCKR')
      
          select case(option%iflowmode)
            case(MPH_MODE,THC_MODE,RICHARDS_MODE)
              do np=1, option%nphase
                call InputReadDouble(input,option,saturation_function%Sr(np))
                call InputErrorMsg(input,option,'Sr','PCKR')
              enddo 
            case default
              call InputReadDouble(input,option,saturation_function%Sr(1))
              call InputErrorMsg(input,option,'Sr','PCKR')
          end select
        
          call InputReadDouble(input,option,saturation_function%m)
          call InputErrorMsg(input,option,'pckrm','PCKR')
          saturation_function%lambda = saturation_function%m

          call InputReadDouble(input,option,saturation_function%alpha)
          call InputErrorMsg(input,option,'alpha','PCKR')

          call InputReadDouble(input,option,saturation_function%pcwmax)
          call InputErrorMsg(input,option,'pcwmax','PCKR')
      
          call InputReadDouble(input,option,saturation_function%betac)
          call InputErrorMsg(input,option,'pbetac','PCKR')
      
          call InputReadDouble(input,option,saturation_function%power)
          call InputErrorMsg(input,option,'pwrprm','PCKR')

          call SaturationFunctionComputeSpline(option,saturation_function)
          
          call SaturationFunctionAddToList(saturation_function, &
                                           realization%saturation_functions)

        enddo
        
        ! allocate dynamic arrays holding saturation function information
        select case(option%iflowmode)
          case(MPH_MODE,THC_MODE,RICHARDS_MODE)
            allocate(option%sir(1:option%nphase,count))
          case default
        end select
  
        ! fill arrays with values from linked list
        saturation_function => realization%saturation_functions
        do 
        
          if (.not.associated(saturation_function)) exit
          
          id = saturation_function%id
          
          if (id > count) then
            call printErrMsg(option,'Saturation function id greater than &
                                    &number of saturation functions')
          endif
          
          select case(option%iflowmode)
            case(MPH_MODE,THC_MODE,RICHARDS_MODE)
              do i=1,option%nphase
                option%sir(i,id) = saturation_function%Sr(i)
              enddo
            case default
          end select
          
          saturation_function => saturation_function%next
          
        enddo
        
        call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                           realization%saturation_function_array)
        
!....................
      
      case ('PHIK','MATERIAL','MATERIALS')

        count = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'PHIK')

          if (InputCheckExit(input,option)) exit
       
          count = count + 1
          material => MaterialCreate()

          call InputReadWord(input,option,material%name,PETSC_TRUE)
          call InputErrorMsg(input,option,'name','PHIK')
                
          call InputReadInt(input,option,material%id)
          call InputErrorMsg(input,option,'id','PHIK')
                
          call InputReadInt(input,option,material%icap)
          call InputErrorMsg(input,option,'icap','PHIK')
  
          call InputReadInt(input,option,material%ithrm)
          call InputErrorMsg(input,option,'ithrm','PHIK')
  
          call InputReadDouble(input,option,material%porosity)
          call InputErrorMsg(input,option,'por','PHIK')
          
          call InputReadDouble(input,option,material%tortuosity)
          call InputErrorMsg(input,option,'tor','PHIK')
  
          call InputReadDouble(input,option,material%permeability(1,1))
          call InputErrorMsg(input,option,'permx','PHIK')
  
          call InputReadDouble(input,option,material%permeability(2,2))
          call InputErrorMsg(input,option,'permy','PHIK')
  
          call InputReadDouble(input,option,material%permeability(3,3))
          call InputErrorMsg(input,option,'permz','PHIK')
  
          call InputReadDouble(input,option,material%permeability_pwr)
          call InputErrorMsg(input,option,'permpwr','PHIK')
          
          material%permeability(1:3,1:3) = material%permeability(1:3,1:3)
          
          call MaterialAddToList(material,realization%materials)
          
        enddo          

        call MaterialConvertListToArray(realization%materials, &
                                        realization%material_array)
                                        
!....................

      case ('USE_TOUCH_OPTIONS')
        option%use_touch_options = PETSC_TRUE

      case ('MPI_IO')
!        call PetscOptionsInsertString('-viewer_binary_mpiio')

      case ('HANDSHAKE_IO')
        call InputReadInt(input,option,option%io_handshake_buffer_size)
        call InputErrorMsg(input,option,'io_handshake_buffer_size','HANDSHAKE_IO')

      case ('OVERWRITE_RESTART_TRANSPORT')
        option%overwrite_restart_transport = PETSC_TRUE

      case ('OVERWRITE_RESTART_FLOW_PARAMS')
        option%overwrite_restart_flow = PETSC_TRUE

      case ('TIME')

        call InputReadStringErrorMsg(input,option,'TIME')
      
        call InputReadWord(input,option,word,PETSC_FALSE)
      
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
          if (OptionPrint(option)) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') realization%output_option%tunit
          endif
          stop
        endif

        call InputReadWord(input,option,word,PETSC_FALSE)
        if (.not.InputError(input)) then
          call StringToUpper(word)
          if (StringCompare(word,'EVERY',FIVE_INTEGER)) then
            periodic_output_flag = PETSC_TRUE
            call InputReadDouble(input,option,periodic_rate)
          endif
        endif

        continuation_flag = PETSC_TRUE
        do
          if (.not.continuation_flag) exit
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          continuation_flag = PETSC_FALSE
          if (index(input%buf,backslash) > 0) continuation_flag = PETSC_TRUE
          input%ierr = 0
          do
            if (InputError(input)) exit
            call InputReadDouble(input,option,temp_real)
            if (.not.InputError(input)) then
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

        call InputReadStringErrorMsg(input,option,'DTST')

        call InputReadDouble(input,option,master_stepper%dt_min)
        call InputDefaultMsg(input,option,'dt_min')
            
        continuation_flag = PETSC_TRUE
        temp_int = 0       
        do
          if (.not.continuation_flag) exit
          call InputReadFlotranString(input,option)
          if (InputError(input)) exit
          continuation_flag = PETSC_FALSE
          if (index(input%buf,backslash) > 0) continuation_flag = PETSC_TRUE
          input%ierr = 0
          do
            if (InputError(input)) exit
            call InputReadDouble(input,option,temp_real)
            if (.not.InputError(input)) then
              waypoint => WaypointCreate()
              waypoint%time = temp_real
              call InputReadDouble(input,option,waypoint%dt_max)
              call InputErrorMsg(input,option,'dt_max','dtst')
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
        call BreakthroughRead(breakthrough,input,option)
        call RealizationAddBreakthrough(realization,breakthrough)        
      
!....................
      case('SDST')
        print *, 'SDST needs to be implemented'
        stop
#if 0
! Needs implementation         
        allocate(master_stepper%steady_eps(option%nflowdof))
        do j=1,option%nflowdof
          call InputReadDouble(input,option,master_stepper%steady_eps(j))
          call InputDefaultMsg(input,option,'steady tol')
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
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time','WALLCLOCK_STOP') 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*3600.d0
      
!....................
      case default
    
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select

  enddo

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
  use String_module

  implicit none 

  type(option_type) :: option
  
  call StringToUpper(option%flowmode)
  select case(option%flowmode)
    case('THC')
      option%iflowmode = THC_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%nflowdof = 3
      option%nflowspec = 2
    case('RICHARDS')
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%nflowdof = 1
      option%nflowspec = 1
    case('MPH','MPHASE')
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 1      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
    case default
      option%io_buffer = 'Mode: '//trim(option%flowmode)//' not recognized.'
      call printErrMsg(option)
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
  use Patch_module
  use Level_module

  implicit none
  
  type(realization_type) :: realization
  
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
        call RegionReadFromFile(region,realization%option, &
                                region%filename)
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
  use Patch_module
  use Logging_module
  use Input_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  character(len=MAXWORDLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  type(input_type), pointer :: input
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  input => realization%input

  if (index(filename,'.h5') > 0) then
    call HDF5ReadMaterialsFromFile(realization,filename)
  else
    call GridCreateNaturalToGhostedHash(grid,option)
    status = 0
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      option%io_buffer = 'File: ' // trim(filename) // ' not found.'
      call printErrMsg(option)
    endif
    call PetscLogEventBegin(logging%event_hash_map, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                            PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
    do
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call InputReadInt(input,option,material_id)
        call InputErrorMsg(input,option,'material id','STRATA')
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
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr, ierr2
  PetscInt :: count, read_count, i
  PetscInt, pointer :: indices(:)
  PetscReal, pointer :: values(:)
  PetscInt, parameter :: block_size = 10000
  Vec :: natural_vec, global_vec

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
      option%io_buffer = 'File: ' // trim(filename) // ' not found.'
      call printErrMsg(option)
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
      if (option%myrank == option%io_rank) &
        read(fid,*,iostat=ierr) values(1:read_count)
      call mpi_bcast(ierr,ONE_INTEGER,MPI_INTEGER,option%io_rank, &
                     option%comm,ierr2)      
      if (ierr /= 0) then
        option%io_buffer = 'Insufficent data in file: ' // filename
        call printErrMsg(option)
      endif
      if (option%myrank == option%io_rank) then
        call VecSetValues(natural_vec,read_count,indices,values,INSERT_VALUES, &
                          ierr)
      endif
      count = count + read_count
    enddo
    call mpi_bcast(count,ONE_INTEGER,MPI_INTEGER,option%io_rank, &
                   option%comm,ierr)      
    if (count /= grid%nmax) then
      write(option%io_buffer, &
            '("Number of data in file (",i8, &
            &") does not match size of vector (", &
            &i8,")")'), count, grid%nlmax
      call printErrMsg(option)
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
