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
#include "finclude/petscts.h"


  public :: Init, InitReadStochasticCardFromInput, InitReadInputFilenames

contains

! ************************************************************************** !
!
! Init: Initializes pflotran
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine Init(simulation)

  use Simulation_module
  use Option_module
  use Grid_module
  use Solver_module
  use Discretization_module
  use Realization_class
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
! use Mass_Balance_module
  use Logging_module  
  use Database_module
  use Database_hpt_module
  use Input_module
  use Condition_Control_module
  
  use Flash2_module
  use MPHASE_module
  use Immis_module
  use Miscible_module
  use Richards_module
  use Richards_MFD_module
  use TH_module
  use THC_module
  use THMC_module
  use General_module
  
  use Reactive_Transport_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  
  use Secondary_Continuum_module, only : SecondaryRTUpdateIterate
  
  use Global_module
  use Variables_module
  use Water_EOS_module
!  use Utility_module
  use Output_module
  use Output_Aux_module
  use Regression_module
    
#ifdef SURFACE_FLOW
  use Surface_Field_module
  use Surface_Flow_Module
  use Surface_Global_module
  use Surface_Init_module
  use Surface_Realization_class
  use Surface_TH_module
  use Unstructured_Grid_module
#endif

#ifdef PROCESS_MODEL
  use Synchronizer_module
  use Process_Model_Coupler_module
  use Process_Model_Richards_class
  use Process_Model_RT_class
  use Process_Model_TH_class
  use Process_Model_THC_class
  use Process_Model_Base_class
  use Process_Model_module
#endif

  implicit none
  
  type(simulation_type) :: simulation
  character(len=MAXSTRINGLENGTH) :: filename, filename_out

  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
!geh: for some reason, changing type(realization) to class(realization) causes
!     the code to crash with gfortran 4.7
!geh  class(realization_type), pointer :: realization
  type(realization_type), pointer :: realization
  class(realization_type), pointer :: realization_class_ptr
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  type(flow_debug_type), pointer :: debug
  type(waypoint_list_type), pointer :: waypoint_list
  type(input_type), pointer :: input
  type(output_variable_type), pointer :: output_variable
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: temp_int
  PetscInt :: flowortranpc    
  PetscErrorCode :: ierr
  PCSide:: pcside
  PetscReal :: r1, r2, r3, r4, r5, r6
  PetscReal :: min_value
  SNESLineSearch :: linesearch

#ifdef SURFACE_FLOW
  type(stepper_type), pointer               :: surf_flow_stepper
  type(solver_type), pointer                :: surf_flow_solver
  type(surface_field_type), pointer         :: surf_field
  type(surface_realization_type), pointer   :: surf_realization
#endif

#ifdef PROCESS_MODEL
  type(process_model_coupler_type), pointer :: cur_process_model_coupler
  type(process_model_coupler_type), pointer :: cur_process_model_coupler_top
  class(process_model_base_type), pointer :: cur_process_model
#endif

  ! popped in TimestepperInitializeRun()
  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr)
  call PetscLogEventBegin(logging%event_init,ierr)
  
  ! set pointers to objects
#ifdef PROCESS_MODEL
  ! set pointers to objects
  flow_stepper => TimestepperCreate()
  tran_stepper => TimestepperCreate()
  realization => RealizationCreate(simulation%option) 
#ifdef SURFACE_FLOW
  surf_realization  => SurfRealizCreate(option)
  surf_flow_stepper => TimestepperCreate()
  surf_field        => surf_realization%surf_field  
#endif  
#else
  flow_stepper => simulation%flow_stepper
  tran_stepper => simulation%tran_stepper
  realization => simulation%realization
#ifdef SURFACE_FLOW
  surf_realization  => simulation%surf_realization
  surf_flow_stepper => simulation%surf_flow_stepper
  surf_field        => surf_realization%surf_field  
#endif
#endif

  discretization => realization%discretization
  option => realization%option
  field => realization%field
  debug => realization%debug
  input => realization%input
  
  option%init_stage = PETSC_TRUE

  nullify(flow_solver)
  nullify(tran_solver)
  
  if (OptionPrintToScreen(option)) then
    temp_int = 6
    call InitPrintPFLOTRANHeader(option,temp_int)
  endif
  
  realization%input => InputCreate(IN_UNIT,option%input_filename,option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  if (OptionPrintToFile(option)) then
    call InitPrintPFLOTRANHeader(option,option%fid_out)
  endif
  
  ! read required cards
  call InitReadRequiredCardsFromInput(realization)
#ifdef SURFACE_FLOW
  surf_realization%input => InputCreate(IN_UNIT,option%input_filename,option)
  surf_realization%subsurf_filename = realization%discretization%filename
  call SurfaceInitReadRequiredCards(surf_realization)
#endif

  patch => realization%patch

  if (associated(patch)) then
     if (associated(patch%grid)) then
        grid => patch%grid
     endif
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
    option%use_isothermal = PETSC_TRUE  ! assume default isothermal when only transport
#ifndef PROCESS_MODEL  
    call TimestepperDestroy(simulation%flow_stepper)
#endif
    nullify(flow_stepper)
  endif
    
  ! initialize transport mode
  if (option%ntrandof > 0) then
    tran_solver => tran_stepper%solver
  else
#ifndef PROCESS_MODEL  
    call TimestepperDestroy(simulation%tran_stepper)
#endif
    nullify(tran_stepper)
  endif

#ifdef SURFACE_FLOW
  ! initialize surface-flow mode
  if (option%nsurfflowdof > 0) then
    surf_flow_solver => surf_flow_stepper%solver
    waypoint_list => WaypointListCreate()
    surf_realization%waypoints => waypoint_list
  else
#ifndef PROCESS_MODEL  
    call TimestepperDestroy(simulation%surf_flow_stepper)
#endif
    nullify(surf_flow_solver)
  endif
#endif

  ! initialize plot variables
  realization%output_option%output_variable_list => OutputVariableListCreate()
  realization%output_option%aveg_output_variable_list => OutputVariableListCreate()
#ifdef SURFACE_FLOW
  ! initialize plot variables
  surf_realization%output_option%output_variable_list => &
    OutputVariableListCreate()
  surf_realization%output_option%aveg_output_variable_list => &
    OutputVariableListCreate()
#endif

  ! read in the remainder of the input file
#ifdef SURFACE_FLOW
  call InitReadInput(realization,flow_stepper,tran_stepper, &
                     simulation%regression,surf_realization, &
                     surf_flow_stepper) 
#else
  call InitReadInput(realization,flow_stepper,tran_stepper, &
                     simulation%regression) 
#endif  
  call InputDestroy(realization%input)

  ! initialize reference density
  if (option%reference_water_density < 1.d-40) then
#ifndef DONT_USE_WATEOS
    call wateos(option%reference_temperature,option%reference_pressure, &
                option%reference_water_density,r1,r2,r3,r4,r5,r6, &
                option%scale,ierr)
#else
    call density(option%reference_temperature,option%reference_pressure, &
                 option%reference_water_density)
#endif                 
  endif
  
  ! read reaction database
  
  if (associated(realization%reaction)) then
    if (realization%reaction%use_full_geochemistry) then
        call DatabaseRead(realization%reaction,option)
        call BasisInit(realization%reaction,option)    
    else
      ! turn off activity coefficients since the database has not been read
      realization%reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(realization%reaction%primary_species_print(option%ntrandof))
      realization%reaction%primary_species_print = PETSC_TRUE
    endif
  endif

  ! Initialize flow databases (e.g. span wagner, etc.)
  select case(option%iflowmode)
    case(MPH_MODE, FLASH2_MODE, IMS_MODE)
      call init_span_wanger(realization)
  end select
  

  ! create grid and allocate vectors
  call RealizationCreateDiscretization(realization)
#ifdef SURFACE_FLOW
  if (option%nsurfflowdof>0) then
    call SurfRealizCreateDiscretization(surf_realization)
  endif
#endif  

  call RegressionCreateMapping(simulation%regression,realization)

  if (realization%discretization%itype == STRUCTURED_GRID .or. &
      realization%discretization%itype == STRUCTURED_GRID_MIMETIC) then
    if (OptionPrintToScreen(option)) then
      write(*,'(/," Requested processors and decomposition = ", &
               & i5,", npx,y,z= ",3i4)') &
          option%mycommsize,grid%structured_grid%npx, &
          grid%structured_grid%npy,grid%structured_grid%npz
      write(*,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
          grid%structured_grid%npx_final,grid%structured_grid%npy_final, &
          grid%structured_grid%npz_final
    endif
    if (OptionPrintToScreen(option)) then
      write(option%fid_out,'(/," Requested processors and decomposition = ", &
                           & i5,", npx,y,z= ",3i4)') &
          option%mycommsize,grid%structured_grid%npx,grid%structured_grid%npy, &
          grid%structured_grid%npz
      write(option%fid_out,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
          grid%structured_grid%npx_final,grid%structured_grid%npy_final, &
          grid%structured_grid%npz_final
    endif
  endif

  ! update flow mode based on optional input
  if (option%nflowdof > 0) then
  
    if (flow_solver%J_mat_type == MATAIJ) then
      select case(option%iflowmode)
        case(MPH_MODE,TH_MODE,THC_MODE,THMC_MODE,IMS_MODE, FLASH2_MODE, G_MODE, MIS_MODE)
          option%io_buffer = 'AIJ matrix not supported for current mode: '// &
                             option%flowmode
          call printErrMsg(option)
      end select
    endif

    if (OptionPrintToScreen(option)) then
      write(*,'(" number of dofs = ",i3,", number of phases = ",i3,i2)') &
        option%nflowdof,option%nphase
      select case(option%iflowmode)
        case(FLASH2_MODE)
          write(*,'(" mode = FLASH2: p, T, s/X")')
        case(MPH_MODE)
          write(*,'(" mode = MPH: p, T, s/X")')
        case(IMS_MODE)
          write(*,'(" mode = IMS: p, T, s")')
        case(MIS_MODE)
          write(*,'(" mode = MIS: p, Xs")')
        case(TH_MODE)
          write(*,'(" mode = THC: p, T")')
        case(THC_MODE)
          write(*,'(" mode = THC: p, T, s/X")')
        case(THMC_MODE)
          write(*,'(" mode = THMC: p, T, s/X")')
        case(RICHARDS_MODE)
          write(*,'(" mode = Richards: p")')  
        case(G_MODE)    
      end select
    endif

    call printMsg(option,"  Beginning setup of FLOW SNES ")

    call SolverCreateSNES(flow_solver,option%mycomm)  
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
      case(TH_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,THResidual, &
                             realization,ierr)
      case(THC_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,THCResidual, &
                             realization,ierr)
      case(THMC_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,THMCResidual, &
                             realization,ierr)
      case(RICHARDS_MODE)
        select case(realization%discretization%itype)
          case(STRUCTURED_GRID_MIMETIC)
            call SNESSetFunction(flow_solver%snes,field%flow_r_faces, &
                                 RichardsResidualMFDLP, &
                                 realization,ierr)
          case default
            call SNESSetFunction(flow_solver%snes,field%flow_r, &
                                 RichardsResidual, &
                                 realization,ierr)
        end select
      case(MPH_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,MphaseResidual, &
                             realization,ierr)
      case(IMS_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,ImmisResidual, &
                             realization,ierr)
      case(MIS_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,MiscibleResidual, &
                             realization,ierr)
      case(FLASH2_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,FLASH2Residual, &
                             realization,ierr)
      case(G_MODE)
        call SNESSetFunction(flow_solver%snes,field%flow_r,GeneralResidual, &
                             realization,ierr)
    end select
    
    if (flow_solver%J_mat_type == MATMFFD) then
      call MatCreateSNESMF(flow_solver%snes,flow_solver%J,ierr)
    endif

    select case(option%iflowmode)
      case(TH_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             THJacobian,realization,ierr)
      case(THC_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             THCJacobian,realization,ierr)
      case(THMC_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             THMCJacobian,realization,ierr)
      case(RICHARDS_MODE)
        select case(realization%discretization%itype)
          case(STRUCTURED_GRID_MIMETIC)
            call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             RichardsJacobianMFDLP,realization,ierr)
          case default !sp 
            call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             RichardsJacobian,realization,ierr)
        end select

      case(MPH_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             MPHASEJacobian,realization,ierr)
      case(IMS_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             ImmisJacobian,realization,ierr)
      case(MIS_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             MiscibleJacobian,realization,ierr)
      case(FLASH2_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             FLASH2Jacobian,realization,ierr)
      case(G_MODE)
        call SNESSetJacobian(flow_solver%snes,flow_solver%J,flow_solver%Jpre, &
                             GeneralJacobian,realization,ierr)
    end select
    
    ! by default turn off line search
    call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
    call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)
    ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
    if (option%verbosity >= 1) then
      string = '-flow_snes_view'
      call PetscOptionsInsertString(string, ierr)
    endif

    call SolverSetSNESOptions(flow_solver)

    ! If we are using a structured grid, set the corresponding flow DA 
    ! as the DA for the PCEXOTIC preconditioner, in case we choose to use it.
    ! The PCSetDA() call is ignored if the PCEXOTIC preconditioner is 
    ! no used.  We need to put this call after SolverCreateSNES() so that 
    ! KSPSetFromOptions() will already have been called.
    ! I also note that this preconditioner is intended only for the flow 
    ! solver.  --RTM
    if ((realization%discretization%itype == STRUCTURED_GRID_MIMETIC).or.&
                (realization%discretization%itype == STRUCTURED_GRID)) then
      call PCSetDM(flow_solver%pc, &
                   realization%discretization%dm_nflowdof,ierr);
    endif

    option%io_buffer = 'Solver: ' // trim(flow_solver%ksp_type)
    call printMsg(option)
    option%io_buffer = 'Preconditioner: ' // trim(flow_solver%pc_type)
    call printMsg(option)

    ! shell for custom convergence test.  The default SNES convergence test  
    ! is call within this function. 
    flow_stepper%convergence_context => &
      ConvergenceContextCreate(flow_solver,option,grid)
    call SNESSetConvergenceTest(flow_solver%snes,ConvergenceTest, &
                                flow_stepper%convergence_context, &
                                PETSC_NULL_FUNCTION,ierr) 

    
 
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
            dabs(option%saturation_change_limit) > 0.d0) then
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPreCheck(linesearch, &
                                         RichardsCheckUpdatePre, &
                                         realization,ierr)
        endif
      case(G_MODE)
        call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
        call SNESLineSearchSetPreCheck(linesearch, &
                                       GeneralCheckUpdatePre, &
                                       realization,ierr)
      case(TH_MODE)
        if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
            dabs(option%pressure_change_limit) > 0.d0 .or. &
            dabs(option%temperature_change_limit) > 0.d0) then
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPreCheck(linesearch, &
                                         THCheckUpdatePre, &
                                         realization,ierr)
        endif
      case(THC_MODE)
        if (dabs(option%pressure_dampening_factor) > 0.d0 .or. &
            dabs(option%pressure_change_limit) > 0.d0 .or. &
            dabs(option%temperature_change_limit) > 0.d0) then
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPreCheck(linesearch, &
                                         THCCheckUpdatePre, &
                                         realization,ierr)
        endif
    end select
    
    
    if (option%check_stomp_norm) then
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPostCheck(linesearch, &
                                          RichardsCheckUpdatePost, &
                                          realization,ierr)
        case(G_MODE)
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPostCheck(linesearch, &
                                          GeneralCheckUpdatePost, &
                                          realization,ierr)
        case(TH_MODE)
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPostCheck(linesearch, &
                                          THCheckUpdatePost, &
                                          realization,ierr)
        case(THC_MODE)
          call SNESGetSNESLineSearch(flow_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPostCheck(linesearch, &
                                          THCCheckUpdatePost, &
                                          realization,ierr)
      end select
    endif
        
    
    call printMsg(option,"  Finished setting up FLOW SNES ")

#ifdef SURFACE_FLOW
    if(option%nsurfflowdof>0) then


      if(option%surf_flow_explicit) then

        ! Setup PETSc TS for explicit surface flow solution
        call printMsg(option,"  Beginning setup of SURF FLOW TS ")

        call SolverCreateTS(surf_flow_solver,option%mycomm)
        call TSSetProblemType(surf_flow_solver%ts,TS_NONLINEAR,ierr)
        select case(option%iflowmode)
          case (RICHARDS_MODE)
            call TSSetRHSFunction(surf_flow_solver%ts,PETSC_NULL_OBJECT, &
                                  SurfaceFlowRHSFunction, &
                                  surf_realization,ierr)
          case (TH_MODE)
            call TSSetRHSFunction(surf_flow_solver%ts,PETSC_NULL_OBJECT, &
                                  SurfaceTHRHSFunction, &
                                  surf_realization,ierr)
        end select
        call TSSetDuration(surf_flow_solver%ts,ONE_INTEGER, &
                           surf_realization%waypoints%last%time,ierr)

      else

        ! Setup PETSc SNES for implicit surface flow solution
        call printMsg(option,"  Beginning setup of SURF FLOW SNES ")

        call SolverCreateSNES(surf_flow_solver,option%mycomm)
        call SNESSetOptionsPrefix(surf_flow_solver%snes, "surf_flow_",ierr)
        call SolverCheckCommandLine(surf_flow_solver)

        if (surf_flow_solver%Jpre_mat_type == '') then
          if (surf_flow_solver%J_mat_type /= MATMFFD) then
            surf_flow_solver%Jpre_mat_type = surf_flow_solver%J_mat_type
          else
            surf_flow_solver%Jpre_mat_type = MATBAIJ
          endif
        endif

        call DiscretizationCreateJacobian( &
                                  surf_realization%discretization, &
                                  NFLOWDOF, &
                                  surf_flow_solver%Jpre_mat_type, &
                                  surf_flow_solver%Jpre, &
                                  option)

        call MatSetOption(surf_flow_solver%Jpre,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr)
        call MatSetOption(surf_flow_solver%Jpre,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)

        call MatSetOptionsPrefix(surf_flow_solver%Jpre,"surf_flow_",ierr)

        if (surf_flow_solver%J_mat_type /= MATMFFD) then
          surf_flow_solver%J = surf_flow_solver%Jpre
        endif

        call SNESSetFunction(surf_flow_solver%snes,surf_field%flow_r, &
                              SurfaceFlowResidual, &
                              surf_realization,ierr)

        call SNESSetJacobian(surf_flow_solver%snes,surf_flow_solver%J, &
                            surf_flow_solver%Jpre, &
                            SurfaceFlowJacobian,surf_realization,ierr)
        ! by default turn off line search
        call SNESGetSNESLineSearch(surf_flow_solver%snes, linesearch, ierr)
        call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)

        ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
        if (option%verbosity >= 1) then
          string = '-surf_flow_snes_view'
          call PetscOptionsInsertString(string, ierr)
        endif

        call SolverSetSNESOptions(surf_flow_solver)

        option%io_buffer = 'Solver: ' // trim(surf_flow_solver%ksp_type)
        call printMsg(option)
        option%io_buffer = 'Preconditioner: ' // trim(surf_flow_solver%pc_type)
        call printMsg(option)

        ! shell for custom convergence test.  The default SNES convergence test
        ! is call within this function.
        surf_flow_stepper%convergence_context => &
          ConvergenceContextCreate(surf_flow_solver,option,grid)
        call SNESSetConvergenceTest(surf_flow_solver%snes,ConvergenceTest, &
                                    surf_flow_stepper%convergence_context, &
                                    PETSC_NULL_FUNCTION,ierr)
      endif ! if(option%surface_flow_explicit)
    endif ! if(option%nsurfflowdof>0)
#endif

  endif

  
  ! update transport mode based on optional input
  if (option%ntrandof > 0) then

    call printMsg(option,"  Beginning setup of TRAN SNES ")
    
    call SolverCreateSNES(tran_solver,option%mycomm)  
    call SNESSetOptionsPrefix(tran_solver%snes, "tran_",ierr)
    call SolverCheckCommandLine(tran_solver)
      
     if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
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
    else
      tran_solver%J_mat_type = MATAIJ
      tran_solver%Jpre_mat_type = MATAIJ

      call DiscretizationCreateJacobian(discretization,ONEDOF, &
                                        tran_solver%Jpre_mat_type, &
                                        tran_solver%Jpre,option)
    endif

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

    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then

      call SNESSetFunction(tran_solver%snes,field%tran_r,RTResidual,&
                           realization,ierr)

      if (tran_solver%J_mat_type == MATMFFD) then
        call MatCreateSNESMF(tran_solver%snes,tran_solver%J,ierr)
      endif
      
      call SNESSetJacobian(tran_solver%snes,tran_solver%J,tran_solver%Jpre, &
                           RTJacobian,realization,ierr)

      ! this could be changed in the future if there is a way to ensure that the linesearch
      ! update does not perturb concentrations negative.
      call SNESGetSNESLineSearch(tran_solver%snes, linesearch, ierr)
      call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)
      
      if (option%use_mc) then
        call SNESLineSearchSetPostCheck(linesearch, &
                                        SecondaryRTUpdateIterate, &
                                        realization,ierr)      
      endif
      
      ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
      if (option%verbosity >= 1) then
        string = '-tran_snes_view'
        call PetscOptionsInsertString(string, ierr)
      endif

    endif

    ! ensure setting of SNES options since they set KSP and PC options too
    call SolverSetSNESOptions(tran_solver)

    option%io_buffer = 'Solver: ' // trim(tran_solver%ksp_type)
    call printMsg(option)
    option%io_buffer = 'Preconditioner: ' // trim(tran_solver%pc_type)
    call printMsg(option)

    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then

      ! shell for custom convergence test.  The default SNES convergence test  
      ! is call within this function. 
      tran_stepper%convergence_context => &
        ConvergenceContextCreate(tran_solver,option,grid)
      call SNESSetConvergenceTest(tran_solver%snes,ConvergenceTest, &
                                  tran_stepper%convergence_context, &
                                  PETSC_NULL_FUNCTION,ierr) 

      ! this update check must be in place, otherwise reactive transport is likely
      ! to fail
      if (associated(realization%reaction)) then
        if (realization%reaction%check_update) then
          call SNESGetSNESLineSearch(tran_solver%snes, linesearch, ierr)
          call SNESLineSearchSetPreCheck(linesearch,RTCheckUpdate, &
                                         realization,ierr)
        endif
      endif
    endif
    
    call printMsg(option,"  Finished setting up TRAN SNES ")
  
  endif

  if (OptionPrintToScreen(option)) write(*,'("++++++++++++++++++++++++++++++++&
                     &++++++++++++++++++++++++++++",/)')

  call PetscLogEventBegin(logging%event_setup,ierr)
  
  ! read any regions provided in external files
  call readRegionFiles(realization)
  ! clip regions and set up boundary connectivity, distance  
  call RealizationLocalizeRegions(realization)
  call RealizatonPassPtrsToPatches(realization)
  ! link conditions with regions through couplers and generate connectivity
  call RealProcessMatPropAndSatFunc(realization)
  call RealizationProcessCouplers(realization)
  call RealizationProcessConditions(realization)
  call RealProcessFluidProperties(realization)
  call assignMaterialPropToRegions(realization)
  if(realization%discretization%lsm_flux_method) &
    call GridComputeMinv(realization%discretization%grid, &
                         realization%discretization%stencil_width,option)


  call RealizationInitAllCouplerAuxVars(realization)
  if (option%ntrandof > 0) then
    call printMsg(option,"  Setting up TRAN Realization ")
    call RealizationInitConstraints(realization)
    call printMsg(option,"  Finished setting up TRAN Realization ")  
  endif
  call RealizationPrintCouplers(realization)
  ! should we still support this
  if (option%use_generalized_grid) then 
    call printMsg(option,'Reading structured grid from hdf5')
    if (.not.associated(patch%imat)) &
      allocate(patch%imat(grid%ngmax))  ! allocate material id array
    call ReadStructuredGridHDF5(realization)
  endif
  call PetscLogEventEnd(logging%event_setup,ierr)
  if (.not.option%steady_state) then
    ! add waypoints associated with boundary conditions, source/sinks etc. to list
    call RealizationAddWaypointsToList(realization)
    ! fill in holes in waypoint data
    call WaypointListFillIn(option,realization%waypoints)
    call WaypointListRemoveExtraWaypnts(option,realization%waypoints)
  ! geh- no longer needed
  !  ! convert times from input time to seconds
  !  call WaypointConvertTimes(realization%waypoints,realization%output_option%tconv)
  endif
  
  if (associated(flow_stepper)) then
    flow_stepper%cur_waypoint => realization%waypoints%first
  endif
  if (associated(tran_stepper)) then
    tran_stepper%cur_waypoint => realization%waypoints%first
  endif
  
  ! initialize global auxiliary variable object
  call GlobalSetup(realization)
  ! initialize FLOW
  ! set up auxillary variable arrays
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(TH_MODE)
        call THSetup(realization)
      case(THC_MODE)
        call THCSetup(realization)
      case(THMC_MODE)
        call THMCSetup(realization)
      case(RICHARDS_MODE)
        call RichardsSetup(realization)
      case(MPH_MODE)
        call MphaseSetup(realization)
      case(IMS_MODE)
        call ImmisSetup(realization)
      case(MIS_MODE)
        call MiscibleSetup(realization)
      case(FLASH2_MODE)
        call Flash2Setup(realization)
      case(G_MODE)
        call GeneralSetup(realization)
    end select
  
    ! assign initial conditionsRealizAssignFlowInitCond
    call CondControlAssignFlowInitCond(realization)

    ! override initial conditions if they are to be read from a file
    if (len_trim(option%initialize_flow_filename) > 1) then
      call readFlowInitialCondition(realization, &
                                    option%initialize_flow_filename)
    endif
  
    select case(option%iflowmode)
      case(TH_MODE)
        call THUpdateAuxVars(realization)
      case(THC_MODE)
        call THCUpdateAuxVars(realization)
      case(THMC_MODE)
        call THMCUpdateAuxVars(realization)
      case(RICHARDS_MODE)
#ifdef DASVYAT
       if (option%mimetic) then
!        call RichardsInitialPressureReconstruction(realization)
!        write(*,*) "RichardsInitialPressureReconstruction"
!        read(*,*)
       end if
#endif 
        call RichardsUpdateAuxVars(realization)
      case(MPH_MODE)
        call MphaseUpdateAuxVars(realization)
      case(IMS_MODE)
        call ImmisUpdateAuxVars(realization)
      case(MIS_MODE)
        call MiscibleUpdateAuxVars(realization)
      case(FLASH2_MODE)
        call Flash2UpdateAuxVars(realization)
      case(G_MODE)
        call GeneralUpdateAuxVars(realization,PETSC_TRUE)
    end select
  else ! no flow mode specified
    if (len_trim(realization%nonuniform_velocity_filename) > 0) then
      call InitReadVelocityField(realization)
    endif
  endif

  if (option%ntrandof > 0) then
    call RTSetup(realization)

    ! initialize densities and saturations
    if (option%nflowdof == 0) then
      call GlobalSetAuxVarScalar(realization,option%reference_pressure, &
                                 LIQUID_PRESSURE)
      call GlobalSetAuxVarScalar(realization,option%reference_temperature, &
                                 TEMPERATURE)
      call GlobalSetAuxVarScalar(realization,option%reference_saturation, &
                                 LIQUID_SATURATION)
      call GlobalSetAuxVarScalar(realization,option%reference_water_density, &
                                 LIQUID_DENSITY)
    endif

    ! initial concentrations must be assigned after densities are set !!!
    call CondControlAssignTranInitCond(realization)
    ! override initial conditions if they are to be read from a file
    if (len_trim(option%initialize_transport_filename) > 1) then
      call readTransportInitialCondition(realization, &
                                         option%initialize_transport_filename)
    endif
    ! PETSC_FALSE = no activity coefficients
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)
    ! at this point the auxvars have been computed with activity coef = 1.d0
    ! to use intitial condition with activity coefs /= 1.d0, must update
    ! activity coefs and recompute auxvars
    if (realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
      !geh: you may ask, why call this twice....  We need to iterate at least
      !     once to ensure that the activity coefficients are more accurate.
      !     Otherwise, the total component concentrations can be quite
      !     different from what is defined in the input file.
      call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
    endif
  endif
  
  ! Add plot variables that are not mode specific
  if (realization%output_option%print_porosity) then
    ! add porosity to header
    call OutputVariableAddToList( &
           realization%output_option%output_variable_list, &
           'Porosity',OUTPUT_GENERIC,'-',POROSITY)  
  endif
  if (realization%output_option%print_permeability) then
    ! add permeability to header
    call OutputVariableAddToList( &
           realization%output_option%output_variable_list, &
           'Permeability X',OUTPUT_GENERIC,'m^2',PERMEABILITY)
  endif
  if (realization%output_option%print_iproc) then
    output_variable => OutputVariableCreate('Processor ID',OUTPUT_DISCRETE,'', &
                                            PROCESSOR_ID)
    output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
    output_variable%iformat = 1 ! integer
    call OutputVariableAddToList( &
           realization%output_option%output_variable_list,output_variable)
  endif

  ! write material ids
  output_variable => OutputVariableCreate('Material ID',OUTPUT_DISCRETE,'', &
                                          MATERIAL_ID)
  output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList( &
         realization%output_option%output_variable_list,output_variable)  

  
  ! print info
  if (associated(flow_stepper)) then
    string = 'Flow Stepper:'
    call TimestepperPrintInfo(flow_stepper,option%fid_out,string,option)
  endif    
  if (associated(tran_stepper)) then
    string = 'Transport Stepper:'
    call TimestepperPrintInfo(tran_stepper,option%fid_out,string,option)
  endif    
#ifdef SURFACE_FLOW
   if (option%nsurfflowdof>0) then
    string = 'Surface Flow Stepper:'
    call TimestepperPrintInfo(surf_flow_stepper,option%fid_out,string,option)
  endif
#endif
  if (associated(flow_solver)) then
    string = 'Flow Newton Solver:'
    call SolverPrintNewtonInfo(flow_solver,OptionPrintToScreen(option), &
                               OptionPrintToFile(option),option%fid_out, &
                               string)
  endif    
  if (associated(tran_solver)) then
    string = 'Transport Newton Solver:'
    call SolverPrintNewtonInfo(tran_solver,OptionPrintToScreen(option), &
                               OptionPrintToFile(option),option%fid_out, &
                               string)
  endif    
#ifdef SURFACE_FLOW
  if (associated(surf_flow_solver).and.(.not.option%surf_flow_explicit)) then
    string = 'Surface Flow Newton Solver:'
    call SolverPrintNewtonInfo(surf_flow_solver,OptionPrintToScreen(option), &
                               OptionPrintToFile(option),option%fid_out, &
                               string)
  endif
#endif
  if (associated(flow_solver)) then
    string = 'Flow Linear Solver:'
    call SolverPrintLinearInfo(flow_solver,string,option)
  endif    
  if (associated(tran_solver)) then
    string = 'Transport Linear Solver'
    call SolverPrintLinearInfo(tran_solver,string,option)
  endif    
#ifdef SURFACE_FLOW
  if (associated(surf_flow_solver).and.(.not.option%surf_flow_explicit)) then
    string = 'Surface Flow Linear Solver:'
    call SolverPrintLinearInfo(surf_flow_solver,string,option)
  endif
  if (associated(surf_flow_solver).and.option%surf_flow_explicit) then
    string = 'Surface Flow TS Solver:'
    if (OptionPrintToScreen(option)) then
      write(*,*),' '
      write(*,*),string
    endif
    call TSView(surf_flow_solver%ts,PETSC_VIEWER_STDOUT_WORLD,ierr)
  endif
#endif

  if (debug%print_couplers) then
    call verifyAllCouplers(realization)
  endif
  if (debug%print_waypoints) then
    call WaypointListPrint(realization%waypoints,option,realization%output_option)
  endif

#ifdef OS_STATISTICS
  call RealizationPrintGridStatistics(realization)
#endif
  
  ! check for non-initialized data sets, e.g. porosity, permeability
  call RealizationNonInitializedData(realization)

#if defined(PETSC_HAVE_HDF5)
#if !defined(HDF5_BROADCAST)
  call printMsg(option,"Default HDF5 method is used in Initialization")
#else
  call printMsg(option,"Glenn's HDF5 broadcast method is used in Initialization")
#endif
#endif
!PETSC_HAVE_HDF5

#ifdef SURFACE_FLOW
  if(option%nsurfflowdof > 0) then
    ! Check if surface-flow is compatible with the given flowmode
    select case(option%iflowmode)
      case(RICHARDS_MODE,TH_MODE)
      case default
        option%io_buffer = 'For surface-flow only RICHARDS and TH mode implemented'
        call printErrMsgByRank(option)
    end select

    call SurfaceInitReadRegionFiles(surf_realization)
    call SurfRealizMapSurfSubsurfGrids(realization,surf_realization)
    call SurfRealizLocalizeRegions(surf_realization)
    call SurfRealizPassFieldPtrToPatches(surf_realization)
    call SurfRealizProcessMatProp(surf_realization)
    call SurfRealizProcessCouplers(surf_realization)
    call SurfRealizProcessConditions(surf_realization)
    !call RealProcessFluidProperties(surf_realization)
    call SurfaceInitMatPropToRegions(surf_realization)
    call SurfRealizInitAllCouplerAuxVars(surf_realization)
    !call SurfaceRealizationPrintCouplers(surf_realization)

    ! add waypoints associated with boundary conditions, source/sinks etc. to list
    call SurfRealizAddWaypointsToList(surf_realization)
    call WaypointListFillIn(option,surf_realization%waypoints)
    call WaypointListRemoveExtraWaypnts(option,surf_realization%waypoints)
    if (associated(flow_stepper)) then
      simulation%surf_flow_stepper%cur_waypoint => surf_realization%waypoints%first
    endif

    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call SurfaceFlowSetup(surf_realization)
      case default
      case(TH_MODE)
        call SurfaceTHSetup(surf_realization)
    end select

    call SurfaceGlobalSetup(surf_realization)
    ! initialize FLOW
    ! set up auxillary variable arrays

    ! assign initial conditionsRealizAssignFlowInitCond
    call CondControlAssignFlowInitCondSurface(surf_realization)

    ! override initial conditions if they are to be read from a file
    if (len_trim(option%surf_initialize_flow_filename) > 1) then
      option%io_buffer = 'For surface-flow initial conditions cannot be read from file'
      call printErrMsgByRank(option)
    endif
  
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call SurfaceFlowUpdateAuxVars(surf_realization)
        if (surf_realization%option%subsurf_surf_coupling == SEQ_COUPLED) then
          call SurfaceFlowCreateSurfSubsurfVec( &
                          simulation%realization, surf_realization)
        endif
      case(TH_MODE)
        call SurfaceTHUpdateAuxVars(surf_realization)
        if (surf_realization%option%subsurf_surf_coupling == SEQ_COUPLED) then
          call SurfaceTHCreateSurfSubsurfVec( &
                          simulation%realization, surf_realization)
        endif
      case default
        option%io_buffer = 'For surface-flow only RICHARDS and TH mode implemented'
        call printErrMsgByRank(option)
    end select
  endif ! option%nsurfflowdof > 0
#endif

  call printMsg(option," ")
  call printMsg(option,"  Finished Initialization")
  call PetscLogEventEnd(logging%event_init,ierr)

#ifdef PROCESS_MODEL
  !----------------------------------------------------------------------------!
  ! This section for setting up new process model approach
  !----------------------------------------------------------------------------!
  simulation%output_option => realization%output_option
  simulation%synchronizer => SynchronizerCreate()
  simulation%synchronizer%option => realization%option
  simulation%synchronizer%output_option => realization%output_option
!  simulation%synchronizer%waypoints => WaypointListCopy(realization%waypoints)
  simulation%synchronizer%waypoints => realization%waypoints
  nullify(cur_process_model)
  cur_process_model_coupler => ProcessModelCouplerCreate()
  cur_process_model_coupler%option => option
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        cur_process_model => PMRichardsCreate()
      case(TH_MODE)
        cur_process_model => PMTHCreate()
      case(THC_MODE)
        cur_process_model => PMTHCCreate()
    end select
    cur_process_model%option => realization%option
    cur_process_model%output_option => realization%output_option
    cur_process_model_coupler%process_model_list => cur_process_model
    cur_process_model_coupler%pm_ptr%ptr => cur_process_model
    cur_process_model_coupler%depth = 0
  endif
  if (option%ntrandof > 0) then
    cur_process_model => PMRTCreate()
    cur_process_model%output_option => realization%output_option
    cur_process_model%option => realization%option
    if (associated(cur_process_model_coupler%process_model_list)) then
      cur_process_model_coupler%below => ProcessModelCouplerCreate()
      cur_process_model_coupler%below%option => option
      cur_process_model_coupler%below%process_model_list => cur_process_model
      cur_process_model_coupler%below%pm_ptr%ptr => cur_process_model
      cur_process_model_coupler%below%depth = 1
    else
      cur_process_model_coupler%process_model_list => cur_process_model
      cur_process_model_coupler%pm_ptr%ptr => cur_process_model
      cur_process_model_coupler%depth = 0
    endif
  endif  
  simulation%process_model_coupler_list => cur_process_model_coupler
  simulation%synchronizer%process_model_coupler_list => &
    cur_process_model_coupler
  
  cur_process_model_coupler_top => simulation%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    cur_process_model_coupler => cur_process_model_coupler_top
    do
      if (.not.associated(cur_process_model_coupler)) exit
      cur_process_model => cur_process_model_coupler%process_model_list
      do
        if (.not.associated(cur_process_model)) exit
        realization_class_ptr => realization
        select type(cur_process_model)
          class is (process_model_richards_type)
            call cur_process_model%PMRichardsSetRealization( &
                                                         realization_class_ptr)
            call cur_process_model_coupler%SetTimestepper(flow_stepper)
            flow_stepper%dt = option%flow_dt
          class is (process_model_rt_type)
            call cur_process_model%PMRTSetRealization(realization_class_ptr)
            call cur_process_model_coupler%SetTimestepper(tran_stepper)
            tran_stepper%dt = option%tran_dt
          class is (process_model_th_type)
            call cur_process_model%PMTHSetRealization(realization_class_ptr)
            call cur_process_model_coupler%SetTimestepper(flow_stepper)
            flow_stepper%dt = option%flow_dt
          class is (process_model_thc_type)
            call cur_process_model%PMTHCSetRealization(realization_class_ptr)
            call cur_process_model_coupler%SetTimestepper(flow_stepper)
            flow_stepper%dt = option%flow_dt
        end select
        call cur_process_model%Init()
        call SNESSetFunction( &
                           cur_process_model_coupler%timestepper%solver%snes, &
                           cur_process_model%residual_vec, &
                           PMResidual, &
                           cur_process_model_coupler%pm_ptr,ierr)
        call SNESSetJacobian( &
                           cur_process_model_coupler%timestepper%solver%snes, &
                           cur_process_model_coupler%timestepper%solver%J, &
                           cur_process_model_coupler%timestepper%solver%Jpre, &
                           PMJacobian, &
                           cur_process_model_coupler%pm_ptr,ierr)
        cur_process_model => cur_process_model%next
      enddo
      cur_process_model_coupler => cur_process_model_coupler%below
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%next
  enddo
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
#endif  

end subroutine Init

! ************************************************************************** !
!
! InitReadStochasticCardFromInput: Reads stochastic card from input file
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine InitReadStochasticCardFromInput(stochastic,option)

  use Option_module
  use Input_module
  use Stochastic_Aux_module

  implicit none
  
  type(stochastic_type), pointer :: stochastic
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  type(input_type), pointer :: input
  PetscBool :: print_warning
  
  input => InputCreate(IN_UNIT,option%input_filename,option)

  ! MODE information
  string = "STOCHASTIC"
  print_warning = PETSC_FALSE
  call InputFindStringInFile(input,option,string,print_warning)

  if (.not.InputError(input)) then
    if (.not.associated(stochastic)) then
      stochastic => StochasticCreate()
    endif
    call StochasticRead(stochastic,input,option)
  endif
  
  call InputDestroy(input)

end subroutine InitReadStochasticCardFromInput

! ************************************************************************** !
!
! InitReadInputFilenames: Reads filenames for multi-simulation runs
! author: Glenn Hammond
! date: 08/11/09
!
! ************************************************************************** !
subroutine InitReadInputFilenames(option,filenames)

  use Option_module
  use Input_module

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: filename_count
  type(input_type), pointer :: input
  PetscBool :: card_found

  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = "FILENAMES"
  call InputFindStringInFile(input,option,string) 

  card_found = PETSC_FALSE
  if (InputError(input)) then
    ! if the FILENAMES card is not included, we will assume that only
    ! filenames exist in the file.
    rewind(input%fid)
  else
    card_found = PETSC_TRUE
  endif
    
  filename_count = 0     
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
  enddo
  
  allocate(filenames(filename_count))
  filenames = ''
  rewind(input%fid) 

  if (card_found) then
    string = "FILENAMES"
    call InputFindStringInFile(input,option,string) 
  endif
  
  filename_count = 0     
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
    filenames(filename_count) = filename
  enddo

  call InputDestroy(input)

end subroutine InitReadInputFilenames

! ************************************************************************** !
!
! InitReadRequiredCardsFromInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine InitReadRequiredCardsFromInput(realization)

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_module
  use String_module
  use Patch_module
  use Level_module
  use Realization_class

  use Reaction_module  
  use Reaction_Aux_module  

  implicit none

  type(realization_type) :: realization

  character(len=MAXSTRINGLENGTH) :: string
  
  type(patch_type), pointer :: patch, patch2 
  type(level_type), pointer :: level
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
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
#if defined(SCORPIO)
  string = "HDF5_WRITE_GROUP_SIZE"
  call InputFindStringInFile(input,option,string)
  if (.not.InputError(input)) then  
    call InputReadInt(input,option,option%hdf5_write_group_size)
    call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
    call InputSkipToEnd(input,option,'HDF5_WRITE_GROUP_SIZE')
  endif

  string = "HDF5_READ_GROUP_SIZE"
  call InputFindStringInFile(input,option,string)
  if (.not.InputError(input)) then  
    call InputReadInt(input,option,option%hdf5_read_group_size)
    call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')
  endif
 rewind(input%fid)

  call Create_IOGroups(option)

#endif

!.........................................................................

  ! GRID information
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call DiscretizationReadRequiredCards(discretization,input,option)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID,STRUCTURED_GRID_MIMETIC)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(realization%level_list)) then
        realization%level_list => LevelCreateList()
      endif
      level => LevelCreate()
      call LevelAddToList(level,realization%level_list)
      call PatchAddToList(patch,level%patch_list)
      realization%patch => patch
  end select
!.........................................................................

  if ((realization%discretization%itype == STRUCTURED_GRID).or. &
        (realization%discretization%itype == STRUCTURED_GRID_MIMETIC)) then  ! look for processor decomposition
    
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
 
      if (option%myrank == option%io_rank .and. option%print_to_screen) then
        option%io_buffer = ' Processor Decomposition:'
        call printMsg(option)
        write(option%io_buffer,'("  npx   = ",3x,i4)') grid%structured_grid%npx
        call printMsg(option)
        write(option%io_buffer,'("  npy   = ",3x,i4)') grid%structured_grid%npy
        call printMsg(option)
        write(option%io_buffer,'("  npz   = ",3x,i4)') grid%structured_grid%npz
        call printMsg(option)
      endif
  
      if (option%mycommsize /= grid%structured_grid%npx * &
                             grid%structured_grid%npy * &
                             grid%structured_grid%npz) then
        write(option%io_buffer,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%npx*grid%structured_grid%npy* &
                       grid%structured_grid%npz,' commsize = ',option%mycommsize
        call printErrMsg(option)
      endif
    endif
  endif
  
!.........................................................................

  ! CHEMISTRY information
  string = "CHEMISTRY"
  call InputFindStringInFile(input,option,string)

  if (.not.InputError(input)) then
    call ReactionInit(realization%reaction,input,option)
  endif
    
end subroutine InitReadRequiredCardsFromInput

! ************************************************************************** !
!
! InitReadInput: Reads pflow input file
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
#ifdef SURFACE_FLOW
subroutine InitReadInput(subsurface_realization,subsurface_flow_stepper, &
                         subsurface_tran_stepper,subsurface_regression, &
                         surface_realization,surface_flow_stepper)
#else
subroutine InitReadInput(subsurface_realization,subsurface_flow_stepper, &
                         subsurface_tran_stepper,subsurface_regression)
#endif

  use Simulation_module
  use Option_module
  use Field_module
  use Grid_module
  use Unstructured_Grid_Aux_module
  use Structured_Grid_module
  use Solver_module
  use Material_module
  use Saturation_Function_module  
  use Dataset_Aux_module
  use Fluid_module
  use Realization_class
  use Timestepper_module
  use Region_module
  use Condition_module
  use Constraint_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Reaction_module
  use Reaction_Aux_module
  use Discretization_module
  use Input_module
  use String_module
  use Units_module
  use Uniform_Velocity_module
  use Mineral_module
  use Regression_module
  use Output_Aux_module
  use Output_Tecplot_module
  
#ifdef SURFACE_FLOW
  use Surface_Flow_module
  use Surface_Init_module
#endif
#ifdef SOLID_SOLUTION
  use Solid_Solution_module, only : SolidSolutionReadFromInputFile
#endif
 
  implicit none
  
  type(realization_type), pointer :: subsurface_realization
  type(stepper_type), pointer :: subsurface_flow_stepper
  type(stepper_type), pointer :: subsurface_tran_stepper
  type(regression_type), pointer :: subsurface_regression
  
#ifdef SURFACE_FLOW
  type(surface_realization_type), pointer :: surf_realization
  type(stepper_type), pointer :: surf_flow_stepper
#endif

  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
    
  PetscBool :: continuation_flag
  
  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal :: units_conversion
  PetscInt :: temp_int
  PetscInt :: count, id
  
  PetscBool :: velocities
  PetscBool :: flux_velocities
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate
  
  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_type), pointer :: sec_tran_constraint
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  
  type(waypoint_type), pointer :: waypoint
  
  type(material_property_type), pointer :: material_property
  type(fluid_property_type), pointer :: fluid_property
  type(saturation_function_type), pointer :: saturation_function

  type(realization_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch   
  type(solver_type), pointer :: flow_solver
  type(solver_type), pointer :: tran_solver
  type(solver_type), pointer :: default_solver
  type(stepper_type), pointer :: flow_stepper
  type(stepper_type), pointer :: tran_stepper
  type(stepper_type), pointer :: default_stepper
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(uniform_velocity_dataset_type), pointer :: uniform_velocity_dataset
  type(dataset_type), pointer :: dataset
  type(input_type), pointer :: input

  nullify(flow_stepper)
  nullify(tran_stepper)
  nullify(flow_solver)
  nullify(tran_solver)
  
  realization => subsurface_realization
  patch => realization%patch

  if (associated(patch)) grid => patch%grid

  option => realization%option
  output_option => realization%output_option
  field => realization%field
  reaction => realization%reaction
  input => realization%input

  tran_stepper => subsurface_tran_stepper
  if (associated(tran_stepper)) then
    tran_solver => tran_stepper%solver
    tran_solver%itype = TRANSPORT_CLASS
  endif
  flow_stepper => subsurface_flow_stepper
  if (associated(flow_stepper)) then
    flow_solver => flow_stepper%solver
    flow_solver%itype = FLOW_CLASS
  endif

  if (associated(flow_stepper)) then
    default_stepper => flow_stepper
    default_solver => flow_solver
  else
    default_stepper => tran_stepper
    default_solver => tran_solver
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
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('CHEMISTRY')
        call ReactionReadPass2(reaction,input,option)

!....................
      case ('NONUNIFORM_VELOCITY')
        call InputReadNChars(input,option, &
                             realization%nonuniform_velocity_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','NONUNIFORM_VELOCITY') 

      case ('UNIFORM_VELOCITY')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        uniform_velocity_dataset%rank = 3
        uniform_velocity_dataset%interpolation_method = 1 ! 1 = STEP
        uniform_velocity_dataset%is_cyclic = PETSC_FALSE
        allocate(uniform_velocity_dataset%times(1))
        uniform_velocity_dataset%times = 0.d0
        allocate(uniform_velocity_dataset%values(3,1))
        uniform_velocity_dataset%values = 0.d0
        call InputReadDouble(input,option,uniform_velocity_dataset%values(1,1))
        call InputErrorMsg(input,option,'velx','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(2,1))
        call InputErrorMsg(input,option,'vely','UNIFORM_VELOCITY')
        call InputReadDouble(input,option,uniform_velocity_dataset%values(3,1))
        call InputErrorMsg(input,option,'velz','UNIFORM_VELOCITY')
        ! read units, if present
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          units_conversion = UnitsConvertToInternal(word,option) 
          uniform_velocity_dataset%values(:,1) = &
            uniform_velocity_dataset%values(:,1) * units_conversion
        endif
        call UniformVelocityDatasetVerify(option,uniform_velocity_dataset)
        realization%uniform_velocity_dataset => uniform_velocity_dataset
      
      case ('VELOCITY_DATASET')
        uniform_velocity_dataset => UniformVelocityDatasetCreate()
        call UniformVelocityDatasetRead(uniform_velocity_dataset,input,option)
        realization%uniform_velocity_dataset => uniform_velocity_dataset

!....................
      case ('DEBUG')
        call DebugRead(realization%debug,input,option)
               
!....................
      case ('PRINT_PRIMAL_GRID')
        option%print_explicit_primal_grid = PETSC_TRUE
        
!....................
      case ('PRINT_DUAL_GRID')
        option%print_explicit_dual_grid = PETSC_TRUE

!....................
      case ('MAX_CHANGE')
        call InputReadDouble(input,option,option%dpmxe)
        call InputErrorMsg(input,option,'dpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dtmpmxe)
        call InputErrorMsg(input,option,'dtmpmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dsmxe)
        call InputErrorMsg(input,option,'dsmxe','MAX_CHANGE')
        call InputReadDouble(input,option,option%dcmxe)
        call InputErrorMsg(input,option,'dcmxe','MAX_CHANGE')
        
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
        if (option%iflowmode == G_MODE) then
          call FlowConditionGeneralRead(flow_condition,input,option)
        else
          call FlowConditionRead(flow_condition,input,option)
        endif
        call FlowConditionAddToList(flow_condition,realization%flow_conditions)
        nullify(flow_condition)
        
!....................
      case ('TRANSPORT_CONDITION')
        if (.not.associated(reaction)) then
          option%io_buffer = 'TRANSPORT_CONDITIONs not supported without ' // &
            'CHEMISTRY.'
          call printErrMsg(option)
        endif
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
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
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
        dataset => DatasetCreate()
        call InputReadWord(input,option,dataset%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Dataset name') 
        call DatasetRead(dataset,input,option)
        call DatasetAddToList(dataset,realization%datasets)
        nullify(dataset)
        
!....................

      case('REFERENCE_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_pressure)
        call InputDefaultMsg(input,option,'Reference Pressure') 

!....................

      case('REFERENCE_DENSITY')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_water_density)
        call InputDefaultMsg(input,option,'Reference Density') 

!....................

      case('MINIMUM_HYDROSTATIC_PRESSURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%minimum_hydrostatic_pressure)
        call InputDefaultMsg(input,option,'Minimum Hydrostatic Pressure') 

!......................

      case('REFERENCE_TEMPERATURE')
        call InputReadStringErrorMsg(input,option,card)
        call InputReadDouble(input,option,option%reference_temperature)
        call InputDefaultMsg(input,option,'Reference Temperature') 

!......................

      case('ANI_RELATIVE_PERMEABILTY')
        option%ani_relative_permeability = PETSC_TRUE

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

      case('NONISOTHERMAL')
        option%use_isothermal = PETSC_FALSE

!......................

      case('ISOTHERMAL')
        option%use_isothermal = PETSC_TRUE
        
!......................

      case('MULTIPLE_CONTINUUM')
        option%use_mc = PETSC_TRUE
        
!......................

      case('UPDATE_FLOW_PERMEABILITY')
        option%update_flow_perm = PETSC_TRUE
        
!......................

      case('SECONDARY_CONTINUUM_SOLVER')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can only be used ' // &
                             'with MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif      
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('KEARST')
            option%secondary_continuum_solver = 1
          case('HINDMARSH')
            option%secondary_continuum_solver = 2
          case('THOMAS')
            option%secondary_continuum_solver = 3
          case default
            option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                               'HINDMARSH or KEARST. For single component'// &
                               'chemistry THOMAS can be used.'
          call printErrMsg(option)    
        end select        
!....................

      case('SECONDARY_CONSTRAINT')
        if (.not.option%use_mc) then
          option%io_buffer = 'SECONDARY_CONSTRAINT can only be used with ' // &
                             'MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif
        if (.not.associated(reaction)) then
          option%io_buffer = 'SECONDARY_CONSTRAINT not supported without' // &
                             'CHEMISTRY.'
          call printErrMsg(option)
        endif
        sec_tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input,option,sec_tran_constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'secondary constraint','name') 
        call printMsg(option,sec_tran_constraint%name)
        call TranConstraintRead(sec_tran_constraint,reaction,input,option)
        realization%sec_transport_constraint => sec_tran_constraint
        nullify(sec_tran_constraint)        

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
         if (OptionPrintToScreen(option)) print *, option%m_nacl
!......................

      case ('RESTART')
        option%restart_flag = PETSC_TRUE
        call InputReadNChars(input,option,option%restart_filename,MAXSTRINGLENGTH, &
                             PETSC_TRUE)
        call InputErrorMsg(input,option,'RESTART','Restart file name') 
        call InputReadDouble(input,option,option%restart_time)
        if (input%ierr == 0) then
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr == 0) then
            option%restart_time = option%restart_time* &
                                  UnitsConvertToInternal(word,option)
          else
            call InputDefaultMsg(input,option,'RESTART, time units')
          endif
        endif

!......................

      case ('CHECKPOINT')
        option%checkpoint_flag = PETSC_TRUE
        call InputReadInt(input,option,option%checkpoint_frequency)
        call InputErrorMsg(input,option,'CHECKPOINT','Checkpoint frequency') 

!......................

      case ('NUMERICAL_JACOBIAN_FLOW')
        option%numerical_derivatives_flow = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_RXN')
        option%numerical_derivatives_rxn = PETSC_TRUE

!......................

      case ('NUMERICAL_JACOBIAN_MULTI_COUPLE')
        option%numerical_derivatives_multi_coupling = PETSC_TRUE

!......................

      case ('COMPUTE_STATISTICS')
        option%compute_statistics = PETSC_TRUE

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
            if (associated(default_stepper)) then
              call TimestepperRead(default_stepper,input,option)
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
            if (associated(default_solver)) then
              call SolverReadLinear(default_solver,input,option)
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
            if (associated(default_solver)) then
              call SolverReadNewton(default_solver,input,option)
            else
              call InputSkipToEnd(input,option,card)
            endif
        end select

!....................

      case ('FLUID_PROPERTY')

        fluid_property => FluidPropertyCreate()
        call FluidPropertyRead(fluid_property,input,option)
        call FluidPropertyAddToList(fluid_property,realization%fluid_properties)
        nullify(fluid_property)
        
!....................

      case ('SATURATION_FUNCTION')
      
        saturation_function => SaturationFunctionCreate(option)
        call InputReadWord(input,option,saturation_function%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','SATURATION_FUNCTION')
        call SaturationFunctionRead(saturation_function,input,option)
        call SaturationFunctionComputeSpline(option,saturation_function)
        call PermFunctionComputeSpline(option,saturation_function)
        call SaturationFunctionAddToList(saturation_function, &
                                         realization%saturation_functions)
        nullify(saturation_function)   

!....................
      
      case ('MATERIAL_PROPERTY')

        material_property => MaterialPropertyCreate()
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')        
        call MaterialPropertyRead(material_property,input,option)
        call MaterialPropertyAddToList(material_property,realization%material_properties)
        nullify(material_property)

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

      case ('INITIALIZE_FLOW_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_flow_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_FLOW_FROM_FILE') 

      case ('INITIALIZE_TRANSPORT_FROM_FILE')
        call InputReadNChars(input,option,option%initialize_transport_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','INITIALIZE_TRANSPORT_FROM_FILE') 

      case ('CENTRAL_DIFFERENCE')
        option%use_upwinding = PETSC_FALSE

!....................
      case ('OBSERVATION')
        observation => ObservationCreate()
        call ObservationRead(observation,input,option)
        call RealizationAddObservation(realization,observation)        
      
!.....................
      case ('WALLCLOCK_STOP')
        option%wallclock_stop_flag = PETSC_TRUE
        call InputReadDouble(input,option,option%wallclock_stop_time)
        call InputErrorMsg(input,option,'stop time','WALLCLOCK_STOP') 

        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) word = 'h'
        call InputDefaultMsg(input,option,'WALLCLOCK_STOP time units')
        units_conversion = UnitsConvertToInternal(word,option) 
        ! convert from hrs to seconds and add to start_time
        option%wallclock_stop_time = option%start_time + &
                                     option%wallclock_stop_time*units_conversion
      
!....................
      case ('OUTPUT')
        velocities = PETSC_FALSE
        flux_velocities = PETSC_FALSE
        mass_flowrate = PETSC_FALSE
        energy_flowrate = PETSC_FALSE
        aveg_mass_flowrate = PETSC_FALSE
        aveg_energy_flowrate = PETSC_FALSE
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','OUTPUT') 
          call StringToUpper(word)
          select case(trim(word))
            case('NO_FINAL','NO_PRINT_FINAL')
              output_option%print_final = PETSC_FALSE
            case('NO_INITIAL','NO_PRINT_INITIAL')
              output_option%print_initial = PETSC_FALSE
            case('PROCESSOR_ID')
              output_option%print_iproc = PETSC_TRUE
            case('PERMEABILITY')
              output_option%print_permeability = PETSC_TRUE
            case('POROSITY')
              output_option%print_porosity = PETSC_TRUE
            case('MASS_BALANCE')
              option%compute_mass_balance_new = PETSC_TRUE
            case('PRINT_COLUMN_IDS')
              output_option%print_column_ids = PETSC_TRUE
            case('TIMES')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'units','OUTPUT')
              units_conversion = UnitsConvertToInternal(word,option) 
              continuation_flag = PETSC_TRUE
              do
                continuation_flag = PETSC_FALSE
                if (index(input%buf,backslash) > 0) &
                  continuation_flag = PETSC_TRUE
                input%ierr = 0
                do
                  if (InputError(input)) exit
                  call InputReadDouble(input,option,temp_real)
                  if (.not.InputError(input)) then
                    waypoint => WaypointCreate()
                    waypoint%time = temp_real*units_conversion
                    waypoint%print_output = PETSC_TRUE    
                    call WaypointInsertInList(waypoint,realization%waypoints)
                  endif
                enddo
                if (.not.continuation_flag) exit
                call InputReadFlotranString(input,option)
                if (InputError(input)) exit
              enddo
            case('OUTPUT_FILE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,OUTPUT_FILE')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_file = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%output_file_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,OUTPUT_FILE')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,OUTPUT_FILE.'
                  call printErrMsg(option)
              end select
            case('SCREEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment','OUTPUT,SCREEN')
              call StringToUpper(word)
              select case(trim(word))
                case('OFF')
                  option%print_to_screen = PETSC_FALSE
                case('PERIODIC')
                  call InputReadInt(input,option,output_option%screen_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,SCREEN')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,SCREEN.'
                  call printErrMsg(option)
              end select
            case('PERIODIC')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                                 'OUTPUT,PERIODIC')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC,TIME')
                  units_conversion = UnitsConvertToInternal(word,option) 
                  output_option%periodic_output_time_incr = temp_real* &
                                                            units_conversion
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (input%ierr == 0) then
                    if (StringCompareIgnoreCase(word,'between')) then

                      call InputReadDouble(input,option,temp_real)
                      call InputErrorMsg(input,option,'start time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'start time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      units_conversion = UnitsConvertToInternal(word,option) 
                      temp_real = temp_real * units_conversion
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      if (.not.StringCompareIgnoreCase(word,'and')) then
                        input%ierr = 1
                      endif
                      call InputErrorMsg(input,option,'and', &
                                          'OUTPUT,PERIODIC,TIME"')
                      call InputReadDouble(input,option,temp_real2)
                      call InputErrorMsg(input,option,'end time', &
                                         'OUTPUT,PERIODIC,TIME')
                      call InputReadWord(input,option,word,PETSC_TRUE)
                      call InputErrorMsg(input,option,'end time units', &
                                         'OUTPUT,PERIODIC,TIME')
                      temp_real2 = temp_real2 * units_conversion
                      do
                        waypoint => WaypointCreate()
                        waypoint%time = temp_real
                        waypoint%print_output = PETSC_TRUE    
                        call WaypointInsertInList(waypoint,realization%waypoints)
                        temp_real = temp_real + output_option%periodic_output_time_incr
                        if (temp_real > temp_real2) exit
                      enddo
                      output_option%periodic_output_time_incr = 0.d0
                    else
                      input%ierr = 1
                      call InputErrorMsg(input,option,'between', &
                                          'OUTPUT,PERIODIC,TIME')
                    endif
                  endif                  
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,PERIODIC,'// &
                                     'TIMESTEP.'
                  call printErrMsg(option)
              end select
            case('PERIODIC_OBSERVATION')
              output_option%print_observation = PETSC_TRUE
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'time increment', &
                'OUTPUT, PERIODIC_OBSERVATION')
              call StringToUpper(word)
              select case(trim(word))
                case('TIME')
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'time increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'time increment units', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIME')
                  units_conversion = UnitsConvertToInternal(word,option) 
                  output_option%periodic_tr_output_time_incr = temp_real* &
                                                               units_conversion
                case('TIMESTEP')
                  call InputReadInt(input,option, &
                                    output_option%periodic_tr_output_ts_imod)
                  call InputErrorMsg(input,option,'timestep increment', &
                                     'OUTPUT,PERIODIC_OBSERVATION,TIMESTEP')
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,'// &
                                     'PERIODIC_OBSERVATION,TIMESTEP.'
                  call printErrMsg(option)
              end select
            case('FORMAT')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','OUTPUT,FORMAT') 
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputDefaultMsg(input,option, &
                                       'OUTPUT,FORMAT,HDF5,# FILES')
                  if (len_trim(word) > 0) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                        output_option%times_per_h5_file = 1
                        call InputReadWord(input,option,word,PETSC_TRUE)
                        if (len_trim(word)>0) then
                          select case(trim(word))
                            case('TIMES_PER_FILE')
                              call InputReadInt(input,option, &
                                              output_option%times_per_h5_file)
                              call InputErrorMsg(input,option,'timestep increment', &
                                        'OUTPUT,FORMAT,MULTIPLE_FILES,TIMES_PER_FILE')
                            case default
                              option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,'// &
                                     'FORMAT,MULTIPLE_FILES,TIMES_PER_FILE.'
                              call printErrMsg(option)
                          end select
                        endif
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recongnized.  Use "SINGLE_FILE" or ' // &
                          '"MULTIPLE_FILES".'
                        call printErrMsg(option)
                    end select
                  endif
                case ('MAD')
                  output_option%print_mad = PETSC_TRUE
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'TECPLOT','OUTPUT,FORMAT') 
                  call StringToUpper(word)
                  select case(trim(word))
                    case('POINT')
                      output_option%tecplot_format = TECPLOT_POINT_FORMAT
                    case('BLOCK')
                      output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                    case('FEBRICK')
                      output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                    case default
                      option%io_buffer = 'TECPLOT format (' // trim(word) // &
                                         ') not recongnized.'
                      call printErrMsg(option)
                  end select
                  if (output_option%tecplot_format == TECPLOT_POINT_FORMAT &
                      .and. option%mycommsize > 1) then
                    output_option%tecplot_format = TECPLOT_BLOCK_FORMAT
                  endif
                  if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
                    output_option%tecplot_format = TECPLOT_FEBRICK_FORMAT
                  endif
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  option%io_buffer = 'Keyword: ' // trim(word) // &
                                     ' not recognized in OUTPUT,FORMAT.'
                  call printErrMsg(option)
              end select
            case('VELOCITIES')
              velocities = PETSC_TRUE
            case('FLUXES_VELOCITIES')
              flux_velocities = PETSC_TRUE
            case('FLOWRATES','FLOWRATE')
              mass_flowrate = PETSC_TRUE
              energy_flowrate = PETSC_TRUE
            case('MASS_FLOWRATE')
              mass_flowrate = PETSC_TRUE
            case('ENERGY_FLOWRATE')
              energy_flowrate = PETSC_TRUE
            case('AVERAGE_FLOWRATES','AVERAGE_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
              aveg_energy_flowrate = PETSC_TRUE
            case('AVERAGE_MASS_FLOWRATE')
              aveg_mass_flowrate = PETSC_TRUE
            case('AVERAGE_ENERGY_FLOWRATE')
              aveg_energy_flowrate = PETSC_TRUE
            case ('HDF5_WRITE_GROUP_SIZE')
              call InputReadInt(input,option,option%hdf5_write_group_size)
              call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')
            case('VARIABLES')
              call OutputVariableRead(input,option,output_option%output_variable_list)
            case('AVERAGE_VARIABLES')
              call OutputVariableRead(input,option,output_option%aveg_output_variable_list)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in OUTPUT.'
              call printErrMsg(option)              
          end select
        enddo
        if (velocities) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_velocities = PETSC_TRUE
          if (output_option%print_hdf5) &
            output_option%print_hdf5_velocities = PETSC_TRUE
          if (output_option%print_vtk) &
            output_option%print_vtk_velocities = PETSC_TRUE
        endif
        if (flux_velocities) then
          if (output_option%print_tecplot) &
            output_option%print_tecplot_flux_velocities = PETSC_TRUE
          if (output_option%print_hdf5) &
           output_option%print_hdf5_flux_velocities = PETSC_TRUE
        endif
        if(output_option%aveg_output_variable_list%nvars>0) then
          if(output_option%periodic_output_time_incr==0.d0) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES defined without' // &
                               ' PERIODIC TIME being set.'
            call printErrMsg(option)
          endif
          if(.not.output_option%print_hdf5) then
            option%io_buffer = 'Keyword: AVERAGE_VARIABLES only defined for FORMAT HDF5'
            call printErrMsg(option)
          endif
        endif
        if (mass_flowrate.or.energy_flowrate.or.aveg_mass_flowrate.or.aveg_energy_flowrate) then
          if (output_option%print_hdf5) then
#ifndef STORE_FLOWRATES
            option%io_buffer='To output FLOWRATES/MASS_FLOWRATE/ENERGY_FLOWRATE, '// &
              'compile with -DSTORE_FLOWRATES'
            call printErrMsg(option)
#endif
            output_option%print_hdf5_mass_flowrate = mass_flowrate
            output_option%print_hdf5_energy_flowrate = energy_flowrate
            output_option%print_hdf5_aveg_mass_flowrate = aveg_mass_flowrate
            output_option%print_hdf5_aveg_energy_flowrate = aveg_energy_flowrate
            if(aveg_mass_flowrate.or.aveg_energy_flowrate) then
              if(output_option%periodic_output_time_incr==0.d0) then
                option%io_buffer = 'Keyword: AVEGRAGE_FLOWRATES/ ' // &
                  'AVEGRAGE_MASS_FLOWRATE/ENERGY_FLOWRATE defined without' // &
                  ' PERIODIC TIME being set.'
                call printErrMsg(option)
              endif
            endif
           option%store_flowrate = PETSC_TRUE
          else
            option%io_buffer='Output FLOWRATES/MASS_FLOWRATE/ENERGY_FLOWRATE ' // &
              'only available in HDF5 format'
            call printErrMsg(option)
          endif
        endif

!.....................
      case ('REGRESSION')
        call RegressionRead(subsurface_regression,input,option)

!.....................
      case ('TIME')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word','TIME') 
          select case(trim(word))
            case('STEADY_STATE')
              option%steady_state = PETSC_TRUE
            case('FINAL_TIME')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units','TIME')
              realization%output_option%tunit = trim(word)
              realization%output_option%tconv = UnitsConvertToInternal(word,option)
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*realization%output_option%tconv
              waypoint%print_output = PETSC_TRUE              
              call WaypointInsertInList(waypoint,realization%waypoints)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Initial Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Initial Timestep Size Time Units','TIME')
              default_stepper%dt_min = temp_real*UnitsConvertToInternal(word,option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Maximum Timestep Size','TIME') 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Maximum Timestep Size Time Units','TIME')
              waypoint => WaypointCreate()
              waypoint%dt_max = temp_real*UnitsConvertToInternal(word,option)
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,temp_real)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time','TIME') 
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  call InputErrorMsg(input,option,'Maximum Timestep Size Update Time Units','TIME')
                  waypoint%time = temp_real*UnitsConvertToInternal(word,option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" after ' // &
                                     'maximum timestep size should be "at".'
                  call printErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif     
              call WaypointInsertInList(waypoint,realization%waypoints)
            case default
              option%io_buffer = 'Keyword: ' // trim(word) // &
                                 ' not recognized in TIME.'
              call printErrMsg(option)              
          end select
        enddo

        if (associated(flow_stepper)) then
          flow_stepper%dt_min = default_stepper%dt_min
        endif
        if (associated(tran_stepper)) then
          tran_stepper%dt_min = default_stepper%dt_min
        endif
        option%flow_dt = default_stepper%dt_min
        option%tran_dt = default_stepper%dt_min
      
#ifdef SURFACE_FLOW
!.....................
      case ('SURFACE_FLOW')
        call SurfaceInitReadInput(surf_realization, &
                              simulation%surf_flow_stepper%solver,input,option)
        simulation%surf_flow_stepper%dt_min = surf_realization%dt_min
        simulation%surf_flow_stepper%dt_max = surf_realization%dt_max
        option%surf_subsurf_coupling_flow_dt = surf_realization%dt_coupling
        option%surf_flow_dt=simulation%surf_flow_stepper%dt_min

        ! Add first waypoint
        waypoint => WaypointCreate()
        waypoint%time = 0.d0
        call WaypointInsertInList(waypoint,surf_realization%waypoints)

        ! Add final_time waypoint to surface_realization
        waypoint => WaypointCreate()
        waypoint%final = PETSC_TRUE
        waypoint%time = realization%waypoints%last%time
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,surf_realization%waypoints)
#endif

!......................
      case ('HDF5_READ_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_read_group_size)
        call InputErrorMsg(input,option,'HDF5_READ_GROUP_SIZE','Group size')

!......................
      case ('HDF5_WRITE_GROUP_SIZE')
        call InputReadInt(input,option,option%hdf5_write_group_size)
        call InputErrorMsg(input,option,'HDF5_WRITE_GROUP_SIZE','Group size')

!....................
      case default
    
        option%io_buffer = 'Keyword ' // trim(word) // ' in input file ' // &
                           'not recognized'
        call printErrMsg(option)

    end select

  enddo
                                      
end subroutine InitReadInput

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
    case('TH')
      option%iflowmode = TH_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2
      option%nflowdof = 2
      option%nflowspec = 1
      option%use_isothermal = PETSC_FALSE
    case('THC')
      option%iflowmode = THC_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%use_isothermal = PETSC_FALSE
   case('THMC')
      option%iflowmode = THMC_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 6
      option%nflowspec = 2
      option%nmechdof = 3
      option%use_isothermal = PETSC_FALSE
    case('MIS','MISCIBLE')
      option%iflowmode = MIS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 2
      option%nflowspec = 2
    case('RICHARDS')
      option%iflowmode = RICHARDS_MODE
      option%nphase = 1
      option%liquid_phase = 1      
      option%nflowdof = 1
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE
    case('MPH','MPHASE')
      option%iflowmode = MPH_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%use_isothermal = PETSC_FALSE
    case('FLA2','FLASH2')
      option%iflowmode = FLASH2_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
      option%use_isothermal = PETSC_FALSE
   case('IMS','IMMIS','THS')
      option%iflowmode = IMS_MODE
      option%nphase = 2
      option%liquid_phase = 1      
      option%gas_phase = 2      
      option%nflowdof = 3
      option%nflowspec = 2
      option%itable = 2
    case('GENERAL')
      option%iflowmode = G_MODE
      option%nphase = 2
      option%liquid_phase = 1  ! liquid_pressure
      option%gas_phase = 2     ! gas_pressure

      option%air_pressure_id = 3
      option%capillary_pressure_id = 4
      option%vapor_pressure_id = 5

      option%water_id = 1
      option%air_id = 2
      option%energy_id = 3

      option%nflowdof = 3
      option%nflowspec = 2
      option%use_isothermal = PETSC_FALSE
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

  use Realization_class
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Level_module
  
  use HDF5_module

  implicit none
  
  type(realization_type) :: realization
  
  PetscReal, pointer :: icap_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: por0_p(:)
  PetscReal, pointer :: tor0_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xz_p(:)
  PetscReal, pointer :: perm_xy_p(:)
  PetscReal, pointer :: perm_yz_p(:)
  PetscReal, pointer :: perm_pow_p(:)
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: patch  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  type(material_property_type), pointer :: material_property, null_material_property
  type(region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  
  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field

  ! loop over all patches and allocation material id arrays
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      if (.not.associated(cur_patch%imat)) then
        allocate(cur_patch%imat(cur_patch%grid%ngmax))
        ! initialize to "unset"
        cur_patch%imat = -999
        ! also allocate saturation function id
        allocate(cur_patch%sat_func_id(cur_patch%grid%ngmax))
        cur_patch%sat_func_id = -999
      endif
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      grid => cur_patch%grid
      strata => cur_patch%strata%first
      do
        if (.not.associated(strata)) exit
        ! Read in cell by cell material ids if they exist
        if (.not.associated(strata%region) .and. strata%active) then
          call readMaterialsFromFile(realization,strata%realization_dependent, &
                                     strata%material_property_filename)
        ! Otherwise, set based on region
        else
          update_ghosted_material_ids = PETSC_TRUE
          region => strata%region
          material_property => strata%material_property
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
            if (strata%active) then
              cur_patch%imat(ghosted_id) = material_property%id
            else
              ! if not active, set material id to zero
              cur_patch%imat(ghosted_id) = 0
            endif
          enddo
        endif
        strata => strata%next
      enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
    
  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call RealLocalToLocalWithArray(realization,MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_material_property => MaterialPropertyCreate()
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
        if (option%mimetic) then
          call GridVecGetArrayF90(grid,field%perm0_xz,perm_xz_p,ierr)
          call GridVecGetArrayF90(grid,field%perm0_xy,perm_xy_p,ierr)
          call GridVecGetArrayF90(grid,field%perm0_yz,perm_yz_p,ierr)
        endif
        call GridVecGetArrayF90(grid,field%perm_pow,perm_pow_p,ierr)
      endif
      call GridVecGetArrayF90(grid,field%porosity0,por0_p,ierr)
      call GridVecGetArrayF90(grid,field%tortuosity0,tor0_p,ierr)
        
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        material_id = cur_patch%imat(ghosted_id)
        if (material_id == 0) then ! accommodate inactive cells
          material_property => null_material_property
        else if (material_id > 0 .and. &
                 material_id <= &
                 size(realization%material_property_array)) then
          material_property => &
            realization%material_property_array(material_id)%ptr
          if (.not.associated(material_property)) then
            write(dataset_name,*) material_id
            option%io_buffer = 'No material property for material id ' // &
                               trim(adjustl(dataset_name)) &
                               //  ' defined in input file.'
            call printErrMsgByRank(option)
          endif
        else if (material_id < -998) then 
          write(dataset_name,*) grid%nG2A(ghosted_id)
          option%io_buffer = 'Uninitialized material id in patch at cell ' // &
                             trim(adjustl(dataset_name))
          call printErrMsgByRank(option)
        else if (material_id > size(realization%material_property_array)) then
          write(option%io_buffer,*) material_id
          option%io_buffer = 'Unmatched material id in patch:' // &
            adjustl(trim(option%io_buffer))
          call printErrMsgByRank(option)
        else
          option%io_buffer = 'Something messed up with material ids. ' // &
            ' Possibly material ids not assigned to all grid cells. ' // &
            ' Contact Glenn!'
          call printErrMsgByRank(option)
        endif
        if (option%nflowdof > 0) then
          patch%sat_func_id(ghosted_id) = material_property%saturation_function_id
          icap_loc_p(ghosted_id) = material_property%saturation_function_id
          ithrm_loc_p(ghosted_id) = material_property%id
          perm_xx_p(local_id) = material_property%permeability(1,1)
          perm_yy_p(local_id) = material_property%permeability(2,2)
          perm_zz_p(local_id) = material_property%permeability(3,3)
          if (option%mimetic) then
            perm_xz_p(local_id) = material_property%permeability(1,3)
            perm_xy_p(local_id) = material_property%permeability(1,2)
            perm_yz_p(local_id) = material_property%permeability(2,3)
          endif
!          perm_pow_p(local_id) = ???
        endif
        por0_p(local_id) = material_property%porosity
        tor0_p(local_id) = material_property%tortuosity
      enddo

      if (option%nflowdof > 0) then
        call GridVecRestoreArrayF90(grid,field%icap_loc,icap_loc_p,ierr)
        call GridVecRestoreArrayF90(grid,field%ithrm_loc,ithrm_loc_p,ierr)
        call GridVecRestoreArrayF90(grid,field%perm0_xx,perm_xx_p,ierr)
        call GridVecRestoreArrayF90(grid,field%perm0_yy,perm_yy_p,ierr)
        call GridVecRestoreArrayF90(grid,field%perm0_zz,perm_zz_p,ierr)
        if (option%mimetic) then
          call GridVecRestoreArrayF90(grid,field%perm0_xz,perm_xz_p,ierr)
          call GridVecRestoreArrayF90(grid,field%perm0_xy,perm_xy_p,ierr)
          call GridVecRestoreArrayF90(grid,field%perm0_yz,perm_yz_p,ierr)
        endif
        call GridVecRestoreArrayF90(grid,field%perm_pow,perm_pow_p,ierr)
      endif
      call GridVecRestoreArrayF90(grid,field%porosity0,por0_p,ierr)
      call GridVecRestoreArrayF90(grid,field%tortuosity0,tor0_p,ierr)
        
      ! read in any user-defined property fields
      do material_id = 1, size(realization%material_property_array)
        material_property => &
               realization%material_property_array(material_id)%ptr
        if (associated(material_property)) then
          if (associated(material_property%permeability_dataset)) then
            call readPermeabilitiesFromFile(realization,material_property)
          endif
          if (associated(material_property%porosity_dataset)) then
            group_name = ''
            dataset_name = 'Porosity'
            call HDF5ReadCellIndexedRealArray(realization,field%work, &
                       material_property%porosity_dataset%filename, &
                       group_name, &
                       dataset_name, &
                       material_property%porosity_dataset%realization_dependent)
            call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
            call GridVecGetArrayF90(grid,field%porosity0,por0_p,ierr)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == &
                  material_property%id) then
                por0_p(local_id) = vec_p(local_id)
              endif
            enddo
            call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
            call GridVecRestoreArrayF90(grid,field%porosity0,por0_p,ierr)
          endif
        endif
      enddo
      
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  call MaterialPropertyDestroy(null_material_property)
  nullify(null_material_property)

  ! update ghosted values
  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%perm_xx_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%perm_yy_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%perm_zz_loc,ONEDOF)   
    
    if (option%mimetic) then
      call DiscretizationGlobalToLocal(discretization,field%perm0_xz, &
                                       field%perm_xz_loc,ONEDOF)  
      call DiscretizationGlobalToLocal(discretization,field%perm0_xy, &
                                       field%perm_xy_loc,ONEDOF)  
      call DiscretizationGlobalToLocal(discretization,field%perm0_yz, &
                                       field%perm_yz_loc,ONEDOF)   
    endif
     
    call DiscretizationLocalToLocal(discretization,field%icap_loc, &
                                    field%icap_loc,ONEDOF)   
    call DiscretizationLocalToLocal(discretization,field%ithrm_loc, &
                                    field%ithrm_loc,ONEDOF)
    call RealLocalToLocalWithArray(realization,SATURATION_FUNCTION_ID_ARRAY)
  endif
  
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%porosity_loc,ONEDOF)
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                   field%tortuosity_loc,ONEDOF)

end subroutine assignMaterialPropToRegions

! ************************************************************************** !
!
! verifyAllCouplers: Verifies the connectivity of a coupler
! author: Glenn Hammond
! date: 1/8/08
!
! ************************************************************************** !
subroutine verifyAllCouplers(realization)

  use Realization_class
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

  use Realization_class
  use Discretization_module
  use Option_module 
  use Coupler_module
  use Condition_module
  use Grid_module
  use Output_module
  use Output_Tecplot_module, only : OutputVectorTecplot  
  use Patch_module

  implicit none

  type(realization_type) :: realization
  type(coupler_list_type), pointer :: coupler_list

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(coupler_type), pointer :: coupler
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: iconn, icell, local_id
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option

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
    else if (len_trim(coupler%tran_condition_name) > 0) then
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

  use Realization_class
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
        if (region%grid_type == STRUCTURED_GRID_REGION) then
          call HDF5ReadRegionFromFile(realization,region,region%filename)
        else
          !geh: Do not skip this subouroutine if PETSC_HAVE_HDF5 is not 
          !     defined.  The subroutine prints an error message if not defined
          !     informing the user of the error.  If you skip the subroutine,
          !     no error message is printed and the user is unaware of the
          !     region not being read.
          call HDF5ReadUnstructuredGridRegionFromFile(realization%option,region, &
                                                      region%filename)
        endif
      else if (index(region%filename,'.ss') > 0) then
        region%sideset => RegionCreateSideset()
        call RegionReadFromFile(region%sideset,region%filename, &
                                realization%option)
      else if (index(region%filename,'.ex') > 0) then
        call RegionReadFromFile(region%explicit_faceset,region%cell_ids, &
                                region%filename,realization%option)
        region%num_cells = size(region%cell_ids)
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
subroutine readMaterialsFromFile(realization,realization_dependent,filename)

  use Realization_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_module

  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  PetscBool :: realization_dependent
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscBool :: append_realization_id
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  Vec :: global_vec
  Vec :: local_vec
  PetscErrorCode :: ierr

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization

  if (index(filename,'.h5') > 0) then
    group_name = 'Materials'
    dataset_name = 'Material Ids'
    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                    option)
    call HDF5ReadCellIndexedIntegerArray(realization,global_vec, &
                                         filename,group_name, &
                                         dataset_name,realization_dependent)
    call DiscretizationGlobalToLocal(discretization,global_vec,local_vec,ONEDOF)

    call GridCopyVecToIntegerArray(grid,patch%imat,local_vec,grid%ngmax)

    call VecDestroy(global_vec,ierr)
    call VecDestroy(local_vec,ierr)
  else
    call PetscLogEventBegin(logging%event_hash_map,ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP,filename,option)
    do
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ! natural ids in hash are zero-based
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        call InputReadInt(input,option,material_id)
        call InputErrorMsg(input,option,'material id','STRATA')
        patch%imat(ghosted_id) = material_id
      endif
    enddo
    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr)
  endif
  
end subroutine readMaterialsFromFile

! ************************************************************************** !
!
! readPermeabilitiesFromFile: Reads in grid cell permeabilities
! author: Glenn Hammond
! date: 01/19/09
!
! ************************************************************************** !
subroutine readPermeabilitiesFromFile(realization,material_property)

  use Realization_class
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Discretization_module
  use Logging_module
  use Input_module
  use Material_module
  use HDF5_module
  
  implicit none
  
  type(realization_type) :: realization
  type(material_property_type) :: material_property

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: local_id, ghosted_id, natural_id
  PetscReal :: permeability
  PetscBool :: append_realization_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscInt :: idirection
  PetscInt :: temp_int
  PetscReal :: ratio, scale
  Vec :: global_vec
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: perm_xx_p(:)
  PetscReal, pointer :: perm_yy_p(:)
  PetscReal, pointer :: perm_zz_p(:)
  PetscReal, pointer :: perm_xyz_p(:)

  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization

  call GridVecGetArrayF90(grid,field%perm0_xx,perm_xx_p,ierr)
  call GridVecGetArrayF90(grid,field%perm0_yy,perm_yy_p,ierr)
  call GridVecGetArrayF90(grid,field%perm0_zz,perm_zz_p,ierr)
  
  if (index(material_property%permeability_dataset%filename,'.h5') > 0) then
    group_name = ''
    if (material_property%permeability_dataset%realization_dependent) then
      append_realization_id = PETSC_TRUE
    else
      append_realization_id = PETSC_FALSE
    endif

    call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                    option)
    if (material_property%isotropic_permeability .or. &
        (.not.material_property%isotropic_permeability .and. &
         material_property%vertical_anisotropy_ratio > 0.d0)) then
      dataset_name = 'Permeability'
      call HDF5ReadCellIndexedRealArray(realization,global_vec, &
                          material_property%permeability_dataset%filename, &
                          group_name,dataset_name,append_realization_id)
      call GridVecGetArrayF90(grid,global_vec,vec_p,ierr)
      ratio = 1.d0
      scale = 1.d0
      !TODO(geh): fix so that ratio and scale work for perms outside
      ! of dataset
      if (material_property%vertical_anisotropy_ratio > 0.d0) then
        ratio = material_property%vertical_anisotropy_ratio
      endif
      if (material_property%permeability_scaling_factor > 0.d0) then
        scale = material_property%permeability_scaling_factor
      endif
      do local_id = 1, grid%nlmax
        if (patch%imat(grid%nL2G(local_id)) == material_property%id) then
          perm_xx_p(local_id) = vec_p(local_id)*scale
          perm_yy_p(local_id) = vec_p(local_id)*scale
          perm_zz_p(local_id) = vec_p(local_id)*ratio*scale
        endif
      enddo
      call GridVecRestoreArrayF90(grid,global_vec,vec_p,ierr)
    else
      temp_int = Z_DIRECTION
      if (grid%itype == STRUCTURED_GRID_MIMETIC) temp_int = YZ_DIRECTION
      do idirection = X_DIRECTION,temp_int
        select case(idirection)
          case(X_DIRECTION)
            dataset_name = 'PermeabilityX'
          case(Y_DIRECTION)
            dataset_name = 'PermeabilityY'
          case(Z_DIRECTION)
            dataset_name = 'PermeabilityZ'
          case(XY_DIRECTION)
            dataset_name = 'PermeabilityXY'
            call GridVecGetArrayF90(grid,field%perm0_xy,perm_xyz_p,ierr)
          case(XZ_DIRECTION)
            dataset_name = 'PermeabilityXZ'
            call GridVecGetArrayF90(grid,field%perm0_xz,perm_xyz_p,ierr)
          case(YZ_DIRECTION)
            dataset_name = 'PermeabilityYZ'
            call GridVecGetArrayF90(grid,field%perm0_yz,perm_xyz_p,ierr)
        end select          
        call HDF5ReadCellIndexedRealArray(realization,global_vec, &
                                          material_property%permeability_dataset%filename, &
                                          group_name, &
                                          dataset_name,append_realization_id)
        call GridVecGetArrayF90(grid,global_vec,vec_p,ierr)
        select case(idirection)
          case(X_DIRECTION)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == material_property%id) then
                perm_xx_p(local_id) = vec_p(local_id)
              endif
            enddo
          case(Y_DIRECTION)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == material_property%id) then
                perm_yy_p(local_id) = vec_p(local_id)
              endif
            enddo
          case(Z_DIRECTION)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == material_property%id) then
                perm_zz_p(local_id) = vec_p(local_id)
              endif
            enddo
          case(XY_DIRECTION,XZ_DIRECTION,YZ_DIRECTION)
            do local_id = 1, grid%nlmax
              if (patch%imat(grid%nL2G(local_id)) == material_property%id) then
                perm_xyz_p(local_id) = vec_p(local_id)
              endif
            enddo
            select case(idirection)
              case(XY_DIRECTION)
                call GridVecRestoreArrayF90(grid,field%perm0_xy,perm_xyz_p, &
                                            ierr)
              case(XZ_DIRECTION)
                call GridVecRestoreArrayF90(grid,field%perm0_xz,perm_xyz_p, &
                                            ierr)
              case(YZ_DIRECTION)
                call GridVecRestoreArrayF90(grid,field%perm0_yz,perm_xyz_p, &
                                            ierr)
            end select
        end select
        call GridVecRestoreArrayF90(grid,global_vec,vec_p,ierr)
      enddo
    endif
    call VecDestroy(global_vec,ierr)
  else

    call PetscLogEventBegin(logging%event_hash_map,ierr)
    call GridCreateNaturalToGhostedHash(grid,option)
    input => InputCreate(IUNIT_TEMP, &
                material_property%permeability_dataset%filename,option)
    do
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,natural_id)
      call InputErrorMsg(input,option,'natural id','STRATA')
      ghosted_id = GridGetLocalGhostedIdFromHash(grid,natural_id)
      if (ghosted_id > 0) then
        if (patch%imat(ghosted_id) /= material_property%id) cycle
        local_id = grid%nG2L(ghosted_id)
        if (local_id > 0) then
          call InputReadDouble(input,option,permeability)
          call InputErrorMsg(input,option,'permeability','STRATA')
          perm_xx_p(local_id) = permeability
          perm_yy_p(local_id) = permeability
          perm_zz_p(local_id) = permeability
        endif
      endif
    enddo

    call InputDestroy(input)
    call GridDestroyHashTable(grid)
    call PetscLogEventEnd(logging%event_hash_map,ierr)
  endif
  
  call GridVecRestoreArrayF90(grid,field%perm0_xx,perm_xx_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm0_yy,perm_yy_p,ierr)
  call GridVecRestoreArrayF90(grid,field%perm0_zz,perm_zz_p,ierr)
  
end subroutine readPermeabilitiesFromFile

! ************************************************************************** !
!
! readVectorFromFile: Reads data from a file into an associated vector
! author: Glenn Hammond
! date: 03/18/08
!
! ************************************************************************** !
subroutine readVectorFromFile(realization,vector,filename,vector_type)

  use Realization_class
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
  PetscErrorCode :: ierr
  PetscInt :: count, read_count, i
  PetscInt :: flag
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
    ! to be taken care of later in readPermeabilitiesFromFile()
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
      flag = ierr
      call MPI_Bcast(flag,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                     option%mycomm,ierr)      
      if (flag /= 0) then
        option%io_buffer = 'Insufficent data in file: ' // filename
        call printErrMsg(option)
      endif
      if (option%myrank == option%io_rank) then
        call VecSetValues(natural_vec,read_count,indices,values,INSERT_VALUES, &
                          ierr)
      endif
      count = count + read_count
    enddo
    call MPI_Bcast(count,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                   option%mycomm,ierr)      
    if (count /= grid%nmax) then
      write(option%io_buffer,'("Number of data in file (",i8, &
      & ") does not match size of vector (",i8,")")') count, grid%nlmax
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

! ************************************************************************** !
!
! readFlowInitialCondition: Assigns flow initial condition from HDF5 file
! author: Glenn Hammond
! date: 03/05/10
!
! ************************************************************************** !
subroutine readFlowInitialCondition(realization,filename)

  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Level_module
  use Patch_module
  use Discretization_module
  use HDF5_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: local_id, idx, offset
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  if (option%iflowmode /= RICHARDS_MODE) then
    option%io_buffer = 'Reading of flow initial conditions from HDF5 ' // &
                       'file (' // trim(filename) // &
                       'not currently not supported for mode: ' // &

                       trim(option%flowmode)
  endif      

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid

       ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid,field%flow_xx, xx_p, ierr); CHKERRQ(ierr)

      ! Pressure for all modes 
      offset = 1
      group_name = ''
      dataset_name = 'Pressure'
      call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                        filename,group_name, &
                                        dataset_name,option%id>0)
      call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
      do local_id=1, grid%nlmax
        if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
        if (dabs(vec_p(local_id)) < 1.d-40) then
          print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
               ': Potential error - zero pressure in Initial Condition read from file.'
        endif
        idx = (local_id-1)*option%nflowdof + offset
        xx_p(idx) = vec_p(local_id)
      enddo
      call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)

      call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
        
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)  
  call VecCopy(field%flow_xx, field%flow_yy, ierr)

end subroutine readFlowInitialCondition


! ************************************************************************** !
!
! readTransportInitialCondition: Assigns transport initial condition from 
!                                HDF5 file
! author: Glenn Hammond
! date: 03/05/10
!
! ************************************************************************** !
subroutine readTransportInitialCondition(realization,filename)

  use Realization_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Level_module
  use Reactive_Transport_module
  use Reaction_Aux_module
  use Discretization_module
  use HDF5_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: local_id, idx, offset, idof
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(reaction_type), pointer :: reaction

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid

       ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid,field%tran_xx,xx_p, ierr); CHKERRQ(ierr)

      ! Primary species concentrations for all modes 
      do idof = 1, option%ntrandof ! primary aqueous concentrations
        offset = idof
        group_name = ''
        dataset_name = reaction%primary_species_names(idof)
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          filename,group_name, &
                                          dataset_name,option%id>0)
        call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
        do local_id=1, grid%nlmax
          if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
          if (vec_p(local_id) < 1.d-40) then
            print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
              ': Zero free-ion concentration in Initial Condition read from file.'
          endif
          idx = (local_id-1)*option%ntrandof + offset
          xx_p(idx) = vec_p(local_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
     
      enddo     

      call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)
        
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr)
  
end subroutine readTransportInitialCondition

! ************************************************************************** !
!
! Create_IOGroups: Create sub-communicators that are used in initialization 
!                  and output HDF5 routines. 
! author: Vamsi Sripathi
! date: 07/14/09
!
! ************************************************************************** !
subroutine Create_IOGroups(option)

  use Option_module
  use Logging_module

#if defined(SCORPIO)
  use hdf5
#endif

  implicit none

  type(option_type) :: option
  PetscErrorCode :: ierr

#if defined(SCORPIO)

  PetscMPIInt :: numiogroups

  call PetscLogEventBegin(logging%event_create_iogroups,ierr)

  ! Initialize HDF interface to define global constants  
  call h5open_f(ierr)

  if (option%hdf5_read_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_READ_GROUP_SIZE & 
            & in the input file (pflotran.in) is either not set or &
            & its value is less than or equal to ZERO. &
            & HDF5_READ_GROUP_SIZE =  ",i6)') &
             option%hdf5_read_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let one process read and broadcast to everyone
    option%hdf5_read_group_size = option%mycommsize
  endif         
 
  if (option%hdf5_write_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_WRITE_GROUP_SIZE & 
            &in the input file (pflotran.in) is either not set or &
            &its value is less than or equal to ZERO. &
            &HDF5_WRITE_GROUP_SIZE =  ",i6)') &
             option%hdf5_write_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let everyone write separately 
    option%hdf5_write_group_size = 1
  endif                    

  ! create read IO groups
  numiogroups = option%mycommsize/option%hdf5_read_group_size
  call scorpio_iogroup_init(numiogroups, option%mycomm, option%ioread_group_id, ierr)

  if ( option%hdf5_read_group_size == option%hdf5_write_group_size ) then
    ! reuse read_group to use for writing too as both groups are same size
    option%iowrite_group_id = option%ioread_group_id
  else   
      ! create write IO groups
      numiogroups = option%mycommsize/option%hdf5_write_group_size
      call scorpio_iogroup_init(numiogroups, option%mycomm, option%iowrite_group_id, ierr)
  end if

    write(option%io_buffer, '(" Read group id :  ", i6)') option%ioread_group_id
    call printMsg(option)      
    write(option%io_buffer, '(" Write group id :  ", i6)') option%iowrite_group_id
    call printMsg(option)      
  call PetscLogEventEnd(logging%event_create_iogroups,ierr)
#endif
! SCORPIO
 
end subroutine Create_IOGroups

! ************************************************************************** !
!
! InitPrintPFLOTRANHeader: Initializes pflotran
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine InitPrintPFLOTRANHeader(option,fid)

  use Option_module
  
  implicit none
  
  PetscInt :: fid
  
  type(option_type) :: option
  
  write(fid,'(" PFLOTRAN Header")') 
  
end subroutine InitPrintPFLOTRANHeader

! ************************************************************************** !
!
! InitReadVelocityField: Reads fluxes in for transport with no flow.
! author: Glenn Hammond
! date: 02/05/13
!
! ************************************************************************** !
subroutine InitReadVelocityField(realization)

  use Realization_class
  use Patch_module
  use Field_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Discretization_module
  use HDF5_module

  implicit none
  
  type(realization_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: idir, iconn, sum_connection
  PetscInt :: ghosted_id_up, local_id
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_loc_p(:)
  PetscReal, pointer :: vec_p(:)
  type(coupler_type), pointer :: boundary_condition  
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization
  
  filename = realization%nonuniform_velocity_filename

  group_name = ''
  do idir = 1, 3
    select case(idir)
      case(1)
        dataset_name = 'Internal Velocity X'
      case(2)
        dataset_name = 'Internal Velocity Y'
      case(3)
        dataset_name = 'Internal Velocity Z'
    end select
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call DiscretizationGlobalToLocal(discretization,field%work,field%work_loc, &
                                     ONEDOF)
    call GridVecGetArrayF90(grid,field%work_loc,vec_loc_p,ierr)
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0  
    do 
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        ghosted_id_up = cur_connection_set%id_up(iconn)
        if (cur_connection_set%dist(idir,iconn) > 0.9d0) then
          patch%internal_velocities(1,sum_connection) = vec_loc_p(ghosted_id_up)
        endif
      enddo
      cur_connection_set => cur_connection_set%next
    enddo
    call GridVecRestoreArrayF90(grid,field%work_loc,vec_loc_p,ierr)
  enddo
  
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    dataset_name = boundary_condition%name
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      patch%boundary_velocities(1,sum_connection) = vec_p(local_id)
    enddo
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine InitReadVelocityField
            
end module Init_module
