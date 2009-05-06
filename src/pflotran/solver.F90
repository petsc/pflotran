module Solver_module
 
  implicit none

  private
 
#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"
#include "finclude/petscmg.h"

  type, public :: solver_type
    PetscInt :: itype            ! type: flow or transport
    PetscReal :: linear_atol       ! absolute tolerance
    PetscReal :: linear_rtol       ! relative tolerance
    PetscReal :: linear_dtol       ! divergence tolerance
    PetscInt :: linear_maxit     ! maximum number of iterations
    
    PetscReal :: newton_atol       ! absolute tolerance
    PetscReal :: newton_rtol       ! relative tolerance
    PetscReal :: newton_stol       ! relative tolerance (relative to previous iteration)
    PetscReal :: newton_dtol       ! divergence tolerance
    PetscReal :: newton_inf_res_tol    ! infinity tolerance for residual
    PetscReal :: newton_inf_upd_tol    ! infinity tolerance for update
    PetscInt :: newton_maxit     ! maximum number of iterations
    PetscInt :: newton_maxf      ! maximum number of function evaluations

    PetscTruth :: use_galerkin_mg  ! If true, precondition linear systems with 
                                   ! Galerkin-type geometric multigrid.
    PetscInt :: galerkin_mg_levels  ! Number of discretization levels for 
                                    ! the Galerkin MG (includes finest level).
    PetscInt :: galerkin_mg_levels_x
    PetscInt :: galerkin_mg_levels_y
    PetscInt :: galerkin_mg_levels_z

    ! Jacobian matrix
    Mat :: J    ! Jacobian
    Mat :: Jpre ! Jacobian to be used in preconditioner
    MatType :: J_mat_type
    MatType :: Jpre_mat_type

    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

    Mat, pointer :: interpolation(:)
      ! Hierarchy of interpolation operators for Galerkin multigrid.

    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType  :: pc_type
    KSP   ::  ksp
    PC    ::  pc
    
    PetscTruth :: inexact_newton

    PetscTruth :: print_convergence
    PetscTruth :: print_detailed_convergence
    PetscTruth :: print_linear_iterations
    PetscTruth :: check_infinity_norm
    PetscTruth :: force_at_least_1_iteration    
            
  end type solver_type
  
  public :: SolverCreate, &
            SolverDestroy, &
            SolverReadLinear, &
            SolverReadNewton, &
            SolverCreateSNES, &
            SolverSetSNESOptions, &
            SolverPrintNewtonInfo, &
            SolverPrintLinearInfo, &
            SolverCheckCommandLine
  
contains

! ************************************************************************** !
!
! SolverCreate: Allocates and initializes a new (empty) Solver object
! Note that this does not create the PETSc solver contexts associated 
! with the Solver.  These contexts are created via a subsequent call to 
! SolverCreateSNES().
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SolverCreate()

  implicit none
  
  type(solver_type), pointer :: SolverCreate
  
  type(solver_type), pointer :: solver
  
  allocate(solver)
  
  ! initialize to default values
  solver%itype = NULL_CLASS
  solver%linear_atol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%linear_rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%linear_dtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%linear_maxit = PETSC_DEFAULT_INTEGER
  
  solver%newton_atol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%newton_rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%newton_stol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%newton_dtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%newton_inf_res_tol = 1.d-50 ! arbitrarily set by geh
  solver%newton_inf_upd_tol = 1.d-50 ! arbitrarily set by geh
  solver%newton_maxit = PETSC_DEFAULT_INTEGER
  solver%newton_maxf = PETSC_DEFAULT_INTEGER

  solver%use_galerkin_mg = PETSC_FALSE
  solver%galerkin_mg_levels = 1
  solver%galerkin_mg_levels_x = 1
  solver%galerkin_mg_levels_y = 1
  solver%galerkin_mg_levels_z = 1
  
  solver%J = 0
  solver%Jpre = 0
  solver%J_mat_type = MATBAIJ
  solver%Jpre_mat_type = ''
!  solver%interpolation = 0
  nullify(solver%interpolation)
  solver%matfdcoloring = 0
  solver%snes = 0
  solver%ksp_type = KSPBCGS
  solver%pc_type = ""
  solver%ksp = 0
  solver%pc = 0
  
  solver%inexact_newton = PETSC_FALSE
  
  solver%print_convergence = PETSC_TRUE
  solver%print_detailed_convergence = PETSC_FALSE
  solver%print_linear_iterations = PETSC_FALSE
  solver%check_infinity_norm = PETSC_TRUE
  solver%force_at_least_1_iteration = PETSC_TRUE
    
  SolverCreate => solver
  
end function SolverCreate

! ************************************************************************** !
!
! SolverCreateSNES: Create PETSc SNES object
! author: Glenn Hammond
! date: 02/12/08
!
! ************************************************************************** !
subroutine SolverCreateSNES(solver,comm)

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call SNESCreate(comm,solver%snes,ierr)
  call SNESSetFromOptions(solver%snes,ierr) 

  ! grab handles for ksp and pc
  call SNESGetKSP(solver%snes,solver%ksp,ierr)
  call KSPGetPC(solver%ksp,solver%pc,ierr)

end subroutine SolverCreateSNES
  
! ************************************************************************** !
!
! SolverSetSNESOptions: Sets options for SNES
! author: Glenn Hammond
! date: 02/12/08
!
! ************************************************************************** !
subroutine SolverSetSNESOptions(solver)

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: myrank
  ! needed for SNESLineSearchGetParams()/SNESLineSearchSetParams()
  PetscReal :: alpha, maxstep, steptol
  PetscErrorCode :: ierr
  PetscInt :: i
  
  ! if ksp_type or pc_type specified in input file, set them here
  if (len_trim(solver%ksp_type) > 1) &
    call KSPSetType(solver%ksp,solver%ksp_type,ierr)
  if (len_trim(solver%pc_type) > 1) &
    call PCSetType(solver%pc,solver%pc_type,ierr)

  call KSPSetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                        solver%linear_dtol,solver%linear_maxit,ierr)

  ! allow override from command line
  call KSPSetFromOptions(solver%ksp,ierr)
  call PCSetFromOptions(solver%pc,ierr)
    
  ! get the ksp_type and pc_type incase of command line override.
  call KSPGetType(solver%ksp,solver%ksp_type,ierr)
  call PCGetType(solver%pc,solver%pc_type,ierr)

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(solver%snes, solver%newton_atol, solver%newton_rtol, &
                         solver%newton_stol,solver%newton_maxit, &
                         solver%newton_maxf,ierr)

  ! set inexact newton, currently applies default settings
  if (solver%inexact_newton) &
    call SNESKSPSetUseEW(solver%snes,PETSC_TRUE,ierr)

!  call SNESLineSearchSet(solver%snes,SNESLineSearchNo,PETSC_NULL)

  ! Setup for n-level Galerkin multigrid.
  if (solver%use_galerkin_mg) then
    call PCSetType(solver%pc, PCMG,ierr)
    call PCMGSetLevels(solver%pc, solver%galerkin_mg_levels, &
                       PETSC_NULL_OBJECT,ierr)
    do i=1,solver%galerkin_mg_levels-1
      call PCMGSetInterpolation(solver%pc, i, solver%interpolation(i),ierr)
      call PCMGSetGalerkin(solver%pc,ierr)
    enddo
  endif
  
  ! allow override from command line; for some reason must come before
  ! LineSearchParams, or they crash
  call SNESSetFromOptions(solver%snes,ierr) 
    
  call SNESLineSearchGetParams(solver%snes, alpha, maxstep, steptol,ierr)  
  call SNESLineSearchSetParams(solver%snes, alpha, maxstep, solver%newton_stol,ierr)  

  call SNESGetTolerances(solver%snes,solver%newton_atol,solver%newton_rtol, &
                         solver%newton_stol,solver%newton_maxit, &
                         solver%newton_maxf,ierr)

  call KSPGetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                         solver%linear_dtol,solver%linear_maxit,ierr)

end subroutine SolverSetSNESOptions
  
! ************************************************************************** !
!
! SolverReadLinear: Reads parameters associated with linear solver
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine SolverReadLinear(solver,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type) :: input
  type(option_type) :: option
  PetscErrorCode :: ierr
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2, prefix
  character(len=MAXSTRINGLENGTH) :: string

  select case(solver%itype)
    case(FLOW_CLASS)
      prefix = '-flow_'
    case(TRANSPORT_CLASS)
      prefix = '-tran_'
  end select

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','LINEAR SOLVER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('SOLVER_TYPE','SOLVER','KRYLOV_TYPE','KRYLOV','KSP','KSP_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'ksp_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PREONLY')
            solver%ksp_type = KSPPREONLY
          case('GMRES')
            solver%ksp_type = KSPGMRES
          case('FGMRES')
            solver%ksp_type = KSPFGMRES
          case('BCGS','BICGSTAB','BI-CGSTAB')
            solver%ksp_type = KSPBCGS
          case('IBCGS','IBICGSTAB','IBI-CGSTAB')
            solver%ksp_type = KSPIBCGS
          case('RICHARDSON')
            solver%ksp_type = KSPRICHARDSON
          case('CG')
            solver%ksp_type = KSPCG
          case default
            option%io_buffer  = 'Krylov solver type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select

      case('PRECONDITIONER_TYPE','PRECONDITIONER','PC','PC_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'pc_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PCNONE')
            solver%pc_type = PCNONE
          case('ILU','PCILU')
            solver%pc_type = PCILU
          case('LU','PCLU')
            solver%pc_type = PCLU
          case('BJACOBI','BLOCK_JACOBI')
            solver%pc_type = PCBJACOBI
          case('ASM','ADDITIVE_SCHWARZ')
            solver%pc_type = PCASM
         case('HYPRE')
            solver%pc_type = PCHYPRE
         case('SHELL')
            solver%pc_type = PCSHELL
          case default
            option%io_buffer  = 'Preconditioner type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select

      case('HYPRE_OPTIONS')
        do
          call InputReadFlotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          call InputReadWord(input,option,keyword,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','LINEAR SOLVER, HYPRE options')   
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'type','LINEAR SOLVER, HYPRE options')  
              call StringToLower(word)
              select case(trim(word))
                case('pilut','parasails','boomeramg','euclid')
                  string = trim(prefix) // 'pc_hypre_type'
                  call PetscOptionsSetValue(trim(string),trim(word),ierr); 
                case default
                  option%io_buffer  = 'HYPRE preconditioner type: ' // &
                                      trim(word) // ' unknown.'
                  call printErrMsg(option)
              end select
            case('BOOMERAMG_CYCLE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG cycle type','LINEAR SOLVER, HYPRE options')  
              call StringToLower(word)
              string = trim(prefix) // 'pc_hypre_boomeramg_cycle_type'
              select case(trim(word))
                case('V')
                  call PetscOptionsSetValue(trim(string),'1',ierr); 
                case('W')
                  call PetscOptionsSetValue(trim(string),'2',ierr); 
                case default
                  option%io_buffer  = 'HYPRE BoomerAMG cycle type: ' &
                                      // trim(word) // ' unknown.'
                  call printErrMsg(option)
              end select
            case('BOOMERAMG_MAX_LEVELS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_levels'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_MAX_ITER')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum iterations', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_iter'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_TOL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG convergence tolerance', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_tol'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_TRUNCFACTOR')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG interpolation truncation factor', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_truncfactor'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_AGG_NL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG # levels aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_nl'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_AGG_NUM_PATHS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG # paths for aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_num_paths'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_STRONG_THRESHOLD')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG threshold for strong connectivity', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_strong_threshold'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_GRID_SWEEPS_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG number of grid sweeps up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_all'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_GRID_SWEEPS_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG number of grid sweeps down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_down'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_GRID_SWEEPS_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG number of grid sweeps up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_up'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_GRID_SWEEPS_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG number of grid sweeps for coarse level', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_coarse'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_TYPE_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation type for up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_all'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_TYPE_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation type for down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_down'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_TYPE_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation type for up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_up'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_TYPE_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation type for coarse grids', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_coarse'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_all'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_level'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG outer relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_outer_relax_weight_all'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG outer relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // 'pc_hypre_boomeramg_outer_relax_weight_level'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_NO_CF')
              string = trim(prefix) // 'pc_hypre_boomeramg_no_CF'
              call PetscOptionsSetValue(trim(string),'',ierr); 
            case('BOOMERAMG_MEASURE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG measure type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_measure_type'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_COARSEN_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG coarsen type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_coarsen_type'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_INTERPOLATION_TYPE','BOOMERAMG_INTERP_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG interpolation type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_interp_type'
              call PetscOptionsSetValue(trim(string),trim(word),ierr); 
            case('BOOMERAMG_NODAL_COARSEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG set nodal coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_coarsen'
              call PetscOptionsSetValue(trim(string),'',ierr); 
            case('BOOMERAMG_NODAL_RELAXATION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG nodal relaxation via Schwarz', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_relaxation'
              call PetscOptionsSetValue(trim(string),'',ierr); 
            case default
              option%io_buffer  = 'HYPRE option: ' // trim(keyword) // ' unknown.'
              call printErrMsg(option)
          end select
        enddo

      case('ATOL')
        call InputReadDouble(input,option,solver%linear_atol)
        call InputDefaultMsg(input,option,'linear_atol')

      case('RTOL')
        call InputReadDouble(input,option,solver%linear_rtol)
        call InputDefaultMsg(input,option,'linear_rtol')

      case('DTOL')
        call InputReadDouble(input,option,solver%linear_dtol)
        call InputDefaultMsg(input,option,'linear_dtol')
   
      case('MAXIT')
        call InputReadInt(input,option,solver%linear_maxit)
        call InputDefaultMsg(input,option,'linear_maxit')

      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in linear solver'    
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine SolverReadLinear

! ************************************************************************** !
!
! SolverReadNewton: Reads parameters associated with linear solver
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine SolverReadNewton(solver,input,option)

  use Input_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','NEWTON SOLVER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case ('INEXACT_NEWTON')
        solver%inexact_newton = PETSC_TRUE

      case ('NO_PRINT_CONVERGENCE')
        solver%print_convergence = PETSC_FALSE

      case ('NO_INF_NORM','NO_INFINITY_NORM')
        solver%check_infinity_norm = PETSC_FALSE

      case ('NO_FORCE_ITERATION')
        solver%force_at_least_1_iteration = PETSC_FALSE

      case ('PRINT_DETAILED_CONVERGENCE')
        solver%print_detailed_convergence = PETSC_TRUE

      case ('PRINT_LINEAR_ITERATIONS')
        solver%print_linear_iterations = PETSC_TRUE

      case('ATOL')
        call InputReadDouble(input,option,solver%newton_atol)
        call InputDefaultMsg(input,option,'newton_atol')

      case('RTOL')
        call InputReadDouble(input,option,solver%newton_rtol)
        call InputDefaultMsg(input,option,'newton_rtol')

      case('STOL')
        call InputReadDouble(input,option,solver%newton_stol)
        call InputDefaultMsg(input,option,'newton_stol')
      
      case('DTOL')
        call InputReadDouble(input,option,solver%newton_dtol)
        call InputDefaultMsg(input,option,'newton_dtol')
   
      case('ITOL', 'INF_TOL', 'ITOL_RES', 'INF_TOL_RES')
        call InputReadDouble(input,option,solver%newton_inf_res_tol)
        call InputDefaultMsg(input,option,'newton_inf_res_tol')
   
      case('ITOL_UPDATE', 'INF_TOL_UPDATE')
        call InputReadDouble(input,option,solver%newton_inf_upd_tol)
        call InputDefaultMsg(input,option,'newton_inf_upd_tol')
   
      case('MAXIT')
        call InputReadInt(input,option,solver%newton_maxit)
        call InputDefaultMsg(input,option,'newton_maxit')

      case('MAXF')
        call InputReadInt(input,option,solver%newton_maxf)
        call InputDefaultMsg(input,option,'newton_maxf')

      case('MATRIX_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('BAIJ')
            solver%J_mat_type = MATBAIJ
          case('AIJ')
            solver%J_mat_type = MATBAIJ
          case('MFFD','MATRIX_FREE')
            solver%J_mat_type = MATMFFD
          case('HYPRESTRUCT')
!            solver%J_mat_type = MATHYPRESTRUCT
          case default
            option%io_buffer = 'Matrix type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select
        
      case('PRECONDITIONER_MATRIX_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('BAIJ')
            solver%Jpre_mat_type = MATBAIJ
          case('AIJ')
            solver%Jpre_mat_type = MATBAIJ
          case('MFFD','MATRIX_FREE')
            solver%Jpre_mat_type = MATMFFD
          case('HYPRESTRUCT')
 !           solver%Jpre_mat_type = MATHYPRESTRUCT
          case default
            option%io_buffer  = 'Preconditioner Matrix type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select
        
      case default
        option%io_buffer = 'Keyword: '//keyword// &
                           ' not recognized in Newton solver'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine SolverReadNewton

! ************************************************************************** !
!
! SolverPrintLinearInfo: Prints information about linear solver
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine SolverPrintLinearInfo(solver,print_to_screen,print_to_file,fid, &
                                 header)

  implicit none
  
  type(solver_type) :: solver
  PetscTruth :: print_to_screen
  PetscTruth :: print_to_file
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header  

  if (print_to_screen) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(*,'(" solver: ",a)') trim(solver%ksp_type)
    write(*,'("precond: ",a)') trim(solver%pc_type)
    write(*,'("   atol:",1pe12.4)') solver%linear_atol
    write(*,'("   rtol:",1pe12.4)') solver%linear_rtol
    write(*,'("   dtol:",1pe12.4)') solver%linear_dtol
    write(*,'("  maxit:",i7)') solver%linear_maxit
  endif
  
  if (print_to_file) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(fid,'(" solver: ",a)') trim(solver%ksp_type)
    write(fid,'("precond: ",a)') trim(solver%pc_type)
    write(fid,'("   atol:",1pe12.4)') solver%linear_atol
    write(fid,'("   rtol:",1pe12.4)') solver%linear_rtol
    write(fid,'("   dtol:",1pe12.4)') solver%linear_dtol
    write(fid,'("  maxit:",i7)') solver%linear_maxit
  endif

end subroutine SolverPrintLinearInfo


! ************************************************************************** !
!
! SolverPrintNewtonInfo: Prints information about Newton solver
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine SolverPrintNewtonInfo(solver,print_to_screen,print_to_file,fid, &
                                 header)    

  implicit none
  
  type(solver_type) :: solver
  PetscTruth :: print_to_screen
  PetscTruth :: print_to_file  
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header  

  if (print_to_screen) then
    write(*,*) 
    write(*,'(a)') trim(header)
    write(*,'("     atol:",1pe12.4)') solver%newton_atol
    write(*,'("     rtol:",1pe12.4)') solver%newton_rtol
    write(*,'("     stol:",1pe12.4)') solver%newton_stol
    write(*,'("     dtol:",1pe12.4)') solver%newton_dtol
    write(*,'("inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(*,'("inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(*,'("    maxit:",i6)') solver%newton_maxit
    write(*,'("     maxf:",i6)') solver%newton_maxf
    if (solver%inexact_newton) then
      write(*,'("inexact newton: on")')
    else
      write(*,'("inexact newton: off")')
    endif
        
    if (solver%print_convergence) then
      write(*,'("print convergence: on")')
    else
      write(*,'("print convergence: off")')
    endif
        
    if (solver%print_detailed_convergence) then
      write(*,'("print detailed convergence: on")')
    else
      write(*,'("print detailed convergence: off")')
    endif
        
    if (solver%check_infinity_norm) then
      write(*,'("check infinity norm: on")')
    else
      write(*,'("check infinity norm: off")')
    endif
        
    if (solver%force_at_least_1_iteration) then
      write(*,'("force at least 1 iteration: on")')
    else
      write(*,'("force at least 1 iteration: off")')
    endif
  endif

  if (print_to_file) then
    write(fid,*) 
    write(fid,'(a)') trim(header)
    write(fid,'("     atol:",1pe12.4)') solver%newton_atol
    write(fid,'("     rtol:",1pe12.4)') solver%newton_rtol
    write(fid,'("     stol:",1pe12.4)') solver%newton_stol
    write(fid,'("     dtol:",1pe12.4)') solver%newton_dtol
    write(fid,'("inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(fid,'("inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(fid,'("    maxit:",i6)') solver%newton_maxit
    write(fid,'("     maxf:",i6)') solver%newton_maxf
    if (solver%inexact_newton) then
      write(fid,'("inexact newton: on")')
    else
      write(fid,'("inexact newton: off")')
    endif
        
    if (solver%print_convergence) then
      write(fid,'("print convergence: on")')
    else
      write(fid,'("print convergence: off")')
    endif
        
    if (solver%print_detailed_convergence) then
      write(fid,'("print detailed convergence: on")')
    else
      write(fid,'("print detailed convergence: off")')
    endif
        
    if (solver%check_infinity_norm) then
      write(fid,'("check infinity norm: on")')
    else
      write(fid,'("check infinity norm: off")')
    endif
        
    if (solver%force_at_least_1_iteration) then
      write(fid,'("force at least 1 iteration: on")')
    else
      write(fid,'("force at least 1 iteration: off")')
    endif
  endif

end subroutine SolverPrintNewtonInfo

! ************************************************************************** !
!
! SolverCheckCommandLine: Parses the command line for various solver 
! options.
! Note: In order to use the PETSc OptionsPrefix associated with 
! solver%snes in parsing the options, the call to SolverCheckCommandLine() 
! should come after the SNESSetOptionsPrefix(solver%snes,...) call.
! author: Richard Tran Mills
! date: 05/09/2008
!
! ************************************************************************** !
subroutine SolverCheckCommandLine(solver)

  implicit none
  
  type(solver_type) :: solver

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: prefix
  character(len=MAXSTRINGLENGTH) :: mat_type
  PetscTruth :: is_present

  if (solver%snes /= 0) then
    call SNESGetOptionsPrefix(solver%snes, prefix, ierr)
  else
    prefix = PETSC_NULL_CHARACTER
  endif

  ! Parse the options to determine if the matrix type has been specified.
  call PetscOptionsGetString(prefix, '-mat_type', mat_type, is_present,ierr)
  if (is_present) solver%J_mat_type = mat_type
  
  call PetscOptionsGetString(prefix, '-pre_mat_type', mat_type, is_present,ierr)
  if (is_present) solver%Jpre_mat_type = mat_type

  ! Parse the options for the Galerkin multigrid solver.
  ! Users can specify the number of levels of coarsening via the 
  ! 'galerkin_mg N' option, which will set the number of levels in the 
  ! x, y, and z directions all to N.  For semi-coarsening, however, 
  ! it is possible to set the number of levels in each direction 
  ! individually via options such as '-galerkin_mg_x N', which would 
  ! override the number of levels in the x direction set by '-galerkin_mg'.
  call PetscOptionsGetInt(prefix, '-galerkin_mg', &
                          solver%galerkin_mg_levels, solver%use_galerkin_mg, &
                          ierr)
  if (solver%use_galerkin_mg) then
    solver%galerkin_mg_levels_x = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_y = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_z = solver%galerkin_mg_levels
  endif

  call PetscOptionsGetInt(prefix, '-galerkin_mg_x', &
                          solver%galerkin_mg_levels_x, is_present,ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(prefix, '-galerkin_mg_y', &
                          solver%galerkin_mg_levels_y, is_present,ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(prefix, '-galerkin_mg_z', &
                          solver%galerkin_mg_levels_z, is_present,ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE

  if (solver%use_galerkin_mg) then
    solver%J_mat_type = MATAIJ
      ! Must use AIJ above, as BAIJ is not supported for Galerkin MG solver.
    solver%galerkin_mg_levels = max(solver%galerkin_mg_levels_x, &
                                    solver%galerkin_mg_levels_y, &
                                    solver%galerkin_mg_levels_z)
  endif
                             

end subroutine SolverCheckCommandLine

! ************************************************************************** !
!
! SolverDestroy: Deallocates a solver
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SolverDestroy(solver)

  implicit none
  
  type(solver_type), pointer :: solver
  
  PetscErrorCode :: ierr
  integer :: i

  if (.not.associated(solver)) return

  if (solver%Jpre == solver%J) then
    solver%Jpre = 0
  elseif (solver%Jpre /= 0) then
    call MatDestroy(solver%Jpre)
  endif
  if (solver%J /= 0) call MatDestroy(solver%J,ierr)
  if (associated(solver%interpolation)) then
    do i=1,solver%galerkin_mg_levels-1
      call MatDestroy(solver%interpolation(i),ierr)
    enddo
    deallocate(solver%interpolation)
  endif
  if (solver%matfdcoloring /= 0) call MatFDColoringDestroy(solver%matfdcoloring,ierr)
! if (solver%snes /= 0) call SNESDestroy(solver%snes)
  solver%ksp = 0
  solver%pc = 0
    
  deallocate(solver)
  nullify(solver)
  
end subroutine SolverDestroy
  
end module Solver_module
