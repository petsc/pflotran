module Solver_module
 
  implicit none

  private
 
#include "definitions.h"

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscmg.h"

  type, public :: solver_type
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

    ! Jacobian matrix
    Mat :: J
    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

    Mat :: interpolation ! Interpolation operator for Galerkin multigrid.

    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType  :: pc_type
    KSP   ::  ksp
    PC    ::  pc
    
    PetscTruth :: inexact_newton

    PetscTruth :: print_convergence
    PetscTruth :: print_detailed_convergence
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
            SolverPrintLinearInfo
  
contains

! ************************************************************************** !
!
! SolverCreate: Allocates and initializes a new Solver object
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
  
  solver%J = 0
  solver%interpolation = 0
  solver%matfdcoloring = 0
  solver%snes = 0
  solver%ksp_type = KSPBCGS
  solver%pc_type = ""
  solver%ksp = 0
  solver%pc = 0
  
  solver%inexact_newton = PETSC_FALSE
  
  solver%print_convergence = PETSC_TRUE
  solver%print_detailed_convergence = PETSC_FALSE
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
subroutine SolverCreateSNES(solver)

  implicit none
  
  type(solver_type) :: solver

  PetscErrorCode :: ierr
  
  call SNESCreate(PETSC_COMM_WORLD, solver%snes, ierr)
  call SNESSetFromOptions(solver%snes, ierr) 

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
subroutine SolverSetSNESOptions(solver, option)

  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type), pointer :: option

  PetscMPIInt :: myrank
  ! needed for SNESLineSearchGetParams()/SNESLineSearchSetParams()
  PetscReal :: alpha, maxstep, steptol
  PetscErrorCode :: ierr
  
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
                         solver%newton_maxf, ierr)

  ! set inexact newton, currently applies default settings
  if (solver%inexact_newton) &
    call SNESKSPSetUseEW(solver%snes,PETSC_TRUE,ierr)

!  call SNESLineSearchSet(solver%snes,SNESLineSearchNo,PETSC_NULL,ierr)

  ! Setup for 2-level Galerkin multigrid.
  if (option%use_galerkin_mg) then
    call PCSetType(solver%pc, PCMG, ierr)
    call PCMGSetLevels(solver%pc, 2, PETSC_NULL_OBJECT, ierr)
    call PCMGSetInterpolation(solver%pc, 1, solver%interpolation, ierr)
    call PCMGSetGalerkin(solver%pc, ierr)
  endif
  
  ! allow override from command line; for some reason must come before
  ! LineSearchParams, or they crash
  call SNESSetFromOptions(solver%snes, ierr) 
    
  call SNESLineSearchGetParams(solver%snes, alpha, maxstep, steptol, ierr)  
  call SNESLineSearchSetParams(solver%snes, alpha, maxstep, solver%newton_stol, ierr)  

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
subroutine SolverReadLinear(solver,fid,myrank)

  use Fileio_module
  
  implicit none

  type(solver_type) :: solver
  PetscInt :: fid
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(myrank,'keyword','LINEAR SOLVER', ierr)
    call fiWordToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('SOLVER_TYPE','SOLVER','KRYLOV_TYPE','KRYLOV','KSP','KSP_TYPE')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(myrank,'ksp_type','SOLVER', ierr)   
        call fiWordToUpper(word)
        select case(trim(word))
          case('NONE','PREONLY')
            solver%ksp_type = KSPPREONLY
          case('GMRES')
            solver%ksp_type = KSPGMRES
          case('BCGS','BICGSTAB','BI-CGSTAB')
            solver%ksp_type = KSPBCGS
          case default
            string  = 'ERROR: Krylov solver type: ' // trim(word) // ' unknown.'
            if (myrank == 0) print *, string
            stop
        end select
      case('PRECONDITIONER_TYPE','PRECONDITIONER','PC','PC_TYPE')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(myrank,'pc_type','SOLVER', ierr)   
        call fiWordToUpper(word)
        select case(trim(word))
          case('ILU','PCILU')
            solver%pc_type = PCILU
          case('LU','PCLU')
            solver%pc_type = PCLU
          case('BJACOBI','BLOCK_JACOBI')
            solver%pc_type = PCBJACOBI
          case('ASM','ADDITIVE_SCHWARTZ')
            solver%pc_type = PCASM
          case default
            string  = 'ERROR: Preconditioner type: ' // trim(word) // ' unknown.'
            if (myrank == 0) print *, string
            stop
        end select

      case('ATOL')
        call fiReadDouble(string,solver%linear_atol,ierr)
        call fiDefaultMsg(myrank,'linear_atol',ierr)

      case('RTOL')
        call fiReadDouble(string,solver%linear_rtol,ierr)
        call fiDefaultMsg(myrank,'linear_rtol',ierr)

      case('DTOL')
        call fiReadDouble(string,solver%linear_dtol,ierr)
        call fiDefaultMsg(myrank,'linear_dtol',ierr)
   
      case('MAXIT')
        call fiReadInt(string,solver%linear_maxit,ierr)
        call fiDefaultMsg(myrank,'linear_maxit',ierr)

      case default
        if (myrank == 0) print *, 'Keyword: '//keyword// &
                                  &' not recognized in linear solver'    

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
subroutine SolverReadNewton(solver,fid,myrank)

  use Fileio_module
  
  implicit none

  type(solver_type) :: solver
  PetscInt :: fid
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)

    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(myrank,'keyword','NEWTON SOLVER', ierr)
    call fiWordToUpper(keyword)   
      
    select case(trim(keyword))
    
      case ('INEXACT_NEWTON')
        solver%inexact_newton = .true.

      case ('NO_PRINT_CONVERGENCE')
        solver%print_convergence = PETSC_FALSE

      case ('NO_INF_NORM','NO_INFINITY_NORM')
        solver%check_infinity_norm = PETSC_FALSE

      case ('NO_FORCE_ITERATION')
        solver%force_at_least_1_iteration = PETSC_FALSE

      case ('PRINT_DETAILED_CONVERGENCE')
        solver%print_detailed_convergence = PETSC_TRUE

      case('ATOL')
        call fiReadDouble(string,solver%newton_atol,ierr)
        call fiDefaultMsg(myrank,'newton_atol',ierr)

      case('RTOL')
        call fiReadDouble(string,solver%newton_rtol,ierr)
        call fiDefaultMsg(myrank,'newton_rtol',ierr)

      case('STOL')
        call fiReadDouble(string,solver%newton_stol,ierr)
        call fiDefaultMsg(myrank,'newton_stol',ierr)
      
      case('DTOL')
        call fiReadDouble(string,solver%newton_dtol,ierr)
        call fiDefaultMsg(myrank,'newton_dtol',ierr)
   
      case('ITOL', 'INF_TOL', 'ITOL_RES', 'INF_TOL_RES')
        call fiReadDouble(string,solver%newton_inf_res_tol,ierr)
        call fiDefaultMsg(myrank,'newton_inf_res_tol',ierr)
   
      case('ITOL_UPDATE', 'INF_TOL_UPDATE')
        call fiReadDouble(string,solver%newton_inf_upd_tol,ierr)
        call fiDefaultMsg(myrank,'newton_inf_upd_tol',ierr)
   
      case('MAXIT')
        call fiReadInt(string,solver%newton_maxit,ierr)
        call fiDefaultMsg(myrank,'newton_maxit',ierr)

      case('MAXF')
        call fiReadInt(string,solver%newton_maxf,ierr)
        call fiDefaultMsg(myrank,'newton_maxf',ierr)

      case default
        if (myrank == 0) print *, 'Keyword: '//keyword// &
                                  &' not recognized in Newton solver'
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
subroutine SolverPrintLinearInfo(solver,fid,header,myrank)

  implicit none
  
  type(solver_type) :: solver
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header  
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string

  if (myrank == 0) then
    write(*,*) 
    write(fid,*) 
    write(*,'(a)') trim(header)
    write(fid,'(a)') trim(header)
    write(*,'(" atol:",1pe12.4)') solver%linear_atol
    write(fid,'(" atol:",1pe12.4)') solver%linear_atol
    write(*,'(" rtol:",1pe12.4)') solver%linear_rtol
    write(fid,'(" rtol:",1pe12.4)') solver%linear_rtol
    write(*,'(" dtol:",1pe12.4)') solver%linear_dtol
    write(fid,'(" dtol:",1pe12.4)') solver%linear_dtol
    write(*,'("maxit:",i7)') solver%linear_maxit
    write(fid,'("maxit:",i7)') solver%linear_maxit
  endif

end subroutine SolverPrintLinearInfo


! ************************************************************************** !
!
! SolverPrintNewtonInfo: Prints information about Newton solver
! author: Glenn Hammond
! date: 02/23/08
!
! ************************************************************************** !
subroutine SolverPrintNewtonInfo(solver,fid,header,myrank)    

  implicit none
  
  type(solver_type) :: solver
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: header  
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string

  if (myrank == 0) then
  
    write(*,*) 
    write(fid,*) 
    write(*,'(a)') trim(header)
    write(fid,'(a)') trim(header)
    write(*,'("     atol:",1pe12.4)') solver%newton_atol
    write(fid,'("     atol:",1pe12.4)') solver%newton_atol
    write(*,'("     rtol:",1pe12.4)') solver%newton_rtol
    write(fid,'("     rtol:",1pe12.4)') solver%newton_rtol
    write(*,'("     stol:",1pe12.4)') solver%newton_stol
    write(fid,'("     stol:",1pe12.4)') solver%newton_stol
    write(*,'("     dtol:",1pe12.4)') solver%newton_dtol
    write(fid,'("     dtol:",1pe12.4)') solver%newton_dtol
    write(*,'("inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(fid,'("inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(*,'("inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(fid,'("inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(*,'("    maxit:",i6)') solver%newton_maxit
    write(fid,'("    maxit:",i6)') solver%newton_maxit
    write(*,'("     maxf:",i6)') solver%newton_maxf
    write(fid,'("     maxf:",i6)') solver%newton_maxf
  
    if (solver%inexact_newton) then
      write(*,'("inexact newton: on")')
      write(fid,'("inexact newton: on")')
    else
      write(*,'("inexact newton: off")')
      write(fid,'("inexact newton: off")')
    endif
        
    if (solver%print_convergence) then
      write(*,'("print convergence: on")')
      write(fid,'("print convergence: on")')
    else
      write(*,'("print convergence: off")')
      write(fid,'("print convergence: off")')
    endif
        
    if (solver%print_detailed_convergence) then
      write(*,'("print detailed convergence: on")')
      write(fid,'("print detailed convergence: on")')
    else
      write(*,'("print detailed convergence: off")')
      write(fid,'("print detailed convergence: off")')
    endif
        
    if (solver%check_infinity_norm) then
      write(*,'("check infinity norm: on")')
      write(fid,'("check infinity norm: on")')
    else
      write(*,'("check infinity norm: off")')
      write(fid,'("check infinity norm: off")')
    endif
        
    if (solver%force_at_least_1_iteration) then
      write(*,'("force at least 1 iteration: on")')
      write(fid,'("force at least 1 iteration: on")')
    else
      write(*,'("force at least 1 iteration: off")')
      write(fid,'("force at least 1 iteration: off")')
    endif
  endif

end subroutine SolverPrintNewtonInfo

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

  if (.not.associated(solver)) return
    
  if (solver%J /= 0) call MatDestroy(solver%J,ierr)
  if (solver%interpolation /= 0) call MatDestroy(solver%interpolation,ierr)
  if (solver%matfdcoloring /= 0) call MatFDColoringDestroy(solver%matfdcoloring,ierr)
  if (solver%snes /= 0) call SNESDestroy(solver%snes,ierr)
  solver%ksp = 0
  solver%pc = 0
    
  deallocate(solver)
  nullify(solver)
  
end subroutine SolverDestroy
  
end module Solver_module
