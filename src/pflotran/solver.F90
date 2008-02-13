module Solver_module
 
  implicit none

  private
 
#include "definitions.h"

#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsnes.h"

  type, public :: solver_type
    PetscReal :: atol       ! absolute tolerance
    PetscReal :: rtol       ! relative tolerance
    PetscReal :: stol       ! relative tolerance (relative to previous iteration)
    PetscReal :: dtol       ! divergence tolerance
    PetscReal :: inf_tol    ! infinity tolerance
    PetscInt :: maxit     ! maximum number of iterations
    PetscInt :: maxf      ! maximum number of function evaluations
    
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
    
    PetscTruth :: inexact_newton
            
  end type solver_type
  
  interface SolverRead
    module procedure SolverReadPflow
  end interface SolverRead

  public :: SolverCreate, &
            SolverDestroy, &
            SolverRead, &
            SolverCreateSNES, &
            SolverSetSNESOptions
  
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
  solver%atol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%stol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%maxit = PETSC_DEFAULT_INTEGER
  solver%maxf = PETSC_DEFAULT_INTEGER
  
  solver%J = 0
  solver%matfdcoloring = 0
  solver%snes = 0
  solver%ksp_type = KSPBCGS
  solver%pc_type = ""
  solver%ksp = 0
  solver%pc = 0
  
  solver%inexact_newton = PETSC_FALSE
  
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

  use Option_module

  implicit none
  
  type(solver_type) :: solver

  PetscErrorCode :: ierr
  
  call SNESCreate(PETSC_COMM_WORLD, solver%snes, ierr)

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
subroutine SolverSetSNESOptions(solver,option)

  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  ! needed for SNESLineSearchGetParams()/SNESLineSearchSetParams()
  PetscReal :: alpha, maxstep, steptol
  PetscErrorCode :: ierr
  
  ! if ksp_type or pc_type specified in input file, set them here
  if (len_trim(solver%ksp_type) > 1) &
    call KSPSetType(solver%ksp,solver%ksp_type,ierr)
  if (len_trim(solver%pc_type) > 1) &
    call PCSetType(solver%pc,solver%pc_type,ierr)

  call KSPSetTolerances(solver%ksp,solver%rtol,solver%atol,solver%dtol, &
                        10000,ierr)

  ! allow override from command line
  call KSPSetFromOptions(solver%ksp,ierr)
  call PCSetFromOptions(solver%pc,ierr)
    
  ! get the ksp_type and pc_type incase of command line override.
  call KSPGetType(solver%ksp,solver%ksp_type,ierr)
  call PCGetType(solver%pc,solver%pc_type,ierr)

  call printMsg(option,'Solver: '//trim(solver%ksp_type))
  call printMsg(option,'Preconditioner: '//trim(solver%pc_type))

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(solver%snes, solver%atol, solver%rtol, solver%stol, & 
                         solver%maxit, solver%maxf, ierr)

  ! set inexact newton, currently applies default settings
  if (solver%inexact_newton == PETSC_TRUE) &
    call SNESKSPSetUseEW(solver%snes,PETSC_TRUE,ierr)

  ! allow override from command line; for some reason must come before
  ! LineSearchParams, or they crash
  call SNESSetFromOptions(solver%snes, ierr) 
    
  call SNESLineSearchGetParams(solver%snes, alpha, maxstep, steptol, ierr)  
  call SNESLineSearchSetParams(solver%snes, alpha, maxstep, solver%stol, ierr)  


end subroutine SolverSetSNESOptions
  
! ************************************************************************** !
!
! SolverReadPflow: Reads debugging data from the input file
! author: Glenn Hammond
! date: 12/21/07
!
! ************************************************************************** !
subroutine SolverReadPflow(solver,fid,myrank)

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
    call fiErrorMsg(myrank,'keyword','SOLVER', ierr)
    call fiWordToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('KRYLOV_TYPE','KRYLOV','KSP','KSP_TYPE')
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
          case('ILU')
            solver%pc_type = PCILU
          case('LU')
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
        
    end select 
  
  enddo  

end subroutine SolverReadPflow

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
  if (solver%matfdcoloring /= 0) call MatFDColoringDestroy(solver%matfdcoloring,ierr)
  if (solver%snes /= 0) call SNESDestroy(solver%snes,ierr)
  solver%ksp = 0
  solver%pc = 0
    
  deallocate(solver)
  nullify(solver)
  
end subroutine SolverDestroy
  
end module Solver_module
