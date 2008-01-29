module Solver_module
 
  implicit none

  private
 
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
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
        
  end type solver_type
  
  interface SolverRead
    module procedure SolverReadPflow
  end interface SolverRead

  public :: SolverCreate, &
            SolverDestroy, &
            SolverComputeMFJacobian, &
            SolverMonitorH, &
            SolverRead
  
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
  
  SolverCreate => solver
  
end function SolverCreate
  
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
! SolverComputeMFJacobian: Sets Jacobian B = J
! author:
! date: 
!
! ************************************************************************** !
subroutine SolverComputeMFJacobian(snes, x, J, B, flag, ctx, ierr)
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: x
  Mat, intent(out) :: J, B
  MatStructure, intent(in) :: flag
  PetscInt, intent(inout) :: ctx(*)
  PetscErrorCode, intent(out) :: ierr

  call MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY, ierr)
  B = J
  
end subroutine SolverComputeMFJacobian

! ************************************************************************** !
!
! SolverMonitorH: Sets Jacobian B = J
! author:
! date: 
!
! ************************************************************************** !
subroutine SolverMonitorH(snes, its, norm, solver, option)
  
  use Option_module
  
  implicit none

  SNES, intent(in) :: snes
  PetscInt, intent(in) :: its
  PetscReal, intent(in) :: norm
  type(solver_type) :: solver
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  PetscMPIInt :: myrank
  PetscReal :: h
  
  call MatMFFDGetH(solver%J, h, ierr)

  if (option%myrank == 0) then
    write(*,*) "#At SNES iteration ", its, "h is ", h
  endif

end subroutine SolverMonitorH

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
