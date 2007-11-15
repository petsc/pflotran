module Solver_module
 
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
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsnes.h"

  type, public :: solver_type
    real*8 :: atol       ! absolute tolerance
    real*8 :: rtol       ! relative tolerance
    real*8 :: stol       ! relative tolerance (relative to previous iteration)
    real*8 :: dtol       ! divergence tolerance
    real*8 :: inf_tol    ! infinity tolerance
    integer :: maxit     ! maximum number of iterations
    integer :: maxf      ! maximum number of function evaluations
    integer :: idt_switch
    
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
  
  public :: SolverCreate, &
            SolverDestroy, &
            SolverComputeMFJacobian, &
            SolverMonitorH
  
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
  solver%idt_switch = 0
  
  solver%J = 0
  solver%matfdcoloring = 0
  solver%snes = 0
  solver%ksp_type = ""
  solver%pc_type = ""
  solver%ksp = 0
  solver%pc = 0
  
  SolverCreate => solver
  
end function SolverCreate
  
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
  integer, intent(inout) :: ctx(*)
  integer, intent(out) :: ierr

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
  integer, intent(in) :: its
  PetscReal, intent(in) :: norm
  type(solver_type) :: solver
  type(option_type) :: option
  
  integer :: ierr
  integer :: myrank
  PetscScalar :: h
  
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
