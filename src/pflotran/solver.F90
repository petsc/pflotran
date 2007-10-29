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
  end type solver_type
  
  public :: createSolver, &
            ComputeMFJacobian, &
            MonitorH
  
contains

! ************************************************************************** !
!
! createSolver: Allocates and initializes a new Solver object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function createSolver()

  implicit none
  
  type(solver_type), pointer :: createSolver
  
  type(solver_type), pointer :: solver
  
  allocate(solver)
  
  ! initialize to default values
  solver%atol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%rtol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%stol = PETSC_DEFAULT_DOUBLE_PRECISION
  solver%maxit = PETSC_DEFAULT_INTEGER
  solver%maxf = PETSC_DEFAULT_INTEGER
  solver%idt_switch = 0
  
  createSolver => solver
  
end function createSolver
  
! ************************************************************************** !
!
! ComputeMFJacobian: Sets Jacobian B = J
! author:
! date: 
!
! ************************************************************************** !
subroutine ComputeMFJacobian(snes, x, J, B, flag, ctx, ierr)
  
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
  
end subroutine ComputeMFJacobian

! ************************************************************************** !
!
! MonitorH: Sets Jacobian B = J
! author:
! date: 
!
! ************************************************************************** !
subroutine MonitorH(snes, its, norm, option)
  
  use Option_module
  
  implicit none

  SNES, intent(in) :: snes
  integer, intent(in) :: its
  PetscReal, intent(in) :: norm
  type(option_type) :: option
  
  integer :: ierr
  integer :: myrank
  PetscScalar :: h
  
  call MatMFFDGetH(option%J, h, ierr)

  if (option%myrank == 0) then
    write(*,*) "#At SNES iteration ", its, "h is ", h
  endif

end subroutine MonitorH
  
end module Solver_module
