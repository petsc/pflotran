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
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

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
  
end module Solver_module
