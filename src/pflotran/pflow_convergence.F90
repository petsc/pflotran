module pflow_convergence_module

  implicit none

  private

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petsclog.h"

  public :: PFLOWConvergenceTest
  
contains

subroutine PFLOWConvergenceTest(snes_,it,xnorm,pnorm,fnorm,reason,ctx,ierr)

  implicit none
  
  SNES :: snes_
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: pnorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscInt :: ctx
  PetscErrorCode :: ierr
  
  Vec :: solution
  Vec :: update
  Vec :: residual
  PetscReal :: inorm_solution  !infinity norm
  PetscReal :: inorm_update  !infinity norm
  PetscReal :: inorm_residual  !infinity norm
  
!typedef enum {/* converged */
!              SNES_CONVERGED_FNORM_ABS         =  2, /* F < F_minabs */
!              SNES_CONVERGED_FNORM_RELATIVE    =  3, /* F < F_mintol*F_initial */
!              SNES_CONVERGED_PNORM_RELATIVE    =  4, /* step size small */
!              SNES_CONVERGED_ITS               =  5, /* maximum iterations reached */
!              SNES_CONVERGED_TR_DELTA          =  7,
!              /* diverged */
!              SNES_DIVERGED_FUNCTION_DOMAIN    = -1,  
!              SNES_DIVERGED_FUNCTION_COUNT     = -2,  
!              SNES_DIVERGED_LINEAR_SOLVE       = -3, 
!              SNES_DIVERGED_FNORM_NAN          = -4, 
!              SNES_DIVERGED_MAX_IT             = -5,
!              SNES_DIVERGED_LS_FAILURE         = -6,
!              SNES_DIVERGED_LOCAL_MIN          = -8,  /* || J^T b || is small, implies converged to local minimum of F() */
!              SNES_CONVERGED_ITERATING         =  0} SNESConvergedReason;

! first call the default convergence test
  call SNESDefaultConverged(snes_,it,xnorm,pnorm,fnorm,reason,ctx,ierr)
 
! if you are not happy, apply some other criteria

!  call SNESGetSolution(snes_,solution)
!  call SNESGetFunction(snes_,residual,PETSC_NULL_INTEGER,PETSC_NULL_OBJECT, &
!                       ierr)
!  call SNESGetSolutionUpdate(snes_,update)
!  
!  ! infinity norms
!  call VecNorm(residual,NORM_INFINITY,inorm_solution,ierr)
!  call VecNorm(residual,NORM_INFINITY,inorm_update,ierr)
!  call VecNorm(residual,NORM_INFINITY,inorm_residual,ierr)


end subroutine PFLOWConvergenceTest

end module pflow_convergence_module
