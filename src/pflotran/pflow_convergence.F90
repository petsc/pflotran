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

subroutine PFLOWConvergenceTest(snes_,it,xnorm,pnorm,fnorm,reason,grid,ierr)

  use pflow_gridtype_module

  implicit none
  
  SNES :: snes_
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: pnorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  type(pflowGrid) :: grid
  PetscErrorCode :: ierr
  
  PetscInt :: ctx = 0

  Vec :: solution
  Vec :: update
  Vec :: residual
  PetscReal :: inorm_solution  !infinity norm
  PetscReal :: inorm_update  !infinity norm
  PetscReal :: inorm_residual  !infinity norm
  
  character(len=128) :: string
  
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

!first call the default convergence test - always take one iteration
#ifdef CHUAN
  call SNESGetIterationNumber(grid%snes, it, ierr)
  if(it == 0) then
    reason = 0
    return
  endif
#endif

  call SNESDefaultConverged(snes_,it,xnorm,pnorm,fnorm,reason,PETSC_NULL_OBJECT,ierr)
 
! if you are not happy, apply some other criteria
#ifdef CHUAN
  call SNESGetFunction(snes_,residual,PETSC_NULL_OBJECT,PETSC_NULL_INTEGER, &
                       ierr)
  if(reason > 0) return
  
  call VecNorm(residual,NORM_INFINITY,inorm_residual,ierr)
  
  if(inorm_residual < grid%inf_tol) then
    if (grid%myrank == 0) print *, 'converged from infinity', inorm_residual
    reason = 1
  endif    
 
  if (grid%myrank == 0) print*, 'snes_default', xnorm,pnorm,fnorm,reason

#endif

#ifdef GLENN
  call SNESGetSolution(snes_,solution,ierr)
                              ! the ctx object should really be PETSC_NULL_OBJECT.  A bug in petsc
  call SNESGetFunction(snes_,residual,PETSC_NULL_OBJECT,PETSC_NULL_INTEGER, &
                       ierr)
  call SNESGetSolutionUpdate(snes_,update,ierr)
  
  ! infinity norms
  call VecNorm(solution,NORM_INFINITY,inorm_solution,ierr)
  call VecNorm(update,NORM_INFINITY,inorm_update,ierr)
  call VecNorm(residual,NORM_INFINITY,inorm_residual,ierr)


  if (grid%myrank == 0) then
    select case(reason)
      case(SNES_CONVERGED_FNORM_ABS)
        string = "SNES_CONVERGED_FNORM_ABS"
      case(SNES_CONVERGED_FNORM_RELATIVE)
        string = "SNES_CONVERGED_FNORM_RELATIVE"
      case(SNES_CONVERGED_PNORM_RELATIVE)
        string = "SNES_CONVERGED_PNORM_RELATIVE"
      case(SNES_CONVERGED_ITS)
        string = "SNES_CONVERGED_ITS"
      case(SNES_CONVERGED_TR_DELTA)
        string = "SNES_CONVERGED_TR_DELTA"
!      case(SNES_DIVERGED_FUNCTION_DOMAIN)
!        string = "SNES_DIVERGED_FUNCTION_DOMAIN"
      case(SNES_DIVERGED_FUNCTION_COUNT)
        string = "SNES_DIVERGED_FUNCTION_COUNT"
      case(SNES_DIVERGED_LINEAR_SOLVE)
        string = "SNES_DIVERGED_LINEAR_SOLVE"
      case(SNES_DIVERGED_FNORM_NAN)
        string = "SNES_DIVERGED_FNORM_NAN"
      case(SNES_DIVERGED_MAX_IT)
        string = "SNES_DIVERGED_MAX_IT"
      case(SNES_DIVERGED_LS_FAILURE)
        string = "SNES_DIVERGED_LS_FAILURE"
      case(SNES_DIVERGED_LOCAL_MIN)
        string = "SNES_DIVERGED_LOCAL_MIN"
      case(SNES_CONVERGED_ITERATING)
        string = "SNES_CONVERGED_ITERATING"
    end select

    print *, 'reason: ', reason, ' - ', trim(string)
    print *, 'its :', it
    print *, 'xnorm2: ', xnorm
    print *, 'pnorm2: ', pnorm
    print *, 'fnorm2: ', fnorm
    print *, 'inf_norm_solution: ', inorm_solution
    print *, 'inf_norm_update: ', inorm_update
    print *, 'inf_norm_residual: ', inorm_residual
    print *
  endif
#endif

end subroutine PFLOWConvergenceTest




end module pflow_convergence_module
