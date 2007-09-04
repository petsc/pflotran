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
  PetscReal :: inorm_solution  !infinity norms
  PetscReal :: inorm_update  
  PetscReal :: inorm_residual  
  
#ifdef GLENN
  integer :: i
  PetscReal, allocatable :: fnorm_solution_stride(:)
  PetscReal, allocatable :: fnorm_update_stride(:)
  PetscReal, allocatable :: fnorm_residual_stride(:)
  PetscReal, allocatable :: inorm_solution_stride(:)
  PetscReal, allocatable :: inorm_update_stride(:)
  PetscReal, allocatable :: inorm_residual_stride(:)
  
  PetscReal :: norm1_solution
  PetscReal :: norm1_update
  PetscReal :: norm1_residual
  PetscReal, allocatable :: norm1_solution_stride(:)
  PetscReal, allocatable :: norm1_update_stride(:)
  PetscReal, allocatable :: norm1_residual_stride(:)
  
  integer, allocatable :: imax_solution(:)
  integer, allocatable :: imax_update(:)
  integer, allocatable :: imax_residual(:)
  PetscReal, allocatable :: max_solution_val(:)
  PetscReal, allocatable :: max_update_val(:)
  PetscReal, allocatable :: max_residual_val(:)
  
  character(len=128) :: string
  
#endif

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

#ifdef CHUAN
  ! always take one iteration
  call SNESGetIterationNumber(grid%snes,it,ierr)
#endif

  call SNESDefaultConverged(snes_,it,xnorm,pnorm,fnorm,reason,PETSC_NULL_OBJECT,ierr)
 
#ifdef CHUAN
  if (reason <= 0 .and. it > 0) then
  
    call SNESGetFunction(snes_,residual,PETSC_NULL_OBJECT,PETSC_NULL_INTEGER, &
                         ierr)

    call VecNorm(residual,NORM_INFINITY,inorm_residual,ierr)
  
    if (inorm_residual < grid%inf_tol) then
      if (grid%myrank == 0) print *, 'converged from infinity', inorm_residual
      reason = 1
    endif

  else
    reason = 0
  endif    
 
  if (grid%myrank == 0) print *, 'snes_default', xnorm,pnorm,fnorm,reason

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

  call VecNorm(solution,NORM_1,norm1_solution,ierr)
  call VecNorm(update,NORM_1,norm1_update,ierr)
  call VecNorm(residual,NORM_1,norm1_residual,ierr)
  
  allocate(fnorm_solution_stride(grid%ndof))
  allocate(fnorm_update_stride(grid%ndof))
  allocate(fnorm_residual_stride(grid%ndof))
  allocate(inorm_solution_stride(grid%ndof))
  allocate(inorm_update_stride(grid%ndof))
  allocate(inorm_residual_stride(grid%ndof))
  allocate(norm1_solution_stride(grid%ndof))
  allocate(norm1_update_stride(grid%ndof))
  allocate(norm1_residual_stride(grid%ndof))
  
  allocate(imax_solution(grid%ndof))
  allocate(imax_update(grid%ndof))
  allocate(imax_residual(grid%ndof))
  allocate(max_solution_val(grid%ndof))
  allocate(max_update_val(grid%ndof))
  allocate(max_residual_val(grid%ndof))

  call VecStrideNormAll(solution,NORM_1,norm1_solution_stride,ierr)
  call VecStrideNormAll(update,NORM_1,norm1_update_stride,ierr)
  call VecStrideNormAll(residual,NORM_1,norm1_residual_stride,ierr)
  call VecStrideNormAll(solution,NORM_2,fnorm_solution_stride,ierr)
  call VecStrideNormAll(update,NORM_2,fnorm_update_stride,ierr)
  call VecStrideNormAll(residual,NORM_2,fnorm_residual_stride,ierr)
  call VecStrideNormAll(solution,NORM_INFINITY,inorm_solution_stride,ierr)
  call VecStrideNormAll(update,NORM_INFINITY,inorm_update_stride,ierr)
  call VecStrideNormAll(residual,NORM_INFINITY,inorm_residual_stride,ierr)
  
  ! can't use VecStrideMaxAll since the index location is not currently supported.
  do i=1,grid%ndof
    call VecStrideMax(solution,i-1,imax_solution(i),max_solution_val(i),ierr)
    call VecStrideMax(update,i-1,imax_update(i),max_update_val(i),ierr)
    call VecStrideMax(residual,i-1,imax_residual(i),max_residual_val(i),ierr)
    ! tweak the index to get the cell id from the mdof vector
    imax_solution(i) = imax_solution(i)/grid%ndof
    imax_update(i) = imax_update(i)/grid%ndof
    imax_residual(i) = imax_residual(i)/grid%ndof
  enddo

  if (grid%myrank == 0) then
    select case(reason)
      case (1)
        string = "CONVERGED_USER_NORM_INF"
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
      case default
        string = "UNKNOWN"
    end select

    print *, 'reason: ', reason, ' - ', trim(string)
    print *, 'its :', it
    print *, 'norm1_solution: ', norm1_solution
    print *, 'norm1_update: ', norm1_update
    print *, 'norm1_residual: ', norm1_residual
    print *, 'norm2_solution: ', xnorm
    print *, 'norm2_update: ', pnorm
    print *, 'norm2_residual: ', fnorm
    print *, 'inf_norm_solution: ', inorm_solution
    print *, 'inf_norm_update: ', inorm_update
    print *, 'inf_norm_residual: ', inorm_residual
    print *, 'max locations by dof:'
    do i=1,grid%ndof
      print *, '  dof: ', i
      print *, '    solution max: ', imax_solution(i), max_solution_val(i)
      print *, '    update max: ', imax_update(i), max_update_val(i)
      print *, '    residual max: ', imax_residual(i), max_residual_val(i)
    enddo
    print *, 'norm by dof:'
    do i=1,grid%ndof
      print *, '  dof: ', i
      print *, '    norm1_solution_stride: ', norm1_solution_stride(i)
      print *, '    norm2_solution_stride: ', fnorm_solution_stride(i)
      print *, '    inf_norm_solution_stride: ', inorm_solution_stride(i)
      print *, '    -'
      print *, '    norm1_update_stride: ', norm1_update_stride(i)
      print *, '    norm2_update_stride: ', fnorm_update_stride(i)
      print *, '    inf_norm_update_stride: ', inorm_update_stride(i)
      print *, '    -'
      print *, '    norm1_residual_stride: ', norm1_residual_stride(i)
      print *, '    norm2_residual_stride: ', fnorm_residual_stride(i)
      print *, '    inf_norm_residual_stride: ', inorm_residual_stride(i)
    enddo
    print *
  endif
  
  deallocate(fnorm_solution_stride)
  deallocate(fnorm_update_stride)
  deallocate(fnorm_residual_stride)
  deallocate(inorm_solution_stride)
  deallocate(inorm_update_stride)
  deallocate(inorm_residual_stride)
  deallocate(norm1_solution_stride)
  deallocate(norm1_update_stride)
  deallocate(norm1_residual_stride)
  
#endif

end subroutine PFLOWConvergenceTest




end module pflow_convergence_module
