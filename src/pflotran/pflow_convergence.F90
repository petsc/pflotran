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

subroutine PFLOWConvergenceTest(snes_,it,xnorm,pnorm,fnorm,reason,option,ierr)

  use Option_module

  implicit none
  
  SNES :: snes_
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: pnorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  type(option_type) :: option
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
  logical :: print_sol_norm_info = .false.
  logical :: print_upd_norm_info = .false.
  logical :: print_res_norm_info = .false.
  logical :: print_norm_by_dof_info = .false.
  logical :: print_max_val_and_loc_info = .false.
  logical :: print_1_norm_info = .false.
  logical :: print_2_norm_info = .false.
  logical :: print_inf_norm_info = .false.
  
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
!              SNES_DIVERGED_LOCAL_MIN          = -8,  /* || J^T b || is small,
!                                        implies converged to local minimum of F() */
!              SNES_CONVERGED_ITERATING         =  0} SNESConvergedReason;

#ifdef CHUAN
  ! always take one iteration
  call SNESGetIterationNumber(option%snes,it,ierr)
  if (it == 0) then
    reason = 0
    return
  endif
#endif

  call SNESDefaultConverged(snes_,it,xnorm,pnorm,fnorm,reason,PETSC_NULL_OBJECT,ierr)
 
#ifdef CHUAN
  if (reason <= 0) then
  
    call SNESGetFunction(snes_,residual,PETSC_NULL_OBJECT,PETSC_NULL_INTEGER, &
                         ierr)

    call VecNorm(residual,NORM_INFINITY,inorm_residual,ierr)
  
    if (inorm_residual < option%inf_tol) then
      if (option%myrank == 0) print *, 'converged from infinity', inorm_residual
      reason = 1
    endif

  endif    
 
  if (option%myrank == 0) print *, 'snes_default', xnorm,pnorm,fnorm,reason

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
  
  allocate(fnorm_solution_stride(option%ndof))
  allocate(fnorm_update_stride(option%ndof))
  allocate(fnorm_residual_stride(option%ndof))
  allocate(inorm_solution_stride(option%ndof))
  allocate(inorm_update_stride(option%ndof))
  allocate(inorm_residual_stride(option%ndof))
  allocate(norm1_solution_stride(option%ndof))
  allocate(norm1_update_stride(option%ndof))
  allocate(norm1_residual_stride(option%ndof))
  
  allocate(imax_solution(option%ndof))
  allocate(imax_update(option%ndof))
  allocate(imax_residual(option%ndof))
  allocate(max_solution_val(option%ndof))
  allocate(max_update_val(option%ndof))
  allocate(max_residual_val(option%ndof))

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
  do i=1,option%ndof
    call VecStrideMax(solution,i-1,imax_solution(i),max_solution_val(i),ierr)
    call VecStrideMax(update,i-1,imax_update(i),max_update_val(i),ierr)
    call VecStrideMax(residual,i-1,imax_residual(i),max_residual_val(i),ierr)
    ! tweak the index to get the cell id from the mdof vector
    imax_solution(i) = imax_solution(i)/option%ndof
    imax_update(i) = imax_update(i)/option%ndof
    imax_residual(i) = imax_residual(i)/option%ndof
  enddo

  if (option%myrank == 0) then
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

    ! uncomment the lines below to determine data printed
    
    !print_sol_norm_info = .true.  ! solution norm information
    !print_upd_norm_info = .true.  ! update norm information
    print_res_norm_info = .true.  ! residual norm information
  
    print_norm_by_dof_info = .true.
    print_max_val_and_loc_info = .true.

    print_1_norm_info = .true.
    print_2_norm_info = .true.
    print_inf_norm_info = .true.

    print *
    print *, 'reason: ', reason, ' - ', trim(string)
    print *, 'its :', it
    if (print_1_norm_info) then
      if (print_sol_norm_info) print *, 'norm_1_solution:   ', norm1_solution
      if (print_upd_norm_info) print *, 'norm_1_update:     ', norm1_update
      if (print_res_norm_info) print *, 'norm_1_residual:   ', norm1_residual
    endif
    if (print_2_norm_info) then
      if (print_sol_norm_info) print *, 'norm_2_solution:   ', xnorm
      if (print_upd_norm_info) print *, 'norm_2_update:     ', pnorm
      if (print_res_norm_info) print *, 'norm_2_residual:   ', fnorm
    endif
    if (print_inf_norm_info) then
      if (print_sol_norm_info) print *, 'norm_inf_solution: ', inorm_solution
      if (print_upd_norm_info) print *, 'norm_inf_update:   ', inorm_update
      if (print_res_norm_info) print *, 'norm_inf_residual: ', inorm_residual
    endif
    if (print_max_val_and_loc_info) then
      print *, 'max locations by dof:'
      do i=1,option%ndof
        print *, '  dof: ', i
        if (print_sol_norm_info) &
          print *, '    solution max: ', imax_solution(i), max_solution_val(i)
        if (print_upd_norm_info) &
          print *, '    update max:   ', imax_update(i), max_update_val(i)
        if (print_res_norm_info) &
          print *, '    residual max: ', imax_residual(i), max_residual_val(i)
      enddo
    endif
    if (print_norm_by_dof_info) then
      print *, 'norm by dof:'
      do i=1,option%ndof
        print *, '  dof: ', i
        if (print_sol_norm_info) then
          if (print_1_norm_info) &
            print *, '    norm_1_solution:   ', norm1_solution_stride(i)
          if (print_2_norm_info) &
            print *, '    norm_2_solution:   ', fnorm_solution_stride(i)
          if (print_inf_norm_info) &
            print *, '    norm_inf_solution: ', inorm_solution_stride(i)
          if (print_1_norm_info .or. print_2_norm_info .or. &
              print_inf_norm_info) print *, '    -'
        endif
        if (print_upd_norm_info) then
          if (print_1_norm_info) &
            print *, '    norm_1_update:   ', norm1_update_stride(i)
          if (print_2_norm_info) &
            print *, '    norm_2_update:   ', fnorm_update_stride(i)
          if (print_inf_norm_info) &
            print *, '    norm_inf_update: ', inorm_update_stride(i)
          if (print_1_norm_info .or. print_2_norm_info .or. &
              print_inf_norm_info) print *, '    -'
        endif
        if (print_res_norm_info) then
          if (print_1_norm_info) &
            print *, '    norm_1_residual:   ', norm1_residual_stride(i)
          if (print_2_norm_info) &
            print *, '    norm_2_residual:   ', fnorm_residual_stride(i)
          if (print_inf_norm_info) &
            print *, '    norm_inf_residual: ', inorm_residual_stride(i)
        endif
      enddo
    endif
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
