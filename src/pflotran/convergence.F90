module Convergence_module

  use Solver_module
  use Option_module

  implicit none

  private

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petsclog.h"

  type, public :: convergence_context_type
    type(solver_type), pointer :: solver
    type(option_type), pointer :: option
  end type convergence_context_type


  public :: ConvergenceContextCreate, ConvergenceTest, &
            ConvergenceContextDestroy
  
contains

! ************************************************************************** !
!
! ConvergenceContextCreate: Creates a context containing pointer
!                           for convergence subroutines
! author: Glenn Hammond
! date: 02/12/08
!
! ************************************************************************** !
function ConvergenceContextCreate(solver,option)

  implicit none
  
  type(convergence_context_type), pointer :: ConvergenceContextCreate
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  
  type(convergence_context_type), pointer :: context
  
  allocate(context)
  context%solver => solver
  context%option => option

  ConvergenceContextCreate => context

end function ConvergenceContextCreate

! ************************************************************************** !
!
! ConvergenceTest: User defined convergence test
! author: Glenn Hammond
! date: 02/12/08
!
! ************************************************************************** !
subroutine ConvergenceTest(snes_,it,xnorm,pnorm,fnorm,reason,context,ierr)

  implicit none
  
  SNES :: snes_
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: pnorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  type(convergence_context_type) :: context
  PetscErrorCode :: ierr
  
  PetscInt :: ctx = 0

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  
  Vec :: solution_vec
  Vec :: update_vec
  Vec :: residual_vec
  PetscReal :: inorm_solution  !infinity norms
  PetscReal :: inorm_update  
  PetscReal :: inorm_residual  
  
  PetscInt :: i, ndof
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

  KSP :: ksp
  
  PetscInt, allocatable :: imax_solution(:)
  PetscInt, allocatable :: imax_update(:)
  PetscInt, allocatable :: imax_residual(:)
  PetscReal, allocatable :: max_solution_val(:)
  PetscReal, allocatable :: max_update_val(:)
  PetscReal, allocatable :: max_residual_val(:)
  
  PetscInt, allocatable :: imin_solution(:)
  PetscInt, allocatable :: imin_update(:)
  PetscInt, allocatable :: imin_residual(:)
  PetscReal, allocatable :: min_solution_val(:)
  PetscReal, allocatable :: min_update_val(:)
  PetscReal, allocatable :: min_residual_val(:)
  PetscReal, pointer :: vec_ptr(:)
  
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  character(len=MAXWORDLENGTH) :: word
  logical :: print_sol_norm_info = .false.
  logical :: print_upd_norm_info = .false.
  logical :: print_res_norm_info = .false.
  logical :: print_norm_by_dof_info = .false.
  logical :: print_max_val_and_loc_info = .false.
  logical :: print_1_norm_info = .false.
  logical :: print_2_norm_info = .false.
  logical :: print_inf_norm_info = .false.

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

  solver => context%solver
  option => context%option

  if (option%use_touch_options) then
    word = 'detailed_convergence'
    if (OptionCheckTouch(option,word)) then
      if (solver%print_detailed_convergence) then
        solver%print_detailed_convergence = PETSC_FALSE
      else
        solver%print_detailed_convergence = PETSC_TRUE
      endif
    endif
  endif
  
  ! always take one iteration
  if (solver%force_at_least_1_iteration) then
!    call SNESGetIterationNumber(snes_,it,ierr)
    if (it == 0) then
      reason = 0
      return
    endif
  endif

  call SNESDefaultConverged(snes_,it,xnorm,pnorm,fnorm,reason, &
                            PETSC_NULL_OBJECT,ierr)

!  if (reason <= 0 .and. solver%check_infinity_norm) then
  if (solver%check_infinity_norm) then
  
    call SNESGetFunction(snes_,residual_vec,PETSC_NULL_OBJECT, &
                         PETSC_NULL_INTEGER,ierr)

    call VecNorm(residual_vec,NORM_INFINITY,inorm_residual,ierr)

    call SNESGetSolutionUpdate(snes_,update_vec,ierr)
    call VecNorm(update_vec,NORM_INFINITY,inorm_update,ierr)

    if (inorm_residual < solver%newton_inf_res_tol) then
!      if (option%myrank == 0) print *, 'converged from infinity', inorm_residual
      reason = 10
    else
!      if (reason > 0 .and. inorm_residual > 100.d0*solver%newton_inf_res_tol) &
!        reason = 0
    endif

    if (inorm_update < solver%newton_inf_upd_tol .and. it > 0) then
!      if (option%myrank == 0) print *, 'converged from infinity', inorm_residual
      reason = 11
    endif

    if (option%myrank == 0 .and. solver%print_convergence) then
      i = int(reason)
      select case(i)
        case(2)
          string = 'atol'
        case(3)
          string = 'rtol'
        case(4)
          string = 'stol'
        case(10)
          string = 'itol_res'
        case(11)
          string = 'itol_upd'
        case default
          write(string,'(i3)') reason
      end select
      write(*,'(i3," fnrm:",es9.2, &
              & " pnrm:",es9.2, &
              & " inrmr:",es9.2, &
              & " inrmu:",es9.2, &
              & " rsn: ",a)') it, fnorm, pnorm, inorm_residual, inorm_update, &
                              trim(string)
    endif
  else
    if (option%myrank == 0 .and. solver%print_convergence) then
      i = int(reason)
      select case(i)
        case(2)
          string = 'atol'
        case(3)
          string = 'rtol'
        case(4)
          string = 'stol'
        case(10)
          string = 'itol_res'
        case(11)
          string = 'itol_upd'
        case default
          write(string,'(i3)') reason
      end select
      write(*,'(i3," fnrm:",es10.2, &
              & " pnrm:",es10.2, &
              & 32x, &
              & " rsn: ",a)') it, fnorm, pnorm, trim(string)
    endif
  endif    

  if (solver%print_detailed_convergence) then

    call SNESGetSolution(snes_,solution_vec,ierr)
    ! the ctx object should really be PETSC_NULL_OBJECT.  A bug in petsc
    call SNESGetFunction(snes_,residual_vec,PETSC_NULL_OBJECT, &
                         PETSC_NULL_INTEGER, &
                         ierr)
    call SNESGetSolutionUpdate(snes_,update_vec,ierr)
    
    ! infinity norms
    call VecNorm(solution_vec,NORM_INFINITY,inorm_solution,ierr)
    call VecNorm(update_vec,NORM_INFINITY,inorm_update,ierr)
    call VecNorm(residual_vec,NORM_INFINITY,inorm_residual,ierr)

    call VecNorm(solution_vec,NORM_1,norm1_solution,ierr)
    call VecNorm(update_vec,NORM_1,norm1_update,ierr)
    call VecNorm(residual_vec,NORM_1,norm1_residual,ierr)
    
    call VecGetBlockSize(solution_vec,ndof,ierr)
    
    allocate(fnorm_solution_stride(ndof))
    allocate(fnorm_update_stride(ndof))
    allocate(fnorm_residual_stride(ndof))
    allocate(inorm_solution_stride(ndof))
    allocate(inorm_update_stride(ndof))
    allocate(inorm_residual_stride(ndof))
    allocate(norm1_solution_stride(ndof))
    allocate(norm1_update_stride(ndof))
    allocate(norm1_residual_stride(ndof))
    
    allocate(imax_solution(ndof))
    allocate(imax_update(ndof))
    allocate(imax_residual(ndof))
    allocate(max_solution_val(ndof))
    allocate(max_update_val(ndof))
    allocate(max_residual_val(ndof))

    allocate(imin_solution(ndof))
    allocate(imin_update(ndof))
    allocate(imin_residual(ndof))
    allocate(min_solution_val(ndof))
    allocate(min_update_val(ndof))
    allocate(min_residual_val(ndof))

    call VecStrideNormAll(solution_vec,NORM_1,norm1_solution_stride,ierr)
    call VecStrideNormAll(update_vec,NORM_1,norm1_update_stride,ierr)
    call VecStrideNormAll(residual_vec,NORM_1,norm1_residual_stride,ierr)
    call VecStrideNormAll(solution_vec,NORM_2,fnorm_solution_stride,ierr)
    call VecStrideNormAll(update_vec,NORM_2,fnorm_update_stride,ierr)
    call VecStrideNormAll(residual_vec,NORM_2,fnorm_residual_stride,ierr)
    call VecStrideNormAll(solution_vec,NORM_INFINITY,inorm_solution_stride,ierr)
    call VecStrideNormAll(update_vec,NORM_INFINITY,inorm_update_stride,ierr)
    call VecStrideNormAll(residual_vec,NORM_INFINITY,inorm_residual_stride,ierr)
    
    ! can't use VecStrideMaxAll since the index location is not currently supported.
    do i=1,ndof
      call VecStrideMax(solution_vec,i-1,imax_solution(i),max_solution_val(i),ierr)
      call VecStrideMax(update_vec,i-1,imax_update(i),max_update_val(i),ierr)
      call VecStrideMax(residual_vec,i-1,imax_residual(i),max_residual_val(i),ierr)
      ! tweak the index to get the cell id from the mdof vector
      imax_solution(i) = imax_solution(i)/ndof
      imax_update(i) = imax_update(i)/ndof
      imax_residual(i) = imax_residual(i)/ndof
    enddo

    do i=1,ndof
      call VecStrideMin(solution_vec,i-1,imin_solution(i),min_solution_val(i),ierr)
      call VecStrideMin(update_vec,i-1,imin_update(i),min_update_val(i),ierr)
      call VecStrideMin(residual_vec,i-1,imin_residual(i),min_residual_val(i),ierr)
      ! tweak the index to get the cell id from the mdof vector
      imin_solution(i) = imin_solution(i)/ndof
      imin_update(i) = imin_update(i)/ndof
      imin_residual(i) = imin_residual(i)/ndof
    enddo

    if (option%myrank == 0) then
      select case(reason)
        case (10)
          string = "CONVERGED_USER_NORM_INF_REL"
        case (11)
          string = "CONVERGED_USER_NORM_INF_UPD"
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
      
      print_sol_norm_info = .true.  ! solution_vec norm information
      print_upd_norm_info = .true.  ! update_vec norm information
      print_res_norm_info = .true.  ! residual_vec norm information
    
      !print_norm_by_dof_info = .true.
      print_max_val_and_loc_info = .true.

      !print_1_norm_info = .true.
      print_2_norm_info = .true.
      print_inf_norm_info = .true.

      print *
      print *, 'reason: ', reason, ' - ', trim(string)
      print *, 'SNES iteration :', it
      call SNESGetKSP(snes_,ksp,ierr)
      call KSPGetIterationNumber(ksp,i,ierr)
      print *, 'KSP iterations :', i
      if (print_1_norm_info) then
        if (print_sol_norm_info) print *, 'norm_1_solution:   ', norm1_solution
        if (print_upd_norm_info) print *, 'norm_1_update:     ', norm1_update
        if (print_res_norm_info) print *, 'norm_1_residual:   ', norm1_residual
      endif
      if (print_2_norm_info) then
        if (print_sol_norm_info) print *, 'norm_2_solution:   ', fnorm_solution_stride
        if (print_upd_norm_info) print *, 'norm_2_update:     ', fnorm_update_stride
        if (print_res_norm_info) print *, 'norm_2_residual:   ', fnorm_residual_stride
      endif
      if (print_inf_norm_info) then
        if (print_sol_norm_info) print *, 'norm_inf_solution: ', inorm_solution
        if (print_upd_norm_info) print *, 'norm_inf_update:   ', inorm_update
        if (print_res_norm_info) print *, 'norm_inf_residual: ', inorm_residual
      endif
      if (print_max_val_and_loc_info) then
        print *, 'max/min locations (zero-based index) by dof:'
        do i=1,ndof
          print *, '  dof: ', i
          if (print_sol_norm_info) then
            print *, '    solution_vec max: ', imax_solution(i), max_solution_val(i)
            print *, '    solution_vec min: ', imin_solution(i), min_solution_val(i)
          endif
          if (print_upd_norm_info) then ! since update is -dx, need to invert
            print *, '    update_vec max:   ', imin_update(i), -1.d0*min_update_val(i)
            print *, '    update_vec min:   ', imax_update(i), -1.d0*max_update_val(i)
          endif
          if (print_res_norm_info) then
            print *, '    residual_vec max: ', imax_residual(i), max_residual_val(i)
            print *, '    residual_vec min: ', imin_residual(i), min_residual_val(i)
          endif
        enddo
      endif
      if (print_norm_by_dof_info) then
        print *, 'norm by dof:'
        do i=1,ndof
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
    
    deallocate(imax_solution)
    deallocate(imax_update)
    deallocate(imax_residual)
    deallocate(max_solution_val)
    deallocate(max_update_val)
    deallocate(max_residual_val)

    deallocate(imin_solution)
    deallocate(imin_update)
    deallocate(imin_residual)
    deallocate(min_solution_val)
    deallocate(min_update_val)
    deallocate(min_residual_val)
    
  endif
  
end subroutine ConvergenceTest

! ************************************************************************** !
!
! ConvergenceContextDestroy: Destroy context
! author: Glenn Hammond
! date: 02/12/08
!
! ************************************************************************** !
subroutine ConvergenceContextDestroy(context)

  implicit none
  
  type(convergence_context_type), pointer :: context
  
  if (associated(context)) deallocate(context)
  nullify(context)

end subroutine ConvergenceContextDestroy

end module Convergence_module
