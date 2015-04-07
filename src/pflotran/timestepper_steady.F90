module Timestepper_Steady_class

  use Timestepper_Base_class
  use Convergence_module
  use Solver_module
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
 
  type, public, extends(timestepper_base_type) :: timestepper_steady_type
  
    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations

    type(solver_type), pointer :: solver
    type(convergence_context_type), pointer :: convergence_context
  
  contains

    procedure, public :: Init => TimestepperSteadyInit
    procedure, public :: StepDT => TimestepperSteadyStepDT
    procedure, public :: UpdateDT => TimestepperSteadyUpdateDT

  end type timestepper_steady_type

  public :: TimestepperSteadyCreate

contains

! ************************************************************************** !

function TimestepperSteadyCreate()
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreate

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  stepper%solver => SolverCreate()

  TimestepperSteadyCreate => stepper

end function TimestepperSteadyCreate

! ************************************************************************** !

subroutine TimestepperSteadyInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(timestepper_steady_type) :: this

  call TimestepperBaseInit(this)

end subroutine TimestepperSteadyInit


! ************************************************************************** !

subroutine TimestepperSteadyUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  ! 

  use PM_Base_class
  use Option_module
  
  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  

end subroutine TimestepperSteadyUpdateDT

! ************************************************************************** !

subroutine TimestepperSteadyStepDT(this, process_model, stop_flag)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use PM_Base_class
  use Option_module
  use PM_Geomechanics_Force_class
  use Output_module, only : Output

  use Solver_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscviewer.h"
#include "finclude/petscsnes.h"

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  PetscBool :: failure

  PetscErrorCode :: ierr
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscInt :: snes_reason
  PetscInt :: icut
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm
  Vec :: residual_vec

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver

  solver => this%solver
  option => process_model%option

  sum_newton_iterations = 0
  sum_linear_iterations = 0
  icut = 0

  call process_model%InitializeTimestep()

  call process_model%PreSolve()

  call SNESSolve(solver%snes, PETSC_NULL_OBJECT, &
                 process_model%solution_vec, ierr);CHKERRQ(ierr)
     
  call SNESGetIterationNumber(solver%snes, num_newton_iterations,  &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes, num_linear_iterations,  &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in steady-solve, reason: ', &
                snes_reason
    endif
    failure = PETSC_TRUE
    return
  endif

  sum_newton_iterations = sum_newton_iterations + num_newton_iterations
  sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
  this%steps = this%steps + 1
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_newton_iterations = num_newton_iterations
  this%num_linear_iterations = num_linear_iterations  

  ! print screen output
  call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_OBJECT, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
    select type(pm => process_model)
      class is(pm_geomech_force_type)
        if (associated(pm%geomech_realization%geomech_discretization%grid)) then
          scaled_fnorm = fnorm/pm%geomech_realization%geomech_discretization% &
                          grid%nmax_node
        else
          scaled_fnorm = fnorm
        endif
    end select
    write(*,*) ''
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  if (option%print_screen_flag) print *, ""
  
  if (option%print_file_flag) then
    write(option%fid_out, '(" STEADY-SOLVE ",i6," snes_conv_reason: ",i4,/, &
      &"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]")') &
      this%steps, &
      snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations
  endif  

  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""

end subroutine TimestepperSteadyStepDT

end module Timestepper_Steady_class
