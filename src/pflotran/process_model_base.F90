module Process_Model_Base_class

  use Option_module
  use Output_Aux_module
  use Realization_Base_Class

  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

#ifdef ABSTRACT
  type, abstract, public :: process_model_base_type
#else
  type, public :: process_model_base_type
#endif
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    Vec :: solution_vec
    Vec :: residual_vec
    class(realization_base_type), pointer :: realization_base
    class(process_model_base_type), pointer :: next
  contains
#ifdef ABSTRACT  
    procedure(PMBaseInit), public, deferred :: Init
    procedure(PMBaseThisOnly), public, deferred :: InitializeRun
    procedure(PMBaseThisOnly), public, deferred :: FinalizeRun
    procedure(PMBaseResidual), public, deferred :: Residual
    procedure(PMBaseJacobian), public, deferred :: Jacobian
    procedure(PMBaseUpdateTimestep), public, deferred :: UpdateTimestep
    procedure(PMBaseThisOnly), public, deferred :: InitializeTimestep
    procedure(PMBaseThisOnly), public, deferred :: PreSolve
    procedure(PMBaseThisOnly), public, deferred :: PostSolve
    procedure(PMBaseThisOnly), public, deferred :: FinalizeTimestep
    procedure(PMBaseFunctionThisOnly), public, deferred :: AcceptSolution
    procedure(PMBaseCheckUpdatePre), public, deferred :: CheckUpdatePre
    procedure(PMBaseCheckUpdatePost), public, deferred :: CheckUpdatePost
    procedure(PMBaseThisOnly), public, deferred :: TimeCut
    procedure(PMBaseThisOnly), public, deferred :: UpdateSolution
    procedure(PMBaseThisOnly), public, deferred :: MaxChange
    procedure(PMBaseComputeMassBalance), public, deferred :: ComputeMassBalance
    procedure(PMBaseThisOnly), public, deferred :: Destroy
#else
    procedure, public :: Init => PMBaseInit
    procedure, public :: InitializeRun => PMBaseThisOnly
    procedure, public :: FinalizeRun => PMBaseThisOnly
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: UpdateTimestep => PMBaseUpdateTimestep
    procedure, public :: InitializeTimestep => PMBaseThisOnly
    procedure, public :: PreSolve => PMBaseThisOnly
    procedure, public :: PostSolve => PMBaseThisOnly
    procedure, public :: FinalizeTimestep => PMBaseThisOnly
    procedure, public :: AcceptSolution => PMBaseFunctionThisOnly
    procedure, public :: CheckUpdatePre => PMBaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMBaseCheckUpdatePost
    procedure, public :: TimeCut => PMBaseThisOnly
    procedure, public :: UpdateSolution => PMBaseThisOnly
    procedure, public :: MaxChange => PMBaseThisOnly
    procedure, public :: ComputeMassBalance => PMBaseComputeMassBalance
    procedure, public :: Destroy => PMBaseThisOnly
#endif
  end type process_model_base_type
  
#ifdef ABSTRACT  
  abstract interface
    subroutine PMBaseInit(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseInit

    subroutine PMBaseResidual(this,snes,xx,r,ierr)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      SNES :: snes
      Vec :: xx
      Vec :: r
      PetscErrorCode :: ierr
    end subroutine PMBaseResidual

    subroutine PMBaseJacobian(this,snes,xx,A,B,flag,ierr)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      SNES :: snes
      Vec :: xx
      Mat :: A, B
      MatStructure flag
      PetscErrorCode :: ierr
    end subroutine PMBaseJacobian
    
    subroutine PMBaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscReal :: dt
      PetscReal :: dt_max
      PetscInt :: iacceleration
      PetscInt :: num_newton_iterations
      PetscReal :: tfac(:)
    end subroutine PMBaseUpdateTimestep
    
    subroutine PMBaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      SNESLineSearch :: line_search
      Vec :: P
      Vec :: dP
      PetscBool :: changed
      PetscErrorCode :: ierr
    end subroutine PMBaseCheckUpdatePre
    
    subroutine PMBaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                      P1_changed,ierr)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      SNESLineSearch :: line_search
      Vec :: P0
      Vec :: dP
      Vec :: P1
      PetscBool :: dP_changed
      PetscBool :: P1_changed
      PetscErrorCode :: ierr
    end subroutine PMBaseCheckUpdatePost
  
    subroutine PMBasePostSolve(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscBool :: solution_accepted
    end subroutine PMBasePostSolve
    
    subroutine PMBaseThisOnly(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseThisOnly
    
    subroutine PMBaseThisTime(this,time)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscReal :: time
    end subroutine PMBaseThisTime
    
    function PMBaseFunctionThisOnly(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscBool ::  PMBaseFunctionThisOnly
    end function PMBaseFunctionThisOnly
    
    subroutine PMBaseComputeMassBalance(this,mass_balance_array)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscReal :: mass_balance_array(:)
    end subroutine PMBaseComputeMassBalance

  end interface
#endif  
  
  public :: PMBaseCreate
  
  public :: PMBaseResidual
  public :: PMBaseJacobian
  
contains

subroutine PMBaseCreate(this)

  implicit none
  
  class(process_model_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%option)
  nullify(this%output_option)
  nullify(this%realization_base)
  this%solution_vec = 0
  this%residual_vec = 0
  nullify(this%next)
  
end subroutine PMBaseCreate

#if 0
! ************************************************************************** !
!
! PMBaseRunTo: Runs the actual simulation.
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMBaseRunTo(this,time)

  implicit none
  
  class(process_model_base_type) :: this  
  PetscReal :: time
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%RunTo(time)
  endif
  
end subroutine PMBaseRunTo
#endif

#ifndef ABSTRACT
subroutine PMBaseInit(this)
  implicit none
  class(process_model_base_type) :: this
end subroutine PMBaseInit

subroutine PMBaseResidual(this,snes,xx,r,ierr)
  implicit none
  class(process_model_base_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
end subroutine PMBaseResidual

subroutine PMBaseJacobian(this,snes,xx,A,B,flag,ierr)
  implicit none
  class(process_model_base_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
end subroutine PMBaseJacobian
    
subroutine PMBaseUpdateTimestep(this,dt,dt_max,iacceleration, &
                                num_newton_iterations,tfac)
  implicit none
  class(process_model_base_type) :: this
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
end subroutine PMBaseUpdateTimestep
    
subroutine PMBaseCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  implicit none
  class(process_model_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
end subroutine PMBaseCheckUpdatePre
    
subroutine PMBaseCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  class(process_model_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
end subroutine PMBaseCheckUpdatePost
  
subroutine PMBasePostSolve(this)
  implicit none
  class(process_model_base_type) :: this
  PetscBool :: solution_accepted
end subroutine PMBasePostSolve
    
subroutine PMBaseThisOnly(this)
  implicit none
  class(process_model_base_type) :: this
end subroutine PMBaseThisOnly
    
subroutine PMBaseThisTime(this,time)
  implicit none
  class(process_model_base_type) :: this
  PetscReal :: time
end subroutine PMBaseThisTime
    
function PMBaseFunctionThisOnly(this)
  implicit none
  class(process_model_base_type) :: this
  PetscBool ::  PMBaseFunctionThisOnly
  PMBaseFunctionThisOnly = PETSC_TRUE
end function PMBaseFunctionThisOnly
    
subroutine PMBaseComputeMassBalance(this,mass_balance_array)
  implicit none
  class(process_model_base_type) :: this
  PetscReal :: mass_balance_array(:)
end subroutine PMBaseComputeMassBalance
#endif

end module Process_Model_Base_class
