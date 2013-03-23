module Process_Model_Base_class

  use Option_module

  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
  
  type, abstract, public :: process_model_base_type
    type(option_type), pointer :: option
    class(process_model_base_type), pointer :: next
  contains
    procedure(PMBaseInit), public, deferred :: Init
    procedure(PMBaseThisOnly), public, deferred :: InitializeRun
    procedure(PMBaseThisOnly), public, deferred :: FinalizeRun
    procedure(PMBaseResidual), public, deferred :: Residual
    procedure(PMBaseJacobian), public, deferred :: Jacobian
    procedure(PMBaseUpdateTimestep), public, deferred :: UpdateTimestep
    procedure(PMBaseThisOnly), public, deferred :: UpdatePreSolve
    procedure(PMBaseCheckUpdatePre), public, deferred :: CheckUpdatePre
    procedure(PMBaseCheckUpdatePost), public, deferred :: CheckUpdatePost
    procedure(PMBaseThisOnly), public, deferred :: TimeCut
    procedure(PMBaseThisOnly), public, deferred :: UpdateSolution
    procedure(PMBaseThisOnly), public, deferred :: MaxChange
    procedure(PMBaseComputeMassBalance), public, deferred :: ComputeMassBalance
    procedure(PMBaseThisOnly), public, deferred :: Destroy
  end type process_model_base_type
  
  interface
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
  
    subroutine PMBaseThisOnly(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseThisOnly
    
    subroutine PMBaseComputeMassBalance(this,mass_balance_array)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscReal :: mass_balance_array(:)
    end subroutine PMBaseComputeMassBalance

  end interface
  
  public :: PMBaseCreate
  
contains

subroutine PMBaseCreate(this)

  implicit none
  
  class(process_model_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%option)
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
end module Process_Model_Base_class
