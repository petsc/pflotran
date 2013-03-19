module Process_Model_Base_class

  use Timestepper_module

  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
  
  type, abstract, public :: process_model_base_type
    class(process_model_base_type), pointer :: next
  contains
    procedure(PMBaseInit), public, deferred :: Init
    procedure(PMBaseInitializeRun), public, deferred :: InitializeRun
    procedure(PMBaseFinalizeRun), public, deferred :: FinalizeRun
    procedure(PMBaseResidual), public, deferred :: Residual
    procedure(PMBaseJacobian), public, deferred :: Jacobian
    procedure(PMBaseCheckUpdatePre), public, deferred :: CheckUpdatePre
    procedure(PMBaseCheckUpdatePost), public, deferred :: CheckUpdatePost
    procedure(PMBaseTimeCut), public, deferred :: TimeCut
    procedure(PMBaseUpdateSolution), public, deferred :: UpdateSolution
    procedure(PMBaseMaxChange), public, deferred :: MaxChange
    procedure(PMBaseComputeMassBalance), public, deferred :: ComputeMassBalance
    procedure(PMBaseDestroy), public, deferred :: Destroy
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
  
    subroutine PMBaseTimeCut(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseTimeCut
    
    subroutine PMBaseUpdateSolution(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseUpdateSolution     

    subroutine PMBaseMaxChange(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseMaxChange
    
    subroutine PMBaseComputeMassBalance(this,mass_balance_array)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
      PetscReal :: mass_balance_array(:)
    end subroutine PMBaseComputeMassBalance
    
    subroutine PMBaseDestroy(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseDestroy
  end interface
  
contains

! ************************************************************************** !
!
! PMBaseRunTo: Runs the actual simulation.
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMBaseRunTo(this,time)

  implicit none
  
  type(process_model_type) :: this  
  PetscReal :: time
  
  ! do something here
  
  if (this%next) then
    this%next%RunTo(time)
  endif
  
end subroutine PMBaseRunTo

end module Process_Model_Base_class
