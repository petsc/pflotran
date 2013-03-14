module Process_Model_Base_class

  ! add accessory modules here

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
    procedure, public :: Init => PMBaseInit
    procedure, public :: InitializeTimeStep => PMBaseInitializeTimeStep
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: CheckUpdatePre => PMBaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMBaseCheckUpdatePost
    procedure, public :: TimeCut => PMBaseTimeCut
    procedure, public :: UpdateSolution => PMBaseUpdateSolution
    procedure, public :: MaxChange => PMBaseMaxChange
    procedure, public :: ComputeMassBalance => PMBaseComputeMassBalance
    procedure, public :: Destroy => PMBaseDestroy
  end type process_model_base_type
  
  interface
    subroutine PMBaseInit(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseInit

    subroutine PMBaseInitializeTimestep(this)
      import process_model_base_type
      implicit none
      class(process_model_base_type) :: this
    end subroutine PMBaseInitializeTimestep 

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

end module Process_Model_Base_class
