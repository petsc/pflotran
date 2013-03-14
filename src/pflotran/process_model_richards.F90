module Process_Model_Richards_class

  use Process_Model_Base_class
  use Richards_module
  use Realization_class
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(process_model_base_type) :: process_model_richards_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMRichardsInit
    procedure, public :: InitializeTimeStep => PMRichardsInitializeTimestep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: Destroy => PMRichardsDestroy
  end type process_model_richards_type
  
contains

! ************************************************************************** !
!
! PMRichardsInit: Initializes variables associated with Richard
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsInit(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsSetup(this%realization)
  
end subroutine PMRichardsInit

! ************************************************************************** !
!
! PMRichardsInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsInitializeTimestep(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsSetup(this%realization)

end subroutine PMRichardsInitializeTimestep 

! ************************************************************************** !
!
! PMRichardsResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsResidual(this,snes,xx,r,ierr)
  implicit none
  class(process_model_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call RichardsResidual(snes,xx,r,this%realization,ierr)
  
end subroutine PMRichardsResidual

! ************************************************************************** !
!
! PMRichardsJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsJacobian(this,snes,xx,A,B,flag,ierr)
  implicit none
  class(process_model_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
  call RichardsJacobian(snes,xx,A,B,flag,this%realization,ierr)
  
end subroutine PMRichardsJacobian
    
! ************************************************************************** !
!
! PMRichardsCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsCheckUpdatePre(this,line_search,P,dP,changed,ierr)
  implicit none
  class(process_model_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call RichardsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
  
end subroutine PMRichardsCheckUpdatePre
    
! ************************************************************************** !
!
! PMRichardsCheckUpdatePost: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  class(process_model_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
  call RichardsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
                               P1_changed,this%realization,ierr)

end subroutine PMRichardsCheckUpdatePost
  
! ************************************************************************** !
!
! PMRichardsTimeCut: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsTimeCut(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsTimeCut(this%realization)

end subroutine PMRichardsTimeCut
    
! ************************************************************************** !
!
! PMRichardsUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsUpdateSolution(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsUpdateSolution(this%realization)

end subroutine PMRichardsUpdateSolution     

! ************************************************************************** !
!
! PMRichardsMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsMaxChange(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsMaxChange(this%realization)

end subroutine PMRichardsMaxChange
    
! ************************************************************************** !
!
! PMRichardsComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsComputeMassBalance(this,mass_balance_array)

  implicit none
  
  class(process_model_richards_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call RichardsComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !
!
! PMRichardsDestroy: Destroys Richards process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRichardsDestroy(this)

  implicit none
  
  class(process_model_richards_type) :: this
  
  call RichardsDestroy(this%realization)
  
end subroutine PMRichardsDestroy
  
end module Process_Model_Richards_class
