module Process_Model_RT_class

  use Process_Model_Base_class
  use Reactive_Transport_module
  use Realization_class
  
  implicit none

  private

#include "definitions.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"

  type, public, extends(process_model_base_type) :: process_model_rt_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMRTInit
    procedure, public :: InitializeTimeStep => PMRTInitializeTimestep
    procedure, public :: Residual => PMRTResidual
    procedure, public :: Jacobian => PMRTJacobian
    procedure, public :: CheckUpdatePre => PMRTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRTCheckUpdatePost
    procedure, public :: TimeCut => PMRTTimeCut
    procedure, public :: UpdateSolution => PMRTUpdateSolution
    procedure, public :: MaxChange => PMRTMaxChange
    procedure, public :: ComputeMassBalance => PMRTComputeMassBalance
    procedure, public :: Destroy => PMRTDestroy
  end type process_model_rt_type
  
contains

! ************************************************************************** !
!
! PMRTInit: Initializes variables associated with Richard
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInit(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTSetup(this%realization)
  
end subroutine PMRTInit

! ************************************************************************** !
!
! PMRTInitializeTimestep: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTInitializeTimestep(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTSetup(this%realization)

end subroutine PMRTInitializeTimestep 

! ************************************************************************** !
!
! PMRTResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTResidual(this,snes,xx,r,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call RTResidual(snes,xx,r,this%realization,ierr)
  
end subroutine PMRTResidual

! ************************************************************************** !
!
! PMRTJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTJacobian(this,snes,xx,A,B,flag,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  PetscErrorCode :: ierr
  
  call RTJacobian(snes,xx,A,B,flag,this%realization,ierr)
  
end subroutine PMRTJacobian
    
! ************************************************************************** !
!
! PMRTCheckUpdatePre: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePre(this,line_search,P,dP,changed,ierr)

  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P
  Vec :: dP
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call RTCheckUpdate(line_search,P,dP,changed,this%realization,ierr)
  
end subroutine PMRTCheckUpdatePre
    
! ************************************************************************** !
!
! PMRTCheckUpdatePost: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
                                  P1_changed,ierr)
  implicit none
  
  class(process_model_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dP_changed
  PetscBool :: P1_changed
  PetscErrorCode :: ierr
  
!  call RTCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
!                               P1_changed,this%realization,ierr)

end subroutine PMRTCheckUpdatePost
  
! ************************************************************************** !
!
! PMRTTimeCut: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTTimeCut(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTTimeCut(this%realization)

end subroutine PMRTTimeCut
    
! ************************************************************************** !
!
! PMRTUpdateSolution: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTUpdateSolution(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTUpdateSolution(this%realization)

end subroutine PMRTUpdateSolution     

! ************************************************************************** !
!
! PMRTMaxChange: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTMaxChange(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTMaxChange(this%realization)

end subroutine PMRTMaxChange
    
! ************************************************************************** !
!
! PMRTComputeMassBalance: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTComputeMassBalance(this,mass_balance_array)

  implicit none
  
  class(process_model_rt_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call RTComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRTComputeMassBalance

! ************************************************************************** !
!
! PMRTDestroy: Destroys RT process model
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMRTDestroy(this)

  implicit none
  
  class(process_model_rt_type) :: this
  
  call RTDestroy(this%realization)
  
end subroutine PMRTDestroy
  
end module Process_Model_RT_class
