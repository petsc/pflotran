module Process_Model_module

  use Process_Model_Base_class
  use Process_Model_Richards_class
  use Process_Model_RT_class
  
  implicit none

  private

#include "definitions.h"

  ! Since the context (ctx) for procedures passed to PETSc must be declared 
  ! as a "type" instead of a "class", object is a workaround for passing the 
  ! process model as context of a procedure where one can pass the
  ! process_model_pointer_type to a procedure, declaring it as e.g.
  !
  ! type(process_model_pointer_type) :: pm_ptr
  !
  ! and use the ptr:
  !
  ! pm_ptr%this%Residual
  !  
  type, public :: process_model_pointer_type
    class(process_model_base_type), pointer :: ptr
  end type process_model_pointer_type

  public :: PMResidual, &
            PMJacobian

contains

! ************************************************************************** !
!
! PMResidual: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMResidual(snes,xx,r,this,ierr)

  use Option_module
  use Realization_class
  
  use Richards_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(process_model_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef PM_TOP_DEBUG    
  print *, 'PMResidual()'
#endif

  call this%ptr%Residual(snes,xx,r,ierr)

end subroutine PMResidual

! ************************************************************************** !
!
! PMJacobian: 
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine PMJacobian(snes,xx,A,B,flag,this,ierr)

  use Option_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  MatStructure flag
  type(process_model_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef PM_TOP_DEBUG    
  print *, 'PMJacobian()'
#endif

  call this%ptr%Jacobian(snes,xx,A,B,flag,ierr)
    
end subroutine PMJacobian
    
end module Process_Model_module
