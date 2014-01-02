module PM_module

  use PM_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  ! Since the context (ctx) for procedures passed to PETSc must be declared 
  ! as a "type" instead of a "class", object is a workaround for passing the 
  ! process model as context of a procedure where one can pass the
  ! pm_pointer_type to a procedure, declaring it as e.g.
  !
  ! type(pm_pointer_type) :: pm_ptr
  !
  ! and use the ptr:
  !
  ! pm_ptr%this%Residual
  !  
  type, public :: pm_pointer_type
    class(pm_base_type), pointer :: ptr
  end type pm_pointer_type

  public :: PMResidual, &
            PMJacobian, &
            PMRHSFunction

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
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(pm_pointer_type) :: this
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
  type(pm_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef PM_TOP_DEBUG    
  print *, 'PMJacobian()'
#endif

  call this%ptr%Jacobian(snes,xx,A,B,flag,ierr)
    
end subroutine PMJacobian
    
! ************************************************************************** !
!
! PMRHSFunction:
! author: Gautam Bisht
! date: 04/12/13
!
! ************************************************************************** !
subroutine PMRHSFunction(ts,time,xx,ff,this,ierr)

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscts.h"

  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  type(pm_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef PM_TOP_DEBUG
  print *, 'PMRHSFunction()'
#endif

  call this%ptr%RHSFunction(ts,time,xx,ff,ierr)

end subroutine PMRHSFunction

end module PM_module
