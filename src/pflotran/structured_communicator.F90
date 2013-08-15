module Structured_Communicator_class

  use Communicator_Base_module
  use Structured_Grid_module  
  
  implicit none

  private

#include "definitions.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"

  type, public, extends(communicator_type) :: structured_communicator_type
    DM :: dm
  contains
    procedure, public :: SetDM => StructuredSetDM
    procedure, public :: GlobalToLocal => StructuredGlobalToLocal
    procedure, public :: LocalToGlobal => StructuredLocalToGlobal
    procedure, public :: LocalToLocal => StructuredLocalToLocal
    procedure, public :: GlobalToNatural => StructuredGlobalToNatural
    procedure, public :: NaturalToGlobal => StructuredNaturalToGlobal
!geh: finalization not yet supported by gfortran.
!    final :: StructuredCommunicatorDestroy
    procedure, public :: Destroy => StructuredCommunicatorDestroy

  end type structured_communicator_type
  
  public :: StructuredCommunicatorCreate
  
contains

! ************************************************************************** !
!
! StructuredCommunicatorCreate: Allocates and initializes a new communicator 
!                               object for structured grids
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
function StructuredCommunicatorCreate()

  implicit none
  
  class(structured_communicator_type), pointer :: StructuredCommunicatorCreate
  
  class(structured_communicator_type), pointer :: communicator
  
  allocate(communicator)
  communicator%dm = 0

  StructuredCommunicatorCreate => communicator  
  
end function StructuredCommunicatorCreate

! ************************************************************************** !
!
! StructuredSetDM: Sets pointer to DM
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine StructuredSetDM(this,dm_ptr)

  use DM_Kludge_module

  implicit none
  
  class(structured_communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm = dm_ptr%dm
  
end subroutine StructuredSetDM

! ************************************************************************** !
!
! StructuredGlobalToLocal: Performs global to local communication with DM
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredGlobalToLocal(this,source,destination)

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMGlobalToLocalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMGlobalToLocalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine StructuredGlobalToLocal

! ************************************************************************** !
!
! StructuredLocalToGlobal: Performs local to global communication with DM
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredLocalToGlobal(this,source,destination)

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMLocalToGlobalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMLocalToGlobalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine StructuredLocalToGlobal

! ************************************************************************** !
!
! StructuredLocalToLocal: Performs local to local communication with DM
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredLocalToLocal(this,source,destination)

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMLocalToLocalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMLocalToLocalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine StructuredLocalToLocal

! ************************************************************************** !
!
! StructuredGlobalToNatural: Performs global to natural communication with DM
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredGlobalToNatural(this,source,destination)

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMDAGlobalToNaturalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMDAGlobalToNaturalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine StructuredGlobalToNatural

! ************************************************************************** !
!
! StructuredNaturalToGlobal: Performs natural to global communication with DM
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredNaturalToGlobal(this,source,destination)

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMDANaturalToGlobalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMDANaturalToGlobalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine StructuredNaturalToGlobal

! ************************************************************************** !
!
! StructuredCommunicatorDestroy: Deallocates a communicator object for 
!                                structured grids
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine StructuredCommunicatorDestroy(this)

  implicit none
  
  class(structured_communicator_type) :: this
  
  PetscErrorCode :: ierr
  
  if (this%dm /= 0) then
    call DMDestroy(this%dm,ierr)
  endif
  this%dm = 0
  
end subroutine StructuredCommunicatorDestroy

end module Structured_Communicator_class
