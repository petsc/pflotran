module Unstructured_Communicator_class

  use Communicator_Base_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Explicit_module  
  
  implicit none

  private

#include "definitions.h"
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmshell.h90"


  type, public, extends(communicator_type) :: unstructured_communicator_type
    DM :: dm
    type(ugdm_type), pointer :: ugdm
  contains
    procedure, public :: SetDM => UnstructuredSetDM
    procedure, public :: GlobalToLocal => UnstructuredGlobalToLocal
    procedure, public :: LocalToGlobal => UnstructuredLocalToGlobal
    procedure, public :: LocalToLocal => UnstructuredLocalToLocal
    procedure, public :: GlobalToNatural => UnstructuredGlobalToNatural
    procedure, public :: NaturalToGlobal => UnstructuredNaturalToGlobal
!geh: finalization not yet supported by gfortran.
!    final :: UnstructuredCommunicatorDestroy
    procedure, public :: Destroy => UnstructuredCommunicatorDestroy
  end type unstructured_communicator_type

  public :: UnstructuredCommunicatorCreate
  
contains

! ************************************************************************** !
!
! UnstructuredCommunicatorCreate: Allocates and initializes a new communicator 
!                               object for unstructured grids
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
function UnstructuredCommunicatorCreate()

  implicit none
  
  class(unstructured_communicator_type), pointer :: &
    UnstructuredCommunicatorCreate
  
  class(unstructured_communicator_type), pointer :: communicator
  
  allocate(communicator)
  nullify(communicator%ugdm)
  communicator%dm = 0

  UnstructuredCommunicatorCreate => communicator  
  
end function UnstructuredCommunicatorCreate

! ************************************************************************** !
!
! UnstructuredSetDM: Sets pointer to DM
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
subroutine UnstructuredSetDM(this,dm_ptr)

  use Discretization_module

  implicit none
  
  class(unstructured_communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm = dm_ptr%dm
  this%ugdm => dm_ptr%ugdm
  
end subroutine UnstructuredSetDM

! ************************************************************************** !
!
! UnstructuredGlobalToLocal: Performs global to local communication 
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredGlobalToLocal(this,source,destination)

!TODO(geh): move to communicator_base.F90

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMGlobalToLocalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
  call DMGlobalToLocalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredGlobalToLocal

! ************************************************************************** !
!
! UnstructuredLocalToGlobal: Performs local to global communication 
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredLocalToGlobal(this,source,destination)

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_ltog,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(this%ugdm%scatter_ltog,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
      
!  call DMLocalToGlobalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
!  call DMLocalToGlobalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredLocalToGlobal

! ************************************************************************** !
!
! UnstructuredLocalToLocal: Performs local to local communication 
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredLocalToLocal(this,source,destination)

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call VecScatterBegin(this%ugdm%scatter_ltol,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(this%ugdm%scatter_ltol,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
!  call DMDALocalToLocalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
!  call DMDALocalToLocalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredLocalToLocal

! ************************************************************************** !
!
! UnstructuredGlobalToNatural: Performs global to natural communication 
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredGlobalToNatural(this,source,destination)

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_gton,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(this%ugdm%scatter_gton,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
!  call DMDAGlobalToNaturalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
!  call DMDAGlobalToNaturalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredGlobalToNatural

! ************************************************************************** !
!
! UnstructuredNaturalToGlobal: Performs natural to global communication 
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredNaturalToGlobal(this,source,destination)

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_ntog,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(this%ugdm%scatter_ntog,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
!  call DMDANaturalToGlobalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
!  call DMDANaturalToGlobalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredNaturalToGlobal

! ************************************************************************** !
!
! UnstructuredCommunicatorDestroy: Deallocates a communicator object for 
!                                  unstructured grids
! author: Glenn Hammond
! date: 03/15/13
!
! ************************************************************************** !
subroutine UnstructuredCommunicatorDestroy(this)

  implicit none
  
  class(unstructured_communicator_type) :: this
  
  PetscErrorCode :: ierr
  
  if (associated(this%ugdm)) then
    call UGridDMDestroy(this%ugdm)
  endif
  nullify(this%ugdm)
  if (this%dm /= 0) then
    call DMDestroy(this%dm,ierr)
  endif
  this%dm = 0  
  
end subroutine UnstructuredCommunicatorDestroy

end module Unstructured_Communicator_class
