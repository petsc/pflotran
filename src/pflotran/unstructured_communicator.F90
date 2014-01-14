module Unstructured_Communicator_class

  use Communicator_Base_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Explicit_module  
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
  
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

function UnstructuredCommunicatorCreate()
  ! 
  ! Allocates and initializes a new communicator
  ! object for unstructured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

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

subroutine UnstructuredSetDM(this,dm_ptr)
  ! 
  ! Sets pointer to DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use DM_Kludge_module

  implicit none
  
  class(unstructured_communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm = dm_ptr%dm
  this%ugdm => dm_ptr%ugdm
  
end subroutine UnstructuredSetDM

! ************************************************************************** !

subroutine UnstructuredGlobalToLocal(this,source,destination)
  ! 
  ! Performs global to local communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

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

subroutine UnstructuredLocalToGlobal(this,source,destination)
  ! 
  ! Performs local to global communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

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

subroutine UnstructuredLocalToLocal(this,source,destination)
  ! 
  ! Performs local to local communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call VecScatterBegin(this%ugdm%scatter_ltol,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(this%ugdm%scatter_ltol,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  
!  call DMLocalToLocalBegin(this%dm,source,INSERT_VALUES,destination,ierr)
!  call DMLocalToLocalEnd(this%dm,source,INSERT_VALUES,destination,ierr)
  
end subroutine UnstructuredLocalToLocal

! ************************************************************************** !

subroutine UnstructuredGlobalToNatural(this,source,destination)
  ! 
  ! Performs global to natural communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

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

subroutine UnstructuredNaturalToGlobal(this,source,destination)
  ! 
  ! Performs natural to global communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

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

subroutine UnstructuredCommunicatorDestroy(this)
  ! 
  ! Deallocates a communicator object for
  ! unstructured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  
  PetscErrorCode :: ierr
  
  if (associated(this%ugdm)) then
    call UGridDMDestroy(this%ugdm)
  endif
  nullify(this%ugdm)
  if (this%dm /= 0) then
    !geh: all DMs are currently destroyed in realization.  This DM is solely
    !     a pointer.  This will need to change, but skip for now.
    !call DMDestroy(this%dm,ierr)
  endif
  this%dm = 0  
  
end subroutine UnstructuredCommunicatorDestroy

end module Unstructured_Communicator_class
