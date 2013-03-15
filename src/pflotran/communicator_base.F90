module Communicator_Base_module

  implicit none

  private

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, abstract, public :: communicator_type
  contains
    procedure(VecToVec), public, deferred :: GlobalToLocal 
    procedure(VecToVec), public, deferred :: LocalToGlobal
    procedure(VecToVec), public, deferred :: LocalToLocal 
    procedure(VecToVec), public, deferred :: GlobalToNatural 
    procedure(VecToVec), public, deferred :: NaturalToGlobal 
    procedure(BaseDestroy), public, deferred :: Destroy 
  end type communicator_type
  
  abstract interface
  
    subroutine VecToVec(this,source,destination)
      import communicator_type
      implicit none
      class(communicator_type) :: this
      Vec :: source
      Vec :: destination
    end subroutine VecToVec

    subroutine BaseDestroy(this)
      import communicator_type
      implicit none
      class(communicator_type) :: this
    end subroutine BaseDestroy

  end interface
  
end module Communicator_Base_module
