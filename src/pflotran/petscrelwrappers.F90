module PetscRelWrappers
  private 
  
  public :: VecScatterBegin_wrap, VecScatterEnd_wrap

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

contains

#if (PETSC_VERSION_RELEASE == 1)

subroutine VecScatterBegin_wrap(inctx,x,y,addv,mode,ierr)
  implicit none
  
  VecScatter, intent(in) :: inctx
  Vec, intent(in) :: x
  Vec, intent(out) :: y
  InsertMode, intent(in) :: addv
  ScatterMode, intent(in) :: mode
  PetscErrorCode, intent(out) :: ierr

  call VecScatterBegin(x,y,addv,mode,inctx,ierr)
end subroutine VecScatterBegin_wrap


subroutine VecScatterEnd_wrap(inctx,x,y,addv,mode,ierr)
  implicit none

  VecScatter, intent(in) :: inctx
  Vec, intent(in) :: x
  Vec, intent(out) :: y
  InsertMode, intent(in) :: addv
  ScatterMode, intent(in) :: mode
  PetscErrorCode, intent(out) :: ierr

  call VecScatterEnd(x,y,addv,mode,inctx,ierr)
end subroutine VecScatterEnd_wrap

#endif

end module PetscRelWrappers
