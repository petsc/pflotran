#ifdef SURFACE_FLOW 

module Surface_Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: surface_field_type 

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_r

  end type surface_field_type

  public :: SurfaceFieldCreate, &
            SurfaceFieldDestroy

contains

! ************************************************************************** !
!
! SurfaceFieldCreate: Allocates and initializes a new surface Field object
! author: Gautam Bisht
! date: 01/17/2012
!
! ************************************************************************** !
function SurfaceFieldCreate()

  implicit none
  
  type(surface_field_type), pointer :: SurfaceFieldCreate
  
  type(surface_field_type), pointer :: surface_field
  
  allocate(surface_field)

  ! nullify PetscVecs
  surface_field%flow_r = 0

  SurfaceFieldCreate => surface_field

end function SurfaceFieldCreate

! ************************************************************************** !
!
! SurfaceFieldDestroy: Deallocates a field object
! author: Gautam Bisht
! date: 01/17/2012
!
! ************************************************************************** !
subroutine SurfaceFieldDestroy(surface_field)

  implicit none
  
  type(surface_field_type), pointer :: surface_field
  
  PetscErrorCode :: ierr

  ! Destroy PetscVecs
  if (surface_field%flow_r /= 0) call VecDestroy(surface_field%flow_r,ierr)

end subroutine SurfaceFieldDestroy

end module Surface_Field_module

#endif