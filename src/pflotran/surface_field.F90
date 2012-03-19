#ifdef SURFACE_FLOW

module Surface_Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: surface_field_type

    Vec :: mannings0, mannings_loc

    Vec :: work, work_loc

    ! residual vectors
    Vec :: flow_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc

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
  surface_field%mannings0 = 0
  surface_field%mannings_loc = 0
  surface_field%work = 0
  surface_field%work_loc = 0
  surface_field%flow_r = 0
  surface_field%flow_xx = 0
  surface_field%flow_xx_loc = 0

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
  if (surface_field%mannings0 /= 0) call VecDestroy(surface_field%mannings0,ierr)
  if (surface_field%work /= 0) call VecDestroy(surface_field%work,ierr)
  if (surface_field%work_loc  /= 0) call VecDestroy(surface_field%work_loc,ierr)
  if (surface_field%mannings_loc /= 0) call VecDestroy(surface_field%mannings_loc,ierr)
  if (surface_field%flow_r /= 0) call VecDestroy(surface_field%flow_r,ierr)
  if (surface_field%flow_xx /= 0) call VecDestroy(surface_field%flow_xx,ierr)
  if (surface_field%flow_xx_loc /= 0) call VecDestroy(surface_field%flow_xx_loc,ierr)

end subroutine SurfaceFieldDestroy

end module Surface_Field_module

#endif