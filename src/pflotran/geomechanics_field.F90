#ifdef GEOMECH

module Geomechanics_Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: geomechanics_field_type
    Vec :: work
    Vec :: work_loc
    ! residual vectors
    Vec :: disp_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: disp_xx, disp_xx_loc, disp_dxx, disp_yy, disp_accum

  end type geomechanics_field_type

  public :: GeomechanicsFieldCreate, &
            GeomechanicsFieldDestroy

contains

! ************************************************************************** !
!
! GeomechanicsFieldCreate: Allocates and initializes a new geomechanics
!                          Field object
! author: Satish Karra, LANL
! date: 06/05/13
!
! ************************************************************************** !
function GeomechanicsFieldCreate()

  implicit none
  
  type(geomechanics_field_type), pointer     :: GeomechanicsFieldCreate
  type(geomechanics_field_type), pointer     :: geomechanics_field
  
  allocate(geomechanics_field)

  ! nullify PetscVecs
  geomechanics_field%work = 0
  geomechanics_field%work_loc = 0
  
  geomechanics_field%disp_r = 0
  geomechanics_field%disp_xx = 0
  geomechanics_field%disp_xx_loc = 0
  geomechanics_field%disp_dxx = 0
  geomechanics_field%disp_yy = 0
  geomechanics_field%disp_accum = 0
  
  GeomechanicsFieldCreate => geomechanics_field

end function GeomechanicsFieldCreate

! ************************************************************************** !
!
! GeomechanicsFieldDestroy: Deallocates a geomechanics field object
! author: Satish Karra, LANL
! date: 06/05/13
!
! ************************************************************************** !
subroutine GeomechanicsFieldDestroy(geomechanics_field)

  implicit none
  
  type(geomechanics_field_type), pointer     :: geomechanics_field
  PetscErrorCode                             :: ierr
  
  ! Destroy PetscVecs
  if (geomechanics_field%work /= 0) call VecDestroy(geomechanics_field%work,ierr)
  if (geomechanics_field%work_loc  /= 0) call VecDestroy(geomechanics_field%work_loc,ierr)

  if (geomechanics_field%disp_r /= 0) call VecDestroy(geomechanics_field%disp_r,ierr)
  if (geomechanics_field%disp_xx /= 0) call VecDestroy(geomechanics_field%disp_xx,ierr)
  if (geomechanics_field%disp_xx_loc /= 0) call VecDestroy(geomechanics_field%disp_xx_loc,ierr)
  if (geomechanics_field%disp_dxx /= 0) call VecDestroy(geomechanics_field%disp_dxx,ierr)
  if (geomechanics_field%disp_yy /= 0) call VecDestroy(geomechanics_field%disp_yy,ierr)
  if (geomechanics_field%disp_accum /= 0) call VecDestroy(geomechanics_field%disp_accum,ierr)
  
end subroutine GeomechanicsFieldDestroy

end module Geomechanics_Field_module

#endif
