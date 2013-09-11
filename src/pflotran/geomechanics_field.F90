#ifdef GEOMECH

module Geomechanics_Field_module

  use PFLOTRAN_Constants_module
! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: geomech_field_type
    Vec :: work
    Vec :: work_loc
    ! residual vectors
    Vec :: disp_r
    Vec :: press      ! store pressure from subsurf
    Vec :: press_loc
    Vec :: temp_loc
    Vec :: temp       ! store temperature from subsurf
    Vec :: subsurf_vec_1dof ! MPI
    Vec :: imech_loc

    ! Solution vectors (xx = current iterate)
    Vec :: disp_xx, disp_xx_loc

  end type geomech_field_type

  public :: GeomechFieldCreate, &
            GeomechFieldDestroy

contains

! ************************************************************************** !
!
! GeomechFieldCreate: Allocates and initializes a new geomechanics
!                     Field object
! author: Satish Karra, LANL
! date: 06/05/13
!
! ************************************************************************** !
function GeomechFieldCreate()

  implicit none
  
  type(geomech_field_type), pointer     :: GeomechFieldCreate
  type(geomech_field_type), pointer     :: geomech_field
  
  allocate(geomech_field)

  ! nullify PetscVecs
  geomech_field%work = 0
  geomech_field%work_loc = 0
  
  geomech_field%disp_r = 0
  geomech_field%disp_xx = 0
  geomech_field%disp_xx_loc = 0
  
  geomech_field%press = 0
  geomech_field%press_loc = 0
  geomech_field%temp = 0
  geomech_field%temp_loc = 0
  geomech_field%subsurf_vec_1dof = 0
  geomech_field%imech_loc = 0
  
  GeomechFieldCreate => geomech_field

end function GeomechFieldCreate

! ************************************************************************** !
!
! GeomechFieldDestroy: Deallocates a geomechanics field object
! author: Satish Karra, LANL
! date: 06/05/13
!
! ************************************************************************** !
subroutine GeomechFieldDestroy(geomech_field)

  implicit none
  
  type(geomech_field_type), pointer     :: geomech_field
  PetscErrorCode                        :: ierr
  
  ! Destroy PetscVecs
  if (geomech_field%work /= 0) call VecDestroy(geomech_field%work,ierr)
  if (geomech_field%work_loc  /= 0) call VecDestroy(geomech_field%work_loc,ierr)

  if (geomech_field%disp_r /= 0) call VecDestroy(geomech_field%disp_r,ierr)
  if (geomech_field%disp_xx /= 0) call VecDestroy(geomech_field%disp_xx,ierr)
  if (geomech_field%disp_xx_loc /= 0) call VecDestroy(geomech_field%disp_xx_loc,ierr)
  
  if (geomech_field%press /= 0) call VecDestroy(geomech_field%press,ierr)
  if (geomech_field%press_loc /= 0) call VecDestroy(geomech_field%press_loc,ierr)
  if (geomech_field%temp /= 0) call VecDestroy(geomech_field%temp,ierr)
  if (geomech_field%temp_loc /= 0) call VecDestroy(geomech_field%temp_loc,ierr)

  if (geomech_field%subsurf_vec_1dof /= 0 ) &
    call VecDestroy(geomech_field%subsurf_vec_1dof,ierr)
  if (geomech_field%imech_loc /= 0) call VecDestroy(geomech_field%imech_loc,ierr)

  
end subroutine GeomechFieldDestroy

end module Geomechanics_Field_module

#endif
