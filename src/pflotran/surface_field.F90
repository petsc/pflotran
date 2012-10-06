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

    Vec :: area
    
    Vec :: vol_subsurf_2_surf    ! MPI +ve value => Flow from subsurface to surface
    Vec :: press_subsurf         ! MPI

    Vec :: Dq                    ! MPI
    Vec :: por                   ! MPI
    Vec :: sat_func_id           ! MPI
    Vec :: perm_xx               ! MPI
    Vec :: perm_yy               ! MPI
    Vec :: perm_zz               ! MPI
    Vec :: subsurf_xx            ! MPI
    Vec :: subsurf_yy            ! MPI
    Vec :: subsurf_zz            ! MPI

    Vec :: subsurf_temp_vec      ! MPI

    ! residual vectors
    Vec :: flow_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum

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

  surface_field%area = 0
  
  surface_field%flow_r = 0
  surface_field%flow_xx = 0
  surface_field%flow_xx_loc = 0
  surface_field%flow_dxx = 0
  surface_field%flow_yy = 0
  surface_field%flow_accum = 0
  
  surface_field%vol_subsurf_2_surf = 0
  surface_field%press_subsurf = 0

  surface_field%Dq = 0
  surface_field%por = 0
  surface_field%sat_func_id = 0
  surface_field%perm_xx = 0
  surface_field%perm_yy = 0
  surface_field%perm_zz = 0
  surface_field%subsurf_xx = 0
  surface_field%subsurf_yy = 0
  surface_field%subsurf_zz = 0
  
  surface_field%subsurf_temp_vec = 0

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
  if (surface_field%mannings_loc /= 0) call VecDestroy(surface_field%mannings_loc,ierr)

  if (surface_field%work /= 0) call VecDestroy(surface_field%work,ierr)
  if (surface_field%work_loc  /= 0) call VecDestroy(surface_field%work_loc,ierr)

  if (surface_field%area  /= 0) call VecDestroy(surface_field%area,ierr)
  
  if (surface_field%vol_subsurf_2_surf /= 0) &
    call VecDestroy(surface_field%vol_subsurf_2_surf,ierr)
  if (surface_field%vol_subsurf_2_surf /= 0) &
    call VecDestroy(surface_field%press_subsurf,ierr)
  if (surface_field%press_subsurf /= 0) &
    call VecDestroy(surface_field%press_subsurf,ierr)

  if (surface_field%Dq /= 0) call VecDestroy(surface_field%Dq,ierr)

  if (surface_field%perm_xx/=0) call VecDestroy(surface_field%perm_xx,ierr)
  if (surface_field%perm_yy/=0) call VecDestroy(surface_field%perm_yy,ierr)
  if (surface_field%perm_zz/=0) call VecDestroy(surface_field%perm_zz,ierr)

  if (surface_field%por/=0) call VecDestroy(surface_field%por,ierr)
  if (surface_field%sat_func_id/=0) call VecDestroy(surface_field%sat_func_id,ierr)

  if (surface_field%subsurf_xx/=0) call VecDestroy(surface_field%subsurf_xx,ierr)
  if (surface_field%subsurf_yy/=0) call VecDestroy(surface_field%subsurf_yy,ierr)
  if (surface_field%subsurf_zz/=0) call VecDestroy(surface_field%subsurf_zz,ierr)
    

  if (surface_field%flow_r /= 0) call VecDestroy(surface_field%flow_r,ierr)
  if (surface_field%flow_xx /= 0) call VecDestroy(surface_field%flow_xx,ierr)
  if (surface_field%flow_xx_loc /= 0) call VecDestroy(surface_field%flow_xx_loc,ierr)
  if (surface_field%flow_dxx /= 0) call VecDestroy(surface_field%flow_dxx,ierr)
  if (surface_field%flow_yy /= 0) call VecDestroy(surface_field%flow_yy,ierr)
  if (surface_field%flow_accum /= 0) call VecDestroy(surface_field%flow_accum,ierr)
  
  if (surface_field%subsurf_temp_vec/=0) call VecDestroy(surface_field%subsurf_temp_vec,ierr)

end subroutine SurfaceFieldDestroy

end module Surface_Field_module

#endif