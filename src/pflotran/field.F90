module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  type, public :: field_type 
    
!geh material id
    PetscInt, pointer :: imat(:)
    
    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)

    ! 1 degree of freedom
    Vec :: porosity0, porosity_loc
    Vec :: tor_loc
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm_xx_loc, perm_yy_loc, perm_zz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    
    ! NDOF degree of freedom
    ! residual vector
    Vec :: r            
    ! Solution vectors
    Vec :: xx, xx_loc, dxx, yy, accum
   
  end type 

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !
!
! FieldCreate: Allocates and initializes a new Field object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function FieldCreate()

  implicit none
  
  type(field_type), pointer :: FieldCreate
  
  type(field_type), pointer :: field
  
  allocate(field)
  
  ! nullify PetscVecs
  field%porosity0 = 0
  field%porosity_loc = 0
  field%tor_loc = 0
  field%ithrm_loc = 0
  field%icap_loc = 0
  field%iphas_loc = 0
  field%iphas_old_loc = 0

  field%perm_xx_loc = 0
  field%perm_yy_loc = 0
  field%perm_zz_loc = 0
  field%perm0_xx = 0
  field%perm0_yy = 0
  field%perm0_zz = 0
  field%perm_pow = 0
  
  field%r = 0
  field%xx = 0
  field%xx_loc = 0
  field%dxx = 0
  field%yy = 0
  field%accum = 0
  
  nullify(field%imat)
  nullify(field%internal_velocities)
  nullify(field%boundary_velocities)
  
  FieldCreate => field
  
end function FieldCreate

! ************************************************************************** !
!
! FieldDestroy: Deallocates a field object
! author: Glenn Hammond
! date: 11/15/07
!
! ************************************************************************** !
subroutine FieldDestroy(field)

  implicit none
  
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr

  ! Destroy PetscVecs
  if (field%porosity0 /= 0) call VecDestroy(field%porosity0,ierr)
  if (field%porosity_loc /= 0) call VecDestroy(field%porosity_loc,ierr)
  if (field%tor_loc /= 0) call VecDestroy(field%tor_loc,ierr)
  if (field%ithrm_loc /= 0) call VecDestroy(field%ithrm_loc,ierr)
  if (field%icap_loc /= 0) call VecDestroy(field%icap_loc,ierr)
  if (field%iphas_loc /= 0) call VecDestroy(field%iphas_loc,ierr)
  if (field%iphas_old_loc /= 0) call VecDestroy(field%iphas_old_loc,ierr)

  if (field%perm_xx_loc /= 0) call VecDestroy(field%perm_xx_loc,ierr)
  if (field%perm_yy_loc /= 0) call VecDestroy(field%perm_yy_loc,ierr)
  if (field%perm_zz_loc /= 0) call VecDestroy(field%perm_zz_loc,ierr)
  if (field%perm0_xx /= 0) call VecDestroy(field%perm0_xx,ierr)
  if (field%perm0_yy /= 0) call VecDestroy(field%perm0_yy,ierr)
  if (field%perm0_zz /= 0) call VecDestroy(field%perm0_zz,ierr)
  if (field%perm_pow /= 0) call VecDestroy(field%perm_pow,ierr)
  
  if (field%r /= 0) call VecDestroy(field%r,ierr)
  if (field%xx /= 0) call VecDestroy(field%xx,ierr)
  if (field%xx_loc /= 0) call VecDestroy(field%xx_loc,ierr)
  if (field%dxx /= 0) call VecDestroy(field%dxx,ierr)
  if (field%yy /= 0) call VecDestroy(field%yy,ierr)
  if (field%accum /= 0) call VecDestroy(field%accum,ierr)
  
  if (associated(field%imat)) deallocate(field%imat)
  nullify(field%imat)
  if (associated(field%internal_velocities)) deallocate(field%internal_velocities)
  nullify(field%internal_velocities)
  if (associated(field%boundary_velocities)) deallocate(field%boundary_velocities)
  nullify(field%boundary_velocities)
  
    
  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
