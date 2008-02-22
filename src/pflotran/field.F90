module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"

  type, public :: field_type 
    
!geh material id
    ! 1 degree of freedom
    Vec :: porosity0, porosity_loc
    Vec :: tor_loc
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm_xx_loc, perm_yy_loc, perm_zz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    
    Vec :: saturation_loc
    
    ! residual vectors
    Vec :: flow_r          
    Vec :: tran_r          
      
    ! Solution vectors
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum
    Vec :: tran_xx, tran_xx_loc, tran_dxx, tran_yy, tran_accum
   
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
  
  field%saturation_loc = 0
  
  field%flow_r = 0
  field%flow_xx = 0
  field%flow_xx_loc = 0
  field%flow_dxx = 0
  field%flow_yy = 0
  field%flow_accum = 0
  
  field%tran_r = 0
  field%tran_xx = 0
  field%tran_xx_loc = 0
  field%tran_dxx = 0
  field%tran_yy = 0
  field%tran_accum = 0
  
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
  
  if (field%saturation_loc /= 0) call VecDestroy(field%saturation_loc,ierr)
  
  if (field%flow_r /= 0) call VecDestroy(field%flow_r,ierr)
  if (field%flow_xx /= 0) call VecDestroy(field%flow_xx,ierr)
  if (field%flow_xx_loc /= 0) call VecDestroy(field%flow_xx_loc,ierr)
  if (field%flow_dxx /= 0) call VecDestroy(field%flow_dxx,ierr)
  if (field%flow_yy /= 0) call VecDestroy(field%flow_yy,ierr)
  if (field%flow_accum /= 0) call VecDestroy(field%flow_accum,ierr)
  
  if (field%tran_r /= 0) call VecDestroy(field%tran_r,ierr)
  if (field%tran_xx /= 0) call VecDestroy(field%tran_xx,ierr)
  if (field%tran_xx_loc /= 0) call VecDestroy(field%tran_xx_loc,ierr)
  if (field%tran_dxx /= 0) call VecDestroy(field%tran_dxx,ierr)
  if (field%tran_yy /= 0) call VecDestroy(field%tran_yy,ierr)
  if (field%tran_accum /= 0) call VecDestroy(field%tran_accum,ierr)
    
  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
