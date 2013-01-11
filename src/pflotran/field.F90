module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: field_type 
    
    !get material id
    ! 1 degree of freedom
    Vec :: porosity0, porosity_loc
    Vec :: tortuosity0, tortuosity_loc
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm_xx_loc, perm_yy_loc, perm_zz_loc
    Vec :: perm_xz_loc, perm_xy_loc, perm_yz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    Vec :: perm0_xz, perm0_xy, perm0_yz
    
    Vec :: work, work_loc

    Vec :: volume, volume0
    
    ! residual vectors
    Vec :: flow_r          
    Vec :: tran_r
    
    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum
    Vec :: tran_xx, tran_xx_loc, tran_dxx, tran_yy, tran_accum

    ! vectors for operator splitting
    Vec :: tran_rhs
    Vec :: tran_rhs_coef
    
    Vec :: tran_log_xx, tran_work_loc
    
    Vec :: flow_ts_mass_balance, flow_total_mass_balance
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

    ! residual vectors for face unknows
    Vec :: flow_r_faces, flow_r_loc_faces          
      
    ! Solution vectors for face unknows
    Vec :: flow_xx_faces, flow_xx_loc_faces, flow_dxx_faces, flow_yy_faces, flow_bc_loc_faces
    Vec :: work_loc_faces   

    ! vector that holds the second layer of ghost cells for tvd
    Vec :: tvd_ghosts

  end type field_type

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
  field%tortuosity0 = 0
  field%tortuosity_loc = 0
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
  field%perm_xz_loc = 0
  field%perm_xy_loc = 0
  field%perm_yz_loc = 0
  field%perm0_xz = 0
  field%perm0_xy = 0
  field%perm0_yz = 0
  field%perm_pow = 0
  
  field%work = 0
  field%work_loc = 0

  field%volume = 0
  field%volume0 = 0
  
  field%flow_r = 0
  field%flow_xx = 0
  field%flow_xx_loc = 0
  field%flow_dxx = 0
  field%flow_yy = 0
  field%flow_accum = 0

  field%tran_r = 0
  field%tran_log_xx = 0
  field%tran_xx = 0
  field%tran_xx_loc = 0
  field%tran_dxx = 0
  field%tran_yy = 0
  field%tran_accum = 0
  field%tran_work_loc = 0

  field%tvd_ghosts = 0

  field%tran_rhs = 0
  field%tran_rhs_coef = 0
  
  field%flow_ts_mass_balance = 0
  field%flow_total_mass_balance = 0
  field%tran_ts_mass_balance = 0
  field%tran_total_mass_balance = 0

  field%flow_r_faces = 0
  field%flow_r_loc_faces = 0
  field%flow_xx_faces = 0
  field%flow_xx_loc_faces = 0
  field%flow_dxx_faces = 0
  field%flow_yy_faces = 0
  field%flow_bc_loc_faces = 0
  field%work_loc_faces = 0
   

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
  if (field%tortuosity0 /= 0) call VecDestroy(field%tortuosity0,ierr)
  if (field%tortuosity_loc /= 0) call VecDestroy(field%tortuosity_loc,ierr)
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
  if (field%perm_xz_loc /= 0) call VecDestroy(field%perm_xz_loc,ierr)
  if (field%perm_xy_loc /= 0) call VecDestroy(field%perm_xy_loc,ierr)
  if (field%perm_yz_loc /= 0) call VecDestroy(field%perm_yz_loc,ierr)
  if (field%perm0_xz /= 0) call VecDestroy(field%perm0_xz,ierr)
  if (field%perm0_xy /= 0) call VecDestroy(field%perm0_xy,ierr)
  if (field%perm0_yz /= 0) call VecDestroy(field%perm0_yz,ierr)
  if (field%perm_pow /= 0) call VecDestroy(field%perm_pow,ierr)
  
  if (field%work /= 0) call VecDestroy(field%work,ierr)
  if (field%work_loc /= 0) call VecDestroy(field%work_loc,ierr)

  if (field%volume /= 0) call VecDestroy(field%volume,ierr)
  if (field%volume0 /= 0) call VecDestroy(field%volume0,ierr)
  
  if (field%flow_r /= 0) call VecDestroy(field%flow_r,ierr)
  if (field%flow_xx /= 0) call VecDestroy(field%flow_xx,ierr)
  if (field%flow_xx_loc /= 0) call VecDestroy(field%flow_xx_loc,ierr)
  if (field%flow_dxx /= 0) call VecDestroy(field%flow_dxx,ierr)
  if (field%flow_yy /= 0) call VecDestroy(field%flow_yy,ierr)
  if (field%flow_accum /= 0) call VecDestroy(field%flow_accum,ierr)
  
  if (field%tran_r /= 0) call VecDestroy(field%tran_r,ierr)
  if (field%tran_log_xx /= 0) call VecDestroy(field%tran_log_xx,ierr)
  if (field%tran_xx /= 0) call VecDestroy(field%tran_xx,ierr)
  if (field%tran_xx_loc /= 0) call VecDestroy(field%tran_xx_loc,ierr)
  if (field%tran_dxx /= 0) call VecDestroy(field%tran_dxx,ierr)
  if (field%tran_yy /= 0) call VecDestroy(field%tran_yy,ierr)
  if (field%tran_accum /= 0) call VecDestroy(field%tran_accum,ierr)
  if (field%tran_work_loc /= 0) call VecDestroy(field%tran_work_loc,ierr)
  
  if (field%tran_rhs /= 0) call VecDestroy(field%tran_rhs,ierr)
  if (field%tran_rhs_coef /= 0) call VecDestroy(field%tran_rhs_coef,ierr)

  if (field%flow_ts_mass_balance /= 0) &
    call VecDestroy(field%flow_ts_mass_balance,ierr)
  if (field%flow_total_mass_balance /= 0) &
    call VecDestroy(field%flow_total_mass_balance,ierr)
  if (field%tran_ts_mass_balance /= 0) &
    call VecDestroy(field%tran_ts_mass_balance,ierr)
  if (field%tran_total_mass_balance /= 0) &
    call VecDestroy(field%tran_total_mass_balance,ierr)
    
  if (field%tvd_ghosts /= 0) &
    call VecDestroy(field%tvd_ghosts,ierr)

  if (field%flow_r_faces/= 0) &
    call VecDestroy(field%flow_r_faces,ierr)

  if (field%flow_r_loc_faces /= 0) &
    call VecDestroy(field%flow_r_loc_faces ,ierr)

  if (field%flow_xx_faces /= 0) &
    call VecDestroy(field%flow_xx_faces ,ierr)

  if (field%flow_xx_loc_faces /= 0) &
    call VecDestroy(field%flow_xx_loc_faces ,ierr)

  if (field%flow_dxx_faces /= 0) &
    call VecDestroy(field%flow_dxx_faces ,ierr)

  if (field%flow_yy_faces /= 0) &
    call VecDestroy(field%flow_yy_faces ,ierr)

  if (field%flow_bc_loc_faces /= 0) &
    call VecDestroy(field%flow_bc_loc_faces ,ierr)

  if (field%work_loc_faces /= 0) &
    call VecDestroy(field%work_loc_faces ,ierr)

  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
