module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type, public :: field_type 
    
    !get material id
    ! 1 degree of freedom
    Vec :: porosity0
    Vec :: porosity_mnrl_loc
    Vec :: tortuosity0
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm0_xx, perm0_yy, perm0_zz
    Vec :: perm0_xz, perm0_xy, perm0_yz
    
    Vec :: work, work_loc

    Vec :: volume0
    
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

    ! vectors to save temporally average quantities
    Vec, pointer :: avg_vars_vec(:)
    PetscInt :: nvars

    ! vectors to save temporally average flowrates
    Vec :: flowrate_inst
    Vec :: flowrate_aveg

    ! vectors to save velocity at face
    Vec :: vx_face_inst
    Vec :: vy_face_inst
    Vec :: vz_face_inst

    Vec, pointer :: max_change_vecs(:)

  end type field_type

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !

function FieldCreate()
  ! 
  ! Allocates and initializes a new Field object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(field_type), pointer :: FieldCreate
  
  type(field_type), pointer :: field
  
  allocate(field)
  
  ! nullify PetscVecs
  field%porosity0 = 0
  field%porosity_mnrl_loc = 0
  field%tortuosity0 = 0
  field%ithrm_loc = 0
  field%icap_loc = 0
  field%iphas_loc = 0
  field%iphas_old_loc = 0

  field%perm0_xx = 0
  field%perm0_yy = 0
  field%perm0_zz = 0
  field%perm0_xz = 0
  field%perm0_xy = 0
  field%perm0_yz = 0
  
  field%work = 0
  field%work_loc = 0

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

  nullify(field%avg_vars_vec)
  field%nvars = 0

  field%flowrate_inst = 0
  field%flowrate_aveg = 0

  field%vx_face_inst = 0
  field%vy_face_inst = 0
  field%vz_face_inst = 0

  nullify(field%max_change_vecs)

  FieldCreate => field
  
end function FieldCreate

! ************************************************************************** !

subroutine FieldDestroy(field)
  ! 
  ! Deallocates a field object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/15/07
  ! 

  implicit none
  
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: ivar

  ! Destroy PetscVecs
  if (field%porosity0 /= 0) then
    call VecDestroy(field%porosity0,ierr)
  endif
  if (field%porosity_mnrl_loc /= 0) then
    call VecDestroy(field%porosity_mnrl_loc,ierr)
  endif
  if (field%tortuosity0 /= 0) then
    call VecDestroy(field%tortuosity0,ierr)
  endif
  if (field%ithrm_loc /= 0) then
    call VecDestroy(field%ithrm_loc,ierr)
  endif
  if (field%icap_loc /= 0) then
    call VecDestroy(field%icap_loc,ierr)
  endif
  if (field%iphas_loc /= 0) then
    call VecDestroy(field%iphas_loc,ierr)
  endif
  if (field%iphas_old_loc /= 0) then
    call VecDestroy(field%iphas_old_loc,ierr)
  endif

  if (field%perm0_xx /= 0) then
    call VecDestroy(field%perm0_xx,ierr)
  endif
  if (field%perm0_yy /= 0) then
    call VecDestroy(field%perm0_yy,ierr)
  endif
  if (field%perm0_zz /= 0) then
    call VecDestroy(field%perm0_zz,ierr)
  endif
  if (field%perm0_xz /= 0) then
    call VecDestroy(field%perm0_xz,ierr)
  endif
  if (field%perm0_xy /= 0) then
    call VecDestroy(field%perm0_xy,ierr)
  endif
  if (field%perm0_yz /= 0) then
    call VecDestroy(field%perm0_yz,ierr)
  endif
  
  if (field%work /= 0) then
    call VecDestroy(field%work,ierr)
  endif
  if (field%work_loc /= 0) then
    call VecDestroy(field%work_loc,ierr)
  endif

  if (field%volume0 /= 0) then
    call VecDestroy(field%volume0,ierr)
  endif
  
  if (field%flow_r /= 0) then
    call VecDestroy(field%flow_r,ierr)
  endif
  if (field%flow_xx /= 0) then
    call VecDestroy(field%flow_xx,ierr)
  endif
  if (field%flow_xx_loc /= 0) then
    call VecDestroy(field%flow_xx_loc,ierr)
  endif
  if (field%flow_dxx /= 0) then
    call VecDestroy(field%flow_dxx,ierr)
  endif
  if (field%flow_yy /= 0) then
    call VecDestroy(field%flow_yy,ierr)
  endif
  if (field%flow_accum /= 0) then
    call VecDestroy(field%flow_accum,ierr)
  endif
  
  if (field%tran_r /= 0) then
    call VecDestroy(field%tran_r,ierr)
  endif
  if (field%tran_log_xx /= 0) then
    call VecDestroy(field%tran_log_xx,ierr)
  endif
  if (field%tran_xx /= 0) then
    call VecDestroy(field%tran_xx,ierr)
  endif
  if (field%tran_xx_loc /= 0) then
    call VecDestroy(field%tran_xx_loc,ierr)
  endif
  if (field%tran_dxx /= 0) then
    call VecDestroy(field%tran_dxx,ierr)
  endif
  if (field%tran_yy /= 0) then
    call VecDestroy(field%tran_yy,ierr)
  endif
  if (field%tran_accum /= 0) then
    call VecDestroy(field%tran_accum,ierr)
  endif
  if (field%tran_work_loc /= 0) then
    call VecDestroy(field%tran_work_loc,ierr)
  endif
  
  if (field%tran_rhs /= 0) then
    call VecDestroy(field%tran_rhs,ierr)
  endif
  if (field%tran_rhs_coef /= 0) then
    call VecDestroy(field%tran_rhs_coef,ierr)
  endif

  if (field%flow_ts_mass_balance /= 0) then
    call VecDestroy(field%flow_ts_mass_balance,ierr)
  endif
  if (field%flow_total_mass_balance /= 0) then
    call VecDestroy(field%flow_total_mass_balance,ierr)
  endif
  if (field%tran_ts_mass_balance /= 0) then
    call VecDestroy(field%tran_ts_mass_balance,ierr)
  endif
  if (field%tran_total_mass_balance /= 0) then
    call VecDestroy(field%tran_total_mass_balance,ierr)
  endif
    
  if (field%tvd_ghosts /= 0) then
    call VecDestroy(field%tvd_ghosts,ierr)
  endif

  if (field%flow_r_faces/= 0) then
    call VecDestroy(field%flow_r_faces,ierr)
  endif

  if (field%flow_r_loc_faces /= 0) then
    call VecDestroy(field%flow_r_loc_faces ,ierr)
  endif

  if (field%flow_xx_faces /= 0) then
    call VecDestroy(field%flow_xx_faces ,ierr)
  endif

  if (field%flow_xx_loc_faces /= 0) then
    call VecDestroy(field%flow_xx_loc_faces ,ierr)
  endif

  if (field%flow_dxx_faces /= 0) then
    call VecDestroy(field%flow_dxx_faces ,ierr)
  endif

  if (field%flow_yy_faces /= 0) then
    call VecDestroy(field%flow_yy_faces ,ierr)
  endif

  if (field%flow_bc_loc_faces /= 0) then
    call VecDestroy(field%flow_bc_loc_faces ,ierr)
  endif

  if (field%work_loc_faces /= 0) then
    call VecDestroy(field%work_loc_faces ,ierr)
  endif

  do ivar = 1,field%nvars
    call VecDestroy(field%avg_vars_vec(ivar),ierr)
  enddo

  if (field%flowrate_inst/=0) then
    call VecDestroy(field%flowrate_inst,ierr)
  endif
  if (field%flowrate_aveg/=0) then
    call VecDestroy(field%flowrate_aveg,ierr)
  endif

  if (field%vx_face_inst/=0) then
    call VecDestroy(field%vx_face_inst,ierr)
  endif
  if (field%vy_face_inst/=0) then
    call VecDestroy(field%vy_face_inst,ierr)
  endif
  if (field%vz_face_inst/=0) then
    call VecDestroy(field%vz_face_inst,ierr)
  endif

  if (associated(field%max_change_vecs)) then
    call VecDestroyVecsF90(size(field%max_change_vecs), &
                           field%max_change_vecs,ierr)
  endif

  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
