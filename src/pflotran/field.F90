module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"
! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

  type, public :: field_type 
  
    real*8, pointer :: density_bc(:),d_p_bc(:),d_t_bc(:), d_s_bc(:),d_c_bc(:),&
                       avgmw_bc(:),avgmw_c_bc(:),&
                       hh_bc(:),h_p_bc(:),h_t_bc(:),h_s_bc(:), h_c_bc(:), &
                       viscosity_bc(:),v_p_bc(:),v_t_bc(:),&
                       uu_bc(:),u_p_bc(:),u_t_bc(:),u_s_bc(:), u_c_bc(:),&    
                       df_bc(:),df_p_bc(:),df_t_bc(:),df_s_bc(:), df_c_bc(:), &
                       hen_bc(:),hen_p_bc(:),hen_t_bc(:),hen_s_bc(:),hen_c_bc(:), &
                       pc_bc(:),pc_p_bc(:),pc_t_bc(:),pc_s_bc(:),pc_c_bc(:), &
                       kvr_bc(:),kvr_p_bc(:),kvr_t_bc(:),kvr_s_bc(:),kvr_c_bc(:)
    real*8, pointer :: xphi_co2(:),xxphi_co2(:),den_co2(:), dden_co2(:)
    
!geh material id
    integer, pointer :: imat(:)
    
    real*8, pointer :: internal_velocities(:,:)
    real*8, pointer :: boundary_velocities(:,:)

    real*8, pointer :: xphi_co2_bc(:), xxphi_co2_bc(:)

    !block BC values read from input
    real*8, pointer :: pressurebc0(:,:)
    real*8, pointer :: velocitybc0(:,:)

    !-------------------------------------------------------------------
    ! Quantities defined at each grid point.
    ! NOTE: I adopt the convention that _loc indicates the local portion
    ! of any global vector.
    !-------------------------------------------------------------------

    ! One degree of freedom: Physical coordinates.
    Vec :: porosity0, porosity_loc, tor_loc
    Vec :: ithrm_loc, icap_loc, iphas_loc, iphas_old_loc
    Vec :: phis

    Vec :: conc
    Vec :: ttemp, ttemp_loc, temp ! 1 dof

    ! Three degrees of freedom:
    Vec :: perm_xx_loc, perm_yy_loc, perm_zz_loc
    Vec :: perm0_xx, perm0_yy, perm0_zz, perm_pow
    ! Multiple degrees of freedom (equal to number of phases present):
    Vec :: var_loc
    Vec :: ppressure, ppressure_loc, pressure, dp
    Vec :: ssat, ssat_loc, sat    ! saturation
    Vec :: xxmol, xxmol_loc, xmol ! mole fraction
    Vec :: density       ! Density at time k
    Vec :: ddensity, ddensity_loc  ! Density at time k+1
    Vec :: d_p, d_p_loc  ! dD/dp at time k+1
    Vec :: d_t, d_t_loc  ! dD/dT at time k+1
    Vec :: d_c, d_c_loc  ! dD/dT at time k+1
    Vec :: d_s, d_s_loc  ! dD/dT at time k+1
    Vec :: avgmw,avgmw_loc  ! Density at time k+1molecular weight at time k+1
    Vec :: avgmw_c,avgmw_c_loc
    Vec :: h             ! H     at time k
    Vec :: hh, hh_loc    ! H     at time k+1
    Vec :: h_p, h_p_loc  ! dH/dp at time k+1
    Vec :: h_t, h_t_loc  ! dH/dT at time k+1
    Vec :: h_c, h_c_loc  ! dD/dT at time k+1
    Vec :: h_s, h_s_loc  ! dD/dT at time k+1
    Vec :: u            ! H     at time k
    Vec :: uu, uu_loc    ! H     at time k+1
    Vec :: u_p, u_p_loc  ! dH/dp at time k+1
    Vec :: u_t, u_t_loc  ! dH/dT at time k+1
    Vec :: u_c, u_c_loc  ! dD/dT at time k+1
    Vec :: u_s, u_s_loc  ! dD/dT at time k+1
    Vec :: hen, hen_loc    ! H     at time k+1
    Vec :: hen_p, hen_p_loc  ! dH/dp at time k+1
    Vec :: hen_t, hen_t_loc  ! dH/dT at time k+1
    Vec :: hen_c, hen_c_loc  ! dD/dT at time k+1
    Vec :: hen_s, hen_s_loc  ! dD/dT at time k+1
    Vec :: df, df_loc    ! H     at time k+1
    Vec :: df_p, df_p_loc  ! dH/dp at time k+1
    Vec :: df_t, df_t_loc  ! dH/dT at time k+1
    Vec :: df_c, df_c_loc  ! dD/dT at time k+1
    Vec :: df_s, df_s_loc  ! dD/dT at time k+1
    Vec :: viscosity, viscosity_loc  !kept for early routine
   
    Vec :: v_p, v_p_loc  ! dv/dp at time k+1
    Vec :: v_t, v_t_loc  ! dv/dT at time k+1
    Vec :: pcw, pcw_loc    ! H     at time k+1
    Vec :: pc_p, pc_p_loc  ! dH/dp at time k+1
    Vec :: pc_t, pc_t_loc  ! dH/dT at time k+1
    Vec :: pc_c, pc_c_loc  ! dD/dT at time k+1
    Vec :: pc_s, pc_s_loc  ! dD/dT at time k+1
    Vec :: kvr, kvr_loc    ! H     at time k+1
    Vec :: kvr_p, kvr_p_loc  ! d/dp at time k+1
    Vec :: kvr_t, kvr_t_loc  ! dm/dT at time k+1
    Vec :: kvr_c, kvr_c_loc  ! d/d at time k+1
    Vec :: kvr_s, kvr_s_loc  ! dD/dT at time k+1
    Vec :: r             ! The residual.  (NOT the negative of the residual.)

    Vec :: vl, vvl, vg, vvg ! phase (liquid and gas) velocities stored at interfaces
 
    real*8, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
    real*8, pointer :: vvlbc(:), vvgbc(:)
    real*8, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)

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

#include "definitions.h"
  
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
  field%phis = 0
  field%conc = 0
  field%ttemp = 0
  field%ttemp_loc = 0
  field%temp = 0
  field%perm_xx_loc = 0
  field%perm_yy_loc = 0
  field%perm_zz_loc = 0
  field%perm0_xx = 0
  field%perm0_yy = 0
  field%perm0_zz = 0
  field%perm_pow = 0
  
  field%ppressure = 0
  field%ppressure_loc = 0
  field%pressure = 0
  field%dp = 0
  field%ssat = 0
  field%ssat_loc = 0
  field%sat = 0
  field%xxmol = 0
  field%xmol = 0
  field%xxmol_loc = 0
  field%xmol = 0
  field%density = 0
  field%ddensity = 0
  field%ddensity_loc = 0
  field%d_p = 0
  field%d_p_loc = 0
  field%d_t = 0
  field%d_t_loc = 0
  field%d_c = 0
  field%d_c_loc = 0
  field%d_s = 0
  field%d_s_loc = 0
  field%avgmw = 0
  field%avgmw_loc = 0
  field%avgmw_c = 0
  field%avgmw_c_loc = 0
  field%h = 0
  field%hh = 0
  field%hh_loc = 0
  field%h_p = 0
  field%h_p_loc = 0
  field%h_t = 0
  field%h_t_loc = 0
  field%h_c = 0
  field%h_c_loc = 0
  field%h_s = 0
  field%h_s_loc = 0
  field%u = 0
  field%uu = 0
  field%uu_loc = 0
  field%u_p = 0
  field%u_p_loc = 0
  field%u_t = 0
  field%u_t_loc = 0
  field%u_c = 0
  field%u_c_loc = 0
  field%u_s = 0
  field%u_s_loc = 0
  field%hen = 0
  field%hen_loc = 0
  field%hen_p = 0
  field%hen_p_loc = 0
  field%hen_t = 0
  field%hen_t_loc = 0
  field%hen_c = 0
  field%hen_c_loc = 0
  field%hen_s = 0
  field%hen_s_loc = 0
  field%df = 0
  field%df_loc = 0
  field%df_s = 0
  field%df_s_loc = 0
  field%viscosity = 0
  field%viscosity_loc = 0
  field%v_p = 0
  field%v_p_loc = 0
  field%v_t = 0
  field%v_t_loc = 0
  field%pcw = 0
  field%pcw_loc = 0
  field%pc_p = 0
  field%pc_p_loc = 0
  field%pc_t = 0
  field%pc_t_loc = 0
  field%pc_c = 0
  field%pc_c_loc = 0
  field%pc_s = 0
  field%pc_s_loc = 0
  field%kvr = 0
  field%kvr_loc = 0
  field%kvr_p = 0
  field%kvr_p_loc = 0
  field%kvr_t = 0
  field%kvr_t_loc = 0
  field%kvr_c = 0
  field%kvr_c_loc = 0
  field%kvr_s = 0
  field%kvr_s_loc = 0
  
  field%vl = 0
  field%vvl = 0
  field%vg = 0
  field%vvg = 0
  
  field%r = 0
  field%xx = 0
  field%xx_loc = 0
  field%dxx = 0
  field%yy = 0
  field%accum = 0
  
  nullify(field%imat)
  nullify(field%density_bc)
  nullify(field%d_p_bc)
  nullify(field%d_t_bc)
  nullify(field%d_s_bc)
  nullify(field%d_c_bc)
  nullify(field%avgmw_bc)
  nullify(field%avgmw_c_bc)
  nullify(field%hh_bc)
  nullify(field%h_p_bc)
  nullify(field%h_t_bc)
  nullify(field%h_s_bc)
  nullify(field%h_c_bc)
  nullify(field%viscosity_bc)
  nullify(field%v_p_bc)
  nullify(field%v_t_bc)
  nullify(field%uu_bc)
  nullify(field%u_p_bc)   
  nullify(field%u_t_bc)   
  nullify(field%u_s_bc)   
  nullify(field%u_c_bc)   
  nullify(field%df_bc)   
  nullify(field%df_p_bc)   
  nullify(field%df_t_bc)   
  nullify(field%df_s_bc)   
  nullify(field%df_c_bc)   
  nullify(field%hen_bc)   
  nullify(field%hen_p_bc)   
  nullify(field%hen_t_bc)   
  nullify(field%hen_s_bc)   
  nullify(field%hen_c_bc)   
  nullify(field%pc_bc)   
  nullify(field%pc_p_bc)   
  nullify(field%pc_t_bc)   
  nullify(field%pc_s_bc)   
  nullify(field%pc_c_bc)   
  nullify(field%kvr_bc)   
  nullify(field%kvr_p_bc)   
  nullify(field%kvr_t_bc)   
  nullify(field%kvr_s_bc)   
  nullify(field%kvr_c_bc) 
    
  nullify(field%xphi_co2)   
  nullify(field%xxphi_co2)   
  nullify(field%den_co2)   
  nullify(field%dden_co2)   
    
!  nullify(field%pressurebc)   
!  nullify(field%velocitybc)   
!  nullify(field%tempbc)   
!  nullify(field%concbc)   
!  nullify(field%sgbc)   
  nullify(field%xphi_co2_bc)   
  nullify(field%xxphi_co2_bc)   
!  nullify(field%xxbc)   
!  nullify(field%varbc)   
!  nullify(field%iphasebc)   
  nullify(field%pressurebc0)   
  nullify(field%velocitybc0) 
  
   
  nullify(field%vl_loc)
  nullify(field%vvl_loc)
  nullify(field%vg_loc)
  nullify(field%vvg_loc)
  nullify(field%vvlbc)
  nullify(field%vvgbc)
  
  nullify(field%rtot)
  nullify(field%rate)
  nullify(field%area_var)
  nullify(field%delx)
  
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
  
  integer :: myrank, ierr
  
  call MPI_Comm_Rank(PETSC_COMM_WORLD,myrank,ierr)
  if (myrank == 0) then
    print *, 'Need to implement FieldDestroy'
  endif
  
  ! all kinds of stuff needs to be added here.
  
  
  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
