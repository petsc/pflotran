!======================================================================

#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*option%nflowdof)
#define PPRESSURE(j,n)     xx_p(j+(n-1)*option%nflowdof)
#define PRESSURE(j,n)      yy_p(j+(n-1)*option%nflowdof)
#define TTEMP_LOC(n)       xx_loc_p(2+(n-1)*option%nflowdof)
#define TTEMP(n)           xx_p(2+(n-1)*option%nflowdof)
#define TEMP(n)            yy_p(2+(n-1)*option%nflowdof)
#define CCONC_LOC(n)       xx_loc_p(3+(n-1)*option%nflowdof)
#define CCONC(n)           xx_p(3+(n-1)*option%nflowdof)
#define CONC(n)            yy_p(3+(n-1)*option%nflowdof)

module thc_field_module
  
  implicit none
 
  private 
#include "definitions.h"
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
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
  
  type, public :: thc_field_type
    Vec :: density
    Vec :: ddensity, ddensity_loc
    Vec :: viscosity, viscosity_loc
    Vec :: d_p, d_p_loc
    Vec :: d_t, d_t_loc
    Vec :: h
    Vec :: hh, hh_loc
    Vec :: h_p, h_p_loc
    Vec :: h_t, h_t_loc
    Vec :: v_p, v_p_loc
    Vec :: v_t, v_t_loc
    Vec :: phis
 
    PetscReal, pointer :: vvl_loc(:)
    PetscReal, pointer :: d_p_bc(:)
    PetscReal, pointer :: d_t_bc(:)
    PetscReal, pointer :: hh_bc(:)
    PetscReal, pointer :: h_p_bc(:)
    PetscReal, pointer :: h_t_bc(:)
    PetscReal, pointer :: viscosity_bc(:)
    PetscReal, pointer :: v_p_bc(:)
    PetscReal, pointer :: v_t_bc(:)
    PetscReal, pointer :: density_bc(:)
    PetscReal, pointer :: vvlbc(:)
    
  end type thc_field_type
  
  type(thc_field_type), pointer, public :: thc_field
  
end module thc_field_module  

module thc_option_module
  
  implicit none
  
  private 
  
#include "definitions.h"
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
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
  
  type, public :: thc_option_type
    PetscInt :: jh2o = 1
    PetscReal, pointer :: rate(:)
    PetscReal, pointer :: area_var(:)
    
!   solid reaction rate
    PetscInt :: ityprxn
    PetscReal :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    PetscReal ::qu_kin, yh2o_in_co2=0.D0
        
  end type thc_option_type
  
  type(thc_option_type), pointer, public :: thc_option
  
end module thc_option_module  

module THC_module

 use Connection_module
 use Realization_module
 use Grid_module
 use Option_module
 use Coupler_module
 use Field_module
 
 use thc_field_module
 use thc_option_module
 
 private

#include "definitions.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "finclude/petscda.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petsclog.h"

  public THCResidual, THCJacobian, THCSetup, THCCheckpointRead, &
          THCCheckpointWrite, THCTimeCut, THCInitializeSolidReaction, &
          THCGetTecplotHeader, THCGetVarFromArray

contains

! ************************************************************************** !
!
! THCCheckpointWrite: Writes vecs to checkpoint file
! author: 
! date: 
!
! ************************************************************************** !
subroutine THCCheckpointWrite(grid,viewer)

  use Grid_module

  implicit none
  
  type(grid_type) grid
  PetscViewer :: viewer
  
  PetscErrorCode :: ierr
  
  call VecView(thc_field%hh, viewer, ierr)
  call VecView(thc_field%ddensity, viewer, ierr)
  ! solid volume fraction
  if (thc_option%rk > 0.d0) then
    call VecView(thc_field%phis, viewer, ierr)
  endif  
  
end subroutine THCCheckpointWrite

! ************************************************************************** !
!
! THCCheckpointRead: Reads vecs from checkpoint file
! author: 
! date: 
!
! ************************************************************************** !
subroutine THCCheckpointRead(grid,viewer)

  use Grid_module

  implicit none
  
  type(grid_type) grid
  PetscViewer :: viewer
  
  PetscErrorCode :: ierr
  
  call VecLoadIntoVector(viewer,thc_field%hh, ierr)
  call VecCopy(thc_field%hh,thc_field%h, ierr)
  call VecLoadIntoVector(viewer,thc_field%ddensity, ierr)
  call VecCopy(thc_field%ddensity,thc_field%density, ierr)
  ! solid volume fraction
  if (thc_option%rk > 0.d0) then
    call VecLoadIntoVector(viewer,thc_field%phis,ierr)
  endif 
  
end subroutine THCCheckpointRead

! ************************************************************************** !
!
! THCSetup: Creates arrays for auxilliary variables, initializes, etc.
! author: 
! date: 
!
! ************************************************************************** !
subroutine THCSetup(realization)

  use Realization_module
 
  implicit none
  
  type(realization_type) :: realization

  call THCInitializeSolidReaction(realization)
  
! these vecs need to be local to THC, not in the field type    
#if 0

    ! declarations formerly in field.F90
    PetscReal, pointer :: density_bc(:),d_p_bc(:),d_t_bc(:), d_s_bc(:),d_c_bc(:),&
                       avgmw_bc(:),avgmw_c_bc(:),&
                       hh_bc(:),h_p_bc(:),h_t_bc(:),h_s_bc(:), h_c_bc(:), &
                       viscosity_bc(:),v_p_bc(:),v_t_bc(:),&
                       uu_bc(:),u_p_bc(:),u_t_bc(:),u_s_bc(:), u_c_bc(:),&    
                       df_bc(:),df_p_bc(:),df_t_bc(:),df_s_bc(:), df_c_bc(:), &
                       hen_bc(:),hen_p_bc(:),hen_t_bc(:),hen_s_bc(:),hen_c_bc(:), &
                       pc_bc(:),pc_p_bc(:),pc_t_bc(:),pc_s_bc(:),pc_c_bc(:), &
                       kvr_bc(:),kvr_p_bc(:),kvr_t_bc(:),kvr_s_bc(:),kvr_c_bc(:)
                       
    Vec :: conc
    Vec :: ttemp, ttemp_loc, temp ! 1 dof
    
! someone will have to sort out whether these belong in thc, from field.F90
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
                           
    PetscReal, pointer :: vl_loc(:), vvl_loc(:), vg_loc(:), vvg_loc(:)
    PetscReal, pointer :: vvlbc(:), vvgbc(:)
    PetscReal, pointer :: rtot(:,:),rate(:),area_var(:), delx(:,:)

! from pflow_init
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)  
      call VecDuplicate(field%porosity_loc, field%ttemp_loc, ierr)
  end select
  
    ! should these be moved to their respective modules
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      ! nphase degrees of freedom
      call GridCreateVector(grid,NPHASEDOF,field%pressure,GLOBAL)
      call VecDuplicate(field%pressure, field%sat, ierr)
      call VecDuplicate(field%pressure, field%xmol, ierr)
      call VecDuplicate(field%pressure, field%ppressure, ierr)
      call VecDuplicate(field%pressure, field%ssat, ierr)
      call VecDuplicate(field%pressure, field%dp, ierr)
      call VecDuplicate(field%pressure, field%density, ierr)
      call VecDuplicate(field%pressure, field%ddensity, ierr)
      call VecDuplicate(field%pressure, field%avgmw, ierr)
      call VecDuplicate(field%pressure, field%d_p, ierr)
      call VecDuplicate(field%pressure, field%d_t, ierr)
      call VecDuplicate(field%pressure, field%h, ierr)
      call VecDuplicate(field%pressure, field%hh, ierr)
      call VecDuplicate(field%pressure, field%h_p, ierr)
      call VecDuplicate(field%pressure, field%h_t, ierr)
      call VecDuplicate(field%pressure, field%viscosity, ierr)
      call VecDuplicate(field%pressure, field%v_p, ierr)
      call VecDuplicate(field%pressure, field%v_t, ierr)
     ! xmol may not be nphase DOF, need change later 
      call VecDuplicate(field%pressure, field%flow_xxmol, ierr)
  end select

  ! should these be moved to their respective modules?
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      call GridCreateVector(grid,NPHASEDOF, field%ppressure_loc, LOCAL)
      call VecDuplicate(field%ppressure_loc, field%ssat_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%flow_xxmol_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%ddensity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%avgmw_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%d_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%hh_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%h_t_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%viscosity_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_p_loc, ierr)
      call VecDuplicate(field%ppressure_loc, field%v_t_loc, ierr)
  end select


  ! move to mph, maybe thc too
! Note: VecAssemblyBegin/End needed to run on the Mac - pcl (11/21/03)!
  if (field%conc /= 0) then
    call VecAssemblyBegin(field%conc,ierr)
    call VecAssemblyEnd(field%conc,ierr)
  endif
  if (field%xmol /= 0) then
    call VecAssemblyBegin(field%xmol,ierr)
    call VecAssemblyEnd(field%xmol,ierr)
  endif


                       
      allocate(field%density_bc(option%nphase))
      allocate(field%d_p_bc(option%nphase))
      allocate(field%d_t_bc(option%nphase))
      allocate(field%d_c_bc(option%nphase))
      allocate(field%d_s_bc(option%nphase))
      allocate(field%avgmw_bc(option%nphase))
      allocate(field%avgmw_c_bc(option%nphase*option%npricomp))
      allocate(field%hh_bc(option%nphase))
      allocate(field%h_p_bc(option%nphase))
      allocate(field%h_t_bc(option%nphase))
      allocate(field%h_c_bc(option%nphase*option%npricomp))
      allocate(field%h_s_bc(option%nphase))
      allocate(field%uu_bc(option%nphase))
      allocate(field%u_p_bc(option%nphase))
      allocate(field%u_t_bc(option%nphase))
      allocate(field%u_c_bc(option%nphase*option%npricomp))
      allocate(field%u_s_bc(option%nphase))
      allocate(field%df_bc(option%nphase*option%nspec))
      allocate(field%df_p_bc(option%nphase*option%nspec))
      allocate(field%df_t_bc(option%nphase*option%nspec))
      allocate(field%df_c_bc(option%nphase*option%nspec*option%npricomp))
      allocate(field%df_s_bc(option%nphase*option%nspec))
      allocate(field%hen_bc(option%nphase*option%nspec))
      allocate(field%hen_p_bc(option%nphase*option%nspec))
      allocate(field%hen_t_bc(option%nphase*option%nspec))
      allocate(field%hen_c_bc(option%nphase*option%nspec*option%npricomp))
      allocate(field%hen_s_bc(option%nphase*option%nspec))
      allocate(field%viscosity_bc(option%nphase))
      allocate(field%v_p_bc(option%nphase))
      allocate(field%v_t_bc(option%nphase))
      allocate(field%pc_bc(option%nphase))
      allocate(field%pc_p_bc(option%nphase))
      allocate(field%pc_t_bc(option%nphase))
      allocate(field%pc_c_bc(option%nphase*option%npricomp))
      allocate(field%pc_s_bc(option%nphase))
      allocate(field%kvr_bc(option%nphase))
      allocate(field%kvr_p_bc(option%nphase))
      allocate(field%kvr_t_bc(option%nphase))
      allocate(field%kvr_c_bc(option%nphase*option%npricomp))
      allocate(field%kvr_s_bc(option%nphase))
      
      temp_int = ConnectionGetNumberInList(grid%internal_connection_list)
      allocate(field%vl_loc(temp_int))
      allocate(field%vvl_loc(temp_int))
      allocate(field%vg_loc(temp_int))
      allocate(field%vvg_loc(temp_int))
      field%vl_loc = 0.D0
      field%vvl_loc = 0.D0
      field%vg_loc = 0.D0
      field%vvg_loc = 0.D0

  ! no longer supported in option, needs to be localized to thc  
  !set specific phase indices
  option%jh2o = 1; option%jgas =1
  select case(option%nphase)
    case(2)
      option%jco2 = 2
      option%jgas =3 
    case(3)
      option%jco2 = 2
      option%jgas =3 
  end select 
  
  ! these vecs need to be stored within this module, not in field
  select case(option%iflowmode)
    case(MPH_MODE,THC_MODE)
      call VecDuplicate(field%porosity0, field%conc, ierr)
      call VecDuplicate(field%porosity0, field%temp, ierr)
      call VecDuplicate(field%porosity0, field%ttemp, ierr)
  end select    
    
        
#endif

end subroutine THCSetup

! ************************************************************************** !
!
! THCTimeCut: Resets arrays for time step cut
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !
subroutine THCTimeCut(realization)
 
  use Realization_module
  use Option_module
  use Grid_module
  use Field_module
 
  implicit none
  
  type(realization_type) :: realization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  
  PetscReal, pointer :: xx_p(:),yy_p(:)
  PetscInt :: dof_offset,re
  PetscErrorCode :: ierr
  PetscInt :: local_id

  grid => realization%grid
  option => realization%option
  field => realization%field
 
#if 0
 ! this can't be write
  call VecCopy(thc_field%pressure, thc_field%ppressure, ierr)
  call VecCopy(thc_field%temp, thc_field%ttemp, ierr)
#endif  
 
end subroutine THCTimeCut

  subroutine THCResidual(snes, xx, r, grid)
  
  use water_eos_module
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Vec, intent(out) :: r
  type(realization_type) :: realization  ! What should 'intent' be for this?
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option

! external VISW, PSAT, WATEOS
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  PetscErrorCode :: ierr
  PetscInt :: n, ng, iconn, nr
  PetscInt :: i, i1, i2, j, jn, jng, jm, jm1, jm2, jmu
  PetscInt :: m, m1, m2, mu, n1, n2, ip1, ip2, p1, p2, t1, t2, c1, c2
  PetscInt :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:), &
               density_p(:), ddensity_p(:), ddensity_loc_p(:), phis_p(:), &
               viscosity_p(:), viscosity_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               vl_p(:), &
               d_p_p(:), d_p_loc_p(:), &
               d_t_p(:), d_t_loc_p(:), &
               h_p(:),  hh_p(:), hh_loc_p(:), &
               h_p_p(:), h_p_loc_p(:), &
               h_t_p(:), h_t_loc_p(:), &
               v_p_p(:), v_p_loc_p(:), &
               v_t_p(:), v_t_loc_p(:), &
               ithrm_loc_p(:)
  PetscReal :: dd1, dd2, diflux, diff, eeng, eng, cond, den, trans, density_ave, &
            hflx, fluxc, fluxe, fluxh, flux, fluxp, gravity, &
            fluxv, fluxbc, q, v_darcy, pvoldt, voldt, accum, eps
  PetscReal :: dd, f1, f2, ff, por1, por2, perm1, perm2
  PetscReal :: Dq, Dk  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: sat_pressure  ! Saturation pressure of water.
  PetscReal :: dw_kg, dw_mol
  PetscReal :: tsrc1, qsrc1, qqsrc, csrc1, hsrc1,enth_src
  PetscReal :: tempbc, concbc
  PetscReal, pointer :: pressurebc(:)
  
  logical*4 :: enthalpy_flag
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_list_type), pointer :: connection_list
  type(connection_type), pointer :: cur_connection_set
  PetscInt :: ibndtyp
  PetscReal :: upweight
  PetscReal :: distance, area, fraction_upwind
  PetscReal :: distance_gravity

  eps = 1.d-6

  grid => realization%grid
  option => realization%option
  field => realization%field
  
!  degrees of freedom:
!  p - 1:option%nphase
!  T - option%nphase + 1
!  C - option%nphase + 2
!  option%nflowdof = option%nphase + 2 = 3
  !---------------------------------------------------------------------------
  ! Calculate the density and viscosity of water at step k+1 at each local
  ! node.
  !---------------------------------------------------------------------------
 !  call recondition_bc(grid)
  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(thc_field%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(thc_field%hh, hh_p, ierr)
  call VecGetArrayF90(thc_field%viscosity, viscosity_p, ierr)

  if (option%ideriv == 1) then
    call VecGetArrayF90(thc_field%d_p, d_p_p, ierr)
    call VecGetArrayF90(thc_field%d_t, d_t_p, ierr)
    call VecGetArrayF90(thc_field%h_p, h_p_p, ierr)
    call VecGetArrayF90(thc_field%h_t, h_t_p, ierr)
    call VecGetArrayF90(thc_field%v_p, v_p_p, ierr)
    call VecGetArrayF90(thc_field%v_t, v_t_p, ierr)
  endif
  
  !get phase properties
  
  do j = 1, option%nphase
    if (j == thc_option%jh2o) then
      !pure liquid water eos
      do n = 1, grid%nlmax
        jn = j + (n-1)*option%nphase
      
!       first arg temp, 2nd pressure, use molar units for density 
!         and other properties

        call PSAT(TTEMP(n), sat_pressure, ierr)
        if (option%ideriv == 1) then
          call wateos(TTEMP(n),PPRESSURE(j,n),dw_kg,dw_mol,d_p_p(jn), &
          d_t_p(jn),hh_p(jn),h_p_p(jn),h_t_p(jn),option%scale,ierr)
          call VISW(TTEMP(n),PPRESSURE(j,n),sat_pressure, &
          viscosity_p(jn),v_t_p(jn),v_p_p(jn),ierr)
          
          !units: TTEMP [C], PPRESSURE [Pa]
          !       dw_mol [mol/dm^3], dw_kg [kg/m^3]
          !       hh_p [MJ/kmol], h_t_p [MJ/C/kmol]
          
!         print *,'THC: ',n,TTEMP(n),PPRESSURE(j,n),hh_p(jn),h_t_p(jn), &
!         dw_mol,dw_kg
        else
          call wateos_noderiv(TTEMP(n),PPRESSURE(j,n),dw_kg,dw_mol, &
          hh_p(jn),option%scale,ierr)
          call VISW_noderiv(TTEMP(n),PPRESSURE(j,n),sat_pressure, &
          viscosity_p(jn),ierr)
        endif
        ddensity_p(jn) = dw_mol ! mol/Liter
      enddo

!   add additional fluid phases here ...

!   else if (j == option%jgas) then
      !pure gas water + air eos
!     do n = 1, grid%nlmax
!       jn = j + (n-1)*option%nphase
!       call steameos(TTEMP(n),PPRESSURE(j,n),pa,dg_kg,dg_mol,d_p_p(jn), &
!       d_t_p(jn),hh_p(jn),h_p_p(jn),h_t_p(jn),option%scale,ierr)
!     enddo
      
!   else if (j == grid%jco2) then

!     do n = 1, option%nphase
!       ...
!     enddo
    endif
  enddo
  
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(thc_field%ddensity, ddensity_p, ierr)
  call VecRestoreArrayF90(thc_field%hh, hh_p, ierr)
  call VecRestoreArrayF90(thc_field%viscosity, viscosity_p, ierr)

  if (option%ideriv == 1) then  
    call VecRestoreArrayF90(thc_field%d_p, d_p_p, ierr)
    call VecRestoreArrayF90(thc_field%d_t, d_t_p, ierr)
    call VecRestoreArrayF90(thc_field%h_p, h_p_p, ierr)
    call VecRestoreArrayF90(thc_field%h_t, h_t_p, ierr)
    call VecRestoreArrayF90(thc_field%v_p, v_p_p, ierr)
    call VecRestoreArrayF90(thc_field%v_t, v_t_p, ierr)
  endif

  !---------------------------------------------------------------------------
  ! Now that we have calculated the density and viscosity for all local 
  ! nodes, we can perform the global-to-local scatters.
  !---------------------------------------------------------------------------
  call GridGlobalToLocal(grid, xx, field%flow_xx_loc, NFLOWDOF)

  call GridGlobalToLocal(grid, thc_field%ddensity, thc_field%ddensity_loc, &
                         NPHASEDOF)

  call GridGlobalToLocal(grid, thc_field%hh, thc_field%hh_loc, NPHASEDOF)

  call GridGlobalToLocal(grid, thc_field%viscosity, thc_field%viscosity_loc, &
                         NPHASEDOF)

  call GridLocalToLocal(grid, field%perm_xx_loc, field%perm_xx_loc, ONEDOF)
  call GridLocalToLocal(grid, field%perm_yy_loc, field%perm_yy_loc, ONEDOF)
  call GridLocalToLocal(grid, field%perm_zz_loc, field%perm_zz_loc, ONEDOF)

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%flow_yy, yy_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(thc_field%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(thc_field%density, density_p, ierr)
  call VecGetArrayF90(thc_field%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(thc_field%h, h_p, ierr)
  call VecGetArrayF90(thc_field%viscosity_loc, viscosity_loc_p, ierr)
  
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)

  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)

  if (option%ideriv == 1) then
    call GridGlobalToLocal(grid, thc_field%d_p, thc_field%d_p_loc, NPHASEDOF)

    call GridGlobalToLocal(grid, thc_field%d_t, thc_field%d_t, NPHASEDOF)

    call GridGlobalToLocal(grid, thc_field%h_p, thc_field%h_p_loc, NPHASEDOF)

    call GridGlobalToLocal(grid, thc_field%h_t, thc_field%h_t_loc, NPHASEDOF)

    call GridGlobalToLocal(grid, thc_field%v_p, thc_field%v_p_loc, NPHASEDOF)

    call GridGlobalToLocal(grid, thc_field%v_t, thc_field%v_t_loc, NPHASEDOF)

    call VecGetArrayF90(thc_field%d_p_loc, d_p_loc_p, ierr)
    call VecGetArrayF90(thc_field%d_t_loc, d_t_loc_p, ierr)
    call VecGetArrayF90(thc_field%h_p_loc, h_p_loc_p, ierr)
    call VecGetArrayF90(thc_field%h_t_loc, h_t_loc_p, ierr)
    call VecGetArrayF90(thc_field%v_p_loc, v_p_loc_p, ierr)
    call VecGetArrayF90(thc_field%v_t_loc, v_t_loc_p, ierr)
  endif

  if (thc_option%rk > 0.d0) then
    call VecGetArrayF90(thc_field%phis,phis_p,ierr)
  endif

  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------

  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    p1 = 1 + (n-1)*option%nflowdof
    t1 = p1 + 1
    c1 = t1 + 1
    voldt = volume_p(n) / option%flow_dt
    pvoldt = porosity_loc_p(ng) * voldt
    accum = 0.d0
    do j = 1, option%nphase
      jn = j + (n-1)*option%nphase
      jng = j + (ng-1)*option%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
      accum = accum + pvoldt * (ddensity_loc_p(jng) - density_p(jn))
    enddo
    
!   pressure residual
    r_p(p1) = accum
    
!   heat residual
    i = ithrm_loc_p(ng)
    j = thc_option%jh2o
    jn = j+(n-1)*option%nphase
    jng = j+(ng-1)*option%nphase
    ! rho U = rho H - p
    eeng = ddensity_loc_p(jng)*hh_loc_p(jng)-option%scale*PPRESSURE_LOC(j,ng)
    eng = density_p(jn) * h_p(jn) - option%scale*PRESSURE(j,n)
    r_p(t1) = pvoldt * (eeng - eng) &
                + (1.d0-porosity_loc_p(ng)) * option%dencpr(i) * &
                (TTEMP_LOC(ng) - TEMP(n)) * voldt

    !tracer
    r_p(c1) = pvoldt * (CCONC_LOC(ng) - CONC(n)) * option%ret
    
    !kinetic rate term
    if (thc_option%rk > 0.d0) then
      thc_option%rate(n) = 0.d0
      if (thc_option%ityprxn == 1) then
        if (phis_p(n) > 0.d0 .or. CCONC_LOC(ng)/thc_option%ceq > 1.d0+eps) then
          thc_option%rate(n) = -thc_option%rk * thc_option%area_var(n) * &
          (1.d0 - CCONC_LOC(ng)/thc_option%ceq)
          r_p(c1) = r_p(c1) + thc_option%rate(n) * volume_p(n)
          r_p(t1) = r_p(t1) - thc_option%delHs * thc_option%rate(n) * &
          volume_p(n)
!         print *,'resTHC: ',n,ng,phis_p(n),option%rate(n),option%area_var(n), &
!         CCONC_LOC(ng),CONC(n)
        endif
      else if (thc_option%ityprxn == 2) then
        if (phis_p(n) > 0.d0 .or. TTEMP_LOC(ng) < 0.d0+eps) then
          thc_option%rate(n) = -thc_option%rk * thc_option%area_var(n) * TTEMP_LOC(ng)
          r_p(p1) = r_p(p1) + thc_option%rate(n) * volume_p(n)
          r_p(t1) = r_p(t1) - thc_option%delHs * thc_option%rate(n) * &
          volume_p(n)
!         print *,'resTHC: ',n,ng,phis_p(n),option%rate(n),option%area_var(n), &
!         CCONC_LOC(ng),CONC(n)
        endif
      endif
    endif
    
  enddo

!---------------------------------------------------------------------------
! Flux terms for interior nodes
!---------------------------------------------------------------------------

  connection_list => grid%internal_connection_list
  cur_connection_set => connection_list%first

  ! Note: Currently, grid%internal_connection_list (which is a list of 
  ! *sets* of connections) only contains one entry: the set of connections 
  ! for the primary continuum.  Therefore I am NOT looping over the 
  ! elements of grid%internal_connection_list for now.  --RTM

  do iconn = 1, cur_connection_set%num_connections !For each interior connection...
    m1 = cur_connection_set%id_up(iconn)  ! m1, m2 are ghosted indices,
    m2 = cur_connection_set%id_dn(iconn)  ! indicating upstream/downstream nodes

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2) ! Ghost to local mapping
    
    p1 = 1 + (n1-1)*option%nflowdof
    p2 = 1 + (n2-1)*option%nflowdof
    t1 = p1 + 1
    t2 = p2 + 1
    c1 = t1 + 1
    c2 = t2 + 1
    
    fraction_upwind = cur_connection_set%dist(-1,iconn)
    distance = cur_connection_set%dist(0,iconn)
    area = cur_connection_set%area(iconn)
    ! The gravity vector can point in any direction, so we must compute the 
    ! dot product of gravity with the direction vector between grid cells to 
    ! get the 'distance_gravity' quantity.
    distance_gravity = distance * OptionDotProduct(option%gravity, &
                       cur_connection_set%dist(1:3,iconn))
    dd1 = distance*fraction_upwind
    dd2 = distance-dd1 ! should avoid truncation error
    upweight = dd2/(dd1+dd2)

    ! for now, just assume diagonal tensor
    perm1 = perm_xx_loc_p(m1)*abs(cur_connection_set%dist(1,iconn))+ &
            perm_yy_loc_p(m1)*abs(cur_connection_set%dist(2,iconn))+ &
            perm_zz_loc_p(m1)*abs(cur_connection_set%dist(3,iconn))

    perm2 = perm_xx_loc_p(m2)*abs(cur_connection_set%dist(1,iconn))+ &
            perm_yy_loc_p(m2)*abs(cur_connection_set%dist(2,iconn))+ &
            perm_zz_loc_p(m2)*abs(cur_connection_set%dist(3,iconn))
    
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    fluxp = 0.d0
    fluxh = 0.d0
    fluxv = 0.d0

    do j = 1, option%nphase
      jm1 = j + (m1-1) * option%nphase
      jm2 = j + (m2-1) * option%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

!     D1 = perm1 / viscosity_loc_p(jm1)
!     D2 = perm2 / viscosity_loc_p(jm2)
!     Dq = (D1 * D2) / (dd2*D1 + dd1*D2)

      D1 = perm1 * viscosity_loc_p(jm2)
      D2 = perm2 * viscosity_loc_p(jm1)
      den = dd2*D1 + dd1*D2
      Dq = (perm1 * perm2) / den

      density_ave = f2 * ddensity_loc_p(jm1) + f1 * ddensity_loc_p(jm2)

      gravity = option%fmwh2o * distance_gravity 

      v_darcy = -Dq * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) & 
                - gravity * density_ave)

      ! store velocities defined at interfaces in PETSc Vec vl at upstream node
      thc_field%vvl_loc(iconn) = v_darcy     ! use for coupling to ptran
      if (n1 > 0) then               ! If the upstream node is not a ghost node...
        vl_p(ip1+3*(n1-1)) = v_darcy ! use for print out of velocity
      endif
      
!     if (grid%t/grid%tconv >= 0.1-1.e-6) &
!     print *,'pflowTHC: ',iconn,n1,n2,v_darcy*365.*24.*60.*60.

      q = v_darcy * area

      flux = density_ave * q

      fluxp = fluxp + flux
      
      !upstream weighting
      if (q > 0.d0) then
        mu = m1
      else
        mu = m2
      endif

      jmu = j + (mu-1) * option%nphase
      fluxh = fluxh + q * density_ave * hh_loc_p(jmu)
      fluxv = fluxv + q * CCONC_LOC(mu)
    enddo

    !heat & tracer residual
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = option%ckwet(i1)
    D2 = option%ckwet(i2)

    Dk = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = Dk * area
    hflx = cond * (TTEMP_LOC(m2) - TTEMP_LOC(m1))
    fluxe = fluxh - hflx

    por1 = porosity_loc_p(m1)
    por2 = porosity_loc_p(m2)
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * option%difaq
    diflux = diff * area * (CCONC_LOC(m2) - CCONC_LOC(m1))
    fluxc = option%fc * fluxv - diflux

    ! interface on rhs
    if (n1 > 0) then  ! If upstream node is not a ghost node...
      r_p(p1) = r_p(p1) + fluxp
      r_p(t1) = r_p(t1) + fluxe         ! heat
      r_p(c1) = r_p(c1) + fluxc         ! tracer
    endif

    ! interface on lhs
    if (n2 > 0) then ! If downstream node is not a ghost node...
      r_p(p2) = r_p(p2) - fluxp
      r_p(t2) = r_p(t2) - fluxe
      r_p(c2) = r_p(c2) - fluxc
    endif
  enddo  ! End loop over interior connections.

!---------------------------------------------------------------------------
! Flux terms for boundary nodes.
!---------------------------------------------------------------------------

  boundary_condition => realization%flow_boundary_conditions%first

  do
    if (.not. associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection

    do iconn = 1, cur_connection_set%num_connections
      
      area = cur_connection_set%area(iconn)
      distance = cur_connection_set%dist(0,iconn)
      distance_gravity = distance * OptionDotProduct(option%gravity, &
                         cur_connection_set%dist(1:3,iconn))
      gravity = option%fmwh2o * distance_gravity

      m = cur_connection_set%id_dn(iconn) ! Note that here, m is NOT ghosted.
      ng = grid%nL2G(m)

      ! this must be zeroed since field%vvlbc references it.  
      v_darcy = 0.d0
      
      p1 = 1 + (m-1) * option%nflowdof
      t1 = p1 + 1
      c1 = t1 + 1

      ! for now, just assume diagonal tensor
      perm1 = perm_xx_loc_p(ng)*abs(cur_connection_set%dist(1,iconn))+ &
              perm_yy_loc_p(ng)*abs(cur_connection_set%dist(2,iconn))+ &
              perm_zz_loc_p(ng)*abs(cur_connection_set%dist(3,iconn))
      ! The below assumes a unit gravity vector of [0,0,1]
      distance_gravity = abs(cur_connection_set%dist(3,iconn)) * &
                         cur_connection_set%dist(0,iconn)
      
      ibndtyp = boundary_condition%condition%itype(1)

      if(ibndtyp == 2) then
        ! solve for pb from Darcy's law given qb /= 0
        boundary_condition%aux_real_var(THC_TEMPERATURE_DOF, iconn) = TTEMP_LOC(ng)
      else if(ibndtyp == 3) then
        boundary_condition%aux_real_var(THC_TEMPERATURE_DOF, iconn) = TTEMP_LOC(ng)
      endif

      tempbc = boundary_condition%aux_real_var(THC_TEMPERATURE_DOF, iconn)
      pressurebc => boundary_condition%aux_real_var(THC_PRESSURE_DOF:option%nphase, iconn)
      concbc = boundary_condition%aux_real_var(THC_CONCENTRATION_DOF, iconn)

      do j = 1, option%nphase
        if (j == thc_option%jh2o) then
          !pure water eos
          call PSAT(tempbc,sat_pressure,ierr)
          if (option%ideriv == 1) then
            call wateos(tempbc,pressurebc(j),dw_kg, &
            dw_mol,thc_field%d_p_bc(j),thc_field%d_t_bc(j),thc_field%hh_bc(j), &
            thc_field%h_p_bc(j),thc_field%h_t_bc(j),option%scale,ierr)
            
            call VISW(tempbc,pressurebc(j), &
            sat_pressure,thc_field%viscosity_bc(j), &
            thc_field%v_t_bc(j),thc_field%v_p_bc(j),ierr)

  !         print *,'thc-bc: ',iconn,j,m,ng,ibc,option%nphase, &
  !         option%tempbc(iconn),option%pressurebc(j),dw_mol, &
  !         dw_kg,field%d_p_bc(j),field%d_t_bc(j),field%hh_bc(j), &
  !         field%h_p_bc(j),field%h_t_bc(j),sat_pressure, &
  !         field%viscosity_bc(j),field%v_t_bc(j),field%v_p_bc(j)
          else
            call wateos_noderiv(tempbc,pressurebc(j), &
            dw_kg,dw_mol,thc_field%hh_bc(j),option%scale,ierr)
            
            call VISW_noderiv(tempbc,pressurebc(j), &
            sat_pressure,thc_field%viscosity_bc(j),ierr)
          endif
          thc_field%density_bc(j) = dw_mol

  !     else if (j == grid%jco2) then

  !       add additional fluid phases here ...

  !       do n = 1, option%nphase
  !         ...
  !       enddo
        else ! for testing purposes only
          jng = j + (ng-1) * option%nphase
          thc_field%viscosity_bc(j) = viscosity_loc_p(jng)
          thc_field%hh_bc(j) = hh_loc_p(jng)
          thc_field%density_bc(j) = ddensity_loc_p(jng)
        endif
      enddo

      if(ibndtyp == 1) then  ! Dirichlet BC for p, T, C ...
      
        fluxp = 0.d0
        fluxh = 0.d0
        fluxc = 0.d0
        do j = 1, option%nphase
          jm = j + (m-1) * option%nphase
          jng = j + (ng-1) * option%nphase

          Dq = perm1 / thc_field%viscosity_bc(j) / distance 
          
          !note: darcy vel. is positive for flow INTO boundary node
          v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - pressurebc(j) &
                    - gravity * thc_field%density_bc(j))
          q = v_darcy * area
        
  !       print *,'pflowTHC-1: ',iconn,m,ng,ibc,v_darcy,ibndtyp

          fluxbc = q * thc_field%density_bc(j)
          fluxp = fluxp - fluxbc

          !upstream weighting
          if (q > 0.d0) then
            fluxh = fluxh + q * thc_field%density_bc(j) * thc_field%hh_bc(j)
            fluxc = fluxc + q * concbc  ! note: need to change 
          else                          ! definition of C for multiple phases
            fluxh = fluxh + q * thc_field%density_bc(j) * hh_loc_p(jng)
            fluxc = fluxc + q * CCONC_LOC(ng) 
          endif
        enddo
          
        r_p(p1) = r_p(p1) + fluxp

        !heat residual
        i1 = ithrm_loc_p(ng)
        cond = option%ckwet(i1) * area / distance
        r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - tempbc) - fluxh
        
        !tracer
        trans = porosity_loc_p(ng) * option%difaq * area / distance
        r_p(c1) = r_p(c1) + trans * (CCONC_LOC(ng) - concbc) &
                  - option%fc * fluxc
                  
  !     print *,'THC: ',iconn,ng,ibc,q,trans,concbc,CCONC_LOC(ng)

      else if(ibndtyp == 2) then ! Constant velocity q, grad T,C = 0

        fluxp = 0.d0
        do j = 1, option%nphase
          jm = j + (m-1) * option%nphase
          jng = j + (ng-1) * option%nphase
          
          v_darcy = pressurebc(j)
          
          q = v_darcy * thc_field%density_bc(j) * area
          fluxp = fluxp - q
        enddo

        r_p(p1) = r_p(p1) + fluxp

        !heat residual: specified temperature
  !     i1 = ithrm_loc_p(ng)
  !     cond = option%ckwet(i1) * area / distance
  !     r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - option%tempbc(ibc)) &
  !                           + fluxh
        
        !tracer: specified concentration
  !     trans = option%difaq * area / distance
        !check for upstreaming weighting and iface even or odd etc.
        !use zero gradient BC
        
  !     r_p(c1) = r_p(c1) + trans * (CCONC_LOC(ng) - concbc) &
  !                               + fluxc

      else if(ibndtyp == 3) then  ! fixed p, grad T, C = 0
      
        fluxp = 0.d0
        fluxh = 0.d0
        fluxc = 0.d0
        do j = 1, option%nphase
          jm = j + (m-1) * option%nphase
          jng = j + (ng-1) * option%nphase
          
          Dq = perm1 / viscosity_loc_p(jng) / distance
          
          !v_darcy is positive for fluid flowing into block
          v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - pressurebc(j) &
                    - gravity * ddensity_loc_p(jng))
          q = v_darcy * area
        
  !       print *,'pflowTHC-3: ',iconn,m,ng,ibc,v_darcy,option%ibndtyp
          
          fluxbc = q * ddensity_loc_p(jng)
          fluxp = fluxp - fluxbc
          
          !upstream weighting not needed: c_N = c_b
          fluxh = fluxh + q * ddensity_loc_p(jng) * hh_loc_p(jng)
          fluxc = fluxc + q * CCONC_LOC(ng)

          !upstream weighting
  !       if (q > 0.d0) then
  !         fluxh = fluxh + q * field%density_bc(j) * field%hh_bc(j)
  !         fluxc = fluxc + q * concbc       ! note: need to change 
  !       else                          ! definition of C for multiple phases
  !         fluxh = fluxh + q * field%density_bc(j) * hh_loc_p(jng)
  !         fluxc = fluxc + q * CCONC_LOC(ng) 
  !       endif
        enddo

        !mass (pressure)
        r_p(p1) = r_p(p1) + fluxp

        !heat
        r_p(t1) = r_p(t1) - fluxh

        !solute (tracer)
        r_p(c1) = r_p(c1) - option%fc * fluxc

      else if(ibndtyp == 4) then ! grad p = 0, fixed T, grad C = 0

       ! fluxp = 0.d0
       ! do j = 1, option%nphase
       !   jm = j + (m-1) * option%nphase
       !   jng = j + (ng-1) * option%nphase
          
       !   v_darcy = option%velocitybc(j,iconn)
       !   
       !   q = v_darcy * field%density_bc(j) * area
       !   fluxp = fluxp - q
       ! enddo

        ! r_p(p1) = r_p(p1) + fluxp
        
        !heat residual: specified temperature
        i1 = ithrm_loc_p(ng)
        cond = option%ckwet(i1) * area / distance
        r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - tempbc) 
                            ! + fluxh
  !     print *, 'thc:res:bc4: ',   
        !tracer: specified concentration
  !     trans = option%difaq * area / distance
        !check for upstreaming weighting and iface even or odd etc.
        !use zero gradient BC
        
  !     r_p(c1) = r_p(c1) + trans * (CCONC_LOC(ng) - concbc) &
  !                               + fluxc
        
      endif

      thc_field%vvlbc(iconn) = v_darcy
    enddo

    boundary_condition => boundary_condition%next

  enddo  ! End loop over boundary connection linked list.

!---------------------------------------------------------------------------
! add source/sink terms
!---------------------------------------------------------------------------

  source_sink => realization%flow_source_sinks%first
  do
    if(.not. associated(source_sink)) exit

#if 0      
    kk1 = grid%k1src(nr) - grid%nzs
    kk2 = grid%k2src(nr) - grid%nzs
    jj1 = grid%j1src(nr) - grid%nys
    jj2 = grid%j2src(nr) - grid%nys
    ii1 = grid%i1src(nr) - grid%nxs
    ii2 = grid%i2src(nr) - grid%nxs
        
    kk1 = max(1,kk1)
    kk2 = min(grid%nlz,kk2)
    jj1 = max(1,jj1)
    jj2 = min(grid%nly,jj2)
    ii1 = max(1,ii1)
    ii2 = min(grid%nlx,ii2)
        
    if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
    do i = 2, grid%ntimsrc
      if (grid%timesrc(i,nr) == grid%t) then
        tsrc1 = grid%tempsrc(i,nr)
        qsrc1 = grid%qsrc(i,nr)
        csrc1 = grid%csrc(i,nr)
        hsrc1 = grid%hsrc(i,nr)
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1 = f1*grid%qsrc(i,nr) + f2*grid%qsrc(i-1,nr)
        csrc1 = f1*grid%csrc(i,nr) + f2*grid%csrc(i-1,nr)
        hsrc1 = f1*grid%hsrc(i,nr) + f2*grid%hsrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
!   print *,'pflowTHC: ', grid%myrank,i,grid%timesrc(i,nr), &
!   grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1

#endif

    ! check whether enthalpy dof is included
    if (source_sink%condition%num_sub_conditions > THC_CONCENTRATION_DOF) then
      enthalpy_flag = .true.
    else
      enthalpy_flag = .false.
    endif

    qsrc1 = source_sink%condition%pressure%dataset%cur_value(1)
    tsrc1 = source_sink%condition%temperature%dataset%cur_value(1)
    csrc1 = source_sink%condition%concentration%dataset%cur_value(1)
    if (enthalpy_flag) hsrc1 = source_sink%condition%enthalpy%dataset%cur_value(1)

    qsrc1 = qsrc1 / option%fmwh2o  ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
    csrc1 = csrc1 / option%fmwco2

    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      p1 = (local_id-1)*option%nflowdof + 1
      t1 = p1 + 1
      c1 = t1 + 1

      if (enthalpy_flag) then
        r_p(t1) = r_p(t1) - hsrc1
      endif

      if (qsrc1 > 0.d0) then ! injection
        call wateos_noderiv(tsrc1,PPRESSURE_LOC(thc_option%jh2o,ghosted_id), &
                            dw_kg,dw_mol,enth_src,option%scale,ierr)

        qqsrc = qsrc1/dw_mol  ! Do we still need to divide this? 
                              ! Or has this changed in the overhaul? --RTM
            
        r_p(p1) = r_p(p1) - qsrc1
        r_p(t1) = r_p(t1) - qsrc1*enth_src
        r_p(c1) = r_p(c1) - qqsrc*csrc1
      
      else if (qsrc1 < 0.d0) then ! withdrawal
        qqsrc = qsrc1/ddensity_loc_p(ng)
        enth_src = hh_loc_p(ng)
                
        r_p(p1) = r_p(p1) - qsrc1
        r_p(t1) = r_p(t1) - qsrc1*enth_src
        r_p(c1) = r_p(c1) - qqsrc*CCONC_LOC(ng)
      end if

    enddo ! End loop over connections in cur_connection_set

    source_sink => source_sink%next
  enddo ! End loop over source-sinks linked list.

  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%density, density_p, ierr)
  call VecRestoreArrayF90(thc_field%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%h, h_p, ierr)
  call VecRestoreArrayF90(thc_field%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)

  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)

  if (option%ideriv == 1) then
    call VecRestoreArrayF90(thc_field%d_p_loc, d_p_loc_p, ierr)
    call VecRestoreArrayF90(thc_field%d_t_loc, d_t_loc_p, ierr)
    call VecRestoreArrayF90(thc_field%h_p_loc, h_p_loc_p, ierr)
    call VecRestoreArrayF90(thc_field%h_t_loc, h_t_loc_p, ierr)
    call VecRestoreArrayF90(thc_field%v_p_loc, v_p_loc_p, ierr)
    call VecRestoreArrayF90(thc_field%v_t_loc, v_t_loc_p, ierr)
  endif
  
  if (thc_option%rk > 0.d0) then
    call VecRestoreArrayF90(thc_field%phis,phis_p,ierr)
  endif

! call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

  end subroutine THCResidual

!======================================================================

  subroutine THCJacobian(snes, xx, A, B, flag, grid, ierr)

  use water_eos_module

  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(grid_type), intent(inout) :: grid
  PetscInt, intent(out) :: flag

! external WATEOS, VISW, PSAT
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  PetscErrorCode :: ierr
  PetscInt :: n, ng, iconn
  PetscInt :: i, i1, i2, j, jn, jng, jm, jm1, jm2
  PetscInt :: m, m1, m2, mu, n1, n2, ip1, ip2
  PetscInt :: p1,p2,t1,t2,c1,c2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
  PetscReal :: elempp, elempt, elem1, elem2, v_darcy, q, qp1, qp2, qt1, qt2, eps
  PetscReal, pointer :: porosity_loc_p(:),volume_p(:),xx_loc_p(:), &
               ddensity_loc_p(:), phis_p(:), &
               d_p_loc_p(:),d_t_loc_p(:), &
               hh_loc_p(:),h_p_loc_p(:),h_t_loc_p(:), &
               viscosity_loc_p(:),v_p_loc_p(:),v_t_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               ithrm_loc_p(:)
  PetscReal :: cond, trans, trans1, trans2, gravity, &
            density_ave, dw_mol, dw_kg, den, voldt, pvoldt
  PetscReal :: daccep, daccet, dtrans, dd1, dd2, dd, f1, f2, u1, u2, qdt1, qdt2
  PetscReal :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2, cupstrm
  PetscReal :: por1, por2, perm1, perm2, diff, hm1, hm2
  PetscReal :: Dk, Dq           ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2, Dk1, Dk2 ! "Diffusion" constants upstream and downstream from a face.
  PetscReal :: sat_pressure  ! Saturation pressure of water.
  
  PetscReal :: blkmat1(3,3),blkmat2(3,3)
  
#ifdef COMPILE_BROKEN
  eps = 1.d-6

  flag = SAME_NONZERO_PATTERN

  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
! Not sure, but snes/examples/tutorials/ex18.c does it.  --RTM
  call GridGlobalToLocal(grid, xx, field%flow_xx_loc, NDOF)
                          
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(thc_field%ddensity_loc, ddensity_loc_p, ierr)
! call VecGetArrayF90(field%density, density_p, ierr)
  call VecGetArrayF90(thc_field%viscosity_loc, viscosity_loc_p, ierr)
  
  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(thc_field%d_p_loc, d_p_loc_p, ierr)
  call VecGetArrayF90(thc_field%d_t_loc, d_t_loc_p, ierr)
  call VecGetArrayF90(thc_field%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(thc_field%h_p_loc, h_p_loc_p, ierr)
  call VecGetArrayF90(thc_field%h_t_loc, h_t_loc_p, ierr)
  call VecGetArrayF90(thc_field%v_p_loc, v_p_loc_p, ierr)
  call VecGetArrayF90(thc_field%v_t_loc, v_t_loc_p, ierr)
  
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  
  if (option%rk > 0.d0) then
    call VecGetArrayF90(field%phis,phis_p,ierr)
  endif
  
  !---------------------------------------------------------------------------
  ! Calculate jacobian for accumulation term for local nodes.
  ! Structure:
  !  do j = 1, ndof
  !    do l = 1, ndof
  !      J(j,l) = ... (p,p) <diagonal>
  !    enddo
  !    J(j, Np+1) = ... (p,T)
  !  enddo
  !  do l = 1, ndof
  !    J(Np+1, l) = ... (T,p)
  !  enddo
  !  J(Np+1, Np+1) = ... (T,T)
  !---------------------------------------------------------------------------
  
  blkmat1 = 0.d0
  blkmat2 = 0.d0

  ! Accumulation terms
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)

    p1 = (ng-1)*option%nflowdof ! = 1 + (ng-1)*option%nflowdof-1
    t1 = p1 + 1           ! = 2 + (ng-1)*option%nflowdof-1
    c1 = t1 + 1           ! = 3 + (ng-1)*option%nflowdof-1

    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    !sum over fluid phases
    elempp = 0.d0
    elempt = 0.d0
    daccep = 0.d0
    daccet = 0.d0
    do j = 1, option%nphase
      jn = j + (n-1)*option%nphase
      jng = j + (ng-1)*option%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.

      !deriv wrt pressure: (p,p)
      elempp = elempp + pvoldt * d_p_loc_p(jng) 
       
      !deriv wrt temperature: (p,T)
      elempt = elempt + pvoldt * d_t_loc_p(jng) 
      
      !deriv of energy accumulation term wrt p: (T,p)
      daccep = daccep + d_p_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_p_loc_p(jng) - option%scale

      !deriv of energy accumulation term wrt T: (T,T)
      daccet = daccet + d_t_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_t_loc_p(jng)

!     print *,'TH: ',d_p_loc_p(jn),d_t_loc_p(jng),ddensity_loc_p(jng), &
!     hh_loc_p(jng),h_p_loc_p(jng),h_t_loc_p(jng),option%scale
    enddo
    
    !(p,p)
    if (option%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,p1,1,p1,elempp,ADD_VALUES,ierr)
    else
      blkmat1(1,1) = elempp
    endif

    !(p,T)
    if (option%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,p1,1,t1,elempt,ADD_VALUES,ierr)
    else
      blkmat1(1,2) = elempt
    endif

    !energy eqn. - cross term: (T,p)
    elem1 = pvoldt * daccep
    if (option%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(2,1) = elem1
    endif

    !energy eqn. - diagonal term: (T,T)
    i = ithrm_loc_p(ng)
    elem1 = pvoldt * daccet + (1.d0-porosity_loc_p(ng))*grid%dencpr(i)*voldt
    if (option%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(2,2) = elem1
    endif
    
    !tracer: (C,C)
    elem1 = pvoldt * grid%ret
    if (option%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(3,3) = elem1
    endif

    !kinetic rate term
    if (option%rk > 0.d0) then
      if (option%ityprxn == 1) then
        if (phis_p(n) > 0.d0 .or. CCONC_LOC(ng)/option%ceq > 1.d0+eps) then
          !(C,C)
          elem1 = option%rk * option%area_var(n) / option%ceq * volume_p(n)
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = blkmat1(3,3) + elem1
          endif

          !(T,C)
          elem1 = -option%rk * option%area_var(n)/option%ceq * option%delHs * &
          volume_p(n)
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,t1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(2,3) = elem1
          endif
        endif
      else if (option%ityprxn == 2) then
        if (phis_p(n) > 0.d0 .or. TTEMP_LOC(ng) < 0.d0+eps) then
          !(p,T)
          elem1 = option%rk * option%area_var(n) * volume_p(n)
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(1,2) = blkmat1(1,2) + elem1
          endif

          !(T,T)
          elem1 = -option%rk * option%area_var(n) * option%delHs * &
          volume_p(n)
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(2,2) = blkmat1(2,2) + elem1
          endif
        endif
      endif
    endif
    if (option%iblkfmt == 1) then
      call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat1, &
      ADD_VALUES,ierr)
    endif
  enddo

  !---------------------------------------------------------------------------
  ! Flux terms for interior nodes
  !---------------------------------------------------------------------------
  do iconn = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(iconn)
    m2 = grid%nd2(iconn)
    blkmat1 =0.D0
    blkmat2 =0.D0
    
    n1 = grid%nG2L(m1)
    n2 = grid%nG2L(m2)

    dd1 = grid%dist1(iconn)
    dd2 = grid%dist2(iconn)
    
    ip1 = option%iperm1(iconn)
    ip2 = option%iperm2(iconn)
    
!   perm1 = perm_loc_p(ip1+3*(m1-1))
!   perm2 = perm_loc_p(ip2+3*(m2-1))
    
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(m1)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(m1)
    else
      perm1 = perm_zz_loc_p(m1)
    endif
    
    if (ip2 == 1) then
      perm2 = perm_xx_loc_p(m2)
    else if (ip2 == 2) then
      perm2 = perm_yy_loc_p(m2)
    else
      perm2 = perm_zz_loc_p(m2)
    endif
    
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    por1 = porosity_loc_p(m1)
    por2 = porosity_loc_p(m2)
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * option%difaq
    dtrans = diff * grid%area(iconn)
      
    gravity = option%fmwh2o * option%gravity * grid%delz(iconn)
    
    dfluxp1 = 0.d0
    dfluxp2 = 0.d0
    dfluxt1 = 0.d0
    dfluxt2 = 0.d0
    
    p1 = (m1-1) * option%nflowdof
    t1 = p1 + 1
    c1 = t1 + 1
    p2 = (m2-1) * option%nflowdof
    t2 = p2 + 1
    c2 = t2 + 1

    do j = 1, option%nphase
      jm1 = j + (m1-1) * option%nphase
      jm2 = j + (m2-1) * option%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

!     D1 = perm1 / viscosity_loc_p(jm1)
!     D2 = perm2 / viscosity_loc_p(jm2)
!     Dq = (D1 * D2) / (dd2*D1 + dd1*D2)

      D1 = perm1 * viscosity_loc_p(jm2)
      D2 = perm2 * viscosity_loc_p(jm1)
      den = dd2*D1 + dd1*D2
      Dq = (perm1 * perm2) / den

      density_ave = f2 * ddensity_loc_p(jm1) + f1* ddensity_loc_p(jm2)

      v_darcy = -Dq * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) & 
                - gravity * density_ave)

      q = v_darcy * grid%area(iconn)
      
      qp1 =  Dq * grid%area(iconn) * (1.d0 + gravity * f2 * d_p_loc_p(jm1))
      qp2 = -Dq * grid%area(iconn) * (1.d0 - gravity * f1 * d_p_loc_p(jm2))

      qt1 = -q * dd1 * perm2 / den * v_t_loc_p(jm1) &
            + Dq * gravity * f2 * d_t_loc_p(jm1) * grid%area(iconn)
      qt2 = -q * dd2 * perm1 / den * v_t_loc_p(jm2) &
            + Dq * gravity * f1 * d_t_loc_p(jm2) * grid%area(iconn)

      trans1 = density_ave * qp1 + q * f2 * d_p_loc_p(jm1)
      trans2 = density_ave * qp2 + q * f1 * d_p_loc_p(jm2)

      qdt1 = density_ave * qt1 + q * f2 * d_t_loc_p(jm1)
      qdt2 = density_ave * qt2 + q * f1 * d_t_loc_p(jm2)

      !upstream weighting 
      if (q > 0.d0) then
        mu = m1
        u1 = q
        u2 = 0.d0
      else
        mu = m2
        u1 = 0.d0
        u2 = q
      endif
      
      ! note: qp1,2; qt1,2 are nonzero for q>0 and q<0
      if (m1 == mu) then
        hm1 = hh_loc_p(jm1)
        dfluxp1 = dfluxp1 + q * (density_ave * h_p_loc_p(jm1) &
                            + f2 * d_p_loc_p(jm1) * hm1) &
                            + qp1 * density_ave * hm1
        dfluxt1 = dfluxt1 + q * (density_ave * h_t_loc_p(jm1) &
                            + f2 * d_t_loc_p(jm1) * hm1) &
                            + qt1 * density_ave * hm1
        dfluxp2 = dfluxp2 + (q * f1 * d_p_loc_p(jm2) + qp2 * density_ave) * hm1
        dfluxt2 = dfluxt2 + (q * f1 * d_t_loc_p(jm2) + qt2 * density_ave) * hm1
      else
        hm2 = hh_loc_p(jm2)
        dfluxp1 = dfluxp1 + (q * f2 * d_p_loc_p(jm1) + qp1 * density_ave) * hm2
        dfluxt1 = dfluxt1 + (q * f2 * d_t_loc_p(jm1) + qt1 * density_ave) * hm2
        dfluxp2 = dfluxp2 + q * (density_ave * h_p_loc_p(jm2) &
                            + f1 * d_p_loc_p(jm2) * hm2) &
                            + qp2 * density_ave * hh_loc_p(jm2)
        dfluxt2 = dfluxt2 + q * (density_ave * h_t_loc_p(jm2) &
                            + f1 * d_t_loc_p(jm2) * hm2) &
                            + qt2 * density_ave * hm2
      endif
    enddo ! end loop phases (not really correct---but nphase=1 so no harm)
    
    !heat flux terms for thermal conduction
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    Dk1 = option%ckwet(i1)
    Dk2 = option%ckwet(i2)

    Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)

    cond = Dk * grid%area(iconn)
    
    if (q > 0.d0) then
      cupstrm = CCONC_LOC(m1)
    else
      cupstrm = CCONC_LOC(m2)
    endif

      ! Now add the flux contributions for this phase.
      ! Note that fluxes through a downstream face should be added to the
      ! residual component at the cell, while fluxes through an upstream  
      ! face should be subtracted.  (The divergence gives the net OUTFLOW 
      ! rate per unit volume.)  Thus, when working with pressure 
      ! differences,(ppressure(jm2) - ppressure(jm1)) should be 
      ! *subtracted* at the upstream node n1 because q = -D*div(P).
      
      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        ! liquid flux terms: (p,p)
        elem1 = trans1
        elem2 = trans2
!       print *,'THC-pp1',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p1,1,p2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
          blkmat2(1,1) = elem2
        endif

        ! liquid flux terms: (p,T)
        elem1 = qdt1
        elem2 = qdt2
!       print *,'THC-pT1',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p1,1,t2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,2) = elem1
          blkmat2(1,2) = elem2
        endif

        ! tracer flux terms: (C,C)
        elem1 = option%fc * u1 + dtrans
        elem2 = option%fc * u2 - dtrans
!       print *,'THC-CC1',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,c1,1,c2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(3,3) = elem1
          blkmat2(3,3) = elem2
        endif
    
      !(T,T)
      elem1 = dfluxt1 + cond
      elem2 = dfluxt2 - cond
!     print *,'THC-TT1',elem1,elem2
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t1,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
        blkmat2(2,2) = elem2
      endif

      ! heat flux terms: (T,p)
      elem1 = dfluxp1
      elem2 = dfluxp2
!     print *,'THC-Tp1',iconn,elem1,elem2,qp1,qp2,tupstrm
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t1,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
        blkmat2(2,1) = elem2
      endif

      ! tracer flux terms: (C,p)
      elem1 = option%fc * qp1 * cupstrm
      elem2 = option%fc * qp2 * cupstrm
!     print *,'THC-Cp1',elem1,elem2
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c1,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,1) = elem1
        blkmat2(3,1) = elem2
      endif

      ! tracer flux terms: (C,T)
      elem1 = option%fc * qt1 * cupstrm
      elem2 = option%fc * qt2 * cupstrm
!     print *,'THC-CT1',elem1,elem2,qt1,qt2,cupstrm
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c1,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c1,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,2) = elem1
        blkmat2(3,2) = elem2
      endif
      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,m1-1,1,m1-1,blkmat1, &
        ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,m1-1,1,m2-1,blkmat2, &
        ADD_VALUES,ierr)
      endif
    endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
      ! liquid flux terms: (p,p)
        elem1 = -trans1
        elem2 = -trans2
!       print *,'THC-pp2',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p2,1,p1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p2,1,p2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
          blkmat2(1,1) = elem2
        endif

        ! liquid flux terms: (p,T)
        elem1 = -qdt1
        elem2 = -qdt2
!       print *,'THC-pT2',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p2,1,t1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p2,1,t2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,2) = elem1
          blkmat2(1,2) = elem2
        endif

      ! tracer flux terms: (C,C)
        elem1 = -(option%fc * u1 + dtrans)
        elem2 = -(option%fc * u2 - dtrans)
!       print *,'THC-CC2',elem1,elem2
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c2,1,c1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,c2,1,c2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(3,3) = elem1
          blkmat2(3,3) = elem2
        endif

      !(T,T)
      elem1 = -(dfluxt1 + cond)
      elem2 = -(dfluxt2 - cond)
!     print *,'THC-TT2',elem1,elem2
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t2,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t2,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
        blkmat2(2,2) = elem2
      endif

      ! heat flux terms: (T,p)
      elem1 = -dfluxp1
      elem2 = -dfluxp2
!     print *,'TH-Tp2',iconn,elem1,elem2,qp1,qp2,tupstrm,dfluxp1,dfluxp2
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t2,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t2,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
        blkmat2(2,1) = elem2
      endif

      ! tracer flux terms: (C,p)
      elem1 = -option%fc * qp1 * cupstrm
      elem2 = -option%fc * qp2 * cupstrm
!     print *,'THC-Cp2',elem1,elem2
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c2,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c2,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,1) = elem1
        blkmat2(3,1) = elem2
      endif

      ! tracer flux terms: (C,T)
      elem1 = -option%fc * qt1 * cupstrm
      elem2 = -option%fc * qt2 * cupstrm
!     print *,'THC-CT2',elem1,elem2,qt1,qt2,cupstrm
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c2,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c2,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,2) = elem1
        blkmat2(3,2) = elem2
      endif
      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,m2-1,1,m2-1,blkmat2, &
        ADD_VALUES,ierr)
        call MatSetValuesBlockedLocal(A,1,m2-1,1,m1-1,blkmat1, &
        ADD_VALUES,ierr)
      endif
    endif

  enddo ! end loop connections

  !------------------------------------------------------------------------
! Flux terms for boundary nodes.
!------------------------------------------------------------------------

  do iconn = 1, grid%nconnbc

    m = grid%mblkbc(iconn)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)
    blkmat1 =0.D0
    p1 = (ng-1)*option%nflowdof
    t1 = p1 + 1
    c1 = t1 + 1

    ibc = grid%ibconn(iconn)
    ip1 = grid%ipermbc(iconn)

!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
        
    gravity = option%fmwh2o * option%gravity * grid%delzbc(iconn)
    
    if(option%ibndtyp(ibc) == 2) then
      ! solve for pb from Darcy's law given qb /= 0
      option%tempbc(iconn) = TTEMP_LOC(ng)
    else if(option%ibndtyp(ibc) == 3) then
      option%tempbc(iconn) = TTEMP_LOC(ng)
    endif

    do j = 1, option%nphase
      if (j == option%jh2o) then
        !pure water eos
        call wateos(option%tempbc(iconn),option%pressurebc(j,iconn),dw_kg, &
        dw_mol,field%d_p_bc(j),field%d_t_bc(j),field%hh_bc(j), &
        field%h_p_bc(j),field%h_t_bc(j),option%scale,ierr)
        field%density_bc(j) = dw_mol
        call PSAT(option%tempbc(iconn),sat_pressure,ierr)
        call VISW_noderiv(option%tempbc(iconn),option%pressurebc(j,iconn), &
        sat_pressure,field%viscosity_bc(j),ierr)

!     else if (j == grid%jco2) then

!       add additional fluid phases here ...

      else ! for testing purposes only
        jng = j + (ng-1) * option%nphase
        field%viscosity_bc(j) = viscosity_loc_p(jng)
        field%hh_bc(j) = hh_loc_p(jng)
        field%density_bc(j) = ddensity_loc_p(jng)
      endif
    enddo

    if(option%ibndtyp(ibc) == 1) then 
    
    ! Dirichlet BC: fixed p, T, C
    
      dfluxp = 0.d0
      dfluxt = 0.d0
      do j=1, option%nphase
        jm = j + (m-1) * option%nphase
        jng = j + (ng-1) * option%nphase
        
        ! For the value of the "diffusion" constant at the interface, I 
        ! just use its value as defined in the boundary cell.  I don't 
        ! know if this is the best thing to do, but it will work for now. 
        ! (E.g. temperature could be different at the boundary compared  
        ! to node center resulting in different density, viscosity etc. 
        ! - pcl)
        Dq = perm1 / field%viscosity_bc(j) / grid%distbc(iconn) * &
        grid%areabc(iconn)

        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - option%pressurebc(j,iconn) &
                  - gravity * field%density_bc(j))
        q = v_darcy
        
        qp1 = -Dq
!       qt1 = -q * v_t_loc_p(jng) / viscosity_loc_p(jng)
        qt1 = 0.d0

        trans = qp1 * field%density_bc(j) ! dR/dp = dq/dp rho

        ! (p,p)
        elem1 = -trans
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
        endif

        ! (p,T) = 0
!       elem1 = -qt1 * field%density_bc(j)
!       if (option%iblkfmt == 0) then
!         call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
!       else
!         blkmat1(1,2) = elem1
!       endif

        if (q > 0.d0) then
          dfluxp = dfluxp + qp1 * field%density_bc(j) * field%hh_bc(j)
          dfluxt = dfluxt + qt1 * field%density_bc(j) * field%hh_bc(j)
        else
          dfluxp = dfluxp + q * (ddensity_loc_p(jng) * h_p_loc_p(jng) &
          + d_p_loc_p(jng) * hh_loc_p(jng)) &
          + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)
          dfluxt = dfluxt + q * (ddensity_loc_p(jng) * h_t_loc_p(jng) &
          + d_t_loc_p(jng) * hh_loc_p(jng)) &
          + qt1 * ddensity_loc_p(jng) * hh_loc_p(jng)
        endif

        ! (C,p)
        if (q > 0.d0) then
          elem1 = -option%fc * qp1 * option%concbc(iconn)
        else
          elem1 = -option%fc * qp1 * CCONC_LOC(ng)
        endif
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif

        ! (C,T)
        if (q > 0.d0) then
          elem1 = -option%fc * qt1 * option%concbc(iconn)
        else
          elem1 = -option%fc * qt1 * CCONC_LOC(ng)
        endif
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif
      
        ! (C,C)
        trans = porosity_loc_p(ng) * option%difaq * grid%areabc(iconn) / &
                grid%distbc(iconn)
        if (q < 0.d0) then
          elem1 = -option%fc * q + trans ! check sign of q
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = elem1
          endif
        else
          elem1 = trans
          if (option%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = elem1
          endif
        endif
      enddo

      ! (T,p)
      elem1 = -dfluxp
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
      endif
      
      ! (T,T)
      i1 = ithrm_loc_p(ng)
      cond = option%ckwet(i1) * grid%areabc(iconn) / grid%distbc(iconn)
      elem1 = -dfluxt + cond
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif
  
   else if(option%ibndtyp(ibc) == 4) then 

 !#if 0
      ! (T,T)
         i1 = ithrm_loc_p(ng)
      cond = option%ckwet(i1) * grid%areabc(iconn) / grid%distbc(iconn)
      elem1 =  cond
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif
!#endif
       
    else if(option%ibndtyp(ibc) == 2) then 

    ! constant velocity q, grad T, C = 0
    
    else if(option%ibndtyp(ibc) == 4) then
        ! (T,T)
         i1 = ithrm_loc_p(ng)
      cond = option%ckwet(i1) * grid%areabc(iconn) / grid%distbc(iconn)
      elem1 =  cond
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif

    else if(option%ibndtyp(ibc) == 3) then 
    
    ! Dirichlet BC: fixed p, grad T, C = 0
    
      dfluxt = 0.d0
      dfluxp = 0.d0
      do j=1, option%nphase
        jm = j + (m-1) * option%nphase
        jng = j + (ng-1) * option%nphase
        
        Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(iconn) * &
        grid%areabc(iconn)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - option%pressurebc(j,iconn) &
                  - gravity * ddensity_loc_p(jng))
        q = v_darcy

!       Note: p = pb at boundary
!       qp1 = -Dq * (1.d0 - gravity * d_p_loc_p(jng))
        qp1 = -Dq
        
        qt1 = -q * v_t_loc_p(jng) / viscosity_loc_p(jng) &
              + Dq * gravity * d_t_loc_p(jng)

!       trans1 = -(qp1 * ddensity_loc_p(jng) + q * d_p_loc_p(jng))
        trans1 = -qp1 * ddensity_loc_p(jng)

        qdt1 = -(qt1 * ddensity_loc_p(jng) + q * d_t_loc_p(jng))

        ! (p,p)
        elem1 = trans1
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
        endif
      
        ! (p,T)
        elem1 = qdt1
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(1,2) = elem1
        endif

!       dfluxp = dfluxp + q * (ddensity_loc_p(jng) * h_p_loc_p(jng) &
!         + d_p_loc_p(jng) * hh_loc_p(jng)) &
!         + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

!       dfluxp = dfluxp + q * ddensity_loc_p(jng) * h_p_loc_p(jng) &
!         + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

        dfluxp = dfluxp + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

        dfluxt = dfluxt + q * (ddensity_loc_p(jng) * h_t_loc_p(jng) &
          + d_t_loc_p(jng) * hh_loc_p(jng)) &
          + qt1 * ddensity_loc_p(jng) * hh_loc_p(jng)

        ! (C,p)
        elem1 = -option%fc * qp1 * CCONC_LOC(ng)
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif

        ! (C,T)
        elem1 = -option%fc * qt1 * CCONC_LOC(ng)
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,t1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,2) = elem1
        endif
      
        ! (C,C), grad C = 0
        elem1 = -option%fc * q
        if (option%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,3) = elem1
        endif
      enddo

      ! (T,p)
      elem1 = -dfluxp
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
      endif
      
      ! (T,T), grad T = 0
      elem1 = -dfluxt
      if (option%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (option%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif
    endif
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%ddensity_loc, ddensity_loc_p, ierr)
! call VecRestoreArrayF90(field%density, density_p, ierr)
  call VecRestoreArrayF90(thc_field%viscosity_loc, viscosity_loc_p, ierr)
  
  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(thc_field%d_p_loc, d_p_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%d_t_loc, d_t_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%h_p_loc, h_p_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%h_t_loc, h_t_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%v_p_loc, v_p_loc_p, ierr)
  call VecRestoreArrayF90(thc_field%v_t_loc, v_t_loc_p, ierr)

  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr)
  
  if (option%rk > 0.d0) then
    call VecRestoreArrayF90(field%phis,phis_p,ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  B = A

! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
#endif  ! End of broken code. 
  end subroutine THCJacobian

! ************************************************************************** !
!
! THCInitializeSolidReaction: Allocates and initializes arrays associated with
!                          mineral reactions
! author: Glenn Hammond
! date: 11/15/07
!
! ************************************************************************** !
subroutine THCInitializeSolidReaction(realization)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  type(realization_type) :: realization
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscInt :: icell
  PetscReal, pointer :: phis_p(:)
  PetscErrorCode :: ierr
  
  grid => realization%grid
  option => realization%option
  field => realization%field
  
  if (thc_option%rk > 0.d0) then
    allocate(thc_option%area_var(grid%nlmax))
    allocate(thc_option%rate(grid%nlmax))
    call VecGetArrayF90(thc_field%phis,phis_p,ierr)
    do icell = 1, grid%nlmax
      phis_p(icell) = thc_option%phis0
      thc_option%area_var(icell) = 1.d0
    enddo
    call VecRestoreArrayF90(thc_field%phis,phis_p,ierr)
  endif
  
end subroutine THCInitializeSolidReaction


! ************************************************************************** !
!
! THCGetTecplotHeader: Returns a Tecplot file header
! author: 
! date: 
!
! ************************************************************************** !
function THCGetTecplotHeader(realization)

  use Realization_module
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THCGetTecplotHeader
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = 'VARIABLES=' // &
           '"X [m]",' // &
           '"Y [m]",' // &
           '"Z [m]",' // &
           '"T [C]",' // &
           '"P [Pa]",' // &
           '"sl",' // &
           '"sg",' // &
           '"Ul",' // &
           '"Ug",'
  do i=1,option%nspec
    write(string2,'(''"Xl('',i2,'')",'')') i
    string = trim(string) // trim(string2)
  enddo
  do i=1,option%nspec
    write(string2,'(''"Xg('',i2,'')",'')') i
    string = trim(string) // trim(string2)
  enddo
  if (thc_option%rk > 0.d0) then
    string = trim(string) // '"Volume Fraction"'
  endif
  string = trim(string) // ',"Phase"'
  if (associated(field%imat)) then
    string = trim(string) // ',"Material_ID"'
  endif
  
  THCGetTecplotHeader = string

end function THCGetTecplotHeader

! ************************************************************************** !
!
! THCGetVarFromArray: Extracts variables indexed by ivar and isubvar
! author: 
! date: 
!
! ************************************************************************** !
subroutine THCGetVarFromArray(realization,vec,ivar,isubvar)

  use Realization_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none
  
  PetscInt, parameter :: TEMPERATURE = 4
  PetscInt, parameter :: PRESSURE = 5
  PetscInt, parameter :: LIQUID_SATURATION = 6
  PetscInt, parameter :: GAS_SATURATION = 7
  PetscInt, parameter :: LIQUID_ENERGY = 8
  PetscInt, parameter :: GAS_ENERGY = 9
  PetscInt, parameter :: LIQUID_MOLE_FRACTION = 10
  PetscInt, parameter :: GAS_MOLE_FRACTION = 11
  PetscInt, parameter :: VOLUME_FRACTION = 12
  PetscInt, parameter :: PHASE = 13
  PetscInt, parameter :: MATERIAL_ID = 14

  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal, pointer :: vec_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%grid
  field => realization%field

  call VecGetArrayF90(vec,vec_ptr,ierr)

  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY)
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)    
        select case(ivar)
          case(TEMPERATURE)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(PRESSURE)
            vec_ptr(local_id) = 0.d0 ! to be provided
          case(LIQUID_ENERGY)
            vec_ptr(local_id) = 0.d0 ! to be provided
        end select
      enddo
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = field%imat(grid%nL2G(local_id))
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr)

end subroutine THCGetVarFromArray

end module THC_module

#undef PPRESSURE_LOC
#undef PPRESSURE
#undef PRESSURE
#undef TTEMP_LOC
#undef TTEMP
#undef TEMP
#undef CCONC_LOC
#undef CCONC
#undef CONC

