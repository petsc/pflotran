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
    

    real*8, pointer :: pressurebc(:,:)
      ! For a Dirichlet BC, pressurebc(j,ibc) gives the partial pressure 
      ! for phase j along the BC block ibc.
    real*8, pointer :: velocitybc(:,:)
      ! For a Neumann BC, velocitybc(j,ibc) gives the velocity q for phase
      ! j along BC block ibc.
    real*8, pointer :: tempbc(:),concbc(:),sgbc(:),xphi_co2_bc(:),xxphi_co2_bc(:)
    real*8, pointer :: xxbc(:,:), varbc(:)
    integer, pointer:: iphasebc(:)

    !block BC values read from input
    real*8, pointer :: pressurebc0(:,:)
    real*8, pointer :: velocitybc0(:,:)
#if 0
! don't believe that these are need anymore
    real*8, pointer :: tempbc0(:),concbc0(:),sgbc0(:)
    real*8, pointer :: xxbc0(:,:)
    integer, pointer:: iphasebc0(:)  
#endif

!   solid reaction rate
    integer*4 :: ityprxn
    real*8 :: rk=0.d0, phis0, areas0, pwrsrf, vbars, ceq, delHs, delEs, wfmts
    real*8 ::qu_kin, yh2o_in_co2=0.D0
    
!   breakthrough curves
    integer :: ibrkcrv = 0
    integer*4, pointer :: ibrktyp(:),ibrkface(:)
    
!   dual continuum
    integer :: idcdm = 0, idcmblk = 0
    real*8, pointer :: fracture_aperture(:), matrix_block(:)
    
    integer*4, pointer :: icap_reg(:),ithrm_reg(:)
    real*8 :: scale
    real*8, pointer :: rock_density(:),cpr(:),dencpr(:),ckdry(:),ckwet(:), &
                       tau(:),cdiff(:),cexp(:)
    real*8, pointer :: swir(:),lambda(:),alpha(:),pckrm(:),pcwmax(:),pcbetac(:), &
                       pwrprm(:),sir(:,:)
    integer, pointer:: icaptype(:)
!geh material id
    integer, pointer :: imat(:)
    real*8 :: m_nacl
    real*8 :: difaq, delhaq, gravity, fmwh2o= 18.0153D0, fmwa=28.96D0, &
              fmwco2=44.0098D0, eqkair, ret=1.d0, fc=1.d0
    
    integer :: ihydrostatic = 0,ideriv = 1
    real*8 :: dTdz,beta,tref,pref,conc0
    real*8 :: hydro_ref_xyz(3)
    
!   table lookup
    integer :: itable=0

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

  !physical constants and defult variables
  field%difaq = 1.d-9 ! m^2/s read from input file
  field%delhaq = 12.6d0 ! kJ/mol read from input file
  field%gravity = 9.8068d0    ! m/s^2
  field%tref   = 50.D0
  field%fmwh2o = 18.01534d0 ! kg H2O/mol H2O
  field%fmwco2 = 44.0098d0
  field%eqkair = 1.d10 ! Henry's constant for air: Xl = eqkair * pa

  ! default brine concentrations
  field%m_nacl = 0.d0
  
  ! nullify PetscVecs
  field%conc = 0
  field%xmol = 0
  
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
  
  ! all kinds of stuff needs to be added here.
  
  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
