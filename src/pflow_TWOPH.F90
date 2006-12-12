!======================================================================


!#include "pflow2PH.F90"

#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*grid%ndof)
#define PPRESSURE(j,n) xx_p(j+(n-1)*grid%ndof)
#define PRESSURE(j,n) yy_p(j+(n-1)*grid%ndof)
#define TTEMP_LOC(n) xx_loc_p(grid%nphase+1+(n-1)*grid%ndof)
#define TTEMP(n) xx_p(grid%nphase+1+(n-1)*grid%ndof)
#define TEMP(n) yy_p(grid%nphase+1+(n-1)*grid%ndof)
#define CCONC_LOC(n) xx_loc_p(grid%nphase+2+(n-1)*grid%ndof)
#define CCONC(n) xx_p(grid%nphase+2+(n-1)*grid%ndof)
#define CONC(n) yy_p(grid%nphase+2+(n-1)*grid%ndof)

module TWOPH_module

  use pflow_gridtype_module

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
!#ifdef USE_PETSC216
!#include "include/finclude/petscsles.h"
!#endif
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petsclog.h"

!#include "pflow_gridtype.h"

public TWOPHResidual, TWOPHJacobian

contains

  subroutine TWOPHResidual(snes, xx, r, grid)

  use water_eos_module
  
  implicit none
  
  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Vec, intent(out) :: r
  type(pflowGrid), intent(inout) :: grid

! external VISW, PSAT, WATEOS
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  integer :: ierr
  integer :: n, ng, nc
  integer :: i, i1, i2, j, j1, j2, j3, jn1, jn2, jng, jm, jm1, jm2
  integer :: m, m1, m2, n1, n2, ip1, ip2
  integer :: ibc  ! Index that specifies a boundary condition block.
  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:), &
               ppressure_p(:), ttemp_p(:), &
               density_p(:), ddensity_p(:), ddensity_loc_p(:), &
               sat_p(:), ssat_loc_p(:), &
               xmol_p(:), xxmol_loc_p(:), &
               viscosity_p(:), viscosity_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               accum_p(:), &
               vl_p(:), iphas_p(:), &
               d_p_p(:), d_p_loc_p(:), &
               d_t_p(:), d_t_loc_p(:), &
               h_p(:), hh_p(:), hh_loc_p(:), &
               h_p_p(:), h_p_loc_p(:), &
               h_t_p(:), h_t_loc_p(:), &
               v_p_p(:), v_p_loc_p(:), &
               v_t_p(:), v_t_loc_p(:), &
               ithrm_p(:), ithrm_loc_p(:), icap_loc_p(:), iphas_loc_p(:)
  real*8 :: dd1, dd2, diflux, diff, eeng, cond, trans, density_ave, &
            hflx, fluxc, fluxe, fluxh, flux, fluxv, fluxbc, &
            fluxp, fluxp1, fluxp2, q, u1, u2, v_darcy, pvoldt, voldt, por, &
            accum1, accum2
  real*8 :: dd, f1, f2, por1, por2, perm1, perm2
  real*8 :: D  ! "Diffusion" constant for a phase
  real*8 :: D1, D2 ! "Diffusion" constants upstream and downstream of a face
  real*8 :: ps,pv,pl,pg,pa,sl,sg,tc,xl1,xl2,xg1,xg2 !,pc
  real*8 :: dw_kg, dw_mol, dg_kg, dg_mol, dencp
  real*8 :: sat_pressure

!----------------------
! degrees of freedom:
!----------------------
! xx   liq   gas   2ph
!----------------------
!  1    p     p     p
!  2    Xl    Xg    pa
!  3    T     T     sg
!----------------------  !---------------------------------------------------------------------------

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%ppressure, ppressure_p, ierr) ; CHKERRQ(ierr)
  call VecGetArrayF90(grid%ttemp, ttemp_p, ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%hh, hh_p, ierr)
  call VecGetArrayF90(grid%viscosity, viscosity_p, ierr)
  call VecGetArrayF90(grid%ithrm, ithrm_p, ierr)
  call VecGetArrayF90(grid%xxmol_loc, xxmol_loc_p, ierr)
! call VecGetArrayF90(grid%xmol, xmol_p, ierr)
  call VecGetArrayF90(grid%ssat_loc, ssat_loc_p, ierr)
! call VecGetArrayF90(grid%sat, sat_p, ierr)
  
  if (grid%ideriv == 1) then
    call VecGetArrayF90(grid%d_p, d_p_p, ierr)
    call VecGetArrayF90(grid%d_t, d_t_p, ierr)
    call VecGetArrayF90(grid%h_p, h_p_p, ierr)
    call VecGetArrayF90(grid%h_t, h_t_p, ierr)
    call VecGetArrayF90(grid%v_p, v_p_p, ierr)
    call VecGetArrayF90(grid%v_t, v_t_p, ierr)
  endif
  
  ! Compute accumulation terms for interior nodes.
  ! Calculate the density, viscosity, and enthalpy of water at step k+1 at 
  ! each local node.
  !---------------------------------------------------------------------------
  
  do n = 1, grid%nlmax
    ng = grid%nL2G(n)
    voldt = volume_p(n) / grid%dt
    por = porosity_loc_p(ng)
    pvoldt = por * voldt
    j1 = 1+(n-1)*grid%ndof
    j2 = 2+(n-1)*grid%ndof
    j3 = 3+(n-1)*grid%ndof
    i = ithrm_p(n)
    dencp = grid%dencpr(i)
    if (iphas_p(n) == 0) then ! gas
      pg = xx_p(j1)
      tc = xx_p(j2)
      xg2 = xx_p(j3)
      xg1 = 1.d0 - xg2
      
      !pure gas steam eos
      jn2 = grid%jgas + (n-1)*grid%nphase
      call steameos(tc,pg,pa,dg_kg,dg_mol,d_p_p(jn2), &
      d_t_p(jn2),hh_p(jn2),h_p_p(jn2),h_t_p(jn2),grid%scale,ierr)
      !call visg(...)
        
!     print *,j,n,grid%nphase,tc,pg,dw_mol,d_p_p(jn2), &
!     d_t_p(jn2),h_p(jn2),h_p_p(jn2),h_t_p(jn2),viscosity_p(jn2)

      ppressure_p(jn2) = pg
      ttemp_p(n) = tc
      xxmol_loc_p(jn2) = xg2
      ddensity_p(jn2) = dg_mol

      accum1 = voldt * (por * dg_mol * xg1 - accum_p(j1))
      accum2 = voldt * (por * dg_mol * xg2 - accum_p(j3))
      
      !heat residual
      eeng = ddensity_p(jn2)*hh_p(jn2)-grid%scale*pg
      ! set residual
      r_p(j1) = accum1
      r_p(j2) = voldt * (por * eeng + (1.d0-por) * dencp * tc - accum_p(j2))
      r_p(j3) = accum2

    else if (iphas_p(n) == 1) then ! liquid
    
      pl = xx_p(j1)
      tc = xx_p(j2)
      xl2 = xx_p(j3)
      xl1 = 1.d0 - xl2
    
      !pure gas water + air eos
      jn1 = grid%jh2o + (n-1)*grid%nphase
      call PSAT(tc, ps, ierr)
      if (grid%ideriv == 1) then
!       first arg temp, 2nd pressure, use molar units for density 
!       and other properties
        call wateos(tc,pl,dw_kg,dw_mol,d_p_p(jn1),d_t_p(jn1), &
        hh_p(jn1),h_p_p(jn1),h_t_p(jn1),grid%scale,ierr)
        call VISW(tc,pl,ps,viscosity_p(jn1),v_t_p(jn1),v_p_p(jn1),ierr)
      else
        call wateos_noderiv(tc,pl,dw_kg,dw_mol,hh_p(jn1),grid%scale,ierr)
        call VISW_noderiv(tc,pl,ps,viscosity_p(jn1),ierr)
      endif
      
      ppressure_p(jn1) = pl
      ttemp_p(n) = tc
      xxmol_loc_p(jn1) = xl2
      ddensity_p(jn1) = dw_mol

      accum1 = voldt * (por * dw_mol * xl1 - accum_p(j1))
      accum2 = voldt * (por * dw_mol * xl2 - accum_p(j3))
    
      !heat residual
      i = ithrm_p(n)
      eeng = dw_mol*hh_p(jn1)-grid%scale*pl
      
      r_p(j1) = accum1
      r_p(j2) = voldt * (por * eeng + (1.d0-por) * dencp * tc - accum_p(j2))
      r_p(j3) = accum2

    else if (iphas_p(n) == 2) then ! 2 phase
    
      pg = xx_p(j1)
      sl = xx_p(j2)
      xg2 = xx_p(j3)
      sg = 1.d0 - sl
      pa = xg2 * pg
      pv = pg - pv
!     tc = tsat(pv)
!     pc = pcap(sl)
!     pl = pg - pc
      xl2 = grid%eqkair * pa
      xl1 = 1.d0 - xl2
      
      jn1 = grid%jh2o + (n-1)*grid%nphase
      jn2 = grid%jgas + (n-1)*grid%nphase
    
      call PSAT(tc, ps, ierr)
      if (grid%ideriv == 1) then
!       first arg temp, 2nd pressure, use molar units for density 
!       and other properties
        call wateos(tc,pl,dw_kg,dw_mol,d_p_p(jn1),d_t_p(jn1), &
        hh_p(jn1),h_p_p(jn1),h_t_p(jn1),grid%scale,ierr)
        call VISW(tc,pl,ps,viscosity_p(jn1),v_t_p(jn1),v_p_p(jn1),ierr)
      else
        call wateos_noderiv(tc,pl,dw_kg,dw_mol,hh_p(jn1),grid%scale,ierr)
        call VISW_noderiv(tc,pl,ps,viscosity_p(jn1),ierr)
      endif
      call steameos(tc,pg,pa,dg_kg,dg_mol,d_p_p(jn2), &
      d_t_p(jn2),hh_p(jn2),h_p_p(jn2),h_t_p(jn2),grid%scale,ierr)
      !call visg(...)

      ppressure_p(jn2) = pg
      ttemp_p(n) = tc
      xxmol_loc_p(jn1) = xl2
      xxmol_loc_p(jn2) = xg2
      ddensity_p(jn1) = dw_mol
      ddensity_p(jn2) = dg_mol
      
      accum1 = voldt * (por * (sl * dw_mol * xl1 + sg * dg_mol * xg1) &
               - accum_p(j1))
      accum2 = voldt * (por * (sl * dw_mol * xl2 + sg * dg_mol * xg2) &
               - accum_p(j3))
    
      !heat residual
      eeng = dw_mol*hh_p(jn1)-grid%scale*pl + dg_mol*hh_p(jn2)-grid%scale*pg

      r_p(j1) = accum1
      r_p(j2) = voldt * (por * eeng + (1.d0-por) * dencp * tc - accum_p(j2))
      r_p(j3) = accum2
      
    endif
  enddo
  
  call VecRestoreArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecRestoreArrayF90(grid%hh, hh_p, ierr)
  call VecRestoreArrayF90(grid%viscosity, viscosity_p, ierr)

  if (grid%ideriv == 1) then  
    call VecRestoreArrayF90(grid%d_p, d_p_p, ierr)
    call VecRestoreArrayF90(grid%d_t, d_t_p, ierr)
    call VecRestoreArrayF90(grid%h_p, h_p_p, ierr)
    call VecRestoreArrayF90(grid%h_t, h_t_p, ierr)
    call VecRestoreArrayF90(grid%v_p, v_p_p, ierr)
    call VecRestoreArrayF90(grid%v_t, v_t_p, ierr)
  endif

  if (grid%flowsteps == 1) then
    call VecCopy(grid%ddensity, grid%density, ierr)
    call VecCopy(grid%hh, grid%h, ierr)
  endif
  !---------------------------------------------------------------------------
  ! Now that we have calculated the density and viscosity for all local 
  ! nodes, we can perform the global-to-local scatters.
  !---------------------------------------------------------------------------
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%ddensity,INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ssat, INSERT_VALUES, &
                            grid%ssat_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ssat, INSERT_VALUES, &
                          grid%ssat_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%xxmol, INSERT_VALUES, &
                            grid%xxmol_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%xxmol, INSERT_VALUES, &
                          grid%xxmol_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                            grid%hh_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                          grid%hh_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof,grid%viscosity,INSERT_VALUES, &
                            grid%viscosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof,grid%viscosity,INSERT_VALUES, &
                          grid%viscosity_loc, ierr)

! call DAGlobalToLocalBegin(grid%da_3_dof, grid%perm, &
!                           INSERT_VALUES, grid%perm_loc, ierr)
! call DAGlobalToLocalEnd(grid%da_3_dof, grid%perm, &
!                         INSERT_VALUES, grid%perm_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_xx, &
                            INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_xx, &
                          INSERT_VALUES, grid%perm_xx_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_yy, &
                            INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_yy, &
                          INSERT_VALUES, grid%perm_yy_loc, ierr)
  call DAGlobalToLocalBegin(grid%da_1_dof, grid%perm_zz, &
                            INSERT_VALUES, grid%perm_zz_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, grid%perm_zz, &
                          INSERT_VALUES, grid%perm_zz_loc, ierr)

  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%yy, yy_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(grid%h, h_p, ierr)
  call VecGetArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  call VecGetArrayF90(grid%accum, accum_p, ierr)
  
  call VecGetArrayF90(grid%icap_loc, icap_loc_p, ierr)
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  call VecGetArrayF90(grid%iphas_loc, iphas_loc_p, ierr)

  call VecGetArrayF90(grid%vl, vl_p, ierr)

  if (grid%ideriv == 1) then
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                            grid%d_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                          grid%d_p_loc, ierr)
                          
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_t, INSERT_VALUES, &
                            grid%d_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_t, INSERT_VALUES, &
                          grid%d_t_loc, ierr)
                          
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%h_p, INSERT_VALUES, &
                            grid%h_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_p, INSERT_VALUES, &
                          grid%h_p_loc, ierr)
                          
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%h_t, INSERT_VALUES, &
                            grid%h_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_t, INSERT_VALUES, &
                          grid%h_t_loc, ierr)
                          
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%v_p, INSERT_VALUES, &
                            grid%v_p_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%v_p, INSERT_VALUES, &
                          grid%v_p_loc, ierr)
                          
    call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%v_t, INSERT_VALUES, &
                            grid%v_t_loc, ierr)
    call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%v_t, INSERT_VALUES, &
                          grid%v_t_loc, ierr)

    call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecGetArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecGetArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecGetArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecGetArrayF90(grid%v_p_loc, v_p_loc_p, ierr)
    call VecGetArrayF90(grid%v_t_loc, v_t_loc_p, ierr)
  endif

!---------------------------------------------------------------------------
! Flux terms for interior nodes
!---------------------------------------------------------------------------
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1)
    n2 = grid%nG2L(m2)
    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)
    ip2 = grid%iperm2(nc)
    
!   perm1 = perm_loc_p(ip1+3*(m1-1))
!   perm2 = perm_loc_p(ip2+3*(m2-1))
    
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(m1)
    else if (ip1==2) then
      perm1 = perm_yy_loc_p(m1)
    else
      perm1 = perm_zz_loc_p(m1)
    endif
    
    if (ip2 == 1) then
      perm2 = perm_xx_loc_p(m2)
    else if (ip2==2) then
      perm2 = perm_yy_loc_p(m2)
    else
      perm2 = perm_zz_loc_p(m2)
    endif
    
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    fluxp1 = 0.d0
    fluxp2 = 0.d0
    fluxh = 0.d0
    fluxv = 0.d0
    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

      D1 = perm1 / viscosity_loc_p(jm1)
      D2 = perm2 / viscosity_loc_p(jm2)

      D = (D1 * D2) / (dd2*D1 + dd1*D2)

      density_ave = f2 * ddensity_loc_p(jm1) + f1 * ddensity_loc_p(jm2)

      v_darcy = -D * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) & 
                - grid%fmwh2o * density_ave * grid%gravity * grid%delz(nc))

      !store velocities defined at interfaces in PETSc Vec vl at upstream node
      if (n1 > 0) vl_p(ip1+3*(n1-1)) = v_darcy
      
      q = v_darcy * grid%area(nc)

      flux = density_ave * q

      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        fluxp1 = fluxp1 + flux ! interface on rhs
      endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
        fluxp2 = fluxp2 + flux ! interface on lhs
      endif
      
      !upstream weighting 
      if (q > 0.d0) then
        u1 = q
        u2 = 0.d0
      else
        u1 = 0.d0
        u2 = q
      endif
      
      fluxh = fluxh + u2 * ddensity_loc_p(jm2) * hh_loc_p(jm2) &
                    + u1 * ddensity_loc_p(jm1) * hh_loc_p(jm1)
                    
      fluxv = fluxv + u2 * CCONC_LOC(m2) + u1 * CCONC_LOC(m1)
    enddo

    !heat & tracer residual
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    D = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = D * grid%area(nc)
    hflx = cond*(TTEMP_LOC(m2) - TTEMP_LOC(m1))
    fluxe = fluxh - hflx

    por1 = porosity_loc_p(m1)
    por2 = porosity_loc_p(m2)
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
    diflux = diff * grid%area(nc) * (CCONC_LOC(m2) - CCONC_LOC(m1))
    fluxc = fluxv - diflux

    if (n1 > 0) then  ! If upstream node is not a ghost node...
      r_p(1+(n1-1)*grid%ndof) = r_p(1+(n1-1)*grid%ndof) + fluxp1
      r_p(2+(n1-1)*grid%ndof) = r_p(2+(n1-1)*grid%ndof) + fluxe   ! heat
      r_p(3+(n1-1)*grid%ndof) = r_p(3+(n1-1)*grid%ndof) + fluxc   ! tracer
    endif

    if (n2 > 0) then ! If downstream node is not a ghost node...
      r_p(1+(n2-1)*grid%ndof) = r_p(1+(n2-1)*grid%ndof) - fluxp2
      r_p(2+(n2-1)*grid%ndof) = r_p(2+(n2-1)*grid%ndof) - fluxe
      r_p(3+(n2-1)*grid%ndof) = r_p(3+(n2-1)*grid%ndof) - fluxc
    endif
  enddo

!---------------------------------------------------------------------------
! Flux terms for boundary nodes.
!---------------------------------------------------------------------------
  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)

!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1==2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
    
    if(grid%ibndtyp(ibc) == 2) then
      ! solve for pb from Darcy's law given qb /= 0
      grid%tempbc(ibc) = TTEMP_LOC(ng)
    else if(grid%ibndtyp(ibc) == 3) then
      grid%tempbc(ibc) = TTEMP_LOC(ng)
    endif

    do j = 1, grid%nphase
      if (j == grid%jh2o) then
        !pure water eos
        call PSAT(grid%tempbc(ibc), sat_pressure, ierr)  
        if (grid%ideriv == 1) then
          call wateos(grid%tempbc(ibc),grid%pressurebc(j,ibc),dw_kg,dw_mol, &
          grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j), &
          grid%h_t_bc(j),grid%scale,ierr)
          call VISW(grid%tempbc(ibc),grid%pressurebc(j,ibc),sat_pressure, &
          grid%viscosity_bc(j),grid%v_t_bc(j),grid%v_p_bc(j),ierr)
        
!         print *,'2ph-bc: ',nc,j,m,ng,ibc,grid%nphase,grid%tempbc(ibc), &
!         grid%pressurebc(j,ibc),dw_mol,dw_kg,grid%d_p_bc(j), &
!         grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j),grid%h_t_bc(j), &
!         sat_pressure,grid%viscosity_bc(j)

        else
          call wateos_noderiv(grid%tempbc(ibc),grid%pressurebc(j,ibc), &
          dw_kg,dw_mol,grid%hh_bc(j),grid%scale,ierr)
          call VISW_noderiv(grid%tempbc(ibc),grid%pressurebc(j,ibc), &
          sat_pressure,grid%viscosity_bc(j),ierr)
        endif
        grid%density_bc(j) = dw_mol
      
!     else if (j == grid%jco2) then

!       add additional fluid phases here ...

!       do n = 1, grid%nphase
!         ...
!       enddo
      else ! for testing purposes only
        jng = j + (ng-1) * grid%nphase
        grid%viscosity_bc(j) = viscosity_loc_p(jng)
        grid%hh_bc(j) = hh_loc_p(jng)
        grid%density_bc(j) = ddensity_loc_p(jng)
      endif
    enddo

    if(grid%ibndtyp(ibc) == 1) then  ! Dirichlet BC for p, T, C ...
    
      fluxp = 0.d0
      fluxh = 0.d0
      fluxc = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase

        D = perm1 / grid%viscosity_bc(j)
        
        v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        - grid%fmwh2o * ddensity_loc_p(jng) * grid%gravity * grid%delzbc(nc) &
                 ) / grid%distbc(nc)
        q = v_darcy * grid%areabc(nc)

        fluxbc = q * grid%density_bc(j)
        fluxp = fluxp - fluxbc

        !upstream weighting
        if (q > 0.d0) then
          fluxh = fluxh + q * grid%density_bc(j) * grid%hh_bc(j)
          fluxc = fluxc + q * grid%concbc(ibc)       ! note: need to change 
        else                          ! definition of C for multiple phases
          fluxh = fluxh + q * ddensity_loc_p(jng) * hh_loc_p(jng)
          fluxc = fluxc + q * CCONC_LOC(ng) 
        endif
        
!       print *,'pflowTH: ',nc,D,fluxh,q,hh_loc_p(jng),r_p(jmp1),fluxbc, &
!       PPRESSURE_LOC(j,ng),grid%pressurebc(j,ibc), &
!       perm_loc_p(ip1+3(ng-1)),viscosity_loc_p(jng)
      enddo
        
      r_p(1+(m-1)*grid%ndof) = r_p(1+(m-1)*grid%ndof) + fluxp

      !heat residual
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      r_p(2+(m-1)*grid%ndof) = r_p(2+(m-1)*grid%ndof) &
                    + cond * (TTEMP_LOC(ng) - grid%tempbc(ibc)) - fluxh
      
      !tracer
      trans = grid%difaq * grid%areabc(nc) / grid%distbc(nc)
      r_p(3+(m-1)*grid%ndof) = r_p(3+(m-1)*grid%ndof) &
                + trans * (CCONC_LOC(ng) - grid%concbc(ibc)) - fluxc

    else if(grid%ibndtyp(ibc) == 2) then  ! Constant velocity q, grad T, C = 0

      fluxp = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        q = grid%velocitybc(j,ibc) * grid%density_bc(j) * grid%areabc(nc)
        fluxp = fluxp - q
      enddo

      r_p(1+(m-1)*grid%ndof) = r_p(1+(m-1)*grid%ndof) + fluxp
      
      !heat residual: specified temperature
!     jmp1 = grid%nphase+1+(m-1)*grid%ndof
!     i1 = ithrm_loc_p(ng)
!     cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
!     r_p(jmp1) = r_p(jmp1) + cond * (TTEMP_LOC(ng) - grid%tempbc(ibc)) &
!                           + fluxh
      
      !tracer: specified concentration
!     trans = grid%difaq * grid%areabc(nc) / grid%distbc(nc)
      !check for upstreaming weighting and iface even or odd etc.
      !use zero gradient BC
      
!     r_p(jmp1+1) = r_p(jmp1+1) + trans * (CCONC_LOC(ng) - grid%concbc(ibc)) &
!                               + fluxc

    else if(grid%ibndtyp(ibc) == 3) then  ! fixed p, grad T, C = 0
    
      fluxp = 0.d0
      fluxh = 0.d0
      fluxc = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        D = perm1 / grid%viscosity_bc(j)
        
        !v_darcy is positive for fluid flowing into block
        v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        - grid%fmwh2o * ddensity_loc_p(jng) * grid%gravity * grid%delzbc(nc) &
                 ) / grid%distbc(nc)
        q = v_darcy * grid%areabc(nc)
        
        fluxbc = q * grid%density_bc(j)
        fluxp = fluxp - fluxbc
        
        if (q > 0.d0) then
          fluxh = fluxh + q * grid%density_bc(j) * grid%hh_bc(j)
          fluxc = fluxc + q * grid%concbc(ibc)
        else
          fluxh = fluxh + q * ddensity_loc_p(jng) * hh_loc_p(jng)
          fluxc = fluxc + q * CCONC_LOC(ng)
        endif
      enddo

      r_p(1+(m-1)*grid%ndof) = r_p(1+(m-1)*grid%ndof) + fluxp

      !heat
      r_p(2+(m-1)*grid%ndof) = r_p(2+(m-1)*grid%ndof) - fluxh

      !solute (tracer)
      r_p(3+(m-1)*grid%ndof) = r_p(3+(m-1)*grid%ndof) - fluxc
    endif
  enddo
  
! call VecGetArrayF90(xx, xx_p, ierr) ; CHKERRQ(ierr)
! do n = 1,grid%nlmax
!   print *,n,PPRESSURE(1,n),TTEMP(n),CONC(n), &
!   (r_p(j+(n-1)*grid%ndof),j=1,grid%ndof)
! enddo
! call VecRestoreArrayF90(xx, xx_p, ierr) ; CHKERRQ(ierr)
  
  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%ssat_loc, ssat_loc_p, ierr)
  call VecRestoreArrayF90(grid%sat, sat_p, ierr)
  call VecRestoreArrayF90(grid%xxmol_loc, xxmol_loc_p, ierr)
  call VecRestoreArrayF90(grid%xmol, xmol_p, ierr)
  call VecRestoreArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(grid%h, h_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)
  call VecRestoreArrayF90(grid%accum, accum_p, ierr)

  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)

  call VecRestoreArrayF90(grid%vl, vl_p, ierr)
  
  if (grid%ideriv == 1) then
    call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
    call VecRestoreArrayF90(grid%v_p_loc, v_p_loc_p, ierr)
    call VecRestoreArrayF90(grid%v_t_loc, v_t_loc_p, ierr)
  endif
  
! call VecScale(1.d-2,r)
  
! call VecView(grid%vl,PETSC_VIEWER_STDOUT_WORLD,ierr)

  end subroutine TWOPHResidual

!======================================================================

  subroutine TWOPHJacobian(snes, xx, A, B, flag, grid, ierr)

  use water_eos_module

  implicit none

!#include "include/finclude/petsc.h"
!#include "include/finclude/petscvec.h"
!#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
!#include "include/finclude/petscmat.h"
!#include "finclude/petscda.h"
!#ifdef USE_PETSC216
!#include "include/finclude/petscsles.h"
!#endif
!#include "include/finclude/petscsnes.h"
!#include "include/finclude/petscviewer.h"
!#include "include/finclude/petscsys.h"
!#include "include/finclude/petscis.h"
!#include "include/finclude/petsclog.h"
!#include "pflow_gridtype.h"
  
  SNES, intent(in) :: snes
  Vec, intent(in) :: xx
  Mat, intent(out) :: A, B
  type(pflowGrid), intent(inout) :: grid
  integer, intent(in) :: flag

! external WATEOS, VISW, PSAT
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  integer :: ierr
  integer :: n, ng, nc
  integer :: i, i1, i2, j, jn, jnp1, jnp2, jng, jngp1, jm, jm1, jm2, jm1p1, jm2p1
  integer :: lm1p1, lm2p1, lnp1, lngp1, m, m1, m2, n1, n2, ip1, ip2
  integer :: ibc  ! Index that specifies a boundary condition block.
  real*8 :: elem1, elem2, v_darcy, q, qdtmp1, qdtmp2
  PetscScalar, pointer :: porosity_loc_p(:),volume_p(:),xx_loc_p(:), &
               ddensity_loc_p(:),density_p(:), &
               d_p_loc_p(:),d_t_loc_p(:), &
               hh_loc_p(:),h_p_loc_p(:),h_t_loc_p(:), &
               viscosity_loc_p(:),v_p_loc_p(:),v_t_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               ithrm_loc_p(:)
  real*8 :: dd1, dd2, cond, trans, density_ave, dw_mol, dw_kg, voldt, pvoldt
  real*8 :: daccep, daccet, dtrans, dd, f1, f2, u1, u2
  real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2 
  real*8 :: por1, por2, perm1, perm2, diff
  real*8 :: D  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.

  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
                          
  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
  call VecGetArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
  call VecGetArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecGetArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
  call VecGetArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
  call VecGetArrayF90(grid%v_p_loc, v_p_loc_p, ierr)
  call VecGetArrayF90(grid%v_t_loc, v_t_loc_p, ierr)
  
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  
  !---------------------------------------------------------------------------
  ! Calculate jacobian for accumulation term for local nodes.
  ! Structure:
  !  do j = 1, nphase
  !    do l = 1, nphase
  !      J(j,l) = ... (p,p) <diagonal>
  !    enddo
  !    J(j, Np+1) = ... (p,T)
  !  enddo
  !  do l = 1, nphase
  !    J(Np+1, l) = ... (T,p)
  !  enddo
  !  J(Np+1, Np+1) = ... (T,T)
  !---------------------------------------------------------------------------
  ! Accumulation terms
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    lnp1 = grid%nphase + (ng-1)*grid%ndof ! = grid%nphase+1+(ng-1)*grid%ndof-1
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    !sum over fluid phases
    daccep = 0.d0
    daccet = 0.d0
    do j = 1, grid%nphase
      jnp1 = j + (ng-1)*grid%ndof-1
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.

      !deriv wrt pressure: (p,p)
      elem1 = pvoldt * d_p_loc_p(jng) 
      call MatSetValuesLocal(A,1,jnp1,1,jnp1,elem1,ADD_VALUES,ierr)
       
      !deriv wrt temperature: (p,T)
      elem1 = pvoldt * d_t_loc_p(jng) 
      call MatSetValuesLocal(A,1,jnp1,1,lnp1,elem1,ADD_VALUES,ierr)
      
      !deriv of energy accumulation term wrt p: (T,p)
      daccep = daccep + d_p_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_p_loc_p(jng) - grid%scale

      !deriv of energy accumulation term wrt T: (T,T)
      daccet = daccet + d_t_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_t_loc_p(jng)

!     print *,'TH: ',d_p_loc_p(jn),d_t_loc_p(jng),ddensity_loc_p(jng), &
!     hh_loc_p(jng),h_p_loc_p(jng),h_t_loc_p(jng),grid%scale
    enddo

    !energy eqn. - cross term: (T,p)
    jnp1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
    elem1 = pvoldt * daccep
    call MatSetValuesLocal(A,1,lnp1,1,jnp1,elem1,ADD_VALUES,ierr)

    !energy eqn. - diagonal term: (T,T)
    i = ithrm_loc_p(ng)
    elem1 = pvoldt * daccet + (1.d0-porosity_loc_p(ng))*grid%dencpr(i)*voldt
    call MatSetValuesLocal(A,1,lnp1,1,lnp1,elem1,ADD_VALUES,ierr)
    
    !tracer: (C,C)
    jnp2 = lnp1 + 1
    elem1 = pvoldt
    call MatSetValuesLocal(A,1,jnp2,1,jnp2,elem1,ADD_VALUES,ierr)

  enddo

  !---------------------------------------------------------------------------
  ! Flux terms for interior nodes
  !---------------------------------------------------------------------------
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc)
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1)
    n2 = grid%nG2L(m2)
    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)
    ip2 = grid%iperm2(nc)
    
!   perm1 = perm_loc_p(ip1+3*(m1-1))
!   perm2 = perm_loc_p(ip2+3*(m2-1))
    
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(m1)
    else if (ip1==2) then
      perm1 = perm_yy_loc_p(m1)
    else
      perm1 = perm_zz_loc_p(m1)
    endif
    
    if (ip2 == 1) then
      perm2 = perm_xx_loc_p(m2)
    else if (ip2==2) then
      perm2 = perm_yy_loc_p(m2)
    else
      perm2 = perm_zz_loc_p(m2)
    endif
    
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd
    
    lm1p1 = grid%nphase + (m1-1) * grid%ndof ! = grid%nphase+1 + (m1-1)*grid%ndof-1
    lm2p1 = grid%nphase + (m2-1) * grid%ndof ! = grid%nphase+1 + (m2-1)*grid%ndof-1

    por1 = porosity_loc_p(m1)
    por2 = porosity_loc_p(m2)
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
    dtrans = diff * grid%area(nc)
    
!   print *,'2PH: ',nc,m1,m2,por1,por2,dtrans,diff,grid%difaq

    dfluxt = 0.d0
    dfluxp = 0.d0

    do j = 1, grid%nphase
      jm1p1 = j + (m1-1) * grid%ndof - 1
      jm2p1 = j + (m2-1) * grid%ndof - 1
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

      D1 = perm1 / viscosity_loc_p(jm1)
      D2 = perm2 / viscosity_loc_p(jm2)

      D = (D1 * D2) / (dd2*D1 + dd1*D2)

      density_ave = f2 * ddensity_loc_p(jm1) + f1* ddensity_loc_p(jm2)

      v_darcy = -D * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) & 
                - grid%fmwh2o * density_ave * grid%gravity * grid%delz(nc))

      q = v_darcy * grid%area(nc)

      trans = density_ave * D * grid%area(nc)
!     trans = density_ave * D * grid%area(nc) * (1.d0 - grid%fmwh2o * &
!     (f2 * d_p_loc_p(jm1) + f1 * d_p_loc_p(jm2)) * grid%gravity * grid%delz(nc))
      
      qdtmp1 = q * (f2 * d_t_loc_p(jm1) - v_t_loc_p(jm1)/viscosity_loc_p(jm1) &
                 * density_ave)
      qdtmp2 = q * (f1 * d_t_loc_p(jm2) - v_t_loc_p(jm2)/viscosity_loc_p(jm2) &
                 * density_ave)

      !upstream weighting 
      if (q > 0.d0) then
        u1 = q
        u2 = 0.d0
      else
        u1 = 0.d0
        u2 = q
      endif

      !ignore changes in viscosity for now
      dfluxp1 = dfluxp1 + u1 * (ddensity_loc_p(jm1) * h_p_loc_p(jm1) &
                            + d_p_loc_p(jm1) * hh_loc_p(jm1))
      dfluxp2 = dfluxp2 + u2 * (ddensity_loc_p(jm2) * h_p_loc_p(jm2) &
                            + d_p_loc_p(jm2) * hh_loc_p(jm2))
      dfluxt1 = dfluxt1 + u1 * (ddensity_loc_p(jm1) * h_t_loc_p(jm1) &
                            + d_t_loc_p(jm1) * hh_loc_p(jm1))
      dfluxt2 = dfluxt2 + u2 * (ddensity_loc_p(jm2) * h_t_loc_p(jm2) &
                            + d_t_loc_p(jm2) * hh_loc_p(jm2))

      ! Now add the flux contributions for this phase.
      ! Note that fluxes through a downstream face should be added to the
      ! residual component at the cell, while fluxes through an upstream face 
      ! should be subtracted.  (The divergence gives the net OUTFLOW rate per
      ! unit volume.)  Thus, when working with pressure differences,
      ! (ppressure(jm2) - ppressure(jm1)) should be *subtracted* at the 
      ! upstream node n1 because q = -D*div(P).
      
      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        ! liquid flux terms: (p,p)
        elem1 =  trans
        elem2 = -trans
        call MatSetValuesLocal(A,1,jm1p1,1,jm1p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm1p1,1,jm2p1,elem2,ADD_VALUES,ierr)

        ! liquid flux terms: (p,T)
        elem1 = qdtmp1
        elem2 = qdtmp2
        call MatSetValuesLocal(A,1,jm1p1,1,lm1p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm1p1,1,lm2p1,elem2,ADD_VALUES,ierr)

        ! tracer flux terms: (C,C) -add upstream weighting-
        elem1 = u1 + dtrans
        elem2 = u2 - dtrans
        elem1 = 1.d0
        elem2 = 1.d0
        call MatSetValuesLocal(A,1,lm1p1+1,1,lm1p1+1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,lm1p1+1,1,lm2p1+1,elem2,ADD_VALUES,ierr)
      endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
      ! liquid flux terms: (p,p)
        elem1 = -trans
        elem2 =  trans
        call MatSetValuesLocal(A,1,jm2p1,1,jm1p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm2p1,1,jm2p1,elem2,ADD_VALUES,ierr)

        ! liquid flux terms: (p,T)
        elem1 = -qdtmp1
        elem2 = -qdtmp2
        call MatSetValuesLocal(A,1,jm2p1,1,lm1p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm2p1,1,lm2p1,elem2,ADD_VALUES,ierr)

      ! tracer flux terms: (C,C) -add upstream weighting-
        elem1 = -u1 - dtrans
        elem2 = -u2 + dtrans
        call MatSetValuesLocal(A,1,lm2p1+1,1,lm1p1+1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,lm2p1+1,1,lm2p1+1,elem2,ADD_VALUES,ierr)
      endif
    enddo
    
    !heat flux terms for thermal conduction
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    D = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = D * grid%area(nc)
    
    jm1p1 = m1*grid%ndof - 1
    jm2p1 = m2*grid%ndof - 1
    if (n1 > 0) then  ! If the upstream node is not a ghost node...
      !(T,T)
      elem1 =  cond + dfluxt1
      elem2 = -cond + dfluxt2
      call MatSetValuesLocal(A,1,jm1p1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,jm1p1,1,jm2p1,elem2,ADD_VALUES,ierr)

      ! heat flux terms: (T,p) -add upstream weighting-
      elem1 =  trans * hh_loc_p(jm1) - dfluxp1
      elem2 = -trans * hh_loc_p(jm2) + dfluxp2
      call MatSetValuesLocal(A,1,lm1p1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,lm1p1,1,jm2p1,elem2,ADD_VALUES,ierr)

      ! tracer flux terms: (C,p) -add upstream weighting-
      elem1 =  trans * CCONC_LOC(m1)
      elem2 = -trans * CCONC_LOC(m1)
      call MatSetValuesLocal(A,1,lm1p1+1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,lm1p1+1,1,jm2p1,elem2,ADD_VALUES,ierr)
    endif

    if (n2 > 0) then ! If the downstream node is not a ghost node...
      !(T,T)
      elem1 = -cond - dfluxt1
      elem2 =  cond - dfluxt2
      call MatSetValuesLocal(A,1,jm2p1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,jm2p1,1,jm2p1,elem2,ADD_VALUES,ierr)

      ! heat flux terms: (T,p) -add upstream weighting-
      elem1 = -trans * hh_loc_p(jm1) + dfluxp1
      elem2 =  trans * hh_loc_p(jm2) - dfluxp2
      call MatSetValuesLocal(A,1,lm2p1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,lm2p1,1,jm2p1,elem2,ADD_VALUES,ierr)

      ! tracer flux terms: (C,p) -add upstream weighting-
      elem1 =  trans * CCONC_LOC(m2)
      elem2 = -trans * CCONC_LOC(m2)
      call MatSetValuesLocal(A,1,lm1p1+1,1,jm1p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,lm1p1+1,1,jm2p1,elem2,ADD_VALUES,ierr)
    endif

  enddo
  
  !---------------------------------------------------------------------------
  ! Flux terms for boundary nodes.
  !---------------------------------------------------------------------------

  do nc=1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)

!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1==2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
    
    if(grid%ibndtyp(ibc) == 2) then
      ! solve for pb from Darcy's law given qb /= 0
      grid%tempbc(ibc) = TTEMP_LOC(ng)
    else if(grid%ibndtyp(ibc) == 3) then
      grid%tempbc(ibc) = TTEMP_LOC(ng)
    endif

    do j = 1, grid%nphase
      if (j == grid%jh2o) then
        !pure water eos: Note boundary derivatives not needed!
        call wateos(grid%tempbc(ibc),grid%pressurebc(j,ibc),dw_kg,dw_mol, &
        grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j), &
        grid%h_t_bc(j),grid%scale,ierr)
        grid%density_bc(j) = dw_mol
        call PSAT(grid%tempbc(ibc),sat_pressure,ierr)
        call VISW_noderiv(grid%tempbc(ibc),grid%pressurebc(j,ibc),sat_pressure, &
        grid%viscosity_bc(j),ierr)

!     else if (j == grid%jco2) then

!       add additional fluid phases here ...

      else ! for testing purposes only
        jng = j + (ng-1) * grid%nphase
        grid%viscosity_bc(j) = viscosity_loc_p(jng)
        grid%hh_bc(j) = hh_loc_p(jng)
        grid%density_bc(j) = ddensity_loc_p(jng)
      endif
    enddo

    lngp1 = grid%nphase + (ng-1)*grid%ndof ! = grid%nphase+1 + (ng-1)*grid%ndof-1

    if(grid%ibndtyp(ibc) == 1) then ! Dirichlet BC: fixed p, T, C
    
      dfluxt = 0.d0
      dfluxp = 0.d0
      do j=1, grid%nphase
        jngp1 = j + (ng-1)*grid%ndof - 1
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        D = perm1 / viscosity_loc_p(jm)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        - grid%fmwh2o * ddensity_loc_p(jng) * grid%gravity * grid%delzbc(nc) &
                 ) / grid%distbc(nc)
        q = v_darcy * grid%areabc(nc)

        trans = ddensity_loc_p(jm) * D / grid%distbc(nc) * grid%areabc(nc)

        ! (p,p)
        elem1 = trans
        call MatSetValuesLocal(A,1,jngp1,1,jngp1,elem1,ADD_VALUES,ierr)
      
        if (q < 0.d0) then
          ! (p,T)
          elem1 = q * d_t_loc_p(jng)
          call MatSetValuesLocal(A,1,jngp1,1,lngp1,elem1,ADD_VALUES,ierr)

          dfluxp = dfluxp + q * (ddensity_loc_p(jng) * h_p_loc_p(jng) &
          + d_p_loc_p(jng) * hh_loc_p(jng))
          dfluxt = dfluxt + q * (ddensity_loc_p(jng) * h_t_loc_p(jng) &
          + d_t_loc_p(jng) * hh_loc_p(jng))
        endif

        ! (T,p)
        elem1 = trans * grid%hh_bc(j) - dfluxp
        call MatSetValuesLocal(A,1,lngp1,1,jngp1,elem1,ADD_VALUES,ierr)

        ! (C,p)
        elem1 = trans * CCONC_LOC(ng)
        call MatSetValuesLocal(A,1,lngp1+1,1,jngp1,elem1,ADD_VALUES,ierr)
      
        ! (C,C)
        if (q > 0.d0) then
          elem1 = q
          call MatSetValuesLocal(A,1,lngp1+1,1,lngp1+1,elem1,ADD_VALUES,ierr)
        endif
      enddo
      
      ! (T,T)
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      elem1 = -dfluxt + cond
      call MatSetValuesLocal(A,1,lngp1,1,lngp1,elem1,ADD_VALUES,ierr)
      
    else if(grid%ibndtyp(ibc) == 2) then ! constant velocity q, grad T, C = 0
      
    else if(grid%ibndtyp(ibc) == 3) then ! Dirichlet BC: fixed p, grad T, C = 0
    
      dfluxt = 0.d0
      dfluxp = 0.d0
      do j=1, grid%nphase
        jngp1 = j + (ng-1)*grid%ndof - 1
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        D = perm1 / viscosity_loc_p(jm)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -D * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,ibc) &
        - grid%fmwh2o * ddensity_loc_p(jng) * grid%gravity * grid%delzbc(nc) &
                 ) / grid%distbc(nc)
        q = v_darcy * grid%areabc(nc)

        trans = ddensity_loc_p(jm) * D / grid%distbc(nc) * grid%areabc(nc)

        ! (p,p)
        elem1 = trans
        call MatSetValuesLocal(A,1,jngp1,1,jngp1,elem1,ADD_VALUES,ierr)
      
        if (q < 0.d0) then
          ! (p,T)
          elem1 = q * d_t_loc_p(jng)
          call MatSetValuesLocal(A,1,jngp1,1,lngp1,elem1,ADD_VALUES,ierr)

          dfluxp = dfluxp + q * (ddensity_loc_p(jng) * h_p_loc_p(jng) &
          + d_p_loc_p(jng) * hh_loc_p(jng))
          dfluxt = dfluxt + q * (ddensity_loc_p(jng) * h_t_loc_p(jng) &
          + d_t_loc_p(jng) * hh_loc_p(jng))
        endif

        ! (T,p)
        elem1 = trans * grid%hh_bc(j) - dfluxp
        call MatSetValuesLocal(A,1,lngp1,1,jngp1,elem1,ADD_VALUES,ierr)

        ! (C,p)
        elem1 = trans * CCONC_LOC(ng)
        call MatSetValuesLocal(A,1,lngp1+1,1,jngp1,elem1,ADD_VALUES,ierr)
      
        ! (C,C)
        if (q > 0.d0) then
          elem1 = q
          call MatSetValuesLocal(A,1,lngp1+1,1,lngp1+1,elem1,ADD_VALUES,ierr)
        endif
      enddo
      
      ! (T,T)
      i1 = ithrm_loc_p(ng)
      elem1 = -dfluxt
      call MatSetValuesLocal(A,1,lngp1,1,lngp1,elem1,ADD_VALUES,ierr)
    endif
  enddo

  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
  call VecRestoreArrayF90(grid%d_t_loc, d_t_loc_p, ierr)
  call VecRestoreArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(grid%h_p_loc, h_p_loc_p, ierr)
  call VecRestoreArrayF90(grid%h_t_loc, h_t_loc_p, ierr)
  call VecRestoreArrayF90(grid%v_p_loc, v_p_loc_p, ierr)
  call VecRestoreArrayF90(grid%v_t_loc, v_t_loc_p, ierr)

  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  B = A

! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  end subroutine TWOPHJacobian

end module TWOPH_module

#undef PPRESSURE_LOC
#undef PPRESSURE
#undef PRESSURE
#undef TTEMP_LOC
#undef TTEMP
#undef TEMP
#undef CCONC_LOC
#undef CCONC
#undef CONC

