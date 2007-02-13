!======================================================================

!#include "pflowTHC.F90"

#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*grid%ndof)
#define PPRESSURE(j,n)     xx_p(j+(n-1)*grid%ndof)
#define PRESSURE(j,n)      yy_p(j+(n-1)*grid%ndof)
#define TTEMP_LOC(n)       xx_loc_p(2+(n-1)*grid%ndof)
#define TTEMP(n)           xx_p(2+(n-1)*grid%ndof)
#define TEMP(n)            yy_p(2+(n-1)*grid%ndof)
#define CCONC_LOC(n)       xx_loc_p(3+(n-1)*grid%ndof)
#define CCONC(n)           xx_p(3+(n-1)*grid%ndof)
#define CONC(n)            yy_p(3+(n-1)*grid%ndof)

module THC_module

 use pflow_gridtype_module
 
 private

#include "include/finclude/petsc.h"
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

!#include "pflow_gridtype.h"

  public THCResidual, THCJacobian

contains

  subroutine THCResidual(snes, xx, r, grid)
  
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
  integer :: n, ng, nc, nr
  integer :: i, i1, i2, j, jn, jng, jm, jm1, jm2, jmu
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2, p1, p2, t1, t2, c1, c2
  integer :: kk1,kk2,jj1,jj2,ii1,ii2, kk, jj, ii
  integer :: ibc  ! Index that specifies a boundary condition block.
  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
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
  real*8 :: dd1, dd2, diflux, diff, eeng, eng, cond, den, trans, density_ave, &
            hflx, fluxc, fluxe, fluxh, flux, fluxp, gravity, &
            fluxv, fluxbc, q, v_darcy, pvoldt, voldt, accum, eps
  real*8 :: dd, f1, f2, ff, por1, por2, perm1, perm2
  real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol
  real*8 :: tsrc1, qsrc1, qqsrc, csrc1, enth_src
  
  eps = 1.d-6

!  degrees of freedom:
!  p - 1:grid%nphase
!  T - grid%nphase + 1
!  C - grid%nphase + 2
!  grid%ndof = grid%nphase + 2 = 3
  !---------------------------------------------------------------------------
  ! Calculate the density and viscosity of water at step k+1 at each local
  ! node.
  !---------------------------------------------------------------------------

  call VecGetArrayF90(xx, xx_p, ierr); CHKERRQ(ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%hh, hh_p, ierr)
  call VecGetArrayF90(grid%viscosity, viscosity_p, ierr)

  if (grid%ideriv == 1) then
    call VecGetArrayF90(grid%d_p, d_p_p, ierr)
    call VecGetArrayF90(grid%d_t, d_t_p, ierr)
    call VecGetArrayF90(grid%h_p, h_p_p, ierr)
    call VecGetArrayF90(grid%h_t, h_t_p, ierr)
    call VecGetArrayF90(grid%v_p, v_p_p, ierr)
    call VecGetArrayF90(grid%v_t, v_t_p, ierr)
  endif
  
  !get phase properties
  
  do j = 1, grid%nphase
    if (j == grid%jh2o) then
      !pure liquid water eos
      do n = 1, grid%nlmax
        jn = j + (n-1)*grid%nphase
      
!       first arg temp, 2nd pressure, use molar units for density 
!         and other properties

        call PSAT(TTEMP(n), sat_pressure, ierr)
        if (grid%ideriv == 1) then
          call wateos(TTEMP(n),PPRESSURE(j,n),dw_kg,dw_mol,d_p_p(jn), &
          d_t_p(jn),hh_p(jn),h_p_p(jn),h_t_p(jn),grid%scale,ierr)
          call VISW(TTEMP(n),PPRESSURE(j,n),sat_pressure, &
          viscosity_p(jn),v_t_p(jn),v_p_p(jn),ierr)
          
          !units: TTEMP [C], PPRESSURE [Pa]
          !       dw_mol [mol/dm^3], dw_kg [kg/m^3]
          !       hh_p [MJ/kmol], h_t_p [MJ/C/kmol]
          
!         print *,'THC: ',n,TTEMP(n),PPRESSURE(j,n),hh_p(jn),h_t_p(jn), &
!         dw_mol,dw_kg
        else
          call wateos_noderiv(TTEMP(n),PPRESSURE(j,n),dw_kg,dw_mol, &
          hh_p(jn),grid%scale,ierr)
          call VISW_noderiv(TTEMP(n),PPRESSURE(j,n),sat_pressure, &
          viscosity_p(jn),ierr)
        endif
        ddensity_p(jn) = dw_mol ! mol/Liter
      enddo

!   add additional fluid phases here ...

!   else if (j == grid%jgas) then
      !pure gas water + air eos
!     do n = 1, grid%nlmax
!       jn = j + (n-1)*grid%nphase
!       call steameos(TTEMP(n),PPRESSURE(j,n),pa,dg_kg,dg_mol,d_p_p(jn), &
!       d_t_p(jn),hh_p(jn),h_p_p(jn),h_t_p(jn),grid%scale,ierr)
!     enddo
      
!   else if (j == grid%jco2) then

!     do n = 1, grid%nphase
!       ...
!     enddo
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

  !---------------------------------------------------------------------------
  ! Now that we have calculated the density and viscosity for all local 
  ! nodes, we can perform the global-to-local scatters.
  !---------------------------------------------------------------------------
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                            grid%hh_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                          grid%hh_loc, ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
                            grid%viscosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
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

  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)

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

  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
  endif

  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    p1 = 1 + (n-1)*grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt
    accum = 0.d0
    do j = 1, grid%nphase
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
      accum = accum + pvoldt * (ddensity_loc_p(jng) - density_p(jn))
    enddo
    
!   pressure residual
    r_p(p1) = accum
    
!   heat residual
    i = ithrm_loc_p(ng)
    j = grid%jh2o
    jn = j+(n-1)*grid%nphase
    jng = j+(ng-1)*grid%nphase
    ! rho U = rho H - p
    eeng = ddensity_loc_p(jng)*hh_loc_p(jng)-grid%scale*PPRESSURE_LOC(j,ng)
    eng = density_p(jn) * h_p(jn) - grid%scale*PRESSURE(j,n)
    r_p(t1) = pvoldt * (eeng - eng) &
                + (1.d0-porosity_loc_p(ng)) * grid%dencpr(i) * &
                (TTEMP_LOC(ng) - TEMP(n)) * voldt

    !tracer
    r_p(c1) = pvoldt * (CCONC_LOC(ng) - CONC(n)) * grid%ret
    
    !kinetic rate term
    if (grid%rk > 0.d0) then
      grid%rate(n) = 0.d0
      if (grid%ityprxn == 1) then
        if (phis_p(n) > 0.d0 .or. CCONC_LOC(ng)/grid%ceq > 1.d0+eps) then
          grid%rate(n) = -grid%rk * grid%area_var(n) * &
          (1.d0 - CCONC_LOC(ng)/grid%ceq)
          r_p(c1) = r_p(c1) + grid%rate(n) * volume_p(n)
          r_p(t1) = r_p(t1) - grid%delHs * grid%rate(n) * &
          volume_p(n)
!         print *,'resTHC: ',n,ng,phis_p(n),grid%rate(n),grid%area_var(n), &
!         CCONC_LOC(ng),CONC(n)
        endif
      else if (grid%ityprxn == 2) then
        if (phis_p(n) > 0.d0 .or. TTEMP_LOC(ng) < 0.d0+eps) then
          grid%rate(n) = -grid%rk * grid%area_var(n) * TTEMP_LOC(ng)
          r_p(p1) = r_p(p1) + grid%rate(n) * volume_p(n)
          r_p(t1) = r_p(t1) - grid%delHs * grid%rate(n) * &
          volume_p(n)
!         print *,'resTHC: ',n,ng,phis_p(n),grid%rate(n),grid%area_var(n), &
!         CCONC_LOC(ng),CONC(n)
        endif
      endif
    endif
    
  enddo

!---------------------------------------------------------------------------
! Flux terms for interior nodes
!---------------------------------------------------------------------------
  do nc = 1, grid%nconn  ! For each interior connection...
    m1 = grid%nd1(nc) ! ghosted
    m2 = grid%nd2(nc)

    n1 = grid%nG2L(m1) ! = zero for ghost nodes
    n2 = grid%nG2L(m2)
    
    p1 = 1 + (n1-1)*grid%ndof
    p2 = 1 + (n2-1)*grid%ndof
    t1 = p1 + 1
    t2 = p2 + 1
    c1 = t1 + 1
    c2 = t2 + 1
    
    dd1 = grid%dist1(nc)
    dd2 = grid%dist2(nc)
    
    ip1 = grid%iperm1(nc)
    ip2 = grid%iperm2(nc)

!   perm1 = perm_loc_p(ip1+3*(m1-1))
!   perm2 = perm_loc_p(ip2+3*(m2-1))
    
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(m1)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(m1)
    else if (ip1 == 3) then
      perm1 = perm_zz_loc_p(m1)
    endif
    
    if (ip2 == 1) then
      perm2 = perm_xx_loc_p(m2)
    else if (ip2 == 2) then
      perm2 = perm_yy_loc_p(m2)
    else if (ip2 == 3) then
      perm2 = perm_zz_loc_p(m2)
    endif
    
    dd = dd1 + dd2
    f1 = dd1/dd
    f2 = dd2/dd

    fluxp = 0.d0
    fluxh = 0.d0
    fluxv = 0.d0
    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
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

      gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)

      v_darcy = -Dq * (PPRESSURE_LOC(j,m2) - PPRESSURE_LOC(j,m1) & 
                - gravity * density_ave)

      ! store velocities defined at interfaces in PETSc Vec vl at upstream node
      grid%vvl_loc(nc) = v_darcy     ! use for coupling to ptran
      if (n1 > 0) then               ! If the upstream node is not a ghost node...
        vl_p(ip1+3*(n1-1)) = v_darcy ! use for print out of velocity
      endif
      
!     if (grid%t/grid%tconv >= 0.1-1.e-6) &
!     print *,'pflowTHC: ',nc,n1,n2,v_darcy*365.*24.*60.*60.

      q = v_darcy * grid%area(nc)

      flux = density_ave * q

      fluxp = fluxp + flux
      
      !upstream weighting
      if (q > 0.d0) then
        mu = m1
      else
        mu = m2
      endif

      jmu = j + (mu-1) * grid%nphase
      fluxh = fluxh + q * density_ave * hh_loc_p(jmu)
      fluxv = fluxv + q * CCONC_LOC(mu)
    enddo

    !heat & tracer residual
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    Dk = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = Dk * grid%area(nc)
    hflx = cond * (TTEMP_LOC(m2) - TTEMP_LOC(m1))
    fluxe = fluxh - hflx

    por1 = porosity_loc_p(m1)
    por2 = porosity_loc_p(m2)
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
    diflux = diff * grid%area(nc) * (CCONC_LOC(m2) - CCONC_LOC(m1))
    fluxc = grid%fc * fluxv - diflux

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
  enddo

!---------------------------------------------------------------------------
! Flux terms for boundary nodes.
!---------------------------------------------------------------------------

  do nc = 1, grid%nconnbc

    gravity = grid%fmwh2o * grid%gravity * grid%delzbc(nc)

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)
    
    p1 = 1 + (m-1) * grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1

    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)
    
!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
    
    if(grid%ibndtyp(ibc) == 2) then
      ! solve for pb from Darcy's law given qb /= 0
      grid%tempbc(nc) = TTEMP_LOC(ng)
    else if(grid%ibndtyp(ibc) == 3) then
      grid%tempbc(nc) = TTEMP_LOC(ng)
    endif

    do j = 1, grid%nphase
      if (j == grid%jh2o) then
        !pure water eos
        call PSAT(grid%tempbc(nc),sat_pressure,ierr)
        if (grid%ideriv == 1) then
          call wateos(grid%tempbc(nc),grid%pressurebc(j,nc),dw_kg, &
          dw_mol,grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j), &
          grid%h_p_bc(j),grid%h_t_bc(j),grid%scale,ierr)
          
          call VISW(grid%tempbc(nc),grid%pressurebc(j,nc), &
          sat_pressure,grid%viscosity_bc(j), &
          grid%v_t_bc(j),grid%v_p_bc(j),ierr)

!         print *,'thc-bc: ',nc,j,m,ng,ibc,grid%nphase, &
!         grid%tempbc(nc),grid%pressurebc(j,ibc),dw_mol, &
!         dw_kg,grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j), &
!         grid%h_p_bc(j),grid%h_t_bc(j),sat_pressure, &
!         grid%viscosity_bc(j),grid%v_t_bc(j),grid%v_p_bc(j)
        else
          call wateos_noderiv(grid%tempbc(nc),grid%pressurebc(j,nc), &
          dw_kg,dw_mol,grid%hh_bc(j),grid%scale,ierr)
          
          call VISW_noderiv(grid%tempbc(nc),grid%pressurebc(j,nc), &
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

        Dq = perm1 / grid%viscosity_bc(j) / grid%distbc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * grid%density_bc(j))
        q = v_darcy * grid%areabc(nc)
      
!       print *,'pflowTHC-1: ',nc,m,ng,ibc,v_darcy,grid%ibndtyp(ibc)

        fluxbc = q * grid%density_bc(j)
        fluxp = fluxp - fluxbc

        !upstream weighting
        if (q > 0.d0) then
          fluxh = fluxh + q * grid%density_bc(j) * grid%hh_bc(j)
          fluxc = fluxc + q * grid%concbc(nc)       ! note: need to change 
        else                          ! definition of C for multiple phases
          fluxh = fluxh + q * grid%density_bc(j) * hh_loc_p(jng)
          fluxc = fluxc + q * CCONC_LOC(ng) 
        endif
      enddo
        
      r_p(p1) = r_p(p1) + fluxp

      !heat residual
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - grid%tempbc(nc)) - fluxh
      
      !tracer
      trans = porosity_loc_p(ng) * grid%difaq * grid%areabc(nc) / grid%distbc(nc)
      r_p(c1) = r_p(c1) + trans * (CCONC_LOC(ng) - grid%concbc(nc)) &
                - grid%fc * fluxc
                
!     print *,'THC: ',nc,ng,ibc,q,trans,grid%concbc(ibc),CCONC_LOC(ng)

    else if(grid%ibndtyp(ibc) == 2) then ! Constant velocity q, grad T,C = 0

      fluxp = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        v_darcy = grid%velocitybc(j,nc)
        
        q = v_darcy * grid%density_bc(j) * grid%areabc(nc)
        fluxp = fluxp - q
      enddo

      r_p(p1) = r_p(p1) + fluxp
      
      !heat residual: specified temperature
!     i1 = ithrm_loc_p(ng)
!     cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
!     r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - grid%tempbc(ibc)) &
!                           + fluxh
      
      !tracer: specified concentration
!     trans = grid%difaq * grid%areabc(nc) / grid%distbc(nc)
      !check for upstreaming weighting and iface even or odd etc.
      !use zero gradient BC
      
!     r_p(c1) = r_p(c1) + trans * (CCONC_LOC(ng) - grid%concbc(ibc)) &
!                               + fluxc

    else if(grid%ibndtyp(ibc) == 3) then  ! fixed p, grad T, C = 0
    
      fluxp = 0.d0
      fluxh = 0.d0
      fluxc = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        !v_darcy is positive for fluid flowing into block
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * ddensity_loc_p(jng))
        q = v_darcy * grid%areabc(nc)
      
!       print *,'pflowTHC-3: ',nc,m,ng,ibc,v_darcy,grid%ibndtyp(ibc)
        
        fluxbc = q * ddensity_loc_p(jng)
        fluxp = fluxp - fluxbc
        
        !upstream weighting not needed: c_N = c_b
        fluxh = fluxh + q * ddensity_loc_p(jng) * hh_loc_p(jng)
        fluxc = fluxc + q * CCONC_LOC(ng)

        !upstream weighting
!       if (q > 0.d0) then
!         fluxh = fluxh + q * grid%density_bc(j) * grid%hh_bc(j)
!         fluxc = fluxc + q * grid%concbc(ibc)       ! note: need to change 
!       else                          ! definition of C for multiple phases
!         fluxh = fluxh + q * grid%density_bc(j) * hh_loc_p(jng)
!         fluxc = fluxc + q * CCONC_LOC(ng) 
!       endif
      enddo

      !mass (pressure)
      r_p(p1) = r_p(p1) + fluxp

      !heat
      r_p(t1) = r_p(t1) - fluxh

      !solute (tracer)
      r_p(c1) = r_p(c1) - grid%fc * fluxc
    endif
    grid%vvlbc(nc) = v_darcy
  enddo
  
! add source/sink terms

  do nr = 1, grid%nblksrc
      
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
        goto 10
      else if (grid%timesrc(i,nr) > grid%t) then
        ff = grid%timesrc(i,nr)-grid%timesrc(i-1,nr)
        f1 = (grid%t - grid%timesrc(i-1,nr))/ff
        f2 = (grid%timesrc(i,nr)-grid%t)/ff
        tsrc1 = f1*grid%tempsrc(i,nr) + f2*grid%tempsrc(i-1,nr)
        qsrc1 = f1*grid%qsrc(i,nr) + f2*grid%qsrc(i-1,nr)
        csrc1 = f1*grid%csrc(i,nr) + f2*grid%csrc(i-1,nr)
        goto 10
      endif
    enddo
 10 continue
    
!   print *,'pflowTHC: ', grid%myrank,i,grid%timesrc(i,nr), &
!   grid%timesrc(i-1,nr),grid%t,f1,f2,ff,qsrc1,csrc1
 
    qsrc1 = qsrc1 / grid%fmwh2o

    if (qsrc1 > 0.d0) then ! injection
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
              n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
              ng = grid%nL2G(n)
              p1 = 1+(n-1)*grid%ndof
              t1 = p1 + 1
              c1 = t1 + 1
              call wateos_noderiv(tsrc1,PPRESSURE_LOC(grid%jh2o,ng), &
              dw_kg,dw_mol,enth_src,grid%scale,ierr)
              qqsrc = qsrc1/dw_mol
              r_p(p1) = r_p(p1) - qsrc1
              r_p(t1) = r_p(t1) - qsrc1*enth_src
              r_p(c1) = r_p(c1) - qqsrc*csrc1

!             print *,'pflowTHC: ',nr,n,ng,qsrc1,dw_mol*grid%fmwh2o, &
!             qqsrc,csrc1,r_p(c1)
          enddo
        enddo
      enddo
        
    else if (qsrc1 < 0.d0) then ! withdrawal
      
      do kk = kk1, kk2
        do jj = jj1, jj2
          do ii = ii1, ii2
              n = ii+(jj-1)*grid%nlx+(kk-1)*grid%nlxy
              ng = grid%nL2G(n)
              p1 = 1+(n-1)*grid%ndof
              t1 = p1 + 1
              c1 = t1 + 1
              qqsrc = qsrc1/ddensity_loc_p(ng)
              enth_src = hh_loc_p(ng)
              r_p(p1) = r_p(p1) - qsrc1
              r_p(t1) = r_p(t1) - qsrc1*enth_src
              r_p(c1) = r_p(c1) - qqsrc*CCONC_LOC(ng)
          enddo
        enddo
      enddo
    endif
  enddo
  
  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%yy, yy_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%hh_loc, hh_loc_p, ierr)
  call VecRestoreArrayF90(grid%h, h_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)

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
  
  if (grid%rk > 0.d0) then
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
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
  type(pflowGrid), intent(inout) :: grid
  integer, intent(out) :: flag

! external WATEOS, VISW, PSAT
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  integer :: ierr
  integer :: n, ng, nc
  integer :: i, i1, i2, j, jn, jng, jm, jm1, jm2
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2
  integer :: p1,p2,t1,t2,c1,c2
  integer :: ibc  ! Index that specifies a boundary condition block.
  real*8 :: elempp, elempt, elem1, elem2, v_darcy, q, qp1, qp2, qt1, qt2, eps
  PetscScalar, pointer :: porosity_loc_p(:),volume_p(:),xx_loc_p(:), &
               ddensity_loc_p(:), phis_p(:), &
               d_p_loc_p(:),d_t_loc_p(:), &
               hh_loc_p(:),h_p_loc_p(:),h_t_loc_p(:), &
               viscosity_loc_p(:),v_p_loc_p(:),v_t_loc_p(:), &
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
               ithrm_loc_p(:)
  real*8 :: cond, trans, trans1, trans2, gravity, &
            density_ave, dw_mol, dw_kg, den, voldt, pvoldt
  real*8 :: daccep, daccet, dtrans, dd1, dd2, dd, f1, f2, u1, u2, qdt1, qdt2
  real*8 :: dfluxp, dfluxt, dfluxp1, dfluxt1, dfluxp2, dfluxt2, cupstrm
  real*8 :: por1, por2, perm1, perm2, diff, hm1, hm2
  real*8 :: Dk, Dq           ! "Diffusion" constant for a phase.
  real*8 :: D1, D2, Dk1, Dk2 ! "Diffusion" constants upstream and downstream from a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  
  real*8 :: blkmat1(3,3),blkmat2(3,3)
  
  eps = 1.d-6

  flag = SAME_NONZERO_PATTERN

  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
                          
  call VecGetArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
! call VecGetArrayF90(grid%density, density_p, ierr)
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
  
  if (grid%rk > 0.d0) then
    call VecGetArrayF90(grid%phis,phis_p,ierr)
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

    p1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
    t1 = p1 + 1           ! = 2 + (ng-1)*grid%ndof-1
    c1 = t1 + 1           ! = 3 + (ng-1)*grid%ndof-1

    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    !sum over fluid phases
    elempp = 0.d0
    elempt = 0.d0
    daccep = 0.d0
    daccet = 0.d0
    do j = 1, grid%nphase
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.

      !deriv wrt pressure: (p,p)
      elempp = elempp + pvoldt * d_p_loc_p(jng) 
       
      !deriv wrt temperature: (p,T)
      elempt = elempt + pvoldt * d_t_loc_p(jng) 
      
      !deriv of energy accumulation term wrt p: (T,p)
      daccep = daccep + d_p_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_p_loc_p(jng) - grid%scale

      !deriv of energy accumulation term wrt T: (T,T)
      daccet = daccet + d_t_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_t_loc_p(jng)

!     print *,'TH: ',d_p_loc_p(jn),d_t_loc_p(jng),ddensity_loc_p(jng), &
!     hh_loc_p(jng),h_p_loc_p(jng),h_t_loc_p(jng),grid%scale
    enddo
    
    !(p,p)
    if (grid%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,p1,1,p1,elempp,ADD_VALUES,ierr)
    else
      blkmat1(1,1) = elempp
    endif

    !(p,T)
    if (grid%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,p1,1,t1,elempt,ADD_VALUES,ierr)
    else
      blkmat1(1,2) = elempt
    endif

    !energy eqn. - cross term: (T,p)
    elem1 = pvoldt * daccep
    if (grid%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(2,1) = elem1
    endif

    !energy eqn. - diagonal term: (T,T)
    i = ithrm_loc_p(ng)
    elem1 = pvoldt * daccet + (1.d0-porosity_loc_p(ng))*grid%dencpr(i)*voldt
    if (grid%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(2,2) = elem1
    endif
    
    !tracer: (C,C)
    elem1 = pvoldt * grid%ret
    if (grid%iblkfmt == 0) then
      call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
    else
      blkmat1(3,3) = elem1
    endif

    !kinetic rate term
    if (grid%rk > 0.d0) then
      if (grid%ityprxn == 1) then
        if (phis_p(n) > 0.d0 .or. CCONC_LOC(ng)/grid%ceq > 1.d0+eps) then
          !(C,C)
          elem1 = grid%rk * grid%area_var(n) / grid%ceq * volume_p(n)
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = blkmat1(3,3) + elem1
          endif

          !(T,C)
          elem1 = -grid%rk * grid%area_var(n)/grid%ceq * grid%delHs * &
          volume_p(n)
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,t1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(2,3) = elem1
          endif
        endif
      else if (grid%ityprxn == 2) then
        if (phis_p(n) > 0.d0 .or. TTEMP_LOC(ng) < 0.d0+eps) then
          !(p,T)
          elem1 = grid%rk * grid%area_var(n) * volume_p(n)
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(1,2) = blkmat1(1,2) + elem1
          endif

          !(T,T)
          elem1 = -grid%rk * grid%area_var(n) * grid%delHs * &
          volume_p(n)
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(2,2) = blkmat1(2,2) + elem1
          endif
        endif
      endif
    endif
    if (grid%iblkfmt == 1) then
      call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat1, &
      ADD_VALUES,ierr)
    endif
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
    diff = (por1 * por2) / (dd2*por1 + dd1*por2) * grid%difaq
    dtrans = diff * grid%area(nc)
      
    gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)
    
    dfluxp1 = 0.d0
    dfluxp2 = 0.d0
    dfluxt1 = 0.d0
    dfluxt2 = 0.d0
    
    p1 = (m1-1) * grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1
    p2 = (m2-1) * grid%ndof
    t2 = p2 + 1
    c2 = t2 + 1

    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
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

      q = v_darcy * grid%area(nc)
      
      qp1 =  Dq * grid%area(nc) * (1.d0 + gravity * f2 * d_p_loc_p(jm1))
      qp2 = -Dq * grid%area(nc) * (1.d0 - gravity * f1 * d_p_loc_p(jm2))

      qt1 = -q * dd1 * perm2 / den * v_t_loc_p(jm1) &
            + Dq * gravity * f2 * d_t_loc_p(jm1) * grid%area(nc)
      qt2 = -q * dd2 * perm1 / den * v_t_loc_p(jm2) &
            + Dq * gravity * f1 * d_t_loc_p(jm2) * grid%area(nc)

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
    Dk1 = grid%ckwet(i1)
    Dk2 = grid%ckwet(i2)

    Dk = (Dk1 * Dk2) / (dd2*Dk1 + dd1*Dk2)

    cond = Dk * grid%area(nc)
    
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
        if (grid%iblkfmt == 0) then
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
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p1,1,t2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,2) = elem1
          blkmat2(1,2) = elem2
        endif

        ! tracer flux terms: (C,C)
        elem1 = grid%fc * u1 + dtrans
        elem2 = grid%fc * u2 - dtrans
!       print *,'THC-CC1',elem1,elem2
        if (grid%iblkfmt == 0) then
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
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t1,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
        blkmat2(2,2) = elem2
      endif

      ! heat flux terms: (T,p)
      elem1 = dfluxp1
      elem2 = dfluxp2
!     print *,'THC-Tp1',nc,elem1,elem2,qp1,qp2,tupstrm
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t1,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
        blkmat2(2,1) = elem2
      endif

      ! tracer flux terms: (C,p)
      elem1 = grid%fc * qp1 * cupstrm
      elem2 = grid%fc * qp2 * cupstrm
!     print *,'THC-Cp1',elem1,elem2
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c1,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,1) = elem1
        blkmat2(3,1) = elem2
      endif

      ! tracer flux terms: (C,T)
      elem1 = grid%fc * qt1 * cupstrm
      elem2 = grid%fc * qt2 * cupstrm
!     print *,'THC-CT1',elem1,elem2,qt1,qt2,cupstrm
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c1,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c1,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,2) = elem1
        blkmat2(3,2) = elem2
      endif
      if (grid%iblkfmt == 1) then
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
        if (grid%iblkfmt == 0) then
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
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p2,1,t1,elem1,ADD_VALUES,ierr)
          call MatSetValuesLocal(A,1,p2,1,t2,elem2,ADD_VALUES,ierr)
        else
          blkmat1(1,2) = elem1
          blkmat2(1,2) = elem2
        endif

      ! tracer flux terms: (C,C)
        elem1 = -(grid%fc * u1 + dtrans)
        elem2 = -(grid%fc * u2 - dtrans)
!       print *,'THC-CC2',elem1,elem2
        if (grid%iblkfmt == 0) then
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
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t2,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t2,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
        blkmat2(2,2) = elem2
      endif

      ! heat flux terms: (T,p)
      elem1 = -dfluxp1
      elem2 = -dfluxp2
!     print *,'TH-Tp2',nc,elem1,elem2,qp1,qp2,tupstrm,dfluxp1,dfluxp2
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t2,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,t2,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
        blkmat2(2,1) = elem2
      endif

      ! tracer flux terms: (C,p)
      elem1 = -grid%fc * qp1 * cupstrm
      elem2 = -grid%fc * qp2 * cupstrm
!     print *,'THC-Cp2',elem1,elem2
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c2,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c2,1,p2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,1) = elem1
        blkmat2(3,1) = elem2
      endif

      ! tracer flux terms: (C,T)
      elem1 = -grid%fc * qt1 * cupstrm
      elem2 = -grid%fc * qt2 * cupstrm
!     print *,'THC-CT2',elem1,elem2,qt1,qt2,cupstrm
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,c2,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,c2,1,t2,elem2,ADD_VALUES,ierr)
      else
        blkmat1(3,2) = elem1
        blkmat2(3,2) = elem2
      endif
      if (grid%iblkfmt == 1) then
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

  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)
    
    p1 = (ng-1)*grid%ndof
    t1 = p1 + 1
    c1 = t1 + 1

    ibc = grid%ibconn(nc)
    ip1 = grid%ipermbc(nc)

!   perm1 = perm_loc_p(ip1+3*(ng-1))
    if (ip1 == 1) then
      perm1 = perm_xx_loc_p(ng)
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
        
    gravity = grid%fmwh2o * grid%gravity * grid%delzbc(nc)
    
    if(grid%ibndtyp(ibc) == 2) then
      ! solve for pb from Darcy's law given qb /= 0
      grid%tempbc(nc) = TTEMP_LOC(ng)
    else if(grid%ibndtyp(ibc) == 3) then
      grid%tempbc(nc) = TTEMP_LOC(ng)
    endif

    do j = 1, grid%nphase
      if (j == grid%jh2o) then
        !pure water eos
        call wateos(grid%tempbc(nc),grid%pressurebc(j,nc),dw_kg, &
        dw_mol,grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j), &
        grid%h_p_bc(j),grid%h_t_bc(j),grid%scale,ierr)
        grid%density_bc(j) = dw_mol
        call PSAT(grid%tempbc(nc),sat_pressure,ierr)
        call VISW_noderiv(grid%tempbc(nc),grid%pressurebc(j,nc), &
        sat_pressure,grid%viscosity_bc(j),ierr)

!     else if (j == grid%jco2) then

!       add additional fluid phases here ...

      else ! for testing purposes only
        jng = j + (ng-1) * grid%nphase
        grid%viscosity_bc(j) = viscosity_loc_p(jng)
        grid%hh_bc(j) = hh_loc_p(jng)
        grid%density_bc(j) = ddensity_loc_p(jng)
      endif
    enddo

    if(grid%ibndtyp(ibc) == 1) then 
    
    ! Dirichlet BC: fixed p, T, C
    
      dfluxp = 0.d0
      dfluxt = 0.d0
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I 
        ! just use its value as defined in the boundary cell.  I don't 
        ! know if this is the best thing to do, but it will work for now. 
        ! (E.g. temperature could be different at the boundary compared  
        ! to node center resulting in different density, viscosity etc. 
        ! - pcl)
        Dq = perm1 / grid%viscosity_bc(j) / grid%distbc(nc) * &
        grid%areabc(nc)

        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * grid%density_bc(j))
        q = v_darcy
        
        qp1 = -Dq
!       qt1 = -q * v_t_loc_p(jng) / viscosity_loc_p(jng)
        qt1 = 0.d0

        trans = qp1 * grid%density_bc(j) ! dR/dp = dq/dp rho

        ! (p,p)
        elem1 = -trans
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
        endif

        ! (p,T) = 0
!       elem1 = -qt1 * grid%density_bc(j)
!       if (grid%iblkfmt == 0) then
!         call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
!       else
!         blkmat1(1,2) = elem1
!       endif

        if (q > 0.d0) then
          dfluxp = dfluxp + qp1 * grid%density_bc(j) * grid%hh_bc(j)
          dfluxt = dfluxt + qt1 * grid%density_bc(j) * grid%hh_bc(j)
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
          elem1 = -grid%fc * qp1 * grid%concbc(nc)
        else
          elem1 = -grid%fc * qp1 * CCONC_LOC(ng)
        endif
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif

        ! (C,T)
        if (q > 0.d0) then
          elem1 = -grid%fc * qt1 * grid%concbc(nc)
        else
          elem1 = -grid%fc * qt1 * CCONC_LOC(ng)
        endif
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif
      
        ! (C,C)
        trans = porosity_loc_p(ng) * grid%difaq * grid%areabc(nc) / &
                grid%distbc(nc)
        if (q < 0.d0) then
          elem1 = -grid%fc * q + trans ! check sign of q
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = elem1
          endif
        else
          elem1 = trans
          if (grid%iblkfmt == 0) then
            call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
          else
            blkmat1(3,3) = elem1
          endif
        endif
      enddo

      ! (T,p)
      elem1 = -dfluxp
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
      endif
      
      ! (T,T)
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      elem1 = -dfluxt + cond
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (grid%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif
      
    else if(grid%ibndtyp(ibc) == 2) then 

    ! constant velocity q, grad T, C = 0
      
    else if(grid%ibndtyp(ibc) == 3) then 
    
    ! Dirichlet BC: fixed p, grad T, C = 0
    
      dfluxt = 0.d0
      dfluxp = 0.d0
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc) * &
        grid%areabc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
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
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(1,1) = elem1
        endif
      
        ! (p,T)
        elem1 = qdt1
        if (grid%iblkfmt == 0) then
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
        elem1 = -grid%fc * qp1 * CCONC_LOC(ng)
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,p1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,1) = elem1
        endif

        ! (C,T)
        elem1 = -grid%fc * qt1 * CCONC_LOC(ng)
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,t1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,2) = elem1
        endif
      
        ! (C,C), grad C = 0
        elem1 = -grid%fc * q
        if (grid%iblkfmt == 0) then
          call MatSetValuesLocal(A,1,c1,1,c1,elem1,ADD_VALUES,ierr)
        else
          blkmat1(3,3) = elem1
        endif
      enddo

      ! (T,p)
      elem1 = -dfluxp
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,1) = elem1
      endif
      
      ! (T,T), grad T = 0
      elem1 = -dfluxt
      if (grid%iblkfmt == 0) then
        call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      else
        blkmat1(2,2) = elem1
      endif

      if (grid%iblkfmt == 1) then
        call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
        blkmat1,ADD_VALUES,ierr)
      endif
    endif
  enddo

  call VecRestoreArrayF90(grid%xx_loc, xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
! call VecRestoreArrayF90(grid%density, density_p, ierr)
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
  
  if (grid%rk > 0.d0) then
    call VecRestoreArrayF90(grid%phis,phis_p,ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  B = A

! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  end subroutine THCJacobian
  
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

