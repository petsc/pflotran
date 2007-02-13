!======================================================================


!#include "pflowTH.F90"

#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*grid%ndof)
#define PPRESSURE(j,n) xx_p(j+(n-1)*grid%ndof)
#define PRESSURE(j,n) yy_p(j+(n-1)*grid%ndof)
#define TTEMP_LOC(n) xx_loc_p(n*grid%ndof)
#define TTEMP(n) xx_p(n*grid%ndof)
#define TEMP(n) yy_p(n*grid%ndof)

module TH_module

use pflow_gridtype_module

private

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
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

!#include "pflow_gridtype.h"

public THResidual, THJacobian

contains

  subroutine THResidual(snes, xx, r, grid)

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
  integer :: i, i1, i2, j, jn, jng, jm, jm1, jm2, jmu
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2, p1, p2, t1, t2
  integer :: ibc  ! Index that specifies a boundary condition block.
  PetscScalar, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:), &
               density_p(:), ddensity_p(:), ddensity_loc_p(:), &
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
  real*8 :: accum, dd1, dd2, eeng, eng, cond, den, density_ave, &
            hflx, fluxh, fluxe, flux, fluxbc, fluxp, gravity, &
            q, v_darcy, voldt, pvoldt
  real*8 :: dd, f1, f2, perm1, perm2
  real*8 :: Dq, Dk  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  real*8 :: sat_pressure  ! Saturation pressure of water.
  real*8 :: dw_kg, dw_mol

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
        else
          call wateos_noderiv(TTEMP(n),PPRESSURE(j,n),dw_kg,dw_mol, &
          hh_p(jn),grid%scale,ierr)
          call VISW_noderiv(TTEMP(n),PPRESSURE(j,n),sat_pressure, &
          viscosity_p(jn),ierr)
        endif
        ddensity_p(jn) = dw_mol
      enddo

!   add additional fluid phases here ...

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

  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    p1 = 1 + (n-1)*grid%ndof
    t1 = p1 + 1
    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt
    accum = 0.d0
    do j = 1, grid%nphase
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
        
!     liquid phase residual
!     r_p(1+(n-1)*grid%ndof) = porosity_loc_p(ng) * voldt &
!      * (d_p_loc_p(jng) * (PPRESSURE_LOC(j,ng) - PRESSURE(j,n)) &
!      + d_t_loc_p(jng) * (TTEMP_LOC(ng) - TEMP(n)))

      accum = accum + pvoldt * (ddensity_loc_p(jng) - density_p(jn))
    enddo
    
!   pressure residual
    r_p(p1) = accum
    
!   heat residual
    i = ithrm_loc_p(ng)
    j = grid%jh2o
    jn = j+(n-1)*grid%nphase
    jng = j+(ng-1)*grid%nphase
    eeng = ddensity_loc_p(jng)*hh_loc_p(jng)-grid%scale*PPRESSURE_LOC(j,ng)
    eng = density_p(jn) * h_p(jn) - grid%scale*PRESSURE(j,n)
    r_p(t1) = (porosity_loc_p(ng) * (eeng - eng) &
                + (1.d0-porosity_loc_p(ng)) * grid%dencpr(i) &
                * (TTEMP_LOC(ng) - TEMP(n))) * voldt
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
    else if (ip2 == 3) then
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

      ! store velocities at interfaces
!     if (n1 > 0) vl_p(ip1+3*(n1-1)) = v_darcy
      
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
    enddo

    !heat residual
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    Dk = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = Dk * grid%area(nc)
    hflx = cond * (TTEMP_LOC(m2) - TTEMP_LOC(m1))
    fluxe = fluxh - hflx
    
    if (n1 > 0) then  ! If upstream node is not a ghost node...
      r_p(p1) = r_p(p1) + fluxp
      r_p(t1) = r_p(t1) + fluxe
    endif

    if (n2 > 0) then ! If downstream node is not a ghost node...
      r_p(p2) = r_p(p2) - fluxp
      r_p(t2) = r_p(t2) - fluxe
    endif
  enddo

!---------------------------------------------------------------------------
! Flux terms for boundary nodes.
!---------------------------------------------------------------------------

  do nc = 1, grid%nconnbc
        
    gravity = grid%fmwh2o * grid%gravity * grid%delzbc(nc)

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)
    
    p1 = 1 + (m-1)*grid%ndof
    t1 = p1 + 1

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
          call wateos(grid%tempbc(nc),grid%pressurebc(j,nc),dw_kg,dw_mol, &
          grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j), &
          grid%h_t_bc(j),grid%scale,ierr)
          call VISW(grid%tempbc(nc),grid%pressurebc(j,nc),sat_pressure, &
          grid%viscosity_bc(j),grid%v_t_bc(j),grid%v_p_bc(j),ierr)

!         print *,'thc-bc: ',nc,j,m,ng,ibc,grid%nphase,grid%tempbc(nc), &
!         grid%pressurebc(j,nc),dw_mol,dw_kg,grid%d_p_bc(j), &
!         grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j),grid%h_t_bc(j), &
!         sat_pressure,grid%viscosity_bc(j),grid%v_t_bc(j),grid%v_p_bc(j)
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

    if(grid%ibndtyp(ibc) == 1) then ! Dirichlet BC for p, T, C ...
    
      fluxh = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase

        Dq = perm1 / grid%viscosity_bc(j) / grid%distbc(nc)
!       Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * grid%density_bc(j))
        q = v_darcy * grid%areabc(nc)

        fluxbc = q * grid%density_bc(j)
        r_p(p1) = r_p(p1) - fluxbc

        !upstream weighting
        if (q > 0.d0) then ! flow into boundary block
          fluxh = fluxh + q * grid%density_bc(j) * grid%hh_bc(j)
        else
          fluxh = fluxh + q * grid%density_bc(j) * hh_loc_p(jng)
        endif
      enddo

      !heat residual
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      r_p(t1) = r_p(t1) + cond * (TTEMP_LOC(ng) - grid%tempbc(nc)) - fluxh

    else if(grid%ibndtyp(ibc) == 2) then ! constant velocity q, grad T,C = 0
    
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        v_darcy = grid%velocitybc(j,nc)
        
        q = v_darcy * grid%density_bc(j) * grid%areabc(nc)
        r_p(p1) = r_p(p1) - q
      enddo

    else if(grid%ibndtyp(ibc) == 3) then ! fixed p = pb, grad T = 0, Tb = T
    
      fluxh = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        ! v_darcy is positive for fluid flowing into block
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * ddensity_loc_p(jng))
        q = v_darcy * grid%areabc(nc)

        fluxbc = q * ddensity_loc_p(jng)
        r_p(p1) = r_p(p1) - fluxbc

        fluxh = fluxh + q * ddensity_loc_p(jng) * hh_loc_p(jng)
      enddo

      !heat
      r_p(t1) = r_p(t1) - fluxh

    endif
    grid%vlbc(nc) = v_darcy
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

! call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

  end subroutine THResidual

!======================================================================

  subroutine THJacobian(snes, xx, A, B, flag, grid, ierr)

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
  integer :: i, i1, i2, j, jn, jng, jm, jm1, jm2
  integer :: m, m1, m2, mu, n1, n2, ip1, ip2
  integer :: ibc  ! Index that specifies a boundary condition block.
  integer :: p1, t1, p2, t2
  real*8 :: daccep, daccet, elem1, elem2
  real*8 :: sat_pressure, dw_kg, dw_mol
  PetscScalar, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), &
                          ddensity_loc_p(:), density_p(:), &
                          viscosity_loc_p(:), &
                      perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
                          d_p_loc_p(:), &
                          d_t_loc_p(:), &
                          hh_loc_p(:), &
                          h_p_loc_p(:), &
                          h_t_loc_p(:), &
                          v_p_loc_p(:), &
                          v_t_loc_p(:), &
                          ithrm_loc_p(:)
  real*8 :: dd1, dd2, den, cond, trans, trans1, trans2, density_ave, voldt, pvoldt
  real*8 :: dfluxp,dfluxt,dfluxp1,dfluxt1,dfluxp2,dfluxt2,gravity
  real*8 :: dd, f1, f2, v_darcy, perm1, perm2, q, qt1, qt2, qp1, qp2, &
            qdt1, qdt2, hm1, hm2
  real*8 :: Dk, Dq  ! "Diffusion" constant for a phase.
  real*8 :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! real*8 :: sat_pressure  ! Saturation pressure of water.

  call MatZeroEntries(A,ierr)

! Is the following necessary-pcl??? We've already done this in residual call.
  call DAGlobalToLocalBegin(grid%da_ndof, xx, INSERT_VALUES, &
                            grid%xx_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_ndof, xx, INSERT_VALUES, &
                          grid%xx_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)
                            
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                            grid%d_p_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                          grid%d_p_loc, ierr)
                            
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_t, INSERT_VALUES, &
                            grid%d_t_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_t, INSERT_VALUES, &
                          grid%d_t_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
                            grid%viscosity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%viscosity, INSERT_VALUES, &
                          grid%viscosity_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%v_p, INSERT_VALUES, &
                            grid%v_p_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%v_p, INSERT_VALUES, &
                          grid%v_p_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%v_t, INSERT_VALUES, &
                            grid%v_t_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%v_t, INSERT_VALUES, &
                          grid%v_t_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                            grid%hh_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%hh, INSERT_VALUES, &
                          grid%hh_loc, ierr)
                            
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%h_p, INSERT_VALUES, &
                            grid%h_p_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_p, INSERT_VALUES, &
                          grid%h_p_loc, ierr)
                            
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%h_t, INSERT_VALUES, &
                            grid%h_t_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%h_t, INSERT_VALUES, &
                          grid%h_t_loc, ierr)
                          
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
  ! Unknowns: (p,s1,s2,...,sNph-1,T,c1,...,cNp)
  ! Structure:
  !  do j = 1, nphase
  !    do l = 1, nphase
  !      J(j,l) = ... (p,p) <diagonal>
  !    enddo
  !    J(j, Nph+1) = ... (p,T)
  !  enddo
  !  do l = 1, nphase
  !    J(Nph+1, l) = ... (T,p)
  !  enddo
  !  J(Nph+1, Nph+1) = ... (T,T)
  !---------------------------------------------------------------------------
  
  !Accumulation terms
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)

    p1 = (ng-1)*grid%ndof ! = 1 + (ng-1)*grid%ndof-1
    t1 = 1 + p1           ! = 2 + (ng-1)*grid%ndof-1

    voldt = volume_p(n) / grid%dt
    pvoldt = porosity_loc_p(ng) * voldt

    !sum over fluid phases
    daccep = 0.d0
    daccet = 0.d0
    do j = 1, grid%nphase
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.

      !deriv wrt pressure: (p,p)
      elem1 = pvoldt * d_p_loc_p(jng) 
      call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)

      !deriv wrt temperature: (p,T)
      elem1 = pvoldt * d_t_loc_p(jng)
!     print *,'TH: Tp ',n,elem1
      call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)

      !deriv of energy accumulation term: (T,p)
      daccep = daccep + d_p_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_p_loc_p(jng) - grid%scale

      !deriv of energy accumulation term: (T,T)
      daccet = daccet + d_t_loc_p(jng) * hh_loc_p(jng) &
      + ddensity_loc_p(jng) * h_t_loc_p(jng)
        
!     print *,'TH: ',d_p_loc_p(jn),d_t_loc_p(jng),ddensity_loc_p(jng), &
!     hh_loc_p(jng),h_p_loc_p(jng),h_t_loc_p(jng),grid%scale
    enddo

    !energy eqn. - cross term: (T,p)
    elem1 = pvoldt * daccep
    call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)

    !energy eqn. - diagonal term: (T,T)
    i = ithrm_loc_p(ng)
    elem1 = pvoldt * daccet + (1.d0-porosity_loc_p(ng))*grid%dencpr(i)*voldt
    call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
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
      
    gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)

    dfluxp1 = 0.d0
    dfluxp2 = 0.d0
    dfluxt1 = 0.d0
    dfluxt2 = 0.d0

    p1 = (m1-1) * grid%ndof
    p2 = (m2-1) * grid%ndof
    t1 = 1 + p1
    t2 = 1 + p2

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
      else
        mu = m2
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

      ! Now add the flux contributions for this phase.
      ! Note that fluxes through a downstream face should be added to the
      ! residual component at the cell, while fluxes through an upstream face 
      ! should be subtracted.  (The divergence gives the net OUTFLOW rate per
      ! unit volume.)  Thus, when working with pressure differences,
      ! (ppressure(jm2) - ppressure(jm1)) should be *subtracted* at the 
      ! upstream node n1 because q = -D*div(P).
      
      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        ! liquid flux terms: (p,p)
        elem1 = trans1
        elem2 = trans2
!       print *,'TH-pp1',elem1,elem2
        call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,p1,1,p2,elem2,ADD_VALUES,ierr)

        ! liquid flux terms: (p,T)
        elem1 = qdt1
        elem2 = qdt2
!       print *,'TH-pT1',n1,elem1,elem2
        call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,p1,1,t2,elem2,ADD_VALUES,ierr)
      endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
        ! liquid flux terms: (p,p)
        elem1 = -trans1
        elem2 = -trans2
!       print *,'TH-pp2',elem1,elem2
        call MatSetValuesLocal(A,1,p2,1,p1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,p2,1,p2,elem2,ADD_VALUES,ierr)

        ! liquid flux terms: (p,T)
        elem1 = -qdt1
        elem2 = -qdt2
!       print *,'TH-pT2',n2,elem1,elem2
        call MatSetValuesLocal(A,1,p2,1,t1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,p2,1,t2,elem2,ADD_VALUES,ierr)
      endif
    enddo

    !heat flux terms for thermal conduction
    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    Dk = (D1 * D2) / (dd2*D1 + dd1*D2)

    cond = Dk * grid%area(nc)

    if (n1 > 0) then  ! If the upstream node is not a ghost node...
      !(T,T)
      elem1 = dfluxt1 + cond
      elem2 = dfluxt2 - cond
!     print *,'TH-TT1',elem1,elem2
      call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,t1,1,t2,elem2,ADD_VALUES,ierr)

      ! heat flux terms: (T,p)
      elem1 = dfluxp1
      elem2 = dfluxp2
!     print *,'THC-Tp1',nc,elem1,elem2,qp1,qp2,tupstrm
      call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,t1,1,p2,elem2,ADD_VALUES,ierr)
    endif

    if (n2 > 0) then ! If the downstream node is not a ghost node...
      !(T,T)
      elem1 = -(dfluxt1 + cond)
      elem2 = -(dfluxt2 - cond)
!     print *,'TH-TT2',elem1,elem2
      call MatSetValuesLocal(A,1,t2,1,t1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,t2,1,t2,elem2,ADD_VALUES,ierr)

      ! heat flux terms: (T,p) -add upstream weighting-
      elem1 = -dfluxp1
      elem2 = -dfluxp2
!     print *,'TH-Tp2',nc,elem1,elem2,qp1,qp2,tupstrm,dfluxp1,dfluxp2
      call MatSetValuesLocal(A,1,t2,1,p1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,t2,1,p2,elem2,ADD_VALUES,ierr)
    endif

  enddo

  !---------------------------------------------------------------------------
  ! Flux terms for boundary nodes.
  !---------------------------------------------------------------------------

  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)
    
    p1 = (ng-1) * grid%ndof
    t1 = p1 + 1
    
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
        call wateos(grid%tempbc(nc),grid%pressurebc(j,nc),dw_kg,dw_mol, &
        grid%d_p_bc(j),grid%d_t_bc(j),grid%hh_bc(j),grid%h_p_bc(j), &
        grid%h_t_bc(j),grid%scale,ierr)
        grid%density_bc(j) = dw_mol
        call PSAT(grid%tempbc(nc),sat_pressure,ierr)
        call VISW_noderiv(grid%tempbc(nc),grid%pressurebc(j,nc),sat_pressure, &
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
    
    if(grid%ibndtyp(ibc) == 1) then ! Dirichlet BC: fixed p, T
    
      dfluxp = 0.d0
      dfluxt = 0.d0
      do j = 1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        
        Dq = perm1 / grid%viscosity_bc(j) / grid%distbc(nc) * grid%areabc(nc)
!       Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc) * grid%areabc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * grid%density_bc(j))
        q = v_darcy
        
        qp1 = -Dq
!       qt1 = -q * v_t_loc_p(jng) / viscosity_loc_p(jng)
        qt1 = 0.d0

        trans = qp1 * grid%density_bc(j)

        ! (p,p)
        elem1 = -trans
        call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)

        ! (p,T) = 0
!       elem1 = -qt1 * grid%density_bc(j)
!       call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)

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
      enddo

      ! (T,p)
      elem1 = -dfluxp
      call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      
      ! (T,T)
      i1 = ithrm_loc_p(ng)
      cond = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      elem1 = -dfluxt + cond
      call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
      
    else if(grid%ibndtyp(ibc) == 2) then ! constant velocity q, grad T = 0

    else if(grid%ibndtyp(ibc) == 3) then ! Dirichlet BC: fixed p, grad T = 0

      dfluxt = 0.d0
      dfluxp = 0.d0
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        Dq = perm1 / viscosity_loc_p(jng) / grid%distbc(nc) * grid%areabc(nc)
        
        !note: darcy vel. is positive for flow INTO boundary node
        v_darcy = -Dq * (PPRESSURE_LOC(j,ng) - grid%pressurebc(j,nc) &
                  - gravity * ddensity_loc_p(jng))
        q = v_darcy

!       qp1 = -Dq * (1.d0 - gravity * d_p_loc_p(jng))
        qp1 = -Dq
        
        qt1 = -q * v_t_loc_p(jng) / viscosity_loc_p(jng) &
              + Dq * gravity * d_t_loc_p(jng)

!       trans1 = -(qp1 * ddensity_loc_p(jng) + q * d_p_loc_p(jng))
        trans1 = -qp1 * ddensity_loc_p(jng)

        qdt1 = -(qt1 * ddensity_loc_p(jng) + q * d_t_loc_p(jng))

        !(p,p)
        elem1 = trans1
        call MatSetValuesLocal(A,1,p1,1,p1,elem1,ADD_VALUES,ierr)
        
        ! (p,T)
        elem1 = qdt1
        call MatSetValuesLocal(A,1,p1,1,t1,elem1,ADD_VALUES,ierr)
      
!       dfluxp = dfluxp + q * (ddensity_loc_p(jng) * h_p_loc_p(jng) &
!         + d_p_loc_p(jng) * hh_loc_p(jng)) &
!         + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

!       dfluxp = dfluxp + q * ddensity_loc_p(jng) * h_p_loc_p(jng) &
!         + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

        dfluxp = dfluxp + qp1 * ddensity_loc_p(jng) * hh_loc_p(jng)

        dfluxt = dfluxt + q * (ddensity_loc_p(jng) * h_t_loc_p(jng) &
          + d_t_loc_p(jng) * hh_loc_p(jng)) &
          + qt1 * ddensity_loc_p(jng) * hh_loc_p(jng)
      enddo

      ! (T,p)
      elem1 = -dfluxp
      call MatSetValuesLocal(A,1,t1,1,p1,elem1,ADD_VALUES,ierr)
      
      ! (T,T), grad T = 0
      elem1 = -dfluxt
      call MatSetValuesLocal(A,1,t1,1,t1,elem1,ADD_VALUES,ierr)
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
  
  end subroutine THJacobian
  
end module TH_module

#undef PPRESSURE_LOC
#undef PPRESSURE
#undef PRESSURE
#undef TTEMP_LOC
#undef TTEMP
#undef TEMP
