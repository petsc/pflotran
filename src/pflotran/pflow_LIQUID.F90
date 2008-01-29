
#define PPRESSURE_LOC(j,n) xx_loc_p(j+(n-1)*grid%ndof)
#define PPRESSURE(j,n) xx_p(j+(n-1)*grid%ndof)
#define PRESSURE(j,n) yy_p(j+(n-1)*grid%ndof)

module LIQUID_module

  use pflow_gridtype_module
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

public LiquidResidual, LiquidJacobian

contains

  subroutine LiquidResidual(snes, ppressure, r, grid)

  use water_eos_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(in) :: ppressure
  Vec, intent(out) :: r
  type(pflowGrid), intent(inout) :: grid
    ! I'm not sure if I can pass a pflowGrid object as a PETSc user context.
    ! If not, I'll have to use some trickery to get this information into
    ! this function.
    ! Also, I'm not sure if this needs to be intent(inout) instead of
    ! intent(in).

! external VISW, PSAT, WATEOS
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  PetscInt :: ierr
  PetscInt :: n, ng, nc
  PetscInt :: j, jn, jng, jm, jm1, jm2
  PetscInt :: m, m1, m2, n1, n2, ip1, ip2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
                      ppressure_p(:), ppressure_loc_p(:), pressure_p(:), &
                      density_p(:), ddensity_p(:), ddensity_loc_p(:), &
                      viscosity_p(:), viscosity_loc_p(:), &
                      vl_p(:), &
                      perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
                      ttemp_p(:), d_p_p(:), d_p_loc_p(:)
  PetscReal :: dd1, dd2, flux, fluxbc, gravity, density_ave, v_darcy, voldt, q
  PetscReal :: dd, f1, f2, perm1, perm2
  PetscReal :: D  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: sat_pressure  ! Saturation pressure of water.
  PetscReal :: dw_kg,dw_mol,d_t_dum,hw_dum,h_p_dum,h_t_dum

  !---------------------------------------------------------------------------
  ! Calculate the density and viscosity of water at step k+1 at each local
  ! node.
  !---------------------------------------------------------------------------
  call VecGetArrayF90(ppressure, ppressure_p, ierr) ; CHKERRQ(ierr)
  call VecGetArrayF90(grid%ttemp, ttemp_p, ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%d_p, d_p_p, ierr)
  call VecGetArrayF90(grid%viscosity, viscosity_p, ierr)

  do n = 1, grid%nlmax

!   use molar density rather than mass density
    call wateos(ttemp_p(n), ppressure_p(n), dw_kg, dw_mol, d_p_p(n), &
    d_t_dum, hw_dum, h_p_dum, h_t_dum, grid%scale, ierr)
    ddensity_p(n) = dw_mol

    call PSAT(ttemp_p(n), sat_pressure, ierr)
    call VISW_noderiv(ttemp_p(n),ppressure_p(n),sat_pressure,viscosity_p(n), &
    ierr)
  enddo
  
  call VecRestoreArrayF90(ppressure, ppressure_p, ierr) ; CHKERRQ(ierr)
  call VecRestoreArrayF90(grid%ttemp, ttemp_p, ierr)
  call VecRestoreArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecRestoreArrayF90(grid%d_p, d_p_p, ierr)
  call VecRestoreArrayF90(grid%viscosity, viscosity_p, ierr)

  !---------------------------------------------------------------------------
  ! Now that we have calculated the density and viscosity for all local 
  ! nodes, we can perform the global-to-local scatters.
  !---------------------------------------------------------------------------
  call DAGlobalToLocalBegin(grid%da_nphase_dof, ppressure, INSERT_VALUES, &
                            grid%ppressure_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, ppressure, INSERT_VALUES, &
                          grid%ppressure_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                            grid%d_p_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                          grid%d_p_loc, ierr)
                          
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

  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ppressure_loc, ppressure_loc_p, ierr)
  call VecGetArrayF90(grid%pressure, pressure_p, ierr)
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
  call VecGetArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)

  call VecGetArrayF90(grid%vl, vl_p, ierr)
    
! call VecView(grid%perm,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    voldt = porosity_loc_p(ng) * volume_p(n) / grid%dt 
    do j = 1, grid%nphase
      jn = j + (n-1)*grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
        
#ifdef USE_COMPRESSIBILITY
      r_p(jn) = (ppressure_loc_p(jng) - pressure_p(jn)) * d_p_loc_p(jng) &
                * voldt
#else
      r_p(jn) = (ddensity_loc_p(jng) - density_p(jn)) * voldt
#endif

!     print *,'pflowLIQ: ',n,jn,r_p(jn),ppressure_loc_p(jng),pressure_p(jn), &
!     d_p_loc_p(jng) !,ddensity_loc_p(jng),density_p(jn)
    enddo
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

    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase

      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

      D1 = perm1 / viscosity_loc_p(jm1)
      D2 = perm2 / viscosity_loc_p(jm2)

      D = (D1 * D2) / (dd2*D1 + dd1*D2)

      density_ave = f2 * ddensity_loc_p(jm1) + f1* ddensity_loc_p(jm2)

      v_darcy = -D * (ppressure_loc_p(jm2) - ppressure_loc_p(jm1) & 
                - gravity * density_ave)

      q = v_darcy * grid%area(nc)
      flux = density_ave * q
    
      ! Now add the flux contributions for this phase.
      ! Note that fluxes through a downstream face should be added to the
      ! residual component at the cell, while fluxes through an upstream face 
      ! should be subtracted.  (The divergence gives the net OUTFLOW rate per
      ! unit volume.)  Thus, when working with pressure differences,
      ! (ppressure(jm2) - ppressure(jm1)) should be *subtracted* at the 
      ! upstream node n1 because q = -D*div(P).
    
      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        jn = j + (n1-1)*grid%nphase
        r_p(jn) = r_p(jn) + flux
          ! Remember, the PETSc SNES solver requires the residual, not its
          ! negative!

        ! store velocities at interfaces
        vl_p(ip1+3*(n1-1)) = v_darcy
      endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
        jn = j + (n2-1)*grid%nphase
        r_p(jn) = r_p(jn) - flux
      endif
    enddo
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
    
    gravity = grid%fmwh2o * grid%gravity * grid%delz(nc)
    
    if(grid%ibndtyp(ibc) == 1) then ! Dirichlet BC: fixed p
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now.
        D = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        v_darcy = -D * (ppressure_loc_p(jng) - grid%pressurebc(j,nc) &
            - gravity * ddensity_loc_p(jng))
        q = v_darcy * grid%areabc(nc)

        fluxbc = q * ddensity_loc_p(jng)
        r_p(jm) = r_p(jm) - fluxbc
      enddo
      
    else if(grid%ibndtyp(ibc) == 2) then ! constant velocity q
    
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        v_darcy = grid%velocitybc(j,nc)
        if( mod(grid%iface(ibc), 2) == 0) then
          ! The cell center is upstream of the boundary. 
          r_p(jm) = r_p(jm) + v_darcy * grid%areabc(nc)
        else
          r_p(jm) = r_p(jm) - v_darcy * grid%areabc(nc)
        endif
      enddo
      
    else if(grid%ibndtyp(ibc) == 3) then  ! Dirichlet BC: fixed p (same as case 1)
    
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now.
        D = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        
        v_darcy = -D * (ppressure_loc_p(jng) - grid%pressurebc(j,nc) &
            - gravity * ddensity_loc_p(jng))
        q = v_darcy * grid%areabc(nc)

        fluxbc = q * ddensity_loc_p(jng)
        r_p(jm) = r_p(jm) - fluxbc
      enddo
    endif
    grid%vlbc(nc) = v_darcy
  enddo

  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ppressure_loc, ppressure_loc_p, ierr)
  call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)

  call VecRestoreArrayF90(grid%vl, vl_p, ierr)
  
! call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

  end subroutine LiquidResidual

!======================================================================

  subroutine LiquidJacobian(snes, ppressure, A, B, flag, grid, ierr)

  use water_eos_module
  
  implicit none


  SNES, intent(in) :: snes
  Vec, intent(in) :: ppressure
  Mat, intent(out) :: A, B
  type(pflowGrid), intent(in) :: grid
  PetscInt, intent(in) :: flag

! external COWAT, VISW, PSAT
    ! Legacy functions for calculating density, internal energy, viscosity,
    ! and saturation pressure of water.

  PetscInt :: ierr
  PetscInt :: n, ng, nc
  PetscInt :: j, jn, jng, jm, jm1, jm2
  PetscInt :: m, m1, m2, n1, n2, ip1, ip2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
  PetscReal :: elem1, elem2
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                      ppressure_loc_p(:), pressure_p(:), &
                      ddensity_p(:), ddensity_loc_p(:), density_p(:), &
                      viscosity_loc_p(:), &
                      perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:), &
                      d_p_loc_p(:)
  PetscReal :: dd1, dd2, trans, density_ave, voldt
  PetscReal :: dd, f1, f2, perm1, perm2
  PetscReal :: D  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.
! PetscReal :: sat_pressure  ! Saturation pressure of water.

  call MatZeroEntries(A,ierr)

  call DAGlobalToLocalBegin(grid%da_nphase_dof, ppressure, INSERT_VALUES, &
                            grid%ppressure_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, ppressure, INSERT_VALUES, &
                          grid%ppressure_loc, ierr)
                          
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                            grid%ddensity_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%ddensity, INSERT_VALUES, &
                          grid%ddensity_loc, ierr)
                            
  call DAGlobalToLocalBegin(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                            grid%d_p_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_nphase_dof, grid%d_p, INSERT_VALUES, &
                          grid%d_p_loc, ierr)
                          
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

  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ppressure_loc, ppressure_loc_p, ierr)
  call VecGetArrayF90(grid%pressure, pressure_p, ierr)
  call VecGetArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecGetArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecGetArrayF90(grid%density, density_p, ierr)
  call VecGetArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecGetArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecGetArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecGetArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecGetArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)

  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(grid%d_p_loc, d_p_loc_p, ierr)

  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    voldt = porosity_loc_p(ng)*volume_p(n)/grid%dt
    do j = 1, grid%nphase
      jng = j + (ng-1)*grid%nphase
        ! jn and jng give the index of the term in local non-ghosted and 
        ! local ghosted vectors, respectively, that corresponds to the jth
        ! phase at node n.
        
      elem1 = d_p_loc_p(jng)*voldt
      call MatSetValuesLocal(A,1,jng-1,1,jng-1,elem1,ADD_VALUES,ierr)
    enddo
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

    do j = 1, grid%nphase
      jm1 = j + (m1-1) * grid%nphase
      jm2 = j + (m2-1) * grid%nphase
      
      ! We need to calculate the "diffusion" constant D at the interface;
      ! it is defined at the cell centers.  We use the harmonic mean of the
      ! values from the two cells.

      D1 = perm1 / viscosity_loc_p(jm1)
      D2 = perm2 / viscosity_loc_p(jm2)

      D = (D1 * D2) / (dd2*D1 + dd1*D2)

      density_ave = f2 * ddensity_loc_p(jm1) + f1* ddensity_loc_p(jm2)

      trans = density_ave * D * grid%area(nc)
    
      ! Now add the flux contributions for this phase.
      ! Note that fluxes through a downstream face should be added to the
      ! residual component at the cell, while fluxes through an upstream face 
      ! should be subtracted.  (The divergence gives the net OUTFLOW rate per
      ! unit volume.)  Thus, when working with pressure differences,
      ! (ppressure(jm2) - ppressure(jm1)) should be *subtracted* at the 
      ! upstream node n1 because q = -D*div(P).
    
      if (n1 > 0) then  ! If the upstream node is not a ghost node...
        jn = j + (n1-1)*grid%nphase
!       r_p(jn) = r_p(jn) - trans*(ppressure_loc_p(jm2) - ppressure_loc_p(jm1))
          ! Remember, the PETSc SNES solver requires the residual, not its
          ! negative!
        elem1 =  trans
        elem2 = -trans
        call MatSetValuesLocal(A,1,jm1-1,1,jm1-1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm1-1,1,jm2-1,elem2,ADD_VALUES,ierr)
      endif

      if (n2 > 0) then ! If the downstream node is not a ghost node...
        jn = j + (n2-1)*grid%nphase
!       r_p(jn) = r_p(jn) + trans*(ppressure_loc_p(jm2) - ppressure_loc_p(jm1))
        elem1 = -trans
        elem2 =  trans
        call MatSetValuesLocal(A,1,jm2-1,1,jm1-1,elem1,ADD_VALUES,ierr)
        call MatSetValuesLocal(A,1,jm2-1,1,jm2-1,elem2,ADD_VALUES,ierr)
      endif
    enddo
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
    else if (ip1 == 2) then
      perm1 = perm_yy_loc_p(ng)
    else
      perm1 = perm_zz_loc_p(ng)
    endif
    
    if (grid%ibndtyp(ibc) == 1) then 
    ! Dirichlet BC: fixed p
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        D = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        trans = ddensity_loc_p(jng) * D * grid%areabc(nc)
        elem1 = trans
        call MatSetValuesLocal(A,1,jng-1,1,jng-1,elem1,ADD_VALUES,ierr)
      enddo
      
    else if (grid%ibndtyp(ibc) == 2) then ! constant velocity q
    
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        if(mod(grid%iface(ibc), 2) == 0) then
          ! The cell center is upstream of the boundary. 
!         r_p(jm) = r_p(jm) + grid%velocitybc(j, ibc) * grid%areabc(nc)
!         elem1 = 
!         call MatSetValuesLocal(A,1,jng-1,1,lng-1,elem1,ADD_VALUES,ierr)
        else
!         r_p(jm) = r_p(jm) - grid%velocitybc(j, ibc) * grid%areabc(nc)
!         elem1 = 
!         call MatSetValuesLocal(A,1,jng-1,1,lng-1,elem1,ADD_VALUES,ierr)
        endif
      enddo

    else if (grid%ibndtyp(ibc) == 3) then

    ! Dirichlet BC: fixed p
      do j=1, grid%nphase
        jm = j + (m-1) * grid%nphase
        jng = j + (ng-1) * grid%nphase
        
        ! For the value of the "diffusion" constant at the interface, I just
        ! use its value as defined in the boundary cell.  I don't know if this
        ! is the best thing to do, but it will work for now. (E.g. temperature
        ! could be different at the boundary compared to node center resulting 
        ! in different density, viscosity etc. - pcl
        D = perm1 / viscosity_loc_p(jng) / grid%distbc(nc)
        trans = ddensity_loc_p(jng) * D * grid%areabc(nc) ! need to add gravity
        elem1 = trans
        call MatSetValuesLocal(A,1,jng-1,1,jng-1,elem1,ADD_VALUES,ierr)
      enddo
    endif
  enddo

  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ppressure_loc, ppressure_loc_p, ierr)
  call VecRestoreArrayF90(grid%pressure, pressure_p, ierr)
  call VecRestoreArrayF90(grid%ddensity_loc, ddensity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ddensity, ddensity_p, ierr)
  call VecRestoreArrayF90(grid%density, density_p, ierr)
  call VecRestoreArrayF90(grid%viscosity_loc, viscosity_loc_p, ierr)
  
! call VecRestoreArrayF90(grid%perm_loc, perm_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_xx_loc, perm_xx_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_yy_loc, perm_yy_loc_p, ierr)
  call VecRestoreArrayF90(grid%perm_zz_loc, perm_zz_loc_p, ierr)
  
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(grid%d_p_loc, d_p_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  B = A

! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  end subroutine LiquidJacobian

end module LIQUID_module

!======================================================================

#undef PPRESSURE_LOC
#undef PPRESSURE
#undef PRESSURE

