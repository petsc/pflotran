
#define TTEMP_LOC(n) xx_loc_p(n*grid%ndof)
#define TTEMP(n) xx_p(n*grid%ndof)
#define TEMP(n) yy_p(n*grid%ndof)

module COND_module

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

public CondResidual, CondJacobian

contains

  subroutine CondResidual(snes, ttemp, r, grid)
  
  implicit none
  
  SNES, intent(in) :: snes
  Vec, intent(in) :: ttemp
  Vec, intent(out) :: r
  type(pflowGrid), intent(in) :: grid
    ! I'm not sure if I can pass a pflowGrid object as a PETSc user context.
    ! If not, I'll have to use some trickery to get this information into
    ! this function.
    ! Also, I'm not sure if this needs to be intent(inout) instead of
    ! intent(in).

  PetscInt :: ierr
  PetscInt :: n, ng, nc
  PetscInt :: i, i1, i2
  PetscInt :: m, m1, m2, n1, n2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
                          ttemp_loc_p(:), temp_p(:), &
                          ithrm_loc_p(:)
  PetscReal :: dd1, dd2, flux, trans, voldt
  PetscReal :: D  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.

  !---------------------------------------------------------------------------
  ! Now that we have calculated the density and viscosity for all local 
  ! nodes, we can perform the global-to-local scatters.
  !---------------------------------------------------------------------------
  call DAGlobalToLocalBegin(grid%da_1_dof, ttemp, INSERT_VALUES, &
                            grid%ttemp_loc, ierr)
  call DAGlobalToLocalEnd(grid%da_1_dof, ttemp, INSERT_VALUES, &
                          grid%ttemp_loc, ierr)

  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%ttemp_loc, ttemp_loc_p, ierr)
  call VecGetArrayF90(grid%temp, temp_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  call VecGetArrayF90(r, r_p, ierr)
  
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  
  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    i = ithrm_loc_p(ng)
    voldt = volume_p(n) / grid%dt
    r_p(n) = (1.d0-porosity_loc_p(ng)) * grid%dencpr(i) * &
             (ttemp_loc_p(ng) - temp_p(n)) * voldt
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

    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    D = (D1 * D2) / (dd2*D1 + dd1*D2)

    trans = D * grid%area(nc)
    flux = trans*(ttemp_loc_p(m2) - ttemp_loc_p(m1))
    
    if (n1 > 0) then  ! If upstream node is not a ghost node...
      r_p(n1) = r_p(n1) - flux
    endif

    if (n2 > 0) then ! If downstream node is not a ghost node...
      r_p(n2) = r_p(n2) + flux
    endif
  enddo
  
  !---------------------------------------------------------------------------
  ! Flux terms for boundary nodes.
  !---------------------------------------------------------------------------

  do nc = 1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    ibc = grid%ibconn(nc)
    
    if(grid%ibndtyp(ibc) == 1) then  ! Dirichlet BC: Tb
      
      i1 = ithrm_loc_p(ng)
      trans = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      
      r_p(m) = r_p(m) + trans * (ttemp_loc_p(ng) - grid%tempbc(ibc))
      
    else if(grid%ibndtyp(ibc) == 2) then  ! constant heat flux q...
      
    else if(grid%ibndtyp(ibc) == 3) then  ! grad T = 0
!     do nothing!
    endif
  enddo

  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%ttemp_loc, ttemp_loc_p, ierr)
  call VecRestoreArrayF90(grid%temp, temp_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  call VecRestoreArrayF90(r, r_p, ierr)

  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  
! call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

  end subroutine CondResidual

!======================================================================

  subroutine CondJacobian(snes, ttemp, A, B, flag, grid, ierr)
  
  implicit none

#include "include/finclude/petscis.h90"

  SNES, intent(in) :: snes
  Vec, intent(in) :: ttemp
  Mat, intent(out) :: A, B
  type(pflowGrid), intent(in) :: grid
  PetscInt, intent(in) :: flag

  PetscInt :: ierr
  PetscInt :: n, ng, nc
  PetscInt :: i, i1, i2
  PetscInt :: m, m1, m2, n1, n2
  PetscInt :: ibc  ! Index that specifies a boundary condition block.
  PetscReal :: elem1, elem2
  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          ithrm_loc_p(:)
  PetscReal :: dd1, dd2, trans
  PetscReal :: D  ! "Diffusion" constant for a phase.
  PetscReal :: D1, D2  ! "Diffusion" constants upstream and downstream of a face.

  call MatZeroEntries(A,ierr)

  call VecGetArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecGetArrayF90(grid%volume, volume_p, ierr)
  
  call VecGetArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)
  
  !---------------------------------------------------------------------------
  ! Calculate accumulation term for interior and exterior (ghosted) nodes.
  !---------------------------------------------------------------------------
  do n = 1, grid%nlmax  ! For each local node do...
    ng = grid%nL2G(n)
    i = ithrm_loc_p(ng)
    elem1 = (1.d0-porosity_loc_p(ng)) * grid%dencpr(i) / grid%dt * volume_p(n)
    call MatSetValuesLocal(A,1,ng-1,1,ng-1,elem1,ADD_VALUES,ierr)
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

    i1 = ithrm_loc_p(m1)
    i2 = ithrm_loc_p(m2)
    D1 = grid%ckwet(i1)
    D2 = grid%ckwet(i2)

    D = (D1 * D2) / (dd2*D1 + dd1*D2)

    trans = D * grid%area(nc)

    if (n1 > 0) then  ! If the upstream node is not a ghost node...
      elem1 =  trans
      elem2 = -trans
      call MatSetValuesLocal(A,1,m1-1,1,m1-1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,m1-1,1,m2-1,elem2,ADD_VALUES,ierr)
    endif

    if (n2 > 0) then ! If the downstream node is not a ghost node...
      elem1 = -trans
      elem2 =  trans
      call MatSetValuesLocal(A,1,m2-1,1,m1-1,elem1,ADD_VALUES,ierr)
      call MatSetValuesLocal(A,1,m2-1,1,m2-1,elem2,ADD_VALUES,ierr)
    endif
  enddo
  
  !---------------------------------------------------------------------------
  ! Flux terms for boundary nodes.
  !---------------------------------------------------------------------------

  do nc=1, grid%nconnbc

    m = grid%mblkbc(nc)  ! Note that here, m is NOT ghosted.
    ng = grid%nL2G(m)

    ibc = grid%ibconn(nc)

    if(grid%ibndtyp(ibc) == 1) then ! Dirichlet BC for T
      
      i1 = ithrm_loc_p(ng)
      trans = grid%ckwet(i1) * grid%areabc(nc) / grid%distbc(nc)
      elem1 = trans
      call MatSetValuesLocal(A,1,ng-1,1,ng-1,elem1,ADD_VALUES,ierr)
      
    else if(grid%ibndtyp(ibc) == 2) then ! constant heat flux q
    
    endif
  enddo

  call VecRestoreArrayF90(grid%porosity_loc, porosity_loc_p, ierr)
  call VecRestoreArrayF90(grid%volume, volume_p, ierr)
  
  call VecRestoreArrayF90(grid%ithrm_loc, ithrm_loc_p, ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  B = A

! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  
  end subroutine CondJacobian

  end module COND_module

#undef TTEMP_LOC
#undef TTEMP
#undef TEMP
