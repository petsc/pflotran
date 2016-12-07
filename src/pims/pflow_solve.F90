! module contains time-steping


 module pflow_solv_module
#include "include/finclude/petscsnes.h"
 use petscsnes
 
 private

 public pflow_solve, pflow_kspsolver_init
 
 contains
 








 
 subroutine pflow_kspsolver_init(grid, pflowsolv)
 
 use pflow_gridtype_module
 
 implicit none
 type(pflowGrid) :: grid
 type(pflow_solver_context), intent(inout) :: pflowsolv 
 integer ierr
 
 call KSPCreate(PETSC_COMM_WORLD,pflowsolv%ksp,ierr)
 call KSPGetPC(pflowsolv%ksp, pflowsolv%pc, ierr)

! pc_type = PCILU
  if (grid%iblkfmt == 0) then
    pflowsolv%pc_type = PCJACOBI
  else
  !  grid%pc_type = PCBJACOBI
    pflowsolv%pc_type = PCILU
  endif
! grid%pc_type = PCSLES
! grid%pc_type = PCASM
! grid%pc_type = PCNONE
  call PCSetType(pflowsolv%pc,pflowsolv%pc_type,ierr)

! call PetscOptionsSetValue('-pc_ilu_damping','1.d-10',ierr)
! call PCILUSetDamping(pc,1.d-14,ierr)

!-------krylov subspace method ----------------
  call KSPSetFromOptions(pflowsolv%ksp,ierr)
  call KSPSetInitialGuessNonzero(pflowsolv%ksp,PETSC_TRUE,ierr)

! ksp_type = KSPGMRES
!  grid%ksp_type = KSPFGMRES
 pflowsolv%ksp_type = KSPFGMRES
  call KSPSetType(pflowsolv%ksp,pflowsolv%ksp_type,ierr)
  
  call KSPSetTolerances(pflowsolv%ksp,grid%rtol,grid%atol,grid%dtol, &
      grid%maxit,ierr)


 end subroutine pflow_kspsolver_init
 
 
 subroutine pflow_solve(grid, pflowsolv,newton,isucc,ierr)
 
 use pflow_gridtype_module
 use translator_ims_module
 use IMS_module
 
 implicit none
 type(pflowGrid) :: grid 
 type(pflow_solver_context) :: pflowsolv
 KSPConvergedReason :: ksp_reason
 integer  newton,isucc,ierr,ichange
 integer icut, its_line
 MatStructure flag
 real*8 rnorm
 type(pflow_localpatch_info),pointer :: locpat
 
 newton=0
 locpat => grid%patchlevel_info(1)%patches(1)%patch_ptr
 do
 
    

    call IMSResidual(pflowsolv%snes,grid%xx,grid%r,grid, ierr)
    print *,' psolve; Get Res'
    
    call VecNorm(grid%r,NORM_INFINITY,rnorm,ierr)
 
  ! note now grid%stol acts as convergence tolerance parameter 
  
  if (newton > 0 .and. rnorm < (grid%stol)) then
     exit ! convergence obtained
  endif
  
  call IMSJacobin(pflowsolv%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
  print *,' psolve; Get Joc'
  
  call VecScale(grid%r,-1D0,ierr)
  call KSPSetOperators(pflowsolv%ksp,grid%J,grid%J,SAME_NONZERO_PATTERN,ierr)
  
  call KSPSolve(pflowsolv%ksp, grid%r,grid%dxx,ierr)
  call KSPGetConvergedReason(pflowsolv%ksp,ksp_reason,ierr)
  call KSPGetIterationNumber(pflowsolv%ksp,its_line,ierr)
  

  newton = newton + 1
  isucc= ksp_reason

  if (newton > grid%newton_max .or. ksp_reason < 0 ) then
    if (grid%myrank==0) &
    print *,'pflowsolv: failed: ',its_line, newton,isucc,rnorm
    isucc=-1
    return    
  endif

     

!---update solution after successful Newton-Raphson iteration
	  call VecAXPY(grid%xx,1.d0,grid%dxx,ierr)
!	  call MPhase_Update(grid%xx,grid,1,ichange)
	enddo
  print *,'Finished ksp', newton
  end subroutine pflow_solve

end module pflow_solv_module
