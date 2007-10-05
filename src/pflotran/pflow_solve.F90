 module pflow_solv_module

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
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"

 public pflow_solve, pflow_kspsolver_init

 contains

 subroutine pflow_kspsolver_init(grid)

 use pflow_gridtype_module

 implicit none
 type(pflowGrid) :: grid 
 integer ierr
 
 call KSPCreate(PETSC_COMM_WORLD,grid%ksp,ierr)
 call KSPGetPC(grid%ksp, grid%pc, ierr)

! pc_type = PCILU
  if (grid%iblkfmt == 0) then
    grid%pc_type = PCJACOBI
  else
  !  grid%pc_type = PCBJACOBI
    grid%pc_type = PCILU
  endif
! grid%pc_type = PCASM
! grid%pc_type = PCNONE
  call PCSetType(grid%pc,grid%pc_type,ierr)

! call PetscOptionsSetValue('-pc_ilu_damping','1.d-10',ierr)
! call PCILUSetDamping(pc,1.d-14,ierr)

!-------krylov subspace method ----------------
  call KSPSetFromOptions(grid%ksp,ierr)
  call KSPSetInitialGuessNonzero(grid%ksp,PETSC_TRUE,ierr)

! ksp_type = KSPGMRES
!  grid%ksp_type = KSPFGMRES
  grid%ksp_type = KSPFGMRES
  call KSPSetType(grid%ksp,grid%ksp_type,ierr)
  
  call KSPSetTolerances(grid%ksp,grid%rtol,grid%atol,grid%dtol, &
      grid%maxit,ierr)



 
 end subroutine pflow_kspsolver_init
 
 
 subroutine pflow_solve(grid, newton,isucc,ierr)
 
 use pflow_gridtype_module
 use translator_mph_module
 use translator_owg_module
 use translator_vad_module
 use translator_Richards_module
 use MPHASE_module
 use Flash_module
 use VADOSE_module
 use OWG_module
 use Richards_module
 
 implicit none
 
#include "include/finclude/petscerror.h"

 type(pflowGrid) :: grid 
 KSPConvergedReason :: ksp_reason
 integer :: newton,isucc,ierr,ichange
 integer :: its_line
!integer icut
 MatStructure flag
 real*8 :: rnorm, epstol
 
 newton=0
 
 do
 
   if(grid%use_mph==PETSc_TRUE)then
 !    call Translator_MPhase_Switching(grid%xx,grid,1,ichange)
     call MPHASEResidual(grid%snes,grid%xx,grid%r,grid,ierr)
   endif
   if(grid%use_vadose==PETSc_TRUE)then
  !   call Translator_vadose_Switching(grid%xx,grid,0,ichange)
     call MPHASEResidual(grid%snes,grid%xx,grid%r,grid,ierr)
   endif

   if(grid%use_richards==PETSc_TRUE)then
  !   call Translator_richards_Switching(grid%xx,grid,0,ichange)
     call RichardsResidual(grid%snes,grid%xx,grid%r,grid,ierr)
   endif

   if(grid%use_flash==PETSc_TRUE) then
   !  call Translator_vadose_Switching(grid%xx,grid,0,ichange)
     call FLashResidual(grid%snes,grid%xx,grid%r,grid,ierr)
   endif

   if(grid%use_owg==PETSc_TRUE) then
     call Translator_OWG_Switching(grid%xx,grid%tref,grid,1,ichange,ierr)
     call OWGResidual(grid%snes,grid%xx,grid%r,grid,ierr)
   endif
   print *,' psolve; Get Res'
 
   if (ierr < 0 .or. ierr ==PETSC_ERR_ARG_DOMAIN) then
      if (grid%myrank==0) &
      print *,'pflowsolv: failed: out of range',its_line, newton,isucc,rnorm
      isucc=-1
      return    
    endif
  
       
           
   call VecNorm(grid%r,NORM_INFINITY,rnorm,ierr)
   
   ! note now grid%stol acts as convergence tolerance parameter 
   
    epstol = grid%inf_tol
   
        
       
    if (grid%myrank == 0) print *,'R2Norm = ', rnorm,epstol
    if (newton > 0 .and. rnorm < epstol) then
!   if (newton > 0 .and. rnorm < grid%stol) then
      exit ! convergence obtained
    endif

    if(grid%use_mph==PETSC_TRUE)then
      call MPHASEJacobian(grid%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
    elseif(grid%use_owg==PETSC_TRUE)then
      call OWGJacobian(grid%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
    elseif(grid%use_vadose==PETSC_TRUE)then
      call VadoseJacobian(grid%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
    elseif(grid%use_flash==PETSC_TRUE)then
      call FlashJacobian(grid%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
    elseif(grid%use_richards==PETSC_TRUE)then
      call RichardsJacobian(grid%snes,grid%xx,grid%J,grid%J,flag,grid,ierr)
    endif
     
    call VecScale(grid%r,-1D0,ierr)
    call KSPSetOperators(grid%ksp,grid%J,grid%J,SAME_NONZERO_PATTERN,ierr)
    
    call KSPSolve(grid%ksp, grid%r,grid%dxx,ierr)
    call KSPGetConvergedReason(grid%ksp,ksp_reason,ierr)
    call KSPGetIterationNumber(grid%ksp,its_line,ierr)


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
  !    call MPhase_Update(grid%xx,grid,1,ichange)
  enddo
  print *,'Finished ksp', newton
  end subroutine pflow_solve

end module pflow_solv_module
