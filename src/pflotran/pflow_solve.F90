 module pflow_solv_module

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

 subroutine pflow_kspsolver_init(realization,solver)

 use Realization_module
 use Option_module
 use Solver_module

 implicit none

 type(realization_type) :: realization
 type(solver_type) :: solver
 PetscInt :: ierr
 
 type(option_type), pointer :: option
 
 option => realization%option
 
 call KSPCreate(PETSC_COMM_WORLD,solver%ksp,ierr)
 call KSPGetPC(solver%ksp, solver%pc, ierr)

! pc_type = PCILU
  if (option%iblkfmt == 0) then
    solver%pc_type = PCJACOBI
  else
  !  solver%pc_type = PCBJACOBI
    solver%pc_type = PCILU
  endif
! solver%pc_type = PCASM
! solver%pc_type = PCNONE
  call PCSetType(solver%pc,solver%pc_type,ierr)

! call PetscOptionsSetValue('-pc_ilu_damping','1.d-10',ierr)
! call PCILUSetDamping(pc,1.d-14,ierr)

!-------krylov subspace method ----------------
  call KSPSetFromOptions(solver%ksp,ierr)
  call KSPSetInitialGuessNonzero(solver%ksp,PETSC_TRUE,ierr)

! ksp_type = KSPGMRES
!  solver%ksp_type = KSPFGMRES
  solver%ksp_type = KSPFGMRES
  call KSPSetType(solver%ksp,solver%ksp_type,ierr)
  
  call KSPSetTolerances(solver%ksp,solver%rtol,solver%atol,solver%dtol, &
      solver%maxit,ierr)



 
 end subroutine pflow_kspsolver_init
 
 
 subroutine pflow_solve(realization,newton,newton_max,isucc,ierr)
 
 use Realization_module
 use Option_module
 use Field_module
 use Solver_module
 
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

 type(realization_type) :: realization
 KSPConvergedReason :: ksp_reason
 PetscInt :: newton,isucc,ichange
 PetscErrorCode :: ierr
 PetscInt :: newton_max
 PetscInt :: its_line
!PetscInt :: icut
 MatStructure flag
 PetscReal :: rnorm, epstol

 type(option_type), pointer :: option
 type(field_type), pointer :: field
 type(solver_type), pointer :: solver 

 option => realization%option
 field => realization%field
  
 newton=0
 
 do
 
   select case(option%imode)
#if 0
! needs to be implemented
   if(option%use_vadose==PETSc_TRUE)then
  !   call Translator_vadose_Switching(field%xx,grid,0,ichange)
     call VadoseResidual(solver%snes,field%xx,field%r,grid,ierr)
   endif
#endif
     case(MPH_MODE)
   !    call Translator_MPhase_Switching(field%xx,grid,1,ichange)
       call MPHASEResidual(solver%snes,field%xx,field%r,realization,ierr)
     case(RICHARDS_MODE)
    !   call Translator_richards_Switching(field%xx,grid,0,ichange)
       call RichardsResidual(solver%snes,field%xx,field%r,realization,ierr)
#if 0
! needs to be implemented
   if(option%use_flash==PETSc_TRUE) then
   !  call Translator_vadose_Switching(field%xx,grid,0,ichange)
     call FLashResidual(solver%snes,field%xx,field%r,grid,ierr)
   endif

   if(option%use_owg==PETSc_TRUE) then
     call Translator_OWG_Switching(field%xx,option%tref,grid,1,ichange,ierr)
     call OWGResidual(solver%snes,field%xx,field%r,grid,ierr)
   endif
!  print *,' psolve; Get Res'
#endif
   end select

   if (ierr < 0 .or. ierr ==PETSC_ERR_ARG_DOMAIN) then
      if (option%myrank==0) &
      print *,'pflowsolv: failed: out of range',its_line, newton,isucc,rnorm
      isucc=-1
      return    
    endif
  
       
           
   call VecNorm(field%r,NORM_INFINITY,rnorm,ierr)
   
   ! note now option%stol acts as convergence tolerance parameter 
   
    epstol = solver%inf_tol
   
        
       
    if (option%myrank == 0) print *,'R2Norm = ', rnorm,epstol
    if (newton > 0 .and. rnorm < epstol) then
!   if (newton > 0 .and. rnorm < option%stol) then
      exit ! convergence obtained
    endif
    
    select case(option%imode)
#if 0
! needs to be implemented
    if(option%use_owg==PETSC_TRUE)then
      call OWGJacobian(solver%snes,field%xx,solver%J,solver%J,flag,grid,ierr)
    elseif(option%use_vadose==PETSC_TRUE)then
      call VadoseJacobian(solver%snes,field%xx,solver%J,solver%J,flag,grid,ierr)
    elseif(option%use_flash==PETSC_TRUE)then
      call FlashJacobian(solver%snes,field%xx,solver%J,solver%J,flag,grid,ierr)
    elseif(option%use_richards==PETSC_TRUE)then
#endif    
      case (MPH_MODE)
        call MPHASEJacobian(solver%snes,field%xx,solver%J,solver%J,flag, &
                            realization,ierr)
      case (RICHARDS_MODE)
        call RichardsJacobian(solver%snes,field%xx,solver%J,solver%J, &
                              flag,realization,ierr)
    end select
     
    call VecScale(field%r,-1D0,ierr)
    call KSPSetOperators(solver%ksp,solver%J,solver%J,SAME_NONZERO_PATTERN,ierr)
    
    call KSPSolve(solver%ksp, field%r,field%dxx,ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
    call KSPGetIterationNumber(solver%ksp,its_line,ierr)


    newton = newton + 1
    isucc= ksp_reason

    if (newton > newton_max .or. ksp_reason < 0 ) then
      if (option%myrank==0) &
      print *,'pflowsolv: failed: ',its_line, newton,isucc,rnorm
      isucc=-1
      return    
    endif

       

  !---update solution after successful Newton-Raphson iteration
      call VecAXPY(field%xx,1.d0,field%dxx,ierr)
  !    call MPhase_Update(field%xx,grid,1,ichange)
  enddo
! print *,'Finished ksp', newton
  end subroutine pflow_solve

end module pflow_solv_module
