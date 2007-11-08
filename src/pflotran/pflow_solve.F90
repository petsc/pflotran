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

 subroutine pflow_kspsolver_init(solution,solver)

 use Solution_module
 use Option_module
 use Solver_module

 implicit none

 type(solution_type) :: solution
 type(solver_type) :: solver
 integer ierr
 
 type(option_type), pointer :: option
 
 option => solution%option
 
 call KSPCreate(PETSC_COMM_WORLD,option%ksp,ierr)
 call KSPGetPC(option%ksp, option%pc, ierr)

! pc_type = PCILU
  if (option%iblkfmt == 0) then
    option%pc_type = PCJACOBI
  else
  !  option%pc_type = PCBJACOBI
    option%pc_type = PCILU
  endif
! option%pc_type = PCASM
! option%pc_type = PCNONE
  call PCSetType(option%pc,option%pc_type,ierr)

! call PetscOptionsSetValue('-pc_ilu_damping','1.d-10',ierr)
! call PCILUSetDamping(pc,1.d-14,ierr)

!-------krylov subspace method ----------------
  call KSPSetFromOptions(option%ksp,ierr)
  call KSPSetInitialGuessNonzero(option%ksp,PETSC_TRUE,ierr)

! ksp_type = KSPGMRES
!  option%ksp_type = KSPFGMRES
  option%ksp_type = KSPFGMRES
  call KSPSetType(option%ksp,option%ksp_type,ierr)
  
  call KSPSetTolerances(option%ksp,solver%rtol,solver%atol,solver%dtol, &
      solver%maxit,ierr)



 
 end subroutine pflow_kspsolver_init
 
 
 subroutine pflow_solve(solution,newton,newton_max,isucc,ierr)
 
 use Solution_module
 use Option_module
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
#include "definitions.h"

 type(solution_type) :: solution
 KSPConvergedReason :: ksp_reason
 integer :: newton,isucc,ierr,ichange
 integer :: newton_max
 integer :: its_line
!integer icut
 MatStructure flag
 real*8 :: rnorm, epstol

 type(option_type), pointer :: option
 type(solver_type), pointer :: solver 

 option => solution%option
 !grid => solution%grid
  
 newton=0
 
 do
 
   select case(option%imode)
#if 0
   if(option%use_mph==PETSc_TRUE)then
 !    call Translator_MPhase_Switching(option%xx,grid,1,ichange)
     call MPHASEResidual(option%snes,option%xx,option%r,grid,ierr)
   endif
   if(option%use_vadose==PETSc_TRUE)then
  !   call Translator_vadose_Switching(option%xx,grid,0,ichange)
     call MPHASEResidual(option%snes,option%xx,option%r,grid,ierr)
   endif
#endif
     case(RICHARDS_MODE)
    !   call Translator_richards_Switching(option%xx,grid,0,ichange)
       call RichardsResidual(option%snes,option%xx,option%r,solution,ierr)
#if 0
   if(option%use_flash==PETSc_TRUE) then
   !  call Translator_vadose_Switching(option%xx,grid,0,ichange)
     call FLashResidual(option%snes,option%xx,option%r,grid,ierr)
   endif

   if(option%use_owg==PETSc_TRUE) then
     call Translator_OWG_Switching(option%xx,option%tref,grid,1,ichange,ierr)
     call OWGResidual(option%snes,option%xx,option%r,grid,ierr)
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
  
       
           
   call VecNorm(option%r,NORM_INFINITY,rnorm,ierr)
   
   ! note now option%stol acts as convergence tolerance parameter 
   
    epstol = solver%inf_tol
   
        
       
    if (option%myrank == 0) print *,'R2Norm = ', rnorm,epstol
    if (newton > 0 .and. rnorm < epstol) then
!   if (newton > 0 .and. rnorm < option%stol) then
      exit ! convergence obtained
    endif
    
    select case(option%imode)
#if 0
    if(option%use_mph==PETSC_TRUE)then
      call MPHASEJacobian(option%snes,option%xx,option%J,option%J,flag,grid,ierr)
    elseif(option%use_owg==PETSC_TRUE)then
      call OWGJacobian(option%snes,option%xx,option%J,option%J,flag,grid,ierr)
    elseif(option%use_vadose==PETSC_TRUE)then
      call VadoseJacobian(option%snes,option%xx,option%J,option%J,flag,grid,ierr)
    elseif(option%use_flash==PETSC_TRUE)then
      call FlashJacobian(option%snes,option%xx,option%J,option%J,flag,grid,ierr)
    elseif(option%use_richards==PETSC_TRUE)then
#endif    
      case (RICHARDS_MODE)
        call RichardsJacobian(option%snes,option%xx,option%J,option%J, &
                              flag,solution,ierr)
    end select
     
    call VecScale(option%r,-1D0,ierr)
    call KSPSetOperators(option%ksp,option%J,option%J,SAME_NONZERO_PATTERN,ierr)
    
    call KSPSolve(option%ksp, option%r,option%dxx,ierr)
    call KSPGetConvergedReason(option%ksp,ksp_reason,ierr)
    call KSPGetIterationNumber(option%ksp,its_line,ierr)


    newton = newton + 1
    isucc= ksp_reason

    if (newton > newton_max .or. ksp_reason < 0 ) then
      if (option%myrank==0) &
      print *,'pflowsolv: failed: ',its_line, newton,isucc,rnorm
      isucc=-1
      return    
    endif

       

  !---update solution after successful Newton-Raphson iteration
      call VecAXPY(option%xx,1.d0,option%dxx,ierr)
  !    call MPhase_Update(option%xx,grid,1,ichange)
  enddo
! print *,'Finished ksp', newton
  end subroutine pflow_solve

end module pflow_solv_module
