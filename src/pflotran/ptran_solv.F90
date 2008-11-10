!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_solv.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_solv.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:52:26  lichtner
! No changes.
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  module ptran_solv_module

  public

  contains

  subroutine ptran_solv (its,da,da_mat,da_1dof,da_kin,ksp)

  use ptran_global_module
  use trdynmem_module
  use ptran_psi_module
  use ptran_multi_module
  use trstdyst_module
  use ptran_destroy_module

  implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
!#include "finclude/petscviewer.h"

  DA    :: da, da_mat, da_1dof, da_kin

  KSPConvergedReason :: ksp_reason
  KSP   ::  ksp

  integer :: iconv, iconv0, ierr, its, j, jn, n !, newton_mx
  real*8 :: delc, res
  real*8, pointer :: cc_p(:), ccloc_p(:), x_p(:), porloc_p(:), &
                     temploc_p(:), ssat_loc_p(:), sat_loc_p(:),pressloc_p(:)
! -----------------------------------------
  ierr = 0

  call DAGlobalToLocalBegin(da_1dof,temp,INSERT_VALUES,temploc,ierr)
  call DAGlobalToLocalEnd(da_1dof,temp,INSERT_VALUES,temploc,ierr)
  call VecGetArrayF90(temploc,temploc_p,ierr)
  
  call DAGlobalToLocalBegin(da_1dof,press,INSERT_VALUES,pressloc,ierr)
  call DAGlobalToLocalEnd(da_1dof,press,INSERT_VALUES,pressloc,ierr)
  call VecGetArrayF90(pressloc,pressloc_p,ierr)

  call DAGlobalToLocalBegin(da_1dof,ssat,INSERT_VALUES,ssat_loc,ierr)
  call DAGlobalToLocalEnd(da_1dof,ssat,INSERT_VALUES,ssat_loc,ierr)
  call DAGlobalToLocalBegin(da_1dof,sat,INSERT_VALUES,sat_loc,ierr)
  call DAGlobalToLocalEnd(da_1dof,sat,INSERT_VALUES,sat_loc,ierr)
  call VecGetArrayF90(ssat_loc,ssat_loc_p,ierr)
  call VecGetArrayF90(sat_loc,sat_loc_p,ierr)

!--begin newton-raphson loop
  icut = 0
  newton = 0
  do

    call DAGlobalToLocalBegin(da,cc,INSERT_VALUES,ccloc,ierr)
    call MatZeroEntries(A,ierr)
    call DAGlobalToLocalEnd(da,cc,INSERT_VALUES,ccloc,ierr)

    call DAGlobalToLocalBegin(da_1dof,porosity,INSERT_VALUES,porloc,ierr)
    call VecSet(b,zero,ierr)
    call VecSet(x,zero,ierr)
    call DAGlobalToLocalEnd(da_1dof,porosity,INSERT_VALUES,porloc,ierr)
    
    call DAGlobalToLocalBegin(da_1dof,tortuosity,INSERT_VALUES,tort_loc,ierr)
    call DAGlobalToLocalEnd(da_1dof,tortuosity,INSERT_VALUES,tort_loc,ierr)

    call VecGetArrayF90(ccloc,ccloc_p,ierr)
    call VecGetArrayF90(porloc,porloc_p,ierr)
    
!   call VecView(ccloc,PETSC_VIEWER_STDOUT_SELF,ierr)
!   call VecView(cc,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   call VecView(ssat,PETSC_VIEWER_STDOUT_WORLD,ierr)
    
    call trpsi (ccloc_p,pressloc_p,temploc_p,ssat_loc_p)
    call trdpsi (ccloc_p,temploc_p,ssat_loc_p)

    call ptran_multi (its,ccloc_p,temploc_p,porloc_p,sat_loc_p,ssat_loc_p)

    call VecRestoreArrayF90(ccloc,ccloc_p,ierr)
    call VecRestoreArrayF90(porloc,porloc_p,ierr)

!---check convergence
    call VecGetArrayF90(b,b_p,ierr)
    iconv = 2
    r2norm = 0.d0
    do n = 1, nlmax
      do j = 1, ncomp
        jn = j+(n-1)*nmat
        if(vb(n) > 1.d0) then
          res = b_p(jn)/vb(n)
        else
          res = b_p(jn)
        endif
        if (abs(res) > eps) iconv = -2
        r2norm = max(r2norm,abs(res))
        
!       print *,'ptran_solv: ',j,n,res,r2norm
      enddo
    enddo
    call VecRestoreArrayF90(b,b_p,ierr)
  
    if(commsize >1)then
      call MPI_ALLREDUCE(iconv,iconv0, 1, MPI_INTEGER,MPI_MIN,PETSC_COMM_WORLD,ierr)
      iconv=iconv0
      call MPI_ALLREDUCE(r2norm,res, 1, MPI_DOUBLE_PRECISION,MPI_MAX,PETSC_COMM_WORLD,ierr)
      r2norm=res
    endif 

    
    if (newton > 0 .and. iconv == 2) exit

    call KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN,ierr)
    
!   call KSPSetRHS(ksp, b, ierr)
!   call KSPSetSolution(ksp, x, ierr)
!   call KSPSolve(ksp,ierr)
    call KSPSolve(ksp,b,x,ierr)
    
    call KSPGetIterationNumber(ksp,its,ierr)
    
!   call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)

!   call KSPGetResidualHistory(ksp,num_history,ierr)
    call KSPGetConvergedReason(ksp,ksp_reason,ierr)

!typedef enum {/* converged */
!             KSP_CONVERGED_RTOL               =  2,
!             KSP_CONVERGED_ATOL               =  3,
!             KSP_CONVERGED_ITS                =  4,
!             KSP_CONVERGED_QCG_NEG_CURVE      =  5,
!             KSP_CONVERGED_QCG_CONSTRAINED    =  6,
!             KSP_CONVERGED_STEP_LENGTH        =  7,
!             /* diverged */
!             KSP_DIVERGED_ITS                 = -3,
!             KSP_DIVERGED_DTOL                = -4,
!             KSP_DIVERGED_BREAKDOWN           = -5,
!             KSP_DIVERGED_BREAKDOWN_BICG      = -6,
!             KSP_DIVERGED_NONSYMMETRIC        = -7,
!             KSP_DIVERGED_INDEFINITE_PC       = -8,
!             KSP_CONVERGED_ITERATING          =  0} KSPConvergedReason;

!   if (myrank==0) &
!   write(*,*) 'ptran_solv-KSP: reason=',ksp_reason,' its=',its, &
!   newton_max,' proc= ',myrank,' step= ',nstep

    newton = newton + 1

!   call MPI_AllReduce(newton,newton_mx,1,MPI_INTEGER, &
!   MPI_MAX,PETSC_COMM_WORLD,ierr)
!   write(*,*) 'ptransolv: ',myrank,nstep,newton,newton_max,newton_mx
!   if (newton < newton_mx) newton = newton_mx

!---reduce time step
 911 continue
    if (newton > newton_max .or. ksp_reason < 0 .or. iflgcut == 1) then
      icut = icut + 1
      
      if (icut > icut_max) then
        if (myrank==0) &
        write(*,*) 'maximum number of cuts exceeded: &
      & stop: ',' nstep= ',nstep,' newt= ',newton,' cuts= ',icut, &
      ' ksp_reason= ',ksp_reason
        call ptran_destroy (da,da_mat,da_1dof,da_kin)
        stop
      endif

      if (myrank==0) &
      write(*,'(''ptransolv-cut time step: '',i6,'' t= '',1pe12.4,'' dt= '',1pe12.4, &
      &'' icut= '',i2,''/'',i2,'' ksp: '',i3,'' newt= '',i3,''/'',i3, &
      &'' iflgcut= '',i2,'' inf_norm= '',1pe12.4)') &
      nstep,t*uyrsec,dt*uyrsec,icut,icut_max,ksp_reason,newton,newton_max,iflgcut,r2norm
!      print *,'ptransolv: cut time step: ',nstep,icut,'/',icut_max, &
!     t*uyrsec,dt*uyrsec,' ksp: ',ksp_reason,' nwtn: ',newton,iflgcut, &
!     ' r2norm= ',r2norm
      
      t = t - dt
      dt = 0.5d0*dt
      t = t + dt
      call VecCopy(c,cc,ierr)
      newton = 0
!     iflgcut = 1
!     call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
      cycle
    endif

    itpetsc(newton) = its

!---update solution after successful Newton-Raphson iteration
    if (loglin == 1) then
      call VecAXPY(cc,1.d0,x,ierr)
    else
      call VecGetArrayF90(cc,cc_p,ierr)
      call VecGetArrayF90(x,x_p,ierr)
      do n = 1, nlmax*ncomp
        delc = x_p(n)
!       write(*,*) 'ptransolv: n/proc/dc/cc=',myrank,nstep,n,x_p(n),cc_p(n)
        if (abs(delc) > abs(tolexp)) then
!         write(*,*) 'ptransolv: n/proc/dc/cc=',nstep,myrank,x_p(n),cc_p(n)
          if (tolexp > 0.d0) then
            x_p(n) = 5.d0*delc/abs(delc)
          else
            iflgcut = 1
!           call VecRestoreArrayF90(cc,cc_p,ierr)
!           call VecRestoreArrayF90(x,x_p,ierr)
!           goto 911
            t = t - dt
            dt = 0.5d0*dt
            t = t + dt
            call VecCopy(c,cc,ierr)
            newton = 0
!           iflgcut = 1
            print *,'ptransolv - cut time step: ',nstep,icut, &
            t*uyrsec,dt*uyrsec,tolexp,delc
!           cycle
            exit
          endif
        endif
      enddo

      cc_p = cc_p*exp(x_p)
      
      call VecRestoreArrayF90(cc,cc_p,ierr)
      call VecRestoreArrayF90(x,x_p,ierr)

!     call MPI_AllReduce(iflgcut,integer_out,1,MPI_INTEGER, &
!     MPI_SUM,PETSC_COMM_WORLD,ierr)
!     write(*,*) 'ptransolv: ',myrank,iflgcut,integer_out
!     if (integer_out > 0) iflgcut = 1
    endif

  enddo

  call VecRestoreArrayF90(temploc,temploc_p,ierr)
  call VecRestoreArrayF90(sat_loc,sat_loc_p,ierr)
  call VecRestoreArrayF90(ssat_loc,ssat_loc_p,ierr)
  call VecRestoreArrayF90(pressloc,pressloc_p,ierr)
!---integrate mineral mass transfer equations
    if (nkin > 0) then
!!    call VecGetArrayF90(phik,phik_p,ierr)
      call trstdyst
!!    call VecRestoreArrayF90(phik,phik_p,ierr)
    endif

  end subroutine ptran_solv

end module ptran_solv_module
