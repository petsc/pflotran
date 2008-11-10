!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_update.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_update.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
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

module ptran_update_module

  public

contains

  subroutine ptran_update (da,da_1dof,da_kin)

    use ptran_global_module
    use trdynmem_module
    use trgamdh_module
   
    implicit none 

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"
#include "finclude/petscda.h"

    DA    :: da, da_1dof, da_kin
  
    integer :: ierr,m,n,nr,lp,llp

    real*8, pointer :: ccloc_p(:), temploc_p(:)
!   real*8 :: dcj,dcmax(commsize)
    real*8 :: sum1,sum2,totrate

    ierr = 0
    
!---compute time derivative
    call VecWAXPY(dc,-1.d0,c,cc,ierr)
    call VecAbs(dc,ierr)
    call VecMax(dc,imax,dcmax0,ierr)
    dcdt = dcmax0/dt*yrsec

!   call VecView(dc,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   call MPI_AllReduce(dc,dcmax0,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
!   MPI_SUM,PETSC_COMM_WORLD,ierr)
    
!---update solution after completed time step
    call VecCopy(cc,c,ierr)
    psi = ppsi

    if (iphase == 2) then
      psig = ppsig
      call VecCopy(ssat,sat,ierr)
    endif

!   do n = 1, ngmax
!     if (ghost_loc_p(n) == -1) cycle
!     psi(:,n) = ppsi(:,n)
!     if (iphase == 2) psig = ppsig
!     print *,'ptran_update: ',n,(psi(j,n),ppsi(j,n),j=1,ncomp)
!   enddo

    if (nkin > 0) then
      do n = 1, nlmax
!       ng = nL2G(n)
        sum1 = zero
        sum2 = zero
        do nr = 1, nkin
          totrate = zero
          do lp = npar1(nr), npar2(nr)
            llp = lp+(n-1)*nkin
!           totrate = totrate + wlam*rrkin_p(lp,n)+wlam1*rkin_p(lp,n)
            totrate = totrate + wlam*rrkin_p(llp)+wlam1*rkin_p(llp)
          enddo
          totrate = vb(n)*totrate

!---------sum contribution of reaction rates producing water and co2
          if (jh2o .gt. 0) then
            sum1 = sum1 + skin(jh2o,nr)*totrate
          endif
          if(jco2 .gt. 0) then
            sum2 = sum2 + skin(jco2,nr)*totrate
          endif
        enddo

!-------compute hydration/dehydration rate for use in FLOW
        rrtoth2o(n) = -sum1
        rrtotco2(n) = -sum2
   !    print *, 'tot rate::  ',n, jh2o,jco2,rrtoth2o(n),rrtotco2(n)

      enddo ! n-loop

!     rkin = rrkin
      call VecCopy(rrkin,rkin,ierr)
    endif

!---compute activity coefficients
    if (iact == 0) then
      gam = one
      gamx = one
    else if (iact == 1) then
      call DAGlobalToLocalBegin(da,cc,INSERT_VALUES,ccloc,ierr)
      call DAGlobalToLocalEnd(da,cc,INSERT_VALUES,ccloc,ierr)
      call VecGetArrayF90(ccloc,ccloc_p,ierr)
      call VecGetArrayF90(temploc,temploc_p,ierr)
      call trgamdh (ccloc_p,temploc_p)
      call VecRestoreArrayF90(ccloc,ccloc_p,ierr)
      call VecRestoreArrayF90(temploc,temploc_p,ierr)
    endif
    
!---update primary mineral surface area
    if (isurf == 1) then
      do n = 1, nlmax
        do m = 1, nkin
          nr = m + (n-1)*nkin
!         print *,'ptranupdate: ',n,m,nr,surf0_p(nr),phik0_p(nr),phik_p(nr),pwrsrf
          if (phik0_p(nr) > 0.d0) &
          surf_p(nr) = surf0_p(nr)*(phik_p(nr)/phik0_p(nr))**pwrsrf
        enddo
      enddo
    endif
  
  end subroutine ptran_update

end module ptran_update_module
