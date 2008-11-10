!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: trstdyst.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: trstdyst.F90,v $
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

!RTM: trstdyst() "transport steady-state" calculates the $I_m$ term of the
!RTM: governing equation.

  module trstdyst_module

  
  private

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  public trstdyst
  contains

  subroutine trstdyst

  use ptran_global_module
  use trdynmem_module

  implicit none

  integer :: ierr,lp,llp,n,nr,nnr
  
  real*8 :: dphik, sum
  real*8, pointer :: porosity_p(:)
  
  ierr = 0
  
!*********************************************************************** 
!     update mineral volume fractions
!***********************************************************************

!-----integrate mineral volume fraction equation using explicit FD
      do n = 1, nlmax
        do nr = 1, nkin
        
          nnr = nr+(n-1)*nkin

!---------compute change in mineral concentration
          dphik = zero
          do lp = npar1(nr), npar2(nr)
            llp = lp+(n-1)*nkin
!           dphik = dphik + wlam*rrkin_p(lp,n)+wlam1*rkin_p(lp,n)
            dphik = dphik + wlam*rrkin_p(llp)+wlam1*rkin_p(llp)
          enddo
          dphik = dt*dphik*vbarkin(nr)
          phik_p(nnr) = phik_p(nnr) + dphik
          if (phik_p(nnr) .gt. zero) then
            if (phik_p(nnr) .lt. 1.e-12) then
              phik_p(nnr) = zero
            endif
          else if (phik_p(nnr) .lt. zero) then
            phik_p(nnr) = zero
          endif
        enddo
      enddo

!-----compute porosity from mineral volume fractions
      do n = 1, nlmax
        sum = zero
        do nr = 1, nkin
          sum=sum+phik_p(nr+(n-1)*nkin)
        enddo
        por_p(n)=one-sum
      enddo

      if (ipor .eq. 1) then
        call VecGetArrayF90(porosity,porosity_p,ierr)
        do n = 1, nlmax
          if (porosity_p(n) .gt. tolpor) then
            porosity_p(n) = por_p(n)
          else
            porosity_p(n) = tolpor
          endif
        enddo
        call VecRestoreArrayF90(porosity,porosity_p,ierr)
      endif

      return
      end subroutine trstdyst

end module trstdyst_module
