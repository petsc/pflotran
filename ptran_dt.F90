!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_dt.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_dt.F90,v $
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

module ptran_dt_module

  public
  
contains

  subroutine ptran_dt (nstep,newton,t,dt,dtmax,tfac,iflgcut,iaccel,myrank)

    implicit none

    integer :: iaccel,nstep,newton,iflgcut,myrank
    real*8 :: t,dt,dtmax,fac,tfac(*)
    
!   write(*,*) 'ptran_dt: ',nstep,myrank,iflgcut,newton,iaccel

    if (iaccel == 0) return ! constant time step

    if (iflgcut == 1) then
      iflgcut = 0
      return
    endif

    fac = 1.d0
    if (newton <= iaccel .and. newton <= 13) then
      fac = tfac(newton)
    else
      fac = 0.5d0
    endif

    dt = fac*dt
    if (nstep > 10 .and. dt > 0.3d0*t) dt = 0.3*t

    if (dt > dtmax) dt = dtmax
    
  end subroutine ptran_dt

end module ptran_dt_module