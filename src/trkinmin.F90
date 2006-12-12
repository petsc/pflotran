!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: trkinmin.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: trkinmin.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.2  2004/01/10 18:32:06  lichtner
! Began work on 2 phase capability.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!
! Revision 1.2  2003/05/13 14:57:03  lichtner
! added header
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module trkinmin_module

  public

contains

  subroutine trkinmin (ccloc_p,temploc_p)

  use ptran_global_module
  use trdynmem_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  integer :: ierr,iflgerr,i,ii,ityprxn,j,jj,l,lp,llp,n,ng,nr,npri,nsec
  
! integer :: irow1(ncomp),icol1(ncomp)

  real*8 :: ajn,alnjn,cdld,cdlod,cxin,den,dprfacdl,dlnqk,dvol, &
            fac1,fac2,fach2o,fl,fk,fkk, &
            prefac,pwrqk,qk,qrb,qrf,rtelim,rgast0,rk,rrkn, &
            satindx,sig,sum1,sum2,totrate,u3,uu1,uu2

  real*8 :: ccloc_p(*), temploc_p(*)

  real*8 :: cdl(nmat,nmat),ajnlog(ncomp),drate(ncomp)

      ierr = 0

      iflgerr = 0
      
      rgast0 = one/(rgaskj*tk0)

      do n = 1, nlmax
      
        cdl = zero

        ng = nL2G(n)

        dvol = wlam*vb(n)

!-------check for zero liquid saturation
!       if ((iphase.eq.2 .or. iphase.eq.0) .and. icksat.eq.1) then
!         if (iwet(n).eq.0) then
!           do nr = 1, nkin
!             do lp = npar1(nr), npar2(nr)
!               llp = lp+(n-1)*nkin
!               rrkin(lp,n) = zero
!               rrkin_p(llp) = zero
!             enddo
!           enddo
!           goto 101
!         endif
!       endif

        if (molal.eq.1) then
          do j = 1, nmass
            ajnlog(j) = log(ccloc_p(j+(ng-1)*ncomp)*gam(j,ng))
          enddo
        else
          fach2o = one/(wh2o*ccloc_p(jh2o+(ng-1)*ncomp))
          do j = 1, nmass
            ajnlog(j) = log(ccloc_p(j+(ng-1)*ncomp)*gam(j,ng)*fach2o)
          enddo
        endif

        if (isothrm .eq. 1) then
          tk = temploc_p(ng)
          uu1 = tk/tk0
          uu2 = one/(tk*rgaskj)-rgast0
          do nr = 1, nkin
            u3 = uu1*exp(-delh(nr)*uu2)
            do lp = npar1(nr), npar2(nr)
              rkf(lp) = rkf0(lp)*u3
              
!             print *,'trkinmin: ',nr,lp,n,ng,' ',namk(nr), &
!             tk,rkf0(lp),rkf(lp),u3
            enddo
            i = ncmplx + ncxkin + ngas + ndxkin(nr)
            eqkin(nr) = -(coef(i,1)*log(tk) &
              + coef(i,2) &
              + coef(i,3)*tk &
              + coef(i,4)/tk &
              + coef(i,5)/(tk*tk))
          enddo
        endif

        do nr = 1, nkin
        
          do lp = npar1(nr), npar2(nr)
            llp = lp+(n-1)*nkin
!           rrkin_p(lp,n) = zero
            rrkin_p(llp) = zero
          enddo

!---------compute forward and backward ion activity product
          qrf = zero
          qrb = zero
          do j = 1, nmass
            ajn = ajnlog(j)*skin(j,nr)
            if(skin(j,nr) .lt. zero) then
              qrb = qrb-ajn
            else 
              qrf = qrf+ajn
            endif
          enddo

!---------add correction for activity of water based on Pitzer model
          if (iact.eq.6) then
            alnjn = skin(jh2o,nr)*ah2o(ng)
            if(skin(jh2o,nr) .lt. zero) then
              qrb = qrb-alnjn
            else 
              qrf = qrf+alnjn
            endif
          endif

          dlnqk = qrb-qrf-eqkin(nr)*aln10

          if (-dlnqk.gt.qkmax) then
            if (iwarn.ge.2 .and. myrank==0) then
              write(*,'("WARNING!--affinity too large in trkinmin: ", &
     &        a12,i4,1p5e11.3)') namk(nr),n,dlnqk,qrb,qrf,qkmax,eqkin(nr)
              write(*,'(3(2x,a8,1pe12.4))') &
              (nam(j),ccloc_p(j+(ng-1)*ncomp),j=1,ncomp)
            endif
            iflgerr = 1
!           return
          endif

!         qk(nr,n) = exp(-dlnqk)
          qk = exp(-dlnqk)
  
!---------rate sign convention chosen so that reaction rate is positive 
!         for precipitation and negative for dissolution with mineral 
!         species on left hand side of reaction.

!         satindx = qk(nr,n)
          satindx = qk

!---------check if mineral is supersaturated or present
          if (satindx.gt.one .or. phik_p(nr+(n-1)*nkin).gt.zero) then

!-----------check if mineral saturation less than threshold 
            if (satindx.gt.one .and. satindx.lt.fkin(nr)) goto 100
            
!-----------reset precipitation threshold
            if (fkin(nr).gt.one) then
              fkin(nr) = one
            endif

!-----------multiply limiting rate by surface area
            rtelim = rlim(nr)*surf_p(nr+(n-1)*nkin)

!-----------compute rate for each parallel reaction and add contribution
!           to jacobian
            do lp = npar1(nr), npar2(nr)
              llp = lp+(n-1)*nkin

!-------------compute prefactor
              npri = nkinpri(lp)
              nsec = nkinsec(lp)
              prefac = one
              if (npri.gt.0 .or. nsec.gt.0) then
                prefac = zero
                if (npri .gt. 0) then
                  do j = 1, npri
                    prefac = prefac + skinpri(j,lp)*ajnlog(jpri(j,lp))
                  enddo
                endif
                if (nsec .gt. 0) then
                  do i = 1, nsec
                    ii = isec(i,lp)
                    cxin = log(gamx(ii,n)*cx(ii,n)*fach2o)
                    prefac = prefac + skinsec(i,lp)*cxin
                  enddo
                endif
                prefac = exp(prefac)
              endif

!-----------------------------------------------------------------------
!  itypkin              reaction rate types
!-----------------------------------------------------------------------
!     20  TST general: rk prod_j a_j^vj prod_i a_i^vi (1-KQ^sigma)
!     25  TST transport limited: 
!           rk prod_j a_j^vj prod_i a_i^vi (1-KQ^sigma)/(1+f*KQ^sigma)
!-----------------------------------------------------------------------

              rk   = rkf(lp)*surf_p(nr+(n-1)*nkin)
              sig  = one/sigma(lp)
              ityprxn = itypkin(lp)
          
!             print *, 'trkinmin: ',lp,n,ng,nr,npri,nsec,qk,prefac, &
!             phik_p(nr+(n-1)*nkin),rk,rkf(lp),surf_p(nr+(n-1)*nkin)

              if (itypkini(ityprxn) == 6) goto 25

   20         continue

!-------------transition state rate law
              pwrqk = satindx**sig
              rrkn = -rk*prefac*(one-pwrqk)
              
!             write(*,*) 'trkinmin: ',n,nr,nbuff,ityprxn, &
!             npar1(nr),npar2(nr),npri,nsec,rk,prefac,pwrqk,sig, &
!             rk,eqkin(nr),satindx,rrkn,phik_p(nr+(n-1)*nkin), &
!             surf_p(nr+(n-1)*nkin),(skin(j,nr),j=1,ncomp)

!-------------jacobian
              fk = rk*prefac*dvol

              fkk = fk*sig*pwrqk
              do l = 1, ncomp
                drate(l) = fkk*skin(l,nr)
              enddo

              do l = 1, ncomp

                if (drate(l).ne.zero) then

!-----------------diagonal terms
                  cdld = skin(l,nr)*drate(l)
                  cdl(l,l) = cdl(l,l) + cdld

!-----------------off diagonal terms: make use of symmetry
                  do j = l+1, ncomp
                    cdlod = skin(j,nr)*drate(l)
                    cdl(j,l) = cdl(j,l) + cdlod
                    cdl(l,j) = cdl(l,j) + cdlod
                  enddo
                endif
              enddo

!-------------derivative of prefactor term
              if (npri.gt.0 .or. nsec.gt.0) then
                do l = 1, ncomp
                  drate(l) = zero
                enddo
                do l = 1, npri
                  jj = jpri(l,lp)
                  dprfacdl = zero
                  if (npri.gt.0) dprfacdl = skinpri(l,lp)
                  if (nsec.gt.0) then
                    do i = 1, nsec
                      ii = isec(i,lp)
                      dprfacdl = dprfacdl + shom(jj,ii)*skinsec(i,lp)
                    enddo
                  endif
                  drate(jj) = dprfacdl*rrkn*dvol
                enddo

!---------------construct jacobian
                do l = 1, ncomp
                  do j = 1, ncomp 
                    cdl(j,l) = cdl(j,l) + skin(j,nr)*drate(l)
                  enddo
                enddo
              endif

              goto 10

   25         continue

!-------------transport-limited transition state rate law
              pwrqk = satindx**sig
              fac1 = one-pwrqk
              fac2 = one+rk*pwrqk/rtelim
              rrkn = -rk*prefac*fac1/fac2

!             write(*,*) 'trkinmin: ',namk(nr),n,rtelim,rk,fac1,fac2,
!    .        rrkn

              if(rrkn.eq.zero) goto 10

!-------------jacobian
              fk = rk*prefac*dvol
              den = one+rk/rtelim*pwrqk
              do l = 1, ncomp
                fl = skin(l,nr)*(one+rk/rtelim) !check
                drate(l) = fk*sig*fl*pwrqk/(den*den)
              enddo
              do l = 1, npri
                jj = jpri(l,lp)
                dprfacdl = zero
                if (npri .gt. 0) dprfacdl = skinpri(l,lp)
                if (nsec .gt. 0) then
                  do i = 1, nsec
                    ii = isec(i,lp)
                    dprfacdl = dprfacdl + shom(jj,ii)*skinsec(i,lp)
                  enddo
                endif
                drate(jj) = drate(jj) - fk*dprfacdl*(one-pwrqk)/den
              enddo
              do l=1,ncomp
                do j=1,ncomp 
                  cdl(j,l) = cdl(j,l) + skin(j,nr)*drate(l)
                enddo
              enddo

   10         continue

!-------------store rate for lpth term
!             rrkin_p(lp,n) = rrkn
              rrkin_p(llp) = rrkn

            enddo ! lp-loop
          endif
  100     continue
        enddo     ! nr-loop
  101   continue
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,cdl, &
          ADD_VALUES,ierr)
        else if (iblkfmt == 2) then
          do j = 1, ncomp
            do l = 1, ncomp
              call MatSetValuesLocal(A,1,j+(ng-1)*nmat-1,1,l+(ng-1)*nmat-1, &
              cdl(j,l),ADD_VALUES,ierr)
            enddo
          enddo
        else if (iblkfmt == 0) then
!         do j = 1, ncomp
!           irow1(j) = j+(ng-1)*ncomp-1
!           icol1(j) = j+(ng-1)*ncomp-1
!         enddo
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,cdl, &
!         ADD_VALUES,ierr)
          do j = 1, ncomp
            do l = 1, ncomp
              if (dfill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,j+(ng-1)*nmat-1,1,l+(ng-1)*nmat-1, &
                cdl(j,l),ADD_VALUES,ierr)
              endif
            enddo
          enddo
        endif
      enddo       ! n-loop

!-----print out rates
!     write(iunit2,*) '*mcyc= ',mcyc,' iter= ',iter,nkinrxn
!     do n = nmax1, nmax2
!       do nr = 1, nkin
!         do lp = npar1(nr), npar2(nr)
!           write(iunit2,'(4i4,1x,a16,1p3e12.4,i3)') n,lp,
!    .      npar1(nr),npar2(nr),namk(nr),rkin(nr,n),rrkin(nr,n)
!         enddo
!       enddo
!     enddo

!-----compute residual
      do n = 1, nlmax
        do nr = 1, nkin
          totrate = zero
          do lp = npar1(nr), npar2(nr)
            llp = lp+(n-1)*nkin
!           totrate = totrate + wlam*rrkin_p(lp,n)+wlam1*rkin_p(lp,n)
            totrate = totrate + wlam*rrkin_p(llp)+wlam1*rkin_p(llp)
          enddo
          totrate = vb(n)*totrate
          do j = 1, ncomp
            b_p(j+(n-1)*nmat) = b_p(j+(n-1)*nmat)-skin(j,nr)*totrate
          enddo
        enddo
      enddo ! n-loop
  
      return
      end subroutine trkinmin

end module trkinmin_module