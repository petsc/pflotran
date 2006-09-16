!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: trgamdh.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: trgamdh.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.2  2004/04/06 17:54:07  lichtner
! Revised check for corner ghost nodes using local index.
!
! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
! initial entry
!
! Revision 1.2  2003/05/13 14:56:54  lichtner
! added header
!

!  pFLOTRAN Version 1.0 LANL
!-----------------------------------------------------------------------
!  Date             Author(s)                Comments/Modifications
!-----------------------------------------------------------------------
!  May  2003        Peter C. Lichtner        Initial Implementation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module trgamdh_module

  public

contains

      subroutine trgamdh (ccloc_p,temploc_p)

      use ptran_global_module
      use trdynmem_module

      implicit none

      integer:: i,it,j,l,n,nnode

      real*8 :: ccloc_p(*),temploc_p(*)

      real*8 :: dcdi,den,dgamdi,didi,f,fach2o,fcomp, &
      si,sum,sumcx,sqrti,ssi, &
      prod,ajn,alogac(ncmx),alogacx(ncxmx)

      do n = 1, ngmax

!-------check for corner ghost nodes
!       if (tk == -999.d0) cycle
        if (ghost_loc_p(n) == -1) cycle
            
        tk = temploc_p(n)

!       if (tk == -999.d0) &
!       write(*,*) 'trgamdh: ',myrank,nstep,n,tk,ngmax,nlmax,nmax, &
!       (ccloc_p(j+(n-1)*ncomp),j=1,ncomp)

        fach2o = one
        if (molal.eq.0) fach2o = one/(wh2o*ccloc_p(jh2o+(n-1)*ncomp))

        fcomp = zero
        do j = 1, ncomp 
          fcomp = fcomp+z(j)*z(j)*ccloc_p(j+(n-1)*ncomp)
        enddo

!**************************************************************** 
!       begin loop in activity coefficient iteration
!**************************************************************** 

        it = 0
!       si = half*fcomp

!       write(*,*) 'trgamdh: ',myrank,nstep,n,sionic(n),tk,ngmax
        
        si = sionic(n)
        ssi = si
        
   20   continue

        it = it+1 
        if (it.gt.50) then
          nnode = n
          goto 130
        endif

!       if (it.gt.175) then
!         write(*,*) 'trgamedh: ',n,it,ssi,si, &
!         (nam(j),gam(j,n),ccloc_p(j+(n-1)*ncomp),j=1,ncomp)
!       endif

        sumcx = zero
        do i = 1, ncmplx
          sumcx = sumcx+zx(i)*zx(i)*cx(i,n)
        enddo
        f = half*(sumcx+fcomp)*fach2o

!**************************************************************** 
!       check convergence of ionic strength 
        if (abs(f-si).lt.1.d-6*si) goto 90

!**************************************************************** 
!       compute derivative of reduced complex concentration 
!       wrt ionic strength
!**************************************************************** 
        if (ncmplx.gt.0) then
          didi = zero
          sqrti = sqrt(si)
          do i = 1, ncmplx
            if (zx(i).ne.zero) then
              sum = half*adebye*zx(i)*zx(i) &
              /(sqrti*(one+bdebye*ax0(i)*sqrti)**2)-bextendx(i)
              do j = 1, ncomp
                if(z(j).ne.zero) then
                  dgamdi = -half*adebye*z(j)**2/(sqrti* &
                  (one+bdebye*a0(j)*sqrti)**2)+bextend(j)
                  sum=sum+shom(j,i)*dgamdi
                endif
              enddo
              dcdi = cx(i,n)*aln10*sum
              didi = didi+half*zx(i)**2*dcdi
            endif
          enddo
          den = one-didi
          if (den.ne.zero) then
            ssi = (f-si*didi)/den
          else
            ssi = f
          endif
        else
          ssi = f
        endif

!-------check for negative ionic strength
        if (ssi .lt. zero .and. myrank==0) then
          write(*,*) 'ionic strength negative: stop!'
          write(*,*) 'node= ',n,' si =',ssi,si,' f= ',f,' didi= ',didi
          write(*,*) 'fcomp =',fcomp,' sumcx =',sumcx,' sqrtsi =',sqrti
          if (myrank == 0) then
            write(iunit2,*) 'primary species'
            do j = 1, ncomp
              write(iunit2,*) nam(j),(ccloc_p(j+(l-1)*ncomp),l=1,ngmax)
            enddo
            write(iunit2,*) 'secondary species'
            do i = 1, ncmplx
              write(iunit2,*) namcx(i),(cx(i,l),l=1,ngmax)
            enddo
          endif
          stop
        endif

        si = ssi
        sqrti = sqrt(si) 
        do j = 1, ncomp 
          if (z(j).ne.zero) then
            alogac(j) = -z(j)**2*adebye*sqrti/(one+bdebye*a0(j)*sqrti) &
            +bextend(j)*si
            gam(j,n) = exp(alogac(j)*aln10) 
          endif
        enddo

        do i = 1, ncmplx
          if (zx(i).ne.zero) then
            alogacx(i) = -zx(i)**2*adebye*sqrti/ &
            (one+bdebye*ax0(i)*sqrti)+bextendx(i)*si
            gamx(i,n) = exp(alogacx(i)*aln10)
          endif
!**************************************************************** 
!         compute concentrations of complexes 
!**************************************************************** 
          cx(i,n) = zero
          if (isothrm .ge. 1) then
              eqhom(i) = -(coef(i,1)*log(tk) &
              + coef(i,2) &
              + coef(i,3)*tk &
              + coef(i,4)/tk &
              + coef(i,5)/(tk*tk))
          endif
          prod = eqhom(i)*aln10  !-log(gamx(i,n)) 
          do j = 1, nmass 
            if (shom(j,i).ne.zero) then 
              ajn = ccloc_p(j+(n-1)*ncomp)*gam(j,n)*fach2o
              if (ajn.gt.zero) then
                prod = prod+log(ajn)*shom(j,i)
              endif
            endif
          enddo
          cx(i,n) = exp(prod)/(fach2o*gamx(i,n))
        enddo

!       write(*,*) 'trgamdh: ',myrank,n,ncomp,(nam(j),gam(j,n),j=1,ncomp)
!       write(*,*) 'trgamdh-cx: ',myrank,n,ncmplx, &
!       (namcx(i),gamx(i,n),i=1,ncmplx)

        goto 20 

   90   continue
   
!       write(*,*) 'trgamdh: ',myrank,n,it,si
        
        if(it > 6) write(*,1000) myrank,it 
 1000   format(' *** excessive iterations in ionic strength loop *** &
     &  proc =',i3,', iter = ',i3)
        sionic(n) = si

      enddo

      return

  130 continue

      if (myrank == 0) then
       ! note includes corner nodes where gam is not defined!
        write(iunit2,1010) it,nnode,ssi,si,(sionic(n),n=1,ngmax)
        write(*     ,1010) it,nnode,ssi,si,(sionic(n),n=1,ngmax)
 1010   format(' ionic strength not converging after ',i4, &
        ' iterations at node: ',i6,/, ' sionic = ',1p2e12.4/ &
        (' ',10(1pg12.4)))
        do j = 1,ncomp
          write(iunit2,'("gam: ",a8,(1p10e12.4))') nam(j), &
          (gam(j,n),n=1,ngmax) ! note includes corner nodes where gam is not defined!
          write(*,'("gam: ",a8,(1p10e12.4))') nam(j), &
          (gam(j,n),n=1,ngmax)
        enddo
      endif
      
      stop
      end subroutine trgamdh

end module trgamdh_module
