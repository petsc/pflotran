!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_psi.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_psi.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:37:21  lichtner
! Removed ghost node corner check using temperature = -999 and replaced with check if m > 0, where m is local node designator.
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

module ptran_psi_module

private 
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscsys.h"
#include "include/finclude/petsclog.h"


public :: trpsi, trdpsi
  contains

      subroutine trpsi (ccloc_p,pressloc_p,temploc_p,ssat_loc_p)

      use ptran_global_module
      use trdynmem_module
      use co2eos_module

      implicit none

      integer :: i,ii,j,k,m,n,ngs, ierr
      real*8 :: alogjn(ncmx),fach2o,prod,sumc,tc
      real*8 :: ccloc_p(*),pressloc_p(*),temploc_p(*),ssat_loc_p(*)
      
      real*8 :: x1m,xmol,henry,poyn,rmix,wmix,eqk,xphico2=0.5d0,vphi
      PetscScalar, pointer :: den_co2_loc_p(:), xphi_co2_loc_p(:)
!***********************************************************************
!     compute total concentrations in liquid phase    
!***********************************************************************

      if (ncmplx.eq.0) then
        do n = 1, ngmax

!         print *,myrank,n,m,rho(n),ccloc_p(1+(n-1)*ncomp),temploc_p(n)
          
          if (ghost_loc_p(n) == -1) cycle
          do j = 1, ncomp
!           ppsi(j,n) = rho(n)*ccloc_p(j+(n-1)*ncomp)
            ppsi(j,n) = ccloc_p(j+(n-1)*ncomp)
!           if (nstep > 20) &
!           write(*,'("trpsi: ",2i3,1x,a8,1x,1p10e20.10)') myrank,n, &
!           nam(j),ccloc_p(j+(n-1)*ncomp),ppsi(j,n),rho(n)
          enddo
        enddo
        goto 100
      endif

      if (icomprs == 0) then
        do n = 1, ngmax
          if (ghost_loc_p(n) == -1) cycle
          if (isothrm .ge. 1) then
            tk = temploc_p(n)+273.15d0
            do i = 1, ncmplx
              eqhom(i) = -(coef(i,1)*log(tk) &
                         + coef(i,2) &
                         + coef(i,3)*tk &
                         + coef(i,4)/tk &
                         + coef(i,5)/(tk*tk))
            enddo
          endif

          if (molal.eq.1) then
            fach2o = one
            do j = 1, nmass
              alogjn(j) = log(ccloc_p(j+(n-1)*ncomp)*gam(j,n))
            enddo
          else
            fach2o = one/(wh2o*ccloc_p(jh2o+(n-1)*ncomp))
            do j = 1, nmass
              alogjn(j) = log(ccloc_p(j+(n-1)*ncomp)*gam(j,n)*fach2o)
            enddo
          endif

          do i = 1, ncmplx
            prod = eqhom(i)*aln10
            do j = 1, nmass 
!             if(shom(j,i).ne.zero) then 
                prod = prod + shom(j,i)*alogjn(j)
!             endif
            enddo

!-----------add correction for activity of water from Pitzer model
            if (iact.eq.6) then
              prod = prod + shom(jh2o,i)*ah2o(n)
            endif
            cx(i,n) = exp(prod)/(fach2o*gamx(i,n))
!           write(*,*) 'trpsi: ',myrank,namcx(i),n,cx(i,n), &
!           prod,fach2o,gamx(i,n)
          enddo

!---------compute psi for primary species
          do j = 1, ncomp 
            sumc = ccloc_p(j+(n-1)*ncomp)
            do i = 1, ncmplx
              sumc = sumc + shom(j,i)*cx(i,n)
            enddo
            ppsi(j,n) = sumc*rho(n)
!           write(*,*) 'trpsi: ',nam(j),n,sumc,ccloc_p(j+(n-1)*ncomp), &
!           ppsi(j,n),rho(n)
          enddo
        enddo

      else ! compressed form

        do n = 1, ngmax
          if (ghost_loc_p(n) == -1) cycle
          if (isothrm .ge. 1) then
            tk = temploc_p(n)+273.15d0
            do i = 1, ncmplx
              eqhom(i) = -(coef(i,1)*log(tk) &
                       + coef(i,2) &
                       + coef(i,3)*tk &
                       + coef(i,4)/tk &
                       + coef(i,5)/(tk*tk))
            enddo
          endif

          do j = 1, nmass
            alogjn(j) = log(ccloc_p(j+(n-1)*ncomp)*gam(j,n))
          enddo

          do i = 1, ncmplx
            prod = eqhom(i)*aln10    !-log(gamx(i,n))
            do k = 1, ncmpr(i)
              j = jcmpr(k,i)
              prod = prod + cshom(k,i)*alogjn(j)
            enddo

!-----------add correction for activity of water from Pitzer model
            if (iact.eq.6) then
              prod = prod + shom(jh2o,i)*ah2o(n)
            endif

            cx(i,n) = exp(prod)/(gamx(i,n))
          enddo

!---------compute psi for primary species
          do j = 1, ncomp 
            sumc = ccloc_p(j+(n-1)*ncomp)
            do k = 1, lc(j)
              sumc = sumc + fshom(j,k)*cx(ki(j,k),n)
            enddo
            ppsi(j,n) = sumc*rho(n)
!           write(*,*) 'trpsi: ',nam(j),n,sumc,ccloc_p(j+(n-1)*ncomp), &
!           ppsi(j,n),rho(n)
          enddo
        enddo
      
      endif
      
  100 continue

      if(idebug.ge.3 .and. myrank==0) then
        do n = ibg1,ibg2
          do j = 1, ncomp
            write (iunit2,510) n,nam(j),ccloc_p(j+(n-1)*ncomp),ppsi(j,n)
          enddo
          write (iunit2,511) (namcx(i),cx(i,n),i=1,ncmplx)
        enddo
  510   format('psi: n cc ppsi ',i5,1x,a8,1p2e12.4)
  511   format('psi: cmplx ',(a12,1pe12.4))
      endif

!==============================================================
!     gas phase: two-phase (liquid + gas) or all gas
!==============================================================

      if (iphase == 1) return

!-----compute total concentrations in gas phase
     call VecGetArrayF90(xphi_co2_loc, xphi_co2_loc_p, ierr)
	  call VecGetArrayF90(den_co2_loc, den_co2_loc_p, ierr)
	   do n = 1, ngmax
        if (ghost_loc_p(n) == -1) cycle
        tc = temploc_p(n)
        tk = tc + 273.15d0

        m = nG2L(n)
       ! if(m<=0) print *, 'ptran_psi:',n,m
        if(ssat_loc_p(n).ge.one) then         ! all liquid block
          do j = 1, ncomp 
            ppsig(j,n) = zero
            do i = 1, ngas
              pgas(i,n) = zero   
            enddo
          enddo

        else

          if (isothrm .ge. 1) then
            do ngs = 1, ngas 
              ii = ncmplx + ncxkin + ngs
              eqgas(ngs) = -(coef(ii,1)*log(tk) &
                           + coef(ii,2) &
                           + coef(ii,3)*tk &
                           + coef(ii,4)/tk &
                           + coef(ii,5)/(tk*tk))
            enddo
          endif
          
          if (jco2g > 0 .and. iphase == 2) then
            xphico2 = xphi_co2_loc_p(n)
            
            if (xphico2 .le. 0.d0) xphico2 = 1.d0
            
         !   call henry_co2_noderiv(xmol,x1m,tc,pgas(jco2g,n)*1.d5*rgasjj*tk, &
          !  xphico2,henry,poyn)
            
         !    call henry_co2_noderiv(xmol,x1m,tc,pressloc_p(n), &
         !    xphico2,henry,poyn)

              call Henry_duan_sun(pressloc_p(n) *1D-5, tc,  henry)  
            
            !note: henry coef. = H/xphico2
            
!          print *,'ptran_psi: ',n,pressloc_p(n),tc,xmol,x1m,xphico2,henry
            
            ! (Vphi [cm^3/mol]: Garcia, 2001)
           ! vphi= 37.51d0 + (-9.585d-2 + (8.74d-4 - 5.044d-7*tc)*tc)*tc
            
            ! fmwh2o, fmwco2: g/mol; rmix, rho: g/cm^3; wmix: g/mol
            ! ccloc_p: mol/L
            
           ! rmix = rho(n) + fmwh2o*ccloc_p(jco2+(n-1)*ncomp)*1.d-3*(fmwco2-rho(n)*vphi)
           ! wmix = rmix + (fmwh2o-fmwco2)*ccloc_p(jco2+(n-1)*ncomp)*1.d-3

!           c = rho x, p = x H/phi, p = K c, K = H/(rho phi)
          !  eqk = henry*1.e-5/(ccloc_p(jco2+(n-1)*ncomp)+ccloc_p(jh2o+(n-1)*ncomp))
!           eqk = henry*1.e-5/56d0 ! this form is only approximate-pcl

!           eqk = henry*1.e-5*fmwh2o*1.d-3/xphico2/(1.d0+ &
!           ccloc_p(jco2+(n-1)*ncomp)*1.d-3*fmwh2o)
!           eqk = henry*1.e-8*fmwh2o/xphico2

!            eqgas(jco2g) = log10(eqk)
            eqgas(jco2g)=  - log(henry*xphico2) / aln10 
          

!           print *,'ptran_psi: ',n,ngs,tc,henry,ssat_loc_p(n), &
!           ccloc_p(jco2+(n-1)*ncomp),ccloc_p(jh2o+(n-1)*ncomp),log10(eqk),eqgas(jco2g) !, &
!           vphi,rmix,wmix,rgasjj*tk,pgas(jco2g,n)
      
		     endif

          fach2o = one
          if (molal.eq.0) fach2o = one/(wh2o*ccloc_p(jh2o+(n-1)*ncomp))
          do j = 1, nmass
            alogjn(j) = log(ccloc_p(j+(n-1)*ncomp)*gam(j,n)*fach2o)
          enddo

!---------pgas_l = K_l prod_j a_j^(nu_jl)/(10 RT) [moles/L]
          do ngs = 1, ngas 
            prod = eqgas(ngs)*aln10
            do j = 1, nmass
              prod = prod+alogjn(j)*sgas(j,ngs)

!             print *,'pgas= ',j,prod,nmass,alogjn(j),sgas(j,ngs),eqgas(ngs),aln10
            enddo
!-----------add correction for activity of water from Pitzer model
            if (iact.eq.6) then
              prod = prod+sgas(jh2o,ngs)*ah2o(n)
            endif
            pgas(ngs,n) = exp(prod)/(rgasjj*tk) ! [moles/L]
            
            pgas(ngs,n)= pgas(ngs,n)*den_co2_loc_p(n)
                
!          print *,'ptranpsi: ',n,pgas(ngs,n),exp(prod),eqgas(ngs), &
!               den_co2(m)
          enddo

!---------compute psig
          do j = 1, ncomp
            ppsig(j,n) = zero
            do i = 1, ngas
              ppsig(j,n) = ppsig(j,n)+sgas(j,i)*pgas(i,n)
              
!             print *,'ptranpsi: ',n,i,j, sgas(j,i),pgas(i,n),ppsig(j,n),wmix
            enddo
          enddo
        endif
      enddo

      call VecRestoreArrayF90(xphi_co2_loc, xphi_co2_loc_p, ierr)
	  call VecRestoreArrayF90(den_co2_loc, den_co2_loc_p, ierr)

      if(idebug.gt.1 .and. myrank==0) then
        do n = 1, ngmax
          if (ghost_loc_p(n) == -1) cycle
          write (iunit2,520) n,(ppsig(j,n),j=1,ncomp)
        enddo
  520   format(' n ppsig ',i5,10e10.3/(13x,10e10.3))
      endif

      return
      end subroutine trpsi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!  PURPOSE:
   
!  This routine computes the derivatives of total concentration in 
!  aqueous and gaseous phases with respect to primary species.
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine trdpsi (ccloc_p,temploc_p,ssat_loc_p)

      use ptran_global_module
      use trdynmem_module

      implicit none

      integer :: i,j,k,l,n

      real*8 :: sum,cln,ucln,cjn,ucjn

      real*8 :: ccloc_p(*),temploc_p(*),ssat_loc_p(*)

!-----------------------------------------------------------------------
!  compute partial derivatives of psij for homogeneous equilibria-liquid
!  dpsi(j,l,n) = c_l dpsi_j/dc_l, that is the second subscript has
!  derivatives with respect to the l_th component (l=1,ncomp) for the 
!  j_th species. compute only the upper half (incl. the diagonal) and 
!  lower half is directly set as dpsi(l,j) = dpsi(j,l)
!-----------------------------------------------------------------------

      if (ncmplx.eq.0) then
        if (loglin.eq.0) then
          do n  = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do l = 1, ncomp
              dpsi(l,l,n) = rho(n)*ccloc_p(l+(n-1)*ncomp)
              do j = l+1, ncomp
                dpsi(j,l,n) = zero
                dpsi(l,j,n) = zero
              enddo
            enddo
          enddo
        else
          do n  = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do l = 1, ncomp
              dpsi(l,l,n) = rho(n)
              do j = l+1, ncomp
                dpsi(j,l,n) = zero
                dpsi(l,j,n) = zero
              enddo
            enddo
          enddo
        endif
        goto 100
      endif

      if (loglin == 0) then ! log

        if (icomprs == 0) then

          do n  = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do l = 1, ncomp
                sum = ccloc_p(l+(n-1)*ncomp)
                do i = 1, ncmplx
                  sum = sum + shom(l,i)*shom(l,i)*cx(i,n)
                enddo
                dpsi(l,l,n) = sum*rho(n)
                do j = l+1, ncomp
                  sum = zero
                  do i = 1, ncmplx
                    sum = sum + shom(j,i)*shom(l,i)*cx(i,n)
                  enddo
                  dpsi(j,l,n) = sum*rho(n)
                  dpsi(l,j,n) = dpsi(j,l,n)
                enddo
            enddo
          enddo

        else ! compressed form
        
          do n = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do j = 1, ncomp
              sum = ccloc_p(j+(n-1)*ncomp)
              do k = 1, lc(j)
                sum = sum + fshom(j,k)*fshom(j,k)*cx(ki(j,k),n)
              enddo
              dpsi(j,j,n) = sum*rho(n)
              do l = j+1, ncomp
                sum = zero
                do k = 1, llc(j,l)
                  sum = sum + ffshom(j,l,k)*cx(kki(j,l,k),n)
                enddo
                dpsi(j,l,n) = sum*rho(n)
                dpsi(l,j,n) = dpsi(j,l,n)
              enddo
            enddo
          enddo
        
        endif
        
      else ! linear

        if (icomprs == 0) then
          do n  = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do l = 1, ncomp
              cln = ccloc_p(l+(n-1)*ncomp)
              ucln = one/cln
              do j = 1, ncomp
                sum = zero
                do i = 1, ncmplx
                  sum = sum + shom(j,i)*shom(l,i)*cx(i,n)
                enddo
                sum = sum*ucln
                if (j == l) sum = one + sum
                dpsi(j,l,n) = sum*rho(n)
              enddo
            enddo
          enddo
          
        else ! compressed form
        
          do n = 1, ngmax
            if (ghost_loc_p(n) == -1) cycle
            do j = 1, ncomp
              cjn = ccloc_p(j+(n-1)*ncomp)
              ucjn = one/cjn
              sum = cjn
              do k = 1, lc(j)
                sum = sum + fshom(j,k)*fshom(j,k)*cx(ki(j,k),n)
              enddo
              dpsi(j,j,n) = sum*rho(n)*ucjn
              do l = j+1, ncomp
                cln = ccloc_p(l+(n-1)*ncomp)
                ucln = one/cln
                sum = zero
                do k = 1, llc(j,l)
                  sum = sum + ffshom(j,l,k)*cx(kki(j,l,k),n)
                enddo
                dpsi(j,l,n) = sum*rho(n)*ucln
                dpsi(l,j,n) = dpsi(j,l,n)
              enddo
            enddo
          enddo
          
        endif
      endif

  100 continue

!-----------------------------------------------------------------------
!     compute partial derivatives of psij for homogeneous equilibria-gas
!-----------------------------------------------------------------------
      
      if (iphase.eq.1) goto 200

      if (ngas.eq.0) then
        do n = 1, ngmax
          if (ghost_loc_p(n) == -1) cycle
          do l = 1, ncomp
            do j = l, ncomp
              dpsig(j,l,n) = zero
              dpsig(l,j,n) = zero
            enddo
          enddo
        enddo
        goto 200
      endif

    if (loglin.eq.0) then

      do n  = 1, ngmax
        if (ghost_loc_p(n) == -1) cycle
        if(ssat_loc_p(n).ge.one) then
          do l = 1, ncomp
            do j = l, ncomp
              dpsig(j,l,n) = zero
              dpsig(l,j,n) = zero
            enddo
          enddo
        else
          do l = 1, ncomp
            sum = zero
            do i = 1, ngas
              sum = sum + sgas(l,i)*sgas(l,i)*pgas(i,n)
              
!             print *,'ptranpsi: ',n,i,l,pgas(i,n),sgas(l,i),sum
            enddo
            dpsig(l,l,n) = sum
              
!           print *,'ptranpsi: ',n,l,l,dpsig(l,l,n)
            
            do j = l+1, ncomp
              sum = zero
              do i = 1, ngas
                sum = sum + sgas(j,i)*sgas(l,i)*pgas(i,n)
              enddo
              dpsig(j,l,n) = sum
              dpsig(l,j,n) = dpsig(j,l,n)
              
!             print *,'ptranpsi: ',n,j,l,dpsig(j,l,n)
            enddo
          enddo
        endif
      enddo

    else ! linear case

      do n  = 1, ngmax
        if (ghost_loc_p(n) == -1) cycle
        if(ssat_loc_p(n).ge.one) then
          do l = 1, ncomp
            do j = l, ncomp
              dpsig(j,l,n) = zero
              dpsig(l,j,n) = zero
            enddo
          enddo
        else
          do l = 1, ncomp
            ucln = one/ccloc_p(l+(n-1)*ncomp)
            do j = 1, ncomp
              sum = zero
              do i = 1, ngas
                sum = sum + sgas(j,i)*sgas(l,i)*pgas(i,n)
              enddo
              dpsig(j,l,n) = sum*ucln
            enddo
          enddo
        endif
      enddo

    endif

  200 continue

      if(idebug.ge.2 .and. myrank==0) then
!     if(nstep > 53) then
        write(iunit2,*) 'proc =',myrank
        write (iunit2,500) t
        ibg1 = 1
        ibg2 = nlmax
        do n = ibg1, ibg2
          do l = 1, ncomp
            write (iunit2,520) n,l,(dpsi(j,l,n),j = 1,ncomp)
          enddo
        enddo
        if (iphase.eq.2 .or. iphase.eq.0) then
          write (iunit2,510) t
          do n = ibg1, ibg2
            do l = 1, ncomp
              write (iunit2,520) n,l,(dpsig(j,l,n),j = 1,ncomp)
            enddo
          enddo
        endif
      endif
 500  format(3x,' n   l   dpsi(j,l,n)  for liquid phase, t =',1pe11.4/)
 510  format(3x,' n   l   dpsig(j,l,n) for gas phase, t =',1pe11.4/)
 520  format(i5,i4,1p7e10.3,/(1p7e10.3))

      return
      end subroutine trdpsi

end module ptran_psi_module
