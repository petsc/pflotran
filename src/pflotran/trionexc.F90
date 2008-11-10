!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: trionexc.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: trionexc.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.6  2004/03/31 23:06:06  lichtner
! Added glam.
!
! Revision 1.5  2003/08/20 16:13:13  lichtner
! minor changes
!
! Revision 1.4  2003/04/19 21:40:52  lichtner
! removed unused code
!
! Revision 1.3  2002/09/01 02:26:55  lichtner
! Added module modflotrn.
!
! Revision 1.2  2002/08/11 17:21:31  lichtner
! Added comment to specify form of ion exchange reaction.
!
! Revision 1.1.1.1  2002/01/18 23:53:37  lichtner
! Initial Entry
!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module trionexc_module

  public

contains

  subroutine trionexc (ccloc_p,porloc_p)


  use ptran_global_module
  use trdynmem_module

  implicit none

#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

!     include 'impl.h'
!     include 'paramtrs.h'
!     include 'flowtran.h'
!     include 'trdebye.h'
!     include 'triounts.h'
!     include 'trcom.h'
!     include 'trscalar.h'

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------
!     nc                ncomp     -number of primary species
!     nexs              nexsite   -number of ion exchange sites
!     nexx              nexmax    -total number of exchange cations on
!                                  different sites (= nc*nexsolid)
!     nsfm              nsrfmin   -number of surface complexation solids
!     nsfx              nsrfmx    -total number of species
!     nsfs              nsrfsit   -number of exchange sites
!-----------------------------------------------------------------------
!-----variables
!-----------------------------------------------------------------------
!     nexsolid
!     nexmax
!     nex1(m), nex2(m)
!     mex(m)
!     jex(j)
!     namex(m), namcat(ncmx)
!     xex0(m,n),xex(m,n),xxex(m,n)
!     cec0(m),cec(m,n)
!     eqiex(m)
!     dpsi(nc,nc,n)
!-----------------------------------------------------------------------
!     colloids
!-----------------------------------------------------------------------
!     namcoll(ncmx)    (not used)
!     ncoll
!     ncollex, ncolsrf
!-----------------------------------------------------------------------

!     dimension cc(nc,*),gam(nc,*),glam(nexx,*),r(*),cdl(ndf,ndf,*),
!    .  xex(nexx,*),xxex(nexx,*),cec(nexs,*),
!    .  ccsorp(nsfx,*),csorp(nsfx,*),csorpf(nsfs,*),siteden(nsfs,*),
!    .  alogpf(nsfm,*),dpsi(ncex,ncex,*)

      real*8 :: ccloc_p(:),porloc_p(:)
      real*8 :: alogjn(ncmx),dkln(nexmx),dchi(nexmx,nexmx),sdpsi(ncmx)
      real*8 :: cdl(ncomp,ncomp)

      real*8 :: fach2o,cj0,cln,dfdk,dkj0,dr,f,fac,facec,facj0, &
                suml,sumx,uz0,uzdk, &
                xj0,xl,zl,zj0
      
      real*8 :: den,dcsc,facp,facv,facvn,flam,flamn,prod, &
                resjn1,resjn2,sum

      integer :: i,iterex,j,jj,j0,jj0,k,l,ll,m,n,ng

      if(nexsolid.eq.0) goto 100
      
! Note: b_p = -r

!********1*********2*********3*********4*********5*********6*********7**
!     compute ion exchange isotherm for the reaction:
!       1/z_j0 A_j0 + 1/z_i Xz_iA_i = 1/z_j0 Xz_j0A_j0 + 1/z_i A_i
!
!     - note: sum_j xex(j,n) = 1, all n, by construction
!     - j0 = pivot element taken as first cation in list
!     - solve for kD_j0 (dkj0)
!***********************************************************************
!-----------------------------------------------------------------------
!     pseudo code do-loop construction
!     do m = 1, nexsolid                ! minerals/colloids
!       do k = nsitex1(m), nsitex2(m)   ! sites
!         do jj = nex1(k), nex2(k)      ! sorbed species
!           j = jex(jj)                 ! cation index
!         enddo
!       enddo
!     enddo
!-----------------------------------------------------------------------

!-----begin loop over nodes
      do n = 1, nlmax

        ng = nL2G(n)

        cdl = zero

        iterex = 0

        do m = 1, nexsolid
          do k = nsitex1(m), nsitex2(m)

            jj0 = nex1(k)
            j0  = jex(jj0)
            zj0 = z(j0)
            uz0 = one/zj0

!-----------conversion factor for molarity to molality
            fach2o = one
!           if (molal.eq.0) fac = one/(wh2o*cc(jh2o,n))

!           cj0  = cc(j0,n)*fach2o
            cj0  = ccloc_p(j0+(ng-1)*ncomp)*fach2o
            dkj0 = xxex(jj0,n)/cj0

            if (dkj0 .le. zero) then
              do jj = nex1(k), nex2(k)
                xxex(jj,n) = zero
              enddo
              goto 99 ! skip node
            endif

   10       iterex = iterex+1

            if(iterex .ge. 50) goto 20 ! quit-too many iterations

            xj0   = dkj0*cj0
            uzdk  = one/(zj0*dkj0)
            facj0 = (glam(jj0,ng)*dkj0/(eqiex(jj0)*gam(j0,ng)))**uz0
            f     = xj0-one
            dfdk  = cj0

            do ll = nex1(k)+1, nex2(k)
              l = jex(ll)
              zl = z(l)
              dkln(ll) = eqiex(ll)*gam(l,ng)/glam(ll,ng)*facj0**zl
!             xl   = dkln(ll)*cc(l,n)*fach2o
              xl   = dkln(ll)*ccloc_p(l+(ng-1)*ncomp)*fach2o
              f    = f + xl
              dfdk = dfdk + zl*xl*uzdk

!             write(*,*) 'trionexc: ',iterex,n,nam(j0),nam(l),cj0, &
!             cc(l,n),gam(l,n),zl,dkln(ll),cc(jh2o,n),dfdk,f,xl,uzdk, &
!             glam(jj0,n),glam(ll,n)
!             write(*,*) 'trionexc: ',m,k,n,iterex,nam(j0),nam(l), &
!             f,tol,glam(jj0,n),glam(ll,n)
            enddo

!           write(*,*) iterex,f,dfdk,dkj0,dkj0/dkln(2)

!-----------check convergence
            if(abs(f) .lt. tol) goto 30 ! convergence obtained
            dkj0 = dkj0 - f/dfdk

!           write(*,*) iterex,f,dfdk,dkj0,tol
!1001       format('ionex: n,m,iter: ',3i3,' f,dfdk,dkj0: ',1p3e12.4)

            goto 10 ! iterate

   20       continue
            write(iunit2,1000) n,m,iterex,f,dfdk,dkj0 
            write(*,1000) n,m,iterex,f,dfdk,dkj0 
 1000       format(' maximum number of iterations in ionex: n,m,iter: ', &
     &      3i3,' f,dfdk,dkj0: ',1p3e12.4)
            stop

   30       continue ! convergence obtained
   
!           print *,'convergence: ',iterex,f,tol

            dkln(jj0) = dkj0

            facec = cec(k,n)*vb(n)/dt
!-----------add porosity factor in colloid accumulation term
            if (m.gt.nexsolid-ncollex) then
              facec = facec*porloc_p(ng)
            endif

!-----------gaines-thomas convention: xex = z_j Cbar_j/sum z_i Cbar_i
            sumx = zero
            do jj = nex1(k), nex2(k) 
              j = jex(jj)
!             xxex(jj,n) = dkln(jj)*cc(j,n)
              xxex(jj,n) = dkln(jj)*ccloc_p(j+(ng-1)*ncomp)
              sumx = sumx+z(j)*xxex(jj,n) 
            enddo

!-----------compute partial derivatives and residuals of sorption 
!           isotherm
            do jj = nex1(k), nex2(k)
              j = jex(jj)
              fac = -z(j)*xxex(jj,n)/sumx
              do ll = nex1(k), nex2(k)
                dpsi(jj,ll,ng) = fac*dkln(ll)
              enddo
              dpsi(jj,jj,ng) = dpsi(jj,jj,ng) + dkln(jj)

!-------------compute chi' for vanselow convention
              if (ionex.eq.1) xxex(jj,n) = z(j)*xxex(jj,n)/sumx

!-------------residual
              dr = facec*(xxex(jj,n)-xex(jj,n))/z(j)
!             r(j+(n-1)*nmat) = r(j+(n-1)*nmat) + dr
              b_p(j+(n-1)*nmat) = b_p(j+(n-1)*nmat) - dr

!             if (ndecay.gt.0) then
!               r(j+(n-1)*nmat) = r(j+(n-1)*nmat)+
!               b_p(j+(n-1)*nmat) = b_p(j+(n-1)*nmat)+ &
!               facec*xlamhalf(j)*xxex(jj,n)/z(j)
!             endif
            enddo

!-----------compute jacobian

            if (ionex .eq. 0) then
!-------------gaines-thomas convention: 
!             Nbar_j = z_j Cbar_j/sum z_i Cbar_i
              do jj = nex1(k), nex2(k)
                j = jex(jj)
                do ll = nex1(k), nex2(k)
                  l = jex(ll)
!                 cln = cc(l,n)/z(j)
                  cln = ccloc_p(l+(ng-1)*ncomp)/z(j)
                  cdl(j,l) = cdl(j,l)+facec*dpsi(jj,ll,ng)*cln
                  if (ndecay.gt.0) then
                    cdl(j,l) = cdl(j,l)+facec*xlamhalf(j)* &
                    dpsi(jj,ll,ng)*cln
                  endif
                enddo
              enddo

            else

!-------------vanselow convention: Nbar_j = Cbar_j/sum Cbar_i

!-------------compute dchi'/dc
              do ll = nex1(k), nex2(k)
                l = jex(ll)
                suml = zero
                do jj = nex1(k), nex2(k)
                  j = jex(jj)
                  suml = suml + z(j)*dpsi(jj,ll,ng)
                enddo
                sdpsi(l) = suml
              enddo
              do jj = nex1(k), nex2(k)
                j = jex(jj)
                do ll = nex1(k), nex2(k)
                  l = jex(ll)
                  dchi(jj,ll) = &
                  (z(j)*dpsi(jj,ll,ng)-xxex(jj,n)*sdpsi(l))/sumx
                enddo
              enddo

!-------------jacobian
              do jj = nex1(k), nex2(k)
                j = jex(jj)
                do ll = nex1(k), nex2(k)
                  dpsi(jj,ll,ng) = dchi(jj,ll)
                  l = jex(ll)
!                 cln = cc(l,n)/z(j)
                  cln = ccloc_p(l+(ng-1)*ncomp)/z(j)
                  cdl(j,l) = cdl(j,l)+facec*dpsi(jj,ll,ng)*cln
                  if (ndecay.gt.0) then
                    cdl(j,l) = cdl(j,l)+facec*xlamhalf(j)* &
                    dpsi(jj,ll,ng)*cln
                  endif
                enddo
              enddo
            endif

!-----------set advective terms for sorption on mobile colloids
!           if (mcyc.gt.0 .and. m.gt.nexsolid-ncollex) then

!-------------set residual and jacobian for species sorbed on colloids
!             call trjaccol (nmat,nc,nexs,nexx,nelmrow,cc,xex,xxex,cec, &
!             dpsi,r,cjac,cdl,vl,trdifl,area,cxm,cxp, &
!             dfrc,ndcon,maxnc,ncdiag,nrow,k,ncex)

!-------------set boundary conditions for species sorbed on colloids
!             call trbccol (nmat,nc,nexs,nexx,xxex,dpsi,dpril,cec, &
!             cc,cdl,r,vlbc,ncc,k,ncex)
!           endif
          enddo ! end loop over sites -k-
        enddo ! end loop over exchange solids including mobile colloids-m-

  99    continue

      enddo ! end loop over nodes -n-

  100 continue

      if(nsrfmin.eq.0) return

!=======================================================================
!     surface complexation model
!=======================================================================

!-----variables                  description
!-----------------------------------------------------------------------
!     nsrfmin                    no. of sorption minerals and colloids
!     ncolsrf                    no. of sorption colloids
!     nsite1, nsite2             range of site indices
!     nsorp1, nsorp2             range of species indices
!     msorp                      namsrf(m) = namkin(msorp(m))
!     csorpf                     unoccupied site concentration
!     csorp, ccsorp              sorbed concentration
!     alogpf                     double layer potential
!     sitede0, siteden           site density
!-----------------------------------------------------------------------
!     pseudo code do-loop construction
!     do m = 1, nsrfmin                ! minerals/colloids
!       do k = nsite1(m), nsite2(m)    ! sites
!         do i = nsorp1(k), nsorp2(k)  ! sorbed species
!         enddo
!       enddo
!     enddo
!-----------------------------------------------------------------------

!-----loop over nodes
      do n = 1, nlmax

        ng = nL2G(n)

        facvn = vb(n)/dt
        facv  = facvn
        flamn = vb(n)
        flam  = flamn

!-------conversion factor for molarity to molality
        fach2o = one
!       if (molal.eq.0) fach2o = one/(wh2o*cc(jh2o,n))

        do j = 1, ncomp
!         alogjn(j) = log(cc(j,n)*gam(j,n)*fach2o)
          alogjn(j) = log(ccloc_p(j+(ng-1)*ncomp)*gam(j,ng)*fach2o)
          if (ncolsrf .gt. 0) then
            do l = 1, ncomp
              dpsi(j,l,ng) = zero
            enddo
          endif
        enddo

!-------loop over mineral and colloid surfaces
        do m = 1, nsrfmin

!         if (m.gt.nsrfmin-ncolsrf) then ! colloids
!           facv = porloc_p(ng)*facvn
!           flam = porloc_p(ng)*flamn
!         endif

          do k = nsite1(m), nsite2(m)

!-----------check for zero sorption sites - loop over sites
            if (siteden(k,n) .le. zero) then
              csorpf(k,n) = zero
              do i = nsorp1(k), nsorp2(k)
                ccsorp(i,n) = zero
                csorp(i,n)  = zero
              enddo
              alogpf(m,n) = zero
              goto 101
            endif

!           if (idblpot.ge.1) then

!-------------compute double layer potential and sorbed concentrations
!             itdblmax = 0
!             call dblpot(ncomp,ncmplx,nsrfmin,nsrfmx,nsrfsit,n,m,cc, &
!             cx,gam,ccsorp,csorpf,siteden,alogpf,alogjn)

!           else

!-------------local equilibrium: compute sorbed concentrations
!-------------first compute reduced sorption concentrations
              sum = zero
              do i = nsorp1(k), nsorp2(k)
                prod = eqsorp(i)*aln10
                do j = 1, ncomp
                  if (ssorp(j,i).ne.zero) then
                    prod = prod + ssorp(j,i)*alogjn(j)
                  endif
                enddo
                facp = exp(prod)
                ccsorp(i,n) = facp
                sum = sum + facp
              enddo

!-------------next compute free surface site density
              den = one + sum
              csorpf(k,n) = siteden(k,n)/den

!-------------finally get actual sorption concentrations
              do i = nsorp1(k), nsorp2(k)
                ccsorp(i,n) = csorpf(k,n)*ccsorp(i,n)
                
!               print *,'trionexc: ',i,k,namscx(i),ccsorp(i,n), &
!               csorpf(k,n)
              enddo
!           endif

!-----------compute residual and jacobian
            do l = 1, ncomp
              sum = zero
              do i = nsorp1(k), nsorp2(k)
                sum = sum + ssorp(l,i)*ccsorp(i,n)
              enddo
              dcsorp(l,k) = sum
            enddo

            do j = 1, ncomp

!-------------residual
              resjn1 = zero
              resjn2 = zero
              do i = nsorp1(k), nsorp2(k)
                resjn2 = resjn2 + ssorp(j,i)*ccsorp(i,n)
                resjn1 = resjn1 + ssorp(j,i)*csorp(i,n)
              enddo
!             r(j+(n-1)*nmat) = r(j+(n-1)*nmat) + facv*(resjn2-resjn1)
              b_p(j+(n-1)*nmat) = b_p(j+(n-1)*nmat) - facv*(resjn2-resjn1)

!-------------radioactive decay of sorbed species
!             if (ndecay.gt.0) then
!               r(j+(n-1)*nmat) = r(j+(n-1)*nmat) +
!               b_p(j+(n-1)*nmat) = b_p(j+(n-1)*nmat) - &
!               flam*xlamhalf(j)*resjn2
!             endif

!-------------jacobian
              if (m.le.nsrfmin-ncolsrf) then ! minerals
                do l = 1, ncomp
                  dcsc = zero
                  do i = nsorp1(k), nsorp2(k)
                    dcsc = dcsc + ssorp(j,i)*ssorp(l,i)*ccsorp(i,n)
                  enddo
                  dcsc = dcsc - dcsorp(j,k)*dcsorp(l,k)/siteden(k,n)
                  cdl(j,l) = cdl(j,l) + facv*dcsc

!-----------------radioactive decay
                  if (ndecay.gt.0) then
                    cdl(j,l) = cdl(j,l) + flam*dcsc*xlamhalf(j)
                  endif
                enddo
              else                           ! colloids
                do l = 1, ncomp
                  dcsc = zero
                  do i = nsorp1(k), nsorp2(k)
                    dcsc = dcsc + ssorp(j,i)*ssorp(l,i)*ccsorp(i,n)
                  enddo
                  dcsc = dcsc - dcsorp(j,k)*dcsorp(l,k)/siteden(k,n)
                  cdl(j,l) = cdl(j,l) + facv*dcsc

!-----------------store dpsi for colloid transport
                  dpsi(j,l,ng) = dpsi(j,l,ng) + dcsc

!-----------------radioactive decay
                  if (ndecay.gt.0) then
                    cdl(j,l) = cdl(j,l) + flam*dcsc*xlamhalf(j)
                  endif
                enddo
              endif
            enddo ! end loop over primary species j-loop
          enddo ! end loop over sites             k-loop
  101     continue
        enddo ! end loop over minerals/colloids   m-loop
      enddo ! end loop over nodes                 n-loop

!-----compute residual and jacobian for colloid contribution to flux
!-----set advective terms for sorption on mobile colloids

!     if (mcyc.gt.0 .and. ncolsrf.gt.0) then
!-------set residual and jacobian for species sorbed on colloids
!       call trjacsrf (nmat,nc,nsfs,nsfx,nelmrow,cc,ccsorp, &
!       dpsi,r,cjac,cdl,vl,trdifl,area,cxm,cxp, &
!       dfrc,ndcon,maxnc,ncdiag,nrow)

!-------set boundary conditions for species sorbed on collides
!       call trbcsrf (nmat,nc,nsfx,ccsorp,dpsi,cc,cdl,r,vlbc)
!     endif

! 500 continue
!-----constant Kd
!     do n = nmax1, nmax2
!       facv = vb(n)/dt
!       do j = 1, ncomp
!         r(j+(n-1)*nmat) = r(j+(n-1)*nmat) + 
!         facv*distcoef(j)*(cc(j,n)-c(j,n))
!         do l = 1, ncomp
!           if (j.eq.l) cdl(j,l) = cdl(j,l) + facv*distcoef(j)
!         enddo
!       enddo
!     enddo

      ibug = 0
      if (ibug .eq. 1) then
        do m = 1, nsrfmin
          do l = nsite1(m), nsite2(m)
            do i = nsorp1(l), nsorp2(l)
              print *,'trionex: ',i,namscx(i),' ',ccsorp(i,1)
            enddo
          enddo
        enddo
      endif

      return
  end subroutine trionexc

end module trionexc_module
