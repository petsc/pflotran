!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_multi.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_multi.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:34:14  lichtner
! Revised boundary conditions for ibndtyp=3 and source/sink term.
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

!RTM: ptran_multi() sets up the linear system that is solved at each 
!RTM: Newton step.  A is the Jacobian, b is the right hand side (b_p is a 
!RTM: pointer to b).  Thus b is -R^{(i)}_{jn} in the discretization.

module ptran_multi_module

contains

  subroutine ptran_multi (its,ccloc_p,temploc_p,porloc_p,sat_loc_p,ssat_loc_p)

  use ptran_global_module
  use trdynmem_module
  use trkinmin_module
  
  use water_eos_module

  implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"
!#include "include/finclude/petscda.h"

  integer :: its,j,jm,jm1,jm2,jn,jng,l,lm1,lm2,lng, &
             m,m1,m2,mm1,mm2,n,n1,n2,nc,ng,nr
  integer :: ierr,i,ii,ii1,ii2,jj,jj1,jj2,kk,kk1,kk2
! integer :: irow1(ncomp),irow2(ncomp),icol1(ncomp),icol2(ncomp)
      
  real*8 :: trans, qflux, val, val1, val2, voldt, wpor, por12, q, u1, u2
  real*8 :: csrc1,tsrc1,qsrc1,qqsrc,sl,ssl,sl1,sl2,sg,ssg,sg1,sg2,pvdt
  real*8 :: ccloc_p(*),temploc_p(*),porloc_p(*),sat_loc_p(*),ssat_loc_p(*)
  real*8 :: ff,f1,f2,dw_kg,dw_mol,enth,scalefac=1.d-6
  real*8 :: por1,por2,tor1,tor2

  real*8 :: blkmat1(ncomp,ncomp),blkmat2(ncomp,ncomp)
  
  real*8, pointer :: press_p(:),tort_loc_p(:)

  ierr = 0

! call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)

  call VecGetArrayF90(b,b_p,ierr) ! note: b = -r

  call VecGetArrayF90(tort_loc,tort_loc_p,ierr)

!---diagonal accumulation term - use ghosted nodes for MatSetValue(s)Local
!   units: vb [m^3]; psi [mol/dm^3]; b [kmol/s]
    if (iphase == 1) then ! pure liquid
      do n = 1, nlmax
        ng = nL2G(n)
        voldt = vb(n)/dt
        pvdt = sat0*porloc_p(ng)*voldt
        do j = 1, ncomp
          jn = j+(n-1)*nmat
          jng = j+(ng-1)*nmat
!         irow1(j) = jng-1
          b_p(jn) = -pvdt*(ppsi(j,ng)-psi(j,ng))
          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1
            val = pvdt*dpsi(j,l,ng)
            blkmat1(j,l) = val
            
!           if (t/yrsec > 48.5) &
!           print *,'ptranmulti-acc: ',myrank,n,ng,porloc_p(ng),vb(n), &
!           psi(j,ng),ppsi(j,ng),dpsi(j,l,ng),b_p(jn)
            
            if (iblkfmt == 0 .and. dfill(j+(l-1)*ncomp).ne.0) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val,ADD_VALUES,ierr)
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val,ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat1, &
          ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,blkmat1, &
!         ADD_VALUES,ierr)
        endif
      enddo

    else if (iphase == 2) then ! two phase gas-liquid

      do n = 1, nlmax
        ng = nL2G(n)
        voldt = vb(n)/dt
        pvdt = porloc_p(ng)*voldt
        ssl = ssat_loc_p(ng)
        ssg = 1.d0 - ssl
        sl = sat_loc_p(ng)
        sg = 1.d0 - sl
        
!       print *,'ptran_multi: ',n,ng,dt,ssl,sl
        
        do j = 1, ncomp
          jn = j+(n-1)*nmat
          jng = j+(ng-1)*nmat
!         irow1(j) = jng-1
          b_p(jn) = -pvdt*(ssl*ppsi(j,ng)-sl*psi(j,ng) + &
          ssg*ppsig(j,ng)-sg*psig(j,ng))
       !  if(n==400)&
       !  print *,'ptranmulti-acc: ',myrank,n,ng,j,porloc_p(ng),dt,pvdt,vb(n), &
       !  ssl,ssg,sl,sg,psi(j,ng),ppsi(j,ng),psig(j,ng),ppsig(j,ng),b_p(jn)

          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1
            val = pvdt*(ssl*dpsi(j,l,ng)+ssg*dpsig(j,l,ng))

!           print *,'ptranmulti-acc: ',myrank,n,ng,j,l,porloc_p(ng),vb(n), &
!           ssl,ssg,sl,sg,psi(j,ng),ppsi(j,ng),dpsi(j,l,ng), &
!           dpsig(j,l,ng),b_p(jn)

            blkmat1(j,l) = val
            if (iblkfmt == 0 .and. dfill(j+(l-1)*ncomp).ne.0) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val,ADD_VALUES,ierr)
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val,ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat1, &
          ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,blkmat1, &
!         ADD_VALUES,ierr)
        endif
      enddo
    endif

!---diagonal/off-diagonal flux terms - interior nodes
!   m = local/ghosted node, n = local node
! RTM: What we do here is loop over all of the interfaces between nodes.
! RTM: For each interface, we calculate the flux across its contribution 
! RTM: to the residual at both of the nodes that share the interface.
    do nc = 1, nconn  !RTM: We loop over the number of node connections.

      m1 = nd1(nc) !RTM: m1, m2 are local ghosted indices.
      m2 = nd2(nc) !RTM: Cell m2 is downstream of m1 (if flow is positive).

      n1 = nG2L(m1) !RTM: n1, n2 are local non-ghosted indices.
      n2 = nG2L(m2)
      
      tor1 = tort_loc_p(m1)
      tor2 = tort_loc_p(m2)
      por1 = porloc_p(m1)
      por2 = porloc_p(m2)

      if (iphase == 1) then
        wpor  = tor1*por1*dist2(nc)+tor2*por2*dist1(nc)
        por12 = sat0*tor1*tor2*por1*por2/wpor
      else
        sl1 = ssat_loc_p(m1)
        sl2 = ssat_loc_p(m2)
        
        if (sl1 >= slcutoff .and. sl2 >=  slcutoff) then
          wpor  = tor1*por1*sl1*dist2(nc)+tor2*por2*sl2*dist1(nc)
          por12 = tor1*tor2*por1*por2*sl1*sl2/wpor
        else
          por12 = 0.D0
        endif
      endif
      trans = por12*difaq*area(nc)

      q = vl(nc)*area(nc)
      
!     print *,'ptran_multi-l: ',nc,q,vl(nc),area(nc),sl1,sl2

      !upstream weighting
      if (q > zero) then
        u1 =  trans+q
        u2 = -trans
      else
        u1 =  trans
        u2 = -trans+q
      endif
      
!     if (t/yrsec >= 1.0) &
!     print *,'ptranmulti-flx: ',myrank,nc,m1,m2,vl(nc)*yrsec
!     print *,'ptranmulti-flx: ',myrank,nc,m1,m2,vl(nc)*yrsec,trans,sat0, &
!              area(nc),porloc_p(m1),porloc_p(m2),ppsi(1,m1),ppsi(1,m2)

      if (n1 > 0) then !RTM: If the upstream node is not a ghost node... 
        
        do j = 1, ncomp
          jm1 = j+(m1-1)*nmat
          jm2 = j+(m2-1)*nmat

          mm1 = jm1-1
          mm2 = jm2-1
!         irow1(j) = mm1
!         irow2(j) = mm2

          qflux = u2*ppsi(j,m2)+u1*ppsi(j,m1)

          jn = j+(n1-1)*nmat
          b_p(jn) = b_p(jn) - qflux ! note: b = -r
          
!         print *,'ptran_multi1: ',myrank,n1,jn,b_p(jn),qflux

          do l = 1, ncomp
            lm1 = l+(m1-1)*nmat
            lm2 = l+(m2-1)*nmat
!           icol1(l) = lm1-1
!           icol2(l) = lm2-1
            val1 = u1*dpsi(j,l,m1)
            val2 = u2*dpsi(j,l,m2)
            if (iblkfmt == 0) then
              if (ofill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,mm1,1,lm1-1,val1, &
                ADD_VALUES,ierr)
                call MatSetValuesLocal(A,1,mm1,1,lm2-1,val2, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val1
              blkmat2(j,l) = val2
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,mm1,1,lm1-1,val1, &
              ADD_VALUES,ierr)
              call MatSetValuesLocal(A,1,mm1,1,lm2-1,val2, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,m1-1,1,m1-1,blkmat1, &
          ADD_VALUES,ierr)
          call MatSetValuesBlockedLocal(A,1,m1-1,1,m2-1,blkmat2, &
          ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,blkmat1, &
!         ADD_VALUES,ierr)
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol2,blkmat2, &
!         ADD_VALUES,ierr)
        endif
      endif

      if (n2 > 0) then !RTM: If the downstream node is not a ghost node...
        do j = 1, ncomp
          jm1 = j+(m1-1)*nmat
          jm2 = j+(m2-1)*nmat

          mm1 = jm1-1
          mm2 = jm2-1
!         irow1(j) = mm1
!         irow2(j) = mm2

          qflux = u2*ppsi(j,m2)+u1*ppsi(j,m1)

          jn = j+(n2-1)*nmat
          b_p(jn) = b_p(jn) + qflux
          
!         print *,'ptran_multi2: ',myrank,n2,jn,b_p(jn),qflux
          
          do l = 1, ncomp
            lm1 = l+(m1-1)*nmat
            lm2 = l+(m2-1)*nmat
!           icol1(l) = lm1-1
!           icol2(l) = lm2-1
            val1 = -u1*dpsi(j,l,m1)
            val2 = -u2*dpsi(j,l,m2)
            if (iblkfmt == 0) then
              if (ofill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,mm2,1,lm2-1,val2, &
                ADD_VALUES,ierr)
                call MatSetValuesLocal(A,1,mm2,1,lm1-1,val1, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val1
              blkmat2(j,l) = val2
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,mm2,1,lm2-1,val2, &
              ADD_VALUES,ierr)
              call MatSetValuesLocal(A,1,mm2,1,lm1-1,val1, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,m2-1,1,m2-1,blkmat2, &
          ADD_VALUES,ierr)
          call MatSetValuesBlockedLocal(A,1,m2-1,1,m1-1,blkmat1, &
          ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow2,ncomp,icol2,blkmat2, &
!         ADD_VALUES,ierr)
!         call MatSetValuesLocal(A,ncomp,irow2,ncomp,icol1,blkmat1, &
!         ADD_VALUES,ierr)
        endif
      endif

      if (iphase == 2) then ! two phase liquid-gas system

        sg1 = 1.d0-ssat_loc_p(m1)
        sg2 = 1.d0-ssat_loc_p(m2)
        wpor = porloc_p(m1)*sg1*dist2(nc)+porloc_p(m2)*sg2*dist1(nc)
 
        if(wpor > 1.D-8)then
          por12 = porloc_p(m1)*sg1*porloc_p(m2)*sg2/wpor
          trans = por12*difgas*area(nc)
        else
          trans=0.D0
        endif
      
        q = vg(nc)*area(nc)

!       if (nc <= 10) &
!       print *,'ptran_multi-g: ',nc,vl(nc)*tconv,vg(nc)*tconv,sg1,sg2

        !upstream weighting
        if (q > zero) then
          u1 =  trans+q
          u2 = -trans
        else
          u1 =  trans
          u2 = -trans+q
        endif

        if (n1 > 0) then !RTM: If the upstream node is not a ghost node... 
        
          do j = 1, ncomp
            jm1 = j+(m1-1)*nmat
            jm2 = j+(m2-1)*nmat

            mm1 = jm1-1
            mm2 = jm2-1
!           irow1(j) = mm1
!           irow2(j) = mm2

            qflux = u2*ppsig(j,m2)+u1*ppsig(j,m1)

            jn = j+(n1-1)*nmat
            b_p(jn) = b_p(jn) - qflux ! note: b = -r
		
	  ! if(j==jco2)&
!       print *,'ptran_multi: ',nc,jn,q,difgas,qflux,b_p(jn), &
!		u2,ppsig(j,m2),u1,ppsig(j,m1)

            do l = 1, ncomp
              lm1 = l+(m1-1)*nmat
              lm2 = l+(m2-1)*nmat
!             icol1(l) = lm1-1
!             icol2(l) = lm2-1
              val1 = u1*dpsig(j,l,m1)
              val2 = u2*dpsig(j,l,m2)
              if (iblkfmt == 0) then
                if (ofill(j+(l-1)*ncomp).ne.0) then
                  call MatSetValuesLocal(A,1,mm1,1,lm1-1,val1, &
                  ADD_VALUES,ierr)
                  call MatSetValuesLocal(A,1,mm1,1,lm2-1,val2, &
                  ADD_VALUES,ierr)
                endif
              else if (iblkfmt == 1) then
                blkmat1(j,l) = val1
                blkmat2(j,l) = val2
              else if (iblkfmt == 2) then
                call MatSetValuesLocal(A,1,mm1,1,lm1-1,val1, &
                ADD_VALUES,ierr)
                call MatSetValuesLocal(A,1,mm1,1,lm2-1,val2, &
                ADD_VALUES,ierr)
              endif
            enddo
          enddo
          if (iblkfmt == 1) then
            call MatSetValuesBlockedLocal(A,1,m1-1,1,m1-1,blkmat1, &
            ADD_VALUES,ierr)
            call MatSetValuesBlockedLocal(A,1,m1-1,1,m2-1,blkmat2, &
            ADD_VALUES,ierr)
!         else
!           call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,blkmat1, &
!           ADD_VALUES,ierr)
!           call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol2,blkmat2, &
!           ADD_VALUES,ierr)
          endif
        endif

        if (n2 > 0) then !RTM: downstream node is not a ghost node...
          do j = 1, ncomp
            jm1 = j+(m1-1)*nmat
            jm2 = j+(m2-1)*nmat

            mm1 = jm1-1
            mm2 = jm2-1
!           irow1(j) = mm1
!           irow2(j) = mm2

            qflux = u2*ppsig(j,m2)+u1*ppsig(j,m1)

            jn = j+(n2-1)*nmat
            b_p(jn) = b_p(jn) + qflux
          
            do l = 1, ncomp
              lm1 = l+(m1-1)*nmat
              lm2 = l+(m2-1)*nmat
!             icol1(l) = lm1-1
!             icol2(l) = lm2-1
              val1 = -u1*dpsig(j,l,m1)
              val2 = -u2*dpsig(j,l,m2)
              if (iblkfmt == 0) then
                if (ofill(j+(l-1)*ncomp).ne.0) then
                  call MatSetValuesLocal(A,1,mm2,1,lm2-1,val2, &
                  ADD_VALUES,ierr)
                  call MatSetValuesLocal(A,1,mm2,1,lm1-1,val1, &
                  ADD_VALUES,ierr)
                endif
              else if (iblkfmt == 1) then
                blkmat1(j,l) = val1
                blkmat2(j,l) = val2
              else if (iblkfmt == 2) then
                call MatSetValuesLocal(A,1,mm2,1,lm2-1,val2, &
                ADD_VALUES,ierr)
                call MatSetValuesLocal(A,1,mm2,1,lm1-1,val1, &
                ADD_VALUES,ierr)
              endif
            enddo
          enddo
          if (iblkfmt == 1) then
            call MatSetValuesBlockedLocal(A,1,m2-1,1,m2-1,blkmat2, &
            ADD_VALUES,ierr)
            call MatSetValuesBlockedLocal(A,1,m2-1,1,m1-1,blkmat1, &
            ADD_VALUES,ierr)
!         else
!           call MatSetValuesLocal(A,ncomp,irow2,ncomp,icol2,blkmat2, &
!           ADD_VALUES,ierr)
!           call MatSetValuesLocal(A,ncomp,irow2,ncomp,icol1,blkmat1, &
!           ADD_VALUES,ierr)
          endif
        endif
      endif
    enddo

!---add boundary conditions

    if (nblkbc == 0) goto 100

    do nc = 1, nconnbc

      m = mblkbc(nc) ! m not ghosted
      ng = nL2G(m)

      if (iphase == 1) then
        trans = porloc_p(ng)*sat0*difaq/distbc(nc)*areabc(nc)
      else
        trans = porloc_p(ng)*ssat_loc_p(ng)*difaq/distbc(nc)*areabc(nc)
      endif
      
      ! q > 0 if flow is INTO block
      q = vlbc(nc)*areabc(nc)

!-----locate boundary face
      ibc = ibconn(nc)
      
!     if (t/yrsec > 48.5) &
!     print *,'ptranmulti-bnd: ',myrank,nc,m,ng,ibc,ibndtyp(ibc), &
!     distbc(nc),areabc(nc)
!     vlbc(nc),trans,sat0,porloc_p(ng),ppsi(1,ng),psibnd(1,ibc)

      if (ibndtyp(ibc) == 1) then
        ! RTM: The if-else below takes care of upwinding.
        if (q > zero) then ! flow into block
          u1 =  trans
          u2 = -trans-q
        else
          u1 =  trans-q
          u2 = -trans
        endif
        do j = 1, ncomp
          jng = j+(ng-1)*nmat
          jm = j+(m-1)*nmat
!         irow1(j) = jng-1

          qflux = u1*ppsi(j,ng)+u2*psibnd(j,ibc)
          
!         print *,'multi-l: ',nc,ng,ibc,q,trans,psibnd(j,ibc),ppsi(j,ng)

          b_p(jm) = b_p(jm) - qflux
          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1
            val = u1*dpsi(j,l,ng)
            if (iblkfmt == 0) then
              if (dfill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
          blkmat1,ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1, &
!         blkmat1,ADD_VALUES,ierr)
        endif
        
      else if (ibndtyp(ibc) == 3 .and. q.ne.zero) then ! zero gradient
      
       if (q > zero) then ! flow into block
         u1 =  zero
         u2 = -q
       else
         u1 = -q
         u2 =  zero
       endif
      
        do j = 1, ncomp
          jng = j+(ng-1)*nmat
          jm = j+(m-1)*nmat
!         irow1(j) = jng-1

         qflux = u1*ppsi(j,ng)+u2*psibnd(j,ibc)
!          qflux = -q*ppsi(j,ng)

          b_p(jm) = b_p(jm) - qflux
          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1

           val = u1*dpsi(j,l,ng)
!            val = -q*dpsi(j,l,ng)
      
!     if (t/yrsec > 48.5) &
!     print *,'ptranmulti-bnd3: ',nc,m,ng,ibc,vlbc(nc),trans,ibndtyp(ibc), &
!     porloc_p(ng),val,b_p(jm),qflux,ppsi(j,ng),dpsi(j,l,ng)
            
            if (iblkfmt == 0) then
              if (dfill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
          
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
          blkmat1,ADD_VALUES,ierr)
!        else
!          call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1, &
!          blkmat1,ADD_VALUES,ierr)
        endif
      endif
    enddo

  if (iphase == 2) then

    do nc = 1, nconnbc

      m = mblkbc(nc) ! m not ghosted
      ng = nL2G(m)

      trans = porloc_p(ng)*(1.d0-ssat_loc_p(ng)) &
      *difgas/distbc(nc)*areabc(nc)

      ! q > 0 if flow is INTO block
      q = vgbc(nc)*areabc(nc)

!-----locate boundary face
      ibc = ibconn(nc)

      if (ibndtyp(ibc) == 1) then
        ! RTM: The if-else below takes care of upwinding.
        if (q > zero) then ! flow into block
          u1 =  trans
          u2 = -trans-q
        else
          u1 =  trans-q
          u2 = -trans
        endif
        do j = 1, ncomp
          jng = j+(ng-1)*nmat
          jm = j+(m-1)*nmat
!         irow1(j) = jng-1
          qflux = u1*ppsig(j,ng)+u2*psigbnd(j,ibc)
	!	  print *,'Ptran-multi: 1 : psaig',q,trans,ng, j,ibc, u1,u2,ppsig(j,ng),psigbnd(j,ibc)
          !print *, xphibc(nc)
          b_p(jm) = b_p(jm) - qflux
          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1
            val = u1*dpsig(j,l,ng)
          
            if (iblkfmt == 0) then
              if (dfill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
          blkmat1,ADD_VALUES,ierr)
!       else
!         call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1, &
!         blkmat1,ADD_VALUES,ierr)
        endif
        
      else if (ibndtyp(ibc) == 3 .and. q.ne.0.D0) then ! zero gradient
      
       if (q > 0.D0) then ! flow into block
         u1 =  0.D0
         u2 = -q
       else
         u1 = -q
         u2 =  0.D0
       endif
      
        do j = 1, ncomp
          jng = j+(ng-1)*nmat
          jm = j+(m-1)*nmat
!         irow1(j) = jng-1

         qflux = u1*ppsig(j,ng)+u2*psigbnd(j,ibc)
	!	 print *,'Ptran-multi: 3 : psaig', q, ng, j,ibc, u1,u2,ppsig(j,ng),psigbnd(j,ibc)
!          qflux = -q*ppsig(j,ng)
          
          b_p(jm) = b_p(jm) - qflux
          do l = 1, ncomp
            lng = l+(ng-1)*nmat
!           icol1(l) = lng-1
           val = u1*dpsig(j,l,ng)
!            val = -q*dpsig(j,l,ng)
            if (iblkfmt == 0) then
              if (dfill(j+(l-1)*ncomp).ne.0) then
                call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
                ADD_VALUES,ierr)
              endif
            else if (iblkfmt == 1) then
              blkmat1(j,l) = val
            else if (iblkfmt == 2) then
              call MatSetValuesLocal(A,1,jng-1,1,lng-1,val, &
              ADD_VALUES,ierr)
            endif
          enddo
        enddo
          
        if (iblkfmt == 1) then
          call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1, &
          blkmat1,ADD_VALUES,ierr)
!        else
!          call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1, &
!          blkmat1,ADD_VALUES,ierr)
        endif
      endif
    enddo
  endif
  
  100 continue
  
!---add source/sink terms

    do nr = 1, nblksrc
      
      kk1 = ks1(nr) - nzs
      kk2 = ks2(nr) - nzs
      jj1 = js1(nr) - nys
      jj2 = js2(nr) - nys
      ii1 = is1(nr) - nxs
      ii2 = is2(nr) - nxs
        
      kk1 = max(1,kk1)
      kk2 = min(nlz,kk2)
      jj1 = max(1,jj1)
      jj2 = min(nly,jj2)
      ii1 = max(1,ii1)
      ii2 = min(nlx,ii2)
        
      if (ii1 > ii2 .or. jj1 > jj2 .or. kk1 > kk2) cycle
      
      do i = 2, nsrc
        if (timesrc(i,nr) == t) then
          qsrc1 = qsrc(i,nr)
          csrc1 = csrc(i,nr)
          goto 10
        else if (timesrc(i,nr) > t) then
          ff = timesrc(i,nr)-timesrc(i-1,nr)
          f1 = (t - timesrc(i-1,nr))/ff
          f2 = (timesrc(i,nr)-t)/ff
          tsrc1 = f1*tempsrc(i,nr) + f2*tempsrc(i-1,nr)
          qsrc1 = f1*qsrc(i,nr) + f2*qsrc(i-1,nr)
          csrc1 = f1*csrc(i,nr) + f2*csrc(i-1,nr)
          goto 10
        endif
      enddo
 10   continue
    
!     print *,'ptranmulti: ', myrank,i,timesrc(i,nr), &
!     timesrc(i-1,nr),t,f1,f2,ff,qsrc1,csrc1

      call VecGetArrayF90(press,press_p,ierr)
        
      if (qsrc1 > zero) then ! injection of H2O
        do kk = kk1,kk2
          do jj = jj1,jj2
            do ii = ii1,ii2
              n = ii+(jj-1)*nlx+(kk-1)*nlxy
              ng = nL2G(n)

              call wateos_noderiv(tsrc1,press_p(n),dw_kg,dw_mol, &
              enth,scalefac,ierr)
              
              !units: press_p [Pa], dw_kg [kg/m^3], dw_mol [mol/dm^3=kmol/m^3]

              qqsrc = qsrc1/dw_kg ! [kg/s / (kg/m^3) = m^3/s]
              do j = 1, ncomp
                jm = j+(n-1)*nmat
                b_p(jm) = b_p(jm) + qqsrc*psisrc(j,nr) ! [m^3/s mol/dm^3 = kmol/s]
              
!               print *,'ptranmulti-h2o: ',nr,n,ng,j,tsrc1,rho(ng), &
!               psisrc(j,nr),qsrc1,qqsrc,press_p(n),dw_kg,dw_mol
              
              enddo
            enddo
          enddo
        enddo
        
      else if (qsrc1 < zero) then ! withdrawal
      
        do kk = kk1,kk2
          do jj = jj1,jj2
            do ii = ii1,ii2
              n = ii+(jj-1)*nlx+(kk-1)*nlxy
              ng = nL2G(n)
              qqsrc = qsrc1/rho(ng)*1.d-3
              do j = 1, ncomp
                jm = j+(n-1)*nmat
                jng = j+(ng-1)*nmat-1
!               irow1(j) = jng
!               icol1(j) = jng
                b_p(jm) = b_p(jm) + qqsrc*ppsi(j,ng)
                do l = 1, ncomp
                  blkmat1(j,l) = -qqsrc*dpsi(j,l,ng)
                  lng = l+(ng-1)*ncomp-1
                  if (iblkfmt == 0 .and. dfill(j+(l-1)*ncomp).ne.0) then
                    call MatSetValuesLocal(A,1,jng,1,lng,blkmat1(j,l), &
                    ADD_VALUES,ierr)
                  else if (iblkfmt == 2) then
                    call MatSetValuesLocal(A,1,jng,1,lng,blkmat1(j,l), &
                    ADD_VALUES,ierr)
                  endif
                enddo
              enddo
              if (iblkfmt == 1) then
                call MatSetValuesBlockedLocal(A,1,ng-1,1,ng-1,blkmat1, &
                ADD_VALUES,ierr)
!             else
!               call MatSetValuesLocal(A,ncomp,irow1,ncomp,icol1,blkmat1, &
!               ADD_VALUES,ierr)
              endif
            enddo
          enddo
        enddo
      endif
        
      if (csrc1 > zero) then ! injection of pure CO2
        csrc1 = csrc1/fmwco2 ! [(kg/s) / (kg/kmol) = kmol/s]
        do kk = kk1,kk2
          do jj = jj1,jj2
            do ii = ii1,ii2
              n = ii+(jj-1)*nlx+(kk-1)*nlxy
              ng = nL2G(n)
              qqsrc = csrc1 ! [kmol/s] *1.d-3 !/rho(ng)*1.d-3
              do j = 1, ncomp
                jm = j+(n-1)*nmat
                if (j == jco2) &
                b_p(jm) = b_p(jm) + qqsrc !*psisrc(j,nr)
              
!               print *,'ptranmulti-co2: ',nr,n,ng,j,fmwco2,csrc1,qqsrc
              
              enddo
            enddo
          enddo
        enddo
      endif
    enddo ! loop nr over sources
    
!---add kinetic reaction rate terms
    if (nkin > 0) then
      if (isolidss == 0) then
        call trkinmin (ccloc_p,temploc_p)
      else if (isolidss > 0) then
        call trkinminss (ccloc_p,temploc_p)
      endif
    endif
    
!---add equilibrium sorption (ion exchange + surface complexation)
!   if (nsrfmx > 0 .or. nexmx > 0) call trionexc (ccloc_p,porloc_p)

    ! RTM: Assemble the Jacobian matrix A.
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

    call VecRestoreArrayF90(b,b_p,ierr)

    call VecRestoreArrayF90(press,press_p,ierr)
    call VecRestoreArrayF90(tortuosity,tort_loc_p,ierr)
    
!   call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
              
  end subroutine ptran_multi

end module ptran_multi_module
