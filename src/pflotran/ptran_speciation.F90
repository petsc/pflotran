!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_speciation.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_speciation.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:53:06  lichtner
! Set rho = 1 for debugging purposes.
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

module ptran_speciation_module

  public

contains

      subroutine trstartup (da,da_1dof,da_kin)

      use ptran_global_module
      use trdynmem_module
      use water_eos_module
      use co2eos_module
	  
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petsclog.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscviewer.h"

  DA    :: da, da_1dof, da_kin
      
      real*8 :: cloctot(ncmx),gamloc(ncmx),pgasloc(ngmx),cloc(ncmx), &
                cxloc(ncxmx),gamxloc(ncxmx),cecloc(nexmx),xexloc(nexmx),&
                sitdnloc(nsitmx),csorpflc(nscxmx),csorplc(nscxmx), &
                ccsorplc(nscxmx)
      
      real*8 :: ah2o0,pgas0(ngmx),profji,dgamdi(ncmx), &
                c0(ncmx),bc0(ncmx),psi0(ncmx),psig0(ncmx), &
                gmolfrac(ncmx,nrgmx)
      
      integer :: ierr,j,jj,i,k,kk,l,m,mcyc,mcyc0,n,nnreg,nruns0,irun, &
                 iflgcon,iireg,iflag,iter
      
      real*8  :: tc,u1, dco2,fc,phi
      
      real*8, pointer :: temp_p(:), nreg_p(:)

!**************************************************************** 
!                    initialize variables
!**************************************************************** 
      ierr = 0
      
      call VecGetArrayF90(nreg,nreg_p,ierr)
      call VecGetArrayF90(temp,temp_p,ierr)
      call VecGetArrayF90(c,c_p,ierr)
      
      if(irestart.le.1) t    = zero
      if(irestart.le.1) mcyc = 0

      call density(tempini,pref0,rho0)
      rho0 = rho0*1.d-3
      rho0 = 1.d0
      
      if (myrank==0) write(iunit2,109) tempini,rho0
  109 format(' temperature =',1pg12.4, &
     &       ' [C], rhoh2o =',1pg12.4,' [g/cm^3]')
     
      if (iact.eq.1 .or. iact.eq.5) then
        if (myrank==0) write(iunit2,110) adebye,bdebye,bextend(1)
      else if (iact.eq.2) then
        if (myrank==0) write(iunit2,111)
      endif
  110 format(' debye-huckel parameters:'/,' a =',1pg12.4,'b =',1pg12.4, &
     &       'bdot =',1pg12.4)
  111 format(' debye-huckel parameters: Davies algorithm')
  
       wm1 = one-w
       wlam1 = one-wlam
       
!***********************************************************************
!     read restart values 
!***********************************************************************
      mcyc0=0 
      if(irestart.ge.2) then
        if (myrank==0) &
        write(iunit2,*) ' restart data read from tape iunit4 '
        mcyc=mcyc0
        stop '>>>FloTran is not set-up for restart!'
      endif

!-----set initial guess for ionic strength <improve this!>
      sionic0 = 1.d-3
      do n = 1, ngmax
        if (ghost_loc_p(n) == -1) cycle
        sionic(n) = sionic0 
      enddo

      if(irestart.le.1) then
        do n = 1, nlmax
!         do j = 1, ncomp
!           psi(j,n) = zero 
!           gam(j,n) = one
!         enddo
          do m = 1, nexsolid
            do k = nsitex1(m), nsitex2(m)
              do j = nex1(k), nex2(k)
                xex(j,n)  = zero
                xxex(j,n) = zero
              enddo
            enddo
          enddo
          if(iphase.eq.2 .or. iphase.eq.0) then
            do j = 1, ncomp
              psig(j,n) = zero
            enddo
          endif

!         do j = 1, ngas
!           pgas(j,n) = zero
!         enddo

!         do i = 1, ncmplx
!           cx(i,n)   = zero
!           gamx(i,n) = one
!         enddo
          if (iact.eq.6) then
            ah2o(n) = zero
          endif
        enddo
      endif

!****************************************************************
!              setup initial conditions 
!****************************************************************

!-----initialize activity coefficients
      do j = 1, ncomp
        gamloc(j) = one
      enddo
      do i = 1, ncmplx
        gamxloc(i) = one
      enddo
      if (iact.eq.6) ah2o0 = zero

!-----determine stoichiometry of conserved quantities
!     if (nconsrv .gt. 0) call trconsrv

!-----initialize cation exchange capacity for all regions
      if (nexsolid .gt. 0) then
        if (myrank == 0) &
        write(iunit2,'(/," region       mineral",14x,"site     cec" &
     &  )')
!       do l = 1, 3
!         ll = 2*l-1
!         write(*,*) 'trstrtup: ',l,ll,ll+1,ibcreg(ll), ibcreg(ll+1)
!         if (ibcreg(ll).gt.0 .and. ibcreg(ll+1).gt.0) then
!           do ireg = ibcreg(ll), ibcreg(ll+1)
            do ireg = 1, initreg
              do m = 1, nexsolid
                do k = nsitex1(m), nsitex2(m)
                  if (idcdm.eq.0) then
                    if (ndloc1.eq.1) then             ! single continuum
                      cec0(ireg,k) = cec0mf(1,k)
                    else                              ! dcm
                      if (mod(ireg,2).ne.0) then
                        cec0(ireg,k) = cec0mf(1,k)    ! matrix
                      else
                        cec0(ireg,k) = cec0mf(2,k)    ! fracture
                      endif
                    endif
                  else
                    if (ireg.eq.1) then
                      cec0(ireg,k) = cec0mf(2,k)      ! fracture
                    elseif (ireg.eq.2) then
                      cec0(ireg,k) = cec0mf(1,k)      ! matrix
                    elseif (ireg.eq.3 .or. ireg.eq.4) then
                      cec0(ireg,k) = cec0mf(2,k)      ! fracture
                    endif
                  endif
                  if (myrank == 0) &
                  write(iunit2,'(i4,14x,a20," site ",i2,1pe12.4)') &
                  ireg,namex(m),k,cec0(ireg,k)
                enddo
              enddo
            enddo
!         endif
!       enddo
      endif

!-----initialize surface site density for all regions
      if (nsrfmin .gt. 0) then
        if (myrank == 0) &
        write(iunit2,'(/," region       mineral",14x,"site     cec" &
        &)')
!       do l = 1, 3
!         ll = 2*l-1
!         write(*,*) 'trstrtup: ',l,ll,ll+1,ibcreg(ll), ibcreg(ll+1)
!         if (ibcreg(ll).gt.0 .and. ibcreg(ll+1).gt.0) then
!           do ireg = ibcreg(ll), ibcreg(ll+1)
            do ireg = 1, initreg
              do m = 1, nsrfmin
                do k = nsite1(m), nsite2(m)
                  if (idcdm.eq.0) then
                    if (ndloc1.eq.1) then             ! single continuum
                      siteden0(ireg,k) = sited0mf(1,k)
                    else                              ! dcm
                      if (mod(ireg,2).ne.0) then
                        siteden0(ireg,k) = sited0mf(1,k)    ! matrix
                      else
                        siteden0(ireg,k) = sited0mf(2,k)    ! fracture
                      endif
                    endif
                  else
                    if (ireg.eq.1) then
                      siteden0(ireg,k) = sited0mf(2,k)      ! fracture
                    elseif (ireg.eq.2) then
                      siteden0(ireg,k) = sited0mf(1,k)      ! matrix
                    elseif (ireg.eq.3 .or. ireg.eq.4) then
                      siteden0(ireg,k) = sited0mf(2,k)      ! fracture
                    endif
                  endif
                  if (myrank == 0) &
                  write(iunit2,'(i4,14x,a20," site ",i2,1pe12.4)') &
                  ireg,namex(m),k,siteden0(ireg,k)
                enddo
              enddo
            enddo
!         endif
!       enddo
      endif
      
      nruns0 = nruns
      if (modetr .eq. 0) then
        if (nruns0 .eq. 0) nruns0 = 1
      else
        nruns0 = 1
      endif

      if (isothrm .ge. 1) goto 100

!-----isothermal system

!     write(*,*) 'trstrtup: ',ibcreg(1),ibcreg(2),ndloc1,isothrm,nlmax

!     do ireg = ibcreg(1), ibcreg(2)
      do ireg = 1, initreg
      
!       write(*,*) 'trstrtup: ',ireg,initreg

!-------PROFile multiple runs
        do irun = 1, nruns0
      
          if (nruns0.gt.1) then
            do j = 1, ncomp
              profji = profile(j,irun)
              if (profji.gt.zero) ctot(j,ireg) = profji
            enddo
          endif
          
        tk = tempini+tkelvin

        if (myrank == 0) then
          write(iunit2,'(/,36("=+"))')
          if (ndloc1.eq.2) then
            if (mod(ireg,2).ne.0) then
              write(iunit2,'(10x,"initial conditions: region",i3, &
     &        " matrx")') ireg
            else
              write(iunit2,'(10x,"initial conditions: region",i3, &
     &        " frac")') ireg
            endif
          else if (idcdm.eq.1) then
            if (ireg.eq.2) then
              write(iunit2,'(10x,"initial conditions: region",i3, &
     &        " matrx")') ireg
            else
              write(iunit2,'(10x,"initial conditions: region",i3, &
     &        " frac")') ireg
            endif
          else
            write(iunit2,'(10x,"initial conditions: region",i3)') ireg
          endif
          write(iunit2,'(36("=+"))')
        endif
      
        iflag = 0
        iireg = ireg

!-------speciate aqueous solution
        call trspeciate (cloc,cxloc,cloctot,pgasloc, &
        cecloc,xexloc,csorplc,ccsorplc,csorpflc,sitdnloc,alogpf, &
        gamloc,gamxloc,dgamdi,rho,tempini,iter)

        if (myrank ==0) then
          if (ndloc1.eq.2) then
            if (mod(ireg,2).ne.0) then
              write(*,'(5x,"--> compute initial conditions: ", &
     &        "region:",i2," matrx  iter =",i4)') ireg,iter
            else
              write(*,'(5x,"--> compute initial conditions: ", &
     &        "region:",i2," frac   iter =",i4)') ireg,iter
            endif
          else if (idcdm.eq.1) then
            if (ireg.eq.2) then
              write(*,'(5x,"--> compute initial conditions: ", &
     &        "region:",i2," matrx  iter =",i4)') ireg,iter
            else
              write(*,'(5x,"--> compute initial conditions: ", &
     &        "region:",i2," frac   iter =",i4)') ireg,iter
            endif
          else
            write(*,'(5x,"--> compute initial conditions: ", &
     &      "region:",i2,"  iter =",i4)') ireg,iter
          endif
        endif
        
        do j = 1, ncomp
          guess(j,ireg) = cloc(j)
          c0(j)         = cloc(j)
          bc0(j)        = cloc(j)
          psi0(j)       = cloctot(j)*rho0
!         if (itype(j,ireg).ne.1) ctot(j,ireg) = cloctot(j)
        enddo

        if (iphase.eq.2 .or. iphase.eq.0) then
             call duanco2(tempini,pgasloc(1) ,dco2,fc,phi)
         
		  u1= (dco2/fmwco2*1D3)
 	  
		  do j = 1, ncomp
            psig0(j) = zero
            do l = 1, ngas
              psig0(j) = psig0(j)+ sgas(j,l)
            enddo
            psig0(j) = psig0(j) *u1!/(rgasjj*(tempini+tkelvin))
       !    print *,"ptran_set_ini: ", l,  tempini, pgasloc(1), psig0(j),u1
		  enddo
          do i = 1, ngas
            pgas0(i) =  pgasloc(1) !*u1
          enddo
        endif

!-------store initial conditions
        do n = 1, nlmax
          
!---------only store if in correct region
          nnreg = nreg_p(n)
          if(nnreg.eq.ireg) then ! controls single and dual contiuum

!           write(*,*) 'ptranspec: ',myrank,n,ireg,nreg_p(n)

!-----------ionic strength
!           sionic(n) = sionic0

!-----------primary species (components)
            do j = 1, ncomp
              kk = j+(n-1)*ncomp
              c_p(kk)   = c0(j)
!             cc(kk)    = c0(j) 
!             cprev(kk) = c0(j)

!             write(*,*) 'trstrtup: ',n,kk,c0(j),nreg_p(n),ireg

!             call VecSetValue(c,kk-1,c0(j),INSERT_VALUES,ierr)
!             call VecSetValue(cc,kk-1,c0(j),INSERT_VALUES,ierr)
!             call VecSetValue(cprev,kk-1,c0(j),INSERT_VALUES,ierr)
              
!             psi(j,n)  = cloctot(j)*rho(n)
!             gam(j,n)  = gamloc(j)
            enddo
            
!           print *,'speciate: ',myrank,n,ireg,(c0(j),j=1,ncomp)
            
!           call VecSetValuesBlocked(c,1,n-1,c0,INSERT_VALUES,ierr)

!-----------activity of H2O
            if (iact.eq.6) then
              ah2o(n) = ah2o0
            endif

!-----------total gaseous concentrations
            if (iphase.eq.2 .or. iphase.eq.0) then
             
			   do j = 1, ncomp
                ppsig(j,n) = psig0(j)
		!		print *, "ptran-spec::", j,n, ppsig(j,n)
              enddo
            endif

!-----------aqueous complexes
!           do i = 1, ncmplx
!             cx(i,n) = cxloc(i)
!             gamx(i,n) = gamxloc(i)
!           enddo

!-----------gaseous species
!           do i = 1, ngas
!             pgas(i,n) = pgas0(i)
!           enddo

!-----------ion exchange
            do m = 1, nexsolid
              do k = nsitex1(m), nsitex2(m)
                cec(k,n) = cecloc(k)
                do jj = nex1(k), nex2(k)
                  xexini(jj) = xexloc(jj)
                  xex(jj,n)  = xexloc(jj)
                  xxex(jj,n) = xexloc(jj)
                enddo
              enddo
            enddo

!-----------surface complexation
            do m = 1, nsrfmin
!             alogpf(m,n) = alogpf(m,1)
              do l = nsite1(m), nsite2(m)
                siteden(l,n) = sitdnloc(l)
                csorpf(l,n)  = csorpflc(l)
                do i = nsorp1(l), nsorp2(l)
                  csorpini(i) = ccsorplc(i)
                  csorp(i,n)  = ccsorplc(i)
                  ccsorp(i,n) = ccsorplc(i)
                enddo
              enddo
            enddo
          endif
        enddo ! end n-loop
        
!       call VecAssemblyBegin(c,ierr)
!       call VecAssemblyEnd(c,ierr)
!       call VecCopy(c,cc,ierr)
        
        if (nruns0>1 .and. myrank==0) then
          write(*,*) 'PROFile Multiple Runs: ',irun
        endif

        enddo
      enddo

!     if (modetr.eq.0) then
!       if (myrank == 0) then
!         write(*,*) &
!         '   --> distribution of species calculation completed!'
!         write(*,'(4x,50("="))')
!         write(*,*)
!         write(iunit2,*) &
!         '   --> distribution of species calculation completed!'
!         write(iunit2,'(4x,50("="))')
!       endif
!       stop
!     endif

      goto 200

  100 continue

!-----nonisothermal system

!     do ireg = ibcreg(1), ibcreg(2)
      do ireg = 1, initreg

        iflgcon = 0

!-------PROFile mulitple runs
        do irun = 1, nruns0
      
          if (nruns0 .gt. 1) then
            do j = 1, ncomp
              profji = profile(j,irun)
              if (profji.gt.zero) ctot(j,ireg) = profji
            enddo
            tempini = profile(ncomp+1,irun)
            tk = tempini+tkelvin
          endif

        do n = 1,nlmax               

          nnreg = nreg_p(n)
          if(nnreg.eq.ireg) then 

            iflgcon = iflgcon+1

!           tk = temp0 + tkelvin
            tc = temp_p(n)
            tk = tc + tkelvin
            
!           write(*,*) 'trstrtup: ',irun,n,nlmax,temp(n),
!    .      profile(ncomp+1,irun)

!-----------set gas constraints
            if (mode .eq. 3) then
!             u1 = press(n)/press(1)
              u1 = one
              do j = 1, ncomp
                if (itype(j,ireg).eq.4) then
                  ctot(j,ireg) = gmolfrac(j,ireg)*u1
                endif
              enddo
            endif

            if (n.gt.1 .and. iflgcon.gt.1) then
              do j = 1, ncomp
!               guess(j,ireg) = cloc(j)
                guess(j,ireg) = zero
             enddo
            endif

            iflag = 0
            iireg = ireg
            call trspeciate (cloc,cxloc,cloctot,pgasloc, &
            cecloc,xexloc,csorplc,ccsorplc,csorpflc,sitdnloc,alogpf, &
            gamloc,gamxloc,dgamdi,rho,tempini,iter)
            
            if (myrank == 0) then
              write(*,'(8x,"--> node: ",i6," region: ",i3," T = ", &
     &        1pg10.3," D = ",1pg10.3," iter = ",i4)') n,ireg, &
     &        tc,rho(n),iter

              write(iunit2,'(8x,"--> node: ",i6," region: ",i3," T = ", &
     &        1pg10.3," D = ",1pg10.3," iter = ",i4)') n,ireg, &
     &        tc,rho(n),iter
            endif
            
!-----------store initial conditions

!           sionic(n) = sionic0

!-----------activity of H2O
            if (iact.eq.6) then
              ah2o(n) = ah2o0
            endif


        if (iphase.eq.2 .or. iphase.eq.0) then
             call duanco2(tempini,pgasloc(1) ,dco2,fc,phi)
         
		  u1= (dco2/fmwco2*1D3) 
 	  
		  do j = 1, ncomp
            psig0(j) = zero
            do l = 1, ngas
              psig0(j) = psig0(j)+ sgas(j,l)
            enddo
            psig0(j) = psig0(j) *u1!/(rgasjj*(tempini+tkelvin))
   !        print *,"ptran_set_ini: ", l,  tempini,pgasloc(1) , psig0(j),u1
		  enddo
          do i = 1, ngas
            pgas0(i) = pgasloc(1) !*u1
          enddo
        endif

    !        if (iphase.eq.2 .or. iphase.eq.0) then
    !          do j = 1, ncomp
    !            psig0(j) = zero
    !            do l = 1, ngas
    !              psig0(j) = psig0(j)+ sgas(j,l)*pgasloc(l)
    !            enddo
    !            psig0(j) = psig0(j)/(rgasjj*tk)
    !          enddo
    !          do i = 1, ngas
    !            pgas0(i) = pgasloc(i)
    !          enddo
    !        endif

            do j = 1, ncomp
              kk = j+(n-1)*ncomp
              guess(j,ireg) = cloc(j)
              c0(j)     = cloc(j)
              bc0(j)    = cloc(j)
              psi0(j)   = cloctot(j)*rho0
              c_p(kk)   = c0(j)
!             cc(kk)    = c0(j) 
!             cprev(kk) = c0(j)
              
!             call VecSetValue(c,kk-1,c0(j),INSERT_VALUES,ierr)
!             call VecSetValue(cc,kk-1,c0(j),INSERT_VALUES,ierr)
!             call VecSetValue(cprev,kk-1,c0(j),INSERT_VALUES,ierr)

!             psi(j,n)  = cloctot(j)*rho(n)
              if(iphase.eq.2 .or. iphase.eq.0) ppsig(j,n) = psig0(j)
			  print *,n,j,ppsig(j,n)
!             gam(j,n)  = gamloc(j)
!             if (itype(j,ireg).ne.1 .and. n.eq.1) ctot(j,ireg) = &
!             cloctot(j)
            enddo
!           call VecSetValuesBlocked(c,1,n-1,c0,INSERT_VALUES,ierr)

!           do i = 1, ncmplx
!             cx(i,n)   = cxloc(i)
!             gamx(i,n) = gamxloc(i)
!           enddo
            do i = 1, ngas
              pgas0(i)  = pgasloc(i)
              pgas(i,n) = pgasloc(i)
            enddo
            do m = 1, nexsolid
              do k = nsitex1(m), nsitex2(m)
                cec(k,n) = cecloc(k)
                do jj = nex1(k), nex2(k)
                  xexini(jj) = xexloc(jj)
                  xex(jj,n)  = xexloc(jj)
                  xxex(jj,n) = xexloc(jj)
                enddo
              enddo
            enddo
            do m = 1, nsrfmin
!             alogpf(m,n) = alogpf(m,1)
              do l = nsite1(m), nsite2(m)
                siteden(l,n) = sitdnloc(l)
                csorpf(l,n)  = csorpflc(l)
                do i = nsorp1(l), nsorp2(l)
                  csorp(i,n)  = ccsorplc(i)
                  ccsorp(i,n) = ccsorplc(i)
                enddo
              enddo
            enddo
          endif
        enddo ! end n loop
        
!       call VecAssemblyBegin(c,ierr)
!       call VecAssemblyEnd(c,ierr)
!       call VecCopy(c,cc,ierr)
        
          if (nruns0>1 .and. myrank==0) then
            write(*,*) 'PROFile Multiple Runs: ',irun
          endif
        enddo ! end PROFile loop
      enddo ! end region loop
      
 200  continue
        
        
!     call VecAssemblyBegin(c,ierr)
!     call VecAssemblyEnd(c,ierr)

      call VecRestoreArrayF90(c,c_p,ierr)
      call VecCopy(c,cc,ierr)

!     call VecView(c,PETSC_VIEWER_STDOUT_WORLD,ierr)
!     call VecView(cc,PETSC_VIEWER_STDOUT_WORLD,ierr)

      call VecRestoreArrayF90(temp,temp_p,ierr)
      call VecRestoreArrayF90(nreg,nreg_p,ierr)
      deallocate(nreg_val)
      
      return
      end subroutine trstartup

!=========================================================================

      subroutine trspeciate (cc,cx,psi,pgas, &
           cec,xex0,csorp,ccsorp,csorpf,siteden,alogpf, &
           gam,gamx,dgamdi,rho,temp,iter)

      use ptran_global_module
      use ptran_dbase_module
            
      implicit none

      real*8 :: cx(*),gam(*),gamx(*),psi(*),cc(*),pgas(*), &
                cec(*),xex0(nexmx),csorp(*),ccsorp(*),csorpf(*),siteden(*)
      real*8 :: rho(*),temp
      real*8 :: alogpf(:,:)

      real*8 :: acratio(ncmx),ret(ncmx),retc(ncmx), &
                tsorb(ncmx),tsorbc(ncmx),phik_reg(nkmx,nrgmx)

      real*8 :: actr,ah2o0,ak,alno2,alogfo2,aph,av,blkden, &
                cjtot,dco2,den,dgamdi(*), &
                distj,dwat,d1m,eh,ehfac, &
                xco2,x1m,xpsi, &
                faraday0,fco2,grainden,tc,fvcoll, &
                ph,pe,percent,perctads,pfac,phico2,phih2o, &
                retard, &
                sj,skdc,skdm,sum,sumpt,sumx,srfpot,surfchrg,skdmj, &
                skdcj,tmpk,tj,qz2,xh2o,xh2oapp,ah2oloc

      integer :: i,ii,j,jj,k,l,m,irun,iter,kk,mm

      integer :: indx(ncxmx)

      character*16 nmisthrm(3)
      character*30 atyp
!     character*20 maspltaq

      data nmisthrm/'Gaines-Thomas','Vanselow','Gapon'/

      data irun/1/

      faraday0 = 23062.3d0
      ehfac = (rgas*tk)/faraday0

!-----compute initial cec and initialize kd
      do m = 1, nexsolid
        mm = mex(m) ! mineral or colloid as primary aqueous species
        do k = nsitex1(m), nsitex2(m)
          if (m.le.nexsolid-ncollex) then
            blkden = phik_reg(mm,ireg)*wtkin(mm)/vbarkin(mm)
            cec(k) = cec0(ireg,k)*blkden
            do jj = nex1(k), nex2(k)
              xex0(jj) = zero
            enddo
            xex0(nex1(k)) = one
          else
            cec(k) = cec0(ireg,k)*wt(mm)*ctot(mm,ireg)
            do jj = nex1(k), nex2(k)
              xex0(jj) = zero
            enddo
            xex0(nex1(k)) = one
          endif

!         write(*,*) 'trspcate-ex: ',k,cec(k)
        enddo
      enddo

!-----compute surface complexation properties
      do mm = 1, nsrfmin
        m = msorp(mm)

        do l = nsite1(mm), nsite2(mm)

!---------set initial guess for double layer potential calculation
!         if (idblpot.eq.0) then
!           alogpf(mm,ireg) = zero
!         else
!           alogpf(mm,ireg) = 6.3748
!         endif

          if (mm.le.nsrfmin-ncolsrf) then
!-----------minerals: Ns/V = Ns/Am . Am/Mm . Mm/Nm . Nm/Vm . Vm/V
            siteden(l) = siteden0(ireg,l)*areamass(mm)*wtkin(m)/ &
            vbarkin(m)*phik_reg(m,ireg)
          else
!-----------colloids: Ns/V = Ns/Ac . Ac/Mc . Mc/Nc . Nc/Vp
            siteden(l) = siteden0(ireg,l)*areamass(mm)*wt(m)* &
            ctot(m,ireg)
          endif

!         write(*,*) 'trspcate-srf: ',l,siteden(l)
        enddo
      enddo

!-----compute distribution of species      
      call treqlib (cx,gamx,gam,psi,cc,pgas,xex0,cec, &
      siteden,csorpf,csorp,ccsorp,alogpf,dgamdi,iter)

      if(iact.eq.0) then
        sum = zero
        do j = 1, ncomp
          if (nam(j) .ne. 'e-') sum = sum+z(j)**2*cc(j)
        enddo
        do i = 1, ncmplx
          sum = sum+zx(i)**2*cx(i)
        enddo
        sionic0 = half*sum
      else if (iact.eq.6) then
        sum = zero
        do j = 1, nmass
          sum = sum + cc(j)
        enddo
        do i = 1, ncmplx
          sum = sum + cx(i)
        enddo
        xh2o = one/(one+wh2o*sum)
        xh2oapp = one-wh2o*sum
!pcl    ah2o = exp(phih2o*log(xh2o))
        ah2oloc = exp(-phih2o*wh2o*sum) ! definition of phih2o
        ah2o0 = -phih2o*wh2o*sum
      endif

!     if (isothrm.ge.1 .and. ireg.le.ibcreg(2)) return
!     if (ireg.le.ibcreg(2)) return

!-----printout results for isothermal case only

!-----order complexes in decending concentration
      indx(1) = 1
      if (ncmplx .gt. 1) call indexx(ncmplx,cx,indx)
      
      if (myrank == 0) &
      write(iunit2,'(" temperature = ",1pe12.4," [C]", &
     &" pressure = ",1pe12.4," [bars]")') temp,pref0*1.e-5
!    &" pressure = ",1pe12.4," [bars]")') tempini,pref0*1.e-5

      if (myrank == 0) write(iunit2,1051) iter,sionic0

!-----compute solution density
      if (jh2o.eq.0) then
        den = rho0
      else
        den = zero
      endif
      do j = 1, ncomp
        den = den + wt(j)*psi(j)
      enddo
      if (myrank == 0) &
      write(iunit2,'(" solution density    =",1pg12.5," [g/cm^3]" &
     &)') den

      if (iact==6 .and. myrank==0) then
        write(iunit2,'( &
     &  " mole fraction H2O   = ",1pe12.4, &
     &  " (approx. ",1pe10.4,")"/ &
     &  " activity of water   = ",1pe12.4, &
     &  " osmotic coefficient = ",1pe12.4)') xh2o,xh2oapp, &
     &  ah2oloc,phih2o
      endif

!-----compute pH
      ph = zero
      if (jph.gt.0) then
        aph = cc(jph)*gam(jph)
        if (jh2o.gt.0) aph = aph/(wh2o*cc(jh2o))
        ph  = -log10(aph)
      else if (iph.gt.0) then
        aph = cx(iph)*gamx(iph)
        if (jh2o.gt.0) aph = aph/(wh2o*cc(jh2o))
        ph  = -log10(aph)
      endif
      if (jph.gt.0 .and. jo2.gt.0 .or. iph.gt.0 .and. jo2.gt.0) then
        if (isothrm .ne. 0) then
!pcl      alnkeh = flogk(coefeh,tempini,tkelvin)
          alnkeh = flogk(coefeh,tk)
        endif

!-------compute fo2 fugacity
        alogfo2 = zero
        alno2 = eqgas(jo2)*aln10
        do kk = 1, ncomp
          ak = cc(kk)*gam(kk)
          if (ak.gt.zero) alno2=alno2+log(ak)*sgas(kk,jo2)
        enddo
        alogfo2 = alno2/aln10
        eh = ehfac*(four*log(aph)+(alogfo2-alnkeh)*aln10)/four
        pe = eh/(ehfac*aln10)
        
!       print *,'ptranspec: ',jo2,alno2,alogfo2,ehfac,alnkeh,aph,pe,eh,four
        
        if (myrank==0) write(iunit2,1090) ph,pe,eh
      else if (jph.gt.0 .or. iph.gt.0) then
        if (myrank==0) write(iunit2,1091) ph
      endif

    if (myrank == 0) then
        
      ! mole fraction
      sumx = 0.d0
      do j = 1, ncomp
        sumx = sumx + cc(j)
      enddo
      do i = 1, ncmplx
        sumx = sumx + cx(i)
      enddo
      
      
      
      write(iunit2,'()') 
      
      write(iunit2,'(72("-"))')
      write(iunit2,'(''species          xmol        Xmol'')')
      write(iunit2,'(72("-"))')
      do j = 1, ncomp
        xpsi = cc(j)
        do i = 1, ncmplx
          xpsi = xpsi + shom(j,i)*cx(i)
        enddo
        if(j/=jco2)&
        write(iunit2,'(a12,1p10e12.4)') nam(j),cc(j)/sumx,xpsi/sumx
       if(j==jco2)&
        write(iunit2,'(a12,1p10e12.4)') nam(j),cc(j)/sumx,xpsi/sumx,&
              10.d0**eqgas(jco2g)*xpsi/(pref0*1D-5)
      enddo
        
      if (jph > 0) then
        write(iunit2,'(72("-"))')
        write(iunit2,1010)
        write(iunit2,'(72("-"))')
        do j = 1, ncomp
          if(cc(jph).gt.0.) acratio(j) = log10(cc(j)) - z(j)*log10(cc(jph)) 
          write(iunit2,1060) nam(j),cc(j),psi(j),gam(j),acratio(j), &
          itype(j,ireg),ncon(j,ireg)
        enddo
        if(ncmplx.ne.0) then
          write(iunit2,1025)
          write(iunit2,'(72("-"))')
          do i = 1, ncmplx
            ii = indx(ncmplx-i+1)
            if(cc(jph).gt.zero .and. cx(ii).gt.zero) &
            actr = log10(cx(ii))-zx(ii)*log10(cc(jph))
            write(iunit2,1040) namcx(ii),cx(ii),gamx(ii),actr,eqhom(ii)
          enddo
        endif
      else
        write(iunit2,'(72("-"))')
        write(iunit2,1015)
        write(iunit2,'(72("-"))')
        do j=1,ncomp
          write(iunit2,1061) nam(j),cc(j),psi(j),gam(j), &
          itype(j,ireg),ncon(j,ireg)
        enddo
        if(ncmplx.ne.0) then
          write(iunit2,1020)
          write(iunit2,'(72("-"))')
          do i=1,ncmplx
            ii = indx(ncmplx-i+1)
            write(iunit2,1040) namcx(ii),cx(ii),gamx(ii),eqhom(ii)
          enddo
        endif
      endif
    endif
    
!      if (myrank==0)&
	   call trsolprd (gam,cc,psi,pgas,ah2o0)
      
      if (iflgco2>0 .and. myrank==0) then
        write(iunit2,'(/," fugacity CO2 = ",1pe12.4," [bars]", &
     &  " fugacity coef. = ",1pe12.4,/, &
     &  " pure CO2 phase density = ",1pe12.4," [kg/m^3]",/, &
     &  " pure H2O phase density = ",1pe12.4," [kg/m^3]",/, &
     &  " mixture density = ",1pe12.4," [kg/m^3]",/, &
     &  " mole fraction CO2(aq) = ",1pe12.4,/, &
     &  " mixture mass fraction CO2(aq) = ",1pe12.4)') &
     &  fco2,phico2,dco2,dwat,d1m,xco2,x1m
      endif

!*********************************************************************** 
!        check electroneutality of initial concentration and flux
!*********************************************************************** 
!     qz1=zero 
!     do j=1,ncomp
!       if (nam(j) .ne. 'e-') qz1=qz1+z(j)*psi(j)
!       write(*,*) 'speciate: ',nam(j),j,qz1,z(j),psi(j)
!     enddo
      qz2=zero 
      do j=1,ncomp
        if (nam(j) .ne. 'e-') qz2=qz2+z(j)*cc(j)
      enddo
      do i=1,ncmplx
        if (nam(j) .ne. 'e-') qz2=qz2+zx(i)*cx(i)
      enddo
      if (myrank==0) write(iunit2,1105) qz2 

!*********************************************************************** 
!-----print out surface complexation properties
!*********************************************************************** 
      if (iprint>=0 .and. nsrfmin>0 .and. myrank==0) then
        write(iunit2,*)
        write(iunit2,*) 'surface complexation properties'
        do m = 1, nsrfmin
          tmpk = temp0
!         srfpot = rgasj*tmpk*alogpf(m,1)/faraday
          srfpot = 0.d0
!         pfac = exp(-alogpf(m,1))

          av = areamass(m)/coverage(m)*1.e-5

          write(iunit2,'(72("-"))')
          if (m .le. nsrfmin-ncolsrf) then            ! minerals
            write(iunit2,'(" mineral: ",a20," area:",1pe10.3, &
     &      " [m^2/g]"," cover:",1pe10.3)') namsrf(m),av,coverage(m)
          else
            write(iunit2,'(" colloid: ",a20," area:",1pe10.3, &
     &      " [m^2/g]"," cover:",1pe10.3)') namsrf(m),av,coverage(m)
          endif
          write(iunit2,'(72("-"))')

          if (srfpot.ne.zero) then
            write(iunit2,'(" pot: ",1pe12.4," [V]"," P: ",1pe12.4, &
     &      /,9x,"  surface charge: ",1pe12.4," [Coul/dm^2]")') &
            srfpot,pfac,surfchrg
            write(iunit2,'(72("-"))')
          endif

          percent = zero
          do l = nsite1(m), nsite2(m)
            sum = csorpf(l)
            do i = nsorp1(l), nsorp2(l)
              sum = sum + ccsorp(i)
            enddo
            if (sum.gt.zero) percent = csorpf(l)/sum
            sumpt = percent
            write(iunit2,'(" site density: ",1pe12.4, &
     &      " [mol/dm^3]")') siteden(l)
            write(iunit2,*) 'species              valence  ', &
     &      '  conc. [mol/dm^3] fraction   log K'
            write(iunit2,'(72("-"))')
            write(iunit2,'(1x,a20,3(1pe12.4,1x))') namsite(l),zsite(l), &
            csorpf(l),percent
            do i = nsorp1(l), nsorp2(l)
              if (sum.gt.zero) percent = ccsorp(i)/sum
              sumpt = sumpt + percent
              write(iunit2,'(1x,a20,4(1pe12.4,1x))') namscx(i),zsrf(i), &
              ccsorp(i),percent,eqsorp(i)
            enddo
            write(iunit2,'(35x,"      total =",1pe11.4,/)') sumpt
          enddo
        enddo

!-------retardation factor
        write(iunit2,'(72("-"))')
        write(iunit2,*) '              distribution coefficient'  
        write(iunit2,'(72("-"))')
        write(iunit2,*) 'component      mineral    colloid', &
        '  retardation  cjtot    % adsorbed'
        write(iunit2,'(6x,"site",4x,"[dm^3/dm^3]",3x,"[L/L]",2x, &
     &  "   [1+Kd]    [mol/L]")')
        write(iunit2,'(72("-"))')
        do j = 1, ncomp
          write(iunit2,'(1x,a12,1p5e12.4)') nam(j)
          skdm = zero
          skdc = zero
          do m = 1, nsrfmin
            if (m .le. nsrfmin-ncolsrf) then            ! minerals
              do k = nsite1(m), nsite2(m)
                do i = nsorp1(k), nsorp2(k)
                  skdm = skdm + ssorp(j,i)*ccsorp(i)
                enddo
              enddo
            else                                        ! colloids
              do k = nsite1(m), nsite2(m)
                do i = nsorp1(k), nsorp2(k)
                  skdc = skdc + ssorp(j,i)*ccsorp(i)
                enddo
              enddo
            endif
            do k = nsite1(m), nsite2(m)
              if (m .le. nsrfmin-ncolsrf) then            ! minerals
                skdmj = zero
                do i = nsorp1(k), nsorp2(k)
                  skdmj = skdmj + ssorp(j,i)*ccsorp(i)
                enddo
                skdmj = skdmj/psi(j)/por_reg(ireg)
                write(iunit2,'(5x,a8,1pe11.4)') namsite(k),skdmj
              else
                skdcj = zero
                do i = nsorp1(k), nsorp2(k)
                  skdcj = skdcj + ssorp(j,i)*ccsorp(i)
                enddo
                skdcj = skdcj/psi(j)
                write(iunit2,'(5x,a8,12x,1pe11.4)') namsite(k),skdcj
              endif
            enddo
          enddo
          cjtot = psi(j) + skdc + skdm/por_reg(ireg)
          skdm = skdm/psi(j)/por_reg(ireg)
          skdc = skdc/psi(j)
!---------note: porosity should correspond to region 
!               (incorrect as stands)
          retard = one+skdm/(one+skdc)
          perctads = skdm/retard*1.d2
          write(iunit2,'(13x,1p6e11.4)') skdm,skdc,retard, &
          cjtot,perctads

!         write(*,*) 'trspcate: ',nam(j),ireg,por_reg(ireg),psi(j)

          write(iunit2,'(72("-"))')
        enddo
        write(iunit2,'()')
      endif

!*********************************************************************** 
!-----print out ion-exchange properties
!*********************************************************************** 
      if (nexmax>0 .and. myrank==0) then
        do j = 1, ncomp
          ret(j) = zero
          retc(j) = zero
          tsorb(j) = zero
          tsorbc(j) = zero
        enddo

        write(iunit2,'()')
        write(iunit2,'(10x,"ion-exchange properties: ",a16, &
     &  " isotherm")') nmisthrm(ionex+1)
        do m = 1, nexsolid-ncollex
          mm = mex(m)
          blkden = phik_reg(mm,ireg)*wtkin(mm)/vbarkin(mm)
          grainden = wtkin(mm)/vbarkin(mm)
          write(iunit2,1081) namex(m),phik_reg(mm,ireg),blkden,grainden
          write(iunit2,'(72("-"))')
          do k = nsitex1(m), nsitex2(m)
            write(iunit2,'(" site: ",i3," cec = ",1pe12.4, &
     &      " [mol/kg], ",1pe12.4," [mol/dm^3]")') &
            k,cec0(ireg,k),cec(k)
            write(iunit2,1082)
            write(iunit2,'(72("-"))')
            do jj = nex1(k), nex2(k)
              j = jex(jj)
              sj = cec(k)*xex0(jj)/z(j)
              distj = sj/psi(j)/por_reg(ireg)
              ret(j) = ret(j) + distj
              tj = psi(j)+sj/por_reg(ireg)
              tsorb(j) = tsorb(j) + sj
              write(iunit2,1092) nam(j),z(j),eqiex(jj),xex0(jj),sj, &
              distj,tj
            enddo
            write(iunit2,'()')
          enddo
        enddo
        do m = nexsolid-ncollex+1, nexsolid
          write(iunit2,1083) namex(m),cc(mex(m))
          write(iunit2,'(72("-"))')
          do k = nsitex1(m), nsitex2(m)
            write(iunit2,'(" site: ",i3," cec = ",1pe12.4, &
     &      " [mol/kg], ",1pe12.4," [mol/dm^3]")') & 
            k,cec0(ireg,k),cec(k)
            write(iunit2,1082)
            write(iunit2,'(72("-"))')
            do jj = nex1(k), nex2(k)
              j = jex(jj)
              sj = cec(k)*xex0(jj)/z(j)
              distj = sj/psi(j)
              retc(j) = retc(j) + distj
              tj = psi(j)+sj
              tsorbc(j) = tsorbc(j) + sj
              write(iunit2,1092) nam(j),z(j),eqiex(jj),xex0(jj),sj, &
              distj,tj
            enddo
            write(iunit2,'()')
          enddo
        enddo

!-------total retardation
        write(iunit2,'(72("-"))')
        write(iunit2,*) '          distribution coef.', &
        '   retardation '
        write(iunit2,*) 'species    Kd_r       Kd_c         R', &
        '        tot sorb   totc sorb'
        write(iunit2,*) '               [L/dm^3] ', &
        '         [1+Kd]'
        write(iunit2,'(72("-"))')
        m = 1
        k = nsitex1(m)
        do j = 1, ncomp
          do jj = nex1(k), nex2(k)
            if (namcat(jj) .eq. nam(j)) then
!             retard = 1+ret(j)/(1+retc(j))
              retard = (1+retc(j)+ret(j))/(1+fvcoll*retc(j))
              write(iunit2,'(1x,a8,1p10e11.4)') nam(j),ret(j),retc(j), &
              retard,tsorb(j),tsorbc(j)
            endif
          enddo
        enddo
        write(iunit2,'(72("-"))')
        write(iunit2,'()')
      endif
 1081 format(' mineral = ',a20,' vol frac =',1pe12.4,/, &
      ' bulk density = ',1pe12.4,' [g/cm^3]', &
      ' grain density = ',1pe12.4,' [g/cm^3]') 
 1082 format(' ion     valence  eqiex       chi         s       ', &
      'Kd    ctot(aq+sorbed)')
 1083 format(' mineral = ',a20,' conc.    =',1pe12.4) 
 1092 format(1x,a12,0pf2.0,1p5e11.4)

 1105 format(/' charge balance - q = ',1p2e12.4)
 1051 format(' iter =',i3,5x,' ionic strength =',1pe12.4)
 1060 format(1x,a12,1pe11.4,1pe12.4,1pe11.4,1pe12.4,1x,i2,1x,a12)
 1061 format(1x,a12,1p3e12.4,1x,i2,1x,a12)
 1070 format(' ',1p10e12.4)
 1072 format(' ion-exchange solid composition',/, &
             ' cation ',5x,'cjbar ',8x,'kd')
 1073 format(' ',a8,1p2e12.4)
 1080 format(/' ',10a12) 
 1090 format(' pH =',1pg12.5,' pe =',1pg12.5,' Eh =',1pe12.4,' [V]')
 1091 format(' pH =',1pg12.5)
 1010 format(' species      molality      tot    act. coef. ', &
             'act. ratio/H+ constraint')
 1015 format(' species       molality       tot     act. coef. ', &
             ' constraint')
 1040 format(1x,a20,1pg12.5,1x,1p5g12.5)
 1020 format(/' complex               molality   act. coef.     log K')
 1025 format(/' complex               molality   act. coef.    act/H+ ', &
              '      log K')
 1030 format(' ',65x,10a8) 

!-----cut and paste initial and boundary conditions
      if (isothrm .eq. 10 .and. myrank==0) then
        write(iunit2,'(/,60("="))')
        write(iunit2,'("cut and paste initial and boundary solution" &
     &  )')
        write(iunit2,'(60("="))')
        write(iunit2,'(":species     itype   guess       ctot", &
     &               "         constraint")')
        do j = 1, ncomp
          write(iunit2,'(a12,1x,i3,3x,1p2e12.4,3x,a20)') nam(j), &
          itype(j,ireg),cc(j),ctot(j,ireg),ncon(j,ireg)
        enddo
        write(iunit2,'(60("="))')
      endif

!-----write out PROFile results to file: prof.xyp
      if (nruns>1 .and. modetr==0 .and. myrank==0) then
        if (irun.eq.1) then
        
!         open(101,file='prof.xyp',status='unknown')
          atyp = 'pf'
!         call openum(101,atyp,maspltaq)
          
          if (isothrm.eq.0) then
            if (nexsolid .gt. 0) then
              write(101,'(a6,(100a12))') '#PROF ', &
              (nam(j),j=1,ncomp), &
              'pH          ','logfO2      ',(namg(l),l=1,ngas), &
              'xco2        ','fco2        ','phico2      ', &
              'denco2      ',(nam(j),j=1,ncomp)
            else
              write(101,'(a6,(100a12))') '#PROF ', &
              (nam(j),j=1,ncomp), &
              'pH          ','logfO2      ',(namg(l),l=1,ngas), &
              'xco2        ','fco2        ','phico2      ', &
              'denco2      '
            endif
          else
            if (nexsolid .gt. 0) then
              write(101,'(a6,(100a12))') '#PROF ','T ', &
              (nam(j),j=1,ncomp), &
              'pH          ','logfO2      ',(namg(l),l=1,ngas), &
              'xco2        ','fco2        ','phico2      ', &
              'denco2      ',(nam(j),j=1,ncomp)
            else
              write(101,'(a6,(100a12))') '#PROF ','T ', &
              (nam(j),j=1,ncomp), &
              'pH          ','logfO2      ',(namg(l),l=1,ngas), &
              'xco2        ','fco2        ','phico2      ', &
              'denco2      '
            endif
          endif
          irun = 0
        endif
        if (isothrm.eq.0) then
          if (nexsolid .gt. 0) then
            write(101,'(1p100e12.4)') (psi(j),j=1,ncomp),ph, &
            alogfo2,(pgas(l),l=1,ngas),xco2,fco2,phico2,dco2, &
            (tsorb(j),j=1,ncomp)
          else
            write(101,'(1p100e12.4)') (psi(j),j=1,ncomp),ph, &
            alogfo2,(pgas(l),l=1,ngas),xco2,fco2,phico2,dco2
          endif
        else
          tc = tk-tkelvin
          if (nexsolid .gt. 0) then
            write(101,'(1p100e12.4)') tc,(psi(j),j=1,ncomp),ph, &
            alogfo2,(pgas(l),l=1,ngas),xco2,fco2,phico2,dco2, &
            (tsorb(j),j=1,ncomp)
          else
            write(101,'(1p100e12.4)') tc,(psi(j),j=1,ncomp),ph, &
            alogfo2,(pgas(l),l=1,ngas),xco2,fco2,phico2,dco2
          endif
        endif
      endif
      
      return
      end subroutine trspeciate

      subroutine indexx(nc,array,indx)

      implicit none
      
      real*8 :: array(*),q
      integer :: i,j,kk,l,nc,indx(*),indxt

      do j=1,nc
        indx(j)=j
      enddo  

      l=nc/2+1
      kk=nc

   10 if(l.gt.1) then
        l=l-1
        indxt=indx(l)
        q=array(indxt)
      else
        indxt=indx(kk)
        q=array(indxt)
        indx(kk)=indx(1)
        kk=kk-1
        if(kk.eq.1) then
          indx(1)=indxt
          goto 30
        endif
      endif
      i=l
      j=l+l
   20 if(j.le.kk) then
        if(j.lt.kk) then
          if(array(indx(j)).lt.array(indx(j+1))) j=j+1
        endif
        if(q.lt.array(indx(j))) then
          indx(i)=indx(j)
          i=j
          j=j+j
        else
          j=kk+1
        endif
        goto 20
      endif
      indx(i)=indxt
      goto 10
   30 continue

      return
      end subroutine indexx

!=========================================================================

      subroutine treqlib (cx,gamx,gam,psi,cc,pgas,xex0,cec, &
      siteden,csorpf,csorp,ccsorp,alogpf,dgamdi,iter)

      use ptran_global_module
            
      implicit none
      
      real*8 cx(*),gamx(*),gam(*),psi(*),cc(*),pgas(*),cec(*),xex0(*), &
      siteden(*),csorpf(*),csorp(*),ccsorp(*),alogpf(:,:)

      real*8 :: zz(ncmx,ncmx),y(ncmx),dgamdi(*)
      
      real*8 :: ah2o0,ctotsav,sumcec
      
      integer :: i,iter,j,k,jchrg,logsav,m

!-----------------------------------------------------------------------
!       constraints        description
!-----------------------------------------------------------------------
!       type  
!                1   t     total concentration
!                2  ex     ion-exchange constraint
!                3   m     mineral phase
!                4   g     gas phase
!                5   -     weight fraction of solvent
!                6   -     mole fraction of solvent
!                7   c     concentration (or activity in the case of ph)
!                8  ph     
!                9   -     activity
!               10   -     conserved quantity  
!               11   -     CO2(g) eos  
!               12   -     CO2(g) eos  
!               -1   e     electroneutrality
!               -2* eh     Eh contraint
!               -3* pe     pe contraint
!-----------------------------------------------------------------------
!                 * -not implemented
!-----------------------------------------------------------------------
!     This routine obtains the initial and boundary fluid compositions
!     through a distribution of species calculation. Several passes are
!     made: first a pre-newton-raphson iteration scheme is used to
!     approximate the initial guess for the newton-raphson scheme;
!     second if charge balance is invoked a solution is found without
!     charge balance and then charge balance is tested to determine if
!     this is possible for the chosen species (unless pH is chosen---
!     which can dramatically alter the distribution of charges!).
!-----------------------------------------------------------------------
   
      logsav = loglin
      loglin = 1
      iter   = 0
      jchrg  = 0
      ah2o0  = zero

!-----transform primary species
!     j0 = 1
!     i0 = 4
!     do j = 1, ncomp
!       psi(j) = ctot(j,ireg)
!     enddo
!     call trfpri(j0,i0,cc,cc,cx,psi,psi,gam,gamx)
      
      do j = 1, ncomp
!       write(*,'("treqlib: ",a12,1x,i3,1p10e12.4)') nam(j), &
!    &  itype(j,ireg),guess(j,ireg),cc(j),psi(j),gam(j),ctot(j,ireg)
        cc(j) = zero
        psi(j) = zero
        gam(j) = one
!       guess(j,ireg) = zero
      enddo

      do i = 1, ncmplx
        gamx(i) = one
      enddo

!-----store initial guess for aqueous concentration

!     write(*,*) 'treqlib: ireg= ',ireg,(guess(j,ireg),j=1,ncomp)
      
      do j=1,ncomp
        if (guess(j,ireg).eq.zero) then
          if      (itype(j,ireg).le.2) then
            cc(j) = abs(ctot(j,ireg))
          else if (itype(j,ireg).eq.3) then
            cc(j) = ctot(j,ireg)
          else if (itype(j,ireg).eq.4 .or. itype(j,ireg).eq.11 .or. &
            itype(j,ireg).eq.12 .or. itype(j,ireg).eq.7) then
            if (ctot(j,ireg) < 0.d0) then
              ctot(j,ireg) = 10.d0**ctot(j,ireg)
            endif
            cc(j) = ctot(j,ireg) ! check if correct initial value-pcl
          else if (itype(j,ireg).eq.5) then
            cc(j) = ctot(j,ireg)
          else if (itype(j,ireg).eq.6) then
            cc(j) = ctot(j,ireg)
          else if (itype(j,ireg).eq.8) then
            if (j.eq.jph) then
              cc(j) = 10.d0**(-ctot(j,ireg))/gam(j)
            else if (j.eq.joh) then
              cc(j) = 10.d0**(-eqhom(iph)-ctot(j,ireg))/gam(j)
            endif
          else if (itype(j,ireg).eq.10) then
            cc(j) = ctot(j,ireg)
          else if (itype(j,ireg).eq.-1) then
            cc(j) = ctot(j,ireg)
          endif
          guess(j,ireg)=cc(j)
        else
          if (itype(j,ireg).ne.7) then
            cc(j) = guess(j,ireg)
          else if (itype(j,ireg).eq.7) then
            cc(j) = ctot(j,ireg)
            guess(j,ireg)=cc(j)
          endif
        endif

!-------reset charge balance constraint on first pass
        if (itype(j,ireg).eq.-1) then
          jchrg = j
          if (j.ne.jph .or. j.ne.joh) then
!           itype(j,ireg) = 7
!pcl        itype(j,ireg) = 1
            cc(j) = ctot(j,ireg)
          else if (j.eq.jph) then
!           itype(j,ireg) = 8
!           cc(j) = 10.d0**(-ctot(j,ireg))/gam(j)
            cc(j) = ctot(j,ireg)
          else if (j.eq.joh) then
!           itype(j,ireg) = 8
!           cc(j) = 10.d0**(-eqhom(iph)-ctot(j,ireg))/gam(j)
            cc(j) = ctot(j,ireg)
          endif
          guess(j,ireg)=cc(j)
        endif
        
!-------reset ctot for sorption
        if (itype(j,ireg).eq.21) then
!         itype(j,ireg) = 2
          sumcec = zero
          do m = 1, nexsolid
            do k = nsitex1(m), nsitex2(m)
              sumcec = sumcec + cec(k)/z(j)
            enddo
          enddo
          ctotsav = ctot(j,ireg)
          ctot(j,ireg) = ctot(j,ireg) + sumcec
          cc(j) = ctot(j,ireg)
!         write(*,*) 'treqlib: ',nam(j),sumcec,ctotsav,ctot(j,ireg)
        endif
      enddo

      if (mnrl.gt.0 .or. ngas.gt.0) &
      call constrnt (y,zz,cc,gam,ah2o0)
      
      call newtonraphson (jchrg,cx,gamx,gam,psi,cc,pgas, &
      xex0,cec,siteden,csorpf,csorp,ccsorp,alogpf,dgamdi,iter)

      loglin = logsav

      return
      end subroutine treqlib
      
!========================================================================

      subroutine newtonraphson (jchrg,cx,gamx,gam,psi,cc,pgas, &
      xex0,cec,siteden,csorpf,csorp,ccsorp,alogpf,dgamdi,iter)
      
      use ptran_global_module
      use ptran_dbase_module

      implicit none
      
      integer :: i,j,k,l,istate,iter,itmax0,iactsav,jchrg,newtconv
      
      real*8 :: alpha,chrg,chrgj,ctchrg,cave,chold,ctotsav,dd,dur,fac,x, &
                tolj,resfac,zchrg,zchrg0

      real*8 :: cx(*),gamx(*),gam(*),psi(*),cc(*),pgas(*), &
                cec(*),xex0(*),siteden(*),csorpf(*),csorp(*),ccsorp(*), &
                alogpf(:,:),dgamdi(*)

      real*8 :: zz(ncmx,ncmx),y(ncmx),dpsi0(ncmx,ncmx)
      
      integer :: indx(ncmx)

      real*8 :: scalefac(ncmx)

!     itmax0 = 5000
      itmax0 = 500
      
!-----store initial guesses
      do j = 1, ncomp
        guess(j,ireg) = cc(j)
        atol(j) = one
      enddo

      iter=0

!-----set up jacobian matrix and residual for newton-raphson eqns.

!-----first obtain convergence without charge balance constraint
      istate = 1
      if (jchrg.gt.0) then
        istate = 0
        j = jchrg
        ctotsav = ctot(jchrg,ireg)
        if (j.ne.jph .or. j.ne.joh) then
!         itype(j,ireg) = 7
!pcl      itype(j,ireg) = 1
          cc(j) = ctot(j,ireg)
        else if (j .eq. jph) then
          itype(j,ireg) = 8
          cc(j) = 10.d0**(-ctot(j,ireg))/gam(j)
!         cc(j) = ctot(j,ireg)
        else if (j .eq. joh) then
          itype(j,ireg) = 8
          cc(j) = 10.d0**(-eqhom(iph)-ctot(j,ireg))/gam(j)
!         cc(j) = ctot(j,ireg)
        endif
      endif

!-----begin newton-raphson iteration without activity coefficient 
!       corrections
      iactsav = iact

!     if (iact.ne.6) iact = 0
      iact = 0

!-----check for negative concentrations
      do j = 1, ncomp
        if (cc(j) .lt. zero) then
          write(*,*) 'eqlib-neg0: ',iter,nam(j),cc(j)
          cc(j) = guess(j,ireg)
        endif
      enddo

   20 continue
      iter = iter+1

!     write(*,*) 'treqlib/newton-1: ',loglin,ireg
!     do j = 1, ncomp
!       write(*,'(i3,1x,a12,1x,i3,1p10e12.4)') iter,nam(j), &
!       itype(j,ireg),cc(j),psi(j),gam(j)
!     enddo

      if(iter.gt.itmax0) goto 333

      call treqres (cc,y,cx,psi,xex0,gam,gamx, &
      pgas,siteden,cec,csorpf,csorp,ccsorp,alogpf,dgamdi)

      call treqjac (cc,cx,psi,gam,dpsi0,pgas, &
      cec,xex0,siteden,ccsorp,zz)

!-----scale jacobian and residual
      if (iter .gt. 100) then
        do j = 1, ncomp
          scalefac(j) = one
          do l = 1, ncomp
            scalefac(j) = max(abs(zz(j,l)),scalefac(j))
          enddo
        enddo
        do j = 1, ncomp
          y(j) = y(j)/scalefac(j)
          do l = 1, ncomp
            zz(j,l) = zz(j,l)/scalefac(j)
          enddo
        enddo
      endif

!-----scale jacobian for logarithmetic concentrations
      if (loglin .eq. 0) then
        do j = 1, ncomp
          do l = 1, ncomp
            zz(j,l) = zz(j,l)*cc(l)
          enddo
        enddo
      endif
      
!       write(*,1010)
!       do j=1,ncomp
!         write(*,1000) iter,nam(j),y(j),cc(j),ctot(j,ireg), &
!         psi(j),gam(j)
!       enddo

      if(iprint.eq.3 .or. iter.eq.itmax0-1 .and. myrank==0) then
        write(iunit2,1010)
        do j=1,ncomp
          write(iunit2,1000) iter,nam(j),y(j),cc(j),ctot(j,ireg), &
          psi(j),gam(j)
          write(*,1000) iter,nam(j),y(j),cc(j),ctot(j,ireg), &
          psi(j),gam(j)
        enddo
        do i=1,ncmplx
          write(iunit2,1020) namcx(i),cx(i),gamx(i)
        enddo
        write(iunit2,1111) (nam(j),j=1,ncomp)
        write(iunit2,*)
        write(iunit2,*) 'jacobian zz(j,l)'
        do j = 1, ncomp
          write(iunit2,1021) nam(j),(zz(j,l), l = 1, ncomp)
        enddo
        write(iunit2,*)
        write(iunit2,*) 'dpsi(j,l)'
        do j = 1, ncomp
          write(iunit2,1021) nam(j),(dpsi0(j,l),l=1,ncomp)
        enddo

!-------check charge balance
        zchrg = zero
        do j = 1, ncomp
          zchrg = zchrg + z(j)*cc(j)
        enddo
        do i = 1, ncmplx
          zchrg = zchrg + zx(i)*cx(i)
        enddo
        zchrg0 = zero
        do j = 1, ncomp
          zchrg0 = zchrg0 + z(j)*psi(j)
        enddo
        write(iunit2,'("  charge balance: ",1p2e12.4)') zchrg,zchrg0
        write(*,'("  charge balance: ",1p2e12.4)') zchrg,zchrg0
      endif

 1000 format(i3,1x,a8,1p5g12.4)
 1010 format(' iter species    residual   molality   tot mol.     psi &
     &    gam')
 1020 format(20x,' complex ',a8,5x,1pg12.4,1pg12.4)
 1021 format(a8,(1p10g12.4))
 1111 format(8x,(16a12))

!-----call solver      
      call ludcmp(zz,ncomp,ncmx,indx,dd)
      call lubksb(zz,ncomp,ncmx,indx,y)

      alpha = one

!-----update concentration
      if (loglin .eq. 1) then ! linear
        do j = 1, ncomp
          fac = 2*alpha*y(j)/cc(j)
          dur = max(one, 2*alpha*y(j)/cc(j))
          chold = cc(j)
          cc(j) = cc(j)-alpha*y(j)/dur
          if (cc(j) .lt. zero) then
            write(*,'("warning neg. conc.: eqlib-",i1,i4,1x,a8, &
     &      1p10e12.4)') k,iter,nam(j),cc(j),alpha,y(j),dur,fac,chold, &
            guess(j,ireg)
            guess(j,ireg) = abs(guess(j,ireg))
!           cc(j) = guess(j,ireg)
!           cc(j) = abs(cc(j))/2
            alpha = half*alpha
!pcl        cc(j) = cc(j) + alpha*y(j)
          endif
        enddo
      else
        do j = 1, ncomp
          if (abs(y(j)).gt.1.e4) then
            write(*,*) 'error in eqlib: ill posed problem!',nam(j), &
            y(j),iter
            write(*,*) 'name       cc      psi      gam      residual'
            do i = 1, ncomp
              write(*,'(a12,1p4e12.4)') nam(i),cc(i),psi(i),gam(i), &
              y(i)
            enddo
            write(*,*) 'STOP!'
            stop
          endif
          resfac = exp(-alpha*y(j))
!         if (resfac .gt. 1.e2) resfac = guess(j,ireg)/cc(j)
          cc(j) = cc(j)*resfac
        enddo
      endif

      call treqres (cc,y,cx,psi,xex0,gam,gamx, &
      pgas,siteden,cec,csorpf,csorp,ccsorp,alogpf,dgamdi)

!-----check convergence
      newtconv = 0
      do j=1,ncomp
!       if (abs(y(j)).gt.atol(j)+rtol(j)*cc(j)) goto 20
!       tolj = min(0.001*cc(j),tol*atol(j))
        tolj = tol*atol(j)

!       if (tolj.lt.tolmin0) tolj=tolmin0
!       if (abs(y(j)).gt.tol*atol(j)) then
        if (abs(y(j)).gt.tolj) then
          newtconv = 1
        endif
        if (iwarn.gt.2) &
        write(*,'("eqlib: ",4i3,1x,a12,1p10e12.4)') iter,newtconv, &
        iactsav,istate,nam(j),cc(j),gam(j),y(j),tolj !,tol,atol(j)
      enddo

      if (newtconv .eq. 1) then
        goto 20
      else if (newtconv.eq.0 .and. iactsav.gt.0 .and. iact.eq.0 &
      .and. istate.eq.1) then
        iact = iactsav
        
!       write(*,*) 'treqlib: convergence with iact=0!',iter
!       do j = 1, ncomp
!         write(*,*) nam(j),cc(j),gam(j),y(j)
!       enddo

        goto 20
      endif

      if (newtconv.eq.0 .and. istate.eq.0) then

        chrg = zero
        chrgj = zero
        do j = 1, ncomp
          chrg = chrg + z(j)*psi(j)
          if (j.ne.jchrg) chrgj = chrgj + z(j)*psi(j)
        enddo
        ctchrg = -chrgj/z(jchrg)

!       x = abs(chrg/chrgj)
!       if (x.eq.1.d0) x=0.5d0
!       cave = x*ctot(jchrg,ireg)+(1.d0-x)*ctchrg
        cave = 0.95d0*ctot(jchrg,ireg)+0.05*ctchrg

        if (iwarn.gt.2) then
          write(*,'("treqlib: chrg ",i3,1x,a12,1x,2i2,1p5e12.4)') &
          iter,nam(jchrg),ireg,istate,x,chrg,cave
          if (myrank == 0) &
          write(iunit2,'("treqlib: chrg ",i3,1x,a12,1x,2i2,1p5e12.4)') & 
          iter,nam(jchrg),ireg,istate,x,chrg,cave
        endif

        ctot(jchrg,ireg) = cave

        iter = 0

!       if (abs(chrg).lt.1.d-3*cave) then
          istate = 1
          itype(jchrg,ireg) = -1
          ctot(jchrg,ireg) = ctotsav
!       endif

        goto 20
      endif

      if (iter.lt.2) goto 20

      if (iprint.ge.3 .and. myrank==0) then
        write(iunit2,'(" iter = ",i4)') iter
        write(*,'(" iter = ",i4)') iter
        write(*,*) 'species    ityp    residual     cc      tot       gam'
        write(iunit2,*) 'species    ityp    residual     cc     tot   gam'
        do j = 1, ncomp
          write(iunit2,1050) nam(j),itype(j,ireg),y(j),cc(j),psi(j),gam(j)
          write(*,1050) nam(j),itype(j,ireg),y(j),cc(j),psi(j),gam(j)
        enddo
      endif

      do j = 1, ncomp
        guess(j,ireg) = cc(j)
      enddo

      return

  333 continue
      if (myrank == 0) then
        write(iunit2,1040) iter,istate
        write(*,1040) iter,istate
        write(*,*) 'species    ityp    residual     cc      tot       gam'
        write(iunit2,*) 'species    ityp    residual     cc     tot   gam'
        do j = 1, ncomp
          write(iunit2,1050) nam(j),itype(j,ireg),y(j),cc(j),psi(j),gam(j)
          write(*,1050) nam(j),itype(j,ireg),y(j),cc(j),psi(j),gam(j)
        enddo
      endif
 1040 format(' no convergence in eqlib after iter =',i4,' istate=',i2) 
 1050 format(' ',a12,i3,1p5g12.4)
      stop
      end subroutine newtonraphson

      subroutine constrnt(y,zz,cc,gam,ah2o)

      use ptran_global_module
      use ptran_dbase_module
            
      implicit real*8 (a-h,o-z)

!     include 'impl.h'
!     include 'paramtrs.h'  
!     include 'trminrl.h'
!     include 'trgas.h'
!     include 'triounts.h'
!     include 'trscalar.h'
!     include 'trcom.h'

      real*8 zz(ncmx,ncmx),y(ncmx),gam(*),cc(*),usvd(ncmx,ncmx), &
      svd(ncmx),vsvd(ncmx,ncmx)

      real*8 :: eqhet(ncmx),shet(ncmx,ncmx)
      integer :: indx(ncmx),jjdx(ncmx)

      character*20 jndx(ncmx),name(ncmx)

      ntot = 0
!-----minerals
      do 10 l = 1, ncomp
        if(itype(l,ireg) .eq. 3) then
          do m = 1, mnrl
            if (ncon(l,ireg) .eq. namrl(m)) then
              ntot = ntot + 1
              jndx(ntot) = nam(l)
              if (isothrm .ge. 1) then
                ii = ncmplx + ncxkin + ngas + m
                alnk(m) = -(coef(ii,1)*log(tk) &
                +           coef(ii,2) &
                +           coef(ii,3)*tk &
                +           coef(ii,4)/tk &
                +           coef(ii,5)/(tk*tk))
              endif
              eqhet(ntot) = alnk(m)
              name(ntot) = namrl(m)
              if (molal.eq.1 .and. jh2o.eq.ncomp+1) then
                do j = 1, ncomp+1
                  shet(j,ntot) = smnrl(j,m)
                enddo
              else
                do j = 1, ncomp
                  shet(j,ntot) = smnrl(j,m)
                enddo
             endif

!-------------check initial guesses
!             spr = eqhet(ntot)*aln10
!             do ll = 1, ncomp
!               spr = spr + shet(ll,ntot)*log(gam(ll)*cc(ll))
!             enddo
!             write(*,'("eqlib: check min ",a8,a12,1p3g12.4)') nam(l),
!    .        namrl(m),spr,cc(l)

            goto 10
          endif
        enddo
        write(*,*) 'constraint mineral not found for species ',nam(l)
        if (myrank == 0) &
        write(iunit2,*) 'constraint mineral not found for species ', &
        nam(l)
      endif

   10 continue

!-----gases
      do 15 l = 1, ncomp
        if(itype(l,ireg).eq.4 .or. itype(l,ireg).eq.12) then
          do m = 1, ngas
            if (ncon(l,ireg) .eq. namg(m)) then
              ntot = ntot + 1
              jndx(ntot) = nam(l)
              if (isothrm .ge. 1) then
                ii = ncmplx + ncxkin + m
                eqgas(m) = -(coef(ii,1)*log(tk) &
                           + coef(ii,2) &
                           + coef(ii,3)*tk &
                           + coef(ii,4)/tk &
                           + coef(ii,5)/(tk*tk))
              endif

!             write(*,*) 'eqlib: ',l,m,k,eqgas(m),ctot(l,ireg)

              eqhet(ntot) = eqgas(m) - log10(ctot(l,ireg))
              name(ntot) = namg(m)
              if (molal.eq.1 .and. jh2o.eq.ncomp+1) then
                do j = 1, ncomp+1
                  shet(j,ntot) = sgas(j,m)
                enddo
              else
                do j = 1, ncomp
                  shet(j,ntot) = sgas(j,m)
                enddo
              endif
              goto 15
            endif
          enddo
          write(*,*) 'constraint gas not found for species ',nam(l)
          if (myrank == 0) &
          write(iunit2,*) 'constraint gas not found for species ', &
          nam(l)            
        endif

   15 continue

!-----solve for new guesses in terms of activities
      fach2o = one
      if (molal.eq.0) fach2o = one/(wh2o*cc(jh2o))
      do m = 1, ntot
        y(m) = -eqhet(m)*aln10
        if (iact.eq.6) then
          y(m) = y(m) - shet(jh2o,m)*ah2o
        endif
        jj = 0
 	do j = 1, ncomp
          iflg = -1
          do i = 1, ntot
	    if (jndx(i) .eq. nam(j)) iflg = 0
          enddo
          if (iflg .eq. 0) then
!-----------set up left hand side
	    jj = jj + 1
            jjdx(jj) = j
	    zz(m,jj) = shet(j,m)
	  else
!-----------set up right hand side
!           write(*,*) 'eqlib: ',m,j,gam(j),cc(j),shet(j,m)
            if (shet(j,m) .ne. zero) then
	      y(m) = y(m) - shet(j,m)*log(gam(j)*cc(j)*fach2o)
            endif
	  endif
        enddo
        if (jj .ne. ntot) then
          write(*,*) 'not enough constraints in solid: stop'
          if (myrank == 0) &
          write(iunit2,*) 'not enough constraints in solid: stop'
          stop
        endif
      enddo

!     write(*,*) 'zone = ',ntot   
!     do m = 1, ntot
!       write(*,100) jndx(m),y(m),(zz(m,l),l=1,ntot)
!     enddo
! 100 format(a8,1pe12.4,(1p10e12.4))

      if (ntot.gt.0) then
!-------check linear independence of mineral reactions
        do i=1,ntot
          do j=1,ntot
            usvd(i,j) = zz(i,j)
          enddo
        enddo
        call svdcmp(usvd,ntot,ntot,ncmx,ncmx,svd,vsvd)
        do j = 1, ntot
          if (svd(j).lt.1.e-12) then
            write(*,*) '        --> equilibrium problem ill posed! STOP'
	    write(*,1313) (svd(l),l=1,ntot)
	    if (myrank==0) write(iunit2,1313) (svd(l),l=1,ntot)
 1313       format(/' singular value decomposition:'/(1p10g12.4))
            stop
          endif
        enddo

        call ludcmp(zz,ntot,ncmx,indx,dd)
        call lubksb(zz,ntot,ncmx,indx,y)
      endif

!-----modify guesses
      do m = 1, ntot
!       write(*,*) 'eqlib: ',nam(jjdx(m)),m,jjdx(m),cc(jjdx(m)),y(m)
	cc(jjdx(m)) = exp(y(m))/gam(jjdx(m))
!       guess(jjdx(m),ireg) = cc(jjdx(m))
!       write(*,*) 'eqlib: ',nam(jjdx(m)),m,jjdx(m),cc(jjdx(m))
      enddo

      return
      end subroutine constrnt
 
      subroutine svdcmp(a,m,n,mp,np,w,v)
!-----------------------------------------------------------------------
!     Given a matrix A, with logical dimensions M by N and physical 
!     dimensions MP by NP, this routine computes its singular value 
!     decomposition, A=U W V^T. The matrix U replaces A on output. The 
!     diagonal matrix of singular values W is output as a vector W. The 
!     matrix V is output as V. M must be greater or equal to N: if it is 
!     smaller, then A should be filled up to square with zero rows.
!-----------------------------------------------------------------------

      parameter (nmax=100)

      implicit real*8 (a-h,o-z)

      dimension a(mp,np),w(np),v(np,np),rv1(nmax)
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if (i.le.m) then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if (scale.ne.0.0) then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            if (i.ne.n) then
              do 15 j=l,n
                s=0.0
                do 13 k=i,m
                  s=s+a(k,i)*a(k,j)
13              continue
                f=s/h
                do 14 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
            endif
            do 16 k= i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if ((i.le.m).and.(i.ne.n)) then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if (scale.ne.0.0) then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            if (i.ne.m) then
              do 23 j=l,m
                s=0.0
                do 21 k=l,n
                  s=s+a(j,k)*a(i,k)
21              continue
                do 22 k=l,n
                  a(j,k)=a(j,k)+s*rv1(k)
22              continue
23            continue
            endif
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if (i.lt.n) then
          if (g.ne.0.0) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=n,1,-1
        l=i+1
        g=w(i)
        if (i.lt.n) then
          do 33 j=l,n
            a(i,j)=0.0
33        continue
        endif
        if (g.ne.0.0) then
          g=1.0/g
          if (i.ne.n) then
            do 36 j=l,n
              s=0.0
              do 34 k=l,m
                s=s+a(k,i)*a(k,j)
34            continue
              f=(s/a(i,i))*g
              do 35 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
          endif
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if ((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if ((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            if ((abs(f)+anorm).ne.anorm) then
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.0/h
              c= (g*h)
              s=-(f*h)
              do 42 j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
42            continue
            endif
43        continue
2         z=w(k)
          if (l.eq.k) then
            if (z.lt.0.0) then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if (its.eq.30) pause 'no convergence in 30 iterations'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=sqrt(f*f+1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 nm=1,n
              x=v(nm,j)
              z=v(nm,i)
              v(nm,j)= (x*c)+(z*s)
              v(nm,i)=-(x*s)+(z*c)
45          continue
            z=sqrt(f*f+h*h)
            w(j)=z
            if (z.ne.0.0) then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 nm=1,m
              y=a(nm,j)
              z=a(nm,i)
              a(nm,j)= (y*c)+(z*s)
              a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      end subroutine svdcmp

!=========================================================================

      subroutine treqres (cc,y,cx,ppsi,xex0,gam,gamx, &
      pgas,siteden,cec,csorpf,csorp,ccsorp,alogpf,dgamdi)

      use ptran_global_module
      use water_eos_module
      use co2eos_module
	  
      
      implicit none

      real*8 :: cc(*),cx(*),ppsi(*),gam(*),gamx(*),alogpf(:,:),dgamdi(*), &
                cec(*),xex0(*),siteden(*),csorpf(*),csorp(*),ccsorp(*), &
                pgas(*),y(*)
      real*8 :: ajln(ncmx)
      
      integer :: i,ii,ierr,j,jj,k,l,m,n,ngs,nns,ncex,iflgam
      
      real*8 :: ah2o0,alogkeq,akeq,cco2,d1m,wmix,uwat,dwat, &
      poynco2,hpco2,x1m,xco2,xmwc,xmww,fco2,dco2,pco2,pco2pa,phico2,tc, &
      frac,fach2o,wfrac,spr,prod,sum,sumc,summ,por00=1.d0

      real*8 :: rmix,eqk,henry,vphi,poyn,xmol,xphico2
      real*8 :: rho(1)
      
      n = 1

      fach2o = one
      if (molal.eq.0) fach2o = one/(wh2o*cc(jh2o))

      do j = 1, nmass
        ajln(j) = log(gam(j)*cc(j)*fach2o)
      enddo

!-----compute concentrations of reversible aqueous complexes
      do i=1,ncmplx
        if (isothrm .ge. 1) then
              eqhom(i) = -(coef(i,1)*log(tk) &
                         + coef(i,2) &
                         + coef(i,3)*tk &
                         + coef(i,4)/tk &
                         + coef(i,5)/(tk*tk))
        endif
 
        cx(i) = zero
        prod  = eqhom(i)*aln10    !- log(gamx(i))
        do j = 1, nmass
          prod = prod+ajln(j)*shom(j,i)
        enddo

!-------add correction for activity of water using Pitzer model
!       ah2o0 = ln[ah2o]
        if (iact.eq.6) then
          prod = prod + ah2o0*shom(jh2o,i)
        endif

        if (prod .lt. 30.d0) then
!         if (gamx(i).gt.1.d-8) then
            cx(i) = exp(prod)/(fach2o*gamx(i))
!         else
!           gamx(i) = one
!           cx(i) = exp(prod)/(fach2o*gamx(i))
!         endif
        else
          cx(i) = one
        endif
        
!       write(*,*) 'treqres: ',namcx(i),ncomp,nmass,prod,cx(i),gamx(i), &
!       ah2o0,fach2o,eqhom(i),(shom(j,i),j=1,ncomp)
      enddo

!-----get activity coefficients
      if      (iact.ge.1 .and. iact.le.5) then
        call trgameq (cc,cx,gam,gamx,dgamdi,iflgam)
!     else if (iact .eq. 2) then
!       call actcoef(cc)
!       call hiact(ncomp,cc)
!     else if (iact .eq. 3) then
!       call actconst
      else if (iact .eq. 6) then
!       call trpitzer (ncomp,nmass,ncmplx,cc,cx,gam,gamx,phih2o, &
!       ah2o0,wh2o,jh2o,molal,sionic0,dbpath,iunit2,isothrm,tk,nam, &
!       namcx,1,1)
      endif

!-----compute total concentrations
      do j=1,ncomp
        if (nam(j) .ne. 'e-') then
          ppsi(j) = cc(j)
        else
          ppsi(j) = zero
        endif
        do i=1,ncmplx
          ppsi(j) = ppsi(j) + shom(j,i)*cx(i)
        enddo
      enddo

      if (nexsolid+nsrfmin.gt.0) then
        ncex = ncomp
        if (nexmax .gt. ncomp) then
          ncex = nexmax
        endif
!       call trionexc (ncomp,nexsite,nexmax,nsrfmin,nsrfmx,nsrfsit, &
!       gam,xex0,xex0,cec,cc,cx, &
!       ccsorp,csorp,csorpf,siteden,r,cdl,alogpf,dpsi,trdifl, &
!       dfrc,nelmrow,nrow,ndcon,maxnc,ncdiag, &
!       cjac,cxm,cxp,vl,area,ncex,n,n)
!       call trionexc (cc,por)
      endif

!-----set up constraint equations for distribution of species
      do j = 1, ncomp

        if (itype(j,ireg) .eq. 1) then

!---------total concentration
          y(j) = ppsi(j) - ctot(j,ireg)

        else if (itype(j,ireg).eq.2 .or. itype(j,ireg).eq.21) then
          
!---------total concentration: sorbed + aqueous
!---------ion exchange
          summ = zero
          sumc = zero
          do m = 1, nexsolid
            if (m.le.nexsolid-ncollex) then
              do k = nsitex1(m), nsitex2(m)
                do jj = nex1(k), nex2(k)
                  if (nam(j) .eq. namcat(jj)) then
                    summ = summ + cec(k)*xex0(jj)/z(j)
                  endif
                enddo
              enddo
            else
              do k = nsitex1(m), nsitex2(m)
                do jj = nex1(k), nex2(k)
                  if (nam(j) .eq. namcat(jj)) then
                    sumc = sumc + cec(k)*xex0(jj)/z(j)
                  endif
                enddo
              enddo
            endif
          enddo

!---------surface complexation
          do m = 1, nsrfmin
            if (m.le.nsrfmin-ncolsrf) then ! minerals
              do k = nsite1(m), nsite2(m)
                do i = nsorp1(k), nsorp2(k)
                  summ = summ + ssorp(j,i)*ccsorp(i)
                enddo
              enddo
            else                           ! colloids
              do k = nsite1(m), nsite2(m)
                do i = nsorp1(k), nsorp2(k)
                  sumc = sumc + ssorp(j,i)*ccsorp(i)
                enddo
              enddo
            endif
          enddo
          summ = summ/por00

          y(j) = ppsi(j) + sumc + summ - ctot(j,ireg)

!       else if (itype(j,ireg) .eq. 9) then
!         tyr = tpath*yrsec
!         ctot(j,ireg) = ci(j)+(one-exp(-alam(j)*tyr))*(cf(j)-ci(j))
!         y(j) = ppsi(j) - ctot(j,ireg)

        else if (itype(j,ireg).eq.3) then

!---------mineral constraint
          do m = 1, mnrl 
            nns = m
            if (namrl(m).eq.ncon(j,ireg)) goto 50
          enddo

          if (myrank==0) &
          write(iunit2,*) 'mineral name not found:',ncon(j,ireg)
          write(*,*) 'mineral name not found:',ncon(j,ireg)

          stop

   50     continue
          m = nns
          if (isothrm .ge. 1) then
              ii = ncmplx + ncxkin + ngas + m
              alnk(m) = -(coef(ii,1)*log(tk) &
                        + coef(ii,2) &
                        + coef(ii,3)*tk &
                        + coef(ii,4)/tk &
                        + coef(ii,5)/(tk*tk))
          endif

          prod = alnk(m)*aln10
          do l = 1, nmass
            prod = prod+ajln(l)*smnrl(l,m)
          enddo

!---------add correction for activity of water using Pitzer model
          if (iact.eq.6) then
            prod = prod + ah2o0*smnrl(jh2o,m)

!           write(*,*) 'treqres: ',iter,prod,ah2o0,smnrl(jh2o,m)
          endif

          if (prod.le.30.d0) then
            spr = exp(prod)
          else 
            spr = 1.e30
          endif
          
          y(j) = spr-one

        else if (itype(j,ireg).eq.4) then

!---------gas constraint
          do m = 1, ngas 
            ngs = m
            if (namg(m).eq.ncon(j,ireg)) goto 90
          enddo
          if (myrank == 0) &
          write(iunit2,*) 'gas name not found'
          write(*,*) 'gas name not found'

          goto 100
   90     continue

          if (isothrm .ge. 1) then
              ii = ncmplx + ncxkin + ngs
              eqgas(ngs) = -(coef(ii,1)*log(tk) &
                           + coef(ii,2) &
                           + coef(ii,3)*tk &
                           + coef(ii,4)/tk &
                           + coef(ii,5)/(tk*tk))
          endif
          
          if (jco2g == -1) then
            call henry_co2_noderiv(xmol,x1m,tc,pgas(jco2g)*1.d5, &
            1.d0,henry,poyn)
            
            !note: henry coef. = H/phi
            
            !xphico2 = ...
            
            ! (Vphi [cm^3/mol]: Garcia, 2001)
            vphi= 37.51 + (-9.585e-2 + 8.74e-4*tc - 5.044e-7*tc*tc)*tc
            
            ! fmwh2o, fmwco2: g/mol; rmix, rho: g/cm^3; wmix: g/mol
            ! ccloc_p: mol/L
            
            rmix = rho(1)+fmwh2o*cc(jco2)*1.d-3*(fmwco2-rho(1)*vphi)
            wmix = rmix+(fmwh2o-fmwco2)*cc(jco2)*1.d-3
            
!           eqk = henry*1.e-8*fmwh2o/xphico2/wmix
!           eqk = henry*1.e-5*fmwh2o*1.d-3/xphico2/(1.d0+cc(jco2)*1.d-3*fmwh2o)
            eqk = henry*1.e-5*fmwh2o*1.d-3!/xphico2

!           print *,'ptran_psi: ',n,ngs,temploc_p(n),henry, &
!           cc(jco2),log10(eqk),eqgas(jco2g), &
!           vphi,rmix,wmix,rgasjj*tk

            eqgas(jco2g) = log10(eqk)
          endif
          prod = eqgas(ngs)*aln10
          do l = 1, nmass
            prod = prod+ajln(l)*sgas(l,ngs)
          enddo

!---------add correction for activity of water using Pitzer model
          if (iact.eq.6) then
            prod = prod + ah2o0*sgas(jh2o,ngs)
          endif

          pgas(ngs) = exp(prod) 
          y(j) = pgas(ngs) - ctot(j,ireg)

  100     continue

        else if (itype(j,ireg).eq.5) then

!---------weight fraction solvent
          sum = zero
          do l = 1, nmass
            sum = sum + wt(l)*cc(l)
          enddo
          do l = 1, ncmplx
            sum = sum + wtx(l)*cx(l)
          enddo
          wfrac = ctot(jh2o,ireg)*1.d-2
          y(j) = cc(jh2o) - wfrac/(one-wfrac)*sum/wh2o

!         write(*,'("treqres: ",a6,1x,1p10e12.4)') nam(j), &
!         y(j),wfrac,sum,fach2o,cc(jh2o)

        else if (itype(j,ireg).eq.6) then

!---------mole fraction solvent
          sum = zero
          do l = 1, nmass
            sum = sum + cc(l)
          enddo
          do l = 1, ncmplx
            sum = sum + cx(l)
          enddo
          frac = ctot(jh2o,ireg)*1.d-2
          y(j) = cc(jh2o) - frac/(one-frac)*sum

!         fac1 = frac/(one-frac)*sum
!         write(*,'("treqres: ",i3,1x,a6,1x,1p10e12.4)') iter,nam(j), &
!         y(j),frac,sum,fach2o,cc(jh2o),fac1

        else if (itype(j,ireg).eq.7) then

!---------fixed concentration constraint
          if (j.eq.jph) then
            y(j) = cc(j)*gam(j)*fach2o-ctot(j,ireg)
!         else if (j.eq.joh) then
!           y(j) = cc(j)*fach2o-ctot(j,ireg)
          else
            y(j) = cc(j)-ctot(j,ireg)
          endif
 
        else if (itype(j,ireg) .eq. 8) then

!---------pH constraint
          if (j.eq.jph) then
            y(j) = cc(j)*gam(j)*fach2o-10.d0**(-ctot(j,ireg))
          else if (j.eq.joh) then

!-----------add correction for activity of water using Pitzer model
            if (iact.eq.6) then
              y(j) = cc(j)*gam(j)*fach2o &
              -10.d0**(eqhom(iph)+ctot(j,ireg))*exp(ah2o0)
            else
              y(j) = cc(j)*gam(j)*fach2o-10.d0**(eqhom(iph)+ctot(j,ireg))
            endif
          else
            y(j) = cc(j)*gam(j)*fach2o-10.d0**(-ctot(j,ireg))
          endif

        else if (itype(j,ireg) .eq. 9) then

!---------activity constraint
          y(j) = cc(j)*gam(j)*fach2o-ctot(j,ireg)
 
!       else if (itype(j,ireg) .eq. 10) then

!!---------conserved quantity
!!         y(j) = ppsi(5)+15./2.*ppsi(2)-45./4.*ppsi(1)

!         y(j) = ppsi(j)
!         do l = 1, nreactcq
!           y(j) = y(j) + sconsrv(j,l)*ppsi(l)
!         enddo

        else if (itype(j,ireg).eq.11 .or. itype(j,ireg).eq.12) then

!---------CO2(g) EOS

          pco2 = ctot(j,ireg)
          pco2pa = pco2*1.e5
          tc = tk-tkelvin

!         print *,'ptran_speciation1: ',tc,pco2,pco2pa
                    
          xmwc = 4.40098D-02
          xmww = 1.801534D-02

!---------cole & mcphearson eos based on mrk eos for co2 (not recommended!)
!         call co2(tc,pco2pa,dc,fc,phi,hc)
!         call sat(tc,ps)
!         call henry(tc,pco2pa,ps,fc,x1m,xco2,hp)

!---------duan, moller, weare (1992) eos with modified henry's law 
!         (recommended)
          call duanco2 (tc,pref0 *1D-5,dco2,fco2,phico2)
          
!         dco2 = [g/cm^3]

!         print *,'ptran_speciation2: ',tc,pco2,dco2,fco2,phico2
          
!          call henry_co2_noderiv (xco2,x1m,tc,pco2*1.e5,phico2,hpco2,poynco2)
          call Henry_duan_sun(pref0 *1D-5, tc,  hpco2)        
          !note: hpco2 = H/phico2 [Pa]
          
 !         call cowat (tc,pco2pa,dwat,uwat,ierr) ! why pco2pa here?
          
!         print *,'ptran-1: ',pco2pa,dwat,uwat,x1m
          
!          call denmix (tc,dwat,x1m,d1m)
          
 !         wmix = xmwc*xmww/(x1m*xmww+(1.d0-x1m)*xmwc) ! [kg/mol]
          
!          cco2 = x1m*d1m/xmwc*1.e-3
           cco2= pco2 * phico2 * hpco2

!         print *,'ptran_speciation3: ',xmww,xmwc,hpco2,poynco2,xco2,x1m,d1m,cco2

 !         akeq = wmix*hpco2*1.e-2/d1m/poynco2 ! units of akeq are bars (/(mol/L))
 !         alogkeq = log10(akeq)
 !         eqgas(jco2g) = alogkeq-log10(phico2)
          akeq = log(hpco2*phico2) / aln10 
          eqgas(jco2g)=  - akeq 
          !1.D0 / hpco2 /(fmwh2o * (akeq/(akeq+ cco2)))

!         print *,'ptran_speciation4: ',jco2,jco2g,wmix,akeq,alogkeq,eqgas(jco2g)
          
!         pstar = 10.d0**eqgas(jco2g)*cco2
          
!         write(*,*) 'treqres: ',pco2,cco2,phico2,tc,alogkeq, &
!         eqgas(jco2g),d1m,x1m,dwat,dco2
          
          if (itype(j,ireg).eq.11) y(j) = ppsi(j)-cco2
          if (itype(j,ireg).eq.12) y(j) = cc(j)-cco2
        
        else if (itype(j,ireg).eq.-1) then

!---------electroneutrality constraint
          y(j) = zero
          do l = 1, ncomp
            if (nam(j) .ne. 'e-') y(j) = y(j) + z(l)*ppsi(l)
          enddo

        endif
      enddo

      return
      end subroutine treqres

!=========================================================================

      subroutine treqjac (cc,cx,cloc,gam,dpsi0,pgas, &
      cec,xex0,siteden,ccsorp,zz)
      
      use ptran_global_module

      implicit none

      integer :: i,j,jj,k,l,ll,m,nns,ngs
      
      real*8 :: ccl,s1,s2,wfrac,sum,fac,fach2o,ffac,frac,por00=1.d0,zz00

      real*8 :: cx(*),pgas(*),cc(*),cloc(*),gam(*),ccsorp(*),siteden(*), &
      cec(*),xex0(*)

      real*8 :: zz(ncmx,ncmx),dpsi0(ncmx,ncmx),didm0(ncmx),corr(ncmx), &
      dxex(nexmx,nexmx),dk(nexmx),sigi(ncmx)

      fach2o = one
      if (molal.eq.0) fach2o = one/(wh2o*cc(jh2o))

      if (iact .ge. 1) then
        do j = 1, ncomp
          corr(j) = zero
!         do i = 1,ncmplx
!           corr(j) = corr(j) + shom(j,i)*dcdi(i)
!         enddo
        enddo
      endif

      do l = 1,ncomp
        ccl = cc(l) 
        didm0(l) = zero

!-------activity coefficient correction to jacobian-not activated
        if(iact .eq. -1) then
          s1 = zero
          s2 = zero
          do i = 1, ncmplx
            if (zx(i).ne.0. .and. shom(l,i).ne.zero) then
              s1 = s1+zx(i)*zx(i)*shom(l,i)*cx(i)
            endif
!           s2 = s2+zx(i)*zx(i)*dcdi(i)
          enddo
          if (ccl .gt. zero) then
            didm0(l) = half*(z(l)*z(l)+s1/ccl)/(one-half*s2)
          endif
        endif

        if (molal.eq.1 .and. jh2o.eq.ncomp+1) then

          do j = 1, ncomp
            dpsi0(j,l) = zero
            if(j.eq.l .and. nam(j).ne.'e-') then
              dpsi0(j,l) = one
            else
              dpsi0(j,l) = zero
            endif
            do i = 1, ncmplx 
              if (ccl.gt.zero .and. shom(j,i).ne.zero) then
                dpsi0(j,l) = dpsi0(j,l)+shom(j,i)*shom(l,i)*cx(i)/ccl
              endif
            enddo
!pcl        if (iact .ge. 1) dpsi0(j,l) = dpsi0(j,l) + corr(j)*didm0(l)
          enddo

        else

          if (l.le.nmass) then

            do j = 1, ncomp
              dpsi0(j,l) = zero
              if(j.eq.l .and. nam(j).ne.'e-') then
                dpsi0(j,l) = one
              else
                dpsi0(j,l) = zero
              endif
              do i = 1, ncmplx 
                if (ccl.gt.zero .and. shom(j,i).ne.zero) then
                  dpsi0(j,l) = dpsi0(j,l)+shom(j,i)*shom(l,i)*cx(i)/ccl
                endif
              enddo
            enddo

          else if (l.eq.jh2o) then

            do i = 1, ncmplx
              sigi(i) = one
              do j = 1, nmass
                sigi(i) = sigi(i) - shom(j,i)
              enddo
            enddo
            do j = 1, ncomp
              dpsi0(j,l) = zero
              if(j.eq.l .and. nam(j).ne.'e-') then
                dpsi0(j,l) = one
              else
                dpsi0(j,l) = zero
              endif
              do i = 1, ncmplx 
                if (ccl.gt.zero .and. shom(j,i).ne.zero) then
                  dpsi0(j,l) = dpsi0(j,l)+shom(j,i)*sigi(i)*cx(i)/ccl
                endif
              enddo
            enddo
          endif
        endif
      enddo

      do j = 1, ncomp

        do l = 1, ncomp
          zz(j,l) = zero
        enddo

        if(itype(j,ireg).eq.1) then

!---------total concentration constraint
          do l = 1, ncomp
            zz(j,l) = dpsi0(j,l)
          enddo

        else if (itype(j,ireg).eq.2 .or. itype(j,ireg).eq.21) then

!---------ion-exchange and surface complexation constraints
          do l = 1, ncomp
            zz(j,l) = dpsi0(j,l)
          enddo

          do m = 1, nexsolid
            if (m.le.nexsolid-ncollex) then
              fac = one/por00
            else
              fac = one
            endif
            do k = nsitex1(m), nsitex2(m)

!-------------gaines-thomas convention: xex = z_j Cbar_j/sum z_i Cbar_i
              sum = zero
              do ll = nex1(k), nex2(k) 
                l = jex(ll)
                sum = sum+z(l)*xex0(ll)
                dk(ll) = xex0(ll)/cc(l)
              enddo
              do jj = nex1(k), nex2(k)
                if (nam(j) .eq. namcat(jj)) then
                  ffac = -z(j)*xex0(jj)/sum
                  do ll = nex1(k), nex2(k)
                    dxex(jj,ll) = ffac*dk(ll)
                  enddo
                  dxex(jj,jj) = dxex(jj,jj) + dk(jj)
                  do ll = nex1(k), nex2(k)
                    l = jex(ll)
                    zz(j,l) = zz(j,l)+cec(k)/z(j)*dxex(jj,ll)*fac 
                  enddo
                endif
              enddo
            enddo
          enddo

!---------surface complexation
          if (nsrfmin .gt. 0) then
            do l = 1, ncomp
              do m = 1, nsrfmin
                if (m.le.nsrfmin-ncolsrf) then ! minerals
                  fac = one/por00
                else
                  fac = one
                endif
                do k = nsite1(m), nsite2(m)
                  do i = nsorp1(k), nsorp2(k)
                    zz(j,l) = zz(j,l) + ssorp(j,i)*ssorp(l,i)*ccsorp(i)/ &
                    cc(l)*fac
                  enddo
                  zz(j,l) = zz(j,l) - dcsorp(j,k)*dcsorp(l,k)/ &
                  (siteden(k)*cc(l))*fac
                enddo
              enddo
            enddo
          endif

        else if (itype(j,ireg).eq.3) then

!---------mineral constraint
          if (j.eq.jh2o) then
            write(*,*) 'illegal constraint for water'
            stop
          endif

          do m = 1, mnrl 
            nns = m
            if(namrl(m).eq.ncon(j,ireg)) goto 50
          enddo
          if (myrank == 0) &
          write(iunit2,*) 'mineral name not found:',ncon(j,ireg)
          write(*,*) 'mineral name not found:',ncon(j,ireg)
          stop

   50     continue
          m = nns
          sum = zero
!         if (iact .ge. 1) then
!           do l = 1, ncomp
!             sum = sum + smnrl(l,m)*dgamdi(l)
!            enddo
!         endif
          do l = 1, ncomp            
            ccl=cc(l) 
            zz(j,l) = zero
            if (ccl.gt.zero) then
              zz(j,l) = smnrl(l,m)/ccl + sum*didm0(l)
            endif
          enddo

        else if (itype(j,ireg).eq.4) then

!---------gaseous constraint
          if (j.eq.jh2o) then
            write(*,*) 'illegal constraint for water'
            stop
          endif

          do m = 1, ngas 
            ngs = m
            if(namg(m).eq.ncon(j,ireg)) goto 90
          enddo
          if (myrank == 0) &
          write(iunit2,*) 'gas name not found'
          write(*,*) 'gas name not found'
          stop

   90     continue
          sum = zero
!         if (iact .ge. 1) then
!           do l = 1, ncomp
!             sum = sum + sgas(l,ngs)*dgamdi(l)
!           enddo
!         endif
          do l = 1, ncomp            
            ccl = cc(l) 
            zz(j,l) = zero
            if (ccl.gt.zero) then
              zz(j,l) = (sgas(l,ngs)/ccl + sum*didm0(l))*pgas(ngs)
            endif
          enddo

  100     continue

!       else if (itype(j,ireg).eq.5) then

!---------concentration buffer constraint
!         do l = 1, ncomp
!           if (j.eq.l) then
!             zz(j,l) = one
!           else
!             zz(j,l) = zero
!           endif
!         enddo

        else if(itype(j,ireg).eq.5) then

!---------weight fraction solvent
          wfrac = ctot(jh2o,ireg)*1.d-2
          fac = wfrac/(one-wfrac)/wh2o
          do l = 1, ncomp
            if (l.eq.jh2o) then
              zz00 = one
            else
              zz00 = zero
            endif
            sum = zero
            do i = 1, nmass
              sum = sum + wt(i)*dpsi0(i,l)
            enddo
            zz(j,l) = zz00-fac*sum
          enddo

        else if(itype(j,ireg).eq.6) then

!---------mole fraction solvent
          frac = ctot(jh2o,ireg)*1.d-2
          do l = 1, ncomp
            if (l.eq.jh2o) then
              zz(j,l) = one
            else
              zz(j,l) = -frac/(one-frac)
            endif
          enddo

        else if(itype(j,ireg).eq.7) then

!---------concentration constraint
          do l = 1, ncomp
            if(j.eq.l) then
              if (j .eq. jph) then
                zz(j,l) = gam(j)*fach2o
              else
                zz(j,l) = one
              endif
            else if(j.ne.l) then
              if (j .eq. jph) then
!               zz(j,l) = gam(j)*cc(j)*dgamdi(j)*didm0(l)
                zz(j,l) = zero
              else
                zz(j,l) = zero
              endif
            endif
          enddo

        else if(itype(j,ireg).eq.8) then

!---------pH constraint
          do l = 1, ncomp
            if(j.eq.l) then
              if (j .eq. jph) then
                zz(j,l) = gam(j)*fach2o
              else
                zz(j,l) = one*fach2o
              endif
            else if(j.ne.l) then
              if (j .eq. jph) then
!               zz(j,l) = gam(j)*cc(j)*dgamdi(j)*didm0(l)
                zz(j,l) = zero
              else
                zz(j,l) = zero
              endif
            endif
          enddo

        else if(itype(j,ireg).eq.9) then
 
          do l = 1, ncomp
            if(j.eq.l) then
              if (j .eq. jph) then
!               zz(j,l) = gam(j)*(one + cc(j)*dgamdi(j)*didm0(j))
                zz(j,l) = one
              else
                zz(j,l) = one
              endif
            else if(j.ne.l) then
              if (j .eq. jph) then
!               zz(j,l) = gam(j)*cc(j)*dgamdi(j)*didm0(l)
                zz(j,l) = zero
              else
                zz(j,l) = zero
              endif
            endif
          enddo

!       else if (itype(j,ireg) .eq. 10) then

!---------conserved quantity
!         do l = 1, ncomp
!           zz(j,l) = dpsi0(j,l)
!           do i = 1, nreactcq
!             zz(j,l) = zz(j,l) + sconsrv(j,i)*dpsi0(i,l)
!           enddo
!         enddo

        else if(itype(j,ireg).eq.11) then

!---------total CO2(g) constraint
          do l = 1, ncomp
            zz(j,l) = dpsi0(j,l)
          enddo

        else if(itype(j,ireg).eq.12) then

!---------CO2(aq) constraint
          do l = 1, ncomp
            if(j.eq.l) then
              zz(j,l) = one
            else if(j.ne.l) then
              zz(j,l) = zero
            endif
          enddo

        else if(itype(j,ireg).eq.-1) then

!---------electroneutrality constraint
          do l = 1, ncomp 
            zz(j,l) = zero 
            do ll = 1, ncomp 
              if (nam(j).ne.'e-') zz(j,l) = zz(j,l)+z(ll)*dpsi0(ll,l)
            enddo
          enddo

        endif
      enddo

      return
      end subroutine treqjac

!=========================================================================

      subroutine trgameq(cc,cx,gam,gamx,dgamdi,iflgam)
      
      use ptran_global_module
      
      implicit none

      integer :: i,iflgam,it,itmax,j

      real*8 :: cc(*),cx(*),gam(*),gamx(*),dgamdi(*)
      real*8 :: alogac(ncmx),alogacx(ncxmx)
      real*8 :: afac,ajn,dcdi,didi,f,fach2o,fcomp,prod,si,sum,sqrti

      fach2o = one
      if (molal.eq.0) fach2o = one/(wh2o*cc(jh2o))

      itmax = 50
      si = sionic0
      
      fcomp = zero
      do j=1,ncomp 
        if (nam(j) .ne. 'e-') fcomp = fcomp + z(j)*z(j)*cc(j)
      enddo

      if (fcomp .eq. zero) return

!     do i = 1, ncxkin
!       fcomp = fcomp + zxk(i)*zxk(i)*ccxkin(i)
!     enddo

      it=0
!**************************************************************** 
!     begin loop in activity coefficient iteration
!**************************************************************** 
   10 continue
      it=it+1 
      if (it .gt. itmax) then
        write(*,1000) it,si
        if (myrank == 0) write(iunit2,1000) it,si
        iflgam = 3
        return
      endif
 1000 format(' *** excessive iterations in ionic strength loop *** &
     &iter =',i3,' ionic strength =',1pg12.4)
      sum = zero
      do i=1,ncmplx
        sum = sum+zx(i)*zx(i)*cx(i)
      enddo
!**************************************************************** 
!     compute new ionic strength
!**************************************************************** 
      f = half*(sum+fcomp)*fach2o

!**************************************************************** 
!     compute derivative of reduced complex concentration 
!                    wrt ionic strength
!**************************************************************** 
      if(ncmplx.gt.0) then
        sqrti=sqrt(si)
        do j = 1, ncomp
          dgamdi(j) = zero
          if(z(j).ne.zero) then
            afac=one+bdebye*a0(j)*sqrti
            dgamdi(j)=-half*adebye*z(j)*z(j)/(sqrti*afac*afac) &
            + bextend(j)
          endif
        enddo
        didi=zero
        do i=1,ncmplx
          if(zx(i).ne.zero) then
            afac = one+bdebye*ax0(i)*sqrti
            sum=half*adebye*zx(i)*zx(i)/(sqrti*afac*afac) - bextendx(i)
            do j=1,ncomp
              if(z(j).ne.zero) then
                sum=sum+shom(j,i)*dgamdi(j)
              endif
            enddo
!           dcdi(i)=cx(i)*aln10*sum
!           didi=didi+half*zx(i)*zx(i)*dcdi(i)
            dcdi=cx(i)*aln10*sum
            didi=didi+half*zx(i)*zx(i)*dcdi
          endif
        enddo
        didi = didi

        si=(f-si*didi)/(one-didi)
        if (si .lt. zero) then
          write(*,1001) it,si
!pcl      if (myrank == 0) write(iunit2,1001) it,si
          iflgam = 3
!         return
          si = -si
        endif
      else
        si = f
      endif

 1001 format(' ionic strength negative: ',i3,1pg12.4)

      sqrti=sqrt(si) 
      do j=1,ncomp 
        if(z(j).ne.zero) then
          alogac(j)=-z(j)*z(j)*adebye*sqrti/(one+bdebye*a0(j)*sqrti) &
          +bextend(j)*si
          gam(j)=exp(alogac(j)*aln10)
        endif
      enddo

      do i=1,ncmplx
        if(zx(i).ne.zero) then
          alogacx(i)=-zx(i)*zx(i)*adebye*sqrti/(one+bdebye*ax0(i) &
          *sqrti)+bextendx(i)*si
          gamx(i)=exp(alogacx(i)*aln10)
        endif
!**************************************************************** 
!       compute concentrations of aqueous complexes
!**************************************************************** 
        cx(i)=zero
        if (isothrm .ge. 1) then
              eqhom(i) = -(coef(i,1)*log(tk) &
              + coef(i,2) &
              + coef(i,3)*tk &
              + coef(i,4)/tk &
              + coef(i,5)/(tk*tk))
        endif
        prod=eqhom(i)*aln10    !-log(gamx(i)) 
        do j=1,nmass 
          if(shom(j,i).ne.zero) then 
            ajn=cc(j)*gam(j)*fach2o
            if(ajn.gt.zero) prod=prod+log(ajn)*shom(j,i)
          endif   
        enddo
!       do k = 1, ncxkin
!         if (shom(ncomp+k,i) .ne. 0.) then
!           ak = gamxk(k)*ccxkin(k)
!           if (ak .gt. 0) then
!             prod = prod + log(ak)*shom(ncomp+k,i)
!           endif
!         endif
!       enddo
        cx(i)=exp(prod)/(fach2o*gamx(i))
      enddo

      if(abs(f-si)/si.lt.1.d-8) goto 20

!**************************************************************** 
      goto 10
!**************************************************************** 

   20 continue

      sionic0 = si

      return
      end subroutine trgameq
      
 !===============================================================

      subroutine trsolprd(gam,cc,psi,pgas,ah2o) 

      use ptran_global_module
      
      implicit none
      
      integer :: ii,j,k,m,mm,ngs,ns,indx(nmmx)

      real*8 :: gam(*),cc(*),psi(*),pgas(*),prod(nmmx)
      
      real*8 :: ak,aj,ah2o,fach2o,prod1,prod10,pot00

!**************************************************************** 
!     compute solubility products for heterogeneous equilibria
!**************************************************************** 
      fach2o = one
      if (molal.eq.0) fach2o=one/(wh2o*cc(jh2o))

      do ns = 1, mnrl
        if (isothrm .ge. 1) then
          ii = ncmplx + ncxkin + ngas + ns
          alnk(ns) = -(coef(ii,1)*log(tk) &
                     + coef(ii,2) &
                     + coef(ii,3)*tk &
                     + coef(ii,4)/tk &
                     + coef(ii,5)/(tk*tk))
        endif

!-------convert log K to base e
        prod(ns) = alnk(ns)*aln10
        do j = 1, nmass
          if(smnrl(j,ns).ne.zero) then
            aj = cc(j)*gam(j)*fach2o
            if(aj.le.zero) goto 60 
            prod(ns) = prod(ns)+log(aj)*smnrl(j,ns) 
          endif
        enddo

!-------add correction for activity of water based on Pitzer model
        if (iact.eq.6) then
          prod(ns) = prod(ns)+ah2o*smnrl(jh2o,ns) 
        endif

        if (ze0(ns).ne.zero) then
          prod(ns) = prod(ns) - ze0(ns)*faraday/(rgasj*tk)*pot00
        endif
        if (abs(prod(ns)) .lt. 1.e-9) prod(ns) = zero
   60   continue
      enddo

!-----order saturation indices in descending order
      indx(1) = 1
      if (mnrl .gt. 1) call indexx(mnrl,prod,indx)
      if (myrank == 0) then
        write(iunit2,'(72("-"))')
        write(iunit2,990) 
  990   format(/' ',15x,'mineral saturation indices')
        write(iunit2,'(" mineral ",13x," SI (Log KQ)",5x, &
     &  " log K")')
        write(iunit2,'(72("-"))')
        do m = 1, mnrl
          mm = indx(mnrl-m+1)
          prod1 = prod(mm)/aln10
          write(iunit2,'(1x,a20,1pe13.5,2x,1pe13.5)') namrl(mm),prod1, &
          alnk(mm)
        enddo
      endif
!**************************************************************** 
!     compute gas partial pressures
!**************************************************************** 

      if (ngas.gt.0 .and. myrank==0) then
        write(iunit2,*)
        write(iunit2,*) ' gas   log partial pressure', &
        '  pressure [bars] ','  log K'
        write(iunit2,'(72("-"))')
      endif

      do ngs=1,ngas 
      
      ! write(*,*) 'trsolprd: ',ngs,isothrm,iflgco2,jco2g
        
        if (isothrm.ge.1 .and. iflgco2.eq.0) then
          ii = ncmplx + ncxkin + ngs
          eqgas(ngs) = -(coef(ii,1)*log(tk) &
          + coef(ii,2) &
          + coef(ii,3)*tk &
          + coef(ii,4)/tk &
          + coef(ii,5)/(tk*tk))
        else if (isothrm.eq.1 .and. iflgco2.eq.1) then
          if (ngs.ne.jco2g) then
            ii = ncmplx + ncxkin + ngs
            eqgas(ngs) = -(coef(ii,1)*log(tk) &
            + coef(ii,2) &
            + coef(ii,3)*tk &
            + coef(ii,4)/tk &
            + coef(ii,5)/(tk*tk))
          endif
        endif
        
!       print *,'ptran_speciation: ',ngs,eqgas(ngs),fach2o,iflgco2
        
        prod(ngs) = eqgas(ngs)*aln10
        do k = 1,nmass
          if (itype(k,ireg).eq.11) then
            ak = psi(k)*gam(k)*fach2o
          else
            ak = cc(k)*gam(k)*fach2o
          endif
          if (ak.gt.zero) prod(ngs) = prod(ngs)+log(ak)*sgas(k,ngs)
          
        ! print *,'ptran_speciation*: ',k,ngs, itype(k,ireg),sgas(k,ngs),ctot(k,ireg),cc(k),psi(k),prod(ngs)
        enddo

!-------add correction for activity of water based on Pitzer model
        if (iact.eq.6) then
          prod(ngs) = prod(ngs)+ah2o*sgas(jh2o,ngs) 
        endif

        pgas(ngs) = exp(prod(ngs))
        prod10 = prod(ngs)/aln10
        
	!	print *,'ptran_speciation**: ',  eqgas(ngs), exp(eqgas(ngs))
		if (myrank==0) &
        write(iunit2,'(1x,a12,1pg12.4,5x,1pg12.4,4x,1pg12.5)') &
        namg(ngs),prod10,pgas(ngs),eqgas(ngs)
      enddo
    ! print *, jh2o,jco2,skin
      return
      end subroutine trsolprd

end module ptran_speciation_module
