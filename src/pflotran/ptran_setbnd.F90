!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_setbnd.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_setbnd.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.2  2004/04/06 17:39:17  lichtner
! Set rho = 1 temporarily for debugging purposes.
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

module ptran_setbnd_module

  public

  contains

      subroutine trsetbnd 

      use ptran_global_module
      use trdynmem_module
      use water_eos_module
      use ptran_speciation_module
      use co2eos_module
    
      implicit none

      real*8 :: cloc(ncmx),cxloc(ncxmx),cloctot(ncmx),pgasloc(ngmx), &
      cecloc(nexmx),xexloc(nexmx), &
      ccsorplc(nscxmx),csorplc(nscxmx),csorpflc(nscxmx), &
      sitdnloc(nscxmx), &
      dgamdi(ncmx),gamloc(ncmx),gamxloc(ncxmx)

      real*8 :: fach2o,rho1,sum,u1, dco2,fc,phi

      integer :: i,idif=0,iireg,iisrc,iter,j,jj,k,l,m,nss

      character :: psinam*7

      ibc = 0
      if(ibcreg(1).eq.0) goto 110 ! source/sink

!-----------------------------------------------------------------------
!-----set boundary condition concentrations
!-----------------------------------------------------------------------
      do 100 ireg = ibcreg(1), ibcreg(2) 

        ibc = ibc+1

!-------skip if zero-flux boundary condition
        if(ibndtyp(ibc).eq.2) goto 100

        if(mode.eq.2) then
          tk = tempbc(ibc)+tkelvin
          rho1 = dwbc(ibc)
        else
!         tk = tbc(ibc)+tkelvin
!         rho1 = dwbc(ibc)*0.0180153d0 ! convert from moles/m3 to gms/cc
!         if (jh2o.gt.0) rho1=one
        endif
        rho1 = 1.d0

!-------initialize activity coefficients to unity
        do j = 1, ncomp
          gamloc(j) = one
        enddo
        do i = 1, ncmplx
          gamxloc(i) = one
        enddo

!-------write to output file
        if (myrank == 0) &
        write (iunit2,1015) ibc,ibndtyp(ibc),iface(ibc),tempbc(ibc)
        
!-------speciate boundary solution
        iireg = ireg
        
!       print *,'ptransetbnd: ',ireg, ibcreg(1), ibcreg(2), initreg
        
        call trspeciate (cloc,cxloc,cloctot,pgasloc, &
        cecloc,xexloc,csorplc,ccsorplc,csorpflc,sitdnloc,alogpf, &
        gamloc,gamxloc,dgamdi,rho,tempbc(ibc),iter)

        print *,"Setbnd:pgas : ", pgasloc
!-------write to screen
        if (myrank ==0) &
        write (*,1010) ibndtyp(ibc),tempbc(ibc),iter
        
        fach2o = one
        if (molal.eq.0) fach2o = one/(wh2o*cloc(jh2o))

        psinam = ' psibnd'
        do j = 1, ncomp
!         psibnd(j,ibc) = cloctot(j) * rho1
          psibnd(j,ibc) = cloctot(j)
          ccbnd (j,ibc) = cloc(j)
          if (myrank==0) &
          write(iunit2,1020) nam(j),psinam,j,ibc,psibnd(j,ibc)
        enddo

        psinam = ' xexbnd'
        do m = 1, nexsolid
          do k = nsitex1(m), nsitex2(m)
            do jj = nex1(k), nex2(k)
              j = jex(jj)
              xexbnd(jj,ibc) = xexloc(jj)
              if (myrank==0) &
              write(iunit2,1020) nam(j),psinam,jj,ibc,xexbnd(jj,ibc)
            enddo
          enddo
        enddo

        if (nsrfmin .gt. 0) then
          psinam = ' psorpbnd'
          do j = 1, ncomp
            sum = zero
            do m = nsrfmin-ncolsrf+1,nsrfmin ! sum over colloids only !
              do k = nsite1(m), nsite2(m)
                do i = nsorp1(k), nsorp2(k)
                  sum = sum + ssorp(j,i)*ccsorplc(i)
                enddo
              enddo
            enddo
            psorpbnd(j,ibc) = sum
            if (myrank==0) &
            write(iunit2,1020) nam(j),psinam,j,ibc,sum
          enddo
        endif

        if (iphase.eq.2 .or. iphase.eq.0) then !  conc for the gas phase
          call duanco2(tempbc(ibc), pgasloc(1),dco2,fc,phi)
         
      u1= (dco2/fmwco2*1D3) /pgasloc(1)
           print *,"ptran_set_BC: ", ibc,  tempbc(ibc), pgasloc(1),dco2,u1
      psinam = 'psigbnd'
          do j = 1, ncomp
            sum = zero
            do l = 1, ngas
              sum = sum + sgas(j,l)*pgasloc(l)
            enddo
            psigbnd(j,ibc) = sum * u1                        
            if (myrank==0) &
            write(iunit2,1020) nam(j),psinam,j,ibc,psigbnd(j,ibc)
          enddo
          do i = 1, ngas
            pgasbnd(i,ibc) = pgasloc(i)!*u1      
      print *,"ptran_setbnd: ",  psigbnd(i,ibc) ,  pgasloc(i)             
          enddo
        endif
  100 continue

  110 if(ibcreg(3).eq.0) return

!-----------------------------------------------------------------------
!-----set source/sinks
!-----------------------------------------------------------------------
      nss = 0
      iisrc = 0
      
      do ireg = ibcreg(3), ibcreg(4)

        tk = tempbc(ireg)+tkelvin
        if (pref0 .gt. zero) then
          call density (tempbc(ireg),pref0,rho1)
          rho1 = rho1*1.d-3
        else
          rho1 = 1.d0
        endif
        rho1 = 1.d0

        do j = 1, ncomp
          gamloc(j) = one
        enddo
        do i = 1, ncmplx
          gamxloc(i) = one
        enddo

        iisrc = iisrc+1
        iireg = ireg
        
!       do j = 1,ncomp
!       print *,'ptran_setbnd1: ',j,ireg,ibcreg(3),ibcreg(4),itype(j,ireg),ctot(j,ireg)
!       enddo
        if (myrank==0) write(iunit2,1000) iisrc,tempbc(ireg)

        call trspeciate (cloc,cxloc,cloctot,pgasloc, &
        cecloc,xexloc,csorplc,ccsorplc,csorpflc,sitdnloc,alogpf, &
        gamloc,gamxloc,dgamdi,rho,tempbc(ireg),iter)
        
!       print *,'ptran_setbnd2: ',ireg,ctot(jco2,ireg)

        nss = nss+1

        if (myrank==0) write(*,1005) nss,tempbc(ireg),iter

        if (myrank==0) write(iunit2,'(/,"source/sink concentrations")')
        psinam = ' psisrc'

!-------set psisrc for the total region
        do k = ks1(iisrc), ks2(iisrc)
          do j = js1(iisrc), js2(iisrc)
            do i = is1(iisrc), is2(iisrc)

!             nss = nss+1
!             m = ndloc1*(i+jl(j)+kl(k))-ndloc2
              m = i+(j-1)*nx+(k-1)*nxy ! global index!
              do l = 1, ncomp
!               psisrc(l,nss) = cloctot(l)*rho(ng)
                psisrc(l,nss) = cloctot(l)*rho1

!               write(*,*) 'trsetbnd: ',myrank,ireg,m,i,j,k,l,nss, &
!               cloctot(l),rho1

                if (myrank == 0) &
                write(iunit2,1020) nam(l),psinam,l,nss,psisrc(l,nss)
              enddo
            enddo
          enddo
        enddo
      enddo

  120 continue

 1000 format (/,70('=')/1x,' source composition: region = ',i3, &
      '  temp = ',1pg10.3/,70('=')/)
 1005 format (/,5x,'--> compute source composition:',i3, &
      '  temp = ',1pg10.3,'  iter =',i4)
 1010 format (5x,'--> boundary condition:  type =',i2, &
      '  temp = ',1pg10.3,'  iter =',i4)   
 1011 format (5x,'--> boundary condition:  type',i2, &
      '  temp = ',1pg10.3,'  iter =',i4,' matrx')   
 1012 format (5x,'--> boundary condition:  type',i2, &
      '  temp = ',1pg10.3,'  iter =',i4,' frac')   
 1015 format (/,70('=')/1x,' boundary condition: region = ',i3, &
      '  type = ',i2,' face = ',i2,'  temp = ',1pg10.3/,70('=')/)   
 1016 format (/,70('=')/1x,' boundary condition:  type',i2, &
      '  temp = ',1pg10.3,' matrx'/,70('=')/)   
 1017 format (/,70('=')/1x,' boundary condition:  type',i2, &
      '  temp = ',1pg10.3,' frac'/,70('=')/)   
 1020 format(10x,'component = ',a12,3x,a7,'(',i2,i3,') =',1p10e12.4)

      return
      end subroutine trsetbnd

end module ptran_setbnd_module
