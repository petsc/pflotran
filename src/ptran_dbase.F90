!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_dbase.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_dbase.F90,v $
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

module ptran_dbase_module

contains

  subroutine ptran_dbase (namc)

  use ptran_global_module
! use fileio_module
  
  implicit none
  
  save
  
      character(len=namlen) :: namc(*)
      
      integer :: i,ii,ix,j,jj,k,l,lhp,lp,m, &
      n,n0,nn0,nbasis,ngs,nngs,nr, &
      iflg,iflgint,iflgck,iflgseck,iflgsec,iflge,ifound, &
      ntemp,nreac,nsec,nsrf,nsite

      real*8 zz,wwt,tmpk,dd,aa0

      real*8 :: alogk0(ntmpmx),tempc(ntmpmx),vec(5,ntmpmx),bvec(5)
      integer :: indx(ndimmx)
      integer :: indxpri(ncmx),indxsec(ncxmx0),indxmin(nmmx0), &
      indxgas(ngmx0),indxirr(nxkmx),indxsrf(nscxmx)
      real*8 s(20)
      real*8 skpri0(ncxmx0,nxkmx)

      real*8 :: alogk(ndimmx),amat(ndimmx,ndimmx), &
      bmat(ndimmx,ncmx+nxkmx),wmat(ndimmx,ndimmx),coef0(ndimmx,5)

  character(len=strlen) :: homepath
  character(len=dbaselen) :: string
  character(len=namlen) :: namsec(ndimmx),namdb(20),namvec(ndimmx)

      n0 = ndimmx

!-----subroutine database orders an independent set of chemical
!     reactions in terms of an arbitrary set of primary species or 
!     components.

!-----------------------------------------------------------------------
!     database structure
!-----------------------------------------------------------------------
!     temperatures at which data is stored
!     primary species data:      species a0 z gfwt
!     aqueous secondary species: species reaction logk a0 z gfwt
!     gaseous secondary species: species molar vol reaction logk z gfwt
!     mineral data:              mineral molar vol reaction logk gfwt
!     surface complexes:         species reaction logk z gfwt
!-----------------------------------------------------------------------
!     definitions
!-----------------------------------------------------------------------
!     secondary species = aqueous complexes + gaseous species
!-----------------------------------------------------------------------

!-----open appropriate database file

      homepath = dbpath

      lhp = len_trim(homepath)
      string(1:lhp) = homepath(1:lhp)
      if (myrank==0) then
        write(*,*) '--> using database: ', &
        homepath(1:len_trim(homepath))
        write(iunit2,*) '--> using database: ', &
        homepath(1:len_trim(homepath))
      endif
      open(iunit5, file=string(1:lhp),action="read",status='old',err=333)

      nbasis  = 5
      iflgint = 0
      iflgck  = 0

!-----set flags for species check
      do j = 1, ncpri
        indxpri(j)= 0
      enddo
      do i = 1, ncmplx
        indxsec(i)= 0
        namvec(i) = namcx(i)
      enddo
      do i = 1, ncxkin
        indxirr(i)= 0
!       namvec(ncpri+i) = namcxk(i)
        namvec(ncmplx+ngas+i) = namcxk(i)
      enddo
      do m = 1, mnrl
        indxmin(m)= 0
      enddo
      do l = 1, ngas
        indxgas(l)= 0
        namvec(ncmplx+l) = namg(l)
      enddo
      do m = 1, nsrfmx
        indxsrf(m)= 0
      enddo

!-----initialize matrix
      do i = 1, ndimmx
        do j = 1, nbasis
          coef(i,j) = zero
        enddo
        do l = 1, ndimmx
          amat(i,l) = zero
        enddo
      enddo

!=======================================================================
!-----begin reading data base
!=======================================================================

      read(iunit5,*,err = 334) namdb(1),ntemp,(tempc(l),l=1,ntemp)

      if (ntemp .gt. ntmpmx) then
        if (myrank==0) write(*,*) &
        'too many temperature points in database!'
        stop
      endif

!-----store interpolation vectors
      if (ntemp .gt. 1) then
        do i = 1, ntemp
          tmpk = tempc(i) + tkelvin
          vec(1,i) = log(tmpk)
          vec(2,i) = one
          vec(3,i) = tmpk
          vec(4,i) = one/tmpk
          vec(5,i) = one/(tmpk*tmpk)
        enddo

        do j = 1, nbasis
          do k = j, nbasis
            wmat(j,k) = zero
            do i = 1, ntemp
              wmat(j,k) = wmat(j,k) + vec(j,i)*vec(k,i)
            enddo
            if (j .ne. k) wmat(k,j) = wmat(j,k)
          enddo
        enddo
        call ludcmp(wmat,nbasis,n0,indx,dd)
        call fit(n0,nbasis,ntemp,iflgint,alogkeh,wmat,coefeh,vec,indx)
        alnkeh = flogk(coefeh,tempini+tkelvin)
      else if (ntemp .eq. 1) then
        alnkeh = alogkeh(2)
      endif

!=======================================================================
!-----read block #1
!-----read and store reversible and irreversible 
!     aqueous species properties
!=======================================================================

      if (iprint==2 .and. myrank==0) then
        write(iunit2,*) 
        write(iunit2,*) 'begin detailed database printout!'
      endif

   10 continue
      read(iunit5,*,err=334) namdb(1),aa0,zz,wwt

      if (iprint==2 .and. myrank==0) write(iunit2,*) &
      'component: ',namdb(1)

      if (namdb(1) .eq. 'null') goto 20

      if (icase .eq. 1) then      ! convert to lower case
        call convrtlc (namdb(1),20)
      else if (icase .eq. 2) then ! convert to upper case
        call convrtuc (namdb(1),20)
      endif

      do j = 1, ncpri
      
!       write(*,*) 'ptrandbase: ',namdb(1),namc(j)

        if (namdb(1) .eq. namc(j)) then
          wt(j) = wwt*1.d-3
          z(j)  = zz
          a0(j) = aa0
          if (iprint==2 .and. myrank==0) &
          write(iunit2,'("species found: ",a12)') namdb(1)
          indxpri(j) = 1
          goto 10
        endif
      enddo

      do i = 1, ncmplx
        if (namdb(1) .eq. namcx(i)) then
          wtx(i) = wwt*1.d-3
          zx(i)  = zz
          ax0(i) = aa0
          if (iprint==2 .and. myrank==0) &
          write(iunit2,'("species found: ",a12)') namdb(1)
          indxsec(i) = 1
          goto 10
        endif
      enddo

      do nsrf = 1, nsrfmin
        do nsite = nsite1(nsrf), nsite2(nsrf)
          if (namdb(1) .eq. namsite(nsite)) then
            zsite(nsite) = zz
            goto 10
          endif
        enddo
      enddo

      goto 10

   20 continue

!=======================================================================
!-----read block #2
!-----read and store reversible and irreversible aqueous reactions
!=======================================================================

      nreac = 0
      
   30 continue

      read(iunit5,*,err=334) namdb(1),n,(s(i+1),namdb(i+1),i=1,n), &
      (alogk0(l),l=1,ntemp),aa0,zz,wwt

      if (iprint==2 .and. myrank==0) write(iunit2,*) &
      namdb(1),(namdb(i+1),i=1,n)

      if (namdb(1) .eq. 'null') goto 70

      s(1) = -one

!-----decide if reaction occurs in system

      iflg = 0
      iflge = 0
      iflgsec = 0
      iflgseck = 0
      do 40 k = 1, n+1

        if (icase .eq. 1) then      ! convert to lower case
          call convrtlc (namdb(k),20)
        else if (icase .eq. 2) then ! convert to upper case
          call convrtuc (namdb(k),20)
        endif

!-------check for water
        if (namdb(k) .eq. 'H2O' .or. namdb(k) .eq. 'h2o') then
          goto 40
        endif

!-------search for primary species
        do j = 1, ncpri
          if (iprint==2 .and. myrank==0) &
          write(iunit2,*) 'comp: ',k,namdb(k),namc(j)
          if (namdb(k) .eq. namc(j)) then
            if (k .eq. 1) then
              iflg = 1
              jj = j
            endif
            goto 40
          endif
        enddo
        
!-------search for aqueous homogeneous kinetic reactions
        do nr = 1, ncxkin
          do lp = nparcxk1(nr),nparcxk2(nr)
            if (iprint==2 .and. myrank==0) &
            write(iunit2,*) 'ncxkin: ',k,namdb(k),namrxnaq(lp)
!           write(*,*) 'trdatbse-ncxkin: ',namdb(k),namrxnaq(lp)
            if (namdb(k) .eq. namrxnaq(lp)) then
!             if (k .eq. 1) then
!               iflg = 1
!               jj = j
!             endif
              goto 40
            endif
          enddo
        enddo 

!-------search for secondary species
        do i = 1, ncmplx
          if (iprint==2 .and. myrank==0) &
          write(iunit2,*) 'cmplx: ', k,namdb(k),namcx(i)
          if (namdb(k) .eq. namcx(i)) then
            iflgsec = 1
            if (k .eq. 1) then
              iflg = 2
              ii = i
            endif
            goto 40
          endif
        enddo

!-------search for gaseous species
        do ngs = 1, ngas
          if (iprint==2 .and. myrank==0) &
          write(iunit2,*) 'gas: ',k,namdb(k),namg(ngs)
          if (namdb(k) .eq. namg(ngs)) then
            nngs = ngs
            goto 40
          endif
        enddo

!-------search for electron
        if (namdb(k) .eq. 'e-') then
          iflge = 1
          goto 40
        endif

!-------species not found - skip reaction
        goto 30 ! get new reaction

   40 continue

      if (iprint==2 .and. myrank==0) &
      write(iunit2,'("species found: ",a12,i3, &
    & 10(f7.3,1x,a8))') namdb(1),n,(s(i+1),namdb(i+1),i=1,n)

!-----check off species found
      k = 1

!-----primary species
      do j = 1, ncpri        
        if (namdb(k) .eq. namc(j)) then
          indxpri(j) = 1
          goto 41
        endif
      enddo
   41 continue
   
!-----kinetic aqueous reactions
      do nr = 1, ncxkin
        do lp = nparcxk1(nr),nparcxk2(nr)
!         write(*,*) 'trdatbse: ',namrxnaq(lp),namdb(k)
          if (namdb(k) .eq. namrxnaq(lp)) then
            do j = 1, ncpri
!             write(*,*) 'trdatbse: ',namcxk(nr),namrxnaq(lp), &
!             namdb(k),namc(j)
              if (namcxk(nr) .eq. namc(j)) then
                indxpri(j) = 1
                goto 44
              endif
            enddo
          endif
        enddo
      enddo
   44 continue

!-----secondary species
      do i = 1, ncmplx
        if (namdb(k) .eq. namcx(i)) then
          indxsec(i) = 1
          goto 42
        endif
      enddo
   42 continue

!-----gaseous species
      do ngs = 1, ngas
        if (namdb(k) .eq. namg(ngs)) then
          indxgas(ngs) = 1
          goto 43
        endif
      enddo
   43 continue

!=======================================================================
!-----check if reaction is reversible (local equilibrium) 
!     or irreversible (kinetic)
!-----irreversible (kinetic)

      do l = 1, ncxkin
        do lp = nparcxk1(l),nparcxk2(l)
        
!         write(*,*) 'trdatbse: ',namcxk(l),namrxnaq(lp),namdb(1)
        
!         if (namdb(1) .eq. namcxk(l)) then
          if (namdb(1) .eq. namrxnaq(lp)) then
            indxirr(l) = 1
!-----------minus sign because lth species on rhs side of reaction
            eqcxk(lp) = -alogk0(1)
            axk0(l)  = aa0
            zxk(l)   = zz
            wtxk(l)  = wwt*1.d-3
          
!           write(*,*) n,lp,namcxk(l),namdb(1),namrxnaq(lp)

            do k = 2, n+1
!-------------primary species
              do j = 1, ncpri
                if (namdb(k) .eq. namc(j)) then
                  skpri(j,lp) = s(k)
                
!                 write(*,*) k,namdb(k),namc(j),s(k)
                
                  goto 133
                endif
              enddo
!-------------secondary species
              do i = 1, ncmplx
                if (namdb(k) .eq. namcx(i)) then
                  skpri0(i,lp) = s(k)
                  goto 133
                endif
              enddo
!-------------gaseous species
              do ngs = 1, ngas
                if (namdb(k) .eq. namg(ngs)) then
                  skpri0(ncmplx+ngs,lp) = s(k)
                  goto 133
                endif
              enddo
              if (namdb(k) .eq. 'e-') then
                zeaq(l) = s(k)
                goto 133
              endif
  133         continue
            enddo
            goto 30 ! read in new species
          endif
        enddo
      enddo

!-----reversible reaction (local equilibrium)

      nreac = nreac + 1

      namsec(nreac) = namdb(1)

      if (ntemp .gt. 1) then
        call fit(ndimmx,nbasis,ntemp,iflgint,alogk0,wmat,bvec,vec,indx)
  
        if (iflgint==1) then
          if (myrank == 0) &
          write(*,'(" -> warning: LogK = 500 encountered for species: ", &
      &   a12,/,4x,1p10e10.2)') namsec(nreac),(alogk0(l),l=1,ntemp)
          if (myrank == 0) &
          write(iunit2,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e10.2)') namsec(nreac),(alogk0(l),l=1,ntemp)
          iflgint = 2
        else if (iflgint .eq. 2) then
          iflgint = 0
        endif

        do j = 1, nbasis
          coef0(nreac,j) = bvec(j)
        enddo

        alogk(nreac) = flogk(bvec,tempini+tkelvin)
        
!       write(*,*) 'ptrandbase: ',nreac,alogk(nreac),(coef0(nreac,j), &
!       j=1,nbasis)

      else if (ntemp .eq. 1) then
        alogk(nreac) = alogk0(1)
      endif

      if (iflg .eq. 1) then
        wt(jj) = wwt*1.d-3
        z(jj)  = zz
        a0(jj) = aa0
      else if (iflg .eq. 2) then
        wtx(ii) = wwt*1.d-3
        zx(ii)  = zz
        ax0(ii) = aa0
      else if (iflg .eq. 3) then
        wtxk(ii) = wwt*1.d-3
        zxk(ii)  = zz
        axk0(ii) = aa0
      endif

!-----rearrange reactions - check for primary species
      do i = 1, n+1
        do j = 1, ncpri
          if (namdb(i) .eq. namc(j)) then
            bmat(nreac,j) = s(i)
            goto 60
          endif
        enddo

!-------check for secondary species
        do k = 1, ncmplx
          if (namdb(i) .eq. namcx(k)) then
            amat(nreac,k) = -s(i)
            goto 60
          endif
        enddo
        do ngs = 1, ngas
          if (namdb(i) .eq. namg(ngs)) then
            amat(nreac,ncmplx+ngs) = -s(i)
            goto 60
          endif
        enddo
   60   continue
      enddo

      goto 30

   70 continue

      if (myrank==0) &
      write(*,*) '        --> finished reading aqueous species'

!------------------end of aqueous complexing reactions------------------

      call dbgas (namc,namsec,ntemp,coef0,alogk,amat,bmat,wmat,bvec,vec, &
      indx,nreac,nbasis,indxgas)
      
      nsec = ncmplx + ngas
      if (myrank == 0) then
        write(iunit2,'(/)')
        write(iunit2,*) '            -->> matrix dimension: ',nsec
        write(*,*)      '            -->> matrix dimension: ',nsec
      endif
      
      if (nsec .gt. ndimmx) then
        if (myrank == 0) &
        write(*,*) 'matrix dimension too large in database! ',ndimmx
        stop
      endif

      if (nreac .ne. nsec) then
        if (myrank == 0) write(iunit2,101) nreac,nsec
        if (myrank == 0) write(*,101) nreac,nsec
        if (myrank == 0) write(iunit2,'("***missing species***")')
        if (myrank == 0) write(*,'("***missing species***")')
  101   format(' nreac =',i4,' nsec = ',i4)

        nn0 = nsec
        if (nreac .gt. nsec) nn0 = nreac
        do 120 ix = 1, ncmplx
          do i = 1, nn0
            if (namcx(ix).eq.namsec(i)) goto 120
          enddo
          if (myrank == 0) write(iunit2,'(i3,2x,a20)') ix,namcx(ix)
          if (myrank == 0) write(*,'(i3,2x,a20)') ix,namcx(ix)
  120   continue

        do 130 ngs = 1, ngas
          do i = 1, nn0
            if (namg(ngs).eq.namsec(i)) goto 130
          enddo
          if (myrank == 0) write(iunit2,'(i3,2x,a20)') ngs,namg(ngs)
          if (myrank == 0) write(*,'(i3,2x,a20)') ngs,namg(ngs)
  130   continue

        do i = 1, nn0
          if (myrank == 0) write(iunit2,'(i3,2x,a20)') i,namsec(i)
        enddo
        goto 190 
      endif

      if (nsec > 0) &
      call inverse (n0,ntemp,nbasis,nsec,namsec,namvec,alogk, &
                                  wmat,bmat,amat,coef0,namc)

      call dbmin (namc,ntemp,wmat,bvec,vec,indx,indxmin,nbasis,nreac)
      
      call dbsrf (namc,ntemp,wmat,bvec,vec,indx,indxsrf,nbasis,nreac)
      
      call dbtrf (namc,namsec,tempc,ntemp,nsec,nreac,nbasis)
      
!-----close data file
      close(iunit5)

      ifound = 0
      do j = 1, ncpri
        if (indxpri(j).ne.1 .and. &
          (namc(j).ne.'h2o' .or. namc(j).ne.'H2O')) then
          ifound = 1
          if (myrank == 0) write(iunit2,*) 'primary species not found: ', &
          namc(j)
          if (myrank == 0) write(*,*) 'primary species not found: ',namc(j)
        endif
      enddo
      do j = 1, ncmplx
        if (indxsec(j) .ne. 1) then
          ifound = 1
          if (myrank == 0) write(iunit2,*) 'secondary species not found: ', &
          namcx(j)
          if (myrank == 0) write(*,*) 'secondary species not found: ',namcx(j)
        endif
      enddo
      do j = 1, mnrl
        if (indxmin(j) .ne. 1) then
          ifound = 1
          if (myrank == 0) write(iunit2,*) 'mineral not found: ',namrl(j)
          if (myrank == 0) write(*,*) 'mineral not found: ',namrl(j)
        endif
      enddo
      do j = 1, nsrfmx
        if (indxsrf(j) .ne. 1) then
          ifound = 1
          if (myrank == 0) write(iunit2,*) 'surface complex not found: ', &
          namscx(j)
          if (myrank == 0) write(*,*) 'surface complex not found: ',namscx(j)
        endif
      enddo
      do j = 1, ngas
        if (indxgas(j) .ne. 1) then
          ifound = 1
          if (myrank == 0) write(iunit2,*) 'gas not found: ',namg(j)
          if (myrank == 0) write(*,*) 'gas not found: ',namg(j)
        endif
      enddo

      if (ifound .eq. 1) stop
      
      return

  333 if (myrank == 0) write(*,*) 'error in opening database file: stop'
      stop
      
  334 if (myrank == 0) write(*,*) 'error in reading database file: ',namdb(1),' stop'
      stop

  190 continue

    if (myrank == 0) then
      write(iunit2,*) 'STOP -- number of reactions not equal to ', &
      'number of secondary species!'
      write(*,*) 'STOP -- number of reactions not equal to ', &
      'number of secondary species!'
     
      write(iunit2,*) ' iflgck = ',iflgck
      do j = 1, ncpri
        if (indxpri(j) .ne. 1) then
          write(iunit2,*) 'primary species not found: ',namc(j)
          write(*,*) 'primary species not found: ',namc(j)
        endif
      enddo
      do j = 1, ncmplx
        if (indxsec(j) .ne. 1) then
          write(iunit2,*) 'secondary species not found: ',namcx(j)
          write(*,*) 'secondary species not found: ',namcx(j)
        endif
      enddo
      do j = 1, mnrl
        if (indxmin(j) .ne. 1) then
          write(iunit2,*) 'mineral not found: ',namrl(j)
          write(*,*) 'mineral not found: ',namrl(j)
        endif
      enddo
      do j = 1, ngas
        if (indxgas(j) .ne. 1) then
          write(iunit2,*) 'gaseous not found: ',namg(j)
          write(*,*) 'gaseous not found: ',namg(j)
        endif
      enddo
      do j = 1, ncxkin
        if (indxirr(j) .ne. 1) then
          write(iunit2,*) 'kinetic aqueous species not found: ',namcxk(j)
          write(*,*) 'kinetic aqueous species not found: ',namcxk(j)
        endif
      enddo
      do j = 1, nsrfmx
        if (indxsrf(j) .ne. 1) then
          write(iunit2,*) 'species not found: ',namscx(j)
          write(*,*) 'species not found: ',namscx(j)
        endif
      enddo
  
      stop
    endif
    
  end subroutine ptran_dbase
  
!========================================================================

  subroutine dbgas (namc,namsec,ntemp,coef0,alogk,amat,bmat,wmat,bvec,vec, &
  indx,nreac,nbasis,indxgas)

  use ptran_global_module
! use fileio_module
  
  implicit none
  
  save
  
      character(len=namlen) :: namc(*)
      
      integer :: i,ix,j,k,l, &
      n,nbasis,ngs,ntemp, &
      iflgint,nreac

      real*8 vbargas,wwt

      real*8 :: alogk0(ntmpmx),vec(5,ntmpmx),bvec(5)
      integer :: indx(ndimmx)
      integer :: indxgas(*)
      real*8 s(20)
      
      real*8 :: alogk(ndimmx),amat(ndimmx,ndimmx), &
                bmat(ndimmx,ncmx+nxkmx), &
                wmat(ndimmx,ndimmx),coef0(ndimmx,5)

  character(len=namlen) :: namsec(ndimmx),namdb(20)

!=======================================================================
!-----read block #3
!-----read and store gas reactions
!=======================================================================

   80 continue

      read(iunit5,*,err=334) namdb(1),vbargas,n,(s(i+1),&
      namdb(i+1),i=1,n),(alogk0(l),l=1,ntemp),wwt

!     write(*,'(a12,1pe12.4,i3,10(1pe12.4,1x,a8))') &
!     namdb(1),vbargas,n,(s(i+1),namdb(i+1),i=1,n)

      if (iprint==2 .and. myrank==0) write(iunit2,*) namdb(1)

      if (namdb(1) .eq. 'null') goto 110

      s(1) = -one

!-----decide if reaction occurs in chosen system

      do 90 i = 1, n+1

        if (icase .eq. 1) then      ! convert to lower case
          call convrtlc (namdb(i),20)
        else if (icase .eq. 2) then ! convert to upper case
          call convrtuc (namdb(i),20)
        endif

!-------check for water
        if (namdb(i) .eq. 'H2O' .or. namdb(i) .eq. 'h2o') then
          goto 90
        endif

        do ngs = 1, ngas
          if (namdb(i) .eq. namg(ngs)) then
            indxgas(ngs) = 1
            wtgas(ngs) = wwt*1.d-3
            goto 90
          endif
        enddo

        do j = 1, ncpri
          if (namdb(i) .eq. namc(j)) then
            goto 90
          endif
        enddo

        do ix = 1, ncmplx
          if (namdb(i) .eq. namcx(ix)) then
            goto 90
          endif
        enddo

        if (namdb(i) .eq. 'e-') then
          goto 90
        endif

!-------species not found - skip reaction
        goto 80

   90 continue

      if (iprint==2 .and. myrank==0) write(iunit2,*) 'species found: ', &
      namdb(1),vbargas,n,(s(i+1),namdb(i+1),i=1,n)

!     write(*,'("species found: ",a8,1pe12.4,i3,10(1pe12.4,1x,a8))') &
!     namdb(1),vbargas,n,(s(i+1),namdb(i+1),i=1,n)

      nreac = nreac + 1

      namsec(nreac) = namdb(1)

      if (ntemp .gt. 1) then
        call fit(ndimmx,nbasis,ntemp,iflgint,alogk0,wmat,bvec,vec,indx)

        if (iflgint .eq. 1) then
          if (myrank == 0) &
          write(*,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e9.2)') namsec(nreac),(alogk0(l),l=1,ntemp)
          if (myrank == 0) &
          write(iunit2,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e9.2)') namsec(nreac),(alogk0(l),l=1,ntemp)
          iflgint = 2
        else if (iflgint .eq. 2) then
          iflgint = 0
        endif

        do j = 1, nbasis
          coef0(nreac,j) = bvec(j)
        enddo
   
        alogk(nreac) = flogk(bvec,tempini+tkelvin)
        
      else if (ntemp .eq. 1) then
        alogk(nreac) = alogk0(1)
      endif

      do 100 i = 1, n+1

!-------rearrange reactions - check for primary species
        do j = 1, ncpri
          if (namdb(i) .eq. namc(j)) then
            bmat(nreac,j) = s(i)
            goto 100
          endif
        enddo

!-------check for secondary species
        do k = 1, ncmplx
          if (namdb(i) .eq. namcx(k)) then
            amat(nreac,k) = -s(i)
            goto 100
          endif
        enddo

        do ngs = 1, ngas
          if (namdb(i) .eq. namg(ngs)) then
            amat(nreac,ncmplx+ngs) = -s(i)
            goto 100
          endif
        enddo

  100 continue
          
      goto 80

  110 continue

      if (myrank == 0) &
      write(*,*) '        --> finished reading gaseous species'
      
      return
      
  334 if (myrank==0) write(*,*) 'error in reading database file: ',namdb(1),' stop'
      stop

!-------------------------end of gas reactions--------------------------

  end subroutine dbgas

!=======================================================================

  subroutine dbmin (namc,ntemp,wmat,bvec,vec,indx,indxmin,nbasis,nreac)
  
  use ptran_global_module
  
  implicit none
  
      character(len=namlen) :: namc(*)
      
      integer :: i,ix,j,k,l,m,mm, &
      n,nbasis,ngs,ndim,ntemp, &
      iflgint,nreac

  real*8 vbar0,wwt

      real*8 :: alogk0(ntmpmx),vec(5,ntmpmx),bvec(5)
      integer :: indx(ndimmx)
      integer :: indxmin(*)
      real*8 s(20)
      
      real*8 :: alogk(ndimmx),wmat(ndimmx,ndimmx)

  character(len=namlen) :: namdb(20)
  
!=======================================================================
!-----read block #4
!-----read and store mineral reactions
!=======================================================================

      ndim = nreac

      do j = 1, ncpri
        do m = 1, mnrl
          smnrl(j,m) = zero
        enddo
      enddo

      if (iprint==2 .and. myrank==0) write(iunit2,*) 'begin mineral block'

  140 continue

      read(iunit5,*,err=334) namdb(1),vbar0,n,(s(i+1),namdb(i+1),i=1,n), &
      (alogk0(l),l=1,ntemp),wwt

      if (iprint==2 .and. myrank==0) write(iunit2,*) namdb(1),vbar0,wwt

      if (namdb(1) .eq. 'null') goto 180

!-----decide if reaction occurs in chosen system

      if (icase .eq. 1) then      ! convert to lower case
        call convrtlc (namdb(1),20)
      else if (icase .eq. 2) then ! convert to upper case
        call convrtuc (namdb(1),20)
      endif
 
      do mm = 1, mnrl
      
!       write(*,*) 'dbmin: ',namrl(mm),namdb(1)
        
        if (namdb(1) .eq. namrl(mm)) then
          indxmin(mm) = 1
          m = mm
          goto 150
        endif
      enddo

      goto 140

  150 continue

      do 160 i = 2, n+1

        if (icase .eq. 1) then      ! convert to lower case
          call convrtlc (namdb(i),20)
        else if (icase .eq. 2) then ! convert to upper case
          call convrtuc (namdb(i),20)
        endif

!-------check for water
        if (namdb(i) .eq. 'H2O' .or. namdb(i) .eq. 'h2o') then
          goto 160
        endif

        do j = 1, ncpri
      
!         write(*,*) 'dbmin: ',namrl(mm),namdb(i),namc(j),ncpri,ncmplx
        
          if (namdb(i) .eq. namc(j)) then
            goto 160
          endif
        enddo

        do ix = 1, ncmplx
          if (namdb(i) .eq. namcx(ix)) then
            goto 160
          endif
        enddo

        do ngs = 1, ngas
          if (namdb(i) .eq. namg(ngs)) then
            goto 160
          endif
        enddo

        if (namdb(i) .eq. 'e-') then
          goto 160
        endif

        do mm = 1, mnrl
          if (namdb(i) .eq. namrl(mm)) then
            goto 160
          endif
        enddo

!-------mineral not found - skip reaction
        if (myrank == 0) &
        write(*,*) 'error reading database: mineral not found: ', &
        namdb(1),' STOP!'
        if (myrank == 0) &
        write(iunit2,*) 'error reading database: mineral not found: ', &
        namdb(1),' STOP!'
        stop

  160 continue

      if (iprint==2 .and. myrank==0) write(iunit2,*) 'mineral found: ',namdb(1)

      wtmin(m) = wwt*1.d-3
      vbar(m) = vbar0*1.d-3

      if (ntemp .gt. 1) then

        call fit(ndimmx,nbasis,ntemp,iflgint,alogk0,wmat,bvec,vec,indx)
          
        if (iflgint .eq. 1) then
          if (myrank == 0) &
          write(*,'(" -> warning: LogK = 500 encountered for species: ", &
      &   a12,/,4x,1p10e9.2)') namdb(1),(alogk0(l),l=1,ntemp)
          if (myrank == 0) &
          write(iunit2,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e9.2)') namdb(1),(alogk0(l),l=1,ntemp)
          iflgint = 2
        else if (iflgint .eq. 2) then
          iflgint = 0
        endif

        alogk(m) = flogk(bvec,tempini+tkelvin)

        do j = 1, nbasis
          coef(ndim+m,j) = bvec(j)
        enddo

      else if (ntemp .eq. 1) then
        alogk(m) = alogk0(1)
      endif

!-----rearrange reactions - check for primary species
      do 172 i = 1, n
        do j = 1, ncpri
          if (iprint==2 .and. myrank==0) write(iunit2,*) 'pri: ', &
          namdb(i+1),namc(j)
          if (namdb(i+1) .eq. namc(j)) then
            smnrl(j,m) = s(i+1)
            goto 172
          endif
        enddo
  172 continue

      do i = 1, n

!-------check for secondary species
        do k = 1, ncmplx
          if (namdb(i+1) .eq. namcx(k)) then
            do j = 1, ncpri
              smnrl(j,m) = smnrl(j,m) + shom(j,k)*s(i+1)
            enddo
            if (ntemp .gt. 1) then
              do l = 1, nbasis
                coef(ndim+m,l) = coef(ndim+m,l) + s(i+1)*coef(k,l)
              enddo
            endif
!...........note minus sign: eqhom = -eqhom as read in from database
            alogk(m) = alogk(m) - s(i+1)*eqhom(k) 
            goto 174
          endif
        enddo

!-------check for gaseous species
        do ngs = 1, ngas
          if (namdb(i+1) .eq. namg(ngs)) then
            do j = 1, ncpri
              smnrl(j,m) = smnrl(j,m) + sgas(j,ngs)*s(i+1)
            enddo
            if (ntemp .gt. 1) then
              do l = 1, nbasis
                coef(ndim+m,l) = coef(ndim+m,l) + s(i+1)*coef(ncmplx+ngs,l)
              enddo
            endif
!...........note minus sign: eqgas = -eqgas as read in from database
            alogk(m) = alogk(m) - s(i+1)*eqgas(ngs)
            goto 174
          endif
        enddo

!-------check for electron
        if (namdb(i+1) .eq. 'e-') then
          ze(m) = s(i+1)
        endif

!-------check for minerals
!       do mm = 1, mnrl
!         if (namdb(i+1) .eq. namrl(mm)) then
!           smin(m,mm) = s(i+1)
!           goto 174
!         endif
!       enddo

  174   continue
      enddo

!-----reverse sign to be consistent with mineral on rhs of reaction
      alnk(m) = -alogk(m)

      goto 140

!-----------------------end of mineral reactions------------------------
  180 continue

      if (myrank == 0) &
      write(*,*) '        --> finished reading minerals'

      return
      
  334 if (myrank == 0) write(*,*) 'error in reading database file: ', &
      namdb(1),' stop'
      stop

  end subroutine dbmin
  
!========================================================================
  
  subroutine dbsrf (namc,ntemp,wmat,bvec,vec,indx,indxsrf,nbasis,nreac)
  
  use ptran_global_module
  
  implicit none
  
! save
  
      character(len=namlen) :: namc(*)
      
      integer :: i,ix,j,k,l,m,mm,n, &
                 nbasis,ngs,ndim,ns,ntemp, &
      iflgint,nreac,nsrf,nsite

  real*8 zz,wwt

      real*8 :: alogk0(ntmpmx),vec(5,ntmpmx),bvec(5)
      integer :: indx(ndimmx)
      integer :: indxsrf(*)
      real*8 s(20)
      
  real*8 :: alogk(ndimmx),wmat(ndimmx,ndimmx)

  character(len=namlen) :: namdb(20)
  
!=======================================================================
!-----read block #5
!-----read and store surface complexation reactions
!=======================================================================

      ndim = nreac + mnrl

!-----initialize stoichiometric coeffs.
      do m = 1, nsrfmin
        do l = nsite1(m), nsite2(m)
          do i = nsorp1(l), nsorp2(l)
            do j = 1, ncpri
              ssorp(j,i) = zero
            enddo
          enddo
        enddo
      enddo

      if (iprint==2 .and. myrank==0) write(iunit2,*) 'begin surface ', &
      'complexation block'

  145 continue

      read(iunit5,*,err=334) namdb(1),n,(s(i+1),namdb(i+1),i=1,n), &
      (alogk0(l),l=1,ntemp),zz,wwt

      if (iprint==2 .and. myrank==0) write(iunit2,*) namdb(1),zz,wwt

      if (namdb(1) .eq. 'null') goto 185

!-----decide if reaction occurs in chosen system

      if (icase .eq. 1) then      ! convert to lower case
        call convrtlc (namdb(1),20)
      else if (icase .eq. 2) then ! convert to upper case
        call convrtuc (namdb(1),20)
      endif
 
      do mm = 1, nsrfmx
        if (namdb(1) .eq. namscx(mm)) then
          indxsrf(mm) = 1
          ns = mm
          goto 155
        endif
      enddo

      goto 145

  155 continue

      do 165 i = 2, n+1

        if (icase .eq. 1) then      ! convert to lower case
          call convrtlc (namdb(i),20)
        else if (icase .eq. 2) then ! convert to upper case
          call convrtuc (namdb(i),20)
        endif

!-------check for water
        if (namdb(i) .eq. 'H2O' .or. namdb(i) .eq. 'h2o') then
          goto 165
        endif

        do j = 1, ncpri
          if (namdb(i) .eq. namc(j)) then
            goto 165
          endif
        enddo

        do ix = 1, ncmplx
          if (namdb(i) .eq. namcx(ix)) then
            goto 165
          endif
        enddo

        do ngs = 1, ngas
          if (namdb(i) .eq. namg(ngs)) then
            goto 165
          endif
        enddo

        if (namdb(i) .eq. 'e-') then
          goto 165
        endif

        do nsrf = 1, nsrfmin
          do nsite = nsite1(nsrf), nsite2(nsrf)
            if (namdb(i) .eq. namsite(nsite)) then
              goto 165
            endif
          enddo
        enddo

!-------surface complex not found - skip reaction
        if (myrank == 0) &
        write(*,*) 'error reading database: surface complex not ', &
        'found: ',namdb(1),' STOP!'
        if (myrank == 0) &
        write(iunit2,*) 'error reading database: surface complex not ', &
        'found: ',namdb(1),' STOP!'
        stop

  165 continue

      if (iprint==2 .and. myrank==0) &
      write(iunit2,*) 'surface complex found: ',namdb(1)

      zsrf(ns) = zz
      wtsrf(ns) = wwt*1.d-3

      if (ntemp .gt. 1) then

        call fit(ndimmx,nbasis,ntemp,iflgint,alogk0,wmat,bvec,vec,indx)
          
        if (iflgint .eq. 1) then
          if (myrank == 0) &
          write(*,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e9.2)') namdb(1),(alogk0(l),l=1,ntemp)
          if (myrank == 0) &
          write(iunit2,'(" -> warning: LogK = 500 encountered for species: ",&
      &   a12,/,4x,1p10e9.2)') namdb(1),(alogk0(l),l=1,ntemp)
          iflgint = 2
        else if (iflgint .eq. 2) then
          iflgint = 0
        endif

        alogk(ns) = flogk(bvec,tempini+tkelvin)

        do j = 1, nbasis
          coef(ndim+ns,j) = bvec(j)
        enddo

      else if (ntemp .eq. 1) then
        alogk(ns) = alogk0(1)
      endif

!-----rearrange reactions - check for primary species
      do 170 i = 1, n
        do j = 1, ncpri
          if (iprint==2 .and. myrank==0) &
          write(iunit2,*) 'pri: ',namdb(i+1),namc(j)

          if (namdb(i+1) .eq. namc(j)) then
            ssorp(j,ns) = s(i+1)
            goto 170
          endif
        enddo
  170 continue

      do i = 1, n

!-------check for secondary species
        do k = 1, ncmplx
          if (namdb(i+1) .eq. namcx(k)) then
            do j = 1, ncpri
              ssorp(j,ns) = ssorp(j,ns) + shom(j,k)*s(i+1)
            enddo
            if (ntemp .gt. 1) then
              do l = 1, nbasis
                coef(ndim+ns,l) = coef(ndim+ns,l) + s(i+1)*coef(k,l)
              enddo
            endif
!...........note minus sign: eqhom = -eqhom as read in from database
            alogk(ns) = alogk(ns) - s(i+1)*eqhom(k) 
            goto 175
          endif
        enddo

!-------check for gaseous species
        do ngs = 1, ngas
          if (namdb(i+1) .eq. namg(ngs)) then
            do j = 1, ncpri
              ssorp(j,ns) = ssorp(j,ns) + sgas(j,ngs)*s(i+1)
            enddo
            if (ntemp .gt. 1) then
              do l = 1, nbasis
                coef(ndim+ns,l) = coef(ndim+ns,l) + &
                                 s(i+1)*coef(ncmplx+ngs,l)
              enddo
            endif
!...........note minus sign: eqgas = -eqgas as read in from database
            alogk(ns) = alogk(ns) - s(i+1)*eqgas(ngs)
            goto 175
          endif
        enddo

!-------check for electron
        if (namdb(i+1) .eq. 'e-') then
          ze(m) = s(i+1)
        endif

  175   continue
      enddo

!-----reverse sign to be consistent with surface complex on rhs of 
!       reaction
      eqsorp(ns) = -alogk(ns)

      goto 145

  185 continue

      if (myrank == 0) &
      write(*,*) '        --> finished reading surface complexation ', &
      'reactions'
      
      return

!-----------------------end of surface complexation reactions-----------
      
  334 if (myrank==0) write(*,*) 'error in reading database file: ',namdb(1),' stop'
      stop

  end subroutine dbsrf
  
!========================================================================

  subroutine dbtrf (namc,namsec,tempc,ntemp,nsec,nreac,nbasis)
  
  use ptran_global_module
  
  implicit none
  
! save
  
      character(len=namlen) :: namc(*)
      
      integer :: i,ii,iflag,j,k,l,ll,lp,m, &
      n,nbasis,ntemp, &
      nreac,nsec

  real*8 epot0,tkl,fac,vbargas

      real*8 :: alogk0(ntmpmx),tempc(ntmpmx)
      real*8 skpri0(ncxmx0,nxkmx)

  character(len=1) :: q
  character(len=namlen) :: varnam
  character(len=namlen) :: namsec(ndimmx)

  q = ' '
!=======================================================================
!-----transform kinetic aqueous reactions: eliminate secondary species
!     By construction skpri(j,k) does not include the kth kinetic 
!     species.
!=======================================================================
      do k = 1, ncxkin
        do lp = nparcxk1(k),nparcxk2(k)

        do i = 1, ncmplx
          if (skpri0(i,lp) .ne. zero) then
            do j = 1, ncpri
              skpri(j,lp) = skpri(j,lp) + shom(j,i)*skpri0(i,lp)
            enddo
!...........note plus sign because all constants now with reversed sign
            eqcxk(lp) = eqcxk(lp) + skpri0(i,lp)*eqhom(i)
          endif
        enddo
        do i = 1, ngas
          if (skpri0(ncmplx+i,lp) .ne. zero) then
            do j = 1, ncpri
              skpri(j,lp) = skpri(j,lp) + sgas(j,i)*skpri0(ncmplx+i,lp)
            enddo
!...........note plus sign because all constants now with reversed sign
            eqcxk(lp) = eqcxk(lp) + skpri0(ncmplx+i,lp)*eqgas(i)
          endif
        enddo

!-------add kinetic primary species representing reaction
        do j = 1, ncpri
          if (namc(j) .eq. namcxk(k)) then
!         if (namc(j) .eq. namrxnaq(lp)) then
            skpri(j,lp) = skpri(j,lp) - one
            goto 111
          endif
        enddo
  111   continue

!-------reverse sign (move all species to rhs: 0 = sum_j skpri(j,k) A_j
        do j = 1, ncpri
          skpri(j,lp) = -skpri(j,lp)
        enddo

        enddo
      enddo
!=======================================================================

!-----write out temperature interpolation coefficients
      if (myrank == 0) then
      if (ntemp>1 .and. (nsec>0 .or. mnrl>0)) then
        write(iunit2,'(/," temperature interpolation coefficients")')
        write(iunit2,7776)
 7776   format(' ',15x,'  ln(T+Tk)        1         T+Tk      (T+Tk)^-1 &
     &  (T+Tk)^-2')
        do i = 1, nsec
          write(iunit2,7777) namsec(i),(coef(i,j),j=1,nbasis)
        enddo
        do i = 1, mnrl
          write(iunit2,7777) namrl(i),(coef(nreac+i,j),j=1,nbasis)
        enddo
 7777   format(' ',a20,5(1pg13.6))
      endif

!-----write out chemical reactions

      write(iunit2,'(/,"-----primary species-----")')
      do j = 1, ncpri
        varnam = q//namc(j)(1:len_trim(namc(j)))//q
        write(iunit2,'(a12,1p3g12.4)') varnam,a0(j),z(j), &
!       write(iunit2,'(a1,a12,a1,1p3g12.4)') q,namc(j),q,a0(j),z(j),
!       write(iunit2,'(a12,1p3g12.4)') namc(j),a0(j),z(j),
        wt(j)*1.d3
      enddo

      if (ncmplx .gt. 0) then
        write(iunit2,'(/)')
        write(iunit2,*) '-----reversible homogeneous reactions-----'
        do i = 1, ncmplx
          n = 0
          do j = 1, ncpri
            if (shom(j,i) .ne. 0) n = n + 1
          enddo
          varnam = q//namcx(i)(1:len_trim(namcx(i)))//q
          write(iunit2,'(/,a20,i4,$)') varnam,n
          do j = 1, ncpri
            if (shom(j,i) .ne. 0) then
              varnam = q//namc(j)(1:len_trim(namc(j)))//q
              write(iunit2,'(f7.3,2x,a10,$)') shom(j,i),varnam
            endif
          enddo
          if (ntemp.eq.1) then
            write(iunit2,'(1p4g12.4)') -eqhom(i),ax0(i),zx(i),wtx(i)*1.d3
          else
            do l = 1, ntemp
              tkl = tempc(l) + tkelvin
              alogk0(l) = -(coef(i,1)*log(tkl) &
                          + coef(i,2)          &
                          + coef(i,3)*tkl      &
                          + coef(i,4)/tkl      &
                          + coef(i,5)/(tkl*tkl))
            enddo
            write(iunit2,'(1p12g12.4)') (-alogk0(l),l=1,ntemp), &
            ax0(i),zx(i),wtx(i)*1.d3
          endif
        enddo
      endif

      if (ncxkin .gt. 0) then
        write(iunit2,'(/,"-----irreversible homogeneous reactions-----")')
        do i = 1, ncxkin
        do lp = nparcxk1(i),nparcxk2(i)
          n = 0
          do j = 1, ncpri
            if (skpri(j,lp) .ne. zero) n = n + 1
          enddo
!         varnam = q//namcxk(i)(1:len_trim(namcxk(i)))//q
          varnam = q//namrxnaq(lp)(1:len_trim(namrxnaq(lp)))//q
          write(iunit2,'(/,a20,i4,$)') varnam,n
          do j = 1, ncpri
            if (skpri(j,lp) .ne. zero) then
              varnam = q//namc(j)(1:len_trim(namc(j)))//q
              write(iunit2,'(f7.3,2x,a10,$)') skpri(j,lp),varnam
            endif
          enddo
          if (zeaq(i) .ne. zero) then
            varnam = q//'e-'//q
            write(iunit2,'(f7.3,2x,a10,$)') zeaq(i),varnam
          endif
          if (ntemp.eq.1) then
            write(iunit2,'(1p4g12.4)') -eqcxk(lp),axk0(i),zxk(i), &
            wtxk(i)*1.d3
          else
            write(iunit2,'(1p4g12.4)') -eqcxk(lp),axk0(i),zxk(i), &
            wtxk(i)*1.d3
          endif
        enddo
        enddo
      endif
        
      if (mnrl.gt.0) then
        write(iunit2,'(/,"-----mineral reactions-----")')
        do m = 1, mnrl
          n = 0
          do j = 1, ncpri
            if (smnrl(j,m) .ne. zero) n = n + 1
          enddo
          varnam = q//namrl(m)(1:len_trim(namrl(m)))//q
          write(iunit2,'(/,a20,1pg12.4,i4,$)') varnam,vbar(m)*1.e3,n
          do j = 1, ncpri
            if (smnrl(j,m).ne.0) then
              varnam = q//namc(j)(1:len_trim(namc(j)))//q
!             write(iunit2,'(f7.3,2x,a10,$)') smnrl(j,m),varnam
              write(iunit2,'(1pg12.4,2x,a10,$)') smnrl(j,m),varnam
            endif
          enddo
          if (ze(m) .ne. zero) then
            varnam = q//'e-'//q
            write(iunit2,'(f7.3,2x,a10,$)') ze(m),varnam
          endif
          if (ntemp.eq.1) then
            write(iunit2,'(1p2g12.4)') -alnk(m),wtmin(m)*1.d3
          else
!           i = ncmplx+ncxkin+ngas+m
            i = ncmplx+ngas+m
            do l = 1, ntemp
              tkl = tempc(l) + tkelvin
              alogk0(l) = -(coef(i,1)*log(tkl) &
                          + coef(i,2)          &
                          + coef(i,3)*tkl      &
                          + coef(i,4)/tkl      &
                          + coef(i,5)/(tkl*tkl))
            enddo
            write(iunit2,'(1p10g12.4)') (-alogk0(l),l=1,ntemp), &
            wtmin(m)*1.d3
          endif

          if (vbar(m).eq.zero) then
            write(*,*) 'warning! -zero molar volume for mineral: ', &
            namrl(m)
            write(iunit2,*) 'warning! -zero molar volume for mineral: ', &
            namrl(m)
          else if (vbar(m).eq.0.5d0) then
            write(*,*) 'warning! -value 500 molar volume for mineral: ', &
            namrl(m)
            write(iunit2,*) 'warning! -value 500 molar volume for ', &
            'mineral: ',namrl(m)
          endif
        enddo
      endif

      if (nsrfmin.gt.0) then
        write(iunit2,'(/,"-----surface complexation reactions-----")')
        do m = 1, nsrfmin
          do ll = nsite1(m), nsite2(m)
            do i = nsorp1(ll), nsorp2(ll)
              n = 1
              do j = 1, ncpri
                if (ssorp(j,i).ne.zero) n = n + 1
              enddo
              varnam = q//namscx(i)(1:len_trim(namscx(i)))//q
              write(iunit2,'(/,a20,i4,$)') varnam,n
              varnam = q//namsite(ll)(1:len_trim(namsite(ll)))//q
              write(iunit2,'(" 1.0",2x,a20,$)') varnam
              do j = 1, ncpri
                if (ssorp(j,i).ne.zero) then
                  varnam = q//namc(j)(1:len_trim(namc(j)))//q
!                 write(iunit2,'(f7.3,2x,a10,$)') ssorp(j,i),varnam
                  write(iunit2,'(1pg12.4,2x,a10,$)') ssorp(j,i),varnam
                endif
              enddo
              if (ze(i) .ne. zero) then
                varnam = q//'e-'//q
                write(iunit2,'(f7.3,2x,a10,$)') ze(i),varnam
              endif
              if (ntemp.eq.1) then
                write(iunit2,'(1p3g12.4)') -eqsorp(i),zsrf(i), &
                wtsrf(i)*1.d3
              else
!               ii = ncmplx+ncxkin+ngas+mnrl+i
                ii = ncmplx+ngas+mnrl+i
                do l = 1, ntemp
                  tkl = tempc(l) + tkelvin
                  alogk0(l) = -(coef(ii,1)*log(tkl) &
                              + coef(ii,2)          &
                              + coef(ii,3)*tkl      &
                              + coef(ii,4)/tkl      &
                              + coef(ii,5)/(tkl*tkl))
                enddo
                write(iunit2,'(1p10g12.4)') (-alogk0(l),l=1,ntemp), &
                zsrf(i),wtsrf(i)*1.d3
              endif
            enddo
          enddo
        enddo
      endif

      if (ngas .gt. 0) then
        vbargas = zero
        write(iunit2, '(/,"-----gas reactions-----")')
        do i = 1, ngas
          n = 0
          do j = 1, ncpri
            if (sgas(j,i) .ne. 0) n = n + 1
          enddo
          varnam = q//namg(i)(1:len_trim(namg(i)))//q
          write(iunit2,'(/,a20,1pg12.4,i4,$)') varnam,vbargas,n
          do j = 1, ncpri
            if (sgas(j,i).ne.0) then
              varnam = q//namc(j)(1:len_trim(namc(j)))//q
              write(iunit2,'(f7.3,2x,a10,$)') sgas(j,i),varnam
            endif
          enddo
          if (ntemp.eq.1) then
            write(iunit2,'(1p2g12.4)') -eqgas(i),wtgas(i)*1.d3
          else
!           ii = ncmplx+ncxkin+i
            ii = ncmplx+i
            do l = 1, ntemp
              tkl = tempc(l) + tkelvin
              alogk0(l) = -(coef(ii,1)*log(tkl) &
                          + coef(ii,2)          &
                          + coef(ii,3)*tkl      &
                          + coef(ii,4)/tkl      &
                          + coef(ii,5)/(tkl*tkl))
            enddo
            write(iunit2,'(1p10g12.4)') (-alogk0(l),l=1,ntemp), &
            wtgas(i)*1.d3
          endif
        enddo
      endif
        
  599 format(' ',26x,12a6)
  600 format(' ',a12,g12.4,(15f6.2))

!-----write out electrochemical standard-state potentials
      iflag = 0
      do m = 1, mnrl
        if (ze(m).ne.0) then
          iflag = 1
          if (iflag.eq.1) then
            write(iunit2,'(/,"name",10x,"logk",8x,"Eh0",9x,"ze")')
          endif
          fac = (ze(m)*faraday)/(rgasj*(tempini+tkelvin))/log(10.d0)
          epot0 = alnk(m)/fac
          write(iunit2,'(a12,1p3e12.4)') namrl(m),alnk(m),epot0,ze(m)
        endif
      enddo
      do i = 1, ncxkin
        if (zeaq(i).ne.0) then
          if (iflag.eq.0) iflag=2
          if (iflag.eq.2) then
            write(iunit2,'(/,"name",10x,"logk",8x,"Eh0",9x,"ze")')
          endif
          fac = (zeaq(i)*faraday)/(rgasj*(tempini+tkelvin))/log(10.d0)
          epot0 = eqcxk(i)/fac
          write(iunit2,'(a12,1p3e12.4)') namcxk(i),eqcxk(i),epot0,zeaq(i)
        endif
      enddo
      
      endif

      return
  end subroutine dbtrf
  
!========================================================================

      subroutine inverse (n0,ntemp,nbasis,nsec,namsec,namvec,alogk, &
                          wmat,bmat,amat,coef0,namc)
      
      use ptran_global_module
      
      implicit none
      
      integer :: n0,nsec,i,ii,iflg,j,k,l,ngss,ngs,ntemp,nbasis
      real*8 :: sum,alogsum,dd
      integer :: indx(ndimmx)
      
  real*8 :: alogk(ndimmx),amat(ndimmx,ndimmx),         &
                           bmat(ndimmx,ncmx+nxkmx),smat(ndimmx,ndimmx),  &
                           ainv(ndimmx,ndimmx),y(ndimmx), &
                           wmat(ndimmx,ndimmx),alogkp(ndimmx), &
                           bb(ndimmx),coef0(ndimmx,5)
  character(len=namlen) :: namsec(ndimmx),namvec(ndimmx)
  character(len=namlen) :: namc(*)

!-----compute inverse matrix to [a]

        do k = 1, nsec
          do l = 1, nsec
            wmat(l,k) = amat(l,k)
          enddo
        enddo
        if (iprint>=2 .and. myrank==0) then
          write(iunit2,'(12x,25(a5,1x))') (namvec(l),l=1,nsec)
          do k = 1, nsec
            write(iunit2,9999) namsec(k),(wmat(k,l),l=1,nsec)
          enddo
          write(iunit2,'(12x,25(a5,1x))') (namc(j),j=1,ncpri)
          do k = 1, nsec
            write(iunit2,9999) namsec(k),(bmat(k,j),j=1,ncpri)
          enddo
 9999     format(' ',a8,(' ',25f6.2))
        endif

        call ludcmp(wmat,nsec,n0,indx,dd)

        do i = 1, nsec
          do k = 1, nsec
            y(k) = zero
            if (k .eq. i) y(k) = one
            bb(k) = y(k)
          enddo

          call lubksb(wmat,nsec,n0,indx,y)

          call mprove(amat,wmat,nsec,n0,indx,bb,y)
          call mprove(amat,wmat,nsec,n0,indx,bb,y)

          do l = 1, nsec
            ainv(l,i) = y(l)
          enddo
          
!         write(*,*) 'ptrandbase: ',i,(ainv(l,i),l=1,nsec)
        enddo

!-------check inverse
        iflg = 0
        do i = 1, nsec
          do j = 1, nsec
            sum = zero
            do k = 1, nsec
              sum = sum + ainv(i,k)*amat(k,j)
            enddo
            if (i.eq.j) then
              if (sum .lt. 1.d0-eps .or. sum .gt. 1.d0+eps) then
                if (myrank==0) write(*,*) 'identity check failed',i,j,sum
                iflg = 1
!               stop
              endif
            else if (i.ne.j) then
              if (abs(sum) .gt. eps) then
                if (myrank==0) write(*,*) 'identity check failed',i,j,sum
                iflg = 1
!               stop
              endif
            endif
            smat(j,i) = sum
          enddo
        enddo

        if (iflg .eq. 1) then
          if (myrank == 0) then
            write(iunit2,105)
            do i = 1, nsec
              write(iunit2,107) i,(smat(l,i), l = 1, nsec)
            enddo
          endif
 105      format(/,' identity check')
 107      format(' ',i3,10(1pg13.6))
          stop
        endif

!-------compute new reaction coefficients and logK's
        do i = 1, nsec
          do j = 1, ncpri
            sum = zero
            do k = 1, nsec
              sum = sum + ainv(i,k)*bmat(k,j)
            enddo
            if (abs(sum) .gt. 1.d-7) then
              smat(j,i) = sum
            else
              smat(j,i) = zero
            endif
          enddo
          alogsum = zero
          do ii = 1, nsec
            alogsum = alogsum + ainv(i,ii)*alogk(ii)
            if (ntemp .gt. 1) then
              do j = 1, nbasis
                coef(i,j) = coef(i,j) + ainv(i,ii)*coef0(ii,j)
              enddo
            endif
          enddo
          alogkp(i) = alogsum
          
!         write(*,*) 'ptrandbase: ',i,alogkp(i),(coef0(i,j),j=1,nbasis)
        enddo

!-------construct sub-matrix shom
        do j = 1, ncpri
          do i = 1, ncmplx
            shom(j,i) = smat(j,i)
          enddo
          do ngs = 1, ngas
            ngss = ncmplx + ngs
            sgas(j,ngs) = smat(j,ngss)
!           print *,'ptrandbase: ',j,ngs,ncpri,ncmplx,ngss,sgas(j,ngs)
          enddo
        enddo

!-------store logK's for aqueous complexes and gases
!             note minus sign to convert from dissociation to
!             association reaction
        do i = 1, ncmplx
          eqhom(i) = -alogkp(i)
        enddo
        do ngs = 1, ngas
          ngss = ncmplx + ngs
          eqgas(ngs) = -alogkp(ngss)
        enddo

  end subroutine inverse
  
!======================================================================

      subroutine fit (n0,nbasis,ntemp,iflgint,alogk0,wmat,bvec,vec,indx)

      use ptran_global_module

      implicit none

      integer :: i,j,k,int(1000),indx(*),nbasis,ntemp,iflgint,n0
      real*8 :: alogk0(*),wmat(n0,*),vec(5,*),bvec(5),dd

      iflgint = 0
      do j = 1, nbasis
        bvec(j) = zero
        do i = 1, ntemp
          if (alogk0(i) .ne. 500.) then
            bvec(j) = bvec(j) + alogk0(i)*vec(j,i)
            int(i) = 1
          else if (alogk0(i) .eq. 500.) then
            iflgint = 1
            int(i) = 0
          else if (alogk0(i) .gt. 500.) then
            if (myrank==0) &
            write(*,*) 'error in fit: log K .gt. 500---stop!'
            stop
          endif
        enddo
      enddo
   
      do j = 1, nbasis
        do k = j, nbasis
          wmat(j,k) = zero
          do i = 1, ntemp
            if (int(i) .eq. 1) then
              wmat(j,k) = wmat(j,k) + vec(j,i)*vec(k,i)
            endif
          enddo
          if (j .ne. k) wmat(k,j) = wmat(j,k)
        enddo
      enddo
      call ludcmp(wmat,nbasis,n0,indx,dd)
      call lubksb(wmat,nbasis,n0,indx,bvec)

      return
      end subroutine fit
      
!========================================================================

      function flogk(b,temp)

      implicit none

      real*8 temp,b(5),flogk

      flogk = b(1)*log(temp) &
            + b(2)           &
            + b(3)*temp      &
            + b(4)/temp      &
            + b(5)/(temp*temp)
     
      return
      end function flogk
      
!========================================================================

      subroutine ludcmp(a,n,np,indx,d)

      implicit none
      
      integer, parameter :: nmax=1000
      
      integer :: i,imax,n,indx(n),j,k,np
      
      real*8 :: a(np,np),vv(nmax),d,aamax,tiny=1.d-20,sum,dum

      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix. [ludcmp]'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0) a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.d0) a(n,n)=tiny
      return
      end subroutine ludcmp
      
!========================================================================

      subroutine lubksb(a,n,np,indx,b)

      implicit none

      integer :: i,ii,ll,j,np,n,indx(n)
      real*8 :: a(np,np),b(n),sum

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end subroutine lubksb
      
!========================================================================

      subroutine mprove (a,alud,n,np,indx,b,x)

      implicit none

      integer, parameter :: nmax=1000

      integer :: i,n,indx(n),j,np
      real*8 :: a(np,np),alud(np,np),b(n),x(n),r(nmax)
      real*8 :: sdp
      
      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+dble(a(i,j))*dble(x(j))
11      continue
        r(i)=sdp
12    continue
      call lubksb(alud,n,np,indx,r)
      do 13 i=1,n
        x(i)=x(i)-r(i)
13    continue
      return
      end subroutine mprove

!=======================================================================

      subroutine convrtlc (char,nchar)

      implicit none

!=======================================================================
!     Convert character string to lower case letters   
!     'char' of length nchar. Nchar is currently set to a maximum of 30, 
!     which can be changed to a larger value if needed.
!=======================================================================
      integer :: i, m, nchar
      
      save lowercase,uppercase
      character*26 lowercase,uppercase
      character(len=*) char

      data lowercase,uppercase &
        /'abcdefghijklmnopqrstuvwxyz', &
         'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do i = 1, nchar
        do m = 1, 26
          if(char(i:i) .eq. uppercase(m:m)) then
            char(i:i) = lowercase(m:m)
            goto 1
          endif
        enddo
    1   continue
      enddo

      return
      end subroutine convrtlc

!=======================================================================

      subroutine convrtuc (char,nchar)

      implicit none
      
!=======================================================================
!     Convert character string to upper case letters   
!     'char' of length nchar. Nchar is currently set to a maximum of 30, 
!     which can be changed to a larger value if needed.
!=======================================================================
      integer :: i, m, nchar
      
      save lowercase,uppercase
      character*26 lowercase,uppercase
      character(len=*) char

      data lowercase,uppercase &
        /'abcdefghijklmnopqrstuvwxyz', &
         'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do i = 1, nchar
        do m = 1, 26
          if(char(i:i) .eq. lowercase(m:m)) then
            char(i:i) = uppercase(m:m)
            goto 1
          endif
        enddo
    1   continue
      enddo

      return
      end subroutine convrtuc

end module ptran_dbase_module
