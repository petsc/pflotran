!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! VERSION/REVISION HISTORY
 
! $Id: ptran_read.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
! $Log: ptran_read.F90,v $
! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
! initial import
!
! Revision 1.3  2004/04/06 17:38:32  lichtner
! Added time unit and conversion.
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

module ptran_read_module

  public

  interface ptran_read
    module procedure ptran_read
  end interface ptran_read

contains

  subroutine ptran_read
  
! keywords: GRID, PROC, OPTS, TOLR, SYST, DXYZ, FLOW, DIFF, DBAS, COMP,
!           BCON, SOUR, SOLV, BRK, TIME, FELD, AQCX, GAS, MNRL, MNIR,
!           IONX, SORP, COUP

  use ptran_global_module
  use fileio_module
  use trdynmem_module
  use ptran_dbase_module

  implicit none

  character(len=strlen) :: string
  character(len=wordlen) :: word, strtim
  character(len=cardlen) :: card
  character(len=namlen) :: name

  integer :: i,ibrk,ierr,iireg,ilog,ir,isrc,ival,j,k,l,lp,m,n,nr
  integer :: n1,ns,ns1,n2,ns2,nss,nsite
  
  real*8 :: fac,val,cec00,eqiex0,area00,cover0,eqsorp00,sited00
! real*8 :: rkfsrf0,rkbsrf0
  
  ierr = 0
  
!-----read input data
  if(myrank==0) write(*,*) '--> read input file: ','ptran.in', &
  ' units(r/w)= ',iunit1,iunit2
  open(iunit1,file='ptran.in',action="read",status='old')
  if(myrank==0) &
  open(iunit2,file='ptran.out',action="write",status='unknown')
  
  do
      
    call fiReadFlotranString(iunit1, string, ierr)

    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)
    
    if (myrank == 0) print *, card
    
    select case (card)
    
!....................

    case ('GRID')
    
      call fiReadStringErrorMsg('GRID',ierr)

      call fiReadInt(string,nx,ierr)
      call fiDefaultMsg('nx',ierr)
      call fiReadInt(string,ny,ierr)
      call fiDefaultMsg('ny',ierr)
      call fiReadInt(string,nz,ierr)
      call fiDefaultMsg('nz',ierr)

      nmax = nx*ny*nz
      nxy = nx*ny
      
      if (myrank == 0) &
      write(iunit2,'(" *GRID",/, &
    & "  nx   = ",3x,i4,/, &
    & "  ny   = ",3x,i4,/, &
    & "  nz   = ",3x,i4,/, &
    & "  nmax = ",i7)') nx,ny,nz,nmax
    
!....................

    case ('PROC')

      call fiReadStringErrorMsg('PROC',ierr)

      call fiReadInt(string,npx,ierr)
      call fiDefaultMsg('npx',ierr)
      call fiReadInt(string,npy,ierr)
      call fiDefaultMsg('npy',ierr)
      call fiReadInt(string,npz,ierr)
      call fiDefaultMsg('npz',ierr)

      if (myrank == 0) &
      write(iunit2,'(/," *PROC",/, &
    & "  npx   = ",3x,i4,/, &
    & "  npy   = ",3x,i4,/, &
    & "  npz   = ",3x,i4)') npx,npy,npz
  
      if (commsize .ne. npx*npy*npz) then
        if (myrank==0) write(*,*) 'Incorrect number of processors specified: ', &
        npx*npy*npz,' commsize = ',commsize
        stop
      endif

!....................

    case ('OPTS')
    
      call fiReadStringErrorMsg('OPTS',ierr)

      call fiReadInt(string,iact,ierr)
      call fiDefaultMsg('iact',ierr)
      
      call fiReadInt(string,loglin,ierr)
      call fiDefaultMsg('loglin',ierr)

      call fiReadInt(string,iaccel,ierr)
      call fiDefaultMsg('iaccel',ierr)

      call fiReadInt(string,icomprs,ierr)
      call fiDefaultMsg('icomprs',ierr)

      call fiReadInt(string,iblkfmt,ierr)
      call fiDefaultMsg('iblkfmt',ierr)

      call fiReadInt(string,isurf,ierr)
      call fiDefaultMsg('isurf',ierr)

      if (myrank==0) &
      write(iunit2,'(/," *OPTS ",/, &
    & "  iact    = ",i2,/, &
    & "  loglin  = ",i2,/, &
    & "  iaccel  = ",i2,/, &
    & "  icomprs = ",i2,/, &
    & "  iblkfmt = ",i2,/, &
    & "  isurf   = ",i2)') &
      iact,loglin,iaccel,icomprs,iblkfmt,isurf

!....................

    case ('TOLR')
    
      call fiReadStringErrorMsg('TOLR',ierr)

      call fiReadInt(string,kmax,ierr)
      call fiDefaultMsg('kmax',ierr)

      call fiReadInt(string,newton_max,ierr)
      call fiDefaultMsg('newton_max',ierr)
      
      call fiReadInt(string,icut_max,ierr)
      call fiDefaultMsg('icut_max',ierr)
      
      call fiReadDouble(string,tolexp,ierr)
      call fiDefaultMsg('tolexp',ierr)
      
      call fiReadInt(string,iwarn,ierr)
      call fiDefaultMsg('iwarn',ierr)
      
      call fiReadInt(string,iprint,ierr)
      call fiDefaultMsg('iprint',ierr)
      
      if (myrank==0) write(iunit2,'(/," *TOLR ",/, &
     &"  kmax       = ",i6,/, &
     &"  newton_max = ",3x,i3,/, &
     &"  icut_max   = ",3x,i3,/, &
     &"  tolexp     = ",3x,1pe12.4,/, &
     &"  iwarn      = ",3x,i3,/, &
     &"  iprint     = ",3x,i3 &
     &   )') &
      kmax,newton_max,icut_max,tolexp,iwarn,iprint

!....................

    case ('SYST')
    
      call fiReadStringErrorMsg('SYST',ierr)

      call fiReadInt(string,iphase,ierr)
      call fiDefaultMsg('iphase',ierr)

      call fiReadInt(string,isothrm,ierr)
      call fiDefaultMsg('isothrm',ierr)
      
      call fiReadInt(string,iread_vel,ierr)
      call fiDefaultMsg('iread_vel',ierr)
      
      call fiReadNChars(string,name,namlen,.true.,ierr)
      call fiErrorMsg('nam','SYST',ierr)
      flowvx = name

      call fiReadNChars(string,name,namlen,.true.,ierr)
      call fiErrorMsg('nam','SYST',ierr)
      flowvy = name

      call fiReadNChars(string,name,namlen,.true.,ierr)
      call fiErrorMsg('nam','SYST',ierr)
      flowvz = name
      
      call fiReadInt(string,iread_sat,ierr)
      call fiDefaultMsg('iread_sat',ierr)
      
      call fiReadNChars(string,name,namlen,.true.,ierr)
      call fiErrorMsg('nam','SYST',ierr)
      fsat = name
      
      if (myrank==0) then
        write(iunit2,'(/," *SYST ",/, &
     &  "  iphase    = ",i6,/, &
     &  "  isothrm   = ",3x,i3,/, &
     &  "  iread_vel = ",3x,i3,/, &
     &  "  iread_sat = ",3x,i3 &
     &  )') iphase,isothrm,iread_vel,iread_sat
        write(iunit2,'(4x,"file_vx:  ",a20)') flowvx
        write(iunit2,'(4x,"file_vy:  ",a20)') flowvy
        write(iunit2,'(4x,"file_vz:  ",a20)') flowvz
        write(iunit2,'(4x,"file_sat: ",a20)') fsat
      endif

!....................

    case ('COUP')
    
      call fiReadStringErrorMsg('COUP',ierr)

      call fiReadInt(string,ipor,ierr)
      call fiDefaultMsg('ipor',ierr)

      call fiReadDouble(string,tolpor,ierr)
      call fiDefaultMsg('tolpor',ierr)
      
      if (myrank==0) then
        write(iunit2,'(/," *COUP ",/, &
     &  "  ipor    = ",i6,/, &
     &  "  tolpor   = ",1pe12.4 &
     &  )') ipor,tolpor
      endif
      
!....................

    case ('DXYZ')

      allocate(dx0(nx))
      allocate(dy0(ny))
      allocate(dz0(nz))

      call readxyz_ptran (dx0,nx,iunit1)
      call readxyz_ptran (dy0,ny,iunit1)
      call readxyz_ptran (dz0,nz,iunit1)

      if (myrank==0) then
        write(IUNIT2,'(/," *DXYZ ")')
        write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (dx0(i),i=1,nx)
        write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (dy0(i),i=1,ny)
        write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (dz0(i),i=1,nz)
      endif
      
!....................

    case ('FLOW')
    
      call fiReadStringErrorMsg('FLOW',ierr)

      call fiReadDouble(string,vlx0,ierr)
      call fiDefaultMsg('vlx',ierr)
      call fiReadDouble(string,vly0,ierr)
      call fiDefaultMsg('vly',ierr)
      call fiReadDouble(string,vlz0,ierr)
      call fiDefaultMsg('vlz',ierr)
      call fiReadDouble(string,vgx0,ierr)
      call fiDefaultMsg('vgx',ierr)
      call fiReadDouble(string,vgy0,ierr)
      call fiDefaultMsg('vgy',ierr)
      call fiReadDouble(string,vgz0,ierr)
      call fiDefaultMsg('vgz',ierr)
      
      if (myrank==0) then
        write(iunit2,'(/," *FLOW ",/ &
    &   "  vl(x, y, z) = ",1p3e12.4," [m/y]")') vlx0,vly0,vlz0
        write(iunit2,'("  vg(x, y, z) = ",1p3e12.4," [m/y]")') &
        vgx0,vgy0,vgz0
      endif      
!....................

    case ('DIFF')
    
      call fiReadStringErrorMsg('DIFF',ierr)
      
      call fiReadDouble(string,difaq,ierr)
      call fiDefaultMsg('difaq',ierr)
      
      call fiReadDouble(string,delhaq,ierr)
      call fiDefaultMsg('delhaq',ierr)
      
      call fiReadDouble(string,difgas,ierr)
      call fiDefaultMsg('difgas',ierr)
      
      call fiReadDouble(string,dgexp,ierr)
      call fiDefaultMsg('dgexp',ierr)
      
      call fiReadDouble(string,por0,ierr)
      call fiDefaultMsg('por0',ierr)

      call fiReadDouble(string,sat0,ierr)
      call fiDefaultMsg('sat0',ierr)

      call fiReadDouble(string,temp0,ierr)
      call fiDefaultMsg('temp0',ierr)

      call fiReadDouble(string,pref0,ierr)
      call fiDefaultMsg('pref0',ierr)

      if (myrank==0) write(iunit2,'(/," *DIFF ",/, &
    & "  difaq  = ",1pe12.4,/, &
    & "  delhq  = ",1pe12.4,/, &
    & "  difgas = ",1pe12.4,/, &
    & "  dgexp  = ",1pe12.4,/, &
    & "  por0   = ",1pe12.4,/, &
    & "  sat0   = ",1pe12.4,/, &
    & "  temp0  = ",1pe12.4,/, &
    & "  pref0  = ",1pe12.4)') difaq,delhaq,difgas,dgexp, &
                              por0,sat0,temp0,pref0

      tempini = temp0

!....................

    case ('DBAS')

      call fiReadFlotranString(iunit1,string,ierr)
      call fiReadStringErrorMsg('DBAS',ierr)

      call fiReadNChars(string,dbpath,strlen,.true.,ierr)
      call fiErrorMsg('dbpath','DBAS',ierr)
      if (myrank == 0) &
      write(iunit2,1050) (dbpath(i:i),i=1,len_trim(dbpath))

1050 format(/,' *DBAS database: ',160a1)
1060 format(160a1)

!....................

    case ('COMP')

      call fiReadStringErrorMsg('COMP',ierr)

      allocate(nreg_val(nmax))
      ireg = 0
      do ! loop over data block: region + composition
        iloop = 0
        ireg = ireg + 1
        do ! loop over regions
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit

          iloop = iloop + 1

          call fiReadInt(string,i1,ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,i2,ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,j1,ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,j2,ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,k1,ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,k2,ierr)
          call fiDefaultMsg('k2',ierr)
          
          do k=k1,k2
            do j=j1,j2
              do i=i1,i2
                n = i+(j-1)*nx+(k-1)*nxy
                nreg_val(n) = ireg
              enddo
            enddo
          enddo
        enddo ! end loop over regions

        if (iloop == 0) exit

        ncomp = 0
        do ! loop over primary species
          
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          
          call fiReadNChars(string,name,namlen,.true.,ierr)
          call fiErrorMsg('nam','COMP',ierr)

          ncomp = ncomp + 1
          
          if (ireg > 1) then
            j = ncomp
            if (nam(j) .ne. name) then
              if (myrank == 0) &
              write(*,*) 'COMP: error primary species out of order-stop!', &
              j,nam(j),' ',name,' block= ',ireg
              stop
            endif
          endif

          if (ncomp > ncmx) then
            if (myrank == 0) &
            write(*,*) 'too many primary species: stop! ',name,ncomp,ncmx
            stop
          endif
          nam(ncomp) = name

          call fiReadInt(string,ival,ierr)
          call fiDefaultMsg('itype',ierr)
          itype(ncomp,ireg) = ival
          
          call fiReadDouble(string,val,ierr)
          call fiDefaultMsg('ctot',ierr)
          ctot(ncomp,ireg) = val

          call fiReadNChars(string,name,namlen,.false.,ierr)
          call fiErrorMsg('ncon','COMP',ierr)
          ncon(ncomp,ireg) = name
        enddo ! end loop over components
      enddo ! end loop over data blocks

      initreg = ireg-1

      if (myrank == 0) then
        write(iunit2,'(/," *COMP: # regions = ",i3, &
     &  " # components = ",i3)') initreg,ncomp
        do iireg = 1, initreg
          write(iunit2,'(/,"region = ",i3)') iireg
          write(iunit2,*) 'species            itype   ctot        ', &
          'constraint'
          do j = 1, ncomp
            write(iunit2,'(a20,i4,1pe12.4,4x,a20)') nam(j), &
            itype(j,iireg),ctot(j,iireg),ncon(j,iireg)
          enddo
        enddo
      endif
      
!....................

!-----read in boundary conditions

    case ('BCON')

!-----------------------------------------------------------------------
!-----boundary conditions:	ibnd:	1-left,		2-right
!					3-top,		4-bottom
!					5-front,	6-back
!
!				ibndtyp:	1-concentration
!				ibndtyp:	2-flux
!				ibndtyp:	3-zero gradient
!-----------------------------------------------------------------------
      ibc = 0
      ir = 0
      ireg = initreg
      ibcreg(1) = initreg + 1
      iregbc1(1) = 1
      do ! loop over blocks
      
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('BCON',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        ibc = ibc + 1  ! RTM: Number of boundary conditions
        ireg = ireg + 1
        
        call fiReadInt(string,ibndtyp(ibc),ierr)
        call fiDefaultMsg('ibndtyp',ierr)

        call fiReadInt(string,iface(ibc),ierr)
        call fiDefaultMsg('iface',ierr)
      
        call fiReadDouble(string,tempbc(ibc),ierr)
        call fiDefaultMsg('tempbc',ierr)

        do ! loop over regions
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ir = ir + 1

          call fiReadInt(string,i1bc(ir),ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,i2bc(ir),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,j1bc(ir),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,j2bc(ir),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,k1bc(ir),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,k2bc(ir),ierr)
          call fiDefaultMsg('k2',ierr)    
        enddo ! end loop over regions
        iregbc2(ibc) = ir
        if (ibc+1 > nrgbcmx) then
          write(*,*) 'dimension of boundary regions too small-stop: ', &
          ibc+1,nrgbcmx
          stop
        else
          iregbc1(ibc+1) = iregbc2(ibc)+1
        endif

        do j = 1, ncomp ! loop over primary species
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
          
          call fiReadNChars(string,name,namlen,.true.,ierr)
          call fiErrorMsg('nam','BCON',ierr)

          if (nam(j) .ne. name) then
            if (myrank == 0) &
            write(*,*) 'BCON: error primary species out of order-stop!', &
            j,nam(j),' ',name,' block= ',ibc
            stop
          endif

          call fiReadInt(string,itype(j,ireg),ierr)
          call fiDefaultMsg('itype',ierr)
          
          call fiReadDouble(string,ctot(j,ireg),ierr)
          call fiDefaultMsg('ctot',ierr)
          
          call fiReadNChars(string,ncon(j,ireg),namlen,.false.,ierr)
          call fiErrorMsg('ncon','BCON',ierr)
        enddo
      enddo ! end loop over blocks
      nblkbc = ibc
      ibcreg(2) = ireg

      if (myrank == 0) then
        write(iunit2,'(/," *BCON: # BC Blocks = ",i3)') nblkbc
        iireg = initreg
        do ibc = 1, nblkbc
          write(iunit2,'(/," Block: ",i3)') ibc
          write(iunit2,*) 'bndtyp  face       temp[C]'
          write(iunit2,'(i4,4x,i4,4x,1pe12.4)') ibndtyp(ibc),iface(ibc), &
          tempbc(ibc)
          write(iunit2,*) '   i1    i2    j1    j2    k1    k2'
          do ir = iregbc1(ibc),iregbc2(ibc)
            write(iunit2,'(6i6)') i1bc(ir),i2bc(ir),j1bc(ir),j2bc(ir), &
            k1bc(ir),k2bc(ir)
          enddo
          write(iunit2,*) 'species            itypbc  ctotbc      ', &
          'constraint'
          iireg = iireg + 1
          do j = 1, ncomp
            write(iunit2,'(a20,i4,1pe12.4,4x,a20)') nam(j), &
            itype(j,iireg),ctot(j,iireg),ncon(j,iireg)
          enddo
        enddo
      endif
      
!....................

    case ('SOUR')
    
      call fiReadStringErrorMsg('SOUR',ierr)

      isrc = 0
      if (ibcreg(2) /= 0) then
        ibcreg(3) = ibcreg(2) + 1
      else
        ibcreg(3) = initreg + 1
        ireg = initreg
      endif
      
      do ! loop over sources
      
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('SOUR',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        isrc = isrc + 1

        ireg = ireg + 1
        
!       tempbc(ireg) = 25.d0
        
!       do ! loop over regions
!         call fiReadFlotranString(iunit1,string,ierr)
!         call fiReadStringErrorMsg('SOUR',ierr)
      
!         if (string(1:1) == '.' .or. string(1:1) == '/') exit
!         ireg = ireg + 1

          call fiReadInt(string,is1(isrc),ierr)
          call fiDefaultMsg('is1',ierr)
          call fiReadInt(string,is2(isrc),ierr)
          call fiDefaultMsg('is2',ierr)
          call fiReadInt(string,js1(isrc),ierr)
          call fiDefaultMsg('js1',ierr)
          call fiReadInt(string,js2(isrc),ierr)
          call fiDefaultMsg('js2',ierr)
          call fiReadInt(string,ks1(isrc),ierr)
          call fiDefaultMsg('ks1',ierr)
          call fiReadInt(string,ks2(isrc),ierr)
          call fiDefaultMsg('ks2',ierr)
          
!       enddo ! end loop over regions
        
        i = 0
        
        do ! loop over time intervals
      
          call fiReadFlotranString(iunit1,string,ierr)
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        
          i = i + 1
        
          call fiReadDouble(string,timesrc(i,isrc),ierr)
          call fiDefaultMsg('timesrc',ierr)

          call fiReadDouble(string,tempsrc(i,isrc),ierr)
          call fiDefaultMsg('tempsrc',ierr)
      
          call fiReadDouble(string,qsrc(i,isrc),ierr)
          call fiDefaultMsg('qsrc',ierr)
      
          call fiReadDouble(string,csrc(i,isrc),ierr)
          call fiDefaultMsg('csrc',ierr)
        
        enddo

        nsrc = i

        do j = 1, ncomp ! loop over primary species
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('SOUR',ierr)
          
          call fiReadNChars(string,name,namlen,.true.,ierr)
          call fiErrorMsg('nam','SOUR',ierr)

          if (nam(j) .ne. name) then
            if (myrank == 0) &
            write(*,*) 'SOUR: error species out of order-stop!', &
            j,nam(j),' ',name,' block= ',ireg
            stop
          endif

          call fiReadInt(string,itype(j,ireg),ierr)
          call fiDefaultMsg('itype',ierr)
          
          call fiReadDouble(string,ctot(j,ireg),ierr)
          call fiDefaultMsg('ctot',ierr)
          
          call fiReadNChars(string,ncon(j,ireg),namlen,.false.,ierr)
          ierr = 0
          call fiErrorMsg('ncon','SOUR',ierr)
        enddo
      enddo
      nblksrc = isrc
      ibcreg(4) = ibcreg(3) + isrc - 1

      if (myrank == 0) then
        write(iunit2,'(/," *SOURce: # sources = ",i3)') nblksrc
        iireg = ibcreg(3) - 1
        do isrc = 1, nblksrc
          iireg = iireg + 1
          tempbc(iireg) = tempsrc(1,isrc) ! source/sink must be isothermal
          write(iunit2,'(/,"source = ",i3)') isrc
          write(iunit2,'("  #     time[s]     temp[C]  QH2O[kg/s]  QCO2[kg/s]")')
          do i = 1, nsrc
            write(iunit2,'(i3,1x,1p10e12.4)') i, &
            timesrc(i,isrc),tempsrc(i,isrc),qsrc(i,isrc),csrc(i,isrc)
          enddo
          write(iunit2,*) 'species            itype   ctot        ', &
          'constraint'
          do j = 1, ncomp
            write(iunit2,'(a20,i4,1pe12.4,4x,a20)') nam(j), &
            itype(j,iireg),ctot(j,iireg),ncon(j,iireg)
          enddo
        enddo
      endif

!....................

    case ('SOLV')
    
      call fiReadStringErrorMsg('SOLV',ierr)

      call fiReadDouble(string,eps,ierr)
      call fiDefaultMsg('eps',ierr)

      call fiReadDouble(string,atol_petsc,ierr)
      call fiDefaultMsg('atol_petsc',ierr)

      call fiReadDouble(string,rtol_petsc,ierr)
      call fiDefaultMsg('rtol_petsc',ierr)

      call fiReadDouble(string,dtol_petsc,ierr)
      call fiDefaultMsg('dtol_petsc',ierr)
      
      call fiReadInt(string,maxits_petsc,ierr)
      call fiDefaultMsg('maxits_petsc',ierr)

      if (myrank==0) write(iunit2,'(/," *SOLV ",/, &
     &"  eps          = ",1pe12.4,/, &
     &"  atol_petsc   = ",1pe12.4,/, &
     &"  rtol_petsc   = ",1pe12.4,/, &
     &"  dtol_petsc   = ",1pe12.4,/, &
     &"  maxits_petsc = ",8x,i4)') &
      eps,atol_petsc,rtol_petsc,dtol_petsc,maxits_petsc

!....................
      
    case ('BRK')

      ibrk = 0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('BRK',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ibrk = ibrk + 1
      
        call fiReadInt(string,i1brk(ibrk),ierr) 
        call fiDefaultMsg('i1',ierr)
        call fiReadInt(string,i2brk(ibrk),ierr)
        call fiDefaultMsg('i2',ierr)
        call fiReadInt(string,j1brk(ibrk),ierr)
        call fiDefaultMsg('j1',ierr)
        call fiReadInt(string,j2brk(ibrk),ierr)
        call fiDefaultMsg('j2',ierr)
        call fiReadInt(string,k1brk(ibrk),ierr)
        call fiDefaultMsg('k1',ierr)
        call fiReadInt(string,k2brk(ibrk),ierr)
        call fiDefaultMsg('k2',ierr)

        call fiReadInt(string,ibrktyp(ibrk),ierr)
        call fiDefaultMsg('ibrktyp',ierr)

        call fiReadInt(string,ibrkface(ibrk),ierr)
        call fiDefaultMsg('ibrkface',ierr)
      enddo
      ibrkcrv = ibrk
            
      if (myrank==0) then
        write(IUNIT2,'(/," *BRK: ibrk = ",i4)') ibrkcrv
        write(IUNIT2,'("  i1  i2  j1  j2  k1  k2  ibrktyp  ibrkface  ")')
        do ibrk = 1, ibrkcrv
          write(IUNIT2,'(6i4,4x,i2,7x,i2)') &
          i1brk(ibrk),i2brk(ibrk), &
          j1brk(ibrk),j2brk(ibrk), &
          k1brk(ibrk),k2brk(ibrk),ibrktyp(ibrk), &
          ibrkface(ibrk)
        enddo
      endif

!....................

    case ('TIME')

      call fiReadStringErrorMsg('TIME',ierr)
      
      call fiReadWord(string,strtim,.false.,ierr)

      tunit = strtim
      
      if (tunit == 's') then
        tconv = 1.d0
      else if (tunit == 'm') then
        tconv = 60.d0
      else if (tunit == 'h') then
        tconv = 60.d0 * 60.d0
      else if (tunit == 'd') then
        tconv = 60.d0 * 60.d0 * 24.d0
      else if (tunit == 'mo') then
        tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
      else if (tunit == 'y') then
        tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
      else
        if (myrank == 0) then
          write(*,'(" Time unit: ",a3,/, &
   &      " Error: time units must be one of ",/, &
   &      "   s -seconds",/,"   m -minutes",/,"   h -hours",/,"   d -days",/, &
   &      "  mo -months",/,"   y -years")') tunit
        endif
        stop
      endif
      
      call fiReadInt(string,kplot,ierr) 
      call fiDefaultMsg('kplot',ierr)

!     do i = 1, kplot
!       call fiReadDouble(string,tplot(i),ierr)
!       call fiDefaultMsg('tplot',ierr)
!     enddo
      
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('TIME',ierr)
      i2 = 0
      do
        i1 = i2 + 1
        i2 = i2+10
        if (i2 > kplot) i2 = kplot
        do i = i1, i2
          call fiReadDouble(string,tplot(i),ierr)
          call fiDefaultMsg('tplot',ierr)
        enddo
        if (i2 == kplot) exit
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('TIME',ierr)
      enddo

      call fiReadFlotranString(iunit1,string,ierr)
      call fiReadStringErrorMsg('TIME',ierr)
      
      call fiReadDouble(string,dt,ierr)
      call fiDefaultMsg('dt',ierr)

      call fiReadDouble(string,dtmax,ierr)
      call fiDefaultMsg('dtmax',ierr)

      if (myrank==0) then
        write(iunit2,'(/," *TIME ",a3,1x,1pe12.4,i4,1p10e12.4)') &
        tunit,tconv,kplot,(tplot(i),i=1,kplot)
        write(iunit2,'("  dt= ",1pe12.4,", dtmax= ",1pe12.4,/)') &
        dt,dtmax
      endif
      
      ! convert time units to seconds
      do i = 1, kplot
        tplot(i) = tconv * tplot(i)
      enddo
      dt = tconv * dt
      dtmax = tconv * dtmax
      
!....................
      
    case ('FELD')

      ireg = 0
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('FIELDS',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ireg = ireg + 1
      
        call fiReadInt(string,i1reg(ireg),ierr) 
        call fiDefaultMsg('i1',ierr)
        call fiReadInt(string,i2reg(ireg),ierr)
        call fiDefaultMsg('i2',ierr)
        call fiReadInt(string,j1reg(ireg),ierr)
        call fiDefaultMsg('j1',ierr)
        call fiReadInt(string,j2reg(ireg),ierr)
        call fiDefaultMsg('j2',ierr)
        call fiReadInt(string,k1reg(ireg),ierr)
        call fiDefaultMsg('k1',ierr)
        call fiReadInt(string,k2reg(ireg),ierr)
        call fiDefaultMsg('k2',ierr)

        call fiReadDouble(string,por_reg(ireg),ierr)
        call fiDefaultMsg('por',ierr)

        call fiReadDouble(string,pref_reg(ireg),ierr)
        call fiDefaultMsg('pref',ierr)

        call fiReadDouble(string,temp_reg(ireg),ierr)
        call fiDefaultMsg('temp',ierr)
      
!       if (myrank==0) write(*,*) 'ptran_read: i1,i2,j1,j2,k1,k2= ', &
!       i1reg(ireg),i2reg(ireg),j1reg(ireg),j2reg(ireg), &
!       k1reg(ireg),k2reg(ireg),por_reg(ireg),pref_reg(ireg),temp_reg(ireg)
      
      enddo
      iregfld = ireg
      
!....................
      
    case ('AQCX')
    
      ncmplx = 0
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('AQCX',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('namcx','AQCX',ierr)
      
        ncmplx = ncmplx + 1
        if (ncmplx > ncxmx) then
          if (myrank == 0) &
          write(*,*) 'too many complexes: stop: ',name,ncmplx,ncxmx
          stop
        endif
        namcx(ncmplx) = name
      enddo
      if (myrank==0) then
        write(iunit2,'(/," *AQCX: ",i3)') ncmplx
        do i = 1, ncmplx
          write(iunit2,*) i,namcx(i)
        enddo
      endif
      
!....................
      
    case ('GAS')
    
      ngas = 0
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('GAS',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('namcx','GAS',ierr)
      
        ngas = ngas + 1
        if (ngas > ngmx) then
          if (myrank == 0) &
          write(*,*) 'too many gaseous species: stop! ',name,ngas,ngmx
          stop
        endif
        namg(ngas) = name
      enddo
      if (myrank==0) then
        write(iunit2,'(/," *GAS: ",i3)') ngas
        do i = 1, ngas
          write(iunit2,*) i,namg(i)
        enddo
      endif
      
!....................
      
    case ('MNRL')
    
      mnrl = 0
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('MNRL',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('namrl','MNRL',ierr)
      
        mnrl = mnrl + 1
        if (mnrl > nmmx) then
          if (myrank == 0) &
          write(*,*) 'too many minerals: stop! ',name,mnrl,nmmx
          stop
        endif
        namrl(mnrl) = name
      enddo
      if (myrank==0) then
        write(iunit2,'(/," *MNRL",i3)') mnrl
        do i = 1, mnrl
          write(iunit2,*) i,namrl(i)
        enddo
      endif

!....................

    case ('MNIR')
    
      allocate(i1kin(nkmx,nrgmx),i2kin(nkmx,nrgmx), &
               j1kin(nkmx,nrgmx),j2kin(nkmx,nrgmx), &
               k1kin(nkmx,nrgmx),k2kin(nkmx,nrgmx))
               
      allocate(phik_reg(nkmx,nrgmx),surf_reg(nkmx,nrgmx))
    
      call fiReadStringErrorMsg('MNIR',ierr)

      nkin = 0
      npar1(1) = 1
      do ! loop over minerals
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('MNIR',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('namk','MNIR',ierr)
      
        nkin = nkin + 1
        if (nkin > nkmx) then
          if (myrank == 0) &
          write(*,*) 'too many minerals: stop: ',name,nkin,nkmx
          stop
        endif
        namk(nkin) = name
        call fiReadInt(string,npar(nkin),ierr)
        call fiDefaultMsg('npar',ierr)
        call fiReadDouble(string,fkin(nkin),ierr)
        call fiDefaultMsg('fkin',ierr)
        call fiReadDouble(string,delh(nkin),ierr)
        call fiDefaultMsg('delh',ierr)
        call fiReadDouble(string,tolpos(nkin),ierr)
        call fiDefaultMsg('tolpos',ierr)

        npar2(nkin) = npar1(nkin) + npar(nkin) - 1
!       write(*,*) 'ptranread1: ',nkin,npar1(nkin),npar2(nkin)
        
        do lp = npar1(nkin), npar2(nkin)
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('MNIR',ierr)
          call fiReadInt(string,itypkin(lp),ierr)
          call fiDefaultMsg('itypkin',ierr)
          call fiReadInt(string,nkinpri(lp),ierr)
          call fiDefaultMsg('nkinpri',ierr)
          call fiReadInt(string,nkinsec(lp),ierr)
          call fiDefaultMsg('nkinsec',ierr)
          call fiReadDouble(string,sigma(lp),ierr)
          call fiDefaultMsg('sigma',ierr)
          call fiReadDouble(string,rkf00(lp),ierr)
          call fiDefaultMsg('rk',ierr)
          do i = 1, nkinpri(lp)
            call fiReadFlotranString(iunit1,string,ierr)
            call fiReadStringErrorMsg('MNIR',ierr)
            call fiReadNChars(string,namprik(i,lp),namlen,.true.,ierr)
            call fiErrorMsg('namk','MNIR',ierr)
            call fiReadDouble(string,skinpri(i,lp),ierr)
            call fiDefaultMsg('skinpri',ierr)
          enddo
          do i = 1, nkinsec(lp)
            call fiReadFlotranString(iunit1,string,ierr)
            call fiReadStringErrorMsg('MNIR',ierr)
            call fiReadNChars(string,namseck(i,lp),namlen,.true.,ierr)
            call fiErrorMsg('namk','MNIR',ierr)
            call fiReadDouble(string,skinsec(i,lp),ierr)
            call fiDefaultMsg('skinsec',ierr)
          enddo
        enddo
        npar1(nkin+1) = npar2(nkin) + 1

        ireg = 0
        do ! loop over regions
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('MNIR',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ireg = ireg + 1
!         write(*,*) 'ptranread2: ',ireg,nkin,npar2(nkin),npar1(nkin)

          call fiReadInt(string,i1kin(nkin,ireg),ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,i2kin(nkin,ireg),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,j1kin(nkin,ireg),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,j2kin(nkin,ireg),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,k1kin(nkin,ireg),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,k2kin(nkin,ireg),ierr)
          call fiDefaultMsg('k2',ierr)

          call fiReadDouble(string,phik_reg(nkin,ireg),ierr)
          call fiDefaultMsg('phik',ierr)

          call fiReadDouble(string,surf_reg(nkin,ireg),ierr)
          call fiDefaultMsg('surf',ierr)
        enddo
        iregkin(nkin) = ireg
        if (ireg > nrgmx) then
          if (myrank == 0) &
          write(iunit2,*) 'too many mineral regions in MNIR: ',ireg
          stop
        endif
      enddo

      if (myrank==0) then
        write(iunit2,'(/," *MNIR: ",i3)') nkin
        do nr = 1, nkin
          write(iunit2,'(/," mineral               fkin        delh", &
    &     "        tolpos")')
          write(iunit2,'(1x,a20,1p3e12.4)') &
          namk(nr),fkin(nr),delh(nr),tolpos(nr)
          write(iunit2,'(" ityp  npri nsec sigma  rk")') 
          do lp = npar1(nr),npar2(nr)
            write(iunit2,'(3i3,1p2e12.4)') itypkin(lp), &
            nkinpri(lp),nkinsec(lp),sigma(lp),rkf00(lp)
            do i = 1, nkinpri(lp)
              write(iunit2,'(1x,a20,1pe12.4)') namprik(i,lp),skinpri(i,lp)
            enddo
            do i = 1, nkinsec(lp)
              write(iunit2,'(1x,a20,1pe12.4)') namseck(i,lp),skinsec(i,lp)
            enddo
          enddo
          write(iunit2,'("  i1  i2  j1  j2  k1  k2    volf        surf ")')
          do ir = 1, iregkin(nr)
            write(iunit2,'(6i4,1p2e12.4)') i1kin(nr,ir),i2kin(nr,ir), &
            j1kin(nr,ir),j2kin(nr,ir),k1kin(nr,ir),k2kin(nr,ir), &
            phik_reg(nr,ir),surf_reg(nr,ir)
          enddo
        enddo
      endif

!....................

    case ('IONX')

!-----------------------------------------------------------------------
!-----ion exchange isotherm
!     ionex        type
!       0          Gaines-Thomas
!       1          Vanselow
!       2          Gapon (not implemented)
!-----------------------------------------------------------------------

!-----ion exchange input data

      call fiReadStringErrorMsg('IONX',ierr)

      call fiReadInt(string,ionex,ierr)
      call fiDefaultMsg('ionex',ierr)
      call fiReadInt(string,ilog,ierr)
      call fiDefaultMsg('ilog',ierr)

!-----initialize counters
      nexsolid = 0 ! nr. minerals
      nexmax = 0   ! total nr. exchange cations
      nsite = 0    ! total nr. exchange sites
      n1 = 1
      ns1 = 1
      j = 0

!-----read mineral
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('IONX',ierr)

        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('nam','name',ierr)

!-------count minerals
        nexsolid = nexsolid + 1

        namex(nexsolid) = name

!-------initialize site counter per mineral
        nss = 0

        do

          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('IONX',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit

!---------count sites
          nss = nss + 1
          nsite = nsite + 1

!---------read site cec [mol/kg]
          call fiReadDouble(string,cec00,ierr)
          call fiDefaultMsg('cec',ierr)

          cec0mf(1,nsite) = cec00 ! matrix

!         if (ndloc1.eq.2 .or. idcdm.eq.1) then
!           cec0mf(1,nsite) = cec00 ! matrix
!           cec0mf(2,nsite) = cec0f ! fracture
!         else
!           cec0mf(1,nsite) = cec00 ! single continuum
!         endif

!---------read exchange ions and selectivity coefficients
          ns = 0
          do

            call fiReadFlotranString(iunit1,string,ierr)
            call fiReadStringErrorMsg('IONX',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit

            call fiReadNChars(string,name,namlen,.true.,ierr)
            call fiErrorMsg('nam','namcat',ierr)
          
            call fiReadDouble(string,eqiex0,ierr)
            call fiDefaultMsg('eqiex0',ierr)

            j = j + 1
            namcat(j) = name
          
            if (ilog == 0) then
              eqiex(j) = eqiex0
            else
              eqiex(j) = 10.d0**eqiex0
            endif
            ns = ns + 1
            nexmax = nexmax + 1
            if (nexmax.gt.nexmx) then
              if (myrank == 0) then
                write(*,*) 'too many ion exchange species-stop!',i,nexmx
                write(iunit2,*) 'too many ion exchange species!',i,nexmx
              endif
              ns = ns-1
              stop
            endif
          enddo
          n2          = n1 + ns - 1
          nex1(nsite) = n1
          nex2(nsite) = n2
          n1          = n2 + 1
        enddo
        ns2               = ns1 + nss - 1
        nsitex1(nexsolid) = ns1
        nsitex2(nexsolid) = ns2
        ns1               = ns2 + 1
      enddo

      nexsite = nsite
      nexmax = j

      if (myrank == 0) then
        write(iunit2,'(/," *IONX: ionex =",i2," ilog = ",i3," nexsolid = ", &
    &   i4," nexsite = ",i4," nexmax = ",i4)') ionex,ilog, &
        nexsolid,nexsite,nexmax
        do m = 1, nexsolid
          write(iunit2,'(i2," ",a20,2i5)') m,namex(m),nsitex1(m),nsitex2(m)
          do l = nsitex1(m), nsitex2(m)
            write(iunit2,'(10x,3i4,1pe12.4)') l,nex1(l),nex2(l),cec0mf(1,l)
            do j = nex1(l), nex2(l)
              write(iunit2,'(20x,a20,1pe12.4)') namcat(j),eqiex(j)
            enddo
          enddo
        enddo
      endif
      
      nexsolid = 0

!....................

    case ('SORP')

!-----initialize counters
      nsrfmin = 0
      nsrfmx = 0
      nsite = 0
      n1 = 1
      ns1 = 1
      do
        call fiReadFlotranString(iunit1,string,ierr)
        call fiReadStringErrorMsg('SORP',ierr)

        if (string(1:1) == '.' .or. string(1:1) == '/') exit

!-------read mineral and specific surface area [m^2/g]
!       call frfmt (ione,20,itwo,iten,izro,izro,iunit1,iunit2)
!       read (image,'(a20,10f10.0)',err=333) name,area00,cover0

!       if (name .eq. blank) goto 994

        call fiReadNChars(string,name,namlen,.true.,ierr)
        call fiErrorMsg('nam','name',ierr)
          
        call fiReadDouble(string,area00,ierr)
        call fiDefaultMsg('area00',ierr)
          
        call fiReadDouble(string,cover0,ierr)
        call fiDefaultMsg('cover0',ierr)

!-------count mineral surfaces
        nsrfmin = nsrfmin + 1
        m = nsrfmin

!-------initialize site counter
        nss = 0

        namsrf(m) = name

        if (cover0.eq.zero) then
          coverage(m) = one
        else
          coverage(m) = cover0
        endif

!-------convert area: m^2/g -> dm^2/kg
        areamass(m) = coverage(m)*area00*1.d5

        do

!---------read site species and site density
          call fiReadFlotranString(iunit1,string,ierr)
          call fiReadStringErrorMsg('SORP',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit

!         call frfmt(ione,20,itwo,20,izro,izro,iunit1,iunit2)
!         read (image,'(a20,2f20.0)',err=333) name,sited00,sited0f
          
          call fiReadNChars(string,name,namlen,.true.,ierr)
          call fiErrorMsg('nam','name',ierr)
          
          call fiReadDouble(string,sited00,ierr)
          call fiDefaultMsg('sited00',ierr)
          
!         call fiReadDouble(string,sited0f,ierr)
!         call fiDefaultMsg('sited0f',ierr)

!---------count sites
          nss = nss + 1
          nsite = nsite + 1

!---------convert site density: # sites/nm^2 --> mol_sites/dm^2: 
!            10^18 nm^2/m^2 / 6.022136736 10^23 mol^-1 * 10^-2 (m/dm)^2
          fac = one/6.022136736*1.d-7
!         siteden0(nsite) = sited00*fac

          sited0mf(1,nsite) = sited00*fac  ! single continuum
!         if (ndloc1.eq.2 .or. idcdm.eq.1) then
!           sited0mf(1,nsite) = sited00*fac  ! matrix
!           sited0mf(2,nsite) = sited0f*fac  ! fracture
!         else
!           sited0mf(1,nsite) = sited00*fac  ! single continuum
!         endif

          namsite(nsite) = name

!---------read surface complexes
          ns = 0
          do
            call fiReadFlotranString(iunit1,string,ierr)
            call fiReadStringErrorMsg('SORP',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit

            call fiReadNChars(string,name,namlen,.true.,ierr)
            call fiErrorMsg('nam','name',ierr)
          
            call fiReadDouble(string,eqsorp00,ierr)
            call fiDefaultMsg('eqsorp00',ierr)
          
!           call fiReadDouble(string,rkfsrf0,ierr)
!           call fiDefaultMsg('rkfsrf0',ierr)
          
!           call fiReadDouble(string,rkbsrf0,ierr)
!           call fiDefaultMsg('rkbsrf0',ierr)

            ns = ns + 1
            nsrfmx = nsrfmx + 1
            if (nsrfmx.gt.nscxmx) then
              write(*,*) 'too many sorption species!',i,nsrfmx,nscxmx
              write(iunit2,*) 'too many sorption species!',i,nsrfmx,nscxmx
              ns = ns-1
              stop
            endif
            namscx(nsrfmx) = name
            eqsorp0(nsrfmx) = eqsorp00
!           rkfsrf(nsrfmx) = rkfsrf0
!           rkbsrf(nsrfmx) = rkbsrf0
          enddo
          n2 = n1 + ns - 1
          nsorp1(nsite) = n1
          nsorp2(nsite) = n2
          n1 = n2 + 1
        enddo
        ns2 = ns1 + nss - 1
        nsite1(m) = ns1
        nsite2(m) = ns2
        ns1 = ns2 + 1
      enddo
      nsrfsit = nsite
      
!         if(iprint.ge.0) then
!           write(iunit2,'(30x,''site  concentration: matx frac'')')
!           write(iunit2,'(14x,a20,1p10e12.4)') name,sited00,sited0f
!    .      siteden0(nsite)
!           write(iunit2,'(27x,''species    eqsorp    kf          kb'')')
!         endif

      if (myrank == 0) then
        write(iunit2,'(/," *SORPtion")')
        do m = 1, nsrfmin
          write(iunit2,'(a20,1p2e12.4,'' [m^2/g] '')') namsrf(m), &
          areamass(m),coverage(m)
          write(iunit2,'(a20,3i5)') namsrf(m),m,nsite1(m),nsite2(m)
          do l = nsite1(m), nsite2(m)
            write(iunit2,'(10x,a20,1x,1pe12.4)') namsite(l), &
            sited0mf(1,l)
            do i = nsorp1(l), nsorp2(l)
              write(iunit2,'(12x,a20)') namscx(i)
            enddo
          enddo
        enddo
      endif
      
      nsrfmin = 0

!....................

    case default
      if (myrank == 0) &
      write(*,*) 'error reading input file: keyword not found-stop!',card
      call PetscFinalize(ierr)
      stop

    end select
  enddo
  
  close (iunit1)
  
  call trinit0
  
!*********************************************************************** 
!     read thermodynamic database
!*********************************************************************** 

!-----initialize electron stoichiometric coefficient
      do m = 1, nkin
        ze(m) = zero
      enddo
      do i = 1, ncxkin
        zeaq(i) = zero
      enddo

!-----if H2O present set nmass one less
!     ncomp: nr. of mass conservation equations consisting of all
!            primary species including H2O if present
!     nmass: primary species appearing in mass action equations,
!            exclusive of H2O if present
!     ncoll: nr. of colloids (must be last in input file)
!     ncpri: primary species including H2O used for writing out reaction
!            stoichiometry in trdatbse.f and incorporating activity of 
!            H2O for iact = 6
!-----------------------------------------------------------------------
      if (jh2o .gt. 0) then
        ncpri = ncomp
        nmass = ncomp - ncoll - 1
      else
        ncpri = ncomp + 1
!       jh2o  = ncpri ! set H2O to end of list if not explicitly present
        nmass = ncomp - ncoll
        if (icase.eq.1) then
          nam(ncpri) = 'h2o'
        else
          nam(ncpri) = 'H2O'
        endif
      endif

      if (ncpri.gt.ncmx) then
        if (myrank == 0) &
        write(*,*) 'too many primary species including h2o:---stop'
        stop
      endif

      if (myrank == 0) write(*,*) '    --> call database'
  
  call ptran_dbase (nam)

      if (myrank == 0) write(*,*) '    --> exit database'

  end subroutine ptran_read

!===================================================================

  subroutine trinit0
  
  use ptran_global_module

  implicit none
  
  integer :: i,j

!-----set species specific indices
      imaster = 1
      if (maspec == ' ') maspec = nam(imaster)
      jh2o  = 0
      jph   = 0
      joh   = 0
      jo2   = 0
      jo2aq = 0
      jco2  = 0
      jhco3 = 0
!     js2o3 = 0
      jfe3  = 0
      jfe2  = 0
      do j=1,ncomp
        if(nam(j).eq.maspec) imaster = j
        if(nam(j).eq.'H+')      jph=j
        if(nam(j).eq.'h+')      jph=j
        if(nam(j).eq.'OH-')     joh=j
        if(nam(j).eq.'oh-')     joh=j
        if(nam(j).eq.'H2O')     jh2o=j
        if(nam(j).eq.'h2o')     jh2o=j
        if(nam(j).eq.'hco3-')   jhco3=j
        if(nam(j).eq.'HCO3-')   jhco3=j
        if(nam(j).eq.'co2(aq)') jco2=j
        if(nam(j).eq.'CO2(aq)') jco2=j
        if(nam(j).eq.'CO2(AQ)') jco2=j
        if(nam(j).eq.'FE+3')    jfe3=j
        if(nam(j).eq.'fe+3')    jfe3=j
        if(nam(j).eq.'FE+2')    jfe2=j
        if(nam(j).eq.'fe+2')    jfe2=j
        if(nam(j).eq.'e-')      jpe=j
        if(nam(j).eq.'o2(aq)')  jo2aq=j
        if(nam(j).eq.'O2(aq)')  jo2aq=j
        if(nam(j).eq.'O2(AQ)')  jo2aq=j
!       if(nam(j).eq.'S2O3-2')  js2o3=j
!       if(nam(j).eq.'s2o3-2')  js2o3=j
!       if(nam(j).eq.'h2s(aq)') jh2s=j
!       if(nam(j).eq.'H2S(aq)') jh2s=j
!       if(nam(j).eq.'H2S(AQ)') jh2s=j
        if (j .le. ngas) then
          if(namg(j).eq.'o2(g)') jo2=j
          if(namg(j).eq.'O2(g)') jo2=j
          if(namg(j).eq.'O2(G)') jo2=j
        endif
      enddo
      
      jo2g  = 0
      jco2g = 0
      jh2g  = 0
      do j = 1, ngas
        if (namg(j).eq.'O2(g)')  jo2g  = j
        if (namg(j).eq.'o2(g)')  jo2g  = j
        if (namg(j).eq.'CO2(g)') jco2g = j
        if (namg(j).eq.'co2(g)') jco2g = j
        if (namg(j).eq.'H2(g)')  jh2g  = j
        if (namg(j).eq.'h2(g)')  jh2g  = j
      enddo

      jstep = imaster
      if (myrank == 0) &
      write(iunit2,'(5x,"--> master species for time-step control: ", &
      & a8)') maspec
      if (myrank == 0) &
      write(*,'(5x,"--> master species for time-step control: ",a8)') &
      maspec

!     if (molal.eq.0 .and. jh2o.eq.0) then
!       write(*,*) 'Error: species H2O missing for molarity input! stop'
!       stop
!     endif

      iph  = 0
      ife3 = 0
      ife2 = 0
      ihco3 = 0
      do i = 1, ncmplx
        if(namcx(i).eq.'H+')   iph = i
        if(namcx(i).eq.'h+')   iph = i
        if(namcx(i).eq.'FE+3') ife3=i
        if(namcx(i).eq.'fe+3') ife3=i
        if(namcx(i).eq.'FE+2') ife2=i
        if(namcx(i).eq.'fe+2') ife2=i
        if(namcx(i).eq.'hco3-') ihco3=i
        if(namcx(i).eq.'HCO3-') ihco3=i
      enddo
      
  end subroutine trinit0

  subroutine readxyz_ptran (a,n,funit)

  use fileio_module
  
  implicit none
  
  integer, intent(in) :: n, funit
  integer :: i, i1, i2, ierr, m, nvalue=10
  real*8, intent(inout) :: a(*)
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 

  save nvalue

    i2 = 0
    do
      i1 = i2+1
      i2 = i2+nvalue
      if(i2.gt.n) i2 = n
      
!     print *,'readxyz: ',i1,i2,n,FUNIT
      
      call fiReadFlotranString(FUNIT,string,ierr)
      call fiReadStringErrorMsg('DXYZ',ierr)
      do i = i1, i2
        call fiReadDouble(string, a(i), ierr)
        if (ierr .eq. 1) a(i) = 0.d0
!       print *,i,i1,i2,nvalue,a(i),n,ierr
!       call fiDefaultMsg("Error reading grid spacing", ierr)
      enddo
      do i = i1,i2
        if(a(i).eq.0.d0) then

!---------if less than nx non-zero values are read, set all the zero
!         values to the last non zero value read. Only for cartesian 
!         system

          do m = i,n
            a(m) = a(i-1)
          enddo
          return
        endif
      enddo
      if(i2.ge.n) exit
    enddo
    
  end subroutine readxyz_ptran

end module ptran_read_module
