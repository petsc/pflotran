module pflow_read_module

public
contains

  subroutine pflow_read_gridsize(inputfile, igeom, nx, ny, nz, npx, npy, npz, &
  nphase, ierr)
  
#include "include/finclude/petscsys.h"
  use petscsys
  use fileio_module
  
  implicit none
  
  character(len=*), intent(in) :: inputfile
  integer, intent(out) :: igeom, nx, ny, nz, npx, npy, npz, nphase
  integer, intent(out) :: ierr

#include "definitions.h"
  integer :: myrank, gridread_flag, commsize
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXCARDLENGTH) :: card

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  open(IUNIT1, file=inputfile, action="read", status="old") 
if (myrank==0)&
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  npx = PETSC_DECIDE; npy = PETSC_DECIDE; npz = PETSC_DECIDE
  ierr = 1
  gridread_flag = 0
  do
!   if(gridread_flag /= 0) exit
    call fiReadFlotranString(IUNIT1, string, ierr)
    if(ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    select case(card)

!....................

    case ('GRID')

      if (myrank == 0) print *, card
      call fiReadStringErrorMsg('GRID',ierr)
      
      call fiReadInt(string,igeom,ierr)
      call fiDefaultMsg('igeom',ierr)
      
      call fiReadInt(string,nx,ierr)
      call fiDefaultMsg('nx',ierr)
      
      call fiReadInt(string,ny,ierr)
      call fiDefaultMsg('ny',ierr)
      
      call fiReadInt(string,nz,ierr)
      call fiDefaultMsg('nz',ierr)
      
      call fiReadInt(string,nphase,ierr)
      call fiDefaultMsg('nphase',ierr)
      print *, "n(x,y,z)=",nx,ny,nz
    
      gridread_flag = 1
      ierr = 0

      if (myrank==0)&
       write(IUNIT2,'(/," *GRID ",/, &
     &"  igeom   = ",i4,/, &
     &"  nx      = ",i4,/, &
     &"  ny      = ",i4,/, &
     &"  nz      = ",i4,/, &
     &"  nphase  = ",i4,/, &
     &"  ndof    = ",i4,/, &
     &"  idcdm   = ",i4,/, &
     &"  itable  = ",i4    &
     &   )') igeom, nx, ny, nz, nphase

    case ('PROC')
    
      if (myrank == 0) print *, card

      call fiReadStringErrorMsg('PROC',ierr)

      call fiReadInt(string,npx,ierr)
      call fiDefaultMsg('npx',ierr)
      call fiReadInt(string,npy,ierr)
      call fiDefaultMsg('npy',ierr)
      call fiReadInt(string,npz,ierr)
      call fiDefaultMsg('npz',ierr)
 
      if (myrank == 0) &
      write(IUNIT2,'(/," *PROC",/, &
    & "  npx   = ",3x,i4,/, &
    & "  npy   = ",3x,i4,/, &
    & "  npz   = ",3x,i4)') npx,npy,npz
  
      call MPI_Comm_size(PETSC_COMM_WORLD, commsize, ierr)
    
      if (commsize .ne. npx*npy*npz) then
        if (myrank==0) write(*,*) 'Incorrect number of processors specified: ', &
        npx*npy*npz,' commsize = ',commsize
        stop
      endif

!....................

    case default
!     print *, "Entered case default"
      ! Skip anything else that isn't a GRID specification.
      exit
    endselect

  enddo
  
  close(IUNIT1)

  end subroutine pflow_read_gridsize


!#include "pflowgrid_readinput.F90"

  subroutine pflowGrid_read_input(grid, timestep,inputfile)
  
  ! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
  !           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR
    use pflow_gridtype_module
  use fileio_module
  
  implicit none

  type(pflowGrid), intent(inout) :: grid
  type(time_stepping_context), intent(inout) :: timestep
 ! type(pflow_localpatch_info) :: locpat
  
  character(len=*), intent(in) :: inputfile
  
  integer :: ierr
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXWORDLENGTH) :: word, strtim
  character(len=MAXCARDLENGTH) :: card
  integer :: i, i1, i2, idum, ireg, isrc, j
  integer :: ibc, ibrk, ir,np
  
  
 

  open(IUNIT1, file=inputfile, action="read", status="old")
  
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if(ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (grid%myrank == 0) print *, card

    select case(card)

!....................

    case ('GRID')

!....................

    case ('PROC')

!....................

    case ('COUP')

      call fiReadStringErrorMsg('COUP',ierr)

      call fiReadInt(string,grid%isync,ierr)
      call fiDefaultMsg('isync',ierr)

      if (grid%myrank == 0) &
      write(IUNIT2,'(/," *COUP",/, &
    & "  isync      = ",3x,i2 &
    & )') grid%isync

!....................

    case ('GRAV')

      call fiReadStringErrorMsg('GRAV',ierr)

      call fiReadDouble(string,grid%gravity,ierr)
      call fiDefaultMsg('gravity',ierr)

      if (grid%myrank == 0) &
      write(IUNIT2,'(/," *GRAV",/, &
    & "  gravity    = "," [m/s^2]",3x,1pe12.4 &
    & )') grid%gravity

!....................

    case ('OPTS')

      call fiReadStringErrorMsg('OPTS',ierr)

      call fiReadInt(string,grid%write_init,ierr)
      call fiDefaultMsg('write_init',ierr)

      call fiReadInt(string,grid%iprint,ierr)
      call fiDefaultMsg('iprint',ierr)

      call fiReadInt(string,grid%imod,ierr)
      call fiDefaultMsg('mod',ierr)

      call fiReadInt(string,grid%itecplot,ierr)
      call fiDefaultMsg('itecplot',ierr)

      call fiReadInt(string,grid%iblkfmt,ierr)
      call fiDefaultMsg('iblkfmt',ierr)

      call fiReadInt(string,grid%ndtcmx,ierr)
      call fiDefaultMsg('ndtcmx',ierr)

      call fiReadInt(string,grid%iran_por,ierr)
      call fiDefaultMsg('iran_por',ierr)

      call fiReadDouble(string,grid%ran_fac,ierr)
      call fiDefaultMsg('iran_por',ierr)
	  
	  call fiReadInt(string,grid%iread_perm,ierr)
      call fiDefaultMsg('iread_perm',ierr)

      if (grid%myrank == 0) &
      write(IUNIT2,'(/," *OPTS",/, &
    & "  write_init = ",3x,i2,/ &
    & "  iprint     = ",3x,i2,/, &
    & "  imod       = ",3x,i2,/, &
    & "  itecplot   = ",3x,i2,/, &
    & "  iblkfmt    = ",3x,i2,/, &
    & "  ndtcmx     = ",3x,i2,/, &
    & "  iran_por   = ",3x,i2,/, &
    & "  ran_fac    = ",3x,1pe12.4 &
    & )') grid%write_init,grid%iprint,grid%imod,grid%itecplot, &
      grid%iblkfmt,grid%ndtcmx,grid%iran_por,grid%ran_fac

!....................

    case ('TOLR')

      call fiReadStringErrorMsg('TOLR',ierr)

      call fiReadInt(string,grid%stepmax,ierr)
      call fiDefaultMsg('stepmax',ierr)

      call fiReadInt(string,timestep%iaccel,ierr)
      call fiDefaultMsg('iaccel',ierr)

      call fiReadInt(string,grid%newton_max,ierr)
      call fiDefaultMsg('newton_max',ierr)

      call fiReadInt(string,grid%icut_max,ierr)
      call fiDefaultMsg('icut_max',ierr)

      call fiReadDouble(string,timestep%dpmxe,ierr)
      call fiDefaultMsg('dpmxe',ierr)

      call fiReadDouble(string,timestep%dsmxe,ierr)
      call fiDefaultMsg('dsmxe',ierr)

      if (grid%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
     &"  flowsteps  = ",i6,/,      &
     &"  iaccel     = ",i3,/,      &
     &"  newtmx     = ",i3,/,      &
     &"  icutmx     = ",i3,/,      &
     &"  dpmxe      = ",1pe12.4,/, &
     &"  dsmxe      = ",1pe12.4)') &
! For commented-out lines to work with the Sun f95 compiler, we have to 
! terminate the string in the line above; otherwise, the compiler tries to
! include the commented-out line as part of the continued string.
      grid%stepmax,timestep%iaccel,grid%newton_max,grid%icut_max, &
      timestep%dpmxe, timestep%dsmxe

!....................

    case ('DXYZ')

	  allocate(grid%dx0(grid%nx))
      allocate(grid%dy0(grid%ny))
      allocate(grid%dz0(grid%nz))
      print *,grid%nx,grid%ny,grid%nz
      call readxyz (grid%dx0,grid%nx)
      call readxyz (grid%dy0,grid%ny)
      call readxyz (grid%dz0,grid%nz)
	  
        
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *DXYZ ")')
        write(IUNIT2,'("  dx  ",/,(1p10e12.4))') (grid%dx0(i),i=1,grid%nx)
        write(IUNIT2,'("  dy  ",/,(1p10e12.4))') (grid%dy0(i),i=1,grid%ny)
        write(IUNIT2,'("  dz  ",/,(1p10e12.4))') (grid%dz0(i),i=1,grid%nz)
      endif

!....................


   case('RAD0')
    
     call fiReadDouble(string,grid%Radius_0,ierr)
      call fiDefaultMsg('R_0',ierr)


   
!....................

    case ('RCTR')

      call fiReadStringErrorMsg('RCTR',ierr)

      call fiReadInt(string,grid%ityprxn,ierr)
      call fiDefaultMsg('ityprxn',ierr)

      call fiReadDouble(string,grid%rk,ierr)
      call fiDefaultMsg('rk',ierr)

      call fiReadDouble(string,grid%phis0,ierr)
      call fiDefaultMsg('phis0',ierr)

      call fiReadDouble(string,grid%areas0,ierr)
      call fiDefaultMsg('areas0',ierr)

      call fiReadDouble(string,grid%pwrsrf,ierr)
      call fiDefaultMsg('pwrsrf',ierr)

      call fiReadDouble(string,grid%vbars,ierr)
      call fiDefaultMsg('vbars',ierr)

      call fiReadDouble(string,grid%ceq,ierr)
      call fiDefaultMsg('ceq',ierr)

      call fiReadDouble(string,grid%delHs,ierr)
      call fiDefaultMsg('delHs',ierr)

      call fiReadDouble(string,grid%delEs,ierr)
      call fiDefaultMsg('delEs',ierr)

      call fiReadDouble(string,grid%wfmts,ierr)
      call fiDefaultMsg('wfmts',ierr)

      if (grid%myrank == 0) &
      write(IUNIT2,'(/," *RCTR",/, &
    & "  ityp   = ",3x,i3,/, &
    & "  rk     = ",3x,1pe12.4," [mol/cm^2/s]",/, &
    & "  phis0  = ",3x,1pe12.4," [-]",/, &
    & "  areas0 = ",3x,1pe12.4," [1/cm]",/, &
    & "  pwrsrf = ",3x,1pe12.4," [-]",/, &
    & "  vbars  = ",3x,1pe12.4," [cm^3/mol]",/, &
    & "  ceq    = ",3x,1pe12.4," [mol/L]",/, &
    & "  delHs  = ",3x,1pe12.4," [J/kg]",/, &
    & "  delEs  = ",3x,1pe12.4," [J/kg]",/, &
    & "  wfmts  = ",3x,1pe12.4," [g/mol]" &
    & )') grid%ityprxn,grid%rk,grid%phis0,grid%areas0,grid%pwrsrf, &
          grid%vbars,grid%ceq,grid%delHs,grid%delEs,grid%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
      grid%rk = grid%rk * grid%areas0 * 1.d3
      grid%vbars = grid%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
      grid%delHs = grid%delHs * grid%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!     grid%delHs = grid%delHs * grid%scale ! convert J/kmol -> MJ/kmol

!....................

       case ('SOLV')
    
      call fiReadStringErrorMsg('SOLV',ierr)

!     call fiReadDouble(string,eps,ierr)
!     call fiDefaultMsg('eps',ierr)

      call fiReadDouble(string,grid%atol,ierr)
      call fiDefaultMsg('atol_petsc',ierr)

      call fiReadDouble(string,grid%rtol,ierr)
      call fiDefaultMsg('rtol_petsc',ierr)

      call fiReadDouble(string,grid%stol,ierr)
      call fiDefaultMsg('stol_petsc',ierr)
      
      grid%dtol=1.D5
!     if(grid%use_ksp == 1)then
	    call fiReadDouble(string,grid%dtol,ierr)
        call fiDefaultMsg('dtol_petsc',ierr)
!     endif
   
	  call fiReadInt(string,grid%maxit,ierr)
      call fiDefaultMsg('maxit',ierr)
      
      call fiReadInt(string,grid%maxf,ierr)
      call fiDefaultMsg('maxf',ierr)

      if (grid%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
     &"  atol_petsc   = ",1pe12.4,/, &
     &"  rtol_petsc   = ",1pe12.4,/, &
     &"  stol_petsc   = ",1pe12.4,/, &
     &"  dtol_petsc   = ",1pe12.4,/, &
     &"  maxit        = ",8x,i5,/, &
     &"  maxf         = ",8x,i5 &
     &    )') &
      grid%atol,grid%rtol,grid%stol,grid%dtol,grid%maxit,grid%maxf

! The line below is a commented-out portion of the format string above.
! We have to put it here because of the stupid Sun compiler.
!    &"  eps          = ",1pe12.4,/, &
 print *,'end pckr' 
!....................

    case ('THRM')

      ireg = 0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('THRM',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ireg = ireg + 1
      
        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg('idum',ierr)

        call fiReadDouble(string,grid%rock_density(ireg),ierr)
        call fiDefaultMsg('rock_density',ierr)

        call fiReadDouble(string,grid%cpr(ireg),ierr)
        call fiDefaultMsg('cpr',ierr)
        
        call fiReadDouble(string,grid%ckdry(ireg),ierr)
        call fiDefaultMsg('ckdry',ierr)
        
        call fiReadDouble(string,grid%ckwet(ireg),ierr)
        call fiDefaultMsg('ckwet',ierr)
        
        call fiReadDouble(string,grid%tau(ireg),ierr)
        call fiDefaultMsg('tau',ierr)

        call fiReadDouble(string,grid%cdiff(ireg),ierr)
        call fiDefaultMsg('cdiff',ierr)

        call fiReadDouble(string,grid%cexp(ireg),ierr)
        call fiDefaultMsg('cexp',ierr)

        !scale thermal properties
        grid%cpr(ireg) = grid%scale * grid%cpr(ireg)
        grid%dencpr(ireg) = grid%rock_density(ireg) * grid%cpr(ireg)
        grid%ckdry(ireg) = grid%scale * grid%ckdry(ireg)
        grid%ckwet(ireg) = grid%scale * grid%ckwet(ireg)
      enddo
      
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *THRM: ",i3)') ireg
        write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
     &                 "     ckwet       tau       cdiff     cexp")')
        write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
     &                 "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
        do i = 1, ireg
          write(IUNIT2,'(i4,1p7e11.4)') i,grid%rock_density(i), &
          grid%cpr(i),grid%ckdry(i),grid%ckwet(i), &
          grid%tau(i),grid%cdiff(i),grid%cexp(i)
        enddo
      endif

!....................

    case ('PCKR')

      ireg = 0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('PCKR',ierr)

        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ireg = ireg + 1
       
        call fiReadInt(string,idum,ierr)
        call fiDefaultMsg('idum',ierr)
        grid%icaptype(ireg)=idum
       	
		   do np=1, grid%nphase
		    call fiReadDouble(string,grid%sir(np,idum),ierr)
            call fiDefaultMsg('sir',ierr)
          enddo 
		        
		call fiReadDouble(string,grid%pckrm(idum),ierr)
        call fiDefaultMsg('lambda',ierr)
        grid%lambda(ireg) = grid%pckrm(ireg)/(-grid%pckrm(ireg) +1.D0)
! Here the lambda is assigned as the same value of m

        call fiReadDouble(string,grid%alpha(idum),ierr)
        call fiDefaultMsg('alpha',ierr)

        call fiReadDouble(string,grid%pcwmax(idum),ierr)
        call fiDefaultMsg('pcwmax',ierr)
      
		call fiReadDouble(string,grid%pcbetac(idum),ierr)
        call fiDefaultMsg('pcbetac',ierr)
      
		call fiReadDouble(string,grid%pwrprm(idum),ierr)
        call fiDefaultMsg('pwrprm',ierr)

     enddo
      
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *PCKR: ",i3)') ireg
        write(IUNIT2,'("  icp swir    lambda         alpha")')
        do j = 1, ireg
		 i=grid%icaptype(j)
		     write(IUNIT2,'(i4,1p8e12.4)') i,(grid%sir(np,i),np=1,grid%nphase), &
               grid%lambda(i),grid%alpha(i),grid%pcwmax(i),grid%pcbetac(i),grid%pwrprm(i)
		enddo
     end if

!....................
      
    case ('PHIK')

      ireg = 0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('PHIK',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ireg = ireg + 1
        
        if (ireg > MAXPERMREGIONS) then
          print *,'Error reading PHIK keyword: too many regions-stop',ireg
          stop
        endif
      
        call fiReadInt(string,grid%i1reg(ireg),ierr) 
        call fiDefaultMsg('i1',ierr)
        call fiReadInt(string,grid%i2reg(ireg),ierr)
        call fiDefaultMsg('i2',ierr)
        call fiReadInt(string,grid%j1reg(ireg),ierr)
        call fiDefaultMsg('j1',ierr)
        call fiReadInt(string,grid%j2reg(ireg),ierr)
        call fiDefaultMsg('j2',ierr)
        call fiReadInt(string,grid%k1reg(ireg),ierr)
        call fiDefaultMsg('k1',ierr)
        call fiReadInt(string,grid%k2reg(ireg),ierr)
        call fiDefaultMsg('k2',ierr)

        call fiReadInt(string,grid%icap_reg(ireg),ierr)
        call fiDefaultMsg('icap',ierr)

        call fiReadInt(string,grid%ithrm_reg(ireg),ierr)
        call fiDefaultMsg('ithrm',ierr)

        call fiReadDouble(string,grid%por_reg(ireg),ierr)
        call fiDefaultMsg('por',ierr)

        call fiReadDouble(string,grid%tor_reg(ireg),ierr)
        call fiDefaultMsg('tor',ierr)

        call fiReadDouble(string,grid%perm_reg(ireg,1),ierr)
        call fiDefaultMsg('permx',ierr)

        call fiReadDouble(string,grid%perm_reg(ireg,2),ierr)
        call fiDefaultMsg('permy',ierr)

        call fiReadDouble(string,grid%perm_reg(ireg,3),ierr)
        call fiDefaultMsg('permz',ierr)

        call fiReadDouble(string,grid%perm_reg(ireg,4),ierr)
        call fiDefaultMsg('permpow',ierr)

!		call fiReadDouble(string,grid%Perm_reg(ireg,5),ierr)
!       call fiDefaultMsg('porokin',ierr)

		
      enddo
      grid%iregperm = ireg
            
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *PHIK: ireg = ",i4)') grid%iregperm
        write(IUNIT2,'("  i1  i2  j1  j2  k1  k2 icap ithrm  por      tor  ", &
      & "     permx      permy      permz [m^2]")')
        do ireg = 1, grid%iregperm
          write(IUNIT2,'(6i4,2i4,1p6e11.4)') &
          grid%i1reg(ireg),grid%i2reg(ireg), &
          grid%j1reg(ireg),grid%j2reg(ireg), &
          grid%k1reg(ireg),grid%k2reg(ireg), &
          grid%icap_reg(ireg),grid%ithrm_reg(ireg), &
          grid%por_reg(ireg),grid%tor_reg(ireg),(grid%perm_reg(ireg,i),i=1,4)
        enddo
      endif

!....................
      
    case ('INIT')
    
      call fiReadInt(string,grid%iread_init,ierr) 
      call fiDefaultMsg('iread_init',ierr)
      
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *INIT: iread = ",i2)') grid%iread_init
      endif
      
      if (grid%iread_init == 0 .or. grid%iread_init == 2) then
      
		  ireg = 0
		  do
			call fiReadFlotranString(IUNIT1,string,ierr)
			call fiReadStringErrorMsg('INIT',ierr)

			if (string(1:1) == '.' .or. string(1:1) == '/') exit
			ireg = ireg + 1

			call fiReadInt(string,grid%i1ini(ireg),ierr) 
			call fiDefaultMsg('i1',ierr)
			call fiReadInt(string,grid%i2ini(ireg),ierr)
			call fiDefaultMsg('i2',ierr)
			call fiReadInt(string,grid%j1ini(ireg),ierr)
			call fiDefaultMsg('j1',ierr)
			call fiReadInt(string,grid%j2ini(ireg),ierr)
			call fiDefaultMsg('j2',ierr)
			call fiReadInt(string,grid%k1ini(ireg),ierr)
			call fiDefaultMsg('k1',ierr)
			call fiReadInt(string,grid%k2ini(ireg),ierr)
			call fiDefaultMsg('k2',ierr)
		   
				 
			 do j=1,grid%ndof
				call fiReadDouble(string,grid%xx_ini(j,ireg),ierr)
				call fiDefaultMsg('xxini',ierr)
			 enddo
		  enddo
		  
		  grid%iregini = ireg
		  
		  if (grid%myrank==0) then
			write(IUNIT2,'("  ireg = ",i4)') grid%iregini
			write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]      ", &
		&   "sl [-]      c [mol/L]")')
			do ireg = 1, grid%iregini
			  write(IUNIT2,'(6i6,1p10e12.4)') &
			  grid%i1ini(ireg),grid%i2ini(ireg), &
			  grid%j1ini(ireg),grid%j2ini(ireg), &
			  grid%k1ini(ireg),grid%k2ini(ireg), &
			  (grid%xx_ini(np,ireg),np =1,grid%ndof)
			 enddo
		  endif

    else if (grid%iread_init == 1) then
    
!     read in initial conditions from file: pflow_init.dat
      if (grid%myrank == 0) then
        write(*,*) '--> read in initial conditions from file: pflow_init.dat'

        open(IUNIT3, file='pflow_init.dat', action="read", status="old")

        ireg = 0
        do
          call fiReadFlotranString(IUNIT3,string,ierr)
!         call fiReadStringErrorMsg('INIT',ierr)

          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ireg = ireg + 1

          call fiReadInt(string,grid%i1ini(ireg),ierr) 
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,grid%i2ini(ireg),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,grid%j1ini(ireg),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,grid%j2ini(ireg),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,grid%k1ini(ireg),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,grid%k2ini(ireg),ierr)
          call fiDefaultMsg('k2',ierr)

                
		 do j=1,grid%ndof
		    call fiReadDouble(string,grid%xx_ini(j,ireg),ierr)
            call fiDefaultMsg('xx_ini',ierr)
  	     enddo
       enddo
        grid%iregini = ireg
        close(IUNIT3)
      endif
    endif

!....................

    case ('TIME')

      call fiReadStringErrorMsg('TIME',ierr)
      
      call fiReadWord(string,strtim,.false.,ierr)
      
      timestep%tunit = strtim

      if (timestep%tunit == 's') then
        grid%tconv = 1.d0
      else if (timestep%tunit == 'm') then
        grid%tconv = 60.d0
      else if (timestep%tunit == 'h') then
        grid%tconv = 60.d0 * 60.d0
      else if (timestep%tunit == 'd') then
        grid%tconv = 60.d0 * 60.d0 * 24.d0
      else if (timestep%tunit == 'mo') then
        grid%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
      else if (timestep%tunit == 'y') then
        grid%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
      else
        if (grid%myrank == 0) then
          write(*,'(" Time unit: ",a3,/, &
   &      " Error: time units must be one of ",/, &
   &      "   s -seconds",/,"   m -minutes",/,"   h -hours",/,"   d -days",/, &
   &      "  mo -months",/,"   y -years")') timestep%tunit
        endif
        stop
      endif

      call fiReadInt(string,timestep%kplot,ierr) 
      call fiDefaultMsg('kplot',ierr)
      
      allocate(timestep%tplot(timestep%kplot))
      
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('TIME',ierr)
      i2 = 0
      do
        i1 = i2 + 1
        i2 = i2+10
        if (i2 > timestep%kplot) i2 = timestep%kplot
        do i = i1, i2
          call fiReadDouble(string,timestep%tplot(i),ierr)
          call fiDefaultMsg('tplot',ierr)
        enddo
        if (i2 == timestep%kplot) exit
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('TIME',ierr)
      enddo

!     call fiReadFlotranString(IUNIT1,string,ierr)
!     call fiReadStringErrorMsg('TIME',ierr)
      
!     call fiReadDouble(string,grid%dt,ierr)
!     call fiDefaultMsg('dt',ierr)

!     call fiReadDouble(string,grid%dt_max,ierr)
!     call fiDefaultMsg('dt_max',ierr)

      if (grid%myrank==0) then
        write(IUNIT2,'(/," *TIME ",a3,1x,i4,/,(1p10e12.4))') timestep%tunit, &
        timestep%kplot,(timestep%tplot(i),i=1,timestep%kplot)
!       write(IUNIT2,'("  dt= ",1pe12.4,", dtmax= ",1pe12.4,/)') &
!       grid%dt,grid%dt_max
      endif
      
      ! convert time units to seconds
      do i = 1, timestep%kplot
        timestep%tplot(i) = grid%tconv * timestep%tplot(i)
      enddo
!     grid%dt = grid%tconv * grid%dt
!     grid%dt_max = grid%tconv * grid%dt_max

!....................

    case ('DTST')

      call fiReadStringErrorMsg('DTST',ierr)

      call fiReadInt(string,timestep%nstpmax,ierr)
      call fiDefaultMsg('nstpmax',ierr)

      allocate(timestep%tstep(timestep%nstpmax))
      allocate(timestep%dtstep(timestep%nstpmax))


      do i = 1, timestep%nstpmax
        call fiReadDouble(string,timestep%tstep(i),ierr)
        call fiDefaultMsg('tstep',ierr)
      enddo

      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('DTST',ierr)
      call fiReadDouble(string,timestep%dt_min,ierr)
      call fiDefaultMsg('dt_min',ierr)
      do i = 1, timestep%nstpmax
        call fiReadDouble(string,timestep%dtstep(i),ierr)
        call fiDefaultMsg('dtstep',ierr)
      enddo
      
      timestep%dt_max = timestep%dtstep(1)
      
      grid%dt = timestep%dt_min
      
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *DTST ",i4,/," tstep= ",(1p10e12.4))') timestep%nstpmax, &
        (timestep%tstep(i),i=1,timestep%nstpmax)
        write(IUNIT2,'(" dtstep= ",1p10e12.4,/)') &
        timestep%dt_min,(timestep%dtstep(i),i=1,timestep%nstpmax)
      endif
      
      ! convert time units to seconds
      do i = 1, timestep%nstpmax
        timestep%tstep(i) = grid%tconv * timestep%tstep(i)
        timestep%dtstep(i) = grid%tconv * timestep%dtstep(i)
      enddo
      grid%dt = grid%tconv * grid%dt
      timestep%dt_min = grid%tconv * timestep%dt_min
      timestep%dt_max = grid%tconv * timestep%dt_max

!....................

    case ('BCON')

!-----------------------------------------------------------------------
!-----boundary conditions:	ibnd:	
!                   1-left,		2-right
!					3-top,		4-bottom
!					5-front,	6-back
!
!	ibndtyp:	1-Dirichlet         (p, T, C)
!	ibndtyp:	2-Neumann/Dirichlet (q, grad T=0, grad C=0)
!	ibndtyp:	3-Dirichlet/Neumann (p, grad T=0, grad C=0)
!-----------------------------------------------------------------------
      ibc = 0
      ir = 0
      grid%iregbc1(1) = 1
      do ! loop over blocks
      
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('BCON',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        ibc = ibc + 1  ! RTM: Number of boundary conditions
        
        call fiReadInt(string,grid%ibndtyp(ibc),ierr)
        call fiDefaultMsg('ibndtyp',ierr)

        call fiReadInt(string,grid%iface(ibc),ierr)
        call fiDefaultMsg('iface',ierr)

        do ! loop over regions
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ir = ir + 1

          call fiReadInt(string,grid%i1bc(ir),ierr)
          call fiDefaultMsg('i1',ierr)
          call fiReadInt(string,grid%i2bc(ir),ierr)
          call fiDefaultMsg('i2',ierr)
          call fiReadInt(string,grid%j1bc(ir),ierr)
          call fiDefaultMsg('j1',ierr)
          call fiReadInt(string,grid%j2bc(ir),ierr)
          call fiDefaultMsg('j2',ierr)
          call fiReadInt(string,grid%k1bc(ir),ierr)
          call fiDefaultMsg('k1',ierr)
          call fiReadInt(string,grid%k2bc(ir),ierr)
          call fiDefaultMsg('k2',ierr)    

          ! Now read the velocities or pressures, depending on the BC type
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BCON',ierr)
   !       do j=1, grid%nphase
   
			 
			 if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then 
				do j=1,grid%ndof
				  call fiReadDouble(string,grid%xxbc0(j,ibc),ierr)
				  call fiDefaultMsg('xxbc',ierr)
				  enddo
			 
			  elseif(grid%ibndtyp(ibc) == 2) then
				   do j=1, grid%nphase       
					 call fiReadDouble(string, grid%velocitybc0(j,ibc), ierr)
					 call fiDefaultMsg("Error reading velocity BCs:", ierr)			  		  	  
				   enddo
				   do j=2,grid%ndof
					 call fiReadDouble(string,grid%xxbc0(j,ibc),ierr)
					 call fiDefaultMsg('xxbc',ierr)
				  enddo

			   endif
					       
        enddo ! End loop over regions.
        
        grid%iregbc2(ibc) = ir
        if (ibc+1 > MAXBCBLOCKS) then
          write(*,*) 'Too many boundary condition blocks specified--stop: ', &
          ibc+1, MAXBCBLOCKS
          stop
        else
          grid%iregbc1(ibc+1) = grid%iregbc2(ibc)+1
        endif
      enddo ! End loop over blocks.
      
      grid%nblkbc = ibc
      
      if (grid%myrank == 0) then
        write(IUNIT2,'(/," *BCON: nblkbc = ",i4)') grid%nblkbc
        do ibc = 1, grid%nblkbc
          write(IUNIT2,'("  ibndtyp = ",i3," iface = ",i2)') grid%ibndtyp(ibc), &
          grid%iface(ibc)
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]    c", &
     &    " [mol/L]")')
          do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
             if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
                write(IUNIT2,'(6i6,1p10e12.4)') &
                grid%i1bc(ireg),grid%i2bc(ireg), &
                grid%j1bc(ireg),grid%j2bc(ireg), &
                grid%k1bc(ireg),grid%k2bc(ireg), &
                (grid%xxbc0(j,ireg),j=1,grid%ndof)
              else if (grid%ibndtyp(ibc) == 2) then
                write(IUNIT2,'(6i4,1p10e12.4)') &
                grid%i1bc(ireg),grid%i2bc(ireg), &
                grid%j1bc(ireg),grid%j2bc(ireg), &
                grid%k1bc(ireg),grid%k2bc(ireg), &
                (grid%velocitybc0(j,ireg),j=1,grid%nphase),&
                (grid%xxbc0(j,ireg),j=2,grid%ndof)
              endif
         enddo
        enddo
      endif

!....................

    case ('SOUR')

      isrc = 0
      ir = 0
      
      do ! loop over sources
      
			call fiReadFlotranString(IUNIT1,string,ierr)
			call fiReadStringErrorMsg('SOUR',ierr)
		  
			if (string(1:1) == '.' .or. string(1:1) == '/') exit

			isrc = isrc + 1  ! Number of sources

			ir = ir + 1

			call fiReadInt(string,grid%i1src(ir),ierr)
			call fiDefaultMsg('i1',ierr)
			call fiReadInt(string,grid%i2src(ir),ierr)
			call fiDefaultMsg('i2',ierr)
			call fiReadInt(string,grid%j1src(ir),ierr)
			call fiDefaultMsg('j1',ierr)
			call fiReadInt(string,grid%j2src(ir),ierr)
			call fiDefaultMsg('j2',ierr)
			call fiReadInt(string,grid%k1src(ir),ierr)
			call fiDefaultMsg('k1',ierr)
			call fiReadInt(string,grid%k2src(ir),ierr)
			call fiDefaultMsg('k2',ierr)    
			print *,'Source', isrc, ir   
			! Read time, temperature, q-source
			i = 0
			do ! loop over time intervals
			
				  call fiReadFlotranString(IUNIT1,string,ierr)
				  call fiReadStringErrorMsg('SOUR',ierr)
				  if (string(1:1) == '.' .or. string(1:1) == '/') exit
				
				  i = i + 1
				
				  if (i+1 > 10) then
					write(*,*) 'Too many times specified in SOURce--stop: ', i+1, 10
					stop
				  endif
				
				  call fiReadDouble(string,grid%timesrc(i,isrc),ierr)
				  call fiDefaultMsg('timesrc',ierr)

				  call fiReadDouble(string,grid%tempsrc(i,isrc),ierr)
				  call fiDefaultMsg('tempsrc',ierr)
			  
				  do j=1, grid%nphase
					call fiReadDouble(string,grid%qsrc(i,isrc,j),ierr)
					call fiDefaultMsg('qsrc',ierr)
				  enddo 
				 
			enddo ! End loop over time.

			grid%ntimsrc = i

			if (grid%ntimsrc > MAXSRCTIMES) then
			  write(*,*) 'Too many source times specified--stop: ', &
			  grid%ntimsrc, MAXSRCTIMES
			  stop
			endif

			if (isrc+1 > MAXSRC) then
			  write(*,*) 'Too many source blocks specified--stop: ', &
			  isrc+1, MAXSRC
			  stop
			endif
      enddo ! End loop over sources.

      grid%nblksrc = isrc

      if (grid%myrank == 0) then
        write(IUNIT2,'(/," *SOURce: nblksrc = ",i4)') grid%nblksrc
        do isrc = 1, grid%nblksrc
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2")')
          write(IUNIT2,'(6i4)') &
          grid%i1src(isrc),grid%i2src(isrc), &
          grid%j1src(isrc),grid%j2src(isrc), &
          grid%k1src(isrc),grid%k2src(isrc)
          write(IUNIT2,'("    t [s]        T [C]    QH2O [kg/s]    QCO2 [kg/s]")')
          do ir = 1, grid%ntimsrc
            write(IUNIT2,'(1p10e12.4)') &
            grid%timesrc(ir,isrc),grid%tempsrc(ir,isrc),grid%qsrc(ir,isrc,:)
            
          enddo
        enddo
      endif

!....................
      
    case ('BRK')

      ibrk = 0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('BRK',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit
        ibrk = ibrk + 1
      
        call fiReadInt(string,grid%i1brk(ibrk),ierr) 
        call fiDefaultMsg('i1',ierr)
        call fiReadInt(string,grid%i2brk(ibrk),ierr)
        call fiDefaultMsg('i2',ierr)
        call fiReadInt(string,grid%j1brk(ibrk),ierr)
        call fiDefaultMsg('j1',ierr)
        call fiReadInt(string,grid%j2brk(ibrk),ierr)
        call fiDefaultMsg('j2',ierr)
        call fiReadInt(string,grid%k1brk(ibrk),ierr)
        call fiDefaultMsg('k1',ierr)
        call fiReadInt(string,grid%k2brk(ibrk),ierr)
        call fiDefaultMsg('k2',ierr)

        call fiReadInt(string,grid%ibrktyp(ibrk),ierr)
        call fiDefaultMsg('ibrktyp',ierr)

        call fiReadInt(string,grid%ibrkface(ibrk),ierr)
        call fiDefaultMsg('ibrkface',ierr)
      enddo
      grid%ibrkcrv = ibrk
            
      if (grid%myrank==0) then
        write(IUNIT2,'(/," *BRK: ibrk = ",i4)') grid%ibrkcrv
        write(IUNIT2,'("  i1  i2  j1  j2  k1  k2  ibrktyp  ibrkface  ")')
        do ibrk = 1, grid%ibrkcrv
          write(IUNIT2,'(6i4,4x,i2,7x,i2)') &
          grid%i1brk(ibrk),grid%i2brk(ibrk), &
          grid%j1brk(ibrk),grid%j2brk(ibrk), &
          grid%k1brk(ibrk),grid%k2brk(ibrk),grid%ibrktyp(ibrk), &
          grid%ibrkface(ibrk)
        enddo
      endif

      if (grid%ndof == 1) grid%ibrkcrv = 0

!....................

    case default
    
      if (grid%myrank == 0) then
        print *, "Error reading input file: keyword not found. Terminating."
      endif
      call PetscFinalize(ierr)
      stop

    end select

  enddo

  close(IUNIT1)
  
  end subroutine pflowgrid_read_input

!======================================================================

  subroutine readxyz (a,n)

  use fileio_module
  
  implicit none
  
  integer*4, intent(in) :: n
  integer*4 :: i, i1, i2, m
  integer ::  ierr, nvalue=10
  real*8, intent(inout) :: a(*)
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 

  save nvalue

!   call fiReadStringErrorMsg('DXYZ',ierr)

  !  call fiReadDouble(string,grid%radius_0,ierr)
   ! call fiDefaultMsg('radius_0',ierr)

    i2 = 0
    do
      i1 = i2+1
      i2 = i2+nvalue
      if(i2.gt.n) i2 = n
      call fiReadFlotranString(IUNIT1,string,ierr)
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
    
  end subroutine readxyz

!======================================================================






end module pflow_read_module
