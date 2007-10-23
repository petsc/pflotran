module Init_module

  use pflow_gridtype_module
  use Grid_module

  implicit none

  private

#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
 
contains

! ************************************************************************** !
!
! initPFLOW: Initializes a pflow grid object
! author: Glenn Hammond
! date: 10/23/07
!
! **************************************************************************
subroutine initPFLOW(solution,grid,filename)

  use Solution_module
  use Grid_module
  use Solver_module
  use Material_module
  use fileio_module
  
  implicit none

  type(solution_type) :: solution
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: filename

  integer :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXNAMELENGTH) :: name
  character(len=MAXCARDLENGTH) :: card
  
  type(solver_type), pointer :: solver
  type(material_type), pointer :: material

  integer :: mcomp, mphas, nx, ny, nz, igeom
  real*8, parameter:: fmwnacl = 58.44277D0, fmwh2o  = 18.01534d0
  integer :: i, i1, i2, idum, ireg, isrc, j
  integer :: ibc, ibrk, ir,np  
  
  call MPI_Comm_rank(PETSC_COMM_WORLD,solution%myrank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,solution%commsize,ierr)

  open(IUNIT1, file=filename, action="read", status="old") 
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  mphas=0
  mcomp = 0

! Read in select required cards
!.........................................................................

  ! GRID information
  string = "GRID"
  call fiFindStringInFile(IUNIT1,string,ierr)
  call fiFindStringErrorMsg(string,ierr)

  ! strip card from front of string
  call fiReadWord(string,word,.false.,ierr)
 
  ! key off igeom for structured vs unstructured 
  call fiReadInt(string,igeom,ierr)
  call fiDefaultMsg('igeom',ierr)
  
  grid => createGrid(igeom) 

  if (grid%is_structured) then ! structured
    call fiReadInt(string,grid%structured_grid%nx,ierr)
    call fiDefaultMsg('nx',ierr)
    
    call fiReadInt(string,grid%structured_grid%ny,ierr)
    call fiDefaultMsg('ny',ierr)
    
    call fiReadInt(string,grid%structured_grid%nz,ierr)
    call fiDefaultMsg('nz',ierr)
    
    grid%structured_grid%nxy = grid%structured_grid%nx*grid%structured_grid%ny
    grid%structured_grid%nmax = grid%structured_grid%nxy * &
                                grid%structured_grid%nz
  else ! unstructured
  endif
      
  call fiReadInt(string,solution%nphase,ierr)
  call fiDefaultMsg('nphase',ierr)

  call fiReadInt(string,solution%nspec,ierr)
  call fiDefaultMsg('nspec',ierr)

  call fiReadInt(string,solution%npricomp,ierr)
  call fiDefaultMsg('npricomp',ierr)

  call fiReadInt(string,solution%ndof,ierr)
  call fiDefaultMsg('ndof',ierr)
      
  call fiReadInt(string,solution%idcdm,ierr)
  call fiDefaultMsg('idcdm',ierr)

  call fiReadInt(string,solution%itable,ierr)
  call fiDefaultMsg('itable',ierr)

  if (solution%myrank==0) then
    if (grid%is_structured) then
      write(IUNIT2,'(/," *GRID ",/, &
        &"  igeom   = ",i4,/, &
        &"  nx      = ",i4,/, &
        &"  ny      = ",i4,/, &
        &"  nz      = ",i4,/, &
        &"  nphase  = ",i4,/, &
        &"  ndof    = ",i4,/, &
        &"  idcdm   = ",i4,/, &
        &"  itable  = ",i4    &
        &   )') solution%igeom, grid%structured_grid%nx, &
                grid%structured_grid%ny, grid%structured_grid%nz, &
                solution%nphase, solution%ndof, solution%idcdm, solution%itable
    else
    endif
  endif

!.........................................................................

  if (grid%is_structured) then  ! look for processor decomposition
    
    ! PROC information
    string = "PROC"
    call fiFindStringInFile(IUNIT1,string,ierr)

    if (ierr /= 0) then

      ! strip card from front of string
      call fiReadWord(string,word,.false.,ierr)
      call fiReadInt(string,grid%structured_grid%nx,ierr)
      call fiDefaultMsg('npx',ierr)
      call fiReadInt(string,grid%structured_grid%ny,ierr)
      call fiDefaultMsg('npy',ierr)
      call fiReadInt(string,grid%structured_grid%nz,ierr)
      call fiDefaultMsg('npz',ierr)
 
      if (solution%myrank == 0) &
        write(IUNIT2,'(/," *PROC",/, &
          & "  npx   = ",3x,i4,/, &
          & "  npy   = ",3x,i4,/, &
          & "  npz   = ",3x,i4)') grid%structured_grid%nx, &
            grid%structured_grid%ny, grid%structured_grid%nz
  
        call MPI_Comm_size(PETSC_COMM_WORLD,solution%commsize,ierr)
        if (solution%commsize /= grid%structured_grid%nx * &
                                 grid%structured_grid%ny * &
                                 grid%structured_grid%nz) then
          if (solution%myrank==0) &
            write(*,*) 'Incorrect number of processors specified: ', &
                       grid%structured_grid%nx*grid%structured_grid%ny* &
                       grid%structured_grid%nz,' commsize = ',solution%commsize
        stop
      endif
    endif
  endif
  
!.........................................................................

  ! COMP information
  string = "COMP"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr /= 0) then

    mcomp = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('COMP',ierr)
      
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','GAS',ierr)
        
      call fiWordToUpper(name) 
      select case(name(1:len_trim(name)))
        case('H2O')
            mcomp = mcomp +1
        case('TRACER')     
            mcomp = mcomp +4
        case('CO2')
            mcomp = mcomp +2   
        case('C10H22')
            mcomp = mcomp +8
        case('AIR')              
            mcomp = mcomp +16
        case('ENG')
            mcomp = mcomp +32
        case default               
           print *,' Wrong comp name::  Stop:' 
           stop
      end select 
    enddo
  endif          

!....................................................................

  ! COMP information
  string = "PHAS"
  call fiFindStringInFile(IUNIT1,string,ierr)

  if (ierr /= 0) then

    mphas = 0
    do
      call fiReadFlotranString(IUNIT1,string,ierr)
      call fiReadStringErrorMsg('phase',ierr)
    
      if (string(1:1) == '.' .or. string(1:1) == '/') exit

      call fiReadWord(string,name,.true.,ierr)
      call fiErrorMsg('namcx','phase',ierr)
        
      call fiWordToUpper(name) 
         
      select case(name(1:len_trim(name)))
        case('ROCK')
          mphas = mphas +1
        case('H2O')     
          mphas = mphas +2
        case('CO2')
          mphas = mphas +4   
        case('GAS')
          mphas = mphas +8
        case('NAPL1')              
          mphas = mphas +16
        case('NAPL2')
          mphas = mphas +32
        case default               
          print *,' Wrong phase name::  Stop:' 
          stop
      end select 
    enddo
  endif
    
!....................................................................


! keywords: GRID, PROC, COUP, GRAV, OPTS, TOLR, DXYZ, DIFF, RADN, HYDR,  
!           SOLV, THRM, PCKR, PHIK, INIT, TIME, DTST, BCON, SOUR, BRK, RCTR
  
#if 0
  type(pflowGrid), intent(inout) :: grid
  character(len=*), intent(in) :: inputfile
  integer :: ierr
#include "definitions.h"
  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXWORDLENGTH) :: word, strtim
  character(len=MAXCARDLENGTH) :: card

  
  open(IUNIT1, file=inputfile, action="read", status="old")
#endif
  
  rewind(IUNIT1)
  do
    call fiReadFlotranString(IUNIT1, string, ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    call fiCharsToUpper(word,len_trim(word))
    call fiReadCard(word,card,ierr)

    if (solution%myrank == 0) print *, 'pflow_read:: ',card

    select case(card)

!....................

      case ('GRID')

!....................

      case ('PROC')
!.....................
      case ('COMP') 
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('COMP',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo
!.....................
      case ('PHAS')
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('PHASE',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
        enddo 
      
!....................

      case ('COUP')

        call fiReadStringErrorMsg('COUP',ierr)

        call fiReadInt(string,solution%isync,ierr)
        call fiDefaultMsg('isync',ierr)

        if (solution%myrank == 0) &
          write(IUNIT2,'(/," *COUP",/, &
            & "  isync      = ",3x,i2 &
            & )') solution%isync

!....................

      case ('GRAV')

        call fiReadStringErrorMsg('GRAV',ierr)

        call fiReadDouble(string,solution%gravity,ierr)
        call fiDefaultMsg('gravity',ierr)

        if (solution%myrank == 0) &
          write(IUNIT2,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1pe12.4 &
            & )') solution%gravity

!....................

      case ('HDF5')
        solution%print_hdf5 = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              solution%print_hdf5_velocities = .true.
            case('FLUX')
              solution%print_hdf5_flux_velocities = .true.
            case default
          end select
            
        enddo

        if (solution%myrank == 0) &
          write(IUNIT2,'(/," *HDF5",10x,i1,/)') solution%print_hdf5

!....................

      case ('TECP')
        solution%print_tecplot = .true.
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          call fiCharsToUpper(word,len_trim(word))
          call fiReadCard(word,card,ierr)

          select case(card)
            case('VELO')
              solution%print_tecplot_velocities = .true.
            case('FLUX')
              solution%print_tecplot_flux_velocities = .true.
            case default
          end select
          
        enddo

        if (solution%myrank == 0) &
          write(IUNIT2,'(/," *TECP",10x,i1,/)') solution%print_tecplot

!....................


      case ('OPTS')

        call fiReadStringErrorMsg('OPTS',ierr)

        call fiReadInt(string,solution%write_init,ierr)
        call fiDefaultMsg('write_init',ierr)

        call fiReadInt(string,solution%iprint,ierr)
        call fiDefaultMsg('iprint',ierr)

        call fiReadInt(string,solution%imod,ierr)
        call fiDefaultMsg('mod',ierr)

        call fiReadInt(string,solution%itecplot,ierr)
        call fiDefaultMsg('itecplot',ierr)

        call fiReadInt(string,solution%iblkfmt,ierr)
        call fiDefaultMsg('iblkfmt',ierr)

        call fiReadInt(string,solution%ndtcmx,ierr)
        call fiDefaultMsg('ndtcmx',ierr)

        call fiReadInt(string,solution%iran_por,ierr)
        call fiDefaultMsg('iran_por',ierr)
  
        call fiReadDouble(string,solution%ran_fac,ierr)
        call fiDefaultMsg('ran_fac',ierr)
    
        call fiReadInt(string,solution%iread_perm,ierr)
        call fiDefaultMsg('iread_perm',ierr)
    
        call fiReadInt(string,solution%iread_geom,ierr)
        call fiDefaultMsg('iread_geom',ierr)


        if (solution%myrank == 0) &
          write(IUNIT2,'(/," *OPTS",/, &
            & "  write_init = ",3x,i2,/ &
            & "  iprint     = ",3x,i2,/, &
            & "  imod       = ",3x,i2,/, &
            & "  itecplot   = ",3x,i2,/, &
            & "  iblkfmt    = ",3x,i2,/, &
            & "  ndtcmx     = ",3x,i2,/, &
            & "  iran_por   = ",3x,i2,/, &
            & "  ran_fac    = ",3x,1pe12.4,/, &
            & "  iread_perm = ",3x,i2,/, &
            & "  iread_geom = ",3x,i2 &
            & )') solution%write_init,solution%iprint,solution%imod,solution%itecplot, &
            solution%iblkfmt,solution%ndtcmx,solution%iran_por,solution%ran_fac, &
            solution%iread_perm,solution%iread_geom

!....................

      case ('TOLR')

        call fiReadStringErrorMsg('TOLR',ierr)

        call fiReadInt(string,solution%stepmax,ierr)
        call fiDefaultMsg('stepmax',ierr)
  
        call fiReadInt(string,solution%iaccel,ierr)
        call fiDefaultMsg('iaccel',ierr)

        call fiReadInt(string,solution%newton_max,ierr)
        call fiDefaultMsg('newton_max',ierr)

        call fiReadInt(string,solution%icut_max,ierr)
        call fiDefaultMsg('icut_max',ierr)

        call fiReadDouble(string,solution%dpmxe,ierr)
        call fiDefaultMsg('dpmxe',ierr)

        call fiReadDouble(string,solution%dtmpmxe,ierr)
        call fiDefaultMsg('dtmpmxe',ierr)
  
        call fiReadDouble(string,solution%dcmxe,ierr)
        call fiDefaultMsg('dcmxe',ierr)

        call fiReadDouble(string,solution%dsmxe,ierr)
        call fiDefaultMsg('dsmxe',ierr)

        if (solution%myrank==0) write(IUNIT2,'(/," *TOLR ",/, &
          &"  flowsteps  = ",i6,/,      &
          &"  iaccel     = ",i3,/,      &
          &"  newtmx     = ",i3,/,      &
          &"  icutmx     = ",i3,/,      &
          &"  dpmxe      = ",1pe12.4,/, &
          &"  dtmpmxe    = ",1pe12.4,/, &
          &     "  dcmxe      = ",1pe12.4,/, &
          &"  dsmxe      = ",1pe12.4)') &
! For commented-out lines to work with the Sun f95 compiler, we have to 
! terminate the string in the line above; otherwise, the compiler tries to
! include the commented-out line as part of the continued string.
          solution%stepmax,solution%iaccel,solution%newton_max,solution%icut_max, &
          solution%dpmxe,solution%dtmpmxe,solution%dcmxe, solution%dsmxe

!....................

      case ('DXYZ')
      
        if (grid%is_structured) then  ! look for processor decomposition
          call readStructuredDXYZ(solution,grid%structured_grid)
        else
          if (solution%myrank == 0) &
            print *, 'ERROR: Keyword "DXYZ" not supported for unstructured grid'
            stop
        endif

!....................


      case('RAD0')
    
        if (grid%is_structured) then  ! look for processor decomposition
          call fiReadDouble(string,grid%structured_grid%Radius_0,ierr)
          call fiDefaultMsg('R_0',ierr)
        else
          if (solution%myrank == 0) &
            print *, 'ERROR: Keyword "RAD0" not supported for unstructured grid'
            stop
        endif


      case ('DIFF')

        call fiReadStringErrorMsg('DIFF',ierr)

        call fiReadDouble(string,solution%difaq,ierr)
        call fiDefaultMsg('difaq',ierr)

        call fiReadDouble(string,solution%delhaq,ierr)
        call fiDefaultMsg('delhaq',ierr)

        if (solution%myrank==0) write(IUNIT2,'(/," *DIFF ",/, &
          &"  difaq       = ",1pe12.4,"[m^2/s]",/, &
          &"  delhaq      = ",1pe12.4,"[kJ/mol]")') &
          solution%difaq,solution%delhaq

!....................

      case ('RCTR')

        call fiReadStringErrorMsg('RCTR',ierr)

        call fiReadInt(string,solution%ityprxn,ierr)
        call fiDefaultMsg('ityprxn',ierr)

        call fiReadDouble(string,solution%rk,ierr)
        call fiDefaultMsg('rk',ierr)

        call fiReadDouble(string,solution%phis0,ierr)
        call fiDefaultMsg('phis0',ierr)

        call fiReadDouble(string,solution%areas0,ierr)
        call fiDefaultMsg('areas0',ierr)

        call fiReadDouble(string,solution%pwrsrf,ierr)
        call fiDefaultMsg('pwrsrf',ierr)

        call fiReadDouble(string,solution%vbars,ierr)
        call fiDefaultMsg('vbars',ierr)

        call fiReadDouble(string,solution%ceq,ierr)
        call fiDefaultMsg('ceq',ierr)

        call fiReadDouble(string,solution%delHs,ierr)
        call fiDefaultMsg('delHs',ierr)

        call fiReadDouble(string,solution%delEs,ierr)
        call fiDefaultMsg('delEs',ierr)

        call fiReadDouble(string,solution%wfmts,ierr)
        call fiDefaultMsg('wfmts',ierr)

        if (solution%myrank == 0) &
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
          & )') solution%ityprxn,solution%rk,solution%phis0,solution%areas0,solution%pwrsrf, &
          solution%vbars,solution%ceq,solution%delHs,solution%delEs,solution%wfmts

 ! convert: mol/cm^2 -> mol/cm^3 -> mol/dm^3 (note area 1/cm)          
        solution%rk = solution%rk * solution%areas0 * 1.d3
        solution%vbars = solution%vbars * 1.d-3 ! convert: cm^3/mol -> L/mol
      
        solution%delHs = solution%delHs * solution%wfmts * 1.d-3 ! convert kJ/kg -> kJ/mol
!        solution%delHs = solution%delHs * solution%scale ! convert J/kmol -> MJ/kmol

!....................

      case ('RADN')

        call fiReadStringErrorMsg('RADN',ierr)

        call fiReadDouble(string,solution%ret,ierr)
        call fiDefaultMsg('ret',ierr)

        call fiReadDouble(string,solution%fc,ierr)
        call fiDefaultMsg('fc',ierr)

        if (solution%myrank==0) write(IUNIT2,'(/," *RADN ",/, &
          &"  ret     = ",1pe12.4,/, &
          &"  fc      = ",1pe12.4)') &
          solution%ret,solution%fc

!....................


      case ('PHAR')

        call fiReadStringErrorMsg('PHAR',ierr)

        call fiReadDouble(string,solution%qu_kin,ierr)
        call fiDefaultMsg('TransReaction',ierr)
        if (solution%myrank==0) write(IUNIT2,'(/," *PHAR ",1pe12.4)')solution%qu_kin
        solution%yh2o_in_co2 = 0.d0
        if (solution%qu_kin > 0.d0) solution%yh2o_in_co2 = 1.d-2 ! check this number!
     
!......................

      case('RICH')
        call fiReadStringErrorMsg('RICH',ierr)
        call fiReadDouble(string,solution%pref,ierr)
        call fiDefaultMsg('Ref. Pressure',ierr) 

!......................

      case('BRIN')
        call fiReadStringErrorMsg('BRIN',ierr)
        call fiReadDouble(string,solution%m_nacl,ierr)
        call fiDefaultMsg('NaCl Concentration',ierr) 

        call fiReadWord(string,word,.false.,ierr)
        call fiWordToUpper(word)
        select case(word(1:len_trim(word)))
          case('MOLAL')
          case('MASS')
            solution%m_nacl = solution%m_nacl /fmwnacl/(1.D0-solution%m_nacl)
          case('MOLE')    
            solution%m_nacl = solution%m_nacl /fmwh2o/(1.D0-solution%m_nacl)
          case default
            print *, 'Wrong unit: ', word(1:len_trim(word))
            stop
         end select 
         if (solution%myrank == 0) print *, solution%m_nacl
!......................

      case ('HYDR')

        call fiReadStringErrorMsg('HYDR',ierr)
  
        call fiReadInt(string,solution%ihydrostatic,ierr)
        call fiDefaultMsg('ihydrostatic',ierr)
       
        call fiReadDouble(string,solution%dTdz,ierr)
        call fiDefaultMsg('dTdz',ierr)

        call fiReadDouble(string,solution%beta,ierr)
        call fiDefaultMsg('beta',ierr)

        call fiReadDouble(string,solution%tref,ierr)
        call fiDefaultMsg('tref',ierr)

        call fiReadDouble(string,solution%pref,ierr)
        call fiDefaultMsg('pref',ierr)

        call fiReadDouble(string,solution%conc0,ierr)
        call fiDefaultMsg('conc0',ierr)

        if (solution%ihydrostatic < 1) solution%ihydrostatic =1
      
        if (solution%myrank==0) write(IUNIT2,'(/," *HYDR ",/, &
          &"  ihydro      = ",i3,/, &
          &"  dT/dz       = ",1pe12.4,/, &
          &"  beta        = ",1pe12.4,/, &
          &"  tref        = ",1pe12.4,/, &
          &"  pref        = ",1pe12.4,/, &
          &"  conc        = ",1pe12.4 &
          &)') &
          solution%ihydrostatic,solution%dTdz,solution%beta,solution%tref,solution%pref, &
          solution%conc0

!....................

      case ('SOLV')
    
        solver => solution%stepper%solver
        call fiReadStringErrorMsg('SOLV',ierr)

!       call fiReadDouble(string,eps,ierr)
!       call fiDefaultMsg('eps',ierr)

        call fiReadDouble(string,solver%atol,ierr)
        call fiDefaultMsg('atol_petsc',ierr)

        call fiReadDouble(string,solver%rtol,ierr)
        call fiDefaultMsg('rtol_petsc',ierr)

        call fiReadDouble(string,solver%stol,ierr)
        call fiDefaultMsg('stol_petsc',ierr)
      
        solver%dtol=1.D5
!       if (solution%use_ksp == 1) then
        call fiReadDouble(string,solver%dtol,ierr)
        call fiDefaultMsg('dtol_petsc',ierr)
!       endif
   
        call fiReadInt(string,solver%maxit,ierr)
        call fiDefaultMsg('maxit',ierr)
      
        call fiReadInt(string,solver%maxf,ierr)
        call fiDefaultMsg('maxf',ierr)
       
        call fiReadInt(string,solver%idt_switch,ierr)
        call fiDefaultMsg('idt',ierr)
        
        solver%inf_tol = solver%atol
        call fiReadDouble(string,solver%inf_tol,ierr)
        call fiDefaultMsg('inf_tol_pflow',ierr)
 
        if (solution%myrank==0) write(IUNIT2,'(/," *SOLV ",/, &
          &"  atol_petsc   = ",1pe12.4,/, &
          &"  rtol_petsc   = ",1pe12.4,/, &
          &"  stol_petsc   = ",1pe12.4,/, &
          &"  dtol_petsc   = ",1pe12.4,/, &
          &"  maxit        = ",8x,i5,/, &
          &"  maxf        = ",8x,i5,/, &
          &"  idt         = ",8x,i5 &
          &    )') &
           solver%atol,solver%rtol,solver%stol,solver%dtol,solver%maxit, &
           solver%maxf,solver%idt_switch

! The line below is a commented-out portion of the format string above.
! We have to put it here because of the stupid Sun compiler.
!    &"  eps          = ",1pe12.4,/, &

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

          call fiReadDouble(string,material%rock_density(ireg),ierr)
          call fiDefaultMsg('rock_density',ierr)

          call fiReadDouble(string,material%cpr(ireg),ierr)
          call fiDefaultMsg('cpr',ierr)
        
          call fiReadDouble(string,material%ckdry(ireg),ierr)
          call fiDefaultMsg('ckdry',ierr)
        
          call fiReadDouble(string,material%ckwet(ireg),ierr)
          call fiDefaultMsg('ckwet',ierr)
        
          call fiReadDouble(string,material%tau(ireg),ierr)
          call fiDefaultMsg('tau',ierr)

          call fiReadDouble(string,material%cdiff(ireg),ierr)
          call fiDefaultMsg('cdiff',ierr)

          call fiReadDouble(string,material%cexp(ireg),ierr)
          call fiDefaultMsg('cexp',ierr)

        !scale thermal properties
          material%cpr(ireg) = solution%scale * material%cpr(ireg)
          material%dencpr(ireg) = material%rock_density(ireg) * material%cpr(ireg)
          material%ckdry(ireg) = solution%scale * material%ckdry(ireg)
          material%ckwet(ireg) = solution%scale * material%ckwet(ireg)
        enddo
      
        if (solution%myrank==0) then
          write(IUNIT2,'(/," *THRM: ",i3)') ireg
          write(IUNIT2,'("  itm rock_density  cpr        ckdry", &
            &                 "     ckwet       tau       cdiff     cexp")')
          write(IUNIT2,'("        [kg/m^3]  [J/kg/K]   [J/m/K/s]", &
            &              "     [J/m/K/s]     [-]        [m^2/s]       [-]")')
          do i = 1, ireg
            write(IUNIT2,'(i4,1p7e11.4)') i,material%rock_density(i), &
            material%cpr(i),material%ckdry(i),material%ckwet(i), &
            material%tau(i),material%cdiff(i),material%cexp(i)
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
          
          call fiReadInt(string,material%icaptype(idum),ierr)
          call fiDefaultMsg('icaptype',ierr)
      
          if (solution%use_mph == PETSC_TRUE .or. solution%use_owg == PETSC_TRUE &
              .or. solution%use_vadose == PETSC_TRUE .or. solution%use_flash == PETSC_TRUE&
              .or. solution%use_richards == PETSC_TRUE) then
            do np=1, grid%nphase
              call fiReadDouble(string,material%sir(np,idum),ierr)
              call fiDefaultMsg('sir',ierr)
            enddo 
          else
            call fiReadDouble(string,material%swir(idum),ierr)
            call fiDefaultMsg('swir',ierr)
          endif
        
          call fiReadDouble(string,material%pckrm(idum),ierr)
          call fiDefaultMsg('lambda',ierr)
          material%lambda(idum) = material%pckrm(idum)/(-material%pckrm(idum) +1.D0)
! Here the lambda is assigned as the same value of m

          call fiReadDouble(string,material%alpha(idum),ierr)
          call fiDefaultMsg('alpha',ierr)

          call fiReadDouble(string,material%pcwmax(idum),ierr)
          call fiDefaultMsg('pcwmax',ierr)
      
          call fiReadDouble(string,material%pcbetac(idum),ierr)
          call fiDefaultMsg('pcbetac',ierr)
      
          call fiReadDouble(string,material%pwrprm(idum),ierr)
          call fiDefaultMsg('pwrprm',ierr)

        enddo

          if (solution%use_mph == PETSC_TRUE .or. &
              solution%use_vadose == PETSC_TRUE .or. &
              solution%use_flash == PETSC_TRUE .or. &
              solution%use_richards == PETSC_TRUE) then
              call pckr_init(solution%nphase,ireg,solution%nlmax, &
                             material%icaptype,material%sir, material%pckrm, &
                             material%lambda,material%alpha,material%pcwmax, &
                             material%pcbetac,material%pwrprm)
          endif 

      
        if (solution%myrank==0) then
          write(IUNIT2,'(/," *PCKR: ",i3)') ireg
          write(IUNIT2,'("  icp swir    lambda         alpha")')
          do j = 1, ireg
            i=material%icaptype(j)
            if (solution%use_mph==PETSC_TRUE .or. &
                solution%use_owg==PETSC_TRUE .or. &
                solution%use_vadose == PETSC_TRUE .or. &
                solution%use_flash == PETSC_TRUE .or. &
                solution%use_richards == PETSC_TRUE) then
              write(IUNIT2,'(i4,1p8e12.4)') i,(material%sir(np,i),np=1, &
                material%nphase),material%lambda(i),material%alpha(i), &
                material%pcwmax(i),material%pcbetac(i),material%pwrprm(i)
            else
              write(IUNIT2,'(i4,1p7e12.4)') i,material%swir(i), &
                material%lambda(i),material%alpha(i),material%pcwmax(i), &
                material%pcbetac(i),material%pwrprm(i)
            endif
          enddo
        end if

        if (solution%use_mph == PETSC_TRUE .or. &
            solution%use_vadose == PETSC_TRUE .or. &
            solution%use_flash == PETSC_TRUE .or. &
            solution%use_richards == PETSC_TRUE) then
          deallocate(material%icaptype, material%pckrm, material%lambda, &
                     material%alpha,material%pcwmax, material%pcbetac, &
                     material%pwrprm)
        endif 
 

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
      
!GEH - Structured Grid Dependence - Begin
          call readRegion(string)
          word = adjustl(string(1:min(31,len_trim(string))))  ! remove leading spaces
          call fiCharsToLower(string,4)
          if (fiStringCompare(string,'file',4)) then ! unstructured
            id = readUnstructuredRegion()
          else
            id = readStructuredRegion(string)
          endif
    
        enddo
        material%iregperm = ireg
        
        if (ireg > MAXINITREGIONS) then
          print *,'error: increase MAXINITREGIONS: regions ',ireg,' > ',MAXINITREGIONS
          stop
        endif

        if (solution%myrank==0) then
          write(IUNIT2,'(/," *PHIK: ireg = ",i4)') material%iregperm
          write(IUNIT2,'("  i1  i2  j1  j2  k1  k2 icap ithrm  por      tor  &
            &",   "     permx      permy      permz [m^2]   permpwr")')
          do ireg = 1, grid%iregperm
!GEH - Structured Grid Dependence - Begin
            write(IUNIT2,'(6i4,2i4,1p6e11.4)') &
                  material%i1reg(ireg),material%i2reg(ireg), &
                  material%j1reg(ireg),material%j2reg(ireg), &
                  material%k1reg(ireg),material%k2reg(ireg), &
                  material%icap_reg(ireg),material%ithrm_reg(ireg), &
                  material%por_reg(ireg),material%tor_reg(ireg), &
                  (material%perm_reg(ireg,i),i=1,4)
!GEH - Structured Grid Dependence - End
          enddo
        endif

!....................
      
      case ('INIT')
    
        call fiReadInt(string,solution%iread_init,ierr) 
        call fiDefaultMsg('iread_init',ierr)
      
        if (solution%myrank==0) then
          write(IUNIT2,'(/," *INIT: iread = ",i2)') grid%iread_init
        endif
      
        if (grid%iread_init == 0 .or. grid%iread_init == 2) then
      
          ireg = 0
          do
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('INIT',ierr)

            if (string(1:1) == '.' .or. string(1:1) == '/') exit
            ireg = ireg + 1
            
!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

            if (solution%use_mph==PETSC_TRUE .or. &
                solution%use_owg==PETSC_TRUE .or. &
                solution%use_vadose == PETSC_TRUE .or. &
                solution%use_flash == PETSC_TRUE .or. &
                solution%use_richards == PETSC_TRUE) then
              call fiReadInt(string,solution%iphas_ini(ireg),ierr)
              call fiDefaultMsg('iphase',ierr)
         
              do j=1,grid%ndof
                call fiReadDouble(string,solution%xx_ini(j,ireg),ierr)
                call fiDefaultMsg('xxini',ierr)
              enddo
            else
              call fiReadDouble(string,solution%pres_ini(ireg),ierr)
              call fiDefaultMsg('pres',ierr)
  
              call fiReadDouble(string,solution%temp_ini(ireg),ierr)
              call fiDefaultMsg('temp',ierr)
  
              call fiReadDouble(string,solution%sat_ini(ireg),ierr)
              call fiDefaultMsg('sat',ierr)
!              grid%sat_ini(ireg)=1.D0 - grid%sat_ini(ireg)
  
              call fiReadDouble(string,solution%conc_ini(ireg),ierr)
              call fiDefaultMsg('conc',ierr)
            endif
          enddo
      
          grid%iregini = ireg
      
          if (solution%myrank==0) then
            write(IUNIT2,'("  ireg = ",i4)') grid%iregini
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]   &
              &   ",    "sl [-]      c [mol/L]")')
            do ireg = 1, grid%iregini
!GEH - Structured Grid Dependence - Begin
              if (solution%use_mph==PETSC_TRUE .or. solution%use_owg==PETSC_TRUE &
                  .or. solution%use_vadose == PETSC_TRUE .or. solution%use_flash == PETSC_TRUE&
                  .or. solution%use_richards == PETSC_TRUE) then
                write(IUNIT2,'(7i4,1p10e12.4)') &
                  grid%i1ini(ireg),grid%i2ini(ireg), &
                  grid%j1ini(ireg),grid%j2ini(ireg), &
                  grid%k1ini(ireg),grid%k2ini(ireg), &
                  grid%iphas_ini(ireg),(grid%xx_ini(np,ireg),np =1,grid%ndof)
              else
                write(IUNIT2,'(6i4,1p10e12.4)') &
                  grid%i1ini(ireg),grid%i2ini(ireg), &
                  grid%j1ini(ireg),grid%j2ini(ireg), &
                  grid%k1ini(ireg),grid%k2ini(ireg), &
                  grid%pres_ini(ireg),grid%temp_ini(ireg),grid%sat_ini(ireg), &
                  grid%conc_ini(ireg)
              endif
!GEH - Structured Grid Dependence - End
            enddo
          endif

        else if (grid%iread_init == 1) then
    
!     read in initial conditions from file: pflow_init.dat
          if (solution%myrank == 0) then
            write(*,*) '--> read in initial conditions from file: &
                        &pflow_init.dat'
  
            open(IUNIT3, file='pflow_init.dat', action="read", status="old")

            ireg = 0
            do
              call fiReadFlotranString(IUNIT3,string,ierr)
!             call fiReadStringErrorMsg('INIT',ierr)

              if (string(1:1) == '.' .or. string(1:1) == '/') exit
              ireg = ireg + 1

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

              if (solution%use_mph==PETSC_TRUE .or. solution%use_owg==PETSC_TRUE &
                  .or. solution%use_vadose == PETSC_TRUE .or. solution%use_flash == PETSC_TRUE&
                  .or. solution%use_richards == PETSC_TRUE) then
                call fiReadInt(string,grid%iphas_ini(ireg),ierr)
                call fiDefaultMsg('iphase_ini',ierr)
            
                do j=1,grid%ndof
                  call fiReadDouble(string,grid%xx_ini(j,ireg),ierr)
                  call fiDefaultMsg('xx_ini',ierr)
                enddo
              else
  
                call fiReadDouble(string,grid%pres_ini(ireg),ierr)
                call fiDefaultMsg('pres',ierr)
  
                call fiReadDouble(string,grid%temp_ini(ireg),ierr)
                call fiDefaultMsg('temp',ierr)
  
                call fiReadDouble(string,grid%sat_ini(ireg),ierr)
                call fiDefaultMsg('sat',ierr)

                call fiReadDouble(string,grid%conc_ini(ireg),ierr)
                call fiDefaultMsg('conc',ierr)
              endif
       
            enddo
            grid%iregini = ireg
            close(IUNIT3)
          endif
        endif

!....................

      case ('TIME')

        call fiReadStringErrorMsg('TIME',ierr)
      
        call fiReadWord(string,strtim,.false.,ierr)
      
        grid%tunit = strtim

        if (grid%tunit == 's') then
          grid%tconv = 1.d0
        else if (grid%tunit == 'm') then
          grid%tconv = 60.d0
        else if (grid%tunit == 'h') then
          grid%tconv = 60.d0 * 60.d0
        else if (grid%tunit == 'd') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0
        else if (grid%tunit == 'mo') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0 * 30.d0
        else if (grid%tunit == 'y') then
          grid%tconv = 60.d0 * 60.d0 * 24.d0 * 365.d0
        else
          if (solution%myrank == 0) then
            write(*,'(" Time unit: ",a3,/, &
              &" Error: time units must be one of ",/, &
              &"   s -seconds",/,"   m -minutes",/,"   h -hours",/, &
              &"   d -days", /, "  mo -months",/,"   y -years")') grid%tunit
          endif
          stop
        endif

        call fiReadInt(string,grid%kplot,ierr) 
        call fiDefaultMsg('kplot',ierr)
      
        allocate(grid%tplot(grid%kplot))
      
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('TIME',ierr)
        i2 = 0
        do
          i1 = i2 + 1
          i2 = i2+10
          if (i2 > grid%kplot) i2 = grid%kplot
          do i = i1, i2
            call fiReadDouble(string,grid%tplot(i),ierr)
            call fiDefaultMsg('tplot',ierr)
          enddo
          if (i2 == grid%kplot) exit
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('TIME',ierr)
        enddo

!       call fiReadFlotranString(IUNIT1,string,ierr)
!       call fiReadStringErrorMsg('TIME',ierr)
      
!       call fiReadDouble(string,grid%dt,ierr)
!       call fiDefaultMsg('dt',ierr)

!       call fiReadDouble(string,grid%dt_max,ierr)
!       call fiDefaultMsg('dt_max',ierr)

        if (solution%myrank==0) then
          write(IUNIT2,'(/," *TIME ",a3,1x,i4,/,(1p10e12.4))') grid%tunit, &
          grid%kplot,(grid%tplot(i),i=1,grid%kplot)
!         write(IUNIT2,'("  dt= ",1pe12.4,", dtmax= ",1pe12.4,/)') &
!         grid%dt,grid%dt_max
        endif
      
        ! convert time units to seconds
        do i = 1, grid%kplot
          grid%tplot(i) = grid%tconv * grid%tplot(i)
        enddo
!       grid%dt = grid%tconv * grid%dt
!       grid%dt_max = grid%tconv * grid%dt_max

!....................

      case ('DTST')

        call fiReadStringErrorMsg('DTST',ierr)
  
        call fiReadInt(string,grid%nstpmax,ierr)
        call fiDefaultMsg('nstpmax',ierr)
  
        allocate(grid%tstep(grid%nstpmax))
        allocate(grid%dtstep(grid%nstpmax))
  
        do i = 1, grid%nstpmax
          call fiReadDouble(string,grid%tstep(i),ierr)
          call fiDefaultMsg('tstep',ierr)
        enddo

        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg('DTST',ierr)
        call fiReadDouble(string,grid%dt_min,ierr)
        call fiDefaultMsg('dt_min',ierr)
        do i = 1, grid%nstpmax
          call fiReadDouble(string,grid%dtstep(i),ierr)
          call fiDefaultMsg('dtstep',ierr)
        enddo
        
        grid%dt_max = grid%dtstep(1)
        
        grid%dt = grid%dt_min
      
        if (solution%myrank==0) then
          write(IUNIT2,'(/," *DTST ",i4,/," tstep= ",(1p10e12.4))')  &
            grid%nstpmax, (grid%tstep(i),i=1,grid%nstpmax)
          write(IUNIT2,'(" dtstep= ",1p10e12.4,/)') &
            grid%dt_min,(grid%dtstep(i),i=1,grid%nstpmax)
        endif
      
        ! convert time units to seconds
        do i = 1, grid%nstpmax
          grid%tstep(i) = grid%tconv * grid%tstep(i)
          grid%dtstep(i) = grid%tconv * grid%dtstep(i)
        enddo
        grid%dt = grid%tconv * grid%dt
        grid%dt_min = grid%tconv * grid%dt_min
        grid%dt_max = grid%tconv * grid%dt_max

!....................

      case ('BCON')

!-----------------------------------------------------------------------
!-----boundary conditions:  ibnd:  
!                   1-left,    2-right
!          3-top,    4-bottom
!          5-front,  6-back
!
!  ibndtyp:  1-Dirichlet         (p, T, C)
!  ibndtyp:  2-Neumann/Dirichlet (q, grad T=0, grad C=0)
!  ibndtyp:  3-Dirichlet/Neumann (p, grad T=0, grad C=0)
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

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

            ! Now read the velocities or pressures, depending on the BC type
            call fiReadFlotranString(IUNIT1,string,ierr)
            call fiReadStringErrorMsg('BCON',ierr)
   
            if (solution%use_mph /= PETSC_TRUE .and.  &
                solution%use_owg /= PETSC_TRUE .and. &
                solution%use_vadose /= PETSC_TRUE .and. &
                solution%use_flash /= PETSC_TRUE .and. &
                solution%use_richards /= PETSC_TRUE) then  
              j=1
              if (grid%nphase>1) j=2
              if (grid%ibndtyp(ibc) == 1) then 
                call fiReadDouble(string, grid%pressurebc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading pressure BCs:", ierr)
              else if (grid%ibndtyp(ibc) == 2) then
                call fiReadDouble(string, grid%velocitybc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading velocity BCs:", ierr)
              else if (grid%ibndtyp(ibc) == 4) then
                call fiReadDouble(string, grid%velocitybc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading velocity BCs:", ierr)  
              else
                call fiReadDouble(string, grid%pressurebc0(j,ibc), ierr)
                call fiDefaultMsg("Error reading pressure BCs:", ierr)
              endif

              if (grid%nphase>1) grid%pressurebc0(1,ibc) = &
                  grid%pressurebc0(2,ibc)
                ! For simple input
              call fiReadDouble(string,grid%tempbc0(ibc),ierr)
              call fiDefaultMsg('tempbc',ierr)

              call fiReadDouble(string,grid%sgbc0(ibc),ierr)
              call fiDefaultMsg('sgbc',ierr)
              grid%sgbc0(ibc) = 1.D0 - grid%sgbc0(ibc) ! read in sl

              call fiReadDouble(string,grid%concbc0(ibc),ierr)
              call fiDefaultMsg('concbc',ierr)
      
            else
     
              call fiReadInt(string,grid%iphasebc0(ir),ierr)
              call fiDefaultMsg('iphase',ierr)
       
              if (grid%ibndtyp(ir) == 1 .or. grid%ibndtyp(ir) == 3 &
                  .or. grid%ibndtyp(ir) == 4) then 
                do j=1,grid%ndof
                  call fiReadDouble(string,grid%xxbc0(j,ir),ierr)
                  call fiDefaultMsg('xxbc',ierr)
                enddo
              elseif (grid%ibndtyp(ir) == 2) then
                do j=1, grid%nphase       
                  call fiReadDouble(string, grid%velocitybc0(j,ir), ierr)
                  call fiDefaultMsg("Error reading velocity BCs:", ierr)
                enddo
                do j=2,grid%ndof
                  call fiReadDouble(string,grid%xxbc0(j,ir),ierr)
                  call fiDefaultMsg('xxbc',ierr)
                enddo
              endif
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

!GEH - Structured Grid Dependence - Begin
        if (solution%myrank == 0) then
          write(IUNIT2,'(/," *BCON: nblkbc = ",i4)') grid%nblkbc
          do ibc = 1, grid%nblkbc
            write(IUNIT2,'("  ibndtyp = ",i3," iface = ",i2)') &
              grid%ibndtyp(ibc), grid%iface(ibc)
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2       p [Pa]     t [C]  &
              &  c",     " [mol/L]")')
            do ireg = grid%iregbc1(ibc), grid%iregbc2(ibc)
              if (solution%use_mph== PETSC_TRUE .or. solution%use_owg== PETSC_TRUE &
                  .or. solution%use_vadose == PETSC_TRUE .or. solution%use_flash == PETSC_TRUE&
                  .or. solution%use_richards == PETSC_TRUE) then
                if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(7i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    grid%iphasebc0(ireg),(grid%xxbc0(j,ireg),j=1,grid%ndof)
                else if (grid%ibndtyp(ibc) == 2 .or. grid%ibndtyp(ibc) == 4) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%velocitybc0(j,ireg),j=1,grid%nphase),&
                    (grid%xxbc0(j,ireg),j=2,grid%ndof)
                endif
              else
                if (grid%ibndtyp(ibc) == 1 .or. grid%ibndtyp(ibc) == 3) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%pressurebc0(j,ireg),j=1,grid%nphase), &
                    grid%tempbc0(ireg), &
                    grid%concbc0(ireg)
                else if (grid%ibndtyp(ibc) == 2 .or. grid%ibndtyp(ibc) == 4) then
                  write(IUNIT2,'(6i4,1p10e12.4)') &
                    grid%i1bc(ireg),grid%i2bc(ireg), &
                    grid%j1bc(ireg),grid%j2bc(ireg), &
                    grid%k1bc(ireg),grid%k2bc(ireg), &
                    (grid%velocitybc0(j,ireg),j=1,grid%nphase), &
                    grid%tempbc0(ireg), &
                    grid%concbc0(ireg)
                endif
              endif
            enddo
          enddo
        endif
!GEH - Structured Grid Dependence - End

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

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

!         print *,'pflowgrid_mod: Source', isrc, ir   
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
      
            call fiReadDouble(string,grid%qsrc(i,isrc),ierr)
            call fiDefaultMsg('qsrc',ierr)
      
            call fiReadDouble(string,grid%csrc(i,isrc),ierr)
            call fiDefaultMsg('csrc',ierr)
         
            call fiReadDouble(string,grid%hsrc(i,isrc),ierr)
            call fiDefaultMsg('hsrc',ierr)


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

!GEH - Structured Grid Dependence - Begin
        if (solution%myrank == 0) then
          write(IUNIT2,'(/," *SOURce: nblksrc = ",i4)') grid%nblksrc
          do isrc = 1, grid%nblksrc
            write(IUNIT2,'("  i1  i2  j1  j2  k1  k2")')
            write(IUNIT2,'(6i4)') &
              grid%i1src(isrc),grid%i2src(isrc), &
              grid%j1src(isrc),grid%j2src(isrc), &
              grid%k1src(isrc),grid%k2src(isrc)
            write(IUNIT2,'("    t [s]        T [C]    QH2O [kg/s]    &
              &QCO2 [kg/s]")')
            do ir = 1, grid%ntimsrc
              write(IUNIT2,'(1p10e12.4)') &
                grid%timesrc(ir,isrc),grid%tempsrc(ir,isrc), &
                grid%qsrc(ir,isrc), &
                grid%csrc(ir,isrc)
            enddo
          enddo
        endif
!GEH - Structured Grid Dependence - End

!....................
      
      case ('BRK')

        ibrk = 0
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('BRK',ierr)
      
          if (string(1:1) == '.' .or. string(1:1) == '/') exit
          ibrk = ibrk + 1

!GEH - Structured Grid Dependence - Begin
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
!GEH - Structured Grid Dependence - End

          call fiReadInt(string,grid%ibrktyp(ibrk),ierr)
          call fiDefaultMsg('ibrktyp',ierr)

          call fiReadInt(string,grid%ibrkface(ibrk),ierr)
          call fiDefaultMsg('ibrkface',ierr)
        enddo
        grid%ibrkcrv = ibrk
            
        if (solution%myrank==0) then
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
      case('SDST')
         
          
          do j=1,grid%ndof
            call fiReadDouble(string,grid%steady_eps(j),ierr)
            call fiDefaultMsg('steady tol',ierr)
          enddo
        if (solution%myrank==0) write(IUNIT2,'(/," *SDST ",/, &
          &"  dpdt        = ",1pe12.4,/, &
          &"  dtmpdt        = ",1pe12.4,/, &
          &"  dcdt        = ",1pe12.4)') &
          grid%steady_eps

!....................
      case default
    
        if (solution%myrank == 0) then
          print *, "Error reading input file: keyword not found. Terminating."
        endif
        call PetscFinalize(ierr)
        stop

    end select

  enddo

  close(IUNIT1)
  
  call createPflowDMs(grid)

end subroutine initPFLOW

subroutine createPflowDMs(grid)

 implicit none
 
 type(pflowGrid) :: grid
 
 if (structured) then
   call createStructuredDMs(grid)
 else
 endif
 
end subroutine createPflowDMs

end module Init_module
