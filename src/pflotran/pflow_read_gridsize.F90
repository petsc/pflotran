module pflow_read_gridsize_module

  implicit none
  private

#include "definitions.h"

  public :: pflow_read_gridsize

contains

  subroutine pflow_read_gridsize(inputfile, igeom, nx, ny, nz, npx, npy, npz, &
  nphase, nspec, npricomp, ndof, idcdm, itable, mcomp,mphas ,ierr)
  
  use fileio_module
  
  implicit none
  
#include "include/finclude/petscdef.h"

  character(len=*), intent(in) :: inputfile

  PetscInt, intent(out) :: igeom, nx, ny, nz, npx, npy, npz, nphase, ndof
  
  PetscInt, intent(out) :: nspec,npricomp,idcdm,itable
  PetscErrorCode :: ierr
  PetscInt, intent(out) :: mcomp, mphas
 
  PetscInt :: gridread_flag, length
  PetscMPIInt :: myrank, commsize
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXCARDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: name
  
  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  open(IUNIT1, file=inputfile, action="read", status="old") 
  open(IUNIT2, file='pflow.out', action="write", status="unknown")

  npx = PETSC_DECIDE; npy = PETSC_DECIDE; npz = PETSC_DECIDE
  ierr = 1
  gridread_flag = 0
  mphas=0; mcomp = 0


  do
!   if(gridread_flag /= 0) exit
    call fiReadFlotranString(IUNIT1, string, ierr)
    if(ierr /= 0) exit

    call fiReadWord(string,word,.false.,ierr)
    length = len_trim(word)
    call fiCharsToUpper(word,length)
    call fiReadCard(word,card,ierr)

    select case(card)

!...............................................................

    case ('GRID')

      if (myrank == 0) print *, card
      call fiReadStringErrorMsg(myrank,'GRID',ierr)
      
      call fiReadInt(string,igeom,ierr)
      call fiDefaultMsg(myrank,'igeom',ierr)
      
      call fiReadInt(string,nx,ierr)
      call fiDefaultMsg(myrank,'nx',ierr)
      
      call fiReadInt(string,ny,ierr)
      call fiDefaultMsg(myrank,'ny',ierr)
      
      call fiReadInt(string,nz,ierr)
      call fiDefaultMsg(myrank,'nz',ierr)
      
      call fiReadInt(string,nphase,ierr)
      call fiDefaultMsg(myrank,'nphase',ierr)

      call fiReadInt(string,nspec,ierr)
      call fiDefaultMsg(myrank,'nspec',ierr)

      call fiReadInt(string,npricomp,ierr)
      call fiDefaultMsg(myrank,'npricomp',ierr)

      call fiReadInt(string,ndof,ierr)
      call fiDefaultMsg(myrank,'ndof',ierr)
      
      call fiReadInt(string,idcdm,ierr)
      call fiDefaultMsg(myrank,'idcdm',ierr)

      call fiReadInt(string,itable,ierr)
      call fiDefaultMsg(myrank,'itable',ierr)

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
     &   )') igeom, nx, ny, nz, nphase, ndof, idcdm, itable

    case ('PROC')
    
      if (myrank == 0) print *, card

      call fiReadStringErrorMsg(myrank,'PROC',ierr)

      call fiReadInt(string,npx,ierr)
      call fiDefaultMsg(myrank,'npx',ierr)
      call fiReadInt(string,npy,ierr)
      call fiDefaultMsg(myrank,'npy',ierr)
      call fiReadInt(string,npz,ierr)
      call fiDefaultMsg(myrank,'npz',ierr)
 
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

!.........................................................................

     case ('COMP')
    
      mcomp =0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg(myrank,'COMP',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadWord(string,name,.true.,ierr)
        call fiErrorMsg(myrank,'namcx','GAS',ierr)
        
        call fiWordToUpper(name) 
        length = len_trim(name)
        select case(name(1:length))
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
!     print *,'comp: ',mcomp   
          
!....................................................................

     case ('PHAS')
    
      mphas =0
      do
        call fiReadFlotranString(IUNIT1,string,ierr)
        call fiReadStringErrorMsg(myrank,'phase',ierr)
      
        if (string(1:1) == '.' .or. string(1:1) == '/') exit

        call fiReadWord(string,name,.true.,ierr)
        call fiErrorMsg(myrank,'namcx','phase',ierr)
        
        call fiWordToUpper(name) 
        length = len_trim(name) 
        select case(name(1:length))
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
!     print *,'phase: ',mphas     
!....................................................................





    case default
!     print *, "Entered case default"
      ! Skip anything else that isn't a GRID specification.
      exit
    endselect

  enddo
  
  close(IUNIT1)

  end subroutine pflow_read_gridsize

end module pflow_read_gridsize_module
