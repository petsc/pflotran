module pflow_read_gridsize_module

contains

  subroutine pflow_read_gridsize(inputfile, igeom, nx, ny, nz, npx, npy, npz, &
  nphase, nspec, npricomp, ndof, idcdm, itable, ierr)
  
  use fileio_module
  
  implicit none
  
  character(len=*), intent(in) :: inputfile
  integer, intent(out) :: igeom, nx, ny, nz, npx, npy, npz, nphase, ndof
  integer, intent(out) :: nspec,npricomp,idcdm,itable,ierr

#include "include/finclude/petsc.h"
#include "include/finclude/petscdef.h"

#include "definitions.h"
  integer :: myrank, gridread_flag, commsize
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXCARDLENGTH) :: card

  call MPI_Comm_rank(PETSC_COMM_WORLD, myrank, ierr)

  open(IUNIT1, file=inputfile, action="read", status="old") 
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

      call fiReadInt(string,nspec,ierr)
      call fiDefaultMsg('nspec',ierr)

      call fiReadInt(string,npricomp,ierr)
      call fiDefaultMsg('npricomp',ierr)

      call fiReadInt(string,ndof,ierr)
      call fiDefaultMsg('ndof',ierr)
      
      call fiReadInt(string,idcdm,ierr)
      call fiDefaultMsg('idcdm',ierr)

      call fiReadInt(string,itable,ierr)
      call fiDefaultMsg('itable',ierr)

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

end module pflow_read_gridsize_module
