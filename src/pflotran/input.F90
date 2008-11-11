module Input_module

  use Option_module

  implicit none

  private

#include "definitions.h"

  type, public :: input_type 
    PetscInt :: fid
    PetscErrorCode :: ierr
    character(len=MAXSTRINGLENGTH) :: buf
    character(len=MAXSTRINGLENGTH) :: err_buf
    character(len=MAXSTRINGLENGTH) :: err_buf2
    PetscTruth :: broadcast_read
  end type input_type
  
  public :: InputReadFlotranString

contains

! ************************************************************************** !
!
! InputCreate: Allocates and initializes a new Input object
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
function InputCreate()

  implicit none
  
  type(input_type), pointer :: InputCreate
  
  type(input_type), pointer :: input
  
  allocate(input)

  input%fid = 0
  input%ierr = 0
  input%buf = ''
  input%err_buf = ''
  input%err_buf2 = ''
  input%broadcast_read = PETSC_FALSE

  InputCreate => input
  
end function InputCreate

#if 0
! ************************************************************************** !
!
! InputDefaultMsg: If ierr /= 0, informs user that default value will be used.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputDefaultMsg(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (input%ierr /= 0) then
    if (option%print_rank) &
      print *, '"', adjustl(trim(input%err_buf)), &
               '" set to default value.'
    input%ierr = 0
  endif

end subroutine InputDefaultMsg

! ************************************************************************** !
!
! InputErrorMsg: If ierr /= 0, If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputErrorMsg(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (input%ierr /= 0) then
    if (option%print_rank) &
      print *, 'Error reading "', adjustl(trim(input%err_buf)), &
               '" under keyword: ',adjustl(trim(input%err_buf2)), '.'
    stop
  endif

end subroutine InputErrorMsg

! ************************************************************************** !
!
! fiReadStringErrorMsg: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadStringErrorMsg(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (input%ierr /= 0) then
    if (option%print_rank) &
      print *, 'Error reading in string in (', &
               adjustl(trim(input%err_buf)), ').'
    stop
  endif

end subroutine InputReadStringErrorMsg

! ************************************************************************** !
!
! fiFindStringErrorMsg: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputFindStringErrorMsg(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (input%ierr /= 0) then
    if (option%print_rank) &
      print *, 'Error: Card (', adjustl(trim(input%err_buf)), ') not ', &
               'found in file.'
    stop
  endif

end subroutine InputFindStringErrorMsg

! ************************************************************************** !
!
! InputReadInt: reads and removes an integer value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadInt(input, option, int)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: int

  character(len=MAXWORDLENGTH) :: word

  call fiReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (input%ierr == 0) then
    read(word,*,iostat=input%ierr) int
  endif

end subroutine InputReadInt

! ************************************************************************** !
!
! InputReadDouble: reads and removes a real value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadDouble(input, option, double)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscReal :: double

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (input%ierr == 0) then
    read(word,*,iostat=input%ierr) double
  endif

end subroutine InputReadDouble
#endif
! ************************************************************************** !
!
! InputReadFlotranString: Reads a string (strlen characters long) from a 
!                           file while avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadFlotranString(option, fid, string, ierr)

  implicit none

  type(option_type) :: option
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr, ierr2

  if (option%broadcast_read) then
    if (option%myrank == 0) then
      call fiReadFlotranString(fid, string, ierr)
    endif
    call mpi_bcast(ierr,ONE_INTEGER,MPI_INTEGER,ZERO_INTEGER, &
                   option%comm,ierr2)
    if (ierr == 0) then  
      call mpi_bcast(string,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     ZERO_INTEGER,option%comm,ierr2)      
    endif
  else
    call fiReadFlotranString(fid, string, ierr)
  endif

end subroutine InputReadFlotranString

#if 0
! ************************************************************************** !
!
! InputReadFlotranString2: Reads a string (strlen characters long) from a 
!                           file while avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadFlotranString2(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  
  PetscErrorCode :: ierr2

  if (input%broadcast_read) then
    if (option%myrank == option%io_rank) then
      call fiReadFlotranString(input, option)
    endif
    call mpi_bcast(input%ierr,ONE_INTEGER,MPI_INTEGER,option%io_rank, &
                   option%comm,ierr2)
    if (input%ierr == 0) then  
      call mpi_bcast(input%buf,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     option%io_rank,option%comm,ierr2)      
    endif
  else
    call fiReadFlotranString(input, option)
  endif

end subroutine InputReadFlotranString2

! ************************************************************************** !
!
! InputReadFlotranStringSlave: Reads a string (strlen characters long) from a 
!                              file while avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadFlotranStringSlave(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) ::  tempstring
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i

  input%ierr = 0

! we initialize the word to blanks to avoid error reported by valgrind
  do i=1,strlen
     word(i:i) = ' '
  enddo
  
  do
    read(fid,'(a256)',iostat=input%ierr) input%buf

    if (input%ierr /= 0) exit

    if (input%buf(1:1) == ':' .or. input%buf(1:1) == '!') cycle

    tempstring = input%buf
    call InputReadWord2(tempstring,word,PETSC_TRUE,input%ierr)
    call InputWordToUpper(word)
    if (word(1:4) == 'SKIP') then
      do 
        read(fid,'(a256)',iostat=input%ierr) tempstring
        if (input%ierr /= 0 .and. option%print_rank) then
          print *, 'End of file reached in fiReadFlotranString.'
          print *, 'SKIP encountered without matching NOSKIP.'
        endif
        call fiReadWord(tempstring,word,PETSC_FALSE,input%ierr)
        call fiWordToUpper(word)
        if (word(1:4) == 'NOSK') exit
      enddo
      if (input%ierr /= 0) exit
    else if (word(1:1) /= ' ' .and. word(1:4) /= 'NOSK') then
      exit
    endif
  enddo
  
  ! Check for comment midway along a string
  if (input%ierr == 0) then
    tempstring = string
    do i=1,strlen
      input%buf(i:i) = ' '
    enddo
    do i=1,len_trim(tempstring)
      if (tempstring(i:i) /= ':' .and. tempstring(i:i) /= '!') then
        input%buf(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine InputReadFlotranStringSlave

! ************************************************************************** !
!
! InputReadWord1: reads and removes a word (consecutive characters) from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadWord1(input, option, word, return_blank_error)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  PetscTruth :: return_blank_error
  
  PetscInt :: i, begins, ends

  if (input%ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,len_trim(word)
    word(i:i) = ' '
  enddo

  input%ierr = len_trim(input%buf)
  
  if (input%ierr == 0) then
    if (return_blank_error) then
      input%ierr = 1
    else
      input%ierr = 0
    endif
    return
  else
    input%ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == ',' .or. &
             input%buf(i:i) == achar(9)) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (input%buf(i:i) /= ' ' .and. input%buf(i:i) /= ',' .and. &
              input%buf(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > 31) ends = begins + 31

    ! Copy (ends-begins) characters to 'word'
    word = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(ends+1:)

  endif

end subroutine InputReadWord1

! ************************************************************************** !
!
! InputReadWord2: reads and removes a word (consecutive characters) from a 
!                 string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadWord2(string, word, return_blank_error, ierr)

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscTruth :: return_blank_error
  PetscErrorCode :: ierr
  
  PetscInt :: i, begins, ends

  if (ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,len_trim(word)
    word(i:i) = ' '
  enddo

  ierr = len_trim(input%buf)
  
  if (ierr == 0) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == ',' .or. &
             input%buf(i:i) == achar(9)) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (input%buf(i:i) /= ' ' .and. input%buf(i:i) /= ',' .and. &
              input%buf(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > 31) ends = begins + 31

    ! Copy (ends-begins) characters to 'word'
    word = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(ends+1:)

  endif

end subroutine InputReadWord2

! ************************************************************************** !
!
! InputReadNChars: reads and removes a specified number of characters from a 
!              string
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine InputReadNChars(input, option, chars, n, return_blank_error)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscTruth :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars

  if (input%ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,n
    chars(i:i) = ' '
  enddo

  input%ierr = len_trim(input%buf)
  if (input%ierr == 0) then
    if (return_blank_error) then
      input%ierr = 1
    else
      input%ierr = 0
    endif
    return
  else
    input%ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == achar(9)) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (input%buf(i:i) /= ' ' .and. input%buf(i:i) /= ',' .and. &
              input%buf(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo

    ends=i-1

    if (ends-begins+1 > n) then ! string read is too large for 'chars'
      ierr = 1
      return
    endif

    ! Copy (ends-begins) characters to 'chars'
    chars = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(ends+1:)

  endif

end subroutine InputReadNChars

! ************************************************************************** !
!
! InputFindStringInFile: Rewinds file and finds the first occurrence of
!                     'string'.  Note that the line must start with 'string'
!                     in order to match and that line is NOT returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
  subroutine InputFindStringInFile(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  
  character(len=strlen) :: string2, string3
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: length1, length2, i

  input%ierr = 0

 ! we initialize the word to blanks to avoid error reported by valgrind
  do i=1,MAXWORDLENGTH
     word(i:i) = ' '
  enddo
 
  length1 = len_trim(string)

  do 
    call fiReadFlotranString(fid,string2,ierr)
    string3 = string2
    call fiReadWord(string2,word,PETSC_TRUE,ierr)
    if (input%ierr /= 0) exit
    length2 = len_trim(word)
    if (length1 == length2 .and. fiStringCompare(string,word,length1)) exit
  enddo
  
  ! if not found, rewind once and try again.  this approach avoids excessive 
  ! reading if successive searches for strings are in descending order in 
  ! the file.
  if (input%ierr /= 0) then
    input%ierr = 0
    rewind(input%fid)
    do 
      call fiReadFlotranString(input%fid,string2,ierr)
      string3 = string2
      call fiReadWord(string2,word,PETSC_TRUE,ierr)
      if (ierr /= 0) exit
      length2 = len_trim(word)
      if (length1 == length2 .and. fiStringCompare(string,word,length1)) exit
    enddo
  endif    
  
  if (ierr == 0) string = trim(string3)

  end subroutine fiFindStringInFile
  
! ************************************************************************** !
!
! fiSkipToEND: Skips to keyword END
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine fiSkipToEND(fid,myrank,card)

  implicit none
  
  PetscInt :: fid
  character(len=MAXWORDLENGTH) :: card
  PetscMPIInt :: myrank
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  do
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(myrank,card,ierr)
    if (fiCheckExit(string)) exit
  enddo

end subroutine fiSkipToEND

! ************************************************************************** !
!
! fiCheckExit: Checks whether an end character (.,/,'END') has been found 
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function fiCheckExit(string)

  implicit none
  
  PetscTruth :: fiCheckExit
  character(len=*) :: string

  if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
      fiStringCompare(string,'END',THREE_INTEGER)) then
    fiCheckExit = PETSC_TRUE
  else
    fiCheckExit = PETSC_FALSE
  endif

end function fiCheckExit

! ************************************************************************** !
!
! InputDestroy: Deallocates an input object
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputDestroy(input)

  implicit none
  
  type(input_type), pointer :: input
  
  deallocate(input)
  nullify(input)
  
end subroutine InputDestroy
#endif

end module Input_module
