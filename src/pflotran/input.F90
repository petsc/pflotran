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

  interface InputReadFlotranString
!    module procedure InputReadFlotranString1
    module procedure InputReadFlotranString2
  end interface
  
  interface InputReadWord
    module procedure InputReadWord1
    module procedure InputReadWord2
  end interface
  
  interface InputReadNChars
    module procedure InputReadNChars1
    module procedure InputReadNChars2
  end interface
  
  interface InputReadInt
    module procedure InputReadInt1
    module procedure InputReadInt2
  end interface
  
  interface InputReadDouble
    module procedure InputReadDouble1
    module procedure InputReadDouble2
  end interface
  
  interface InputError
    module procedure InputError1
    module procedure InputError2
  end interface
  
  interface InputErrorMsg
    module procedure InputErrorMsg1
    module procedure InputErrorMsg2
  end interface
  
  interface InputDefaultMsg
    module procedure InputDefaultMsg1
    module procedure InputDefaultMsg2
  end interface
  
  interface InputReadStringErrorMsg
    module procedure InputReadStringErrorMsg1
    module procedure InputReadStringErrorMsg2
  end interface
  
  interface InputFindStringInFile
    module procedure InputFindStringInFile1
    module procedure InputFindStringInFile2
  end interface
  
  public :: InputCreate, InputDestroy, InputReadFlotranString, &
            InputReadWord, InputReadDouble, InputReadInt, InputCheckExit, &
            InputSkipToEND, InputFindStringInFile, InputErrorMsg, &
            InputDefaultMsg, InputReadStringErrorMsg, &
            InputFindStringErrorMsg, InputError, &
            InputReadNChars, InputReadQuotedWord, &
            InputGetCommandLineInt, &
            InputGetCommandLineReal, &
            InputGetCommandLineTruth, &
            InputGetCommandLineString

contains

! ************************************************************************** !
!
! InputCreate: Allocates and initializes a new Input object
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
function InputCreate(fid,filename)

  implicit none
  
  PetscInt :: fid
  character(len=*) :: filename
  
  type(input_type), pointer :: InputCreate
  PetscInt :: status  
  type(input_type), pointer :: input
  
  allocate(input)

  input%fid = fid
  input%ierr = 0
  input%buf = ''
  input%err_buf = ''
  input%err_buf2 = ''
  input%broadcast_read = PETSC_FALSE

  open(unit=input%fid,file=filename,status="old",iostat=status)
  if (status /= 0) then
    print *, 'file: ', trim(filename), ' not found'
    stop
  endif
  
  InputCreate => input
  
end function InputCreate

#if 1
! ************************************************************************** !
!
! InputDefaultMsg1: If ierr /= 0, informs user that default value will be used.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputDefaultMsg1(input,option,buffer)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer

  if (InputError(input)) then
    input%err_buf = buffer
    call InputDefaultMsg(input,option)
  endif

end subroutine InputDefaultMsg1

! ************************************************************************** !
!
! InputDefaultMsg2: If ierr /= 0, informs user that default value will be used.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputDefaultMsg2(input,option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer =  '"' // trim(input%err_buf) // &
                        '" set to default value.'
    call printMsg(option)
    input%ierr = 0
  endif

end subroutine InputDefaultMsg2

! ************************************************************************** !
!
! InputErrorMsg1: If ierr /= 0, If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputErrorMsg1(input,option,buffer1,buffer2)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer1, buffer2

  if (InputError(input)) then
    input%err_buf = buffer1
    input%err_buf2 = buffer2
    call InputErrorMsg(input,option)
  endif

end subroutine InputErrorMsg1

! ************************************************************************** !
!
! InputErrorMsg: If ierr /= 0, If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputErrorMsg2(input,option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer = 'While reading "' // trim(input%err_buf) // &
                       '" under keyword: ' // trim(input%err_buf2) // '.'
    call printErrMsg(option)
  endif

end subroutine InputErrorMsg2

! ************************************************************************** !
!
! InputReadStringErrorMsg1: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadStringErrorMsg1(input, option, buffer)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer

  if (InputError(input)) then
    input%err_buf = buffer
    call InputReadStringErrorMsg(input, option)
  endif

end subroutine InputReadStringErrorMsg1

! ************************************************************************** !
!
! InputReadStringErrorMsg2: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadStringErrorMsg2(input, option)

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer = 'While reading in string in "' // &
                       trim(input%err_buf) // '.'
    call printErrMsg(option)
  endif

end subroutine InputReadStringErrorMsg2

! ************************************************************************** !
!
! InputFindStringErrorMsg: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputFindStringErrorMsg(input, option, string)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: string

  if (InputError(input)) then
    option%io_buffer = 'Card (' // trim(string) // ') not &
                       &found in file.'
    call printErrMsg(option)    
  endif

end subroutine InputFindStringErrorMsg

! ************************************************************************** !
!
! InputReadInt1: reads and removes an integer value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadInt1(input, option, int)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: int

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (.not.InputError(input)) then
    read(word,*,iostat=input%ierr) int
  endif

end subroutine InputReadInt1

! ************************************************************************** !
!
! InputReadInt2: reads and removes an integer value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadInt2(string, option, int, ierr)

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscInt :: int
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  ierr = 0
  call InputReadWord(string,word,PETSC_TRUE,ierr)
  
  if (.not.InputError(ierr)) then
    read(word,*,iostat=ierr) int
  endif

end subroutine InputReadInt2

! ************************************************************************** !
!
! InputReadDouble1: reads and removes a real value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadDouble1(input, option, double)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscReal :: double

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (.not.InputError(input)) then
    read(word,*,iostat=input%ierr) double
  endif

end subroutine InputReadDouble1

! ************************************************************************** !
!
! InputReadDouble2: reads and removes a real value from a string
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadDouble2(string, option, double, ierr)

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscReal :: double
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  ierr = 0
  call InputReadWord(string,word,PETSC_TRUE,ierr)
  
  if (.not.InputError(ierr)) then
    read(word,*,iostat=ierr) double
  endif

end subroutine InputReadDouble2

#endif

#if 0
! ************************************************************************** !
!
! InputReadFlotranString: Reads a string (strlen characters long) from a 
!                           file while avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine InputReadFlotranString1(option, fid, string, ierr)

 
  implicit none

  type(option_type) :: option
  PetscInt :: fid
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr, ierr2

  if (option%broadcast_read) then
    if (option%myrank == option%io_rank) then
      call fiReadFlotranString(fid, string, ierr)
    endif
    call mpi_bcast(ierr,ONE_INTEGER,MPI_INTEGER,option%io_rank, &
                   option%mycomm,ierr2)
    if (ierr == 0) then  
      call mpi_bcast(string,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     option%io_rank,option%mycomm,ierr2)      
    endif
  else
    call fiReadFlotranString(fid, string, ierr)
  endif

end subroutine InputReadFlotranString1
#endif

#if 1
! ************************************************************************** !
!
! InputReadFlotranString2: Reads a string (strlen characters long) from a 
!                          file while avoiding commented or skipped lines.
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
      call InputReadFlotranStringSlave(input, option)
    endif
    call mpi_bcast(input%ierr,ONE_INTEGER,MPI_INTEGER,option%io_rank, &
                   option%mycomm,ierr2)
    if (.not.InputError(input)) then  
      call mpi_bcast(input%buf,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     option%io_rank,option%mycomm,ierr2)      
    endif
  else
    call InputReadFlotranStringSlave(input, option)
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

  use String_module
  
  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) ::  tempstring
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i

  input%ierr = 0

! we initialize the word to blanks to avoid error reported by valgrind
!  do i=1,MAXWORDLENGTH
!     word(i:i) = ' '
!  enddo
  word = ''
  
  do
    read(input%fid,'(a512)',iostat=input%ierr) input%buf

    if (InputError(input)) exit

    if (input%buf(1:1) == ':' .or. input%buf(1:1) == '!') cycle

    tempstring = input%buf
    call InputReadWord(tempstring,word,PETSC_TRUE,input%ierr)
    call StringToUpper(word)
    if (word(1:4) == 'SKIP') then
      do 
        read(input%fid,'(a512)',iostat=input%ierr) tempstring
       if (InputError(input) .and. OptionPrintToScreen(option)) then
          print *, 'End of file reached in InputReadFlotranStringSlave.'
          print *, 'SKIP encountered without matching NOSKIP.'
        endif
        call InputReadWord(tempstring,word,PETSC_FALSE,input%ierr)
        call StringToUpper(word)
        if (word(1:4) == 'NOSK') exit
      enddo
      if (InputError(input)) exit
    else if (word(1:1) /= ' ' .and. word(1:4) /= 'NOSK') then
      exit
    endif
  enddo
  
  ! Check for comment midway along a string
  if (.not.InputError(input)) then
    tempstring = input%buf
    do i=1,MAXSTRINGLENGTH
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

  if (InputError(input)) return

  ! Initialize character string to blank.
  do i=1,len_trim(word)
    word(i:i) = ' '
  enddo

  input%ierr = len_trim(input%buf)
  
  if (.not.InputError(input)) then
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
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

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

  ierr = len_trim(string)
  
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
    do while(string(i:i) == ' ' .or. string(i:i) == ',' .or. &
             string(i:i) == achar(9)) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine InputReadWord2

! ************************************************************************** !
!
! InputReadNChars1: reads and removes a specified number of characters from a 
!              string
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine InputReadNChars1(input, option, chars, n, return_blank_error)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscTruth :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars

  if (InputError(input)) return

  ! Initialize character string to blank.
  do i=1,n
    chars(i:i) = ' '
  enddo

  input%ierr = len_trim(input%buf)
  if (.not.InputError(input)) then
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
      input%ierr = 1
      return
    endif

    ! Copy (ends-begins) characters to 'chars'
    chars = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(ends+1:)

  endif

end subroutine InputReadNChars1

! ************************************************************************** !
!
! InputReadNChars2: reads and removes a specified number of characters from a 
!              string
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine InputReadNChars2(string, chars, n, return_blank_error, ierr)

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars
  PetscErrorCode :: ierr

  if (InputError(ierr)) return

  ! Initialize character string to blank.
  do i=1,n
    chars(i:i) = ' '
  enddo

  ierr = len_trim(string)
  if (.not.InputError(ierr)) then
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
    do while(string(i:i) == ' ' .or. string(i:i) == achar(9)) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo

    ends=i-1

    if (ends-begins+1 > n) then ! string read is too large for 'chars'
      ierr = 1
      return
    endif

    ! Copy (ends-begins) characters to 'chars'
    chars = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine InputReadNChars2

! ************************************************************************** !
!
! InputReadQuotedWord: reads and removes a word from a string, that is
!                      delimited by "'".
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine InputReadQuotedWord(input, option, word, return_blank_error)

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: i, begins, ends, realends
  PetscTruth :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: word
  PetscTruth :: openquotefound

  if (InputError(input)) return

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  do i=1,len_trim(word)
    word(i:i) = ' '
  enddo
  
!  if(len_trim(word) == 0) then
!    call InputReadFlotranString(input,option)
!    call InputReadStringErrorMsg(input,option,'trdatbse')
!  endif

  if (len_trim(input%buf) == 0) then
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

    if (input%buf(i:i) == "'") then
      openquotefound = PETSC_TRUE
      i=i+1
    endif

    begins=i

    if (openquotefound) then
      do while (input%buf(i:i) /= "'")
        if (i > (MAXWORDLENGTH-1)) exit
        i=i+1
      enddo
    else
    ! Count # of continuous characters (no blanks, commas, etc. in between)
      do while (input%buf(i:i) /= ' ' .and. input%buf(i:i) /= ',' .and. &
                input%buf(i:i) /= achar(9)) ! 9 = tab
        i=i+1
      enddo
    endif

    realends = i
    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'chars'
    word = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(realends+1:)
  endif

end subroutine InputReadQuotedWord

! ************************************************************************** !
!
! InputFindStringInFile1: Rewinds file and finds the first occurrence of
!                     'string'.  Note that the line must start with 'string'
!                     in order to match and that line is NOT returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
  subroutine InputFindStringInFile1(input, option, string)

  use String_module

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  call InputFindStringInFile2(input, option, string, PETSC_TRUE)
  
  end subroutine InputFindStringInFile1

! ************************************************************************** !
!
! InputFindStringInFile2: Rewinds file and finds the first occurrence of
!                     'string'.  Note that the line must start with 'string'
!                     in order to match and that line is NOT returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
  subroutine InputFindStringInFile2(input, option, string, print_warning)

  use String_module

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscTruth :: print_warning
  
  character(len=MAXWORDLENGTH) :: word
  PetscTruth :: found = PETSC_FALSE
  PetscInt :: length1, length2, i

  input%ierr = 0

  length1 = len_trim(string)

  do 
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) exit
    length2 = len_trim(word)
    if (length1 == length2 .and. StringCompare(string,word,length1)) then
      found = PETSC_TRUE
      exit
    endif
  enddo
  
  ! if not found, rewind once and try again.  this approach avoids excessive 
  ! reading if successive searches for strings are in descending order in 
  ! the file.
  if (InputError(input)) then
    input%ierr = 0
    rewind(input%fid)
    do 
      call InputReadFlotranString(input,option)
      if (InputError(input)) exit
      call InputReadWord(input,option,word,PETSC_TRUE)
      if (InputError(input)) exit
      length2 = len_trim(word)
      if (length1 == length2 .and. StringCompare(string,word,length1)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
  endif    
  
  if (.not.found .and. print_warning) then
    option%io_buffer = 'Card (' // trim(string) // ') not found in input file.'
    call printWrnMsg(option)
    input%ierr = 1
  endif
  
  end subroutine InputFindStringInFile2

! ************************************************************************** !
!
! InputSkipToEND: Skips to keyword END
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine InputSkipToEND(input,option,string)

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: string

  do
    call InputReadFlotranString(input,option)
    input%err_buf = 'End of file found before end of card ' // trim(string)
    call InputReadStringErrorMsg(input,option)
    if (InputCheckExit(input,option)) exit
  enddo

end subroutine InputSkipToEND

! ************************************************************************** !
!
! InputCheckExit: Checks whether an end character (.,/,'END') has been found 
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function InputCheckExit(input,option)

  use String_module
  
  implicit none

  type(input_type) :: input
  type(option_type) :: option  
  
  PetscTruth :: InputCheckExit

  if (input%buf(1:1) == '.' .or. input%buf(1:1) == '/' .or. &
      StringCompare(input%buf,'END',THREE_INTEGER)) then
    InputCheckExit = PETSC_TRUE
  else
    InputCheckExit = PETSC_FALSE
  endif

end function InputCheckExit

! ************************************************************************** !
!
! InputError1: Returns true if an error has occurred 
! author: Glenn Hammond
! date: 12/10/08
!
! ************************************************************************** !
function InputError1(input)

  implicit none

  type(input_type) :: input
  
  PetscTruth :: InputError1

  if (input%ierr == 0) then
    InputError1 = PETSC_FALSE
  else
    InputError1 = PETSC_TRUE
  endif

end function InputError1

! ************************************************************************** !
!
! InputError2: Returns true if an error has occurred 
! author: Glenn Hammond
! date: 12/10/08
!
! ************************************************************************** !
function InputError2(ierr)

  implicit none

  PetscErrorCode :: ierr
  
  PetscTruth :: InputError2

  if (ierr == 0) then
    InputError2 = PETSC_FALSE
  else
    InputError2 = PETSC_TRUE
  endif

end function InputError2

! ************************************************************************** !
!
! InputGetCommandLineInt: Returns integer value associated with a command 
!                          line argument
! author: Glenn Hammond
! date: 02/05/09
!
! ************************************************************************** !
subroutine InputGetCommandLineInt(string,int_value,found,option)

  use String_module
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscTruth :: found
  PetscInt :: int_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = iargc()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getarg(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getarg(iarg+1,string2)
        call InputReadInt(string2,option,int_value,ierr)
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'Integer argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineInt

! ************************************************************************** !
!
! InputGetCommandLineReal: Returns real*8 value associated with a command 
!                          line argument
! author: Glenn Hammond
! date: 02/05/09
!
! ************************************************************************** !
subroutine InputGetCommandLineReal(string,double_value,found,option)

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscTruth :: found
  PetscReal :: double_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = iargc()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getarg(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getarg(iarg+1,string2)
        call InputReadDouble(string2,option,double_value,ierr)
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'Real argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineReal


! ************************************************************************** !
!
! InputGetCommandLineString: Returns a string associated with a command 
!                          line argument
! author: Glenn Hammond
! date: 02/05/09
!
! ************************************************************************** !
subroutine InputGetCommandLineString(string,string_value,found,option)

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscTruth :: found
  character(len=MAXSTRINGLENGTH) :: string_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = iargc()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getarg(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getarg(iarg+1,string2)
        call InputReadNChars(string2,string_value,MAXSTRINGLENGTH, &
                             PETSC_TRUE,ierr)
        if (string_value(1:1) == '-') then
          ! no argument exists
          option%io_buffer = 'String argument (' // &
                             trim(adjustl(string_value)) // & 
                             ') for command line argument "' // &
                             trim(adjustl(string)) // '" not recognized.'
          call printErrMsg(option)
        endif
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'String argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineString

! ************************************************************************** !
!
! InputGetCommandLineTruth: Returns logical associated with a command 
!                           line argument
! author: Glenn Hammond
! date: 02/05/09
!
! ************************************************************************** !
subroutine InputGetCommandLineTruth(string,truth_value,found,option)

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscTruth :: found
  PetscTruth :: truth_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = iargc()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getarg(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getarg(iarg+1,string2)
        call InputReadWord(string2,word,PETSC_TRUE,ierr)
      else
        ! check if no argument exists, which is valid and means 'true'
        truth_value = PETSC_TRUE
        exit
      endif    
      if (word(1:1) == '-') then
        ! no argument exists, which is valid and means 'true'
        truth_value = PETSC_TRUE
        exit
      endif
      call StringToLower(word)
      select case(trim(word))
        case('yes','true','1','on')
          truth_value = PETSC_TRUE
        case('no','false','0','off')
          truth_value = PETSC_FALSE
        case default
          option%io_buffer = 'Truth argument for command line argument "' // &
                             trim(adjustl(string)) // '" not recognized.'
          call printErrMsg(option)
      end select
    endif
  enddo
  
end subroutine InputGetCommandLineTruth

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
  
  if (input%fid /= 0) close(input%fid)
  input%fid = 0
  deallocate(input)
  nullify(input)
  
end subroutine InputDestroy
#endif

end module Input_module
