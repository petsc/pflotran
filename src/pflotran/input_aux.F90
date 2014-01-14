module Input_Aux_module

  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: input_type 
    PetscInt :: fid
    PetscErrorCode :: ierr
    character(len=MAXWORDLENGTH) :: filename
    character(len=MAXSTRINGLENGTH) :: buf
    character(len=MAXSTRINGLENGTH) :: err_buf
    character(len=MAXSTRINGLENGTH) :: err_buf2
    PetscBool :: broadcast_read
  end type input_type

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
#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64)
    ! If PetscInt and PetscMPIInt have different sizes (occurs for some builds
    ! with 64 bit indices), then we need to have additional routines for the 
    ! InputReadInt() generic subroutine.  (We use the above check instead of 
    ! directly checking to see if PetscInt and PetscMPIInt have the same size
    ! because the size of PetscInt is not included in the 
    ! $PETSC_DIR/$PETSC_ARCH/include/petscconf.h file.) If the two types have
    ! the same size, then these additional routines for type PetscMPIInt must
    ! *not* be defined, because then the interface becomes ambiguous, since 
    ! Fortran doesn't know the difference between PetscInt and PetscMPIInt if
    ! they are identically sized integers.  --RTM
    module procedure InputReadInt3
    module procedure InputReadInt4
#endif
  end interface
  
  interface InputReadDouble
    module procedure InputReadDouble1
    module procedure InputReadDouble2
  end interface
  
  interface InputReadNDoubles
    module procedure InputReadNDoubles1
    module procedure InputReadNDoubles2
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
  
  public :: InputCreate, InputDestroy, InputReadPflotranString, &
            InputReadWord, InputReadDouble, InputReadInt, InputCheckExit, &
            InputReadNDoubles, &
            InputSkipToEND, InputFindStringInFile, InputErrorMsg, &
            InputDefaultMsg, InputReadStringErrorMsg, &
            InputFindStringErrorMsg, InputError, &
            InputReadNChars, InputReadQuotedWord, &
            InputReadPath, &
            InputGetCommandLineInt, &
            InputGetCommandLineReal, &
            InputGetCommandLineTruth, &
            InputGetCommandLineString, &
            InputReadFilenames

contains

! ************************************************************************** !

function InputCreate(fid,filename,option)
  ! 
  ! Allocates and initializes a new Input object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  use Option_module

  implicit none
  
  PetscInt :: fid
  character(len=*) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: InputCreate
  PetscInt :: status  
  type(input_type), pointer :: input
  
  allocate(input)

  input%fid = fid
  input%filename = filename
  input%ierr = 0
  input%buf = ''
  input%err_buf = ''
  input%err_buf2 = ''
  input%broadcast_read = PETSC_FALSE
  
  if (fid == MAX_IN_UNIT) then
    option%io_buffer = 'MAX_IN_UNIT in definitions.h must be increased to' // &
      ' accommodate a larger number of embedded files.'
    call printErrMsg(option)
  endif

  open(unit=input%fid,file=filename,status="old",iostat=status)
  if (status /= 0) then
    if (len_trim(filename) == 0) filename = '<blank>'
    option%io_buffer = 'File: "' // trim(filename) // '" not found.'
    call printErrMsg(option)
  endif
  
  InputCreate => input
  
end function InputCreate

! ************************************************************************** !

subroutine InputDefaultMsg1(input,option,buffer)
  ! 
  ! If ierr /= 0, informs user that default value will be used.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputDefaultMsg2(input,option)
  ! 
  ! If ierr /= 0, informs user that default value will be used.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputErrorMsg1(input,option,buffer1,buffer2)
  ! 
  ! If ierr /= 0, If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputErrorMsg2(input,option)
  ! 
  ! InputErrorMsg: If ierr /= 0, If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputReadStringErrorMsg1(input, option, buffer)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputReadStringErrorMsg2(input, option)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer = 'While reading in string in "' // &
                       trim(input%err_buf) // '".'
    call printErrMsg(option)
  endif

end subroutine InputReadStringErrorMsg2

! ************************************************************************** !

subroutine InputFindStringErrorMsg(input, option, string)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputReadInt1(input, option, int)
  ! 
  ! reads and removes an integer value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputReadInt2(string, option, int, ierr)
  ! 
  ! reads and removes an integer value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64)

! ************************************************************************** !

subroutine InputReadInt3(input, option, int)
  ! 
  ! InputReadInt3() and InputReadInt4() must only be defined if PetscInt and
  ! PetscMPIInt differ in size.  See notes above in the interface definition.
  ! --RTM
  ! reads and removes an integer value from a string
  ! authors: Glenn Hammond, Richard Mills
  ! 
  ! Date: 2/3/2012
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscMPIInt :: int

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (.not.InputError(input)) then
    read(word,*,iostat=input%ierr) int
  endif

end subroutine InputReadInt3

! ************************************************************************** !

subroutine InputReadInt4(string, option, int, ierr)
  ! 
  ! reads and removes an integer value from a string
  ! authors: Glenn Hammond, Richard Mills
  ! 
  ! Date: 2/3/2012
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscMPIInt :: int
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  ierr = 0
  call InputReadWord(string,word,PETSC_TRUE,ierr)
  
  if (.not.InputError(ierr)) then
    read(word,*,iostat=ierr) int
  endif

end subroutine InputReadInt4

#endif

! ************************************************************************** !

subroutine InputReadDouble1(input, option, double)
  ! 
  ! End of defined(PETSC_USE_64BIT_INDICES) &&
  ! (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64) conditional
  ! reads and removes a real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

subroutine InputReadDouble2(string, option, double, ierr)
  ! 
  ! reads and removes a real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

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

! ************************************************************************** !

subroutine InputReadNDoubles1(input, option, double, n)
  ! 
  ! reads and removes "n" real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/11
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: n
  PetscReal :: double(n)

  PetscInt :: i

  do i = 1, n
    call InputReadDouble(input,option,double(i))
    if (InputError(input)) return
  enddo

end subroutine InputReadNDoubles1

! ************************************************************************** !

subroutine InputReadNDoubles2(string, option, double, n, ierr)
  ! 
  ! reads and removes "n" real values from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/11
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscInt :: n
  PetscReal :: double(n)
  PetscErrorCode :: ierr

  PetscInt :: i

  do i = 1, n
    call InputReadDouble(string,option,double(i),ierr)
    if (InputError(ierr)) return
  enddo

end subroutine InputReadNDoubles2

! ************************************************************************** !

subroutine InputReadPflotranString(input, option)
  ! 
  ! Reads a string (strlen characters long) from a
  ! file while avoiding commented or skipped lines.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  PetscInt :: flag

  if (input%broadcast_read) then
    if (option%myrank == option%io_rank) then
      call InputReadPflotranStringSlave(input, option)
    endif
    flag = input%ierr
    call MPI_Bcast(flag,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                   option%mycomm,ierr)
    input%ierr = flag
    if (.not.InputError(input)) then  
      call MPI_Bcast(input%buf,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     option%io_rank,option%mycomm,ierr)      
    endif
  else
    call InputReadPflotranStringSlave(input, option)
  endif

end subroutine InputReadPflotranString

! ************************************************************************** !

subroutine InputReadPflotranStringSlave(input, option)
  ! 
  ! Reads a string (strlen characters long) from a
  ! file while avoiding commented or skipped lines.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  use String_module
  
  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) ::  tempstring
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  PetscInt :: skip_count

  input%ierr = 0

! we initialize the word to blanks to avoid error reported by valgrind
!  do i=1,MAXWORDLENGTH
!     word(i:i) = ' '
!  enddo
  word = ''
  
  do
    read(input%fid,'(a512)',iostat=input%ierr) input%buf
    call StringAdjustl(input%buf)

    if (InputError(input)) exit

    if (input%buf(1:1) == ':' .or. input%buf(1:1) == '!' .or. &
        input%buf(1:1) == '#' ) cycle

    tempstring = input%buf
    call InputReadWord(tempstring,word,PETSC_TRUE,input%ierr)
    call StringToUpper(word)
    if (word(1:4) == 'SKIP') then
      skip_count = 1
      do 
        read(input%fid,'(a512)',iostat=input%ierr) tempstring
        if (InputError(input)) then
          option%io_buffer = 'End of file reached in ' // &
              'InputReadPflotranStringSlave.  SKIP encountered ' // &
              'without a matching NOSKIP.'
          call printErrMsg(option)              
        endif
        call InputReadWord(tempstring,word,PETSC_FALSE,input%ierr)
        call StringToUpper(word)
        if (word(1:4) == 'SKIP') skip_count = skip_count + 1
        if (word(1:4) == 'NOSK') then
          skip_count = skip_count - 1
          if (skip_count == 0) exit
        endif
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
      if (tempstring(i:i) /= ':' .and. tempstring(i:i) /= '!' .and. &
          tempstring(i:i) /= '#' ) then
        input%buf(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine InputReadPflotranStringSlave

! ************************************************************************** !

subroutine InputReadWord1(input, option, word, return_blank_error)
  ! 
  ! reads and removes a word (consecutive characters) from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: return_blank_error
  
  PetscInt :: i, begins, ends, lenword
  character(len=1) :: tab, backslash

  if (InputError(input)) return
  
  call InputReadWord2(input%buf, word, return_blank_error, input%ierr)

end subroutine InputReadWord1

! ************************************************************************** !

subroutine InputReadWord2(string, word, return_blank_error, ierr)
  ! 
  ! reads and removes a word (consecutive characters) from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscBool :: return_blank_error
  PetscErrorCode :: ierr
  
  PetscInt :: i, begins, ends
  character(len=1) :: tab, backslash  

  if (ierr /= 0) return

  tab = achar(9)
  backslash = achar(92)
  
  ! Initialize character string to blank.
  ! Initialize character string to blank.  len_trim(word) is not
  ! defined if word is allocated but not initialized.  This works on
  ! most compilers, but may not work on some?  Holler if it
  ! errors... - etc
  word = ''
  ! do i=1,len_trim(word)
  !   word(i:i) = ' '
  ! enddo

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
             string(i:i) == tab) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab .and. &
              (i == begins .or. string(i:i) /= backslash))
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

subroutine InputReadNChars1(input, option, chars, n, return_blank_error)
  ! 
  ! reads and removes a specified number of characters from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/00
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscBool :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars
  character(len=1) :: tab, backslash    

  if (InputError(input)) return

  call InputReadNChars2(input%buf, chars, n, return_blank_error, input%ierr)
  
end subroutine InputReadNChars1

! ************************************************************************** !

subroutine InputReadNChars2(string, chars, n, return_blank_error, ierr)
  ! 
  ! reads and removes a specified number of characters from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/00
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars
  PetscErrorCode :: ierr
  character(len=1) :: tab, backslash    

  if (InputError(ierr)) return

  tab = achar(9)
  backslash = achar(92)

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
    do while(string(i:i) == ' ' .or. string(i:i) == tab) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab  .and. &
              (i == begins .or. string(i:i) /= backslash))
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

subroutine InputReadQuotedWord(input, option, word, return_blank_error)
  ! 
  ! reads and removes a word from a string, that is
  ! delimited by "'".
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/00
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: i, begins, ends, realends
  PetscBool :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: word
  PetscBool :: openquotefound
  character(len=1) :: tab, backslash    

  if (InputError(input)) return

  tab = achar(9)
  backslash = achar(92)

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  do i=1,len_trim(word)
    word(i:i) = ' '
  enddo
  
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
    do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == tab) 
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
                input%buf(i:i) /= tab .and. &
                (i == begins .or. input%buf(i:i) /= backslash))
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

subroutine InputReadPath(string, word, return_blank_error, ierr)
  ! 
  ! reads and removes a words from a path
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/10
  ! 

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscBool :: return_blank_error
  PetscErrorCode :: ierr
  
  PetscInt :: i, begins, ends
  character(len=1) :: slash, backslash  

  if (ierr /= 0) return

  slash = achar(47)
  backslash = achar(92)

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
    do while(string(i:i) == ' ' .and. string(i:i) == slash) 
      i=i+1
    enddo

    begins=i

    ! Count # of characters (no slashes in between)
    do while (string(i:i) /= slash .and. &
              (i == begins .or. string(i:i) /= backslash))
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
  
end subroutine InputReadPath

! ************************************************************************** !

subroutine InputFindStringInFile1(input, option, string)
  ! 
  ! Rewinds file and finds the first occurrence of
  ! 'string'.  Note that the line must start with 'string'
  ! in order to match and that line is NOT returned
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  use String_module

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  call InputFindStringInFile2(input, option, string, PETSC_TRUE)
  
end subroutine InputFindStringInFile1

! ************************************************************************** !

subroutine InputFindStringInFile2(input, option, string, print_warning)
  ! 
  ! Rewinds file and finds the first occurrence of
  ! 'string'.  Note that the line must start with 'string'
  ! in order to match and that line is NOT returned
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  use String_module

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: print_warning
  
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found = PETSC_FALSE
  PetscInt :: length1, length2, i

  input%ierr = 0

  length1 = len_trim(string)

  do 
    call InputReadPflotranString(input,option)
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
      call InputReadPflotranString(input,option)
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

subroutine InputSkipToEND(input,option,string)
  ! 
  ! Skips to keyword END
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: string

  do
    call InputReadPflotranString(input,option)
    input%err_buf = 'End of file found before end of card ' // trim(string)
    call InputReadStringErrorMsg(input,option)
    if (InputCheckExit(input,option)) exit
  enddo

end subroutine InputSkipToEND

! ************************************************************************** !

function InputCheckExit(input,option)
  ! 
  ! Checks whether an end character (.,/,'END') has been found
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use String_module
  
  implicit none

  type(input_type) :: input
  type(option_type) :: option  
  PetscInt :: i
  character(len=1) :: tab
  
  PetscBool :: InputCheckExit

  ! We must remove leading blanks and tabs. --RTM
  tab = achar(9)
  i=1
  do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == tab) 
    i=i+1
  enddo

  if (input%buf(i:i) == '/' .or. &
      StringCompare(input%buf(i:),'END',THREE_INTEGER)) then
    InputCheckExit = PETSC_TRUE
  else
    InputCheckExit = PETSC_FALSE
  endif

end function InputCheckExit

! ************************************************************************** !

function InputError1(input)
  ! 
  ! Returns true if an error has occurred
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/08
  ! 

  implicit none

  type(input_type) :: input
  
  PetscBool :: InputError1

  if (input%ierr == 0) then
    InputError1 = PETSC_FALSE
  else
    InputError1 = PETSC_TRUE
  endif

end function InputError1

! ************************************************************************** !

function InputError2(ierr)
  ! 
  ! Returns true if an error has occurred
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/08
  ! 

  implicit none

  PetscErrorCode :: ierr
  
  PetscBool :: InputError2

  if (ierr == 0) then
    InputError2 = PETSC_FALSE
  else
    InputError2 = PETSC_TRUE
  endif

end function InputError2

! ************************************************************************** !

subroutine InputGetCommandLineInt(string,int_value,found,option)
  ! 
  ! Returns integer value associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscInt :: int_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
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

subroutine InputGetCommandLineReal(string,double_value,found,option)
  ! 
  ! Returns real*8 value associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscReal :: double_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
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

subroutine InputGetCommandLineString(string,string_value,found,option)
  ! 
  ! Returns a string associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: string_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
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

subroutine InputGetCommandLineTruth(string,truth_value,found,option)
  ! 
  ! Returns logical associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscBool :: truth_value

  PetscInt :: iarg, narg
  PetscInt :: len
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  len = len_trim(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2,len)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
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

function getCommandLineArgumentCount()
  ! 
  ! Returns the number of command line arguments
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/10
  ! 

  implicit none
  
  integer :: iargc
  
  PetscInt :: getCommandLineArgumentCount
  
  ! initialize to zero
  getCommandLineArgumentCount = 0
  
#if defined(PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT)
  getCommandLineArgumentCount = command_argument_count()
#elif defined(PETSC_HAVE_GETARG)
  getCommandLineArgumentCount = iargc()
#endif

end function getCommandLineArgumentCount

! ************************************************************************** !

subroutine getCommandLineArgument(i,arg)
  ! 
  ! Returns the ith command line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/10
  ! 

  implicit none
  
  PetscInt :: i
  character(len=*) :: arg

  integer*4 :: fortran_int

  fortran_int = i
#if defined(PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT)
  call get_command_argument(fortran_int,arg)
#elif defined(PETSC_HAVE_GETARG)
  call getarg(fortran_int,arg)
#endif

end subroutine getCommandLineArgument

! ************************************************************************** !

subroutine InputReadFilenames(option,filenames)
  ! 
  ! Reads filenames for multi-simulation runs
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 

  use Option_module

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: filename_count
  type(input_type), pointer :: input
  PetscBool :: card_found

  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = "FILENAMES"
  call InputFindStringInFile(input,option,string) 

  card_found = PETSC_FALSE
  if (InputError(input)) then
    ! if the FILENAMES card is not included, we will assume that only
    ! filenames exist in the file.
    rewind(input%fid)
  else
    card_found = PETSC_TRUE
  endif
    
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
  enddo
  
  allocate(filenames(filename_count))
  filenames = ''
  rewind(input%fid) 

  if (card_found) then
    string = "FILENAMES"
    call InputFindStringInFile(input,option,string) 
  endif
  
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
    filenames(filename_count) = filename
  enddo

  call InputDestroy(input)

end subroutine InputReadFilenames

! ************************************************************************** !

subroutine InputDestroy(input)
  ! 
  ! Deallocates an input object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none
  
  type(input_type), pointer :: input
  
  if (input%fid /= 0) close(input%fid)
  input%fid = 0
  deallocate(input)
  nullify(input)
  
end subroutine InputDestroy

end module Input_Aux_module
