  module fileio_module

! use paramtrsmod

  implicit none
  integer, parameter :: dbaselen = 256, strlen = 256
  save

! by default, all variables are private
  private

! specify public variables
  public :: fiToUpper, fiToLower, fiWordToLower, fiWordToUpper,   &
            fiCharsToUpper, fiCharsToLower, &
            fiReadString, fiReadWord, fiReadCard, fiReadNChars,   &
            fiReadFlotranString, fiIsAlpha, &
            fiReadInt, fiReadDouble, fiReadMultDouble, &
            fiDefaultMsg, fiErrorMsg, fiReadStringErrorMsg

  public :: fiReadDBaseString, fiReadDBaseName, fiReadDBaseInt, &
            fiReadDBaseDouble, fiReadDBaseMultDouble

! __________________________________________________________________________ !
! __________________________________________________________________________ !

contains

! ************************************************************************** !
!
! fiDefaultMsg: If ierr /= 0, informs user that default value will be used.
! author: Glenn Hammond
! date: 11/06/00
!
! ************************************************************************** !
subroutine fiDefaultMsg(string, ierr)

  implicit none

  character(len=*) :: string
  integer :: i, ierr

  if (ierr /= 0) then
    print *, '"', (string(i:i),i=1,len_trim(string)), '" set to default value.'
    ierr = 0
  endif

end subroutine fiDefaultMsg

! ************************************************************************** !
!
! fiErrorMsg: If ierr /= 0, If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/06/00
!
! ************************************************************************** !
subroutine fiErrorMsg(string1, string2, ierr)

  implicit none

  character(len=*) :: string1, string2
  integer :: i, ierr

  if (ierr /= 0) then
    print *, 'Error reading "', (string1(i:i),i=1,len_trim(string1)), &
             '" under keyword: ',(string2(i:i),i=1,len_trim(string2)), '.'
    stop
  endif

end subroutine fiErrorMsg

! ************************************************************************** !
!
! fiReadStringErrorMsg: If ierr /= 0, informs user of error and stops.
! author: Glenn Hammond
! date: 11/06/00
!
! ************************************************************************** !
subroutine fiReadStringErrorMsg(string, ierr)

  implicit none

  character(len=*) :: string
  integer :: i, ierr

  if (ierr /= 0) then
    print *, 'Error reading in string in (', &
             (string(i:i),i=1,len_trim(string)), ').'
    stop
  endif

end subroutine fiReadStringErrorMsg

! ************************************************************************** !
!
! fiReadInt: reads and removes an integer value from a string
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
subroutine fiReadInt(string, int, ierr)

  implicit none

  character(len=strlen) :: string
  character(len=32) :: word
  integer :: int, ierr

  call fiReadWord(string,word,.true.,ierr)
  
  if (ierr == 0) then
    read(word,*,iostat=ierr) int
  endif

end subroutine fiReadInt

! ************************************************************************** !
!
! fiReadDouble: reads and removes a real value from a string
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
subroutine fiReadDouble(string, double, ierr)

  implicit none

  character(len=strlen) :: string
  character(len=32) :: word
  real*8 :: double
  integer :: ierr

  call fiReadWord(string,word,.true.,ierr)
  if (ierr == 0) then
    read(word,*,iostat=ierr) double
  endif

end subroutine fiReadDouble

! ************************************************************************** !
!
! fiReadMultDouble: reads n real numbers from an input file
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadMultDouble(fid, string, doubles, n, variable_name, keyword, &
                            subroutine_name)

  implicit none

  character(len=strlen) :: string
  character(len=*) :: variable_name, keyword, subroutine_name
  character(len=32) :: short_string
  integer :: i, n, ierr, fid
  logical :: newline
  real*8 :: doubles(n)

  i = 1
  ierr = 0
  newline = .false.
  do
    call fiReadDouble(string,doubles(i),ierr)
    if (ierr /= 0) then
      if (newline) then
        i = len_trim(variable_name)
        short_string(1:i) = variable_name(1:i)
        short_string(i+1:) = ' (too few values)'
        call fiErrorMsg(short_string,keyword,ierr)
      endif
      ierr = 0
      call fiReadFlotranString(fid,string,ierr)
      call fiReadStringErrorMsg(subroutine_name,ierr)
      newline = .true.
    else
      newline = .false.
      i = i + 1
      if (i > n) exit
    endif
  enddo

end subroutine fiReadMultDouble

! ************************************************************************** !
!
! fiReadString: Reads a string (strlen characters long) from a file
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
subroutine fiReadString(fid, string, ierr)

  implicit none

  character(len=strlen) :: string
  integer :: fid, ierr

!  if (ierr /= 0) return
  ierr = 0
  read(fid,'(a128)',iostat=ierr) string

end subroutine fiReadString

! ************************************************************************** !
!
! fiReadFlotranString: Reads a string (strlen characters long) from a file while
!                      avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine fiReadFlotranString(fid, string, ierr)

  implicit none

  character(len=strlen) :: string, tempstring, word
  integer :: i, fid, ierr

!  if (ierr /= 0) return
  ierr = 0
  
  do
    !read(fid,'(a128)',iostat=ierr) string
    read(fid,'(a256)',iostat=ierr) string
	!print *,string
	if (ierr /= 0) exit
    tempstring = string
    call fiReadWord(tempstring,word,.true.,ierr)
    call fiWordToUpper(word)
    if (word(1:4) == 'SKIP') then
      do 
        read(fid,'(a256)',iostat=ierr) tempstring
        if (ierr /= 0) then
          print *, 'End of file reached in fiReadFlotranString.'
          print *, 'SKIP encountered without matching NOSKIP.'
          stop
        endif
        call fiReadWord(tempstring,word,.false.,ierr)
        if (word(1:4) == 'NOSK') exit
      enddo
    else if (word(1:1) /= ':' .and. word(1:1) /= ' ' .and. &
             word(1:4) /= 'NOSK') then
      exit
    endif
  enddo

  ! Check for comment midway along a string
  if (ierr == 0) then
    tempstring = string
    do i=1,strlen
      string(i:i) = ' '
    enddo
    do i=1,len_trim(tempstring)
      if (tempstring(i:i) /= ':' .and. tempstring(i:i) /= '!') then
        string(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine fiReadFlotranString

! ************************************************************************** !
!
! fiReadCard: reads a four-letter card from a string.
! author: Glenn Hammond
! date: 6/27/00
!
! ************************************************************************** !
subroutine fiReadCard(word, card, ierr)

  implicit none

  integer :: i, ierr, length
  character(len=32) :: word
  character(len=4) :: card

  card = '    '
  if (ierr /= 0) then
    return
  else
    length = len_trim(word)
    if (length > 4) length = 4

    do i=1,length
      card(i:i) = word(i:i)
    enddo
  endif

end subroutine fiReadCard

! ************************************************************************** !
!
! fiReadWord: reads and removes a word (consecutive characters) from a string
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
subroutine fiReadWord(string, word, return_blank_error, ierr)

  implicit none

  integer :: i, ierr, begins, ends
  logical :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=strlen) :: string
  character(len=32) :: word

  if (ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,32
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

    ! Remove leading blanks
    
    i=1
    do while(string(i:i) == ' ') 
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
    if (ends-begins > 31) ends = begins + 31

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine fiReadWord

! ************************************************************************** !
!
! fiReadNChar: reads and removes a specified number of characters from a 
!              string
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine fiReadNChars(string, chars, n, return_blank_error, ierr)

  implicit none

  integer :: i, n, ierr, begins, ends
  logical :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=strlen) :: string
  character(len=n) :: chars

  if (ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,n
    chars(i:i) = ' '
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

    ! Remove leading blanks
    i=1
    do while(string(i:i) == ' ') 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',')
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

end subroutine fiReadNChars

! ************************************************************************** !
!
! fiWordToUpper: converts lowercase characters in a card to uppercase
! author: Glenn Hammond
! date: 6/27/00
!
! ************************************************************************** !
subroutine fiWordToUpper(word)
      
  implicit none

  integer :: i
  character(len=32) :: word

  do i=1,32
    word(i:i) = fiToUpper(word(i:i))
  enddo

end subroutine fiWordToUpper

! ************************************************************************** !
!
! fiWordToLower: converts uppercase characters in a card to lowercase
! author: Glenn Hammond
! date: 6/27/00
!
! ************************************************************************** !
subroutine fiWordToLower(word)
      
  implicit none

  integer :: i
  character(len=32) :: word

  do i=1,32
    word(i:i) = fiToLower(word(i:i))
  enddo

end subroutine fiWordToLower

! ************************************************************************** !
!
! fiCharsToUpper: converts lowercase characters in a card to uppercase
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine fiCharsToUpper(word, n)
      
  implicit none

  integer :: i, n
  character(len=n) :: word

  do i=1,n
    word(i:i) = fiToUpper(word(i:i))
  enddo

end subroutine fiCharsToUpper

! ************************************************************************** !
!
! fiCharsToLower: converts uppercase characters in a card to lowercase
! author: Glenn Hammond
! date: 11/02/00
!
! ************************************************************************** !
subroutine fiCharsToLower(word, n)
      
  implicit none

  integer :: i, n
  character(len=n) :: word

  do i=1,n
    word(i:i) = fiToLower(word(i:i))
  enddo

end subroutine fiCharsToLower

! ************************************************************************** !
!
! fiToUpper: converts a lowercase character to uppercase
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
character function fiToUpper(c)
      
  implicit none

  character :: c

  if (c >= 'a' .and. c <= 'z') then
    fiToUpper = achar(iachar(c) - 32)
  else
    fiToUpper = c
  endif

end function fiToUpper

! ************************************************************************** !
!
! fiToUpper: converts a uppercase character to lowercase
! author: Glenn Hammond
! date: 3/21/00
!
! ************************************************************************** !
character function fiToLower(c)
      
  implicit none

  character :: c

  if (c >= 'A' .and. c <= 'Z') then
    fiToLower = achar(iachar(c) + 32)
  else
    fiToLower = c
  endif

end function fiToLower

! ! ************************************************************************** !
!
! fiIsAlpha: returns true if 'c' is an Alphabet character
! author: Glenn Hammond
! date: 4/04/01
!
! ************************************************************************** !
logical function fiIsAlpha(c)
      
  implicit none

  character :: c

  if ((c >= 'A' .and. c <= 'Z') .or. &
      (c >= 'a' .and. c <= 'z')) then
    fiIsAlpha = .true.
  else
    fiIsAlpha = .false.
  endif

end function fiIsAlpha


! DATABASE SECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ************************************************************************** !
!
! fiReadDBaseString: Reads a long string (dbaselen characters long) from the
!                    database file while avoiding commented or skipped lines.
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadDBaseString(fid, string, ierr)

  implicit none

  character(len=dbaselen) :: string, tempstring, word
  integer :: i, fid, ierr

!  if (ierr /= 0) return
  ierr = 0
  
  do
    read(fid,'(a256)',iostat=ierr) string
    if (ierr /= 0) exit
    tempstring = string
    call fiReadDBaseWord(fid,tempstring,word,.true.,ierr)
    call fiWordToUpper(word)
    if (word(1:4) == 'SKIP') then
      do 
        read(fid,'(a256)',iostat=ierr) tempstring
        if (ierr /= 0) then
          print *, 'End of file reached in fiReadFlotranString.'
          print *, 'SKIP encountered without matching NOSKIP.'
          stop
        endif
        call fiReadDBaseWord(fid,tempstring,word,.false.,ierr)
        if (word(1:4) == 'NOSK') exit
      enddo
    else if (word(1:1) /= ':' .and. word(1:1) /= ' ' .and. &
             word(1:4) /= 'NOSK') then
      exit
    endif
  enddo

  ! Check for comment midway along a string
  if (ierr == 0) then
    tempstring = string
    do i=1,dbaselen
      string(i:i) = ' '
    enddo
    do i=1,len_trim(tempstring)
!      if (tempstring(i:i) /= ':' .and. tempstring(i:i) /= '!') then
     if (tempstring(i:i) /= '!') then
        string(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine fiReadDBaseString

! ************************************************************************** !
!
! fiReadDBaseName: reads and removes a name from a string read from the
!                  database.  "'" are used as delimiters.
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadDBaseName(fid, string, name, return_blank_error, ierr)

  implicit none

  integer :: i, ierr, begins, ends, realends, fid, length
  logical :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  integer, parameter :: n = 20
  character(len=dbaselen) :: string
  character(len=n) :: name
  logical :: openquotefound

  if (ierr /= 0) return

  openquotefound = .false.
  ! Initialize character string to blank.
  do i=1,n
    name(i:i) = ' '
  enddo

  
  if(len_trim(string) == 0) then
    call fiReadDBaseString(fid,string,ierr)
    call fiReadStringErrorMsg('trdatbse',ierr)
  endif

  ierr = 0
  length = len_trim(string)

  ! Remove leading blanks
  i=1
  do while(string(i:i) == ' ') 
    i=i+1
  enddo

  if (string(i:i) == "'") then
    openquotefound = .true.
    i=i+1
  endif

  begins=i

  if (openquotefound) then
    do while (string(i:i) /= "'")
      if (i > length) exit
      i=i+1
    enddo
  else
  ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',')
      i=i+1
    enddo
  endif

  realends = i
  ends=i-1

  ! Avoid copying beyond the end of the word (32 characters).
  if (ends-begins > n - 1) ends = begins + n - 1

  ! Copy (ends-begins) characters to 'chars'
  name = string(begins:ends)
  ! Remove chars from string
  string = string(realends+1:)

end subroutine fiReadDBaseName

! ************************************************************************** !
!
! fiReadDBaseWord: reads and removes a word (consecutive characters) from a 
!                  long string
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadDBaseWord(fid, string, word, return_blank_error, ierr)

  implicit none

  integer :: i, ierr, begins, ends, fid
  logical :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=dbaselen) :: string
  character(len=32) :: word

  if (ierr /= 0) return

  ! Initialize character string to blank.
  do i=1,32
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

    ! Remove leading blanks
    i=1
    do while(string(i:i) == ' ') 
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
    if (ends-begins > 31) ends = begins + 31

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine fiReadDBaseWord

! ************************************************************************** !
!
! fiReadDBaseInt: reads and removes an integer value from a long string
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadDBaseInt(fid, string, int, ierr)

  implicit none

  character(len=dbaselen) :: string
  character(len=32) :: word
  integer :: int, ierr, fid

  if(len_trim(string) == 0) then
    call fiReadDBaseString(fid,string,ierr)
    call fiReadStringErrorMsg('trdatbse',ierr)
  endif

  call fiReadDBaseWord(fid,string,word,.false.,ierr)
  if (ierr == 0) then
    read(word,*,iostat=ierr) int
  endif

end subroutine fiReadDBaseInt

! ************************************************************************** !
!
! fiReadDBaseDouble: reads and removes a real value from a long string
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
subroutine fiReadDBaseDouble(fid, string, double, ierr)

  implicit none

  character(len=dbaselen) :: string
  character(len=32) :: word
  real*8 :: double
  integer :: ierr, fid

  if(len_trim(string) == 0) then
    call fiReadDBaseString(fid,string,ierr)
    call fiReadStringErrorMsg('trdatbse',ierr)
  endif

  call fiReadDBaseWord(fid,string,word,.false.,ierr)
  if (ierr == 0) then
    read(word,*,iostat=ierr) double
  endif

end subroutine fiReadDBaseDouble

! ************************************************************************** !
!
! fiReadDBaseMultDouble: reads n real numbers from an input file
! author: Glenn Hammond
! date: 11/07/00
!
! ************************************************************************** !
  subroutine fiReadDBaseMultDouble(fid, string, doubles, n, ierr)

  implicit none

  character(len=dbaselen) :: string
  integer :: i, n, ierr, fid
  real*8 :: doubles(n)

  ierr = 0
  do i=1,n
    call fiReadDBaseDouble(fid,string,doubles(i),ierr)
    if (ierr /= 0) exit
  enddo

  end subroutine fiReadDBaseMultDouble

  end module fileio_module
