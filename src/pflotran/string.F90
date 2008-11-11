module String_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "definitions.h"

  public :: StringCompare, &
            StringToUpper, &
            StringToLower, &
            StringReadQuotedWord

contains

! ************************************************************************** !
!
! StringCompare: compares two strings
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
PetscTruth function StringCompare(string1,string2,n)

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2
  
  do i=1,n
    if (string1(i:i) /= string2(i:i)) then
      StringCompare = PETSC_FALSE
      return
    endif
  enddo

  StringCompare = PETSC_TRUE
  return

end function StringCompare

! ************************************************************************** !
!
! StringToUpper: converts lowercase characters in a card to uppercase
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringToUpper(string)
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
      string(i:i) = achar(iachar(string(i:i)) - 32)
    endif
  enddo

end subroutine StringToUpper

! ************************************************************************** !
!
! StringToLower: converts uppercase characters in a card to lowercase
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringToLower(string)
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
      string(i:i) = achar(iachar(string(i:i)) + 32)
    endif
  enddo

end subroutine StringToLower

! ************************************************************************** !
!
! StringReadQuotedWord: reads and removes a name from a string read from the
!                       database.  "'" are used as delimiters.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringReadQuotedWord(string, name, return_blank_error, ierr)

  implicit none

  PetscInt :: i, begins, ends, realends, length
  PetscTruth :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: string
  character(len=*) :: name
  logical :: openquotefound
  PetscErrorCode :: ierr

  if (ierr /= 0) return

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  do i=1,len_trim(name)
    name(i:i) = ' '
  enddo

  ierr = 0
  length = len_trim(string)

  ! Remove leading blanks and tabs
  i=1
  do while(string(i:i) == ' ' .or. string(i:i) == achar(9)) 
    i=i+1
  enddo

  if (string(i:i) == "'") then
    openquotefound = PETSC_TRUE
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
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo
  endif

  realends = i
  ends=i-1

  ! Avoid copying beyond the end of the word (32 characters).
  if (ends-begins > MAXWORDLENGTH - 1) ends = begins + MAXWORDLENGTH - 1

  ! Copy (ends-begins) characters to 'chars'
  name = string(begins:ends)
  ! Remove chars from string
  string = string(realends+1:)

end subroutine StringReadQuotedWord

end module String_module
