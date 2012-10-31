module Units_module

  implicit none
  
  private
  
#include "definitions.h"
  
  public :: UnitsConvertToInternal, UnitsConvertToExternal
  
contains

! ************************************************************************** !
!
! UnitsConvert: Converts units to pflotran internal units
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
function UnitsConvertToInternal(units,option)

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToInternal
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: char
  PetscReal :: numerator, denominator
  PetscInt :: istart, iend
  PetscInt :: length
  PetscBool :: set_numerator
  
  UnitsConvertToInternal = 1.d0
  
  numerator = 1.d0
  denominator = 1.d0
  
  length = len_trim(units)
  istart = 1
  iend = 1
  set_numerator = PETSC_TRUE
  do 
    if (istart > length) exit
    do
      char = units(iend:iend)
      if (char == '-' .or. char == '/' .or. char == ' ' .or. iend == length) then
        if (iend == length) then
          word = units(istart:max(istart,iend))
        else
          word = units(istart:max(istart,iend-1))
        endif
        word = adjustl(word)
        if (set_numerator) then
          numerator = numerator * UnitsConvert(word,option)
        else
          denominator = denominator * UnitsConvert(word,option)
        endif
        if (char == '/') then
          set_numerator = PETSC_FALSE
        endif
        istart = iend + 1
        iend = istart
        exit
      endif
      iend = iend + 1
    enddo
  enddo
  
  UnitsConvertToInternal = numerator / denominator

end function UnitsConvertToInternal

! ************************************************************************** !
!
! UnitsConvert: Converts units to pflotran internal units
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
function UnitsConvertToExternal(units,option)

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToExternal
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units,option)

end function UnitsConvertToExternal

! ************************************************************************** !
!
! UnitsConvert: Converts units to pflotran internal units
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
function UnitsConvert(unit,option)

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  type(option_type) :: option
  
  PetscReal :: UnitsConvert
    
  select case(trim(unit))
    ! convert volumes to m^3
    case('cm^3')
      UnitsConvert = 1.d-6
    case('l','dm^3')
      UnitsConvert = 1.d-3
    case('m^3')
      UnitsConvert = 1.d0
    ! convert areas to m^2
    case('cm^2')
      UnitsConvert = 1.d-4
    case('dm^2')
      UnitsConvert = 1.d-2
    case('m^2')
      UnitsConvert = 1.d0
    ! convert times to seconds
    case('s','sec','second')
      UnitsConvert = 1.d0
    case('min','minute')
      UnitsConvert = 60.d0
    case('h','hr','hour')
      UnitsConvert = 3600.d0
    case('d','day')
      UnitsConvert = 24.d0*3600.d0 
    case('w','week')
      UnitsConvert = 7.d0*24.d0*3600.d0 
    case('mo','month')
      UnitsConvert = 30.d0*24.d0*3600.d0 
    case('y','yr')
      UnitsConvert = 365.d0*24.d0*3600.d0
    ! convert length to meters
    case('km')
      UnitsConvert = 1000.d0
    case('m')
      UnitsConvert = 1.d0
    case('dm')
      UnitsConvert = 1.d-1
    case('cm')
      UnitsConvert = 1.d-2
    case('mm')
      UnitsConvert = 1.d-3
    ! convert mass to kg
    case('mol','mole','moles')
      UnitsConvert = 1.d0
    case('ug')
      UnitsConvert = 1.d-9
    case('mg')
      UnitsConvert = 1.d-6
    case('g')
      UnitsConvert = 1.d-3
    case('kg')
      UnitsConvert = 1.d0
    case('1') ! one
      UnitsConvert = 1.d0
    case('C') ! one
      UnitsConvert = 1.d0
    case('M') ! one
      UnitsConvert = 1.d0
    case('gpm') ! convert from gallons per minute to m^3/sec
                !    l/gal     m^3/l   sec/min
      UnitsConvert = 3.785d0 * 1.d-3 / 60.d0
    case default
       option%io_buffer = 'Unit "' // trim(unit) // '" not recognized.'
       call printErrMsg(option)
  end select

end function UnitsConvert

end module Units_module
