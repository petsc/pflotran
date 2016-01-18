module Units_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"
  
  public :: UnitsConvertToInternal, UnitsConvertToExternal
  
contains

! ************************************************************************** !

function UnitsConvertToInternal(units,unit_category,option)
  ! 
  ! UnitsConvert: Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH), dimension(2) :: unit_category
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
          numerator = numerator * UnitsConvert(word,unit_category(1),option)
        else
          denominator = denominator * UnitsConvert(word,unit_category(2),option)
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

! jmf: Must do a search on this to correctly assign unit_category pass
function UnitsConvertToExternal(units,unit_category,option)
  ! 
  ! UnitsConvert: Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH), dimension(2) :: unit_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToExternal
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units,unit_category,option)

end function UnitsConvertToExternal

! ************************************************************************** !

function UnitsConvert(unit,unit_category,option)
  ! 
  ! Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  character(len=MAXWORDLENGTH) :: unit_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvert
    
  select case(trim(unit_category))
    ! convert volumes to cubic meters (m^3)
    case('volume')
      select case(trim(unit))
        case('cm^3')
          UnitsConvert = 1.d-6
        case('l','dm^3')
          UnitsConvert = 1.d-3
        case('m^3')
          UnitsConvert = 1.d0
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of VOLUME. Recognized&
                              units include: cm^3, l, dm^3, m^3 only.'
          call printErrMsg(option)
      end select
    ! convert areas to square meters (m^2)
    case('area') 
      select case(trim(unit))
        case('cm^2')
          UnitsConvert = 1.d-4
        case('dm^2')
          UnitsConvert = 1.d-2
        case('m^2')
          UnitsConvert = 1.d0 
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of AREA. Recognized&
                              units include: cm^2, dm^2, m^2 only.'
          call printErrMsg(option)
      end select
    ! convert lengths to meters (m)
    case('length')
      select case(trim(unit))
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
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of LENGTH. Recognized&
                              units include: km, m, dm, cm, mm only.'
          call printErrMsg(option)
      end select
    ! convert times to seconds (s)
    case('time')
      select case(trim(unit))
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
          UnitsConvert = 365.d0/12.d0*24.d0*3600.d0 
        case('y','yr','year')
          UnitsConvert = 365.d0*24.d0*3600.d0
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of TIME. Recognized&
                              units include: s, sec, second, min,&
                              minute, h, hr, hour, d, day, w, week,&
                              mo, month, y, yr, year only.'
          call printErrMsg(option)
      end select
    ! convert energy
    case('energy')
      select case(time(unit))
        case('J')   
          UnitsConvert = 1.d-6
        case('KJ')   
          UnitsConvert = 1.d-3
        case('MJ')   
          UnitsConvert = 1.d0
        case('W')   
          UnitsConvert = 1.d-6
        case('KW')   
          UnitsConvert = 1.d-3
        case('MW')   
          UnitsConvert = 1.d0
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of ENERGY. Recognized&
                              units include: J, KJ, MJ, W, KW, MW only.'
          call printErrMsg(option)
      end select
    ! convert mass to kilogram (kg)
    case('mass')
      select case(trim(unit))
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
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of MASS. Recognized&
                              units include: mol, mole, moles, ug, mg&
                              g, kg, 1, C, M only.'
          call printErrMsg(option)
      end select
    ! convert volumetric flow rate to cubic meters per second (m^3/s)
    case('flow_rate_volumetric')
      select case(trim(unit))
        case('gpm') !    l/gal     m^3/l   sec/min
          UnitsConvert = 3.785d0 * 1.d-3 / 60.d0
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of VOLUMETRIC FLOW RATE.& 
                              Recognized units include: gpm only.'
          call printErrMsg(option)
      end select
    ! convert mass flow rate to kilograms per second (kg/s)
    case('flow_rate_mass')
      select case(trim(unit))
        case default
          option%io_buffer = 'Unit "' // trim(unit) // '" is not a&
                              recognized unit of MASS FLOW RATE.& 
                              Recognized units include: ? only.'
          call printErrMsg(option)
      end select
    ! no unit category has been assigned during pass
    case('not_assigned')
      ! jmf: not sure what to do here yet!
      !
    ! cannot convert user-supplied units:
    case default
      option%io_buffer = 'Unit "' // trim(unit) // '" is not a recognized&
                          unit under any unit category.'
      call printErrMsg(option)
    end select
  end select


end function UnitsConvert

end module Units_module
