module Units_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"
  
  public :: UnitsConvertToInternal, UnitsConvertToExternal
  
contains

! ************************************************************************** !

function UnitsConvertToInternal(units,units_category,option)
  ! 
  ! UnitsConvert: Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! Notes: Updated/modified by Jenn Frederick 1/18/2016
  ! 

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXSTRINGLENGTH) :: numerator_buffer
  character(len=MAXSTRINGLENGTH) :: denominator_buffer
  character(len=*) :: units_category
  ! A maximum of 3 numerators/denominators are accepted 
  character(len=MAXWORDLENGTH) :: numerator_unit_cats(3)
  character(len=MAXWORDLENGTH) :: denominator_unit_cats(3)
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToInternal
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: char
  PetscReal :: numerator, denominator
  PetscInt :: istart, iend, k
  PetscInt :: length, num_length, denom_length
  PetscInt :: ind, num_denom_ind
  PetscBool :: set_numerator
  
  UnitsConvertToInternal = 1.d0
  numerator = 1.d0
  denominator = 1.d0

  numerator_unit_cats(:) = 'not_assigned'
  denominator_unit_cats(:) = 'not_assigned'

  length = len(units_category)
  ind = index(units_category,"/")
  
  if (ind == 0) then
  ! numerator(s) only
    numerator_buffer = units_category
    num_length = len_trim(numerator_buffer)
    denominator_buffer = ''
    denom_length = 0
  else
  ! denominator(s) exists        
    numerator_buffer = units_category(1:(ind-1))
    num_length = len_trim(numerator_buffer)
    denominator_buffer = units_category((ind+1):length)
    denom_length = len_trim(denominator_buffer)
  endif

  ! separate out numerator(s) unit categories
  ind = 1
  k = 1
  do while (ind /= 0)
    ind = index(numerator_buffer,"-")
    if (ind == 0) then
      numerator_unit_cats(k) = trim(numerator_buffer)
    else
      numerator_unit_cats(k) = numerator_buffer(1:(ind-1))
      numerator_buffer = numerator_buffer((ind+1):num_length)
      num_length = len_trim(numerator_buffer)
    endif
    k = k + 1
    if (k .gt. 3) then
      option%io_buffer = 'Maximum number of numerator units exceeded.&
                         & Please reduce to no more than 3 numerator units.' 
      call printErrMsg(option)
    endif
  enddo
  ! separate out denominator(s) unit categories
  if (denom_length /= 0) then
    ind = 1
    k = 1
    do while (ind /= 0)
      ind = index(denominator_buffer,"-")
      if (ind == 0) then
        denominator_unit_cats(k) = trim(denominator_buffer)
      else
        denominator_unit_cats(k) = denominator_buffer(1:(ind-1))
        denominator_buffer = denominator_buffer((ind+1):denom_length)
        denom_length = len_trim(denominator_buffer)
      endif
      k = k + 1
      if (k .gt. 3) then
        option%io_buffer = 'Maximum number of denominator units exceeded.&
                           & Please reduce to no more than 3 denominator units.' 
        call printErrMsg(option)
      endif
    enddo
  endif

  
  length = len_trim(units)
! if length == 0
  istart = 1
  iend = 1
  k = 1
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
          numerator = numerator * UnitsConvert(word,numerator_unit_cats(k),option)
          k = k + 1
        else
          denominator = denominator * UnitsConvert(word,denominator_unit_cats(k),option)
          k = k + 1
        endif
        if (char == '/') then
          set_numerator = PETSC_FALSE
          k = 1
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

function UnitsConvertToExternal(units,units_category,option)
  ! 
  ! UnitsConvert: Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  character(len=*) :: units_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToExternal
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units,units_category,option)

end function UnitsConvertToExternal

! ************************************************************************** !

function UnitsConvert(unit,units_category,option)
  ! 
  ! Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  character(len=MAXSTRINGLENGTH) :: units_string
  character(len=*) :: units_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvert
    
  select case(trim(units_category))
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
          units_string = 'cm^3, l, dm^3, m^3'
          call UnitsError(unit,units_category,units_string,option)
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
          units_string = 'cm^2, dm^2, m^2'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert lengths to meters (m)
    case('length')
      select case(trim(unit))
        case('km')
          UnitsConvert = 1000.d0
        case('m','met','meter')
          UnitsConvert = 1.d0
        case('dm')
          UnitsConvert = 1.d-1
        case('cm')
          UnitsConvert = 1.d-2
        case('mm')
          UnitsConvert = 1.d-3
        case default
          units_string = 'km, m, dm, cm, mm'
          call UnitsError(unit,units_category,units_string,option)
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
          units_string = 's, sec, second, min, minute, h, hr, hour,&
                         & d, day, w, week, mo, month, y, yr, year'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert energy
    case('energy')
      select case(trim(unit))
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
          units_string = 'J, KJ, MJ, W, KW, MW'
          call UnitsError(unit,units_category,units_string,option)
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
        case default
          units_string = 'mol, mole, moles, ug, mg, g, kg'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert temperature to Celsius (C)
    case('temperature')
      select case(trim(unit))
        case('C') 
          UnitsConvert = 1.d0
        case default
          units_string = 'C (celsius)'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert pressure to Pascal (Pa)
    case('pressure')
      select case(trim(unit))
        case('Pa') 
          UnitsConvert = 1.d0
        case('KPa')
          UnitsConvert = 1.d3
        case('MPa')
          UnitsConvert = 1.d6  
        case default
          units_string = 'Pa, KPa, MPa'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert concentrations to (?)
    ! jmf: check this! (two defaults possible, but they will always be consistent within the input file)
    case('concentration')
      select case(trim(unit))
        case('M') 
          UnitsConvert = 1.d0
        case('m')
          UnitsConvert = 1.d0  
        case default
          units_string = 'M (molarity), m (molality)'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! convert saturations to (?)
    ! jfm: check this!
    case('saturation')
    ! convert forces to Newtons (N)
    case('force')
      select case(trim(unit))
        case('N') 
          UnitsConvert = 1.d0
        case default
          units_string = 'N (Newton)'
          call UnitsError(unit,units_category,units_string,option)
      end select
    ! the parameter is unitless (1)
    case('unitless')
      UnitsConvert = 1.d0
    ! no unit category has been assigned during pass
    case('not_assigned')
      option%io_buffer = 'Unit category has not been assigned.'
      call printErrMsg(option)
    case('unknown')
      option%io_buffer = 'Unit category set to unknown.'
      call printErrMsg(option)
    ! cannot convert user-supplied units:
    case default
      option%io_buffer = 'Unit "' // trim(unit) // '" is not a recognized&
                          & unit under any unit category.'
      call printErrMsg(option)
  end select

end function UnitsConvert

! ************************************************************************** !

subroutine UnitsError(unit,units_category,units_string,option)
  ! 
  ! Forms and prints a units conversion error message.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/18/2016
  ! 

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  character(len=*) :: units_category
  character(len=MAXSTRINGLENGTH) :: units_string
  type(option_type) :: option

  option%io_buffer = 'Unit "' // trim(unit) // '" is not a recognized unit&
                     & of "' // units_category // '". Recognized units&
                     & include: "' // trim(units_string) // '" only.'
  call printErrMsg(option)
  
end subroutine UnitsError

end module Units_module
