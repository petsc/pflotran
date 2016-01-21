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
  ! Notes: Updated/modified by Jenn Frederick 1/21/2016
  ! 

  use Option_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH) :: numerator_units
  character(len=MAXWORDLENGTH) :: denominator_units
  character(len=*) :: units_category
  ! A maximum of 3 numerators/denominators are accepted 
  character(len=MAXWORDLENGTH) :: numerator_unit_cats
  character(len=MAXWORDLENGTH) :: denominator_unit_cats
  type(option_type) :: option
  PetscReal :: numerator, denominator
  PetscBool :: set_units_denom, set_cats_denom
  PetscInt :: length
  PetscInt :: ind
  PetscReal :: UnitsConvertToInternal
  
  UnitsConvertToInternal = 1.d0
  set_units_denom = PETSC_FALSE
  set_cats_denom = PETSC_FALSE
  numerator = 1.d0
  denominator = 1.d0
  numerator_units = ''
  denominator_units = ''
  numerator_unit_cats = 'not_assigned'
  denominator_unit_cats = 'not_assigned'

  length = len(units_category)
  ind = index(units_category,"/")   
  if (ind == 0) then
  ! numerator(s) only
    numerator_unit_cats = trim(units_category)
  else
  ! denominator(s) also exists
    set_cats_denom = PETSC_TRUE
    numerator_unit_cats = units_category(1:(ind-1))
    denominator_unit_cats = units_category((ind+1):length)
  endif

  length = len_trim(units)
  ind = index(units,"/")
  if (ind == 0) then
  ! numerator(s) only
    numerator_units = trim(units)
  else
  ! denominator(s) also exists
    set_units_denom = PETSC_TRUE
    numerator_units = units(1:(ind-1))
    denominator_units = units((ind+1):length)
  endif

  ! check if provided units structure matches units category structure
  if (set_units_denom .neqv. set_cats_denom) then
    if ((set_cats_denom .eqv. PETSC_TRUE) .and. &
        (set_units_denom .eqv. PETSC_FALSE)) then
      option%io_buffer = 'Units provided to not match the expected &
                         &unit category structure (numerator/denominator). &
                         &A unit denominator was expected, but not given.'
      call printErrMsg(option)
    endif
    if ((set_cats_denom .eqv. PETSC_FALSE) .and. &
        (set_units_denom .eqv. PETSC_TRUE)) then
      option%io_buffer = 'Units provided to not match the expected &
                         &unit category structure (numerator/denominator). &
                         &A unit denominator was not expected, but was given.'
      call printErrMsg(option)
    endif
  endif

  numerator = numerator * UnitsConvert(numerator_units, &
                                       numerator_unit_cats, &
                                       option)
  if (set_cats_denom) then
    denominator = denominator * UnitsConvert(denominator_units, &
                                             denominator_unit_cats, &
                                             option)
  endif

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
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units, &
                                                       units_category,option)

end function UnitsConvertToExternal

! ************************************************************************** !

function UnitsConvert(units,units_category,option)
  ! 
  ! Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! Notes: Updated/modified by Jenn Frederick 1/21/2016
  ! 

  use Option_module
  
  implicit none
   
  ! a maximum of 3 unit categories are allowed
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH) :: unit(3)
  character(len=*) :: units_category
  character(len=MAXSTRINGLENGTH) :: unit_category_string
  character(len=MAXWORDLENGTH) :: unit_category(3)
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: units_string
  PetscInt :: length, ind_dash, ind_space, ind 
  PetscInt :: k, j
  PetscInt :: num_units, num_unit_cats
  PetscBool :: category_assigned(3), successful
  PetscReal :: UnitsConvert

  UnitsConvert = 1d0
  unit_category_string = ''
  unit_category(:) = 'not_assigned'
  unit(:) = 'not_assigned'
  category_assigned(:) = PETSC_FALSE
  num_units = 0
  num_unit_cats = 0

  ind_dash = 1
  ! separate out unit categories (separated by "-")
  k = 1
  do while (ind_dash /= 0)
    length = len_trim(units_category)
    ind_dash = index(trim(units_category),"-")
    if (ind_dash == 0) then
      unit_category(k) = trim(units_category)
    else
      unit_category(k) = units_category(1:(ind_dash-1))
      units_category = units_category((ind_dash+1):length)
    endif
    k = k + 1
    if (k .gt. 3) then
      option%io_buffer = 'Maximum number of unit categories exceeded.&
                         & Reduce to no more than 3 units.' 
      call printErrMsg(option)
    endif
  enddo
  num_unit_cats = k - 1

  ind_dash = 1
  ind_space = 1
  ind = 1
  ! separate out units (can be separated by "-" or " ")
  k = 1
  do while (ind /= 0)
    length = len_trim(units)
    ind_dash = index(units,"-")
    ind_space = index(units," ")
    if ((ind_dash == 0) .or. (ind_space == 0)) then
      ind = max(ind_dash,ind_space)
    else
      ind = min(ind_dash,ind_space)
    endif
    if (ind == 0) then
      unit(k) = trim(units)
    else
      unit(k) = units(1:(ind-1))
      units = units((ind+1):length)
    endif
    k = k + 1
    if (k .gt. 3) then
      option%io_buffer = 'Maximum number of units exceeded.&
                         & Reduce to no more than 3 units.' 
      call printErrMsg(option)
    endif
  enddo
  num_units = k - 1

  ! check if expected # of units matches given # of units
  if (num_unit_cats /= num_units) then
    option%io_buffer = 'Mismatch between the number of units expected and the &
                       &number of units given.'
    call printErrMsg(option)
  endif

  k = 1
  do while (k < (num_units + 1))

    select case(trim(unit(k))) !-------------------------------------------------

      !---> VOLUME ---> (m^3)
      case('cm^3','l','dm^3','m^3')
        unit_category_string = 'volume'
        units_string = 'cm^3, l, dm^3, m^3'
        select case(trim(unit(k)))
          case('cm^3')
            UnitsConvert = 1.d-6
          case('l','dm^3')
            UnitsConvert = 1.d-3
          case('m^3')
            UnitsConvert = 1.d0
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      !---> AREA ---> (m^2)
      case('cm^2','dm^2','m^2')
        unit_category_string = 'area'
        units_string = 'cm^2, dm^2, m^2'
        select case(trim(unit(k)))
          case('cm^2')
            UnitsConvert = 1.d-4
          case('dm^2')
            UnitsConvert = 1.d-2
          case('m^2')
            UnitsConvert = 1.d0
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> LENGTH ---> (m)
      case('km','m','met','meter','dm','cm','mm')
        unit_category_string = 'length'
        units_string = 'km, m, met, meter, dm, cm, mm'
        select case(trim(unit(k)))
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
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> TIME ---> (sec)
      case('s','sec','second','min','minute','h','hr','hour','d','day','w', &
           'week','mo','month','y','yr','year')
        unit_category_string = 'time'
        units_string = 's, sec, second, min, minute, h, hr, hour, d, day, w, &
                       &week, mo, month, y, yr, year'
        select case(trim(unit(k)))
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
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> ENERGY ---> (MJ, MW)
      case('J','KJ','MJ','W','KW','MW')
        unit_category_string = 'energy'
        units_string = 'J, KJ, MJ, W, KW, MW'
        select case(trim(unit(k)))
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
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> MASS ---> (kg, mol)
      case('mol','mole','moles','ug','mg','g','kg')
        unit_category_string = 'mass'
        units_string = 'mol, mole, moles, ug, mg, g, kg'
        select case(trim(unit(k)))
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
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> TEMPERATURE ---> (C)
      case('C')
        unit_category_string = 'temperature'
        units_string = 'C (Celcius)'
        select case(trim(unit(k)))
          case('C') 
            UnitsConvert = 1.d0
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> PRESSURE ---> (Pa)
      case('Pa','KPa','MPa')
        unit_category_string = 'pressure'
        units_string = 'Pa, KPa, MPa'
        select case(trim(unit(k)))
          case('Pa') 
            UnitsConvert = 1.d0
          case('KPa')
            UnitsConvert = 1.d3
          case('MPa')
            UnitsConvert = 1.d6  
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
        end select
      ! ---> CONCENTRATION ---> (M)
      case('M')
        unit_category_string = 'concentration'
        units_string = 'M'
        select case(trim(unit(k)))
          case('M') 
            UnitsConvert = 1.d0 
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
      end select
      ! ---> FORCE ---> (N)
      case('N')
        unit_category_string = 'force'
        units_string = 'N (Newton)'
        select case(trim(unit(k)))
          case('N') 
            UnitsConvert = 1.d0
          case default
            call UnitsError(unit(k),unit_category_string,units_string,option)
      end select
      ! ---> SATURATION ---> (1)
      case('saturation')
        unit_category_string = 'saturation'
        UnitsConvert = 1.d0
      ! ---> UNITLESS ---> (1)
      case('unitless')
        unit_category_string = 'unitless'
        UnitsConvert = 1.d0
      ! ---> ERROR No unit category has been assigned during pass!
      case('not_assigned')
        unit_category_string = 'not_assigned'
        option%io_buffer = 'Unit category has not been assigned.'
        call printErrMsg(option)
      ! ---> ERROR Unit category unknown!
      case('unknown')
        unit_category_string = 'unknown'
        option%io_buffer = 'Unit category set to unknown.'
        call printErrMsg(option)
      ! ---> ERROR Cannot convert user-supplied units!
      case default
        option%io_buffer = 'Unit "' // trim(unit(k)) // '" is not a recognized &
                            &unit under any unit category.'
        call printErrMsg(option)

    end select ! unit(k) --------------------------------------------------------
                                                                                
    successful = PETSC_FALSE
    j = 1
    do while (j < (num_unit_cats + 1))
      if ((unit_category(j) == trim(unit_category_string)) .and. &
          (.not. category_assigned(j))) then
        category_assigned(j) = PETSC_TRUE
        successful = PETSC_TRUE
      endif
      if (successful) exit ! after first successful assignment
      j = j + 1
    enddo
    if (.not. successful) then
      option%io_buffer = 'Units of "' // trim(unit_category_string) // '" were &
                          &given, but not expected.'
      call printErrMsg(option)
    endif
    k = k + 1
  enddo

end function UnitsConvert

! ************************************************************************** !

subroutine UnitsError(unit,units_category,units_string,option)
  ! 
  ! Forms and prints a units conversion error message.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/21/2016
  ! 

  use Option_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: unit
  character(len=MAXWORDLENGTH) :: units_category
  character(len=MAXSTRINGLENGTH) :: units_string
  type(option_type) :: option

  option%io_buffer = 'Unit "' // trim(unit) // '" is not a recognized unit&
                     & of "' // trim(units_category) // '". Recognized units&
                     & include: "' // trim(units_string) // '" only.'
  call printErrMsg(option)
  
end subroutine UnitsError

end module Units_module
