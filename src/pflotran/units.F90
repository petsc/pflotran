module Units_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"
  
  public :: UnitsConvertToInternal, UnitsConvertToExternal
  
contains

! ************************************************************************** !

function UnitsConvertToInternal(units_input,units_category,option)
  ! 
  ! Converts units to pflotran internal units
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! Notes: Updated/modified by Jenn Frederick 1/25/2016
  ! 

  use Option_module

  implicit none
  
  character(len=*) :: units_input
  character(len=*) :: units_category
  character(len=MAXSTRINGLENGTH) :: ierr_msg
  character(len=MAXSTRINGLENGTH) :: units_input_buff
  character(len=MAXSTRINGLENGTH) :: units_category_buff1
  character(len=MAXSTRINGLENGTH) :: units_category_buff2
  type(option_type) :: option
  PetscBool :: ierr, multi_option
  PetscInt :: length, ind_or, num_options
  PetscReal :: conversion_factor, UnitsConvertToInternal

  units_input_buff = trim(units_input)
  units_category_buff1 = trim(units_category)
  units_category_buff2 = trim(units_category)
  ierr = PETSC_FALSE
  ierr_msg = ''
  multi_option = PETSC_FALSE
  conversion_factor = 1d0
  UnitsConvertToInternal = 1d0
  num_options = 0
  ind_or = 1
  length = 0

  do while(ind_or /= 0)
    length = len_trim(units_category_buff1)
    ind_or = index(trim(units_category_buff1),"|")   
    if (ind_or == 0) then
      call UnitsConvertParse(units_input_buff,units_category_buff1, &
                             conversion_factor,ierr,ierr_msg)
    else
      multi_option = PETSC_TRUE
      num_options = num_options + 1
      units_category_buff2 = units_category_buff1(1:(ind_or-1))
      call UnitsConvertParse(units_input_buff,units_category_buff2, &
                             conversion_factor,ierr,ierr_msg)
      units_category_buff1 = units_category_buff1((ind_or+1):length)
    endif
    if (.not. ierr) exit
  enddo

  if (ierr) then
    option%io_buffer = 'While converting unit "' // trim(units_input) // &
                       '", under unit category "' // trim(units_category) // &
                       '" :: "' // trim(ierr_msg) // '" '
    call printErrMsg(option)
  endif

  UnitsConvertToInternal = conversion_factor

end function UnitsConvertToInternal

! ************************************************************************** !

subroutine UnitsConvertParse(units,units_category,conversion_factor,ierr, &
                             ierr_msg)
  ! 
  ! Parses unit and unit category numerators from denominators (if they exist)
  ! and gets the conversion factor to internal Pflotran units
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/25/2016
  ! 

  use Option_module

  implicit none

  character(len=*) :: units, units_category
  character(len=MAXWORDLENGTH) :: numerator_units
  character(len=MAXWORDLENGTH) :: denominator_units
  ! A maximum of 3 numerators/denominators are accepted 
  character(len=MAXSTRINGLENGTH) :: numerator_unit_cats
  character(len=MAXSTRINGLENGTH) :: denominator_unit_cats
  character(len=MAXSTRINGLENGTH) :: ierr_msg
  PetscReal :: numerator_conv_factor, denominator_conv_factor
  PetscReal :: conversion_factor
  PetscBool :: set_units_denom, set_cats_denom, ierr
  PetscInt :: length, ind
  
  conversion_factor = 1.d0
  ierr = PETSC_FALSE
  set_units_denom = PETSC_FALSE
  set_cats_denom = PETSC_FALSE
  numerator_conv_factor = 1.d0
  denominator_conv_factor = 1.d0
  ierr_msg = ''
  numerator_units = ''
  denominator_units = ''
  numerator_unit_cats = 'not_assigned'
  denominator_unit_cats = 'not_assigned'

  length = len(units_category)
  ind = index(trim(units_category),"/")   
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
  ind = index(trim(units),"/")
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
    ierr = PETSC_TRUE
    if ((set_cats_denom .eqv. PETSC_TRUE) .and. &
        (set_units_denom .eqv. PETSC_FALSE)) then
      ierr_msg = 'Units provided do not match the expected &
                 &unit category structure (numerator/denominator). &
                 &A unit denominator was expected, but was not given.'
    endif
    if ((set_cats_denom .eqv. PETSC_FALSE) .and. &
        (set_units_denom .eqv. PETSC_TRUE)) then
      ierr_msg = 'Units provided do not match the expected &
                 &unit category structure (numerator/denominator). &
                 &A unit denominator was not expected, but was given.'
    endif
  endif

  if (.not. ierr) then
  call UnitsConvert(numerator_units, numerator_unit_cats, &
                    numerator_conv_factor, ierr, ierr_msg)
  endif

  if (.not. ierr) then
    if (set_cats_denom) then
       call UnitsConvert(denominator_units, denominator_unit_cats, &
                         denominator_conv_factor, ierr, ierr_msg)
    endif
  endif

  if (.not. ierr) then
    conversion_factor = numerator_conv_factor/denominator_conv_factor
  endif

end subroutine UnitsConvertParse

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
  
  character(len=MAXSTRINGLENGTH) :: units
  character(len=*) :: units_category
  type(option_type) :: option
  
  PetscReal :: UnitsConvertToExternal
  
  UnitsConvertToExternal = 1.d0/UnitsConvertToInternal(units, &
                                                       units_category,option)

end function UnitsConvertToExternal

! ************************************************************************** !

subroutine UnitsConvert(units,units_category,conversion_factor,ierr,ierr_msg)
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
  character(len=*) :: units
  character(len=*) :: units_category
  character(len=MAXWORDLENGTH) :: unit(3)
  character(len=MAXWORDLENGTH) :: unit_category(3)
  character(len=MAXSTRINGLENGTH) :: units_buff
  character(len=MAXSTRINGLENGTH) :: units_category_buff
  character(len=MAXSTRINGLENGTH) :: unit_category_string
  character(len=MAXSTRINGLENGTH) :: units_string, ierr_msg
  type(option_type) :: option
  PetscInt :: length, ind_dash, ind_space, ind 
  PetscInt :: k, j
  PetscInt :: num_units, num_unit_cats
  PetscBool :: category_assigned(3), successful, ierr
  PetscReal :: conversion_factor

  conversion_factor = 1d0
  ierr_msg = ''
  ierr = PETSC_FALSE
  unit_category_string = ''
  units_buff = trim(units)
  units_category_buff = trim(units_category)
  unit_category(:) = 'not_assigned'
  unit(:) = 'not_assigned'
  category_assigned(:) = PETSC_FALSE
  num_units = 0
  num_unit_cats = 0

  ind_dash = 1
  ! separate out unit categories (separated by "-")
  k = 1
  do while (ind_dash /= 0)
    length = len_trim(units_category_buff)
    ind_dash = index(trim(units_category_buff),"-")
    if (ind_dash == 0) then
      unit_category(k) = trim(units_category_buff)
    else
      unit_category(k) = units_category_buff(1:(ind_dash-1))
      units_category_buff = units_category_buff((ind_dash+1):length)
    endif
    k = k + 1
    if (k .gt. 3) then
      ierr_msg = 'Maximum number of unit categories exceeded.&
                 & Reduce to no more than 3 units.' 
      ierr = PETSC_TRUE
      exit
    endif
  enddo
  num_unit_cats = k - 1

  ind_dash = 1
  ind_space = 1
  ind = 1
  ! separate out units (can be separated by "-" or " ")
  k = 1
  do while (ind /= 0)
    length = len_trim(units_buff)
    ind_dash = index(trim(units_buff),"-")
    ind_space = index(trim(units_buff)," ")
    if ((ind_dash == 0) .or. (ind_space == 0)) then
      ind = max(ind_dash,ind_space)
    else
      ind = min(ind_dash,ind_space)
    endif
    if (ind == 0) then
      unit(k) = trim(units_buff)
    else
      unit(k) = units_buff(1:(ind-1))
      units_buff = units_buff((ind+1):length)
    endif
    k = k + 1
    if (k .gt. 3) then
      ierr_msg = 'Maximum number of units exceeded.&
                 & Reduce to no more than 3 units.' 
      ierr = PETSC_TRUE
      exit
    endif
  enddo
  num_units = k - 1

  ! check if expected # of units matches given # of units
  if (num_unit_cats /= num_units) then
    ierr_msg = 'Mismatch between the number of units expected and the &
               &number of units given.'
    ierr = PETSC_TRUE
  endif

  if (.not. ierr) then

    k = 1
    do while (k < (num_units + 1))

      select case(trim(unit(k))) !-------------------------------------------------

      !---> VOLUME ---> (m^3)
      case('cm^3','l','dm^3','m^3')
        unit_category_string = 'volume'
        units_string = 'cm^3, l, dm^3, m^3'
        select case(trim(unit(k)))
        case('cm^3')
          conversion_factor = 1.d-6
        case('l','dm^3')
          conversion_factor = 1.d-3
        case('m^3')
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      !---> AREA ---> (m^2)
      case('cm^2','dm^2','m^2')
        unit_category_string = 'area'
        units_string = 'cm^2, dm^2, m^2'
        select case(trim(unit(k)))
        case('cm^2')
          conversion_factor = 1.d-4
        case('dm^2')
          conversion_factor = 1.d-2
        case('m^2')
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> LENGTH ---> (m)
      case('km','m','met','meter','dm','cm','mm')
        unit_category_string = 'length'
        units_string = 'km, m, met, meter, dm, cm, mm'
        select case(trim(unit(k)))
        case('km')
          conversion_factor = 1000.d0
        case('m','met','meter')
          conversion_factor = 1.d0
        case('dm')
          conversion_factor = 1.d-1
        case('cm')
          conversion_factor = 1.d-2
        case('mm')
          conversion_factor = 1.d-3
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> TIME ---> (sec)
      case('s','sec','second','min','minute','h','hr','hour','d','day','w', &
           'week','mo','month','y','yr','year')
        unit_category_string = 'time'
        units_string = 's, sec, second, min, minute, h, hr, hour, d, day, w, &
                       &week, mo, month, y, yr, year'
        select case(trim(unit(k)))
        case('s','sec','second')
          conversion_factor = 1.d0
        case('min','minute')
          conversion_factor = 60.d0
        case('h','hr','hour')
          conversion_factor = 3600.d0
        case('d','day')
          conversion_factor = 24.d0*3600.d0 
        case('w','week')
          conversion_factor = 7.d0*24.d0*3600.d0 
        case('mo','month')
          conversion_factor = 365.d0/12.d0*24.d0*3600.d0 
        case('y','yr','year')
          conversion_factor = 365.d0*24.d0*3600.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> ENERGY ---> (MJ)
      case('J','KJ','MJ')
        unit_category_string = 'energy'
        units_string = 'J, KJ, MJ'
        select case(trim(unit(k)))
        case('J')   
          conversion_factor = 1.d-6
        case('KJ')   
          conversion_factor = 1.d-3
        case('MJ')   
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> ENERGY FLUX or POWER ---> (MW)
      case('W','KW','MW')
        unit_category_string = 'power'
        units_string = 'W, KW, MW'
        select case(trim(unit(k)))
        case('W')   
          conversion_factor = 1.d-6
        case('KW')   
          conversion_factor = 1.d-3
        case('MW')   
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> MASS ---> (kg, mol)
      case('mol','mole','moles','ug','mg','g','kg')
        unit_category_string = 'mass'
        units_string = 'mol, mole, moles, ug, mg, g, kg'
        select case(trim(unit(k)))
        case('mol','mole','moles')
          conversion_factor = 1.d0
        case('ug')
          conversion_factor = 1.d-9
        case('mg')
          conversion_factor = 1.d-6
        case('g')
          conversion_factor = 1.d-3
        case('kg')
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> TEMPERATURE ---> (C)
      case('C')
        unit_category_string = 'temperature'
        units_string = 'C (Celcius)'
        select case(trim(unit(k)))
        case('C') 
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> PRESSURE ---> (Pa)
      case('Pa','KPa','MPa')
        unit_category_string = 'pressure'
        units_string = 'Pa, KPa, MPa'
        select case(trim(unit(k)))
        case('Pa') 
          conversion_factor = 1.d0
        case('KPa')
          conversion_factor = 1.d3
        case('MPa')
          conversion_factor = 1.d6  
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> CONCENTRATION ---> (M)
      case('M')
        unit_category_string = 'concentration'
        units_string = 'M'
        select case(trim(unit(k)))
        case('M') 
          conversion_factor = 1.d0 
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> FORCE ---> (N)
      case('N')
        unit_category_string = 'force'
        units_string = 'N (Newton)'
        select case(trim(unit(k)))
        case('N') 
          conversion_factor = 1.d0
        case default
          call UnitsError(unit(k),unit_category_string,units_string,ierr,ierr_msg)
        end select
      ! ---> SATURATION ---> (1)
      case('saturation')
        unit_category_string = 'saturation'
        conversion_factor = 1.d0
      ! ---> UNITLESS ---> (1)
      case('unitless')
        unit_category_string = 'unitless'
        conversion_factor = 1.d0
      ! ---> ERROR No unit category has been assigned during pass!
      case('not_assigned')
        unit_category_string = 'not_assigned'
        ierr_msg = 'Unit category has not been assigned.'
        ierr = PETSC_TRUE
      ! ---> ERROR Unit category unknown!
      case('unknown')
        unit_category_string = 'unknown'
        ierr_msg = 'Unit category set to unknown.'
        ierr = PETSC_TRUE
      ! ---> ERROR Cannot convert user-supplied units!
      case default
        ierr_msg = 'Unit "' // trim(unit(k)) // '" is not a recognized &
                    &unit under any unit category.'
        ierr = PETSC_TRUE

      end select ! unit(k) --------------------------------------------------------

      if (.not. ierr) then
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
          ierr_msg = 'Units of "' // trim(unit_category_string) // '" were &
                     &given, but not expected.'
          ierr = PETSC_TRUE
        endif
        k = k + 1
      endif ! (.not. ierr)

    enddo ! while (k < (num_units + 1))

  endif ! (.not. ierr)

end subroutine UnitsConvert

! ************************************************************************** !

subroutine UnitsError(unit,units_category,units_string,ierr,ierr_msg)
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
  character(len=MAXSTRINGLENGTH) :: units_string, ierr_msg
  PetscBool :: ierr

  ierr = PETSC_TRUE

  ierr_msg = 'Unit "' // trim(unit) // '" is not a recognized unit&
             & of "' // trim(units_category) // '". Recognized units&
             & include: "' // trim(units_string) // '" only.'
  
end subroutine UnitsError

end module Units_module
