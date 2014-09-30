module Characteristic_Curves_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  PetscReal, parameter :: DEFAULT_PCMAX = 1.d9

  type :: polynomial_type
    PetscReal :: low
    PetscReal :: high
    PetscReal :: coefficients(4)
  end type polynomial_type

  ! Saturation Function
  type :: sat_func_base_type
    type(polynomial_type), pointer :: sat_poly
    type(polynomial_type), pointer :: pres_poly
    PetscReal :: Sr
    PetscReal :: pcmax
  contains
    procedure, public :: Init => SFBaseInit
    procedure, public :: Verify => SFBaseVerify
    procedure, public :: Test => SFBaseTest
    procedure, public :: SetupPolynomials => SFBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFBaseCapillaryPressure
    procedure, public :: Saturation => SFBaseSaturation
  end type sat_func_base_type
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_VG_Init
    procedure, public :: Verify => SF_VG_Verify
    procedure, public :: CapillaryPressure => SF_VG_CapillaryPressure
    procedure, public :: Saturation => SF_VG_Saturation
  end type sat_func_VG_type  
  type, public, extends(sat_func_base_type) :: sat_func_BC_type
    PetscReal :: alpha
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_BC_Init
    procedure, public :: Verify => SF_BC_Verify
    procedure, public :: SetupPolynomials => SF_BC_SetupPolynomials
    procedure, public :: CapillaryPressure => SF_BC_CapillaryPressure
    procedure, public :: Saturation => SF_BC_Saturation
  end type sat_func_BC_type  
  
  ! Relative Permeability Function
  type :: rel_perm_func_base_type
    type(polynomial_type), pointer :: poly
    PetscReal :: Sr
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: Verify => RPFBaseVerify
    procedure, public :: Test => RPF_BAse_Test
    procedure, public :: SetupPolynomials => RPFBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Mualem_Init
    procedure, public :: Verify => RPF_Mualem_Verify
    procedure, public :: SetupPolynomials => RPF_Mualem_SetupPolynomials
    procedure, public :: RelativePermeability => RPF_Mualem_RelPerm
  end type rpf_Mualem_type
  type, public, extends(rpf_Mualem_type) :: rpf_Mualem_VG_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_VG_Gas_Init
    procedure, public :: Verify => RPF_Mualem_VG_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Gas_RelPerm
  end type rpf_Mualem_VG_gas_type
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_Burdine_Init
    procedure, public :: Verify => RPF_Burdine_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_RelPerm
  end type rpf_Burdine_type
  type, public, extends(rpf_Burdine_type) :: rpf_Burdine_BC_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_BC_Gas_Init
    procedure, public :: Verify => RPF_Burdine_BC_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_BC_Gas_RelPerm
  end type rpf_Burdine_BC_gas_type
  ! since the TOUGH2_Corey relative permeability function (IRP=7 in TOUGH2
  ! manual) calculates relative perm as a function of the Mualem-based liquid
  ! relative permeability when Srg = 0., we extend the rpf_Mualem_type to
  ! save code
  type, public, extends(rpf_Mualem_type) :: rpf_TOUGH2_IRP7_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_TOUGH2_IRP7_Gas_Init
    procedure, public :: Verify => RPF_TOUGH2_IRP7_Gas_Verify
    procedure, public :: RelativePermeability => RPF_TOUGH2_IRP7_Gas_RelPerm
  end type rpf_TOUGH2_IRP7_gas_type

  type, public :: characteristic_curves_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(sat_func_base_type), pointer :: saturation_function
    class(rel_perm_func_base_type), pointer :: liq_rel_perm_function
    class(rel_perm_func_base_type), pointer :: gas_rel_perm_function
    class(characteristic_curves_type), pointer :: next
  end type characteristic_curves_type
  
  type, public :: characteristic_curves_ptr_type
    class(characteristic_curves_type), pointer :: ptr
  end type characteristic_curves_ptr_type 
  
  public :: CharacteristicCurvesCreate, &
            CharacteristicCurvesRead, &
            CharacteristicCurvesAddToList, &
            CharCurvesConvertListToArray, &
            CharacteristicCurvesGetID, &
            CharCurvesGetGetResidualSats, &
            CharacteristicCurvesDestroy
  
contains

! ************************************************************************** !

function CharacteristicCurvesCreate()
  ! 
  ! Creates a characteristic curve object that holds parameters and pointers
  ! to functions for calculating saturation, capillary pressure, relative
  ! permeability, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/23/14
  ! 

  implicit none

  class(characteristic_curves_type), pointer :: CharacteristicCurvesCreate
  
  class(characteristic_curves_type), pointer :: characteristic_curves
  
  allocate(characteristic_curves)
  characteristic_curves%name = ''
  characteristic_curves%print_me = PETSC_FALSE
  characteristic_curves%test = PETSC_FALSE
  nullify(characteristic_curves%saturation_function)
  nullify(characteristic_curves%liq_rel_perm_function)
  nullify(characteristic_curves%gas_rel_perm_function)
  nullify(characteristic_curves%next)

  CharacteristicCurvesCreate => characteristic_curves

end function CharacteristicCurvesCreate

! ************************************************************************** !

subroutine CharacteristicCurvesRead(this,input,option)
  ! 
  ! Reads in contents of a saturation_function card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(characteristic_curves_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: iphase
  
  character(len=MAXWORDLENGTH) :: keyword, word, phase_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  class(rel_perm_func_base_type), pointer :: rel_perm_function_ptr

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES'  
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('SATURATION_FUNCTION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation_function_type', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('VAN_GENUCHTEN')
            this%saturation_function => SF_VG_Create()
          case('BROOKS_COREY')
            this%saturation_function => SF_BC_Create()
          case default
            option%io_buffer = 'Keyword: ' // trim(word) // &
              ' not a recognized in saturation function type.'    
            call printErrMsg(option)            
        end select
        call SaturationFunctionRead(this%saturation_function,input,option)
      case('PERMEABILITY_FUNCTION')
        nullify(rel_perm_function_ptr)
        phase_keyword = 'NONE'
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'permeability_function_type', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MUALEM')
            rel_perm_function_ptr => RPF_Mualem_Create()
          case('MUALEM_VG_GAS')
            rel_perm_function_ptr => RPF_Mualem_VG_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE')
            rel_perm_function_ptr => RPF_Burdine_Create()
          case('BURDINE_BC_GAS')
            rel_perm_function_ptr => RPF_Burdine_BC_Gas_Create()
            phase_keyword = 'GAS'
          case('TOUGH2_IRP7_GAS')
            rel_perm_function_ptr => RPF_TOUGH2_IRP7_Gas_Create()
            phase_keyword = 'GAS'
          case default
            option%io_buffer = 'Keyword: ' // trim(word) // &
              ' not a recognized in relative permeability function type.'
            call printErrMsg(option)            
        end select
        call PermeabilityFunctionRead(rel_perm_function_ptr,phase_keyword, &
                                      input,option)
        ! if PHASE is specified, have to align correct pointer
        select case(phase_keyword)
          case('GAS')
            this%gas_rel_perm_function => rel_perm_function_ptr
          case('LIQUID')
            this%liq_rel_perm_function => rel_perm_function_ptr
          case('NONE')
            this%gas_rel_perm_function => rel_perm_function_ptr
            this%liq_rel_perm_function => rel_perm_function_ptr
          case default
            option%io_buffer = 'PHASE keyword: ' // trim(word) // &
              ' not a recognized in relative permeability function type.'
            call printErrMsg(option)            
        end select
      case('TEST') 
        this%test = PETSC_TRUE
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in charateristic_curves.'    
        call printErrMsg(option)
    end select 
  enddo
  
  call CharacteristicCurvesVerify(this,option)

end subroutine CharacteristicCurvesRead

! ************************************************************************** !

subroutine SaturationFunctionRead(saturation_function,input,option)

  ! Reads in contents of a SATURATION_FUNCTION block

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(sat_func_base_type) :: saturation_function
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION,'
  select type(sf => saturation_function)
    class is(sat_func_VG_type)
      error_string = trim(error_string) // 'VAN_GENUCHTEN'
    class is(sat_func_BC_type)
      error_string = trim(error_string) // 'BROOKS_COREY'
  end select
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,saturation_function%Sr)
        call InputErrorMsg(input,option,'liquid residual saturation', &
                           error_string)
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcmax)
        call InputErrorMsg(input,option,'maximum capillary pressure', &
                            error_string)
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    
    if (found) cycle
    
    select type(sf => saturation_function)
      class is(sat_func_VG_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in van Genuchten saturation function.'
            call printErrMsg(option)
        end select
      class is(sat_func_BC_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Brooks-Corey saturation function.'
            call printErrMsg(option)
        end select
      class default
        option%io_buffer = 'Read routine not implemented for saturation ' // &
                           'function class.'
        call printErrMsg(option)
    end select
  enddo
  
  if (smooth) then
    call saturation_function%SetupPolynomials(option,error_string)
  endif

end subroutine SaturationFunctionRead

! ************************************************************************** !

subroutine PermeabilityFunctionRead(permeability_function,phase_keyword, &
                                    input,option)

  ! Reads in contents of a PERMEABILITY_FUNCTION block

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(rel_perm_func_base_type) :: permeability_function
  character(len=MAXWORDLENGTH) :: phase_keyword
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, new_phase_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  new_phase_keyword = 'NONE'
  error_string = 'CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION,'
  select type(rpf => permeability_function)
    class is(rpf_Mualem_type)
      error_string = trim(error_string) // 'MUALEM'
    class is(rpf_Mualem_VG_gas_type)
      error_string = trim(error_string) // 'MUALEM_VG_GAS'
    class is(rpf_Burdine_type)
      error_string = trim(error_string) // 'BURDINE'
    class is(rpf_Burdine_BC_gas_type)
      error_string = trim(error_string) // 'BURDINE_BC_GAS'
  end select
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,permeability_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('PHASE')
        call InputReadWord(input,option,new_phase_keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'phase',error_string)
        call StringToUpper(phase_keyword) 
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    
    if (found) cycle
    
    ! we assume liquid phase if PHASE keyword is not present.
    select type(rpf => permeability_function)
      class is(rpf_Mualem_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Mualem van Genuchten relative ' // &
              'permeability function.'
            call printErrMsg(option)
        end select
      class is(rpf_Mualem_VG_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Mualem van Genuchten gas relative ' // &
              'permeability function.'
            call printErrMsg(option)
        end select
      class is(rpf_Burdine_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Burdine Brooks-Corey relative ' // &
              'permeability function.'
            call printErrMsg(option)
        end select
      class is(rpf_Burdine_BC_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Burdine Brooks-Corey gas relative ' // &
              'permeability function.'
            call printErrMsg(option)
        end select
      class is(rpf_TOUGH2_IRP7_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in TOUGH2 IRP7 gas relative ' // &
              'permeability function.'
            call printErrMsg(option)
        end select
      class default
        option%io_buffer = 'Read routine not implemented for relative ' // &
                           'permeability function class.'
        call printErrMsg(option)
    end select
  enddo
  
  ! check to ensure that the phase is correct if phase_keyword was set to 
  ! something other than 'NONE' prior to the call of this subroutine
  if (StringCompare('NONE',phase_keyword)) then
    phase_keyword = new_phase_keyword
  else if (.not.StringCompare('NONE',new_phase_keyword)) then
    if (.not.StringCompare(phase_keyword,new_phase_keyword)) then
      option%io_buffer = 'Relative permeability function has been set ' // &
        'for the wrong phase (' // trim(phase_keyword) // ' vs ' // &
        trim(new_phase_keyword) // ').'
      call printErrMsg(option)
    endif
  endif

  if (smooth) then
    call permeability_function%SetupPolynomials(option,error_string)
  endif
  
end subroutine PermeabilityFunctionRead

! ************************************************************************** !

subroutine CharacteristicCurvesAddToList(new_characteristic_curves,list)
  ! 
  ! Adds a characteristic curves object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: new_characteristic_curves
  class(characteristic_curves_type), pointer :: list

  class(characteristic_curves_type), pointer :: cur_characteristic_curves
  
  if (associated(list)) then
    cur_characteristic_curves => list
    ! loop to end of list
    do
      if (.not.associated(cur_characteristic_curves%next)) exit
      cur_characteristic_curves => cur_characteristic_curves%next
    enddo
    cur_characteristic_curves%next => new_characteristic_curves
  else
    list => new_characteristic_curves
  endif
  
end subroutine CharacteristicCurvesAddToList

! ************************************************************************** !

subroutine CharCurvesConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the characteristic curves objects in the 
  ! list
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/11/07
  ! 

  use String_module
  use Option_module
  
  implicit none
  
  class(characteristic_curves_type), pointer :: list
  type(characteristic_curves_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  class(characteristic_curves_type), pointer :: cur_characteristic_curves
  PetscInt :: count

  count = 0
  cur_characteristic_curves => list
  do 
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo
  
  if(associated(array)) deallocate(array)
  allocate(array(count))
  
  count = 0
  cur_characteristic_curves => list
  do 
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    array(count)%ptr => cur_characteristic_curves
    if (cur_characteristic_curves%test .and. &
        option%myrank == option%io_rank) then
      call CharacteristicCurvesTest(cur_characteristic_curves,option)
    endif
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

end subroutine CharCurvesConvertListToArray

! ************************************************************************** !

function CharCurvesGetGetResidualSats(characteristic_curves,option)
  ! 
  ! Returns the residual saturations associated with a characteristic curves
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  ! 

  use Option_module
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  PetscReal :: CharCurvesGetGetResidualSats(option%nphase)

  CharCurvesGetGetResidualSats(1) = &
    characteristic_curves%liq_rel_perm_function%Sr
  if (option%nphase > 1) then
    select type(rpf=>characteristic_curves%gas_rel_perm_function)
      class is(rpf_Mualem_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Burdine_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Mualem_VG_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Burdine_BC_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_TOUGH2_IRP7_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class default
        option%io_buffer = 'Relative permeability class not supported in ' // &
          'CharCurvesGetGetResidualSats.'
        call printErrMsg(option)
    end select
  endif

end function CharCurvesGetGetResidualSats

! ************************************************************************** !

function CharacteristicCurvesGetID(characteristic_curves_array, &
                                   characteristic_curves_name, &
                                   material_property_name, option)
  ! 
  ! Returns the ID of the characteristic curves object named
  ! "characteristic_curves_name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  type(characteristic_curves_ptr_type), pointer :: &
    characteristic_curves_array(:)
  character(len=MAXWORDLENGTH) :: characteristic_curves_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CharacteristicCurvesGetID

  CharacteristicCurvesGetID = 0
  do CharacteristicCurvesGetID = 1, size(characteristic_curves_array)
    if (StringCompare(characteristic_curves_name, &
                      characteristic_curves_array( &
                        CharacteristicCurvesGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Characteristic curves "' // &
           trim(characteristic_curves_name) // &
           '" in material property "' // &
           trim(material_property_name) // &
           '" not found among available characteristic curves.'
  call printErrMsg(option)    

end function CharacteristicCurvesGetID

! ************************************************************************** !

subroutine CharacteristicCurvesTest(characteristic_curves,option)
  ! 
  ! Outputs values of characteristic curves over a range of values
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: phase

  call characteristic_curves%saturation_function%Test( &
                                                 characteristic_curves%name, &
                                                 option)
  phase = 'liquid'
  call characteristic_curves%liq_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  phase = 'gas'
  call characteristic_curves%gas_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  
end subroutine CharacteristicCurvesTest

! ************************************************************************** !

subroutine CharacteristicCurvesVerify(characteristic_curves,option)
  ! 
  ! Outputs values of characteristic curves over a range of values
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  string = 'CHARACTERISTIC_CURVES(' // trim(characteristic_curves%name) // &
           '),'

  call characteristic_curves%saturation_function%Verify(string,option)
  call characteristic_curves%liq_rel_perm_function%Verify(string,option)
  call characteristic_curves%gas_rel_perm_function%Verify(string,option)
  
end subroutine CharacteristicCurvesVerify

!!!BASE ROUTINES
! ************************************************************************** !

function PolynomialCreate()

  implicit none
  
  type(polynomial_type), pointer :: PolynomialCreate  

  allocate(PolynomialCreate)
  PolynomialCreate%low = 0.d0
  PolynomialCreate%high = 0.d0
  PolynomialCreate%coefficients(:) = 0.d0
  
end function PolynomialCreate

! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%sat_poly)
  nullify(this%pres_poly)
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine SFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  
  
  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
end subroutine SFBaseVerify

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%poly)
  this%Sr = UNINITIALIZED_DOUBLE
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine RPFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
end subroutine RPFBaseVerify
  
! ************************************************************************** !

subroutine SFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine SFBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing relative permeability functions

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine RPFBaseSetupPolynomials

! ************************************************************************** !

subroutine SFBaseCapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,option)
  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseCapillaryPressure must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseCapillaryPressure

! ************************************************************************** !

subroutine SFBaseSaturation(this,capillary_pressure,liquid_saturation, &
                              dsat_pres,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseSaturation must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseSaturation

! ************************************************************************** !

subroutine SFBaseTest(this,cc_name,option)

  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: pc, pc_increment
  PetscReal :: capillary_pressure(101)
  PetscReal :: liquid_saturation(101)
  PetscReal :: dummy_real
  PetscInt :: count, i

  ! calculate saturation as a function of capillary pressure
  ! start at 1 Pa up to maximum capillary pressure
  pc = 1.d0
  pc_increment = 1.d0
  count = 0
  do
    if (pc > this%pcmax) exit
    count = count + 1
    call this%Saturation(pc,liquid_saturation(count),dummy_real,option)
    capillary_pressure(count) = pc
    if (pc > 0.99*pc_increment*10.d0) pc_increment = pc_increment*10.d0
    pc = pc + pc_increment
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation"'
  do i = 1, count
    write(86,'(2es14.6)') capillary_pressure(i), liquid_saturation(i)
  enddo
  close(86)

 ! calculate capillary pressure as a function of saturation
  do i = 1, 101
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%CapillaryPressure(liquid_saturation(i),capillary_pressure(i), &
                                option)
  enddo
  count = 101

  write(string,*) cc_name
  string = trim(cc_name) // '_sat_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure"'
  do i = 1, count
    write(86,'(2es14.6)') liquid_saturation(i), capillary_pressure(i)
  enddo
  close(86)

end subroutine SFBaseTest

! ************************************************************************** !

subroutine RPF_Base_RelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_Se,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPF_Base_RelPerm must be extended.'
  call printErrMsg(option)
  
end subroutine RPF_Base_RelPerm

! ************************************************************************** !

subroutine RPF_Base_Test(this,cc_name,phase,option)

  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: dummy_real
  PetscInt :: i
  PetscReal :: liquid_saturation(101), kr(101)

  do i = 1, 101
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%RelativePermeability(liquid_saturation(i),kr(i),dummy_real, &
                                   option)
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "' // trim(phase) // ' relative permeability"'
  do i = 1, size(liquid_saturation)
    write(86,'(2es14.6)') liquid_saturation(i), kr(i)
  enddo
  close(86)

end subroutine RPF_Base_Test

! ************************************************************************** !

function SF_VG_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_VG_type), pointer :: SF_VG_Create
  
  allocate(SF_VG_Create)
  call SF_VG_Create%Init()
  
end function SF_VG_Create

! ************************************************************************** !

subroutine SF_VG_Init(this)

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_VG_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine SF_VG_Init

! ************************************************************************** !

subroutine SF_VG_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,VAN_GENUCHTEN'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif   
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   

end subroutine SF_VG_Verify

! ************************************************************************** !

subroutine SF_VG_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  
  implicit none
  
  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: n
  PetscReal :: Se
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: pc_alpha_n
  PetscReal :: pc_alpha
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  n = 1.d0/(1.d0-this%m)
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  one_plus_pc_alpha_n = Se**(-1.d0/this%m)
  pc_alpha_n = one_plus_pc_alpha_n - 1.d0
  pc_alpha = pc_alpha_n**(1.d0/n)
  capillary_pressure = pc_alpha/this%alpha

  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_VG_CapillaryPressure

! ************************************************************************** !

subroutine SF_VG_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_pres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  PetscReal, parameter :: pc_alpha_n_epsilon = 1.d-15
  PetscReal :: n
  PetscReal :: pc_alpha
  PetscReal :: pc_alpha_n
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: Se
  PetscReal :: dSe_pc
  PetscReal :: dsat_pc
  
  dsat_pres = 0.d0
  
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_pc)
      liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
      dsat_pc = (1.d0-this%Sr)*dSe_pc
    endif
  endif

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    n = 1.d0/(1.d0-this%m)
    pc_alpha = capillary_pressure*this%alpha
    pc_alpha_n = pc_alpha**n
    !geh:  This conditional does not catch potential cancelation in 
    !      the dkr_Se deriviative calculation.  Therefore, I am setting
    !      an epsilon here
!        if (1.d0 + pc_alpha_n == 1.d0) then ! check for zero perturbation
    if (pc_alpha_n < pc_alpha_n_epsilon) then 
      liquid_saturation = 1.d0
!      switch_to_saturated = PETSC_TRUE
      return
    endif
    one_plus_pc_alpha_n = 1.d0+pc_alpha_n
    Se = one_plus_pc_alpha_n**(-this%m)
    dSe_pc = -this%m*n*this%alpha*pc_alpha_n/ &
            (pc_alpha*one_plus_pc_alpha_n**(this%m+1.d0))
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_pc = (1.d0-this%Sr)*dSe_pc
  endif
  
end subroutine SF_VG_Saturation

! ************************************************************************** !

function RPF_Mualem_Create()

  ! Creates the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_type), pointer :: RPF_Mualem_Create
  
  allocate(RPF_Mualem_Create)
  call RPF_Mualem_Create%Init()
  
end function RPF_Mualem_Create

! ************************************************************************** !

subroutine RPF_Mualem_Init(this)

  ! Initializes the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_Init

! ************************************************************************** !

subroutine RPF_Mualem_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   
  
end subroutine RPF_Mualem_Verify

! ************************************************************************** !

subroutine RPF_Mualem_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Mualem permeability function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(rpf_Mualem_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)
  PetscReal :: one_over_m, Se_one_over_m, m

  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  this%poly%high = 1.d0   ! saturated
  
  m = this%m
  one_over_m = 1.d0/m
  Se_one_over_m = this%poly%low**one_over_m
  b(1) = 1.d0
  b(2) = sqrt(this%poly%low)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
  b(3) = 0.d0
  b(4) = 0.5d0*b(2)/this%poly%low+ &
          2.d0*this%poly%low**(one_over_m-0.5d0)* &
          (1.d0-Se_one_over_m)**(m-1.d0)* &
          (1.d0-(1.d0-Se_one_over_m)**m)
  
  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  
  this%poly%coefficients(1:4) = b(1:4)
  
  
end subroutine RPF_Mualem_SetupPolynomials

! ************************************************************************** !

subroutine RPF_Mualem_RelPerm(this,liquid_saturation,relative_permeability, &
                              dkr_Se,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mualem_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m

  relative_permeability = 0.d0
  dkr_Se = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      return
    endif
  endif
  
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**this%m)**2.d0
  dkr_Se = 0.5d0*relative_permeability/Se+ &
            2.d0*Se**(one_over_m-0.5d0)* &
                (1.d0-Se_one_over_m)**(this%m-1.d0)* &
                (1.d0-(1.d0-Se_one_over_m)**this%m)
  
end subroutine RPF_Mualem_RelPerm

! ************************************************************************** !

function RPF_Mualem_VG_Gas_Create()

  ! Creates the van Genutchten Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_VG_gas_type), pointer :: RPF_Mualem_VG_Gas_Create
  
  allocate(RPF_Mualem_VG_Gas_Create)
  call RPF_Mualem_VG_Gas_Create%Init()
  
end function RPF_Mualem_VG_Gas_Create

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_Init(this)

  ! Initializes the van Genutchten Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_VG_gas_type) :: this

  call RPF_Mualem_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_VG_Gas_Init

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rpf_Mualem_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_VG_GAS'
  endif  
  call RPF_Mualem_Verify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif 
  
end subroutine RPF_Mualem_VG_Gas_Verify

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_Se,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_VG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_Se = UNINITIALIZED_DOUBLE
  if (Se >= 1.d0) then
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = sqrt(Seg)*(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)
  
end subroutine RPF_Mualem_VG_Gas_RelPerm

! ************************************************************************** !

function SF_BC_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BC_type), pointer :: SF_BC_Create
  
  allocate(SF_BC_Create)
  call SF_BC_Create%Init()
  
end function SF_BC_Create

! ************************************************************************** !

subroutine SF_BC_Init(this)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this
  character(len=MAXWORDLENGTH) :: name
  type(option_type) :: option

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
end subroutine SF_BC_Init

! ************************************************************************** !

subroutine SF_BC_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BROOKS_COREY'
  endif  
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif 
  
end subroutine SF_BC_Verify

! ************************************************************************** !

subroutine SF_BC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Brooks-Corey saturation function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)

  ! polynomial fitting pc as a function of saturation
  ! 1.05 is essentially pc*alpha (i.e. pc = 1.05/alpha)
  this%sat_poly => PolynomialCreate()
  this%sat_poly%low = 1.05d0**(-this%lambda)
  this%sat_poly%high = 1.d0
  
  b = 0.d0
  ! fill right hand side
  ! capillary pressure at 1
  b(1) = 1.05d0/this%alpha 
  ! capillary pressure at 2
  b(2) = 0.d0
  ! derivative of pressure at saturation_1
  ! pc = Se**(-1/lambda)/alpha
  ! dpc_dSe = -1/lambda*Se**(-1/lambda-1)/alpha
  b(3) = -1.d0/this%lambda* &
          this%sat_poly%low**(-1.d0/this%lambda-1.d0)/ &
          this%alpha

  call QuadraticPolynomialSetup(this%sat_poly%low,this%sat_poly%high,b(1:3), &
                                ! indicates derivative given at 1
                                PETSC_TRUE) 
      
  this%sat_poly%coefficients(1:3) = b(1:3)

  ! polynomial fitting saturation as a function of pc
  !geh: cannot invert the pressure/saturation relationship above
  !     since it can result in saturations > 1 with both
  !     quadratic and cubic polynomials
  ! fill matix with values
  this%pres_poly => PolynomialCreate()
  this%pres_poly%low = 0.95/this%alpha
  this%pres_poly%high = 1.05/this%alpha
  
  b = 0.d0
  ! Se at 1
  b(1) = 1.d0
  ! Se at 2
  b(2) = (this%pres_poly%high*this%alpha)** &
          (-this%lambda)
  ! derivative of Se at 1
  b(3) = 0.d0 
  ! derivative of Se at 2
  b(4) = -this%lambda/this%pres_poly%high* &
            (this%pres_poly%high*this%alpha)** &
              (-this%lambda)

  call CubicPolynomialSetup(this%pres_poly%low,this%pres_poly%high,b)

  this%pres_poly%coefficients(1:4) = b(1:4)
  
  
end subroutine SF_BC_SetupPolynomials

! ************************************************************************** !

subroutine SF_BC_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dummy_real
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dummy_real)
      return
    endif
  endif
  capillary_pressure = (Se**(-1.d0/this%lambda))/this%alpha

  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_BC_CapillaryPressure

! ************************************************************************** !

subroutine SF_BC_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_pres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: Se
  PetscReal :: dSe_pc
  PetscReal :: dsat_pc
  
  dsat_pres = 0.d0
  
  ! reference #1
  if (capillary_pressure < this%pres_poly%low) then
    liquid_saturation = 1.d0
    return
  else if (capillary_pressure < this%pres_poly%high) then
    call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                 capillary_pressure,Se,dSe_pc)
  else
    pc_alpha_neg_lambda = (capillary_pressure*this%alpha)**(-this%lambda)
    Se = pc_alpha_neg_lambda
    dSe_pc = -this%lambda/capillary_pressure*pc_alpha_neg_lambda
  endif
  liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
  dsat_pc = (1.d0-this%Sr)*dSe_pc
  
end subroutine SF_BC_Saturation

! ************************************************************************** !

function RPF_Burdine_Create()

  ! Creates the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_type), pointer :: RPF_Burdine_Create
  
  allocate(RPF_Burdine_Create)
  call RPF_Burdine_Create%Init()
  
end function RPF_Burdine_Create

! ************************************************************************** !

subroutine RPF_Burdine_Init(this)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_Init

! ************************************************************************** !

subroutine RPF_Burdine_Verify(this,name,option)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  use Option_module
  
  implicit none
  
  class(rpf_Burdine_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE'
  endif    
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif
  
end subroutine RPF_Burdine_Verify

! ************************************************************************** !

subroutine RPF_Burdine_RelPerm(this,liquid_saturation,relative_permeability, &
                              dkr_Se,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Burdine_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power

  relative_permeability = 0.d0
  dkr_Se = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    return
  endif
  
  ! reference #1
  power = 3.d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se          
  
end subroutine RPF_Burdine_RelPerm

! ************************************************************************** !

function RPF_Burdine_BC_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_gas_type), pointer :: RPF_Burdine_BC_Gas_Create
  
  allocate(RPF_Burdine_BC_Gas_Create)
  call RPF_Burdine_BC_Gas_Create%Init()
  
end function RPF_Burdine_BC_Gas_Create

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this

  call RPF_Burdine_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_BC_Gas_Init

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_BC_GAS'
  endif    
  call RPF_Burdine_Verify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_Burdine_BC_Gas_Verify

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_Se,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Burdine_BC_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_Se = UNINITIALIZED_DOUBLE
  if (Se >= 1.d0) then
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
          ! reference #1
  relative_permeability = Seg*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
  
end subroutine RPF_Burdine_BC_Gas_RelPerm

! ************************************************************************** !

function RPF_TOUGH2_IRP7_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type), pointer :: RPF_TOUGH2_IRP7_Gas_Create
  
  allocate(RPF_TOUGH2_IRP7_Gas_Create)
  call RPF_TOUGH2_IRP7_Gas_Create%Init()
  
end function RPF_TOUGH2_IRP7_Gas_Create

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this

  call RPF_Mualem_Init(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_TOUGH2_IRP7_Gas_Init

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,TOUGH2_IRP7_GAS'
  endif    
  call RPF_Mualem_Verify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_TOUGH2_IRP7_Gas_Verify

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_RelPerm(this,liquid_saturation, &
                                       relative_permeability,dkr_Se,option)
  ! 
  ! TOUGH2 IRP(7) equations from Appendix G of TOUGH2 user manual
  !
  use Option_module
  
  implicit none

  class(rpf_TOUGH2_IRP7_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_Se
  type(option_type), intent(inout) :: option
  
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_Se
  PetscReal :: Se
  PetscReal :: Seg

  relative_permeability = 0.d0
  dkr_Se = UNINITIALIZED_DOUBLE
  
                 ! essentially zero
  if (this%Srg < 1.d-40) then
    call RPF_Mualem_RelPerm(this,liquid_saturation, &
                            liquid_relative_permeability, &
                            liquid_dkr_Se,option)
    relative_permeability = 1.d0 - liquid_relative_permeability
    return
  endif  
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  if (Se >= 1.d0) then
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = Seg*Seg*(1.d0-Se*Se)
  
end subroutine RPF_TOUGH2_IRP7_Gas_RelPerm

! ************************************************************************** !

subroutine PolynomialDestroy(poly)
  ! 
  ! Destroys a polynomial smoother
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  type(polynomial_type), pointer :: poly
  
  if (.not.associated(poly)) return
  
  deallocate(poly)
  nullify(poly)

end subroutine PolynomialDestroy

! ************************************************************************** !

subroutine SaturationFunctionDestroy(sf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(sat_func_base_type), pointer :: sf
  
  if (.not.associated(sf)) return
  
  call PolynomialDestroy(sf%sat_poly)
  call PolynomialDestroy(sf%sat_poly)
  deallocate(sf)
  nullify(sf)

end subroutine SaturationFunctionDestroy

! ************************************************************************** !

subroutine PermeabilityFunctionDestroy(rpf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(rel_perm_func_base_type), pointer :: rpf
  
  if (.not.associated(rpf)) return
  
  call PolynomialDestroy(rpf%poly)
  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionDestroy

! ************************************************************************** !

recursive subroutine CharacteristicCurvesDestroy(cc)
  ! 
  ! Destroys a characteristic curve
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: cc
  
  if (.not.associated(cc)) return
  
  call CharacteristicCurvesDestroy(cc%next)
  
  call SaturationFunctionDestroy(cc%saturation_function)
  call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
  call PermeabilityFunctionDestroy(cc%gas_rel_perm_function)
  
  deallocate(cc)
  nullify(cc)
  
end subroutine CharacteristicCurvesDestroy

end module Characteristic_Curves_module
