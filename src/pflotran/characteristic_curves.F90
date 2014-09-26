module Characteristic_Curves_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  ! Saturation Function
  type :: sat_func_base_type
    PetscReal :: Sr
    PetscReal :: pcmax
  contains
    procedure, public :: Init => SFBaseInit
    procedure, public :: CapillaryPressure => SF_Base_CapillaryPressure
    procedure, public :: Saturation => SF_Base_Saturation
  end type sat_func_base_type
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_VG_Init
    procedure, public :: CapillaryPressure => SF_VG_CapillaryPressure
    procedure, public :: Saturation => SF_VG_Saturation
  end type sat_func_VG_type  
  
  ! Relative Permeability Function
  type :: rel_perm_func_base_type
    PetscReal :: Sr
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_Mualem_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Mualem_Init
    procedure, public :: RelativePermeability => RPF_Mualem_RelPerm
  end type rel_perm_func_Mualem_type
  type, public, extends(rel_perm_func_Mualem_type) :: rpf_Mualem_VG_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_VG_Gas_Init
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Gas_RelPerm
  end type rpf_Mualem_VG_gas_type

  type, public :: characteristic_curves_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
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
            print *, 'Brooks-Corey not implemented'
            stop
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
            print *, 'Burdine not implemented'
            stop
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
      case('VERIFY') 
        this%print_me = PETSC_TRUE
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in charateristic_curves.'    
        call printErrMsg(option)
    end select 
  enddo 

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

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION'
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,saturation_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcmax)
        call InputErrorMsg(input,option,'maximum capillary pressure', &
                            error_string)
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
      class default
        option%io_buffer = 'Read routine not implemented for saturation ' // &
                           'function class.'
        call printErrMsg(option)
    end select
  enddo 

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

  input%ierr = 0
  new_phase_keyword = 'NONE'
  error_string = 'CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION'
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('RESIDUAL_SATURATION','LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,permeability_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case default
        found = PETSC_FALSE
    end select
    
    if (found) cycle
    
    ! we assume liquid phase if PHASE keyword is not present.
    select type(rpf => permeability_function)
      class is(rel_perm_func_Mualem_type)
        select case(keyword)
          case('PHASE')
            call InputReadWord(input,option,new_phase_keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'phase',error_string)
            call StringToUpper(phase_keyword) 
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            option%io_buffer = 'Keyword: ' // trim(keyword) // &
              ' not recognized in Mualem relative permeability function.'
            call printErrMsg(option)
        end select
      class is(rpf_Mualem_VG_gas_type)
        select case(keyword)
          case('PHASE')
            call InputReadWord(input,option,new_phase_keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'phase',error_string)
            call StringToUpper(phase_keyword) 
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
    if (cur_characteristic_curves%print_me .and. &
        option%myrank == option%io_rank) then
      print *, 'CharacteristicCurvesVerify not implemented'
      stop
!      call SaturationFunctionVerify(cur_saturation_function,option)
    endif
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

end subroutine CharCurvesConvertListToArray


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

!!!BASE ROUTINES

! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%Sr = 0.d0
  this%pcmax = 1.d9
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%Sr = 0.d0
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine SF_Base_CapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,option)
  
  use Option_module
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_Base_CapillaryPressure must be extended.'
  call printErrMsg(option)
  
end subroutine SF_Base_CapillaryPressure

! ************************************************************************** !

subroutine SF_Base_Saturation(this,capillary_pressure,liquid_saturation, &
                              dsat_pres,option)

  use Option_module
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_pres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SF_Base_Saturation must be extended.'
  call printErrMsg(option)
  
end subroutine SF_Base_Saturation

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
  this%alpha = 0.d0
  this%m = 0.d0
  
end subroutine SF_VG_Init

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
  
#if 0
  if (capillary_pressure < this%spline%low) then
    liquid_saturation = 1.d0
    relative_perm = 1.d0
    switch_to_saturated = PETSC_TRUE
    return
  else if (capillary_pressure < this%spline%high) then
    call CubicPolynomialEvaluate(this%spline%coefficients, &
                                 capillary_pressure,Se,dSe_pc)
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_pc = (1.d0-this%Sr)*dSe_pc
#else
  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
!    switch_to_saturated = PETSC_TRUE
    return
#endif        
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
  
  class(rel_perm_func_Mualem_type), pointer :: RPF_Mualem_Create
  
  allocate(RPF_Mualem_Create)
  call RPF_Mualem_Create%Init()
  
end function RPF_Mualem_Create

! ************************************************************************** !

subroutine RPF_Mualem_Init(this)

  ! Initializes the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rel_perm_func_Mualem_type) :: this

  call RPFBaseInit(this)
  this%m = 0.d0
  
end subroutine RPF_Mualem_Init

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
  
  implicit none

  class(rel_perm_func_Mualem_type) :: this
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
  
#ifdef MUALEM_SPLINE
  if (Se > this%spline%low) then
    call CubicPolynomialEvaluate(this%spline%coefficients, &
                                 Se,relative_permeability,dkr_Se)
  else
#endif          
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**this%m)**2.d0
  dkr_Se = 0.5d0*relative_permeability/Se+ &
            2.d0*Se**(one_over_m-0.5d0)* &
                (1.d0-Se_one_over_m)**(this%m-1.d0)* &
                (1.d0-(1.d0-Se_one_over_m)**this%m)
#ifdef MUALEM_SPLINE
  endif
#endif          
  
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
  this%Srg = 0.d0
  
end subroutine RPF_Mualem_VG_Gas_Init

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
  dkr_Se = -999.d0
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = sqrt(Seg)*(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)
  
end subroutine RPF_Mualem_VG_Gas_RelPerm

! ************************************************************************** !

recursive subroutine CharacteristicCurvesDestroy(characteristic_curves)
  ! 
  ! Destroys a characteristic curve
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: characteristic_curves
  
  if (.not.associated(characteristic_curves)) return
  
  call CharacteristicCurvesDestroy(characteristic_curves%next)
    
  deallocate(characteristic_curves)
  nullify(characteristic_curves)
  
end subroutine CharacteristicCurvesDestroy

end module Characteristic_Curves_module