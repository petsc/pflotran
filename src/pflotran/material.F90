module Material_module
 
  use Region_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: material_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: permeability(3,3)
    PetscReal :: permeability_pwr
    character(len=MAXWORDLENGTH) :: permeability_filename
    PetscReal :: porosity
    PetscReal :: tortuosity
    PetscInt :: ithrm
    PetscInt :: icap
    type(material_type), pointer :: next
  end type material_type
  
  type, public :: fluid_property_type
    PetscReal, pointer :: diff_base(:)
    PetscReal, pointer :: diff_exp(:)
  end type fluid_property_type
  
  type, public :: material_ptr_type
    type(material_type), pointer :: ptr
  end type material_ptr_type
  
  type, public :: thermal_property_type
    PetscInt :: id
    PetscReal :: rock_density
    PetscReal :: spec_heat
    PetscReal :: therm_cond_dry
    PetscReal :: therm_cond_wet
    PetscReal :: pore_compress
    PetscReal :: pore_expansivity
    PetscReal :: tort_bin_diff
    PetscReal :: vap_air_diff_coef
    PetscReal :: exp_binary_diff
    PetscReal :: enh_binary_diff_coef
    type(thermal_property_type), pointer :: next
  end type thermal_property_type
  
  type, public :: saturation_function_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: saturation_function_ctype
    PetscInt :: saturation_function_itype
    character(len=MAXWORDLENGTH) :: permeability_function_ctype
    PetscInt :: permeability_function_itype
    PetscReal, pointer :: Sr(:)
    PetscReal :: m
    PetscReal :: lambda
    PetscReal :: alpha
    PetscReal :: pcwmax
    PetscReal :: betac
    PetscReal :: power
    PetscInt :: hysteresis_id
    PetscInt :: hysteresis_params(6)
    PetscReal :: BC_pressure_low
    PetscReal :: BC_pressure_high
    PetscReal :: BC_spline_coefficients(4)
    type(saturation_function_type), pointer :: next
  end type saturation_function_type
  
  type, public :: saturation_function_ptr_type
    type(saturation_function_type), pointer :: ptr
  end type saturation_function_ptr_type
  
  public :: MaterialCreate, ThermalPropertyCreate, SaturationFunctionCreate, &
            MaterialDestroy, ThermalPropertyDestroy, &
            SaturationFunctionDestroy, &
            MaterialAddToList, ThermalAddPropertyToList, &
            SaturationFunctionAddToList, &
            MaterialGetPtrFromList, &
            MaterialGetPtrFromArray, &
            SaturationFunctionCompute, &
            SaturatFuncConvertListToArray, &
            MaterialConvertListToArray, &
            SaturationFunctionComputeSpline, &
            FluidPropertyCreate, &
            FluidPropertyDestroy

  PetscInt, parameter :: VAN_GENUCHTEN = 1
  PetscInt, parameter :: BROOKS_COREY = 2
  PetscInt, parameter :: THOMEER_COREY = 3

  PetscInt, parameter :: DEFAULT = 0
  PetscInt, parameter :: BURDINE = 1
  PetscInt, parameter :: MUALEM = 2
  
contains

! ************************************************************************** !
!
! MaterialCreate: Creates a material
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialCreate()
  
  implicit none

  type(material_type), pointer :: MaterialCreate
  
  type(material_type), pointer :: material
  
  allocate(material)
  material%id = 0
  material%name = ''
  material%permeability = 0.d0
  material%permeability_pwr = 0.d0
  material%permeability_filename = ''
  material%porosity = 0.d0
  material%tortuosity = 0.d0
  material%ithrm = 0
  material%icap = 0
  nullify(material%next)
  MaterialCreate => material

end function MaterialCreate

! ************************************************************************** !
!
! MaterialCreate: Creates a thermal property
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function ThermalPropertyCreate()
  
  implicit none

  type(thermal_property_type), pointer :: ThermalPropertyCreate
  
  type(thermal_property_type), pointer :: thermal_property
  
  allocate(thermal_property)
  thermal_property%id = 0
  thermal_property%rock_density = 0.d0
  thermal_property%spec_heat = 0.d0
  thermal_property%therm_cond_dry = 0.d0
  thermal_property%therm_cond_wet = 0.d0
  thermal_property%pore_compress = 0.d0
  thermal_property%pore_expansivity = 0.d0
  thermal_property%tort_bin_diff = 0.d0
  thermal_property%vap_air_diff_coef = 0.d0
  thermal_property%exp_binary_diff = 0.d0
  thermal_property%enh_binary_diff_coef = 0.d0
  nullify(thermal_property%next)
  ThermalPropertyCreate => thermal_property

end function ThermalPropertyCreate

! ************************************************************************** !
!
! FluidPropertyCreate: Creates a fluid property object
! author: Glenn Hammond
! date: 05/07/08
!
! ************************************************************************** !
function FluidPropertyCreate(nphase)
  
  implicit none

  type(fluid_property_type), pointer :: FluidPropertyCreate
  integer :: nphase
  
  type(fluid_property_type), pointer :: fluid_property
  
  allocate(fluid_property)
  allocate(fluid_property%diff_base(nphase))
  allocate(fluid_property%diff_exp(nphase))
  fluid_property%diff_base = 0.d0
  fluid_property%diff_exp = 0.d0
  FluidPropertyCreate => fluid_property

end function FluidPropertyCreate

! ************************************************************************** !
!
! SaturationFunctionCreate: Creates a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function SaturationFunctionCreate(option)
  
  use Option_module
  
  implicit none

  type(saturation_function_type), pointer :: SaturationFunctionCreate
  type(option_type) :: option
  
  type(saturation_function_type), pointer :: saturation_function
  
  allocate(saturation_function)
  saturation_function%id = 0
  saturation_function%saturation_function_ctype = ""
  saturation_function%saturation_function_itype = VAN_GENUCHTEN
  saturation_function%permeability_function_ctype = ""
  saturation_function%permeability_function_itype = BURDINE
  allocate(saturation_function%Sr(option%nphase))
  saturation_function%Sr = 0.d0
  saturation_function%m = 0.d0
  saturation_function%lambda = 0.d0
  saturation_function%alpha = 0.d0
  saturation_function%pcwmax = 0.d0
  saturation_function%betac = 0.d0
  saturation_function%power = 0.d0
  saturation_function%hysteresis_id = 0
  saturation_function%hysteresis_params = 0
  saturation_function%BC_pressure_low = 0.d0
  saturation_function%BC_pressure_high = 0.d0
  saturation_function%BC_spline_coefficients = 0.d0
  nullify(saturation_function%next)
  SaturationFunctionCreate => saturation_function

end function SaturationFunctionCreate

! ************************************************************************** !
!
! MaterialRead: Reads in contents of a material card
! author: Glenn Hammond
! date: 01/13/09
! 
! ************************************************************************** !
subroutine MaterialRead(material,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(material_type) :: material
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','MATERIAL')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('ID') 
        call InputReadInt(input,option,material%id)
        call InputDefaultMsg(input,option,'material id')
      case('SATURATION_FUNCTION') 
        call InputReadInt(input,option,material%icap)
        call InputDefaultMsg(input,option,'material saturation function id')
      case('THERMAL_PROPERTY')
        call InputReadInt(input,option,material%ithrm)
        call InputDefaultMsg(input,option,'material thermal property id')
      case('POROSITY')
        call InputReadDouble(input,option,material%porosity)
        call InputDefaultMsg(input,option,'porosity')
      case('TORTUOSITY')
        call InputReadDouble(input,option,material%tortuosity)
        call InputDefaultMsg(input,option,'tortuosity')
      case('PERMEABILITY')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'MATERIAL,PERMEABILITY')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','MATERIAL,PERMEABILITY')   
          select case(trim(word))
            case('PERM_X')
              call InputReadDouble(input,option,material%permeability(1,1))
              call InputDefaultMsg(input,option,'x permeability')
            case('PERM_Y')
              call InputReadDouble(input,option,material%permeability(2,2))
              call InputDefaultMsg(input,option,'y permeability')
            case('PERM_Z')
              call InputReadDouble(input,option,material%permeability(3,3))
              call InputDefaultMsg(input,option,'z permeability')
            case('PERM_POWER')
              call InputReadDouble(input,option,material%permeability_pwr)
              call InputDefaultMsg(input,option,'permeability power')
            case('RANDOM_DATASET')
              call InputReadWord(input,option,material%permeability_filename,PETSC_TRUE)
              call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
!            case('ISOTROPIC')
            case default
              option%io_buffer = 'keyword not recognized in material,permeability'
              call printErrMsg(option)
          end select
        enddo

      case default
        option%io_buffer = 'Keyword: ' // keyword // &
                           ' not recognized in material'    
        call printErrMsg(option)
    end select 
  
  enddo  


end subroutine 

! ************************************************************************** !
!
! SaturationFunctionComputeSpline: Computes a spline spanning the 
!                                  discontinuity in Brooks Corey
! author: Glenn Hammond
! date: 02/27/08
!
! ************************************************************************** !
subroutine SaturationFunctionComputeSpline(option,saturation_function)
  
  use Option_module
  use Utility_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  
  PetscReal :: A(4,4), x(4), b(4)
  PetscInt :: indx(4)
  PetscInt :: d
  PetscReal :: pressure_high, pressure_low
  
  PetscReal :: alpha

  if (saturation_function%saturation_function_itype /= BROOKS_COREY) &
    return
  
  alpha = saturation_function%alpha
  
  ! fill matix with values
  pressure_high = 1.d0/alpha*2.d0
  pressure_low = 1.d0/alpha*0.5d0
  
  saturation_function%BC_pressure_low = pressure_low
  saturation_function%BC_pressure_high = pressure_high
  
  
  A(1,1) = 1.d0
  A(2,1) = 1.d0
  A(3,1) = 0.d0
  A(4,1) = 0.d0
  
  A(1,2) = pressure_high
  A(2,2) = pressure_low
  A(3,2) = 1.d0
  A(4,2) = 1.d0
  
  A(1,3) = pressure_high**2.d0
  A(2,3) = pressure_low**2.d0
  A(3,3) = 2.d0*pressure_high
  A(4,3) = 2.d0*pressure_low
  
  A(1,4) = pressure_high**3.d0
  A(2,4) = pressure_low**3.d0
  A(3,4) = 3.d0*pressure_high**2.d0
  A(4,4) = 3.d0*pressure_low**2.d0
  
  b(1) = (pressure_high*alpha)**(-saturation_function%lambda)
  b(2) = 1.d0
  b(3) = -saturation_function%lambda/pressure_high* &
        (pressure_high*alpha)**(-saturation_function%lambda)
  b(4) = 0.d0
  
  call ludcmp(A,4,indx,d)
  call lubksb(A,4,indx,b)
  
  saturation_function%BC_spline_coefficients(1:4) = b(1:4)
  
end subroutine SaturationFunctionComputeSpline

! ************************************************************************** !
!
! MaterialAddToList: Adds a thermal property to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine MaterialAddToList(material,list)

  implicit none
  
  type(material_type), pointer :: material
  type(material_type), pointer :: list

  type(material_type), pointer :: cur_material
  
  if (associated(list)) then
    cur_material => list
    ! loop to end of list
    do
      if (.not.associated(cur_material%next)) exit
      cur_material => cur_material%next
    enddo
    cur_material%next => material
  else
    list => material
  endif
  
end subroutine MaterialAddToList

! ************************************************************************** !
!
! MaterialConvertListToArray: Creates an array of pointers to the 
!                                materials in the list
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine MaterialConvertListToArray(list,array)

  implicit none
  
  type(material_type), pointer :: list
  type(material_ptr_type), pointer :: array(:)
    
  type(material_type), pointer :: cur_material
  PetscInt :: max_id

  max_id = 0
  cur_material => list
  do 
    if (.not.associated(cur_material)) exit
    max_id = max(max_id,cur_material%id)
    cur_material => cur_material%next
  enddo
  
  allocate(array(max_id))
  
  cur_material => list
  do 
    if (.not.associated(cur_material)) exit
    array(cur_material%id)%ptr => cur_material
    cur_material => cur_material%next
  enddo

end subroutine MaterialConvertListToArray

! ************************************************************************** !
!
! ThermalAddPropertyToList: Adds a thermal property to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine ThermalAddPropertyToList(thermal_property,list)

  implicit none
  
  type(thermal_property_type), pointer :: thermal_property
  type(thermal_property_type), pointer :: list

  type(thermal_property_type), pointer :: cur_thermal_property
  
  if (associated(list)) then
    cur_thermal_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_thermal_property%next)) exit
      cur_thermal_property => cur_thermal_property%next
    enddo
    cur_thermal_property%next => thermal_property
  else
    list => thermal_property
  endif
  
end subroutine ThermalAddPropertyToList

! ************************************************************************** !
!
! SaturationFunctionAddToList: Adds a saturation function to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine SaturationFunctionAddToList(saturation_function,list)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  type(saturation_function_type), pointer :: list

  type(saturation_function_type), pointer :: cur_saturation_function
  
  if (associated(list)) then
    cur_saturation_function => list
    ! loop to end of list
    do
      if (.not.associated(cur_saturation_function%next)) exit
      cur_saturation_function => cur_saturation_function%next
    enddo
    cur_saturation_function%next => saturation_function
  else
    list => saturation_function
  endif
  
end subroutine SaturationFunctionAddToList

! ************************************************************************** !
!
! SaturatFuncConvertListToArray: Creates an array of pointers to the 
!                                saturation functions in the list
! author: Glenn Hammond
! date: 12/11/07
!
! ************************************************************************** !
subroutine SaturatFuncConvertListToArray(list,array)

  implicit none
  
  type(saturation_function_type), pointer :: list
  type(saturation_function_ptr_type), pointer :: array(:)
    
  type(saturation_function_type), pointer :: cur_saturation_function
  PetscInt :: count

  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    cur_saturation_function => cur_saturation_function%next
  enddo
  
  allocate(array(count))
  
  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    array(count)%ptr => cur_saturation_function
    cur_saturation_function => cur_saturation_function%next
  enddo

end subroutine SaturatFuncConvertListToArray

! ************************************************************************** !
!
! SaturationFunctionCompute: Computes the saturation and relative permeability
!                            (and associated derivatives) as a function of 
!                            capillary pressure
! author: Glenn Hammond
! date: 12/11/07
!
! ************************************************************************** !
subroutine SaturationFunctionCompute(pressure,saturation,relative_perm, &
                                     dsat_pres,dkr_pres, &
                                     saturation_function, &
                                     auxvar1,auxvar2,option)

  use Option_module
  
  implicit none

  PetscReal :: pressure, saturation, relative_perm, dsat_pres, dkr_pres
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1,auxvar2
  type(option_type) :: option

  PetscInt :: iphase
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_pc, dsat_pc, dkr_pc
  PetscReal :: dkr_Se, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: por, perm
  PetscReal :: Fg, a, Pd, PHg
  
  iphase = 1
  dsat_pres = 0.d0
  dkr_pres = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      if (pressure >= option%reference_pressure) then
        saturation = 1.d0
        relative_perm = 1.d0
        return
      else
        alpha = saturation_function%alpha
        pc = option%reference_pressure-pressure
        m = saturation_function%m
        n = 1.d0/(1.d0-m)
        pc_alpha = pc*alpha
        pc_alpha_n = pc_alpha**n
        one_plus_pc_alpha_n = 1.d0+pc_alpha_n
        Sr = saturation_function%Sr(iphase)
        Se = one_plus_pc_alpha_n**(-m)
        dSe_pc = -m*n*alpha*pc_alpha_n/(pc_alpha*one_plus_pc_alpha_n**(m+1))
        saturation = Sr + (1.d0-Sr)*Se
        dsat_pc = (1.d0-Sr)*dSe_pc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'vG Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'vG Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = Se*Se*(1.d0-(1.d0-Se_one_over_m)**m)
          dkr_Se = 2.d0*relative_perm/Se + &
                   Se*Se_one_over_m*(1.d0-Se_one_over_m)**(m-1.d0)
          dkr_pc = dkr_Se*dSe_pc
        case(MUALEM)
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
          dkr_Se = 0.5d0*relative_perm/Se+ &
                   2.d0*Se**(one_over_m-0.5d0)* &
                        (1.d0-Se_one_over_m)**(m-1.d0)* &
                        (1.d0-(1.d0-Se_one_over_m)**m)
          dkr_pc = dkr_Se*dSe_pc
        case default
          option%io_buffer = 'Unknown relative permeabilty function' 
          call printErrMsg(option)
      end select
    case(BROOKS_COREY)
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      pc = option%reference_pressure-pressure
#if 0      
      if (pc < saturation_function%BC_pressure_low) then
        saturation = 1.d0
        relative_perm = 1.d0
        return
      else if (pc < saturation_function%BC_pressure_high) then
        Sr = saturation_function%Sr(iphase)
        Se = saturation_function%BC_spline_coefficients(1)+ &
             saturation_function%BC_spline_coefficients(2)*pc+ &
             saturation_function%BC_spline_coefficients(3)*pc**2.d0+ &
             saturation_function%BC_spline_coefficients(4)*pc**3.d0
        dSe_pc = saturation_function%BC_spline_coefficients(2)+ &
                  saturation_function%BC_spline_coefficients(3)*2.d0*pc+ &
                  saturation_function%BC_spline_coefficients(4)*3.d0*pc**2.d0
        saturation = Sr + (1.d0-Sr)*Se
        dsat_pc = (1.d0-Sr)*dSe_pc
#else
      if (pc < one_over_alpha) then
        saturation = 1.d0
        relative_perm = 1.d0
        return
#endif        
      else
        lambda = saturation_function%lambda
        Sr = saturation_function%Sr(iphase)
        pc_alpha_neg_lambda = (pc*alpha)**(-lambda)
        Se = pc_alpha_neg_lambda
        dSe_pc = -lambda/pc*pc_alpha_neg_lambda
        saturation = Sr + (1.d0-Sr)*Se
!        dsat_pc = -lambda*(1.d0-Sr)/pc*pc_alpha_neg_lambda
        dsat_pc = (1.d0-Sr)*dSe_pc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'BC Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'BC Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          power = 3.d0+2.d0/lambda
          relative_perm = Se**power
          dkr_Se = power*relative_perm/Se
          dkr_pc = dkr_Se*dSe_pc
        case(MUALEM)
          power = 2.5d0+2.d0/lambda
          relative_perm = Se**power
          dkr_Se = power*relative_perm/Se
          dkr_pc = dkr_Se*dSe_pc
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(THOMEER_COREY)
      pc = option%reference_pressure-pressure
      por = auxvar1
      perm = auxvar2*1.013202d15 ! convert from m^2 to mD
      Fg = saturation_function%alpha
      a = saturation_function%m
      Pd = 100.d0*por/sqrt(perm/(3.8068d0*(Fg**(-1.334d0)))) ! psi
      PHg = 9.63051d-4*pc
      if (PHg > Pd) then
        saturation = 1.d0-exp(-Fg/log10(PHg/Pd))
#if 0
        alpha = pc*(1.d0+1.d-8)
        m = 9.63051d-4*alpha
        n = 1.d0-exp(-Fg/log10(m/Pd))
        n = (n-saturation)/(alpha-pc)
#endif        
        dsat_pc = (saturation-1.d0)*Fg/(log10(PHg/Pd)**2.d0)/(pc*2.30258509d0)
        ! Sr assumed to be zero
        relative_perm = saturation**a
        dkr_pc = a*saturation**(a-1.d0)*dsat_pc
      else
        saturation = 1.d0
        relative_perm = 1.d0
        return
      endif
    case default
      option%io_buffer = 'Unknown saturation function'
      call printErrMsg(option)
  end select
  dsat_pres = -dsat_pc 
  dkr_pres = -dkr_pc

end subroutine SaturationFunctionCompute

! ************************************************************************** !
!
! MaterialGetPtrFromList: Returns a pointer to the material matching 
!                         material_name
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialGetPtrFromList(material_name,material_list)

  use String_module
  
  implicit none
  
  type(material_type), pointer :: MaterialGetPtrFromList
  character(len=MAXWORDLENGTH) :: material_name
  type(material_type), pointer :: material_list
  PetscInt :: length
  type(material_type), pointer :: material
    
  nullify(MaterialGetPtrFromList)
  material => material_list
  
  do 
    if (.not.associated(material)) exit
    length = len_trim(material_name)
    if (length == len_trim(material%name) .and. &
        StringCompare(material%name,material_name,length)) then
      MaterialGetPtrFromList => material
      return
    endif
    material => material%next
  enddo
  
end function MaterialGetPtrFromList

! ************************************************************************** !
!
! MaterialGetPtrFromArray: Returns a pointer to the material matching 
!                          material_name
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialGetPtrFromArray(material_name,material_array)

  use String_module

  implicit none
  
  type(material_type), pointer :: MaterialGetPtrFromArray
  character(len=MAXWORDLENGTH) :: material_name
  type(material_ptr_type), pointer :: material_array(:)
  PetscInt :: length
  PetscInt :: imaterial
    
  nullify(MaterialGetPtrFromArray)
  
  do imaterial = 1, size(material_array)
    length = len_trim(material_name)
    if (length == len_trim(material_array(imaterial)%ptr%name) .and. &
        StringCompare(material_array(imaterial)%ptr%name, &
                        material_name,length)) then
      MaterialGetPtrFromArray => material_array(imaterial)%ptr
      return
    endif
  enddo
  
end function MaterialGetPtrFromArray

! ************************************************************************** !
!
! FluidPropertyDestroy: Destroys a fluid property object
! author: Glenn Hammond
! date: 05/07/08
!
! ************************************************************************** !
subroutine FluidPropertyDestroy(fluid_property)

  implicit none
  
  type(fluid_property_type), pointer :: fluid_property
  
  if (.not.associated(fluid_property)) return
  
  if (associated(fluid_property%diff_base)) &
    deallocate(fluid_property%diff_base)
  nullify(fluid_property%diff_base)
  if (associated(fluid_property%diff_exp)) &
    deallocate(fluid_property%diff_exp)
  nullify(fluid_property%diff_exp)
  
  deallocate(fluid_property)
  nullify(fluid_property)
  
end subroutine FluidPropertyDestroy

! ************************************************************************** !
!
! MaterialDestroy: Destroys a material
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine MaterialDestroy(material)

  implicit none
  
  type(material_type), pointer :: material
  
  if (.not.associated(material)) return
  
  call MaterialDestroy(material%next)
  
  deallocate(material)
  nullify(material)
  
end subroutine MaterialDestroy

! ************************************************************************** !
!
! SaturationFunctionDestroy: Destroys a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine ThermalPropertyDestroy(thermal_property)

  implicit none
  
  type(thermal_property_type), pointer :: thermal_property
  
  if (.not.associated(thermal_property)) return
  
  call ThermalPropertyDestroy(thermal_property%next)

  deallocate(thermal_property)
  nullify(thermal_property)
  
end subroutine ThermalPropertyDestroy

! ************************************************************************** !
!
! SaturationFunctionDestroy: Destroys a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine SaturationFunctionDestroy(saturation_function)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  
  if (.not.associated(saturation_function)) return
  
  call SaturationFunctionDestroy(saturation_function%next)
    
  if (associated(saturation_function%Sr)) deallocate(saturation_function%Sr)
  nullify(saturation_function%Sr)
    
  deallocate(saturation_function)
  nullify(saturation_function)
  
end subroutine SaturationFunctionDestroy

end module Material_module
