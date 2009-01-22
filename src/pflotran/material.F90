module Material_module
 

  implicit none

  private

#include "definitions.h"
 
  type, public :: material_property_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: permeability(3,3)
    PetscReal :: permeability_pwr
    character(len=MAXWORDLENGTH) :: permeability_filename
    PetscReal :: porosity
    PetscReal :: tortuosity
    PetscInt :: saturation_function_id
    character(len=MAXWORDLENGTH) :: saturation_function_name
    PetscReal :: rock_density
    PetscReal :: specific_heat
    PetscReal :: thermal_conductivity_dry
    PetscReal :: thermal_conductivity_wet
    PetscReal :: pore_compressibility
    PetscReal :: pore_expansivity    
    type(material_property_type), pointer :: next
  end type material_property_type
  
  type, public :: material_property_ptr_type
    type(material_property_type), pointer :: ptr
  end type material_property_ptr_type
  
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
  
  public :: MaterialPropertyCreate, &
            MaterialPropertyDestroy, &
            MaterialPropertyAddToList, &
            MaterialPropGetPtrFromList, &
            MaterialPropGetPtrFromArray, &
            MaterialPropConvertListToArray, &
            MaterialPropertyRead, &
            SaturationFunctionCreate, &
            SaturationFunctionDestroy, &
            SaturationFunctionAddToList, &
            SaturationFunctionCompute, &
            SaturatFuncConvertListToArray, &
            SaturationFunctionComputeSpline

  PetscInt, parameter :: VAN_GENUCHTEN = 1
  PetscInt, parameter :: BROOKS_COREY = 2
  PetscInt, parameter :: THOMEER_COREY = 3

  PetscInt, parameter :: DEFAULT = 0
  PetscInt, parameter :: BURDINE = 1
  PetscInt, parameter :: MUALEM = 2
  
contains

! ************************************************************************** !
!
! MaterialPropertyCreate: Creates a material property
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialPropertyCreate()
  
  implicit none

  type(material_property_type), pointer :: MaterialPropertyCreate
  
  type(material_property_type), pointer :: material_property
  
  allocate(material_property)
  material_property%id = 0
  material_property%name = ''
  material_property%permeability = 0.d0
  material_property%permeability_pwr = 0.d0
  material_property%permeability_filename = ''
  material_property%porosity = 0.d0
  material_property%tortuosity = 0.d0
  material_property%saturation_function_id = 0
  material_property%saturation_function_name = ''
  material_property%rock_density = 0.d0
  material_property%specific_heat = 0.d0
  material_property%thermal_conductivity_dry = 0.d0
  material_property%thermal_conductivity_wet = 0.d0
  material_property%pore_compressibility = 0.d0
  material_property%pore_expansivity = 0.d0  
  nullify(material_property%next)
  MaterialPropertyCreate => material_property

end function MaterialPropertyCreate

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
! MaterialPropertyRead: Reads in contents of a material_property card
! author: Glenn Hammond
! date: 01/13/09
! 
! ************************************************************************** !
subroutine MaterialPropertyRead(material_property,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(material_property_type) :: material_property
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','MATERIAL_PROPERTY')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('NAME') 
        call InputReadWord(input,option,material_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'material_property name','MATERIAL_PROPERTY')
      case('ID') 
        call InputReadInt(input,option,material_property%id)
        call InputErrorMsg(input,option,'material_property id','MATERIAL_PROPERTY')
      case('SATURATION_FUNCTION') 
        call InputReadWord(input,option,material_property%saturation_function_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'material_property saturation function name','MATERIAL_PROPERTY')
      case('POROSITY')
        call InputReadDouble(input,option,material_property%porosity)
        call InputErrorMsg(input,option,'porosity','MATERIAL_PROPERTY')
      case('TORTUOSITY')
        call InputReadDouble(input,option,material_property%tortuosity)
        call InputErrorMsg(input,option,'tortuosity','MATERIAL_PROPERTY')
      case('PERMEABILITY')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'MATERIAL_PROPERTY,PERMEABILITY')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','MATERIAL_PROPERTY,PERMEABILITY')   
          select case(trim(word))
            case('PERM_X')
              call InputReadDouble(input,option,material_property%permeability(1,1))
              call InputErrorMsg(input,option,'x permeability','MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Y')
              call InputReadDouble(input,option,material_property%permeability(2,2))
              call InputErrorMsg(input,option,'y permeability','MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Z')
              call InputReadDouble(input,option,material_property%permeability(3,3))
              call InputErrorMsg(input,option,'z permeability','MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_ISO')
              call InputReadDouble(input,option,material_property%permeability(1,1))
              call InputErrorMsg(input,option,'isotropic permeability','MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(2,2) = material_property%permeability(1,1)
              material_property%permeability(3,3) = material_property%permeability(1,1)
            case('PERM_POWER')
              call InputReadDouble(input,option,material_property%permeability_pwr)
              call InputErrorMsg(input,option,'permeability power','MATERIAL_PROPERTY,PERMEABILITY')
            case('RANDOM_DATASET')
              call InputReadWord(input,option,material_property%permeability_filename,PETSC_TRUE)
              call InputErrorMsg(input,option,'RANDOM_DATASET,FILENAME','MATERIAL_PROPERTY,PERMEABILITY')   
!            case('ISOTROPIC')
            case default
              option%io_buffer = 'keyword not recognized in material_property,permeability'
              call printErrMsg(option)
          end select
        enddo

      case default
        option%io_buffer = 'Keyword: ' // keyword // &
                           ' not recognized in material_property'    
        call printErrMsg(option)
    end select 
  
  enddo  


end subroutine MaterialPropertyRead

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
! MaterialPropertyAddToList: Adds a material property to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine MaterialPropertyAddToList(material_property,list)

  implicit none
  
  type(material_property_type), pointer :: material_property
  type(material_property_type), pointer :: list

  type(material_property_type), pointer :: cur_material_property
  
  if (associated(list)) then
    cur_material_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_material_property%next)) exit
      cur_material_property => cur_material_property%next
    enddo
    cur_material_property%next => material_property
  else
    list => material_property
  endif
  
end subroutine MaterialPropertyAddToList

! ************************************************************************** !
!
! MaterialPropConvertListToArray: Creates an array of pointers to the 
!                                material_properties in the list
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine MaterialPropConvertListToArray(list,array)

  implicit none
  
  type(material_property_type), pointer :: list
  type(material_property_ptr_type), pointer :: array(:)
    
  type(material_property_type), pointer :: cur_material_property
  PetscInt :: max_id

  max_id = 0
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    max_id = max(max_id,cur_material_property%id)
    cur_material_property => cur_material_property%next
  enddo
  
  allocate(array(max_id))
  
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    array(cur_material_property%id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo

end subroutine MaterialPropConvertListToArray

! ************************************************************************** !
!
! SaturationFunctionAddToList: Adds a saturation function to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine SaturationFunctionAddToList(saturation_function,list)

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
! MaterialPropGetPtrFromList: Returns a pointer to the material property
!                             matching material_name
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialPropGetPtrFromList(material_property_name, &
                                    material_property_list)

  use String_module
  
  implicit none
  
  type(material_property_type), pointer :: MaterialPropGetPtrFromList
  character(len=MAXWORDLENGTH) :: material_property_name
  type(material_property_type), pointer :: material_property_list
  PetscInt :: length
  type(material_property_type), pointer :: material_property
    
  nullify(MaterialPropGetPtrFromList)
  material_property => material_property_list
  
  do 
    if (.not.associated(material_property)) exit
    length = len_trim(material_property_name)
    if (length == len_trim(material_property%name) .and. &
        StringCompare(material_property%name,material_property_name,length)) then
      MaterialPropGetPtrFromList => material_property
      return
    endif
    material_property => material_property%next
  enddo
  
end function MaterialPropGetPtrFromList

! ************************************************************************** !
!
! MaterialPropGetPtrFromArray: Returns a pointer to the material property
!                              matching material_name
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialPropGetPtrFromArray(material_property_name, &
                                     material_property_array)

  use String_module

  implicit none
  
  type(material_property_type), pointer :: MaterialPropGetPtrFromArray
  character(len=MAXWORDLENGTH) :: material_property_name
  type(material_property_ptr_type), pointer :: material_property_array(:)
  PetscInt :: length
  PetscInt :: imaterial_property
    
  nullify(MaterialPropGetPtrFromArray)
  
  do imaterial_property = 1, size(material_property_array)
    length = len_trim(material_property_name)
    if (length == len_trim(material_property_array(imaterial_property)%ptr%name) .and. &
        StringCompare(material_property_array(imaterial_property)%ptr%name, &
                        material_property_name,length)) then
      MaterialPropGetPtrFromArray => material_property_array(imaterial_property)%ptr
      return
    endif
  enddo
  
end function MaterialPropGetPtrFromArray

! ************************************************************************** !
!
! MaterialPropertyDestroy: Destroys a material_property
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine MaterialPropertyDestroy(material_property)

  implicit none
  
  type(material_property_type), pointer :: material_property
  
  if (.not.associated(material_property)) return
  
  call MaterialPropertyDestroy(material_property%next)
  
  deallocate(material_property)
  nullify(material_property)
  
end subroutine MaterialPropertyDestroy

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
