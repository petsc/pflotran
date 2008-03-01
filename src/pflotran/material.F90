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
    PetscReal :: porosity
    PetscReal :: tortuosity
    PetscInt :: ithrm
    PetscInt :: icap
    type(material_type), pointer :: next
  end type material_type
  
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
    PetscInt :: ihist 
    PetscReal :: BC_pressure_offset
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
            SaturationFunctionCompute, &
            SaturatFuncConvertListToArray, &
            MaterialConvertListToArray, &
            SaturationFunctionComputeSpline

  PetscInt, parameter :: VAN_GENUCHTEN = 1
  PetscInt, parameter :: BROOKS_COREY = 2

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
  material%name = ""
  material%permeability = 0.d0
  material%permeability_pwr = 0.d0
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
  saturation_function%ihist = 0
  saturation_function%BC_pressure_offset = 0.d0
  saturation_function%BC_spline_coefficients = 0.d0
  nullify(saturation_function%next)
  SaturationFunctionCreate => saturation_function

end function SaturationFunctionCreate

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
  PetscReal :: pressure_a, pressure_b
  
  PetscReal :: alpha

  if (saturation_function%saturation_function_itype /= BROOKS_COREY) &
    return
  
  alpha = saturation_function%alpha
  
  ! fill matix with values
  saturation_function%BC_pressure_offset = 500.d0
  pressure_a = 1.d0/alpha+saturation_function%BC_pressure_offset
  pressure_b = 1.d0/alpha-saturation_function%BC_pressure_offset
  
  A(1,1) = 1.d0
  A(2,1) = 1.d0
  A(3,1) = 0.d0
  A(4,1) = 0.d0
  
  A(1,2) = pressure_a
  A(2,2) = pressure_b
  A(3,2) = 1.d0
  A(4,2) = 1.d0
  
  A(1,3) = pressure_a**2.d0
  A(2,3) = pressure_b**2.d0
  A(3,3) = 2.d0*pressure_a
  A(4,3) = 2.d0*pressure_b
  
  A(1,4) = pressure_a**3.d0
  A(2,4) = pressure_b**3.d0
  A(3,4) = 3.d0*pressure_a**2.d0
  A(4,4) = 3.d0*pressure_b**2.d0
  
  b(1) = (pressure_a*alpha)**(-1.d0*saturation_function%lambda)
  b(2) = 1.d0
  b(3) = -1.d0*saturation_function%lambda/pressure_a* &
        (pressure_a*alpha)**(-1.d0*saturation_function%lambda)
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
  PetscInt :: max_id = 0

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
                                     dsat_pres,dkr_pres,saturation_function, &
                                     option)

  use Option_module
  
  implicit none

  PetscReal :: pressure, saturation, relative_perm, dsat_pres, dkr_pres
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscInt :: iphase = 1
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_pc, dsat_pc, dkr_pc
  PetscReal :: dkr_Se, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  
  dsat_pres = 0.d0
  dkr_pres = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      if (pressure >= option%pref) then
        saturation = 1.d0
        relative_perm = 1.d0
        return
      else
        alpha = saturation_function%alpha
        pc = option%pref-pressure
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
          call printErrMsg(option,"Unknown relative permeabilty function")
      end select
    case(BROOKS_COREY)
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      pc = option%pref-pressure
#if 1      
      if (pc < one_over_alpha-saturation_function%BC_pressure_offset) then
        saturation = 1.d0
        relative_perm = 1.d0
        return
      else if (pc < one_over_alpha+saturation_function%BC_pressure_offset) then
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
          call printErrMsg(option,"Unknown relative permeabilty function")
      end select
    case default
      call printErrMsg(option,"Unknown saturation function")
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

  use Fileio_module

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
        fiStringCompare(material%name,material_name,length)) then
      MaterialGetPtrFromList => material
      return
    endif
    material => material%next
  enddo
  
end function MaterialGetPtrFromList

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
