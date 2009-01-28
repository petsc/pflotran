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
    PetscReal :: thermal_expansitivity    
    type(material_property_type), pointer :: next
  end type material_property_type
  
  type, public :: material_property_ptr_type
    type(material_property_type), pointer :: ptr
  end type material_property_ptr_type
  

  public :: MaterialPropertyCreate, &
            MaterialPropertyDestroy, &
            MaterialPropertyAddToList, &
            MaterialPropGetPtrFromList, &
            MaterialPropGetPtrFromArray, &
            MaterialPropConvertListToArray, &
            MaterialPropertyRead
  
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
  material_property%thermal_expansitivity = 0.d0  
  nullify(material_property%next)
  MaterialPropertyCreate => material_property

end function MaterialPropertyCreate

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
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
      case('ID') 
        call InputReadInt(input,option,material_property%id)
        call InputErrorMsg(input,option,'id','MATERIAL_PROPERTY')
      case('SATURATION_FUNCTION') 
        call InputReadWord(input,option,material_property%saturation_function_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation function name','MATERIAL_PROPERTY')
      case('ROCK_DENSITY') 
        call InputReadDouble(input,option,material_property%rock_density)
        call InputErrorMsg(input,option,'rock density','MATERIAL_PROPERTY')
      case('SPECIFIC_HEAT') 
        call InputReadDouble(input,option,material_property%specific_heat)
        call InputErrorMsg(input,option,'specific heat','MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_DRY') 
        call InputReadDouble(input,option,material_property%thermal_conductivity_dry)
        call InputErrorMsg(input,option,'dry thermal conductivity','MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_WET') 
        call InputReadDouble(input,option,material_property%thermal_conductivity_wet)
        call InputErrorMsg(input,option,'wet thermal conductivity','MATERIAL_PROPERTY')
      case('PORE_COMPRESSIBILITY') 
        call InputReadDouble(input,option,material_property%pore_compressibility)
        call InputErrorMsg(input,option,'pore compressibility','MATERIAL_PROPERTY')
      case('THERMAL_EXPANSITIVITY') 
        call InputReadDouble(input,option,material_property%thermal_expansitivity)
        call InputErrorMsg(input,option,'thermal expansitivity','MATERIAL_PROPERTY')
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
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in material_property'    
        call printErrMsg(option)
    end select 
  
  enddo  


end subroutine MaterialPropertyRead

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
  type(material_property_type), pointer :: prev_material_property
  type(material_property_type), pointer :: next_material_property
  PetscInt :: i, max_id

#if 0
! don't necessary need right now, but maybe in future
  ! reorder into ascending order
  swapped = PETSC_FALSE
  do
    if (.not.swapped) exit
    cur_material_property => list
    do 
      if (.not.associated(cur_material_property)) exit
      next_material_property => cur_material_property%next
      if (associated(next_material_property)) then
        if (cur_material_property%id > next_material_property%id) then
          ! swap
          if (associated(prev_material_property)) then
            prev_material_property%next => next_material_property
          else
            list => next_material_property
          endif
          cur_material_property%next => next_material_property%next
          next_material_property%next => cur_material_property
          swapped = PETSC_TRUE
        endif
      endif
      prev_material_property => cur_material_property
      cur_material_property => next_material_property
    enddo
  enddo
#endif

  max_id = 0
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    max_id = max(max_id,cur_material_property%id)
    cur_material_property => cur_material_property%next
  enddo
  
  allocate(array(max_id))
  do i = 1, max_id
    nullify(array(i)%ptr)
  enddo
  
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    array(cur_material_property%id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo

end subroutine MaterialPropConvertListToArray

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

end module Material_module
