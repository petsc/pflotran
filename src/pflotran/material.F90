module Material_module
 
  use Dataset_Aux_module

  implicit none

  private

#include "definitions.h"
 
  type, public :: material_property_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: permeability(3,3)
    PetscBool :: isotropic_permeability
    PetscReal :: vertical_anisotropy_ratio ! (vertical / horizontal)
    PetscReal :: permeability_scaling_factor
    character(len=MAXWORDLENGTH) :: permeability_dataset_name
    type(dataset_type), pointer :: permeability_dataset
    PetscReal :: porosity
    character(len=MAXWORDLENGTH) :: porosity_dataset_name
    type(dataset_type), pointer :: porosity_dataset
    PetscReal :: tortuosity
    PetscInt :: saturation_function_id
    character(len=MAXWORDLENGTH) :: saturation_function_name
    PetscReal :: rock_density ! kg/m^3
    PetscReal :: specific_heat ! J/kg-K
    PetscReal :: thermal_conductivity_dry
    PetscReal :: thermal_conductivity_wet
    PetscReal :: alpha    ! conductivity saturation relation exponent
    PetscReal :: youngs_modulus
    PetscReal :: poissons_ratio
#ifdef ICE
    PetscReal :: thermal_conductivity_frozen
    PetscReal :: alpha_fr
#endif
    PetscReal :: pore_compressibility
    PetscReal :: thermal_expansitivity   
    PetscReal :: longitudinal_dispersivity 
    PetscReal :: tortuosity_pwr
    !geh: minral surface area power functions must be defined on a per
    !     mineral basis, look in reaction_aux.F90
    !PetscReal :: mnrl_surf_area_volfrac_pwr
    !PetscReal :: mnrl_surf_area_porosity_pwr
    PetscReal :: permeability_pwr
    PetscReal :: permeability_crit_por
    PetscReal :: permeability_min_scale_fac
    character(len=MAXWORDLENGTH) :: secondary_continuum_name
    PetscReal :: secondary_continuum_length
    PetscReal :: secondary_continuum_matrix_block_size
    PetscReal :: secondary_continuum_fracture_spacing
    PetscReal :: secondary_continuum_radius
    PetscReal :: secondary_continuum_area
    PetscInt :: secondary_continuum_ncells
    PetscReal :: secondary_continuum_epsilon
    PetscReal :: secondary_continuum_aperture
    PetscReal :: secondary_continuum_init_temp
    PetscReal :: secondary_continuum_init_conc
    PetscReal :: secondary_continuum_porosity
    PetscReal :: secondary_continuum_diff_coeff
    PetscReal :: secondary_continuum_mnrl_volfrac
    PetscReal :: secondary_continuum_mnrl_area 
    PetscBool :: secondary_continuum_log_spacing
    PetscReal :: secondary_continuum_outer_spacing
    PetscReal :: secondary_continuum_area_scaling
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
  material_property%isotropic_permeability = PETSC_TRUE
  material_property%vertical_anisotropy_ratio = 0.d0
  material_property%permeability_scaling_factor = 0.d0
  material_property%permeability_pwr = 0.d0
  material_property%permeability_crit_por = 0.d0
  material_property%permeability_min_scale_fac = 1.d0
  material_property%permeability_dataset_name = ''
  nullify(material_property%permeability_dataset)
  material_property%porosity = 0.d0
  material_property%porosity_dataset_name = ''
  nullify(material_property%porosity_dataset)
  material_property%tortuosity = 1.d0
  material_property%tortuosity_pwr = 0.d0
  material_property%saturation_function_id = 0
  material_property%saturation_function_name = ''
  material_property%rock_density = 0.d0
  material_property%specific_heat = 0.d0
  material_property%thermal_conductivity_dry = 0.d0
  material_property%thermal_conductivity_wet = 0.d0
  material_property%alpha = 0.45d0
  material_property%youngs_modulus = 2.d11 ! in Pa
  material_property%poissons_ratio = 0.3
#ifdef ICE
  material_property%thermal_conductivity_frozen = 0.d0
  material_property%alpha_fr = 0.95d0
#endif
  material_property%pore_compressibility = 0.d0
  material_property%thermal_expansitivity = 0.d0  
  material_property%longitudinal_dispersivity = 0.d0
  material_property%secondary_continuum_name = ''
  material_property%secondary_continuum_length = 0.d0
  material_property%secondary_continuum_matrix_block_size = 0.d0
  material_property%secondary_continuum_fracture_spacing = 0.d0
  material_property%secondary_continuum_radius = 0.d0
  material_property%secondary_continuum_area = 0.d0
  material_property%secondary_continuum_epsilon = 1.d0
  material_property%secondary_continuum_aperture = 0.d0
  material_property%secondary_continuum_init_temp = 100.d0
  material_property%secondary_continuum_init_conc = 0.d0
  material_property%secondary_continuum_porosity = 0.5d0
  material_property%secondary_continuum_diff_coeff = 1.d-9
  material_property%secondary_continuum_mnrl_volfrac = 0.d0
  material_property%secondary_continuum_mnrl_area = 0.d0
  material_property%secondary_continuum_ncells = 0
  material_property%secondary_continuum_log_spacing = PETSC_FALSE
  material_property%secondary_continuum_outer_spacing = 1.d-3
  material_property%secondary_continuum_area_scaling = 1.d0
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
  character(len=MAXSTRINGLENGTH) :: string
  type(dataset_type), pointer :: dataset

  PetscInt :: length

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
        call InputReadWord(input,option, &
                           material_property%saturation_function_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation function name', &
                           'MATERIAL_PROPERTY')
      case('ROCK_DENSITY') 
        call InputReadDouble(input,option,material_property%rock_density)
        call InputErrorMsg(input,option,'rock density','MATERIAL_PROPERTY')
      case('SPECIFIC_HEAT') 
        call InputReadDouble(input,option,material_property%specific_heat)
        call InputErrorMsg(input,option,'specific heat','MATERIAL_PROPERTY')
      case('LONGITUDINAL_DISPERSIVITY') 
        call InputReadDouble(input,option,material_property%longitudinal_dispersivity)
        call InputErrorMsg(input,option,'longitudinal_dispersivity','MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_DRY') 
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_dry)
        call InputErrorMsg(input,option,'dry thermal conductivity', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_CONDUCTIVITY_WET') 
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_wet)
        call InputErrorMsg(input,option,'wet thermal conductivity', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_COND_EXPONENT') 
        call InputReadDouble(input,option, &
                             material_property%alpha)
        call InputErrorMsg(input,option,'thermal conductivity exponent', &
                           'MATERIAL_PROPERTY')
#ifdef ICE
      case('THERMAL_CONDUCTIVITY_FROZEN') 
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_frozen)
        call InputErrorMsg(input,option,'frozen thermal conductivity', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_COND_EXPONENT_FROZEN') 
        call InputReadDouble(input,option, &
                             material_property%alpha_fr)
        call InputErrorMsg(input,option,'thermal conductivity frozen exponent', &
                           'MATERIAL_PROPERTY')
#endif
      case('PORE_COMPRESSIBILITY') 
        call InputReadDouble(input,option, &
                             material_property%pore_compressibility)
        call InputErrorMsg(input,option,'pore compressibility', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_EXPANSITIVITY') 
        call InputReadDouble(input,option, &
                             material_property%thermal_expansitivity)
        call InputErrorMsg(input,option,'thermal expansitivity', &
                           'MATERIAL_PROPERTY')
      case('POROSITY')
        call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'porosity','MATERIAL_PROPERTY')
        call StringToUpper(string)
        if (StringCompare(string,'DATASET',SEVEN_INTEGER)) then
          call InputReadNChars(input,option, &
                               material_property%porosity_dataset_name,&
                               MAXWORDLENGTH,PETSC_TRUE)
          call InputErrorMsg(input,option,'DATASET,NAME', &
                             'MATERIAL_PROPERTY,POROSITY')   
        else
          call InputReadDouble(string,option,material_property%porosity, &
                               input%ierr)
          call InputErrorMsg(input,option,'porosity','MATERIAL_PROPERTY')
        endif
      case('TORTUOSITY')
        call InputReadDouble(input,option,material_property%tortuosity)
        call InputErrorMsg(input,option,'tortuosity','MATERIAL_PROPERTY')
      case('PERMEABILITY')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,PERMEABILITY')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,PERMEABILITY')   
          select case(trim(word))
            case('ANISOTROPIC')
              material_property%isotropic_permeability = PETSC_FALSE
            case('VERTICAL_ANISOTROPY_RATIO')
              material_property%isotropic_permeability = PETSC_FALSE
              call InputReadDouble(input,option, &
                                   material_property%vertical_anisotropy_ratio)
              call InputErrorMsg(input,option,'vertical anisotropy ratio', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('ISOTROPIC')
              material_property%isotropic_permeability = PETSC_TRUE
            case('PERMEABILITY_SCALING_FACTOR')
              call InputReadDouble(input,option, &
                                   material_property%permeability_scaling_factor)
              call InputErrorMsg(input,option,'permeability scaling factor', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_X')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,1))
              call InputErrorMsg(input,option,'x permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Y')
              call InputReadDouble(input,option, &
                                   material_property%permeability(2,2))
              call InputErrorMsg(input,option,'y permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_Z')
              call InputReadDouble(input,option, &
                                   material_property%permeability(3,3))
              call InputErrorMsg(input,option,'z permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_XZ')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,3))
              call InputErrorMsg(input,option,'xz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_XY')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,2))
              call InputErrorMsg(input,option,'xy permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_YZ')
              call InputReadDouble(input,option, &
                                   material_property%permeability(2,3))
              call InputErrorMsg(input,option,'yz permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
            case('PERM_ISO')
              call InputReadDouble(input,option, &
                                   material_property%permeability(1,1))
              call InputErrorMsg(input,option,'isotropic permeability', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')
              material_property%permeability(2,2) = &
                material_property%permeability(1,1)
              material_property%permeability(3,3) = &
                material_property%permeability(1,1)
            case('RANDOM_DATASET')
              option%io_buffer = 'RANDOM_DATASET is no longer supported.  ' // &
                'Please use the new DATASET object in the input file and ' // &
                'reference that dataset through "DATASET name" within ' // &
                'the PERMEABILITY card.'
              call printErrMsg(option)
            case('DATASET')
              call InputReadNChars(input,option, &
                                   material_property%permeability_dataset_name,&
                                   MAXWORDLENGTH,PETSC_TRUE)
              call InputErrorMsg(input,option,'DATASET,NAME', &
                                 'MATERIAL_PROPERTY,PERMEABILITY')   
            case default
              option%io_buffer = 'Keyword (' // trim(word) // &
                                 ') not recognized in MATERIAL_PROPERTY,' // &
                                 'PERMEABILITY'
              call printErrMsg(option)
          end select
        enddo
      case('YOUNGS_MODULUS') 
        call InputReadDouble(input,option, &
                             material_property%youngs_modulus)
        call InputErrorMsg(input,option,'youngs modulus', &
                           'MATERIAL_PROPERTY')
      case('POISSONS_RATIO') 
        call InputReadDouble(input,option, &
                             material_property%poissons_ratio)
        call InputErrorMsg(input,option,'poissons_ratio', &
                           'MATERIAL_PROPERTY')
      case('PERMEABILITY_POWER')
        call InputReadDouble(input,option, &
                             material_property%permeability_pwr)
        call InputErrorMsg(input,option,'permeability power','MATERIAL_PROPERTY')
      case('PERMEABILITY_CRIT_POROSITY')
        call InputReadDouble(input,option, &
                             material_property%permeability_crit_por)
        call InputErrorMsg(input,option,'permeability critical porosity','MATERIAL_PROPERTY')
      case('PERMEABILITY_MIN_SCALE_FAC')
        call InputReadDouble(input,option, &
                             material_property%permeability_min_scale_fac)
        call InputErrorMsg(input,option,'permeability min scale factor','MATERIAL_PROPERTY')
      case('TORTUOSITY_POWER')
        call InputReadDouble(input,option, &
                             material_property%tortuosity_pwr)
        call InputErrorMsg(input,option,'tortuosity power','MATERIAL_PROPERTY')
      case('MINERAL_SURFACE_AREA_POWER')
        option%io_buffer = 'Adjustment of mineral surface area based on ' // &
          'mineral volume fraction or porosity must be performed on a ' // &
          'per mineral basis under the MINERAL_KINETICS card.  See ' // &
          'reaction_aux.F90.'
          call printErrMsg(option)
      case('SECONDARY_CONTINUUM')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,SECONDARY_CONTINUUM')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,PERMEABILITY')   
          select case(trim(word))
            case('TYPE')
              call InputReadNChars(input,option, &
                                   material_property%secondary_continuum_name,&
                                   MAXWORDLENGTH,PETSC_TRUE)
              call InputErrorMsg(input,option,'type', &
                                'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('MATRIX_BLOCK_SIZE')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_matrix_block_size)
              call InputErrorMsg(input,option,'matrix_block_size', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('FRACTURE_SPACING')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_fracture_spacing)
              call InputErrorMsg(input,option,'fracture_spacing', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('RADIUS')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_radius)
              call InputErrorMsg(input,option,'radius', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('LENGTH')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_length)
              call InputErrorMsg(input,option,'length', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('AREA')
              call InputReadDouble(input,option, &
                                   material_property%secondary_continuum_area)
              call InputErrorMsg(input,option,'area', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('NUM_CELLS')
              call InputReadInt(input,option, &
                                   material_property%secondary_continuum_ncells)
              call InputErrorMsg(input,option,'number of cells', &
                                 'MATERIAL_PROPERTY, SECONDARY_CONTINUUM')
            case('EPSILON')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_epsilon)
              call InputErrorMsg(input,option,'epsilon', &
                           'MATERIAL_PROPERTY')
            case('APERTURE')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_aperture)
              call InputErrorMsg(input,option,'aperture', &
                           'MATERIAL_PROPERTY')
            case('TEMPERATURE')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_init_temp)
              call InputErrorMsg(input,option,'secondary continuum init temp', &
                           'MATERIAL_PROPERTY')
              option%set_secondary_init_temp = PETSC_TRUE
            case('CONCENTRATION')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_init_conc)
              call InputErrorMsg(input,option,'secondary continuum init conc', &
                           'MATERIAL_PROPERTY')
              option%set_secondary_init_conc = PETSC_TRUE
            case('POROSITY')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_porosity)
              call InputErrorMsg(input,option,'secondary continuum porosity', &
                           'MATERIAL_PROPERTY')
            case('DIFFUSION_COEFFICIENT')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_diff_coeff)
              call InputErrorMsg(input,option,'secondary continuum diff coeff', &
                           'MATERIAL_PROPERTY')
            case('MINERAL_VOLFRAC')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_mnrl_volfrac)
              call InputErrorMsg(input,option,'secondary cont. mnrl volfrac.', &
                           'MATERIAL_PROPERTY')  
            case('MINERAL_AREA')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_mnrl_area)
              call InputErrorMsg(input,option,'secondary cont. mnrl area', &
                           'MATERIAL_PROPERTY')
            case('LOG_GRID_SPACING')
              material_property%secondary_continuum_log_spacing = PETSC_TRUE
            case('OUTER_SPACING')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_outer_spacing)
              call InputErrorMsg(input,option,'secondary cont. outer spacing', &
                           'MATERIAL_PROPERTY')
            case('AREA_SCALING_FACTOR')
              call InputReadDouble(input,option, &
                             material_property%secondary_continuum_area_scaling)
              call InputErrorMsg(input,option,'secondary area scaling factor', &
                           'MATERIAL_PROPERTY')
            case default
              option%io_buffer = 'Keyword (' // trim(word) // &
                                 ') not recognized in MATERIAL_PROPERTY,' // &
                                 'SECONDARY_CONTINUUM'
              call printErrMsg(option)
          end select
        enddo

      case default
        option%io_buffer = 'Keyword (' // trim(keyword) // &
                           ') not recognized in material_property'    
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
subroutine MaterialPropConvertListToArray(list,array,option)

  use Option_module
  use String_module

  implicit none
  
  type(material_property_type), pointer :: list
  type(material_property_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  type(material_property_type), pointer :: cur_material_property
  type(material_property_type), pointer :: prev_material_property
  type(material_property_type), pointer :: next_material_property
  PetscInt :: i, j, length1,length2, max_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

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
  
  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_id))
  id_count = 0
  
  cur_material_property => list
  do 
    if (.not.associated(cur_material_property)) exit
    id_count(cur_material_property%id) = &
      id_count(cur_material_property%id) + 1
    array(cur_material_property%id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo
  
  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Material ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call printMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Material IDs.'
    call printErrMsg(option)
  endif
  
  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Material name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call printMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Material names.'
    call printErrMsg(option)
  endif
  
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
    if (.not.associated(material_property_array(imaterial_property)%ptr)) cycle
    if (length == &
        len_trim(material_property_array(imaterial_property)%ptr%name) .and. &
        StringCompare(material_property_array(imaterial_property)%ptr%name, &
                        material_property_name,length)) then
      MaterialPropGetPtrFromArray => &
        material_property_array(imaterial_property)%ptr
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
  
  ! simply nullify since the datasets reside in a list within realization
  nullify(material_property%permeability_dataset)
  nullify(material_property%porosity_dataset)
    
  deallocate(material_property)
  nullify(material_property)
  
end subroutine MaterialPropertyDestroy

end module Material_module
