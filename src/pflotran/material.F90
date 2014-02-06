module Material_module
 
  use Dataset_Common_HDF5_class

  use PFLOTRAN_Constants_module
  use Material_Aux_class

  implicit none

  private

#include "finclude/petscsys.h"
 
  type, public :: material_property_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: permeability(3,3)
    PetscBool :: isotropic_permeability
    PetscReal :: vertical_anisotropy_ratio ! (vertical / horizontal)
    PetscReal :: permeability_scaling_factor
    character(len=MAXWORDLENGTH) :: permeability_dataset_name
    class(dataset_common_hdf5_type), pointer :: permeability_dataset
    PetscReal :: porosity
    character(len=MAXWORDLENGTH) :: porosity_dataset_name
    class(dataset_common_hdf5_type), pointer :: porosity_dataset
    PetscReal :: tortuosity
    PetscInt :: saturation_function_id
    character(len=MAXWORDLENGTH) :: saturation_function_name
    PetscReal :: rock_density ! kg/m^3
    PetscReal :: specific_heat ! J/kg-K
    PetscReal :: thermal_conductivity_dry
    PetscReal :: thermal_conductivity_wet
    PetscReal :: alpha    ! conductivity saturation relation exponent

    character(len=MAXWORDLENGTH) :: soil_compressibility_function
    PetscReal :: soil_compressibility
    PetscReal :: soil_reference_pressure

    ! ice properties
    PetscReal :: thermal_conductivity_frozen
    PetscReal :: alpha_fr

    PetscReal :: pore_compressibility
    PetscReal :: thermal_expansitivity   
    PetscReal :: dispersivity(3)
    PetscReal :: tortuosity_pwr
    PetscReal :: min_pressure
    PetscReal :: max_pressure
    PetscReal :: max_permfactor
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
  
  ! procedure pointer declarations
  procedure(MaterialCompressSoilDummy), pointer :: &
    MaterialCompressSoilPtr => null()
 
  ! interface blocks
  interface
    subroutine MaterialCompressSoilDummy(auxvar,pressure,compressed_porosity, &
                                         dcompressed_porosity_dp)
    use Material_Aux_class
    implicit none
    class(material_auxvar_type), intent(in) :: auxvar
    PetscReal, intent(in) :: pressure
    PetscReal, intent(out) :: compressed_porosity
    PetscReal, intent(out) :: dcompressed_porosity_dp
    end subroutine MaterialCompressSoilDummy
  end interface 
  
  interface MaterialCompressSoil
    procedure MaterialCompressSoilPtr
  end interface
  
  public :: MaterialPropertyCreate, &
            MaterialPropertyDestroy, &
            MaterialPropertyAddToList, &
            MaterialPropGetPtrFromList, &
            MaterialPropGetPtrFromArray, &
            MaterialPropConvertListToArray, &
            MaterialAnisotropyExists, &
            MaterialSetAuxVarScalar, &
            MaterialSetAuxVarVecLoc, &
            MaterialGetAuxVarVecLoc, &
            MaterialAuxVarCommunicate, &
            MaterialCompressSoil, &
            MaterialPropertyRead, &
            MaterialInitAuxIndices, &
            MaterialAssignPropertyToAux, &
            MaterialSetup
  
contains

! ************************************************************************** !

function MaterialPropertyCreate()
  ! 
  ! Creates a material property
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 
  
  implicit none

  type(material_property_type), pointer :: MaterialPropertyCreate
  
  type(material_property_type), pointer :: material_property
  
  allocate(material_property)
  material_property%id = 0
  material_property%name = ''
  ! initialize to -999.d0 to catch bugs
  material_property%permeability = -999.d0
  material_property%isotropic_permeability = PETSC_TRUE
  material_property%vertical_anisotropy_ratio = 0.d0
  material_property%permeability_scaling_factor = 0.d0
  material_property%permeability_pwr = 0.d0
  material_property%permeability_crit_por = 0.d0
  material_property%permeability_min_scale_fac = 1.d0
  material_property%permeability_dataset_name = ''
  nullify(material_property%permeability_dataset)
  ! initialize to -999.d0 to catch bugs
  material_property%porosity = -999.d0
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

  material_property%soil_compressibility_function = ''
  material_property%soil_compressibility = -999.d0
  material_property%soil_reference_pressure = -999.d0

  material_property%thermal_conductivity_frozen = 0.d0
  material_property%alpha_fr = 0.95d0

  material_property%pore_compressibility = -999.d0
  material_property%thermal_expansitivity = 0.d0  
  material_property%dispersivity = 0.d0
  material_property%min_pressure = 0.d0
  material_property%max_pressure = 1.d6
  material_property%max_permfactor = 1.d0
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

subroutine MaterialPropertyRead(material_property,input,option)
  ! 
  ! Reads in contents of a material_property card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  type(material_property_type) :: material_property
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: length
  PetscBool :: therm_k_frz
  PetscBool :: therm_k_exp_frz

  therm_k_frz = PETSC_FALSE
  therm_k_exp_frz = PETSC_FALSE

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

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
      case('SPECIFIC_HEAT','HEAT_CAPACITY') 
        call InputReadDouble(input,option,material_property%specific_heat)
        call InputErrorMsg(input,option,'specific heat','MATERIAL_PROPERTY')
      case('LONGITUDINAL_DISPERSIVITY') 
        call InputReadDouble(input,option,material_property%dispersivity(1))
        call InputErrorMsg(input,option,'longitudinal_dispersivity','MATERIAL_PROPERTY')
      case('TRANSVERSE_DISPERSIVITY_H') 
        call InputReadDouble(input,option,material_property%dispersivity(2))
        call InputErrorMsg(input,option,'transverse_dispersivity_h','MATERIAL_PROPERTY')
      case('TRANSVERSE_DISPERSIVITY_V') 
        call InputReadDouble(input,option,material_property%dispersivity(3))
        call InputErrorMsg(input,option,'transverse_dispersivity_v','MATERIAL_PROPERTY')
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
      case('THERMAL_CONDUCTIVITY_FROZEN') 
        therm_k_frz = PETSC_TRUE
        call InputReadDouble(input,option, &
                             material_property%thermal_conductivity_frozen)
        call InputErrorMsg(input,option,'frozen thermal conductivity', &
                           'MATERIAL_PROPERTY')
      case('THERMAL_COND_EXPONENT_FROZEN') 
        therm_k_exp_frz = PETSC_TRUE
        call InputReadDouble(input,option, &
                             material_property%alpha_fr)
        call InputErrorMsg(input,option,'thermal conductivity frozen exponent', &
                           'MATERIAL_PROPERTY')
      case('PORE_COMPRESSIBILITY') 
        call InputReadDouble(input,option, &
                             material_property%pore_compressibility)
        call InputErrorMsg(input,option,'pore compressibility', &
                           'MATERIAL_PROPERTY')
      case('SOIL_COMPRESSIBILITY_FUNCTION') 
        call InputReadWord(input,option, &
                           material_property%soil_compressibility_function, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'soil compressibility function', &
                           'MATERIAL_PROPERTY')
      case('SOIL_COMPRESSIBILITY') 
        call InputReadDouble(input,option, &
                             material_property%soil_compressibility)
        call InputErrorMsg(input,option,'soil compressibility', &
                           'MATERIAL_PROPERTY')
      case('SOIL_REFERENCE_PRESSURE') 
        call InputReadDouble(input,option, &
                             material_property%soil_reference_pressure)
        call InputErrorMsg(input,option,'soil reference pressure', &
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
          call InputReadPflotranString(input,option)
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
      case('PERM_FACTOR') 
      ! Permfactor is the multiplier to permeability to increase perm
      ! The perm increase could be due to pressure or other variable
      ! Added by Satish Karra, LANL, 1/8/12
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,PERM_FACTOR')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,PERM_FACTOR')   
          select case(trim(word))
          ! Assuming only ramp function for now
          ! The permfactor ramps from 1 to max_permfactor at max_pressure
          ! and remains same
            case('MIN_PRESSURE')       
              call InputReadDouble(input,option,material_property%min_pressure)
              call InputErrorMsg(input,option,'min pressure','PERM_FACTOR')  
            case('MAX_PRESSURE')       
              call InputReadDouble(input,option,material_property%max_pressure)
              call InputErrorMsg(input,option,'max pressure','PERM_FACTOR')
            case('MAX_PERMFACTOR')       
              call InputReadDouble(input,option,material_property%max_permfactor)
              call InputErrorMsg(input,option,'max permfactor','PERM_FACTOR')
            case default
              option%io_buffer = 'Keyword (' // trim(word) // &
                                 ') not recognized in MATERIAL_PROPERTY,' // &
                                 'PERM_FACTOR'
              call printErrMsg(option)
          end select
        enddo
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
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'MATERIAL_PROPERTY,SECONDARY_CONTINUUM')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'MATERIAL_PROPERTY,SECONDARY_CONTINUUM')   
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

  if ((option%iflowmode == TH_MODE) .or. (option%iflowmode == THC_MODE)) then
     if (option%use_th_freezing .eqv. PETSC_TRUE) then
        if (.not. therm_k_frz) then
           option%io_buffer = 'THERMAL_CONDUCTIVITY_FROZEN must be set ' // &
             'in inputdeck for MODE TH(C) ICE'
           call printErrMsg(option)
        endif
        if (.not. therm_k_exp_frz) then
           option%io_buffer = 'THERMAL_COND_EXPONENT_FROZEN must be set ' // &
             'in inputdeck for MODE TH(C) ICE'
           call printErrMsg(option)
        endif
     endif
  endif

end subroutine MaterialPropertyRead

! ************************************************************************** !

subroutine MaterialPropertyAddToList(material_property,list)
  ! 
  ! Adds a material property to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

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

subroutine MaterialPropConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the
  ! material_properties in the list
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/18/07
  ! 

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

subroutine MaterialSetup(material_parameter, material_property_array, &
                         saturation_function_array, option)
  ! 
  ! Creates arrays for material parameter boject
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/14
  !
  use Option_module
  use Saturation_Function_module
  
  implicit none
  
  type(material_parameter_type) :: material_parameter
  type(material_property_ptr_type) :: material_property_array(:)
  type(saturation_function_ptr_type) :: saturation_function_array(:)
  type(option_type), pointer :: option
  
  PetscInt :: num_sat_func
  PetscInt :: num_mat_prop
  PetscInt :: i
  
  num_mat_prop = size(material_property_array)
  num_sat_func = size(saturation_function_array)
  
  allocate(material_parameter%soil_residual_saturation(option%nphase, &
                                                       num_sat_func))
  material_parameter%soil_residual_saturation = -999.d0
  do i = 1, num_sat_func
    if (associated(saturation_function_array(i)%ptr)) then
      material_parameter%soil_residual_saturation(:, &
                         saturation_function_array(i)%ptr%id) = &
        saturation_function_array(i)%ptr%Sr(:)
    endif
  enddo

  allocate(material_parameter%soil_heat_capacity(num_mat_prop))
  allocate(material_parameter%soil_thermal_conductivity(2,num_mat_prop))
  material_parameter%soil_heat_capacity = -999.d0
  material_parameter%soil_thermal_conductivity = -999.d0
  do i = 1, num_mat_prop
    if (associated(material_property_array(i)%ptr)) then
      ! kg rock/m^3 rock * J/kg rock-K * 1.e-6 MJ/J
      material_parameter%soil_heat_capacity(i) = &
        material_property_array(i)%ptr%specific_heat * option%scale ! J -> MJ
      material_parameter%soil_thermal_conductivity(1,i) = &
        material_property_array(i)%ptr%thermal_conductivity_dry
      material_parameter%soil_thermal_conductivity(2,i) = &
        material_property_array(i)%ptr%thermal_conductivity_wet
    endif
  enddo
  
end subroutine MaterialSetup
  
! ************************************************************************** !

function MaterialPropGetPtrFromList(material_property_name, &
                                    material_property_list)
  ! 
  ! Returns a pointer to the material property
  ! matching material_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

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

function MaterialPropGetPtrFromArray(material_property_name, &
                                     material_property_array)
  ! 
  ! Returns a pointer to the material property
  ! matching material_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

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

function MaterialAnisotropyExists(material_property_list)
  ! 
  ! Determines whether any of the material
  ! properties are anisotropic
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/13
  ! 

  implicit none
  
  type(material_property_type), pointer :: material_property_list

  PetscBool :: MaterialAnisotropyExists
  
  type(material_property_type), pointer :: cur_material_property
    
  MaterialAnisotropyExists = PETSC_FALSE
  
  cur_material_property => material_property_list
  do 
    if (.not.associated(cur_material_property)) exit
    if (.not. cur_material_property%isotropic_permeability) then
      MaterialAnisotropyExists = PETSC_TRUE
      return
    endif
    cur_material_property => cur_material_property%next
  enddo
  
end function MaterialAnisotropyExists


! ************************************************************************** !

subroutine MaterialInitAuxIndices(material_property_ptrs,option)
  !
  ! Initializes the pointer used to index material property arrays
  !
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  !
  use Material_Aux_class
  use String_module
  use Option_module
  
  implicit none
  
  type(material_property_ptr_type), pointer :: material_property_ptrs(:)
  type(option_type) :: option

  PetscInt :: i
  PetscInt :: icount = 0
  
  procedure(MaterialCompressSoilDummy), pointer :: &
    MaterialCompressSoilPtrTmp 
  
!  soil_thermal_conductivity_index = 0
!  soil_heat_capacity_index = 0
  soil_compressibility_index = 0
  soil_reference_pressure_index = 0
  max_material_index = 0
  
  do i = 1, size(material_property_ptrs)
    ! if gaps exist between material ids in the input file, those gaps will
    ! be null and need to be skipped
    if (.not.associated(material_property_ptrs(i)%ptr)) cycle
    MaterialCompressSoilPtrTmp => null()
    if (len_trim(material_property_ptrs(i)%ptr% &
                   soil_compressibility_function) > 1) then
      call StringToUpper(material_property_ptrs(i)%ptr% &
                           soil_compressibility_function)
      select case(material_property_ptrs(i)%ptr%soil_compressibility_function)
        case('BRAGFLO')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilBRAGFLO
        case('LEIJNSE','DEFAULT')
          MaterialCompressSoilPtrTmp => MaterialCompressSoilLeijnse
        case default
          option%io_buffer = 'Soil compressibility function "' // &
            trim(material_property_ptrs(i)%ptr% &
                   soil_compressibility_function) // &
            '" not recognized.'
          call printErrMsg(option)
      end select
      if (.not.associated(MaterialCompressSoilPtr)) then
        MaterialCompressSoilPtr => MaterialCompressSoilPtrTmp
      else if (.not.associated(MaterialCompressSoilPtr, &
                               MaterialCompressSoilPtrTmp)) then
        option%io_buffer = 'All MATERIAL_PROPERTIES must specify the ' // &
          'same soil compressibility function.'
        call printErrMsg(option)
      endif
    endif  
    if (material_property_ptrs(i)%ptr%soil_compressibility > -998.d0 .and. &
        soil_compressibility_index == 0) then
      icount = icount + 1
      soil_compressibility_index = icount
    endif
    if (material_property_ptrs(i)%ptr%soil_reference_pressure > -998.d0 .and. &
        soil_reference_pressure_index == 0) then
      icount = icount + 1
      soil_reference_pressure_index = icount
    endif
!    if (material_property_ptrs(i)%ptr%specific_heat > 0.d0 .and. &
!        soil_heat_capacity_index == 0) then
!      icount = icount + 1
!      soil_heat_capacity_index = icount
!    endif
!    if (material_property_ptrs(i)%ptr%thermal_conductivity_wet > 0.d0 .and. &
!        soil_thermal_conductivity_index == 0) then
!      icount = icount + 1
!      soil_thermal_conductivity_index = icount
!    endif
  enddo
  max_material_index = icount
  
  if (.not.associated(MaterialCompressSoilPtr)) then
    MaterialCompressSoilPtr => MaterialCompressSoilLeijnse
  endif
  
end subroutine MaterialInitAuxIndices

! ************************************************************************** !

subroutine MaterialAssignPropertyToAux(material_auxvar,material_property, &
                                       option)
  !
  ! Initializes the pointer used to index material property arrays
  !
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  !
  use Material_Aux_class
  use Option_module
  
  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(material_property_type) :: material_property
  type(option_type) :: option

  if (material_property%rock_density > -998.d0) then
    material_auxvar%soil_particle_density = &
      material_property%rock_density
  endif
  if (soil_compressibility_index > 0) then
    material_auxvar%soil_properties(soil_compressibility_index) = &
      material_property%soil_compressibility
  endif
  if (soil_reference_pressure_index > 0) then
    material_auxvar%soil_properties(soil_reference_pressure_index) = &
      material_property%soil_reference_pressure
  endif
!  if (soil_heat_capacity_index > 0) then
!    material_auxvar%soil_properties(soil_heat_capacity_index) = &
!      material_property%specific_heat
!  endif
!  if (soil_thermal_conductivity_index > 0) then
!    material_auxvar%soil_properties(soil_thermal_conductivity_index) = &
!      material_property%thermal_conductivity_wet
!  endif
  
end subroutine MaterialAssignPropertyToAux

! ************************************************************************** !

subroutine MaterialSetAuxVarScalar(Material,value,ivar)
  ! 
  ! Sets values of a material auxvar data using a scalar value.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Variables_module
  
  implicit none

  type(material_type) :: Material ! from realization%patch%aux%Material
  PetscReal :: value
  PetscInt :: ivar

  PetscInt :: i
  
  select case(ivar)
    case(VOLUME)
      do i=1, Material%num_aux
        Material%auxvars(i)%volume = value
      enddo
    case(POROSITY)
      do i=1, Material%num_aux
        Material%auxvars(i)%porosity = value
      enddo
    case(TORTUOSITY)
      do i=1, Material%num_aux
        Material%auxvars(i)%tortuosity = value
      enddo
    case(PERMEABILITY_X)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xx_index) = value
      enddo
    case(PERMEABILITY_Y)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_yy_index) = value
      enddo
    case(PERMEABILITY_Z)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_zz_index) = value
      enddo
    case(PERMEABILITY_XY)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xy_index) = value
      enddo
    case(PERMEABILITY_YZ)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_yz_index) = value
      enddo
    case(PERMEABILITY_XZ)
      do i=1, Material%num_aux
        Material%auxvars(i)%permeability(perm_xz_index) = value
      enddo
  end select
  
end subroutine MaterialSetAuxVarScalar

! ************************************************************************** !

subroutine MaterialSetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(VOLUME)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%volume = vec_loc_p(ghosted_id)
      enddo
    case(POROSITY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%porosity = vec_loc_p(ghosted_id)
      enddo
    case(TORTUOSITY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%tortuosity = vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_X)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xx_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_Y)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_yy_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_Z)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_zz_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_XY)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xy_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_YZ)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_yz_index) = &
          vec_loc_p(ghosted_id)
      enddo
    case(PERMEABILITY_XZ)
      do ghosted_id=1, Material%num_aux
        Material%auxvars(ghosted_id)%permeability(perm_xz_index) = &
          vec_loc_p(ghosted_id)
      enddo
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr)

end subroutine MaterialSetAuxVarVecLoc

! ************************************************************************** !

subroutine MaterialGetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(VOLUME)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%volume
      enddo
    case(POROSITY)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%porosity
      enddo
    case(TORTUOSITY)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = Material%auxvars(ghosted_id)%tortuosity
      enddo
    case(PERMEABILITY_X)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xx_index)
      enddo
    case(PERMEABILITY_Y)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_yy_index)
      enddo
    case(PERMEABILITY_Z)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_zz_index)
      enddo
    case(PERMEABILITY_XY)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xy_index)
      enddo
    case(PERMEABILITY_YZ)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_yz_index)
      enddo
    case(PERMEABILITY_XZ)
      do ghosted_id=1, Material%num_aux
        vec_loc_p(ghosted_id) = &
          Material%auxvars(ghosted_id)%permeability(perm_xz_index)
      enddo
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr)

end subroutine MaterialGetAuxVarVecLoc

! ************************************************************************** !

subroutine MaterialAuxVarCommunicate(comm,Material,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of material auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Communicator_Base_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  class(communicator_type), pointer :: comm
  type(material_type) :: Material ! from realization%patch%aux%Material
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  call MaterialGetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)
  call comm%LocalToLocal(vec_loc,vec_loc)
  call MaterialSetAuxVarVecLoc(Material,vec_loc,ivar,isubvar)

end subroutine MaterialAuxVarCommunicate

! ************************************************************************** !

subroutine MaterialCompressSoilLeijnse(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Leijnse, 1992.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  PetscReal :: compression
  PetscReal :: tempreal
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compression = &
    exp(-1.d0 * compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  tempreal = (1.d0 - auxvar%porosity) * compression
  compressed_porosity = 1.d0 - tempreal
  dcompressed_porosity_dp = tempreal * compressibility
  
end subroutine MaterialCompressSoilLeijnse

! ************************************************************************** !

subroutine MaterialCompressSoilBRAGFLO(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Eq. 9.6.9 of BRAGFLO
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compressed_porosity = auxvar%porosity * &
    exp(compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  dcompressed_porosity_dp = compressibility * compressed_porosity
  
end subroutine MaterialCompressSoilBRAGFLO

! ************************************************************************** !

recursive subroutine MaterialPropertyDestroy(material_property)
  ! 
  ! Destroys a material_property
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

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
