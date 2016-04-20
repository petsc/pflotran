module PM_Waste_Form_class

  use PM_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Geometry_module
  use Data_Mediator_Vec_class
  use Dataset_Base_class
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscBool, public :: bypass_warning_message = PETSC_FALSE

! --------------- waste form species packages ---------------------------------
  type, public :: rad_species_type
   PetscReal :: formula_weight
   PetscReal :: decay_constant
   PetscReal :: mass_fraction
   PetscReal :: inst_release_fraction
   PetscInt :: parent_id
   character(len=MAXWORDLENGTH) :: parent
   PetscInt :: column_id
   PetscInt :: ispecies
   character(len=MAXWORDLENGTH) :: name
  end type rad_species_type

! --------------- waste form mechanism types ----------------------------------
  type :: wf_mechanism_base_type
    type(rad_species_type), pointer :: rad_species_list(:)
    PetscInt :: num_species
    PetscBool :: canister_degradation_model
    PetscReal :: vitality_rate_mean
    PetscReal :: vitality_rate_stdev
    PetscReal :: vitality_rate_trunc
    PetscReal :: canister_material_constant
    PetscReal :: matrix_density ! kg/m^3
    character(len=MAXWORDLENGTH) :: name
    class(wf_mechanism_base_type), pointer :: next
  contains
    procedure, public :: Dissolution => WFMechBaseDissolution
  end type wf_mechanism_base_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_glass_type
    PetscReal :: specific_surface_area   
    PetscReal :: dissolution_rate         ! kg-glass/m^2/day
  contains
    procedure, public :: Dissolution => WFMechGlassDissolution
  end type wf_mechanism_glass_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_dsnf_type
    PetscReal :: frac_dissolution_rate    ! 1/day
  contains
    procedure, public :: Dissolution => WFMechDSNFDissolution
  end type wf_mechanism_dsnf_type

  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_fmdm_type
    PetscReal :: specific_surface_area 
    PetscReal :: dissolution_rate         ! kg-matrix/m^2/day
    PetscReal :: frac_dissolution_rate    ! 1/day
  contains
    procedure, public :: Dissolution => WFMechFMDMDissolution
  end type wf_mechanism_fmdm_type
  
  type, public, extends(wf_mechanism_base_type) :: wf_mechanism_custom_type
    PetscReal :: specific_surface_area 
    PetscReal :: dissolution_rate         ! kg-matrix/m^2/day
    PetscReal :: frac_dissolution_rate    ! 1/day
  contains
    procedure, public :: Dissolution => WFMechCustomDissolution
  end type wf_mechanism_custom_type

! --------------- waste form types --------------------------------------------
  type :: waste_form_base_type
    PetscInt :: id
    PetscInt :: local_cell_id
    type(point3d_type) :: coordinate
    PetscReal :: volume
    PetscReal :: exposure_factor
    PetscReal :: eff_dissolution_rate                   ! kg-matrix/sec
    PetscReal, pointer :: instantaneous_mass_rate(:)    ! mol/sec
    PetscReal, pointer :: cumulative_mass(:)            ! mol
    PetscReal, pointer :: rad_mass_fraction(:)          ! g-rad/g-matrix
    PetscReal, pointer :: rad_concentration(:)          ! mol-rad/g-matrix
    PetscReal, pointer :: inst_release_amount(:)        ! of rad
    PetscBool :: canister_degradation_flag
    PetscReal :: canister_vitality
    PetscReal :: canister_vitality_rate
    PetscReal :: eff_canister_vit_rate
    PetscBool :: breached
    character(len=MAXWORDLENGTH) :: mech_name
    class(wf_mechanism_base_type), pointer :: mechanism
    class(waste_form_base_type), pointer :: next
  end type waste_form_base_type

! --------------- waste form process model ------------------------------------
  type, public, extends(pm_base_type) :: pm_waste_form_type
    class(realization_subsurface_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: data_mediator_species
    class(data_mediator_vec_type), pointer :: data_mediator
    class(waste_form_base_type), pointer :: waste_form_list
    class(wf_mechanism_base_type), pointer :: mechanism_list
    PetscBool :: print_mass_balance
  contains
    procedure, public :: PMWFSetRealization
    procedure, public :: Setup => PMWFSetup
    procedure, public :: Read => PMWFRead
    procedure, public :: InitializeRun => PMWFInitializeRun
    procedure, public :: InitializeTimestep => PMWFInitializeTimestep
    procedure, public :: FinalizeTimestep => PMWFFinalizeTimestep
    procedure, public :: UpdateSolution => PMWFUpdateSolution
    procedure, public :: Solve => PMWFSolve
    procedure, public :: Checkpoint => PMWFCheckpoint    
    procedure, public :: Restart => PMWFRestart  
    procedure, public :: InputRecord => PMWFInputRecord
    procedure, public :: Destroy => PMWFDestroy
  end type pm_waste_form_type
  
  public :: PMWFCreate, &
            PMWFSetup, &
            MechanismGlassCreate, &
            MechanismDSNFCreate, &
            MechanismCustomCreate, &
            MechanismFMDMCreate, &
            RadSpeciesCreate
  
contains

! ************************************************************************** !

subroutine MechanismInit(this)
  ! 
  ! Initializes the base waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

  class(wf_mechanism_base_type) :: this

  nullify(this%next)
  nullify(this%rad_species_list)
  this%num_species = 0
  this%matrix_density = UNINITIALIZED_DOUBLE
  this%name = ''
 !---- canister degradation model ----------------------
  this%canister_degradation_model = PETSC_FALSE
  this%vitality_rate_mean = UNINITIALIZED_DOUBLE
  this%vitality_rate_stdev = UNINITIALIZED_DOUBLE
  this%vitality_rate_trunc = UNINITIALIZED_DOUBLE
  this%canister_material_constant = UNINITIALIZED_DOUBLE
 !------------------------------------------------------

end subroutine MechanismInit

! ************************************************************************** !

function MechanismGlassCreate()
  ! 
  ! Creates the glass waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_glass_type), pointer :: MechanismGlassCreate
  
  allocate(MechanismGlassCreate)
  call MechanismInit(MechanismGlassCreate)
  MechanismGlassCreate%specific_surface_area = UNINITIALIZED_DOUBLE  ! m^2/m^3
  MechanismGlassCreate%dissolution_rate = 0.d0  ! kg / sec

end function MechanismGlassCreate

! ************************************************************************** !

function MechanismDSNFCreate()
  ! 
  ! Creates the DSNF (DOE Spent Nuclear Fuel) waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_dsnf_type), pointer :: MechanismDSNFCreate
  
  allocate(MechanismDSNFCreate)
  call MechanismInit(MechanismDSNFCreate)
  MechanismDSNFCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day

end function MechanismDSNFCreate

! ************************************************************************** !

function MechanismFMDMCreate()
  ! 
  ! Creates the FMDM waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_fmdm_type), pointer :: MechanismFMDMCreate
  
  allocate(MechanismFMDMCreate)
  call MechanismInit(MechanismFMDMCreate)
  MechanismFMDMCreate%specific_surface_area = UNINITIALIZED_DOUBLE  ! m^2/m^3
  MechanismFMDMCreate%dissolution_rate = UNINITIALIZED_DOUBLE  ! kg/sec
  MechanismFMDMCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day

end function MechanismFMDMCreate

! ************************************************************************** !

function MechanismCustomCreate()
  ! 
  ! Creates the 'custom' waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none
  
  class(wf_mechanism_custom_type), pointer :: MechanismCustomCreate
  
  allocate(MechanismCustomCreate)
  call MechanismInit(MechanismCustomCreate)
  MechanismCustomCreate%specific_surface_area = UNINITIALIZED_DOUBLE  ! m^2/m^3
  MechanismCustomCreate%dissolution_rate = UNINITIALIZED_DOUBLE  ! kg/sec
  MechanismCustomCreate%frac_dissolution_rate = UNINITIALIZED_DOUBLE  ! 1/day

end function MechanismCustomCreate

! ************************************************************************** !

function RadSpeciesCreate()
  ! 
  ! Creates a radioactive species in the waste form mechanism package
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/09/16

  implicit none
  
  type(rad_species_type) :: RadSpeciesCreate

  RadSpeciesCreate%name = ''
  RadSpeciesCreate%parent = ''
  RadSpeciesCreate%parent_id = UNINITIALIZED_INTEGER
  RadSpeciesCreate%formula_weight = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%decay_constant = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%mass_fraction = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%inst_release_fraction = UNINITIALIZED_DOUBLE
  RadSpeciesCreate%column_id = UNINITIALIZED_INTEGER
  RadSpeciesCreate%ispecies = UNINITIALIZED_INTEGER

end function RadSpeciesCreate

! ************************************************************************** !

function WasteFormCreate()
  ! 
  ! Creates a waste form and initializes all parameters
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016

  implicit none

  type(waste_form_base_type), pointer :: WasteFormCreate

  allocate(WasteFormCreate)
  WasteFormCreate%id = UNINITIALIZED_INTEGER
  WasteFormCreate%local_cell_id = UNINITIALIZED_INTEGER
  WasteFormCreate%coordinate%x = UNINITIALIZED_DOUBLE
  WasteFormCreate%coordinate%y = UNINITIALIZED_DOUBLE
  WasteFormCreate%coordinate%z = UNINITIALIZED_DOUBLE
  WasteFormCreate%volume = UNINITIALIZED_DOUBLE
  WasteFormCreate%exposure_factor = 1.0d0
  WasteFormCreate%eff_dissolution_rate = UNINITIALIZED_DOUBLE
  WasteFormCreate%mech_name = ''
  nullify(WasteFormCreate%instantaneous_mass_rate) ! mol-rad/sec
  nullify(WasteFormCreate%cumulative_mass) ! mol-rad
  nullify(WasteFormCreate%rad_mass_fraction) ! g-rad/g-matrix
  nullify(WasteFormCreate%rad_concentration) ! mol-rad/g-matrix
  nullify(WasteFormCreate%inst_release_amount) ! of rad
  nullify(WasteFormCreate%mechanism)
  nullify(WasteFormCreate%next)
 !------- canister degradation model -----------------
  WasteFormCreate%canister_degradation_flag = PETSC_FALSE
  WasteFormCreate%breached = PETSC_FALSE
  WasteFormCreate%canister_vitality = 0.d0
  WasteFormCreate%canister_vitality_rate = UNINITIALIZED_DOUBLE
  WasteFormCreate%eff_canister_vit_rate = UNINITIALIZED_DOUBLE
 !----------------------------------------------------

end function WasteFormCreate

! ************************************************************************** !

function PMWFCreate()
  ! 
  ! Creates and initializes the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/15/15, 07/20/15
  ! Notes: Modified by Jenn Frederick 03/24/2016

  implicit none
  
  class(pm_waste_form_type), pointer :: PMWFCreate
  
  allocate(PMWFCreate)
  nullify(PMWFCreate%realization)
  nullify(PMWFCreate%data_mediator)
  nullify(PMWFCreate%waste_form_list)
  nullify(PMWFCreate%mechanism_list)  
  PMWFCreate%print_mass_balance = PETSC_FALSE

  call PMBaseInit(PMWFCreate)

end function PMWFCreate

! ************************************************************************** !

subroutine PMWFRead(this,input)
  ! 
  ! Reads input file parameters associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/24/2016

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  
  class(waste_form_base_type), pointer :: cur_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: matched

  option => this%option
  input%ierr = 0
  error_string = 'WASTE_FORM_GENERAL'

  option%io_buffer = 'pflotran card:: ' // trim(error_string)
  call printMsg(option)

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadMechanism(this,input,option,word,error_string,found)
    if (found) cycle
    
    error_string = 'WASTE_FORM_GENERAL'
    call PMWFReadWasteForm(this,input,option,word,error_string,found)
    if (found) cycle
   
    select case(trim(word))
    !-------------------------------------
      case('PRINT_MASS_BALANCE')
        this%print_mass_balance = PETSC_TRUE
    !-------------------------------------
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    !-------------------------------------
    end select
  enddo

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    matched = PETSC_FALSE
    cur_mechanism => this%mechanism_list
    do
      if (.not.associated(cur_mechanism)) exit
      if (StringCompare(trim(cur_waste_form%mech_name), &
                        trim(cur_mechanism%name))) then
        cur_waste_form%mechanism => cur_mechanism
        matched = PETSC_TRUE
      endif
      if (matched) exit
      cur_mechanism => cur_mechanism%next
    enddo
    if (.not.associated(cur_waste_form%mechanism)) then
      option%io_buffer = 'WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mech_name) // &
                         ' not found amoung given mechanism names.'
      call printErrMsg(option)
    endif
    if (initialized(cur_waste_form%canister_vitality_rate) .and. &
        .not. cur_waste_form%mechanism%canister_degradation_model) then
      option%io_buffer = 'WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mech_name) // &
                         ' does not have the canister degradation model turned &
                         &on, but at least one of the waste forms assigned to &
                         &this mechanism specifies a canister vitality rate.'
      call printErrMsg(option)
    endif
    ! both waste form and mechanism canister vitality rate parameters 
    ! are specified:
    if (initialized(cur_waste_form%canister_vitality_rate) .and. &
        ( initialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
          initialized(cur_waste_form%mechanism%vitality_rate_trunc) )) then
      option%io_buffer = 'Either CANISTER_VITALITY_RATE within the &
                         &WASTE_FORM blocks -or- the VITALITY_LOG10_MEAN, &
                         &VITALITY_LOG10_STDEV, and VITALITY_UPPER_TRUNCATION & 
                         &within the WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mechanism%name) // &
                         ' block should be specified, but not both.'
      call printErrMsg(option)
    endif
    ! the canister degradation model is on, but neither canister vitality
    ! rate parameters were given:
    if (cur_waste_form%mechanism%canister_degradation_model) then 
      if ( (uninitialized(cur_waste_form%mechanism%vitality_rate_mean) .or. &
            uninitialized(cur_waste_form%mechanism%vitality_rate_stdev) .or. &
            uninitialized(cur_waste_form%mechanism%vitality_rate_trunc) ) .and. &
          uninitialized(cur_waste_form%canister_vitality_rate) )  then 
        option%io_buffer = 'CANISTER_VITALITY_RATE within the &
                         &WASTE_FORM blocks -or- the VITALITY_LOG10_MEAN, &
                         &VITALITY_LOG10_STDEV, and VITALITY_UPPER_TRUNCATION & 
                         &within the WASTE_FORM MECHANISM ' // &
                         trim(cur_waste_form%mechanism%name) // &
                         ' block should be specified (but not both).'
        call printErrMsg(option)
      endif
    endif
    ! the canister degradation model is on, but neither canister vitality
    ! rate parameters were given:
    if (uninitialized(cur_waste_form%canister_vitality_rate) .and. &
        cur_waste_form%canister_degradation_flag) then
      if (uninitialized(cur_waste_form%mechanism%vitality_rate_mean)) then
        option%io_buffer = 'VITALITY_LOG10_MEAN must be given in the '&
                            // trim(error_string) // ' ' // &
                            trim(cur_waste_form%mechanism%name) // &
                            ', CANISTER_DEGRADATION_MODEL block.'
        call printErrMsg(option)
      endif
      if (uninitialized(cur_waste_form%mechanism%vitality_rate_stdev)) then
        option%io_buffer = 'VITALITY_LOG10_STDEV must be given in the '&
                           // trim(error_string) // ' ' // &
                           trim(cur_waste_form%mechanism%name) // &
                           ', CANISTER_DEGRADATION_MODEL block.'
        call printErrMsg(option)
      endif
      if (uninitialized(cur_waste_form%mechanism%vitality_rate_trunc)) then
        option%io_buffer = 'VITALITY_UPPER_TRUNCATION must be given in the '&
                           // trim(error_string) // ' ' // &
                           trim(cur_waste_form%mechanism%name) // &
                           ', CANISTER_DEGRADATION_MODEL block.'
        call printErrMsg(option)
      endif 
    endif
    cur_waste_form => cur_waste_form%next
  enddo
    
end subroutine PMWFRead

! ************************************************************************** !

subroutine PMWFReadMechanism(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the waste form mechanism
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use String_module
  use Units_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscBool :: added
  character(len=MAXWORDLENGTH) :: word, units, internal_units
  character(len=MAXSTRINGLENGTH) :: temp_buf, string
  type(rad_species_type), pointer :: temp_species_array(:)
  class(wf_mechanism_base_type), pointer :: new_mechanism, cur_mechanism
  PetscInt :: k, j, icol
  PetscReal :: double

  error_string = trim(error_string) // ',MECHANISM'
  found = PETSC_TRUE
  added = PETSC_FALSE
  input%ierr = 0
  allocate(temp_species_array(50))
  k = 0

  select case(trim(keyword))
  !-------------------------------------
    case('MECHANISM')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'mechanism type',error_string)
      call StringToUpper(word)
      select case(trim(word))
      !---------------------------------
        case('GLASS')
          error_string = trim(error_string) // ' GLASS'
          allocate(new_mechanism)
          new_mechanism => MechanismGlassCreate()
      !---------------------------------
        case('DSNF')
          error_string = trim(error_string) // ' DSNF'
          allocate(new_mechanism)
          new_mechanism => MechanismDSNFCreate()
      !---------------------------------
        case('FMDM')
          option%io_buffer = 'FMDM waste form not yet implemented. Sorry! ' &
                             // trim(error_string)
          call printErrMsg(option)
          !error_string = trim(error_string) // ' FMDM'
          !allocate(new_mechanism)
          !new_mechanism => MechanismFMDMCreate()
      !---------------------------------
        case('CUSTOM')
          error_string = trim(error_string) // ' CUSTOM'
          allocate(new_mechanism)
          new_mechanism => MechanismCustomCreate()
      !---------------------------------
        case default
          option%io_buffer = 'Unrecognized mechanism type &
                             &in the ' // trim(error_string) // ' block.'
          call printErrMsg(option)
      !---------------------------------
      end select
      
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !--------------------------
          case('NAME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'mechanism name',error_string)
            call StringToUpper(word)
            new_mechanism%name = trim(word)
        !--------------------------
          case('SPECIFIC_SURFACE_AREA')
            call InputReadDouble(input,option,double)
            call InputErrorMsg(input,option,'specific surface area', &
                               error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr == 0) then
              internal_units = 'm^2/m^3'
              double = UnitsConvertToInternal(word,internal_units,option) * &
                       double
            endif
            select type(new_mechanism)
              type is(wf_mechanism_dsnf_type)
                option%io_buffer = 'SPECIFIC_SURFACE_AREA cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
              type is(wf_mechanism_custom_type)
                new_mechanism%specific_surface_area = double
              type is(wf_mechanism_glass_type)
                new_mechanism%specific_surface_area = double
              type is(wf_mechanism_fmdm_type)
                new_mechanism%specific_surface_area = double
            end select
        !--------------------------
          case('MATRIX_DENSITY')
            call InputReadDouble(input,option,new_mechanism%matrix_density)
            call InputErrorMsg(input,option,'matrix density',error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr == 0) then
              internal_units = 'kg/m^3'
              new_mechanism%matrix_density = UnitsConvertToInternal(word, &
                   internal_units,option) * new_mechanism%matrix_density
            endif
        !--------------------------
          case('FRACTIONAL_DISSOLUTION_RATE')
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                call InputReadDouble(input,option, &
                     new_mechanism%frac_dissolution_rate)
                call InputErrorMsg(input,option,'fractional dissolution rate', &
                                   error_string)
                call InputReadWord(input,option,word,PETSC_TRUE)
                if (input%ierr == 0) then
                  internal_units = 'unitless/sec'
                  new_mechanism%frac_dissolution_rate = &
                    UnitsConvertToInternal(word,internal_units,option) * &
                    new_mechanism%frac_dissolution_rate
                endif
              type is(wf_mechanism_glass_type)
                option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
              type is(wf_mechanism_dsnf_type)
                option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE cannot be &
                                   &specified for ' // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('DISSOLUTION_RATE')
            select type(new_mechanism)
              type is(wf_mechanism_custom_type)
                call InputReadDouble(input,option, &
                     new_mechanism%dissolution_rate)
                call InputErrorMsg(input,option,'dissolution rate',error_string)
                call InputReadWord(input,option,word,PETSC_TRUE)
                if (input%ierr == 0) then
                  internal_units = 'kg/m^2-sec'
                  new_mechanism%dissolution_rate = &
                    UnitsConvertToInternal(word,internal_units,option) * &
                    new_mechanism%dissolution_rate
                endif
              type is(wf_mechanism_glass_type)
                option%io_buffer = 'DISSOLUTION_RATE cannot be specified for ' &
                                   // trim(error_string)
                call printErrMsg(option)
              type is(wf_mechanism_dsnf_type)
                option%io_buffer = 'DISSOLUTION_RATE cannot be specified for ' &
                                   // trim(error_string)
                call printErrMsg(option)
            end select
        !--------------------------
          case('SPECIES')
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              k = k + 1
              temp_species_array(k) = RadSpeciesCreate() 
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'SPECIES name',error_string)
              temp_species_array(k)%name = trim(word)
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES formula weight', &
                                 error_string)
              temp_species_array(k)%formula_weight = double
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES decay rate constant', &
                                 error_string)
              temp_species_array(k)%decay_constant = double
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES initial mass fraction', &
                                 error_string)
              temp_species_array(k)%mass_fraction = double
              call InputReadDouble(input,option,double)
              call InputErrorMsg(input,option,'SPECIES instant release &
                                 &fraction',error_string)
              temp_species_array(k)%inst_release_fraction = double
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                temp_species_array(k)%parent = trim(word)
              else
                temp_species_array(k)%parent = 'no_parent'
              endif
              new_mechanism%num_species = k
            enddo
            if (k == 0) then
              option%io_buffer = 'At least one radionuclide species must be &
                                 &provided in the ' // trim(error_string) // &
                                 ', SPECIES block.'
              call printErrMsg(option)
            endif
            allocate(new_mechanism%rad_species_list(k))
            new_mechanism%rad_species_list(1:k) = temp_species_array(1:k)
            deallocate(temp_species_array)
            k = 0
            do while (k < new_mechanism%num_species)
              k = k + 1
              if (trim(new_mechanism%rad_species_list(k)%parent) == &
                  'no_parent') then
                new_mechanism%rad_species_list(k)%parent_id = 0
              else
                j = 0
                do while (j < new_mechanism%num_species)
                  j = j + 1
                  if (trim(new_mechanism%rad_species_list(k)%parent) == &
                       trim(new_mechanism%rad_species_list(j)%name)) then
                    new_mechanism%rad_species_list(k)%parent_id = j
                    exit
                  endif
                enddo
              endif
            enddo
        !--------------------------
          case('CANISTER_DEGRADATION_MODEL')
            new_mechanism%canister_degradation_model = PETSC_TRUE
            do
              call InputReadPflotranString(input,option)
              if (InputCheckExit(input,option)) exit
              call InputReadWord(input,option,word,PETSC_TRUE)
              call StringToUpper(word)
              select case(trim(word))
              case('VITALITY_LOG10_MEAN')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_mean)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &mean value',error_string)
              case('VITALITY_LOG10_STDEV')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_stdev)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &st. dev. value',error_string)
              case('VITALITY_UPPER_TRUNCATION')
                call InputReadDouble(input,option, &
                                     new_mechanism%vitality_rate_trunc)
                call InputErrorMsg(input,option,'canister vitality log-10 &
                                   &upper truncation value',error_string)
              case('CANISTER_MATERIAL_CONSTANT')
                call InputReadDouble(input,option, &
                                     new_mechanism%canister_material_constant)
                call InputErrorMsg(input,option,'canister material constant', &
                                   error_string)
              case default
                option%io_buffer = 'Keyword ' // trim(word) // &
                                   ' not recognized in the ' // &
                                   trim(error_string) // &
                                   ' CANISTER_DEGRADATION_MODEL block.'
                call printErrMsg(option)
              end select
            enddo
        !--------------------------
          case default
            call InputKeywordUnrecognized(word,error_string,option)
        !--------------------------
        end select
      enddo

     !----------- error messaging ----------------------------------------------
      if (new_mechanism%name == '') then
        option%io_buffer = 'NAME must be specified in ' // trim(error_string) &
                           // ' block.'
        call printErrMsg(option)
      endif
      select type(new_mechanism)
        type is(wf_mechanism_glass_type)
          if (uninitialized(new_mechanism%specific_surface_area)) then
            option%io_buffer = 'SPECIFIC_SURFACE_AREA must be specified in ' &
                               // trim(error_string) // ' ' // &
                               trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
        type is(wf_mechanism_custom_type)
          if (uninitialized(new_mechanism%specific_surface_area) .and. &
              uninitialized(new_mechanism%dissolution_rate) .and. &
              uninitialized(new_mechanism%frac_dissolution_rate)) then
            option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
          if ( (initialized(new_mechanism%frac_dissolution_rate) .and. &
                initialized(new_mechanism%dissolution_rate)    ) .or. &
               (uninitialized(new_mechanism%frac_dissolution_rate) .and. &
                uninitialized(new_mechanism%dissolution_rate)    ) ) then
            option%io_buffer = 'Either FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block. &
                               &Both types of dissolution rates cannot be &
                               &specified.'
            call printErrMsg(option)
          endif
          if ( (initialized(new_mechanism%specific_surface_area) .and. &
                uninitialized(new_mechanism%dissolution_rate)  ) .or. &
               (uninitialized(new_mechanism%specific_surface_area) .and. &
                initialized(new_mechanism%dissolution_rate)      ) ) then
            option%io_buffer = 'FRACTIONAL_DISSOLUTION_RATE or &
                               &DISSOLUTION_RATE with SPECIFIC_SURFACE_AREA &
                               &must be specified in ' // trim(error_string) &
                               // ' ' // trim(new_mechanism%name) // ' block.'
            call printErrMsg(option)
          endif
      end select
      if (uninitialized(new_mechanism%matrix_density)) then
        option%io_buffer = 'MATRIX_DENSITY must be specified in ' // &
                           trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // ' block.'
        call printErrMsg(option)
      endif

      if (new_mechanism%canister_degradation_model .and. &
          uninitialized(new_mechanism%canister_material_constant)) then
        option%io_buffer = 'CANISTER_MATERIAL_CONSTANT must be given in the '&
                           // trim(error_string) // ' ' // &
                           trim(new_mechanism%name) // &
                           ', CANISTER_DEGRADATION_MODEL block.'
        call printErrMsg(option)
      endif

      if (.not.associated(new_mechanism%rad_species_list)) then
        option%io_buffer = 'At least one SPECIES must be specified in the ' // &
          trim(error_string) // ' ' // trim(new_mechanism%name) // ' block.'
        call printErrMsg(option)
      endif

      if (.not.associated(this%mechanism_list)) then
        this%mechanism_list => new_mechanism
      else
        cur_mechanism => this%mechanism_list
        do
          if (.not.associated(cur_mechanism)) exit
          if (.not.associated(cur_mechanism%next)) then
            cur_mechanism%next => new_mechanism
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_mechanism => cur_mechanism%next
        enddo
      endif
      nullify(new_mechanism)
  !-------------------------------------    
    case default !(MECHANISM keyword not found)
      found = PETSC_FALSE
  !-------------------------------------
  end select

end subroutine PMWFReadMechanism

! ************************************************************************** !

subroutine PMWFReadWasteForm(this,input,option,keyword,error_string,found)
  ! 
  ! Reads input file parameters associated with the waste form 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/24/2016
  !
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class 
  use String_module
  use Units_module
  
  implicit none
  
  class(pm_waste_form_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  PetscBool :: added
  character(len=MAXWORDLENGTH) :: word, internal_units
  class(waste_form_base_type), pointer :: new_waste_form, cur_waste_form
  class(wf_mechanism_base_type), pointer :: cur_mechanism

  error_string = trim(error_string) // ',WASTE_FORM'
  found = PETSC_TRUE
  added = PETSC_FALSE

  select case(trim(keyword))
  !-------------------------------------
    case('WASTE_FORM')
      allocate(new_waste_form)
      new_waste_form => WasteFormCreate()
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(trim(word))
        !-----------------------------
          case('EXPOSURE_FACTOR')
            call InputReadDouble(input,option,new_waste_form%exposure_factor)
            call InputErrorMsg(input,option,'exposure factor',error_string)
        !-----------------------------
          case('VOLUME')
            call InputReadDouble(input,option,new_waste_form%volume)
            call InputErrorMsg(input,option,'volume',error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr == 0) then
              internal_units = 'm^3'
              new_waste_form%volume = UnitsConvertToInternal(word, &
                   internal_units,option) * new_waste_form%volume
            endif
        !-----------------------------
          case('WF_COORDINATE')
            call GeometryReadCoordinate(input,option, &
                                        new_waste_form%coordinate,error_string)
        !-----------------------------
          case('MECHANISM_NAME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'mechanism assignment',error_string)
            call StringToUpper(word)
            new_waste_form%mech_name = trim(word)
        !-----------------------------
          case('CANISTER_VITALITY_RATE')
            call InputReadDouble(input,option, &
                                 new_waste_form%canister_vitality_rate)
            call InputErrorMsg(input,option,'canister vitality rate',error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr == 0) then
              internal_units = 'unitless/sec'
              new_waste_form%canister_vitality_rate = UnitsConvertToInternal(word, &
                   internal_units,option) * &
                   new_waste_form%canister_vitality_rate
            endif
        !-----------------------------    
          case default
            call InputKeywordUnrecognized(word,error_string,option)
        !-----------------------------
        end select
      enddo
      
     ! ----------------- error messaging -------------------------------------
      if (Uninitialized(new_waste_form%volume)) then
        option%io_buffer = 'VOLUME must be specified for all waste forms.'
        call printErrMsg(option)
      endif
      if (Uninitialized(new_waste_form%coordinate%z)) then
        option%io_buffer = 'WF_COORDINATE must be specified for all waste forms.'
        call printErrMsg(option)
      endif
      if (new_waste_form%mech_name == '') then
        option%io_buffer = 'MECHANISM_NAME must be specified for &
                           &all waste forms.'
        call printErrMsg(option)
      endif
      !note: do not throw error if EXPOSURE_FACTOR isn't specified (default = 1)
      
      if (.not.associated(this%waste_form_list)) then
        this%waste_form_list => new_waste_form
      else
        cur_waste_form => this%waste_form_list
        do
          if (.not.associated(cur_waste_form)) exit
          if (.not.associated(cur_waste_form%next)) then
            cur_waste_form%next => new_waste_form
            added = PETSC_TRUE
          endif
          if (added) exit
          cur_waste_form => cur_waste_form%next
        enddo
      endif
      nullify(new_waste_form)
  !-------------------------------------
    case default
      found = PETSC_FALSE
  !-------------------------------------
  end select

  if (.not.associated(this%waste_form_list)) then
    option%io_buffer = 'At least one WASTE_FORM must be specified in the &
                       &WASTE_FORM_GENERAL block.'
    call printErrMsg(option)
  endif

end subroutine PMWFReadWasteForm

! ************************************************************************** !

subroutine PMWFSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Realization_Subsurface_class

  implicit none
  
  class(pm_waste_form_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMWFSetRealization

! ************************************************************************** !

subroutine PMWFSetup(this)
  ! 
  ! Maps waste forms to grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Option_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  class(waste_form_base_type), pointer :: cur_waste_form, prev_waste_form
  class(waste_form_base_type), pointer :: next_waste_form
  PetscInt :: i, j, k, local_id
  PetscInt :: waste_form_id
  PetscErrorCode :: ierr
  
  grid => this%realization%patch%grid
  option => this%realization%option
  
  waste_form_id = 0
  nullify(prev_waste_form)
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    waste_form_id = waste_form_id + 1
    local_id = -1
    select case(grid%itype)
      case(STRUCTURED_GRID)
        call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                            cur_waste_form%coordinate%x, &
                                            cur_waste_form%coordinate%y, &
                                            cur_waste_form%coordinate%z, &
                                            i,j,k)
        if (i > 0 .and. j > 0 .and. k > 0) then
          local_id = i + (j-1)*grid%structured_grid%nlx + &
                      (k-1)*grid%structured_grid%nlxy
        endif
      case(IMPLICIT_UNSTRUCTURED_GRID)
        call UGridGetCellFromPoint(cur_waste_form%coordinate%x, &
                                   cur_waste_form%coordinate%y, &
                                   cur_waste_form%coordinate%z, &
                                   grid%unstructured_grid,option,local_id)
      case default
          option%io_buffer = 'Only STRUCTURED_GRID and ' // &
            'IMPLICIT_UNSTRUCTURED_GRID types supported in PMGlass.'
          call printErrMsg(option)
    end select
    if (local_id > 0) then
      cur_waste_form%id = waste_form_id
      cur_waste_form%local_cell_id = local_id
      prev_waste_form => cur_waste_form
      cur_waste_form => cur_waste_form%next
    else
      ! remove waste form
      next_waste_form => cur_waste_form%next
      if (associated(prev_waste_form)) then
        prev_waste_form%next => next_waste_form
      else
        this%waste_form_list => next_waste_form
      endif
      deallocate(cur_waste_form)
      cur_waste_form => next_waste_form
    endif
  enddo
  
end subroutine PMWFSetup

! ************************************************************************** !

 subroutine PMWFInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/25/15
  use Reaction_Aux_module
  use Realization_Base_class
  use Utility_module, only : GetRndNumFromNormalDist
  
  implicit none

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_type) :: this
  
  IS :: is
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: num_waste_form_cells
  PetscInt :: num_species
  PetscInt :: size_of_vec
  PetscInt :: i, j
  PetscInt :: data_mediator_species_id
  PetscInt, allocatable :: species_indices_in_residual(:)
  PetscErrorCode :: ierr

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species
    allocate(cur_waste_form%instantaneous_mass_rate(num_species))
    allocate(cur_waste_form%cumulative_mass(num_species))
    cur_waste_form%instantaneous_mass_rate = 0.d0
    cur_waste_form%cumulative_mass = 0.d0
    allocate(cur_waste_form%rad_mass_fraction(num_species))
    allocate(cur_waste_form%rad_concentration(num_species))
    allocate(cur_waste_form%inst_release_amount(num_species))
    cur_waste_form%rad_mass_fraction = &
      cur_waste_form%mechanism%rad_species_list%mass_fraction
    cur_waste_form%rad_concentration = 0.d0
    cur_waste_form%inst_release_amount = 0.d0
    do j = 1, num_species
      cur_waste_form%mechanism%rad_species_list(j)%ispecies = &
        GetPrimarySpeciesIDFromName( &
        cur_waste_form%mechanism%rad_species_list(j)%name, &
        this%realization%reaction,this%option)
    enddo
   !--------- canister degradation model --------------------
    if (cur_waste_form%mechanism%canister_degradation_model) then
      cur_waste_form%canister_degradation_flag = PETSC_TRUE
      cur_waste_form%canister_vitality = 1.d0
      if (Uninitialized(cur_waste_form%canister_vitality_rate)) then
        call GetRndNumFromNormalDist(cur_waste_form%mechanism%vitality_rate_mean, &
                                     cur_waste_form%mechanism%vitality_rate_stdev,&
                                     cur_waste_form%canister_vitality_rate)
        if (cur_waste_form%canister_vitality_rate > &
            cur_waste_form%mechanism%vitality_rate_trunc) then
          cur_waste_form%canister_vitality_rate = &
            cur_waste_form%mechanism%vitality_rate_trunc
        endif
        ! Given rates are in units of log-10/yr, so convert to 1/yr:
        cur_waste_form%canister_vitality_rate = &
          10.0**(cur_waste_form%canister_vitality_rate)
      endif
    endif
   !----------------------------------------------------------
    cur_waste_form => cur_waste_form%next
  enddo
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
  endif
  
  if (.not.this%option%restart_flag .and. this%print_mass_balance) then
    call PMWFOutputHeader(this)
    call PMWFOutput(this)
  endif

  ! set up mass transfer
  call RealizCreateTranMassTransferVec(this%realization)
  this%data_mediator => DataMediatorVecCreate()
  call this%data_mediator%AddToList(this%realization%tran_data_mediator_list)
  ! create a Vec sized by # waste packages * # primary dofs influenced by 
  ! waste package
  ! count of waste form cells
  cur_waste_form => this%waste_form_list
  num_waste_form_cells = 0
  size_of_vec = 0
  do
    if (.not.associated(cur_waste_form)) exit
    size_of_vec = size_of_vec + cur_waste_form%mechanism%num_species
    num_waste_form_cells = num_waste_form_cells + 1
    cur_waste_form => cur_waste_form%next
  enddo
  call VecCreateSeq(PETSC_COMM_SELF,size_of_vec, &
                    this%data_mediator%vec,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(this%data_mediator%vec,ierr);CHKERRQ(ierr)

  if (num_waste_form_cells > 0) then
    allocate(species_indices_in_residual(size_of_vec))
    species_indices_in_residual = 0
    cur_waste_form => this%waste_form_list
    i = 0
    do
      if (.not.associated(cur_waste_form)) exit
      do j = 1,cur_waste_form%mechanism%num_species
        i = i + 1
        species_indices_in_residual(i) = &
          (cur_waste_form%local_cell_id-1)*this%option%ntrandof + &
          cur_waste_form%mechanism%rad_species_list(j)%ispecies
      enddo
      cur_waste_form => cur_waste_form%next
    enddo                             ! zero-based indexing
    !write(*,*) species_indices_in_residual(:)
    !stop
    species_indices_in_residual(:) = species_indices_in_residual(:) - 1
    ! set to global petsc index
    species_indices_in_residual(:) = species_indices_in_residual(:) + &
      this%realization%patch%grid%global_offset*this%option%ntrandof
  endif
  call ISCreateGeneral(this%option%mycomm,size_of_vec, &
                       species_indices_in_residual, &
                       PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
  if (allocated(species_indices_in_residual)) &
    deallocate(species_indices_in_residual)
  call VecScatterCreate(this%data_mediator%vec,PETSC_NULL_OBJECT, &
                        this%realization%field%tran_r,is, &
                        this%data_mediator%scatter_ctx,ierr);CHKERRQ(ierr)
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  call PMWFSolve(this,0.d0,ierr)
  
end subroutine PMWFInitializeRun

! ************************************************************************** !

subroutine PMWFInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick 03/28/2016

  use Global_Aux_module
  use Material_Aux_class
  use Field_module
  use Option_module
  use Grid_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  class(pm_waste_form_type) :: this
  
  class(waste_form_base_type), pointer :: cur_waste_form
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  PetscReal :: rate
  PetscReal :: dV
  PetscReal :: dt
  PetscInt :: k, p
  PetscInt :: num_species
  PetscErrorCode :: ierr
  PetscInt :: cell_id, idof
  PetscReal :: parent_concentration_old
  PetscReal :: inst_release_molality
  PetscReal, parameter :: conversion = 1.d0/(24.d0*3600.d0)
  PetscReal, pointer :: xx_p(:)

  global_auxvars => this%realization%patch%aux%Global%auxvars
  material_auxvars => this%realization%patch%aux%Material%auxvars
  field => this%realization%field
  option => this%option
  grid => this%realization%patch%grid
  dt = option%tran_dt
  
  if (option%print_screen_flag) then
    write(*,'(/,2("=")," WASTE FORM MODEL ",60("="))')
  endif

  cur_waste_form => this%waste_form_list
  do 
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species
    ! ------ update mass balances after transport step ---------------------
    cur_waste_form%cumulative_mass = cur_waste_form%cumulative_mass + &
                                     cur_waste_form%instantaneous_mass_rate * dt
    ! ------ update matrix volume ------------------------------------------
    dV = cur_waste_form%eff_dissolution_rate / &      ! kg-matrix/sec
         cur_waste_form%mechanism%matrix_density * &  ! kg-matrix/m^3-matrix
         dt                                           ! sec
    cur_waste_form%volume = cur_waste_form%volume - dV
    if (cur_waste_form%volume <= 1.d-10) then
      cur_waste_form%volume = 0.d0
    endif

    k = 0
    ! ------ get species concentrations from mass fractions ----------------
    do while (k < num_species)
      k = k + 1
      if (cur_waste_form%volume <= 0.d0) then
        cur_waste_form%rad_concentration(k) = 0.d0
        cur_waste_form%rad_mass_fraction(k) = 0.d0
      else
        cur_waste_form%rad_concentration(k) = &
          cur_waste_form%rad_mass_fraction(k) / &
          cur_waste_form%mechanism%rad_species_list(k)%formula_weight
      endif
    enddo

    !---------------- vitality degradation function ------------------------
    if (cur_waste_form%canister_degradation_flag) then
      if (cur_waste_form%canister_vitality < 1.d-3) then
        cur_waste_form%canister_vitality = 0.d0
        cur_waste_form%eff_canister_vit_rate = 0.d0
      else
        cur_waste_form%eff_canister_vit_rate = &
          cur_waste_form%canister_vitality_rate * &
          exp( cur_waste_form%mechanism%canister_material_constant * &
          ( (1.d0/333.15d0) - &
          (1.d0/(global_auxvars(grid%nL2G(cur_waste_form%local_cell_id))% &
           temp+273.15d0))) )
        cur_waste_form%canister_vitality = cur_waste_form%canister_vitality &
                      - ( cur_waste_form%eff_canister_vit_rate * &   ! [1/yr]
                          dt * &                                     ! [sec]
                          (1.0/(365.0*24.0*3600.0)) )                ! [yr/sec]
        if (cur_waste_form%canister_vitality < 1.d-3) then
          cur_waste_form%canister_vitality = 0.d0
        endif
      endif
    endif

    k = 0
    !------- instantaneous release ----------------------------------------- 
    if (.not.cur_waste_form%breached .and. &
           cur_waste_form%canister_vitality == 0.d0) then
      call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
      do while (k < num_species)
        k = k + 1
        cur_waste_form%inst_release_amount(k) = &
           (cur_waste_form%mechanism%rad_species_list(k)%inst_release_fraction *&
            cur_waste_form%rad_concentration(k))
        cur_waste_form%rad_concentration(k) = &
           cur_waste_form%rad_concentration(k) - &
           cur_waste_form%inst_release_amount(k)
        ! update mass fractions after instantaneous release
        cur_waste_form%rad_mass_fraction(k) = &
           cur_waste_form%rad_concentration(k) * &
           cur_waste_form%mechanism%rad_species_list(k)%formula_weight
        ! update transport solution vector with mass injection molality
        ! as an alternative to a source term (issue with tran_dt changing)
        idof = cur_waste_form%mechanism%rad_species_list(k)%ispecies + &
               ((cur_waste_form%local_cell_id - 1) * option%ntrandof) 
        cell_id = cur_waste_form%local_cell_id
        inst_release_molality = &                      ! [mol-rad/kg-water]
           ! [mol-rad]
          (cur_waste_form%inst_release_amount(k) * &   ! [mol-rad/g-matrix]
           cur_waste_form%volume * &                   ! [m^3-matrix]
           cur_waste_form%mechanism%matrix_density * & ! [kg-matrix/m^3-matrix]
           1.d3) / &                                  ! [kg-matrix] -> [g-matrix]
           ! [kg-water]
          (material_auxvars(cell_id)%porosity * &         ! [-]
           global_auxvars(cell_id)%sat(LIQUID_PHASE) * &  ! [-]
           material_auxvars(cell_id)%volume * &           ! [m^3]
           global_auxvars(cell_id)%den_kg(LIQUID_PHASE))  ! [kg/m^3-water]
        xx_p(idof) = xx_p(idof) + inst_release_molality
      enddo
      cur_waste_form%breached = PETSC_TRUE 
      call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    endif

    if (cur_waste_form%volume >= 0.d0) then
      k = 0
      !------- decay the radionuclide species --------------------------------
      do while (k < num_species)
        k = k + 1
        ! If species has a parent, save the parent's initial concentration
        if (cur_waste_form%mechanism%rad_species_list(k)%parent_id /= 0) then
          p = cur_waste_form%mechanism%rad_species_list(k)%parent_id
          parent_concentration_old = cur_waste_form%rad_concentration(p)
        endif
        ! Decay the species
        cur_waste_form%rad_concentration(k) = &
           cur_waste_form%rad_concentration(k) * exp(-1.0d0*dt* &
           cur_waste_form%mechanism%rad_species_list(k)%decay_constant)
        ! If species has a parent, adjust new concentration due to ingrowth
        if (cur_waste_form%mechanism%rad_species_list(k)%parent_id /= 0) then
          p = cur_waste_form%mechanism%rad_species_list(k)%parent_id
          cur_waste_form%rad_concentration(k) = &
             cur_waste_form%rad_concentration(k) + &
             ((cur_waste_form%mechanism%rad_species_list(p)%decay_constant* &
             parent_concentration_old)/ &
             (cur_waste_form%mechanism%rad_species_list(k)%decay_constant &
             - cur_waste_form%mechanism%rad_species_list(p)%decay_constant)) * &
             (exp(-1.0d0*dt* &
                  cur_waste_form%mechanism%rad_species_list(p)%decay_constant) &
             - exp(-1.0d0*dt* &
                   cur_waste_form%mechanism%rad_species_list(k)%decay_constant))
        endif
      enddo
      k = 0
      ! ------ update species mass fractions ---------------------------------
      do while (k < num_species)
        k = k + 1
        cur_waste_form%rad_mass_fraction(k) = &
        cur_waste_form%rad_concentration(k) * &
          cur_waste_form%mechanism%rad_species_list(k)%formula_weight
        ! to avoid errors in plotting data when conc is very very low:  
        if (cur_waste_form%rad_mass_fraction(k) <= 1e-40) then
          cur_waste_form%rad_mass_fraction(k) = 0.d0
        endif
      enddo
    endif
    cur_waste_form => cur_waste_form%next
  enddo

  if (this%print_mass_balance) then
    call PMWFOutput(this)
  endif

end subroutine PMWFInitializeTimestep

! ************************************************************************** !

subroutine PMWFSolve(this,time,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_waste_form_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: i, j
  PetscInt :: num_species
  PetscReal, pointer :: vec_p(:)  

  call VecGetArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  cur_waste_form => this%waste_form_list
  i = 0
  do 
    if (.not.associated(cur_waste_form)) exit
    num_species = cur_waste_form%mechanism%num_species    
    if ((cur_waste_form%volume > 0.d0) .and. &
        (cur_waste_form%canister_vitality <= 1.d-40)) then
      ! calculate the mechanism-specific eff_dissolution_rate [kg-matrix/sec]
      call cur_waste_form%mechanism%Dissolution(cur_waste_form,this)
      ! mol/sec
      do j = 1,num_species
        i = i + 1
        cur_waste_form%instantaneous_mass_rate(j) = &
          (cur_waste_form%eff_dissolution_rate * &            ! kg-matrix/sec
           cur_waste_form%mechanism%rad_species_list(j)%formula_weight * &! kmol-rad/kg-rad
           cur_waste_form%rad_mass_fraction(j) * &            ! kg-rad/kg-matrix
           1.d3)                                              ! kmol -> mol
        vec_p(i) = cur_waste_form%instantaneous_mass_rate(j)  ! mol/sec
      enddo
    else ! (canister not breached, or all waste form has dissolved already)
      i = i + num_species
      cur_waste_form%eff_dissolution_rate = 0.d0
      cur_waste_form%instantaneous_mass_rate = 0.d0
    endif
    cur_waste_form => cur_waste_form%next
  enddo
  
  call VecRestoreArrayF90(this%data_mediator%vec,vec_p,ierr);CHKERRQ(ierr)
  
end subroutine PMWFSolve

! ************************************************************************** !

subroutine WFMechBaseDissolution(this,waste_form,pm) 
  ! 
  ! Calculates the waste form dissolution rate; must be extended
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  implicit none
  
  class(wf_mechanism_base_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm

  ! This routine must be extended.

end subroutine WFMechBaseDissolution

! ************************************************************************** !

subroutine WFMechGlassDissolution(this,waste_form,pm) 
  ! 
  ! Calculates the glass waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module  
  use Global_Aux_module

  implicit none
  
  class(wf_mechanism_glass_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm

  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)                                                                    ! 1/day -> 1/sec
  PetscReal, parameter :: time_conversion = 1.d0/(24.d0*3600.d0)

  grid => pm%realization%patch%grid
  global_auxvars => pm%realization%patch%aux%Global%auxvars

  ! kg glass/m^2/day
  this%dissolution_rate = 560.d0*exp(-7397.d0/ &
    (global_auxvars(grid%nL2G(waste_form%local_cell_id))%temp+273.15d0))
  
  ! kg glass / sec
  waste_form%eff_dissolution_rate = &
    this%dissolution_rate * &          ! kg-glass/m^2/day
    this%specific_surface_area * &     ! m^2/kg glass
    this%matrix_density * &            ! kg-glass/m^3-glass
    waste_form%volume * &              ! m^3-glass
    waste_form%exposure_factor * &     ! [-]
    time_conversion                    ! day/sec

end subroutine WFMechGlassDissolution

! ************************************************************************** !

subroutine WFMechDSNFDissolution(this,waste_form,pm) 
  ! 
  ! Calculates the DSNF waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module  
  use Global_Aux_module

  implicit none

  class(wf_mechanism_dsnf_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm  
                                           ! day/sec
  PetscReal, parameter :: time_conversion = 1.d0/(24.d0*3600.d0)
  
  this%frac_dissolution_rate = 1.d0 / (1.1d0*pm%realization%option%tran_dt) 

  ! kg matrix/sec
  waste_form%eff_dissolution_rate = &
    this%frac_dissolution_rate * &           ! 1/sec
    this%matrix_density * &                  ! kg matrix/m^3 matrix
    waste_form%volume * &                    ! m^3 matrix
    waste_form%exposure_factor               ! [-]

end subroutine WFMechDSNFDissolution

! ************************************************************************** !

subroutine WFMechFMDMDissolution(this,waste_form,pm)
  !
  ! Calculates the FMDM waste form dissolution rate
  !
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module
  use Global_Aux_module

  implicit none

  class(wf_mechanism_fmdm_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm

  ! This is a placeholder routine!!!

  

end subroutine WFMechFMDMDissolution

! ************************************************************************** !

subroutine WFMechCustomDissolution(this,waste_form,pm) 
  ! 
  ! Calculates the "custom" waste form dissolution rate
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016

  use Grid_module  
  use Global_Aux_module

  implicit none

  class(wf_mechanism_custom_type) :: this
  class(waste_form_base_type) :: waste_form
  class(pm_waste_form_type) :: pm

  ! Note: Units for dissolution rates have already been converted to
  ! internal units within the PMWFRead routine.

  if (uninitialized(this%frac_dissolution_rate)) then
    ! kg glass / sec
    waste_form%eff_dissolution_rate = &
       this%dissolution_rate * &         ! kg-matrix/m^2/sec
       this%specific_surface_area * &    ! m^2/kg-matrix
       this%matrix_density * &           ! kg-matrix/m^3-matrix
       waste_form%volume * &             ! m^3-matrix
       waste_form%exposure_factor        ! [-]
  else
    ! kg matrix/sec
    waste_form%eff_dissolution_rate = &
       this%frac_dissolution_rate * &     ! [-]/sec
       this%matrix_density * &            ! kg matrix/m^3 matrix
       waste_form%volume * &              ! m^3 matrix
       waste_form%exposure_factor         ! [-]
  endif

end subroutine WFMechCustomDissolution

! ************************************************************************** !

subroutine PMWFFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
end subroutine PMWFFinalizeTimestep

! ************************************************************************** !

subroutine PMWFUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update glass mass here?

end subroutine PMWFUpdateSolution

! ************************************************************************** !

recursive subroutine PMWFFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMWFFinalizeRun

! ************************************************************************** !

subroutine PMWFOutput(this)
  ! 
  ! Sets up output for a waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module
  use Output_Aux_module
  use Global_Aux_module
  use Grid_module

  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  class(waste_form_base_type), pointer :: cur_waste_form
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: i
  
  if (.not.associated(this%waste_form_list)) return
  
100 format(100es18.8)

  option => this%realization%option
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  global_auxvars => this%realization%patch%aux%Global%auxvars
  
  fid = 86
  filename = PMWFOutputFilename(option)
  open(unit=fid,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(fid,100,advance="no") option%time / output_option%tconv
  
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    do i = 1, cur_waste_form%mechanism%num_species
      write(fid,100,advance="no") cur_waste_form%cumulative_mass(i), &
                                  cur_waste_form%instantaneous_mass_rate(i), &
                                  cur_waste_form%rad_mass_fraction(i)
    enddo
    write(fid,100,advance="no") cur_waste_form%eff_dissolution_rate, &
                                cur_waste_form%volume, &
                                cur_waste_form%eff_canister_vit_rate, &
                                cur_waste_form%canister_vitality*100.0
    cur_waste_form => cur_waste_form%next
  enddo
  close(fid)
  
end subroutine PMWFOutput

! ************************************************************************** !

function PMWFOutputFilename(option)
  ! 
  ! Generates filename for waste form output
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Option_module

  implicit none
  
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: PMWFOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  PMWFOutputFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       '-wf_mass-' // trim(adjustl(word)) // '.dat'
  
end function PMWFOutputFilename  

! ************************************************************************** !

subroutine PMWFOutputHeader(this)
  ! 
  ! Writes header for waste form output file
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  use Output_Aux_module
  use Grid_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  class(pm_waste_form_type) :: this
  
  type(output_option_type), pointer :: output_option
  type(grid_type), pointer :: grid
  class(waste_form_base_type), pointer :: cur_waste_form
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  
  if (.not.associated(this%waste_form_list)) return
  
  output_option => this%realization%output_option
  grid => this%realization%patch%grid
  
  fid = 86
  filename = PMWFOutputFilename(this%option)
  open(unit=fid,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // ']"'

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    ! cell natural id
    write(cell_string,*) grid%nG2A(grid%nL2G(cur_waste_form%local_cell_id))
    cell_string = ' (' // trim(adjustl(cell_string)) // ')'
    ! coordinate of waste form
    x_string = BestFloat(cur_waste_form%coordinate%x,1.d4,1.d-2)
    y_string = BestFloat(cur_waste_form%coordinate%y,1.d4,1.d-2)
    z_string = BestFloat(cur_waste_form%coordinate%z,1.d4,1.d-2)
    cell_string = trim(cell_string) // &
             ' (' // trim(adjustl(x_string)) // &
             ' ' // trim(adjustl(y_string)) // &
             ' ' // trim(adjustl(z_string)) // ')'
    do i = 1, cur_waste_form%mechanism%num_species
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Cum. Mass Flux'
      ! cumulative
      units_string = 'mol'
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Inst. Mass Flux'
      ! instantaneous
      units_string = 'mol/s' !// trim(adjustl(output_option%tunit))
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)       
      variable_string = trim(cur_waste_form%mechanism%rad_species_list(i)%name) &
                        // ' Mass Frac.'
      units_string = 'g-rad/g-matrix' 
      call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                               icolumn)
    enddo
    variable_string = 'WF Dissolution Rate'
    units_string = 'kg/s' !// trim(adjustl(output_option%tunit))
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Volume'
    units_string = 'm^3'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Vitality Degradation Rate'
    units_string = '1/yr'
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)
    variable_string = 'WF Canister Vitality'
    units_string = '%' 
    call OutputWriteToHeader(fid,variable_string,units_string,cell_string, &
                             icolumn)

    cur_waste_form => cur_waste_form%next
  enddo
  
  close(fid)
  
end subroutine PMWFOutputHeader

! ***************************************************************************** !

subroutine PMWFCheckpoint(this,viewer)
  ! 
  ! Checkpoints data associated with the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  !
  use Option_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_type) :: this
  PetscViewer :: viewer
  
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: maximum_waste_form_id
  PetscInt :: local_waste_form_count
  PetscInt :: temp_int
  
  Vec :: local, global
  PetscErrorCode :: ierr
  
  this%option%io_buffer = 'PMWFCheckpoint not implemented.'
  call printErrMsg(this%option)
  
  ! calculate maximum waste form id
  maximum_waste_form_id = 0
  local_waste_form_count = 0
  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    local_waste_form_count = local_waste_form_count + 1
    maximum_waste_form_id = max(maximum_waste_form_id,cur_waste_form%id)
    cur_waste_form => cur_waste_form%next
  enddo
  call MPI_Allreduce(maximum_waste_form_id,temp_int,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,this%option%mycomm,ierr)
!  call VecCreateMPI(this%option%mycomm,local_waste_form_count,PETSC_DETERMINE,ierr)
  
                     
end subroutine PMWFCheckpoint

! ***************************************************************************** !

subroutine PMWFRestart(this,viewer)
  ! 
  ! Restarts data associated with Subsurface PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_waste_form_type) :: this
  PetscViewer :: viewer
  
  this%option%io_buffer = 'PMWFRestart not implemented.'
  call printErrMsg(this%option)
!  call RestartFlowProcessModel(viewer,this%realization)
!  call this%UpdateAuxVars()
!  call this%UpdateSolution()
  
end subroutine PMWFRestart

! ************************************************************************** !

subroutine PMWFInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_waste_form_type) :: this

  character(len=MAXWORDLENGTH) :: word
  class(waste_form_base_type), pointer :: cur_waste_form
  PetscInt :: id
  PetscInt :: k

  id = INPUT_RECORD_UNIT

  
end subroutine PMWFInputRecord

! ************************************************************************** !

subroutine PMWFStrip(this)
  ! 
  ! Strips the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15
  ! Notes: Modified by Jenn Frederick, 03/28/2016

  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_waste_form_type) :: this
  
  class(waste_form_base_type), pointer :: cur_waste_form, prev_waste_form

  nullify(this%realization)
  nullify(this%data_mediator)

  cur_waste_form => this%waste_form_list
  do
    if (.not.associated(cur_waste_form)) exit
    prev_waste_form => cur_waste_form
    cur_waste_form => cur_waste_form%next
    call DeallocateArray(prev_waste_form%rad_mass_fraction)
    call DeallocateArray(prev_waste_form%rad_concentration)
    call DeallocateArray(prev_waste_form%inst_release_amount)
    call DeallocateArray(prev_waste_form%instantaneous_mass_rate)
    call DeallocateArray(prev_waste_form%cumulative_mass)
    nullify(prev_waste_form%mechanism)
    deallocate(prev_waste_form)
    nullify(prev_waste_form)
  enddo
  nullify(this%waste_form_list)
  call PMWFMechanismStrip(this)

end subroutine PMWFStrip

! ************************************************************************** !

subroutine PMWFMechanismStrip(this)
  ! 
  ! Strips the waste form mechanisms in the waste form process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/28/2016
  !
  implicit none
  
  class(pm_waste_form_type) :: this
  
  class(wf_mechanism_base_type), pointer :: cur_mechanism, prev_mechanism

  cur_mechanism => this%mechanism_list
  do
    if (.not.associated(cur_mechanism)) exit
    prev_mechanism => cur_mechanism
    cur_mechanism => cur_mechanism%next
    deallocate(prev_mechanism%rad_species_list)
    nullify(prev_mechanism%rad_species_list)
    deallocate(prev_mechanism)
    nullify(prev_mechanism)
  enddo
  nullify(this%mechanism_list)

end subroutine PMWFMechanismStrip
  
! ************************************************************************** !

subroutine PMWFDestroy(this)
  ! 
  ! Destroys the waste form process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/26/15

  implicit none
  
  class(pm_waste_form_type) :: this
  
  call PMWFStrip(this)
  
end subroutine PMWFDestroy


  
end module PM_Waste_Form_class
