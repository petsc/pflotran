module Condition_module
 
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Dataset_Aux_module
  use Time_Series_module
  
  use Surface_Complexation_Aux_module  
  use Mineral_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2

  type, public :: flow_condition_dataset_type
    type(time_series_type), pointer :: time_series
    type(dataset_type), pointer ::  dataset
  end type flow_condition_dataset_type
  
  type, public :: flow_condition_type
    PetscInt :: id                          ! id from which condition can be referenced
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name    ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: saturation
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: well
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: displacement_x
    type(flow_sub_condition_type), pointer :: displacement_y
    type(flow_sub_condition_type), pointer :: displacement_z
    type(flow_general_condition_type), pointer :: general
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type flow_condition_type

  ! data structure for general phase
  type, public :: flow_general_condition_type
    !TODO(geh): check to ensure that general condition is considerered 
    !           wherever sub_condition_ptr is used.
    type(flow_sub_condition_type), pointer :: liquid_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: mole_fraction
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: flux
  end type flow_general_condition_type
    
  type, public :: flow_sub_condition_type
    PetscInt :: itype                  ! integer describing type of condition
    PetscInt :: isubtype
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    type(flow_condition_dataset_type) :: datum
    type(flow_condition_dataset_type) :: gradient
    type(flow_condition_dataset_type) :: flow_dataset
  end type flow_sub_condition_type
  
  type, public :: sub_condition_ptr_type
    type(flow_sub_condition_type), pointer :: ptr
  end type sub_condition_ptr_type
    
  type, public :: condition_ptr_type
    type(flow_condition_type), pointer :: ptr
  end type condition_ptr_type
  
  type, public :: condition_list_type
    PetscInt :: num_conditions
    type(flow_condition_type), pointer :: first
    type(flow_condition_type), pointer :: last
    type(flow_condition_type), pointer :: array(:)    
  end type condition_list_type
  
  type, public :: tran_condition_type
    PetscInt :: id                     ! id from which condition can be referenced
    PetscInt :: itype                  ! integer describing type of condition
    PetscBool :: is_transient
    character(len=MAXWORDLENGTH) :: name  ! name of condition (e.g. initial, recharge)
    type(tran_constraint_coupler_type), pointer :: constraint_coupler_list
    type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
    type(tran_condition_type), pointer :: next
  end type tran_condition_type
  
  type, public :: tran_condition_ptr_type
    type(tran_condition_type), pointer :: ptr
  end type tran_condition_ptr_type
  
  type, public :: tran_condition_list_type
    PetscInt :: num_conditions
    type(tran_condition_type), pointer :: first
    type(tran_condition_type), pointer :: last
    type(tran_condition_ptr_type), pointer :: array(:)    
  end type tran_condition_list_type
  
  type, public :: tran_constraint_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name         
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(mineral_constraint_type), pointer :: minerals
    type(srfcplx_constraint_type), pointer :: surface_complexes
    type(colloid_constraint_type), pointer :: colloids
    PetscBool :: requires_equilibration
    type(tran_constraint_type), pointer :: next    
  end type tran_constraint_type
  
  type, public :: tran_constraint_ptr_type
    type(tran_constraint_type), pointer :: ptr
  end type tran_constraint_ptr_type
  
  type, public :: tran_constraint_list_type
    PetscInt :: num_constraints
    type(tran_constraint_type), pointer :: first
    type(tran_constraint_type), pointer :: last
    type(tran_constraint_ptr_type), pointer :: array(:)    
  end type tran_constraint_list_type
  
  type, public :: tran_constraint_coupler_type
    character(len=MAXWORDLENGTH) :: constraint_name   
    PetscReal :: time
    PetscInt :: num_iterations
    character(len=MAXWORDLENGTH) :: time_units
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(mineral_constraint_type), pointer :: minerals
    type(srfcplx_constraint_type), pointer :: surface_complexes
    type(colloid_constraint_type), pointer :: colloids
    type(global_auxvar_type), pointer :: global_auxvar
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar
    type(tran_constraint_coupler_type), pointer :: next   
  end type tran_constraint_coupler_type
      
  public :: FlowConditionCreate, FlowConditionDestroy, FlowConditionRead, &
            FlowConditionGeneralRead, &
            FlowConditionAddToList, FlowConditionInitList, FlowConditionDestroyList, &
            FlowConditionGetPtrFromList, FlowConditionUpdate, &
            FlowConditionPrint, &
            TranConditionCreate, TranConstraintCreate, &
            TranConditionAddToList, TranConditionInitList, &
            TranConditionDestroyList, TranConditionGetPtrFromList, &
            TranConstraintAddToList, TranConstraintInitList, &
            TranConstraintDestroyList, TranConstraintGetPtrFromList, &
            TranConditionRead, TranConstraintRead, &
            TranConditionUpdate, &
            FlowConditionIsTransient, FlowConditionDatasetGetTimes
    
contains

! ************************************************************************** !
!
! FlowConditionCreate: Creates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function FlowConditionCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_condition_type), pointer :: FlowConditionCreate
  
  type(flow_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%saturation)
  nullify(condition%rate)
  nullify(condition%well)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)
  nullify(condition%sub_condition_ptr)
  nullify(condition%general)
  nullify(condition%itype)
  nullify(condition%next)
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%iphase = 0
  condition%num_sub_conditions = 0
  condition%name = ''
  
  FlowConditionCreate => condition

end function FlowConditionCreate

! ************************************************************************** !
!
! TranConditionCreate: Creates a transport condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function TranConditionCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_condition_type), pointer :: TranConditionCreate
  
  type(tran_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%constraint_coupler_list)
  nullify(condition%cur_constraint_coupler)
  nullify(condition%next)
  condition%id = 0
  condition%itype = 0
  condition%name = ''

  TranConditionCreate => condition

end function TranConditionCreate

! ************************************************************************** !
!
! TranConstraintCreate: Creates a transport constraint (set of concentrations
!                       and constraints for setting boundary or initial 
!                       condition).
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function TranConstraintCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_constraint_type), pointer :: TranConstraintCreate
  
  type(tran_constraint_type), pointer :: constraint
  
  allocate(constraint)
  nullify(constraint%aqueous_species)
  nullify(constraint%minerals)
  nullify(constraint%surface_complexes)
  nullify(constraint%colloids)
  nullify(constraint%next)
  constraint%id = 0
  constraint%name = ''
  constraint%requires_equilibration = PETSC_FALSE
  
  TranConstraintCreate => constraint

end function TranConstraintCreate

! ************************************************************************** !
!
! TranConstraintCouplerCreate: Creates a coupler that ties a constraint to a
!                              transport condition
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
function TranConstraintCouplerCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type), pointer :: TranConstraintCouplerCreate
  
  type(tran_constraint_coupler_type), pointer :: coupler
  
  allocate(coupler)
  nullify(coupler%aqueous_species)
  nullify(coupler%minerals)
  nullify(coupler%surface_complexes)
  nullify(coupler%colloids)
  
  coupler%num_iterations = 0
  nullify(coupler%rt_auxvar)
  nullify(coupler%global_auxvar)
  
  nullify(coupler%next)
  coupler%constraint_name = ''
  coupler%time = 0.d0
  coupler%time_units = ''
  
  TranConstraintCouplerCreate => coupler

end function TranConstraintCouplerCreate

! ************************************************************************** !
!
! FlowGeneralConditionCreate: Creates a condition for general mode
! author: Glenn Hammond
! date: 05/26/11
!
! ************************************************************************** !
function FlowGeneralConditionCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_general_condition_type), pointer :: FlowGeneralConditionCreate
  
  type(flow_general_condition_type), pointer :: general_condition
  
  allocate(general_condition)
  nullify(general_condition%liquid_pressure)
  nullify(general_condition%gas_pressure)
  nullify(general_condition%gas_saturation)
  nullify(general_condition%mole_fraction)
  nullify(general_condition%temperature)
  nullify(general_condition%flux)

  FlowGeneralConditionCreate => general_condition

end function FlowGeneralConditionCreate

! ************************************************************************** !
!
! FlowGeneralSubConditionPtr: Returns a pointer to a subcondition, creating
!                             them if necessary
! author: Glenn Hammond
! date: 06/09/11
!
! ************************************************************************** !
function FlowGeneralSubConditionPtr(sub_condition_name,general, &
                                    option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_general_condition_type) :: general
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowGeneralSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('LIQUID_PRESSURE')
      if (associated(general%liquid_pressure)) then
        sub_condition_ptr => general%liquid_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%liquid_pressure => sub_condition_ptr
      endif
    case('GAS_PRESSURE')
      if (associated(general%gas_pressure)) then
        sub_condition_ptr => general%gas_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','GAS_SATURATION')
      if (associated(general%gas_saturation)) then
        sub_condition_ptr => general%gas_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_saturation => sub_condition_ptr
      endif
    case('TEMPERATURE')
      if (associated(general%temperature)) then
        sub_condition_ptr => general%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%temperature => sub_condition_ptr
      endif
    case('MOLE_FRACTION')
      if (associated(general%mole_fraction)) then
        sub_condition_ptr => general%mole_fraction
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%mole_fraction => sub_condition_ptr
      endif
    case('FLUX')
      if (associated(general%flux)) then
        sub_condition_ptr => general%flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%flux => sub_condition_ptr
      endif
    case default
      option%io_buffer = 'keyword (' // trim(sub_condition_name) // &
                          ') not recognized in general condition,type'
      call printErrMsg(option)
  end select

  FlowGeneralSubConditionPtr => sub_condition_ptr

end function FlowGeneralSubConditionPtr

! ************************************************************************** !
!
! FlowSubConditionCreate: Creates a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
function FlowSubConditionCreate(ndof)

  use Option_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: FlowSubConditionCreate
  
  PetscInt :: ndof
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = 0
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''

  ! by default, create time series
  call FlowConditionDatasetInit(sub_condition%flow_dataset)
  sub_condition%flow_dataset%time_series => TimeSeriesCreate()
  sub_condition%flow_dataset%time_series%rank = ndof
  call FlowConditionDatasetInit(sub_condition%gradient)
  sub_condition%gradient%time_series => TimeSeriesCreate()
  sub_condition%gradient%time_series%rank = 3
  call FlowConditionDatasetInit(sub_condition%datum)
  sub_condition%datum%time_series => TimeSeriesCreate()
  sub_condition%datum%time_series%rank = 3

  FlowSubConditionCreate => sub_condition

end function FlowSubConditionCreate

! ************************************************************************** !
!
! GetFlowSubCondFromArrayByName: returns a pointer to a subcondition with
!                                 matching name
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetFlowSubCondFromArrayByName(sub_condition_ptr_list,name)

  use Input_module
  use String_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: GetFlowSubCondFromArrayByName
  type(sub_condition_ptr_type), pointer :: sub_condition_ptr_list(:)
  character(len=MAXWORDLENGTH) :: name
  
  PetscInt :: idof
  PetscInt :: length
  
  nullify(GetFlowSubCondFromArrayByName)
  length = len_trim(name)
  do idof = 1, size(sub_condition_ptr_list)
    if (length == len_trim(sub_condition_ptr_list(idof)%ptr%name) .and. &
        StringCompare(name,sub_condition_ptr_list(idof)%ptr%name,length)) then
      GetFlowSubCondFromArrayByName => sub_condition_ptr_list(idof)%ptr
      return
    endif
  enddo
  
end function GetFlowSubCondFromArrayByName
          
! ************************************************************************** !
!
! FlowConditionDatasetInit: Initializes a dataset
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine FlowConditionDatasetInit(flow_condition_dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: flow_condition_dataset

  nullify(flow_condition_dataset%time_series)
  nullify(flow_condition_dataset%dataset)
    
end subroutine FlowConditionDatasetInit

! ************************************************************************** !
!
! FlowSubConditionVerify: Verifies the data in a subcondition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine FlowSubConditionVerify(option, condition, sub_condition_name, &
                                  sub_condition, &
                                  default_time, &
                                  default_ctype, default_itype, &
                                  default_flow_dataset, &
                                  default_datum, default_gradient, &
                                  destroy_if_null)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(flow_condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_sub_condition_type), pointer :: sub_condition
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscBool :: default_cyclic
  PetscInt :: default_interpolation
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_flow_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  PetscBool :: destroy_if_null

  PetscInt :: array_size

  if (.not.associated(sub_condition)) return
  
  if (.not. (associated(sub_condition%flow_dataset%time_series%values) .or. &
             associated(sub_condition%flow_dataset%dataset))) then
    if (destroy_if_null) call FlowSubConditionDestroy(sub_condition)
    return
  endif
  
  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call printErrMsg(option)
  endif
  
  call FlowConditionDatasetVerify(option,condition%name,sub_condition_name, &
                                  default_time,sub_condition%flow_dataset, &
                                  default_flow_dataset)
  call FlowConditionDatasetVerify(option,condition%name,sub_condition_name, &
                                  default_time,sub_condition%datum,&
                                  default_datum)
  call FlowConditionDatasetVerify(option,condition%name,sub_condition_name, &
                                  default_time,sub_condition%gradient, &
                                  default_gradient)

end subroutine FlowSubConditionVerify

! ************************************************************************** !
!
! FlowConditionDatasetVerify: Verifies the data in a dataset
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine FlowConditionDatasetVerify(option, condition_name, &
                                      sub_condition_name, &
                                      default_time, &
                                      dataset, default_dataset)
  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: condition_name
  character(len=MAXWORDLENGTH) :: sub_condition_name
  character(len=MAXWORDLENGTH) :: size1, size2
  PetscReal :: default_time
  type(flow_condition_dataset_type) :: dataset
  type(flow_condition_dataset_type) :: default_dataset
  
  if (associated(dataset%dataset) .or. &
      associated(default_dataset%dataset)) then
    call TimeSeriesDestroy(dataset%time_series)
    call TimeSeriesDestroy(default_dataset%time_series)
    if (associated(default_dataset%dataset) .and. &
        .not.associated(dataset%dataset)) then
      dataset%dataset => default_dataset%dataset
    endif
    if (associated(dataset%dataset) .and. &
        .not.associated(default_dataset%dataset)) then
      default_dataset%dataset => dataset%dataset
    endif
  endif
  
  if (associated(dataset%time_series)) then
    call TimeSeriesVerify(option, default_time, dataset%time_series, &
                          default_dataset%time_series)
  endif
  
end subroutine FlowConditionDatasetVerify

! ************************************************************************** !
!
! FlowConditionDatasetGetTimes: Fills an array of times based on dataset
! author: Glenn Hammond
! date: 05/19/11
!
! ************************************************************************** !
subroutine FlowConditionDatasetGetTimes(option, sub_condition, &
                                        max_sim_time, times)
  use Option_module

  implicit none
  
  type(option_type) :: option
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscReal :: max_sim_time
  PetscReal, pointer :: times(:)
  
  type(flow_condition_dataset_type), pointer :: flow_dataset
  
  flow_dataset => sub_condition%flow_dataset

  if (associated(flow_dataset%time_series) .and. &
      associated(flow_dataset%dataset)) then
    option%io_buffer = 'FlowConditionDatasetGetTimes() currently does not ' // &
                       'support both time_series and datasets.'
    call printErrMsg(option)
  endif
  if (associated(flow_dataset%time_series)) then
    call TimeSeriesGetTimes(option, flow_dataset%time_series, max_sim_time, &
                            times)
  endif
  if (associated(flow_dataset%dataset)) then
    call DatasetGetTimes(option, flow_dataset%dataset, max_sim_time, &
                            times)
  endif
 
end subroutine FlowConditionDatasetGetTimes

! ************************************************************************** !
!
! FlowConditionRead: Reads a condition from the input file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine FlowConditionRead(condition,input,option)

  use Option_module
  use Input_module
  use String_module
  use Logging_module  
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(flow_sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, rate, well,&
                                       sub_condition_ptr, saturation, &
                                       displacement_x, displacement_y, &
                                       displacement_z
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_flow_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  type(flow_condition_dataset_type) :: default_well
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, idof
  PetscBool :: found
  PetscBool :: destroy_if_null
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read,ierr)

  default_time = 0.d0
  default_iphase = 0
  call FlowConditionDatasetInit(default_flow_dataset)
  default_flow_dataset%time_series => TimeSeriesCreate()
  default_flow_dataset%time_series%rank = 1
  default_flow_dataset%time_series%interpolation_method = STEP
  default_flow_dataset%time_series%is_cyclic = PETSC_FALSE
  call FlowConditionDatasetInit(default_datum)
  default_datum%time_series => TimeSeriesCreate()
  default_datum%time_series%rank = 3
  call FlowConditionDatasetInit(default_gradient)
  default_gradient%time_series => TimeSeriesCreate()
  default_gradient%time_series%rank = 3
  call FlowConditionDatasetInit(default_well)
  default_well%time_series => TimeSeriesCreate()
  default_well%time_series%rank = 7 + option%nflowspec

  pressure => FlowSubConditionCreate(option%nphase)
  pressure%name = 'pressure'
  flux => pressure
  rate => FlowSubConditionCreate(option%nflowspec)
  rate%name = 'rate'
  well => FlowSubConditionCreate(default_well%time_series%rank)
  well%name = 'well'
  saturation => FlowSubConditionCreate(option%nphase)
  saturation%name = 'saturation'
  temperature => FlowSubConditionCreate(ONE_INTEGER)
  temperature%name = 'temperature'
  concentration => FlowSubConditionCreate(ONE_INTEGER)
  concentration%name = 'concentration'
  enthalpy => FlowSubConditionCreate(option%nphase)
  enthalpy%name = 'enthalpy'
  displacement_x => FlowSubConditionCreate(ONE_INTEGER)
  displacement_y => FlowSubConditionCreate(ONE_INTEGER)
  displacement_z => FlowSubConditionCreate(ONE_INTEGER)
  displacement_x%name = 'displacement_x'
  displacement_y%name = 'displacement_y'
  displacement_z%name = 'displacement_z'


  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  rate%units = 'kg/s'
  well%units = 'Pa'
  saturation%units = ' '
  temperature%units = 'C'
  concentration%units = 'M'
  enthalpy%units = 'KJ/mol'
  displacement_x%units = 'm'
  displacement_y%units = 'm'
  displacement_z%units = 'm'

  
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
            case('Pa','KPa')
              pressure%units = trim(word)
            case('kg/s','kg/yr')
              rate%units = trim(word)
            case('m/s','m/yr')
              flux%units = trim(word)
            case('C','K')
              temperature%units = trim(word)
            case('M','mol/L')
              concentration%units = trim(word)
            case('KJ/mol')
              enthalpy%units = trim(word)
          end select
        enddo
      case('CYCLIC')
        default_flow_dataset%time_series%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_flow_dataset%time_series%interpolation_method = STEP
          case('linear') 
            default_flow_dataset%time_series%interpolation_method = LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(trim(word))
            case('PRESSURE')
              sub_condition_ptr => pressure
            case('RATE')
              sub_condition_ptr => rate
            case('WELL')
              sub_condition_ptr => well
            case('FLUX')
              sub_condition_ptr => flux
            case('SATURATION')
              sub_condition_ptr => saturation
            case('TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONCENTRATION')
              sub_condition_ptr => concentration
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
            case('DISPLACEMENT_X')
              sub_condition_ptr => displacement_x
            case('DISPLACEMENT_Y')
              sub_condition_ptr => displacement_y
            case('DISPLACEMENT_Z')
              sub_condition_ptr => displacement_z
            case default
              option%io_buffer = 'keyword (' // trim(word) // &
                                 ') not recognized in condition,type'
              call printErrMsg(option)
          end select
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'TYPE','CONDITION')   
          call StringToLower(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
            case('scaled_mass_rate')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    option%io_buffer = 'scaled_mass_rate type: ' // &
                      trim(word) // &
                      'not recognized for flow condition "' // &
                      trim(condition%name) // '".'
                    call printErrMsg(option)
                end select
              else
                sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
              endif
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('well','production_well', 'injection_well')
              sub_condition_ptr%itype = WELL_SS
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
            case('scaled_volumetric_rate')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    option%io_buffer = 'scaled_volumetric_rate type: ' // &
                      trim(word) // &
                      'not recognized for flow condition "' // &
                      trim(condition%name) // '".'
                    call printErrMsg(option)
                end select
              else
                sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
              endif
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('unit_gradient')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may ' // &
                  'only be associated with a PRESSURE flow condition.'
                call printErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC
            case('distributed_volumetric_rate')
              sub_condition_ptr%itype = DISTRIBUTED_VOLUMETRIC_RATE_SS
            case('distributed_mass_rate')
              sub_condition_ptr%itype = DISTRIBUTED_MASS_RATE_SS
            case default
              option%io_buffer = 'bc type "' // trim(word) // &
                                 '" not recognized in condition,type'
              call printErrMsg(option)
          end select
        enddo
      case('TIME','TIMES')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION')   
      case('IPHASE')
        call InputReadInt(input,option,default_iphase)
        call InputErrorMsg(input,option,'IPHASE','CONDITION')   
      case('DATUM')
        call FlowConditionReadValues(input,option,word,string,default_datum,word)
      case('GRADIENT','GRAD')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
            case('RATE')
              sub_condition_ptr => rate
            case('WELL')
              sub_condition_ptr => well
            case('FLUX')
              sub_condition_ptr => flux
            case('SATURATION')
              sub_condition_ptr => saturation
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONC','CONCENTRATION')
              sub_condition_ptr => concentration
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
            case default
              option%io_buffer = 'keyword not recognized in condition,type'
              call printErrMsg(option)
          end select
          call FlowConditionReadValues(input,option,word,string, &
                                       sub_condition_ptr%gradient,word)
          nullify(sub_condition_ptr)
        enddo
      case('TEMPERATURE','TEMP')
        call FlowConditionReadValues(input,option,word,string, &
                                     temperature%flow_dataset, &
                                     temperature%units)
      case('ENTHALPY','H')
        call FlowConditionReadValues(input,option,word,string, &
                                     enthalpy%flow_dataset, &
                                     enthalpy%units)
      case('PRESSURE','PRES','PRESS')
        call FlowConditionReadValues(input,option,word,string, &
                                     pressure%flow_dataset, &
                                     pressure%units)
      case('RATE')
        call FlowConditionReadValues(input,option,word,string, &
                                     rate%flow_dataset, &
                                     rate%units)
      case('WELL')
        call FlowConditionReadValues(input,option,word,string, &
                                     well%flow_dataset, &
                                     well%units)
      case('FLUX','VELOCITY','VEL')
        call FlowConditionReadValues(input,option,word,string, &
                                     pressure%flow_dataset, &
                                     pressure%units)
      case('CONC','CONCENTRATION')
        call FlowConditionReadValues(input,option,word,string, &
                                     concentration%flow_dataset, &
                                     concentration%units)
      case('SAT','SATURATION')
        call FlowConditionReadValues(input,option,word,string, &
                                     saturation%flow_dataset, &
                                     saturation%units)
      case('DISPLACEMENT_X')
        call FlowConditionReadValues(input,option,word,string, &
                                     displacement_x%flow_dataset, &
                                     displacement_x%units)
      case('DISPLACEMENT_Y')
        call FlowConditionReadValues(input,option,word,string, &
                                     displacement_y%flow_dataset, &
                                     displacement_y%units) 
      case('DISPLACEMENT_Z')
        call FlowConditionReadValues(input,option,word,string, &
                                     displacement_z%flow_dataset, &
                                     displacement_z%units)
      case('CONDUCTANCE')
        call InputReadDouble(input,option,pressure%flow_dataset%time_series%lame_aux_variable_remove_me)
        call InputErrorMsg(input,option,'CONDUCTANCE','CONDITION')   
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized in flow condition'
        call printErrMsg(option)                                 
    end select 
  
  enddo  
  
  ! check whether
  if (default_iphase == 0) then
    option%io_buffer = '"iphase" not set in condition; set to 1'
    call printWrnMsg(option)
    condition%iphase = 1
  else
    condition%iphase = default_iphase    
  endif
  
  ! update datum and gradient defaults, if null, based on flow_dataset default
  if (default_flow_dataset%time_series%is_cyclic) then
    default_datum%time_series%is_cyclic = PETSC_TRUE
    default_gradient%time_series%is_cyclic = PETSC_TRUE
  endif
  if (default_datum%time_series%interpolation_method == NULL) &
    default_datum%time_series%interpolation_method = &
    default_flow_dataset%time_series%interpolation_method
  if (default_gradient%time_series%interpolation_method == NULL) &
    default_gradient%time_series%interpolation_method = &
    default_flow_dataset%time_series%interpolation_method

  ! check to ensure that a rate condition is not of type pressure   
  if (associated(rate)) then
    select case(rate%itype)
      case(DIRICHLET_BC,NEUMANN_BC,HYDROSTATIC_BC,UNIT_GRADIENT_BC, &
           CONDUCTANCE_BC,ZERO_GRADIENT_BC,SEEPAGE_BC)
        option%io_buffer = 'RATE condition must not be of type: dirichlet, ' // &
          'neumann, zero_gradient, dirichlet_zero_gradient, hydrostatic, ' // &
          'seepage, or conductance".'
        call printErrMsg(option)
    end select
  endif
  ! check to ensure that a pressure condition is not of type rate   
  if (associated(pressure)) then                          
    select case(pressure%itype)
      case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS)
        option%io_buffer = 'PRESSURE or FLUX condition must not be of type: ' // &
          'mass_rate, scaled_mass_rate, volumetric_rate, ' // &
          'scaled_volumetric_rate, equilibrium, or production_well.'
        call printErrMsg(option)
    end select
  endif

  ! verify the datasets
  word = 'pressure/flux'
  call FlowSubConditionVerify(option,condition,word,pressure,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient, PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'well'
  call FlowSubConditionVerify(option,condition,word,well,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'saturation'
  call FlowSubConditionVerify(option,condition,word,saturation,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)

  word = 'concentration'
  call FlowSubConditionVerify(option,condition,word,concentration,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'displacement_x'
  call FlowSubConditionVerify(option,condition,word,displacement_x,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'displacement_y'
  call FlowSubConditionVerify(option,condition,word,displacement_y,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'displacement_z'
  call FlowSubConditionVerify(option,condition,word,displacement_z,default_time, &
                              default_ctype, default_itype, &
                              default_flow_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)

  select case(option%iflowmode)
    case(G_MODE)
      option%io_buffer = 'General mode not supported in original FlowConditionRead.'
      call printMsg(option)
    case(THC_MODE,MPH_MODE,IMS_MODE,FLASH2_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
     
      
      if (.not.associated(temperature)) then
        option%io_buffer = 'temperature condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%temperature => temperature
      
      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%concentration => concentration
      
      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%enthalpy => enthalpy
      
      condition%num_sub_conditions = 4
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 4
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
      if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
        
      allocate(condition%itype(FIVE_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      condition%itype(TWO_INTEGER) = temperature%itype
      condition%itype(THREE_INTEGER) = concentration%itype
      if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype

!#if 0      
    case(MIS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well)) then
        option%io_buffer = 'pressure and rate condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      
      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%concentration => concentration

#if 0
      if (.not.associated(temperature)) then
        option%io_buffer = 'temperature condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%temperature => temperature
      
      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%enthalpy => enthalpy
#endif

      condition%num_sub_conditions = 2
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 2
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs in problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
!     condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => concentration
!     if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy

      allocate(condition%itype(TWO_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      condition%itype(TWO_INTEGER) = concentration%itype
!#endif
    
    case(RICHARDS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate) .and. &
          .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)      
      endif
      
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      
      condition%num_sub_conditions = 1
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      if (associated(pressure)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      elseif (associated(saturation)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => saturation
      elseif (associated(rate)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      endif                         

      allocate(condition%itype(ONE_INTEGER))
      if (associated(pressure)) then 
        condition%itype(ONE_INTEGER) = pressure%itype
      else if (associated(saturation)) then
        condition%itype(ONE_INTEGER) = saturation%itype
      else if (associated(rate)) then
        condition%itype(ONE_INTEGER) = rate%itype
      endif
      
      ! these are not used with richards
      if (associated(temperature)) call FlowSubConditionDestroy(temperature)
      if (associated(enthalpy)) call FlowSubConditionDestroy(enthalpy)

   case(THMC_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
     
      
      if (.not.associated(temperature)) then
        option%io_buffer = 'temperature condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%temperature => temperature
      
      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%concentration => concentration
      
!       if (.not.associated(enthalpy)) then
!         option%io_buffer = 'enthalpy condition null in condition: ' // &
!                             trim(condition%name)      
!         call printErrMsg(option)
!       endif                         
!       condition%enthalpy => enthalpy
      
      if (.not.associated(displacement_x)) then
        option%io_buffer = 'displacement_x condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%displacement_x => displacement_x

      if (.not.associated(displacement_y)) then
        option%io_buffer = 'displacement_y condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%displacement_y => displacement_y

      if (.not.associated(displacement_z)) then
        option%io_buffer = 'displacement_z condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%displacement_z => displacement_z

      condition%num_sub_conditions = 6
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 6
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
!      if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
      condition%sub_condition_ptr(FOUR_INTEGER)%ptr => displacement_x
      condition%sub_condition_ptr(FIVE_INTEGER)%ptr => displacement_y
      condition%sub_condition_ptr(SIX_INTEGER)%ptr => displacement_z

              
      allocate(condition%itype(SIX_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      condition%itype(TWO_INTEGER) = temperature%itype
      condition%itype(THREE_INTEGER) = concentration%itype
!      if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = enthalpy%itype
      condition%itype(FOUR_INTEGER) = displacement_x%itype
      condition%itype(FIVE_INTEGER) = displacement_y%itype
      condition%itype(SIX_INTEGER) = displacement_z%itype

      
  end select
  
  call FlowConditionDatasetDestroy(default_flow_dataset)
  call FlowConditionDatasetDestroy(default_datum)
  call FlowConditionDatasetDestroy(default_gradient)
  call FlowConditionDatasetDestroy(default_well)
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr)

end subroutine FlowConditionRead

! ************************************************************************** !
!
! FlowConditionGeneralRead: Reads a condition from the input file for
!                           general mode
! author: Glenn Hammond
! date: 09/14/11
!
! ************************************************************************** !
subroutine FlowConditionGeneralRead(condition,input,option)

  use Option_module
  use Input_module
  use String_module
  use Logging_module  
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(flow_general_condition_type), pointer :: general
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_flow_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, idof
  PetscInt :: i
  PetscBool :: found
  PetscBool :: destroy_if_null
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read,ierr)


  default_time = 0.d0
  default_iphase = 0
  call FlowConditionDatasetInit(default_flow_dataset)
  default_flow_dataset%time_series => TimeSeriesCreate()
  default_flow_dataset%time_series%rank = 1
  default_flow_dataset%time_series%interpolation_method = STEP
  default_flow_dataset%time_series%is_cyclic = PETSC_FALSE
  call FlowConditionDatasetInit(default_datum)
  default_datum%time_series => TimeSeriesCreate()
  default_datum%time_series%rank = 3
  call FlowConditionDatasetInit(default_gradient)
  default_gradient%time_series => TimeSeriesCreate()
  default_gradient%time_series%rank = 3
  
  select case(option%iflowmode)
    case(G_MODE)
      general => FlowGeneralConditionCreate(option)
      condition%general => general
  end select
  
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('CYCLIC')
        default_flow_dataset%time_series%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_flow_dataset%time_series%interpolation_method = STEP
          case('linear') 
            default_flow_dataset%time_series%interpolation_method = LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(option%iflowmode)
            case(G_MODE)
              sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                              option)
          end select
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'TYPE','CONDITION')   
          call StringToLower(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
            case('scaled_mass_rate')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
            case('scaled_volumetric_rate')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
            case('distributed_volumetric_rate')
              sub_condition_ptr%itype = DISTRIBUTED_VOLUMETRIC_RATE_SS
            case('distributed_mass_rate')
              sub_condition_ptr%itype = DISTRIBUTED_MASS_RATE_SS
            case default
              option%io_buffer = 'bc type "' // trim(word) // &
                                 '" not recognized in condition,type'
              call printErrMsg(option)
          end select
        enddo
      case('DATUM')
        call FlowConditionReadValues(input,option,word,string,default_datum,word)
      case('GRADIENT')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(option%iflowmode)
            case(G_MODE)
              sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                              option)
          end select
          call FlowConditionReadValues(input,option,word,string, &
                                       sub_condition_ptr%gradient,word)
          nullify(sub_condition_ptr)
        enddo
      case('LIQUID_PRESSURE','GAS_PRESSURE','LIQUID_SATURATION', &
           'GAS_SATURATION','TEMPERATURE','MOLE_FRACTION','RATE','FLUX')
        select case(option%iflowmode)
          case(G_MODE)
            sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                            option)
        end select
        call FlowConditionReadValues(input,option,word,string, &
                                     sub_condition_ptr%flow_dataset, &
                                     sub_condition_ptr%units)
        select case(word)
          case('LIQUID_SATURATION') ! convert to gas saturation
            do i = 1, size(sub_condition_ptr%flow_dataset%time_series%values)
              sub_condition_ptr%flow_dataset%time_series%values(i,:) = &
                1.d0 - sub_condition_ptr%flow_dataset%time_series%values(i,:)
            enddo
        end select
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized in flow condition'
        call printErrMsg(option)                                 
    end select 
  
  enddo  
  
  ! update datum and gradient defaults, if null, based on flow_dataset default
  if (default_flow_dataset%time_series%is_cyclic) then
    default_datum%time_series%is_cyclic = PETSC_TRUE
    default_gradient%time_series%is_cyclic = PETSC_TRUE
  endif
  if (default_datum%time_series%interpolation_method == NULL) &
    default_datum%time_series%interpolation_method = &
      default_flow_dataset%time_series%interpolation_method
  if (default_gradient%time_series%interpolation_method == NULL) &
    default_gradient%time_series%interpolation_method = &
      default_flow_dataset%time_series%interpolation_method

  ! need mole fraction and some sort of saturation
  if (associated(general%flux)) then
    ! neumann or mass/volumetric flux
    ! need temperature

    condition%sub_condition_ptr(GENERAL_FLUX_DOF)%ptr => general%flux
    if (.not.associated(general%mole_fraction)) then
      option%io_buffer = 'General Phase rate condition must include ' // &
        'mole fraction '
      call printErrMsg(option)
    endif
    if (.not.associated(general%gas_saturation)) then
      option%io_buffer = 'General Phase rate condition must include ' // &
        'gas or liquid saturation'
      call printErrMsg(option)
    endif
    if (.not.associated(general%temperature)) then
      option%io_buffer = 'General Phase rate condition must include ' // &
        'temperature'
      call printErrMsg(option)
    endif
  else
  
!geh    condition%num_sub_conditions = THREE_INTEGER
!geh    allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
!geh    do idof = 1, condition%num_sub_conditions
!geh      nullify(condition%sub_condition_ptr(idof)%ptr)
!geh    enddo

    ! some sort of dirichlet-based pressure, temperature, etc.
    if (.not.associated(general%liquid_pressure) .and. &
        .not.associated(general%gas_pressure)) then
      option%io_buffer = 'General Phase non-rate condition must include ' // &
        'a liquid or gas pressure'
      call printErrMsg(option)
    endif
    if (.not.associated(general%mole_fraction) .and. &
        .not.associated(general%gas_saturation)) then
      option%io_buffer = 'General Phase non-rate condition must include ' // &
        'mole fraction or gas/liquid saturation'
      call printErrMsg(option)
    endif
    if (.not.associated(general%temperature)) then
      option%io_buffer = 'General Phase non-rate condition must include ' // &
        'temperature'
      call printErrMsg(option)
    endif
    if (associated(general%gas_pressure) .and. &
        associated(general%gas_saturation)) then
      ! two phase condition
      condition%iphase = TWO_PHASE_STATE
    else if (associated(general%liquid_pressure) .and. &
              associated(general%mole_fraction)) then
      ! liquid phase condition
      condition%iphase = LIQUID_STATE
    else if (associated(general%gas_pressure) .and. &
              associated(general%mole_fraction)) then
      ! gas phase condition
      condition%iphase = GAS_STATE
    else 
      option%io_buffer = 'General Phase non-rate condition contains an ' // &
        'unsupported combination of primary dependent variables.'
      call printErrMsg(option)
    endif
  endif
    
  ! verify the datasets
  word = 'liquid pressure'
  call FlowSubConditionVerify(option,condition,word,general%liquid_pressure, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)
  word = 'gas pressure'
  call FlowSubConditionVerify(option,condition,word,general%gas_pressure, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)
  word = 'gas saturation'
  call FlowSubConditionVerify(option,condition,word,general%gas_saturation, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)
  word = 'mole fraction'
  call FlowSubConditionVerify(option,condition,word,general%mole_fraction, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,general%temperature, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)
  word = 'flux'
  call FlowSubConditionVerify(option,condition,word,general%flux, &
                              default_time, default_ctype, default_itype, &
                              default_flow_dataset, default_datum, &
                              default_gradient, PETSC_TRUE)

  call FlowConditionDatasetDestroy(default_flow_dataset)
  call FlowConditionDatasetDestroy(default_datum)
  call FlowConditionDatasetDestroy(default_gradient)
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr)

end subroutine FlowConditionGeneralRead

! ************************************************************************** !
!
! TranConditionRead: Reads a transport condition from the input file
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConditionRead(condition,constraint_list,reaction,input,option)

  use Option_module
  use Input_module
  use String_module
  use Logging_module  
  use Units_module
  
  implicit none
  
  type(tran_condition_type) :: condition
  type(tran_constraint_list_type) :: constraint_list
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  type(tran_constraint_type), pointer :: constraint
  type(tran_constraint_coupler_type), pointer :: constraint_coupler, cur_coupler
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: default_time
  character(len=MAXWORDLENGTH) :: default_time_units
  PetscInt :: default_iphase
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscBool :: found
  PetscInt :: icomp
  PetscBool :: minerals_exist
  PetscErrorCode :: ierr
  PetscReal :: conversion

  call PetscLogEventBegin(logging%event_tran_condition_read,ierr)

  default_time = 0.d0
  default_iphase = 0
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC
  default_time_units = ''

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'TYPE','CONDITION')   
        call StringToLower(word)
        select case(word)
            case('dirichlet')
              condition%itype = DIRICHLET_BC
            case('dirichlet_zero_gradient')
              condition%itype = DIRICHLET_ZERO_GRADIENT_BC
            case('equilibrium')
              condition%itype = EQUILIBRIUM_SS
            case('neumann')
              condition%itype = NEUMANN_BC
            case('mole','mole_rate')
              condition%itype = MASS_RATE_SS
            case('zero_gradient')
              condition%itype = ZERO_GRADIENT_BC
            case default
              option%io_buffer = 'Keyword ' // trim(word) // &
                                 ' not recognized in condition,type'
              call printErrMsg(option)
        end select
      case('TIME')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION') 
      case('TIME_UNITS') 
        call InputReadWord(input,option,word,PETSC_TRUE) 
        call InputErrorMsg(input,option,'UNITS','CONDITION')   
        call StringToLower(word)
        select case(trim(word))     
          case('s','sec','min','m','hr','h','d','day','y','yr')
            default_time_units = trim(word)         
        end select          
      case('CONSTRAINT_LIST')
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT')
              
          if (InputCheckExit(input,option)) exit  

          constraint_coupler => TranConstraintCouplerCreate(option)
          call InputReadDouble(input,option,constraint_coupler%time)
          call InputErrorMsg(input,option,'time','CONSTRAINT_LIST') 
          ! time units are optional  
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'constraint name','CONSTRAINT_LIST') 
          ! read constraint name
          call InputReadWord(input,option,constraint_coupler%constraint_name,PETSC_TRUE)
          if (InputError(input)) then
            constraint_coupler%time_units = default_time_units
            constraint_coupler%constraint_name = trim(word)
          else
            constraint_coupler%time_units = word
          endif
          ! convert time units
          if (len_trim(constraint_coupler%time_units) > 0) then
            constraint_coupler%time = constraint_coupler%time* &
              UnitsConvertToInternal(constraint_coupler%time_units,option)
          endif
          ! add to end of list
          if (.not.associated(condition%constraint_coupler_list)) then
            condition%constraint_coupler_list => constraint_coupler
          else
            cur_coupler => condition%constraint_coupler_list
            do
              if (.not.associated(cur_coupler%next)) exit
              cur_coupler => cur_coupler%next
            enddo
            cur_coupler%next => constraint_coupler
          endif
        enddo
      case('CONSTRAINT')
        constraint => TranConstraintCreate(option)
        constraint_coupler => TranConstraintCouplerCreate(option)
        call InputReadWord(input,option,constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name') 
        option%io_buffer = 'Constraint: ' // trim(constraint%name)
        call printMsg(option)
        call TranConstraintRead(constraint,reaction,input,option)
        call TranConstraintAddToList(constraint,constraint_list)
        constraint_coupler%aqueous_species => constraint%aqueous_species
        constraint_coupler%minerals => constraint%minerals
        constraint_coupler%surface_complexes => constraint%surface_complexes
        constraint_coupler%colloids => constraint%colloids
        constraint_coupler%time = default_time
        ! add to end of coupler list
        if (.not.associated(condition%constraint_coupler_list)) then
          condition%constraint_coupler_list => constraint_coupler
        else
          cur_coupler => condition%constraint_coupler_list
          do
            if (.not.associated(cur_coupler%next)) exit
            cur_coupler => cur_coupler%next
          enddo
          cur_coupler%next => constraint_coupler
        endif        
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                 ' not recognized in transport condition'
        call printErrMsg(option)
    end select 
  
  enddo  
  
  if (len_trim(default_time_units) > 0) then
    conversion = UnitsConvertToInternal(default_time_units,option)
    cur_coupler => condition%constraint_coupler_list
    do
      if (.not.associated(cur_coupler)) exit
      if (len_trim(cur_coupler%time_units) == 0) then
        cur_coupler%time = cur_coupler%time*conversion
      endif
      cur_coupler => cur_coupler%next
    enddo
  endif

  call PetscLogEventEnd(logging%event_tran_condition_read,ierr)

end subroutine TranConditionRead

! ************************************************************************** !
!
! TranConstraintRead: Reads a transport constraint from the input file
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintRead(constraint,reaction,input,option)

  use Option_module
  use Input_module
  use String_module
  use Logging_module
 
  implicit none
  
  type(tran_constraint_type) :: constraint
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: icomp, imnrl
  PetscInt :: isrfcplx
  PetscInt :: length
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(srfcplx_constraint_type), pointer :: srfcplx_constraint
  type(colloid_constraint_type), pointer :: colloid_constraint
  PetscErrorCode :: ierr
  PetscReal :: tempreal

  call PetscLogEventBegin(logging%event_tran_constraint_read,ierr)

  ! read the constraint
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONSTRAINT')
        
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONSTRAINT')   
      
    select case(trim(word))

      case('CONC','CONCENTRATIONS')

        aq_species_constraint => &
          AqueousSpeciesConstraintCreate(reaction,option)

        icomp = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'CONSTRAINT, CONCENTRATIONS')
          
          if (InputCheckExit(input,option)) exit  
          
          icomp = icomp + 1        
          
          if (icomp > reaction%naqcomp) then
            option%io_buffer = 'Number of concentration constraints ' // &
                               'exceeds number of primary chemical ' // &
                               'components in constraint: ' // &
                                trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,aq_species_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'aqueous species name', &
                          'CONSTRAINT, CONCENTRATIONS') 
          option%io_buffer = 'Constraint Species: ' // &
                             trim(aq_species_constraint%names(icomp))
          call printMsg(option)
          
          call InputReadDouble(input,option, &
                               aq_species_constraint%constraint_conc(icomp))
          call InputErrorMsg(input,option,'concentration', &
                          'CONSTRAINT, CONCENTRATIONS')          
          
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputDefaultMsg(input,option, &
                            'CONSTRAINT, CONCENTRATION, constraint_type')
          length = len_trim(word)
          if (length > 0) then
            call StringToUpper(word)
            select case(word)
              case('F','FREE')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_FREE
              case('T','TOTAL')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
              case('TOTAL_SORB')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_TOTAL_SORB
              case('S')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_TOTAL_SORB_AQ_BASED
              case('P','PH')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_PH
              case('L','LOG')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_LOG
              case('M','MINERAL','MNRL') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_MINERAL
              case('G','GAS') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_GAS
              case('SC','CONSTRAINT_SUPERCRIT_CO2') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_SUPERCRIT_CO2
              case('Z','CHG') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_CHARGE_BAL
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                         ' not recognized in constraint,concentration'
                call printErrMsg(option)
            end select 
            
            if (aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_MINERAL .or. &
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_GAS .or.&
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_SUPERCRIT_CO2) then
              call InputReadWord(input,option,aq_species_constraint% &
                                 constraint_aux_string(icomp), &
                                 PETSC_TRUE)
              call InputErrorMsg(input,option,'constraint name', &
                              'CONSTRAINT, CONCENTRATIONS') 
            else
              call InputReadWord(input,option,word,PETSC_FALSE)
              if (input%ierr == 0) then
                call StringToUpper(word)
                select case(word)
                  case('DATASET')
                    call InputReadWord(input,option,aq_species_constraint% &
                                       constraint_aux_string(icomp),PETSC_TRUE)
                    call InputErrorMsg(input,option,'dataset name', &
                                    'CONSTRAINT, CONCENTRATIONS,')
                    aq_species_constraint%external_dataset(icomp) = PETSC_TRUE
                end select
              endif
            endif
          else
            aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
          endif  
        
        enddo  
        
        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of primary species in aqueous constraint.'
          call printErrMsg(option)        
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of primary species in aqueous constraint.'
          call printWrnMsg(option)        
        endif
        
        if (associated(constraint%aqueous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint 
        
      case('MNRL','MINERALS')

        mineral_constraint => MineralConstraintCreate(reaction%mineral,option)

        imnrl = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT, MINERALS')
          
          if (InputCheckExit(input,option)) exit          
          
          imnrl = imnrl + 1

          if (imnrl > reaction%mineral%nkinmnrl) then
            option%io_buffer = &
                     'Number of mineral constraints exceeds number of ' // &
                     'kinetic minerals in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,mineral_constraint%names(imnrl), &
                             PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral name', &
                             'CONSTRAINT, MINERALS')  
          option%io_buffer = 'Constraint Minerals: ' // &
                             trim(mineral_constraint%names(imnrl))
          call printMsg(option)

          ! volume fraction
          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            input%buf = trim(string)
            call InputReadWord(input,option,mineral_constraint% &
                                constraint_aux_string(imnrl),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                            'CONSTRAINT, MINERALS, VOL FRAC')
            mineral_constraint%external_dataset(imnrl) = PETSC_TRUE
            ! set vol frac to NaN to catch bugs
            tempreal = -1.d0
            mineral_constraint%constraint_vol_frac(imnrl) = sqrt(tempreal)
          else
            call InputReadDouble(input,option, &
                                 mineral_constraint%constraint_vol_frac(imnrl))
            call InputErrorMsg(input,option,'volume fraction', &
                               'CONSTRAINT, MINERALS')   
          endif

          ! specific surface area
          call InputReadDouble(input,option, &
                               mineral_constraint%constraint_area(imnrl))
          call InputErrorMsg(input,option,'area', &
                          'CONSTRAINT, MINERALS')          
        
        enddo  
        
        if (imnrl < reaction%mineral%nkinmnrl) then
          option%io_buffer = &
                   'Mineral lists in constraints must provide a volume ' // &
                   'fraction and surface area for all kinetic minerals ' // &
                   '(listed under MINERAL_KINETICS card in CHEMISTRY), ' // &
                   'regardless of whether or not they are present (just ' // &
                   'assign a zero volume fraction if not present).'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%minerals)) then
          call MineralConstraintDestroy(constraint%minerals)
        endif
        constraint%minerals => mineral_constraint 
                            
      case('SURFACE_COMPLEXES')
      
        srfcplx_constraint => &
          SurfaceComplexConstraintCreate(reaction%surface_complexation,option)

        isrfcplx = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option, &
                                       'CONSTRAINT, SURFACE_COMPLEXES')
          
          if (InputCheckExit(input,option)) exit          
          
          isrfcplx = isrfcplx + 1

          if (isrfcplx > reaction%surface_complexation%nkinsrfcplx) then
            option%io_buffer = &
                     'Number of surface complex constraints exceeds ' // &
                     'number of kinetic surface complexes in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,srfcplx_constraint%names(isrfcplx), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'surface complex name', &
                          'CONSTRAINT, SURFACE COMPLEX')  
          option%io_buffer = 'Constraint Surface Complex: ' // &
                             trim(srfcplx_constraint%names(isrfcplx))
          call printMsg(option)
          call InputReadDouble(input,option, &
                               srfcplx_constraint%constraint_conc(isrfcplx))
          call InputErrorMsg(input,option,'concentration', &
                          'CONSTRAINT, SURFACE COMPLEX')          
        enddo  
        
        if (isrfcplx < reaction%surface_complexation%nkinsrfcplx) then
          option%io_buffer = &
                   'Number of surface complex constraints is less than ' // &
                   'number of kinetic surface complexes in surface ' // &
                   'complex constraint.'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%surface_complexes)) then
          call SurfaceComplexConstraintDestroy(constraint%surface_complexes)
        endif
        constraint%surface_complexes => srfcplx_constraint
         
      case('COLL','COLLOIDS')

        colloid_constraint => ColloidConstraintCreate(reaction,option)

        icomp = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT, COLLOIDS')
          
          if (InputCheckExit(input,option)) exit          
          
          icomp = icomp + 1

          if (icomp > reaction%ncoll) then
            option%io_buffer = &
                     'Number of colloid constraints exceeds number of ' // &
                     'colloids in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,colloid_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'colloid name', &
                          'CONSTRAINT, COLLOIDS')  
          option%io_buffer = 'Constraint Colloids: ' // &
                             trim(colloid_constraint%names(icomp))
          call printMsg(option)
          call InputReadDouble(input,option, &
                               colloid_constraint%constraint_conc_mob(icomp))
          call InputErrorMsg(input,option,'mobile concentration', &
                          'CONSTRAINT, COLLOIDS')          
          call InputReadDouble(input,option, &
                               colloid_constraint%constraint_conc_imb(icomp))
          call InputErrorMsg(input,option,'immobile concentration', &
                          'CONSTRAINT, COLLOIDS')          
        
        enddo  
        
        if (icomp < reaction%ncoll) then
          option%io_buffer = &
                   'Colloid lists in constraints must provide mobile ' // &
                   'and immobile concentrations for all colloids ' // &
                   '(listed under the COLLOIDS card in CHEMISTRY), ' // &
                   'regardless of whether or not they are present (just ' // &
                   'assign a small value (e.g. 1.d-40) if not present).'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%colloids)) then
          call ColloidConstraintDestroy(constraint%colloids)
        endif
        constraint%colloids => colloid_constraint 
                                         
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                 ' not recognized in transport constraint'
        call printErrMsg(option)
    end select 
  
  enddo  
  
  call PetscLogEventEnd(logging%event_tran_constraint_read,ierr)

end subroutine TranConstraintRead

! ************************************************************************** !
!
! FlowConditionReadValues: Read the value(s) of a condition variable
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine FlowConditionReadValues(input,option,keyword,string,flow_dataset, &
                                   units)

  use Input_module
  use String_module
  use Option_module
  use Logging_module
  use HDF5_aux_module
  use Units_module
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword
  type(flow_condition_dataset_type) :: flow_dataset
  character(len=MAXWORDLENGTH) :: units
  
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i, icount
  PetscInt :: irank
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_buffer(:)
  type(input_type), pointer :: input2
  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err
#endif

  call PetscLogEventBegin(logging%event_flow_condition_read_values,ierr)    

  nullify(input2)
  filename = ''
  realization_word = ''
  hdf5_path = ''
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToLower(word)
  length = len_trim(word)
  if (length == FOUR_INTEGER .and. StringCompare(word,'file',FOUR_INTEGER)) then 
    input%err_buf2 = trim(keyword) // ', FILE'
    input%err_buf = 'keyword'
    call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
    if (input%ierr == 0) then
      filename = string2
    else
      option%io_buffer = 'The ability to read realization dependent datasets outside the DATASET block is no longer supported'
      call printErrMsg(option)
    endif
    
    if (len_trim(filename) < 2) then
      option%io_buffer = 'No filename listed under Flow_Condition: ' // &
                         trim(keyword)
      call printErrMsg(option)
    endif

    if (index(filename,'.h5') > 0) then
#if !defined(PETSC_HAVE_HDF5)
      write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                               &"read HDF5 formatted flow conditions.")')
      call printErrMsg(option)
#else   
      if (len_trim(hdf5_path) < 1) then
        option%io_buffer = 'No hdf5 path listed under Flow_Condition: ' // &
                           trim(keyword)
        call printErrMsg(option)
      endif

      call h5open_f(hdf5_err)
      option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
      call printMsg(option)
      call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
      call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf5_err,prop_id)
      call h5pclose_f(prop_id,hdf5_err)

      hdf5_path = trim(hdf5_path) // trim(realization_word)
      call HDF5ReadNDimRealArray(option,file_id,hdf5_path,ndims,dims, &
                                 real_buffer)
      option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
      call printMsg(option)  
      call h5fclose_f(file_id,hdf5_err)
      call h5close_f(hdf5_err)
      
      ! dims(1) = size of array
      ! dims(2) = number of data point in time
      if (dims(1)-1 == flow_dataset%time_series%rank) then
        ! alright, the 2d data is layed out in C-style.  now place it in
        ! the appropriate arrays
        allocate(flow_dataset%time_series%times(dims(2)))
        flow_dataset%time_series%times = -999.d0
        allocate(flow_dataset%time_series%values(flow_dataset%time_series%rank,dims(2))) 
        flow_dataset%time_series%values = -999.d0
        icount = 1
        do i = 1, dims(2)
          flow_dataset%time_series%times(i) = real_buffer(icount)
          icount = icount + 1
          do irank = 1, flow_dataset%time_series%rank
            flow_dataset%time_series%values(irank,i) = real_buffer(icount)
            icount = icount + 1
          enddo
        enddo  
      else
        option%io_buffer = 'HDF condition data set rank does not match' // &
          'rank of internal data set.  Email Glenn for additions'
        call printErrMsg(option)
      endif
      if (associated(dims)) deallocate(dims)
      nullify(dims)
      if (associated(real_buffer)) deallocate(real_buffer)
      nullify(real_buffer)
#endif    
    else
      i = index(filename,'.',PETSC_TRUE)
      if (i > 2) then
        filename = filename(1:i-1) // trim(realization_word) // filename(i:)
      else
        filename = trim(filename) // trim(realization_word)
      endif
      input2 => InputCreate(IUNIT_TEMP,filename,option)
      if (flow_dataset%time_series%rank<=3) then
        call FlowConditionReadValuesFromFile(input2,flow_dataset,option)
      else
        call FlowConditionReadValuesFromFile2(input2,flow_dataset,option)
      endif  
      call InputDestroy(input2)
    endif
  else if (StringCompare(word,'dataset')) then
    call InputReadWord(input,option,word,PETSC_TRUE)    
    input%err_buf2 = trim(keyword) // ', DATASET'
    input%err_buf = 'dataset name'
    call InputErrorMsg(input,option)
    flow_dataset%dataset => DatasetCreate()
    flow_dataset%dataset%name = word
  else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then  !sp 
    if (flow_dataset%time_series%rank <= 3) then
      call FlowConditionReadValuesFromFile(input,flow_dataset,option)
    else
      call FlowConditionReadValuesFromFile2(input,flow_dataset,option)
    endif  
  else
    input%buf = trim(string2)
    allocate(flow_dataset%time_series%values(flow_dataset%time_series%rank,1))
    do irank=1,flow_dataset%time_series%rank
      call InputReadDouble(input,option,flow_dataset%time_series%values(irank,1))
      write(input%err_buf,'(a,i2)') trim(keyword) // ' dataset_values, irank = ', irank
      input%err_buf2 = 'CONDITION'
      call InputErrorMsg(input,option) 
    enddo
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) then
      word = trim(keyword) // ' UNITS'
      call InputDefaultMsg(input,option,word)
    else
      units = trim(word)
      flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1) = &
        UnitsConvertToInternal(units,option) * &
        flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1)
    endif
  endif
  call PetscLogEventEnd(logging%event_flow_condition_read_values,ierr)    

end subroutine FlowConditionReadValues

! ************************************************************************** !
!
! FlowConditionReadValuesFromFile: Read values from a external file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine FlowConditionReadValuesFromFile(input,flow_dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module

  implicit none
  
  type(input_type) :: input
  type(flow_condition_dataset_type) :: flow_dataset
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: time_units, data_units
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_times(:), temp_array1(:), temp_array2(:), &
                        temp_array3(:)
  PetscReal :: temp_time
  PetscReal :: conversion
  PetscInt :: max_size
  PetscInt :: temp_max_size
  PetscInt :: count, i, status
  PetscErrorCode :: ierr
  
  time_units = ''
  data_units = ''
  max_size = 1000
  allocate(temp_times(max_size))
  allocate(temp_array1(max_size))
  temp_times = 0.d0
  temp_array1 = 0.d0

  if (flow_dataset%time_series%rank > 1) then
    allocate(temp_array2(max_size))
    temp_array2 = 0.d0
  endif
  
  if (flow_dataset%time_series%rank > 2) then
    allocate(temp_array3(max_size))
    temp_array3 = 0.d0
  endif

  count = 0
  ierr = 0
  do
    call InputReadFlotranString(input,option)
    ! reach the end of file or close out block
    if (InputError(input)) exit  ! check for end of file
    if (InputCheckExit(input,option)) exit  ! check for end of list
    ! check for units on first line
    if (count == 0) then
      string = input%buf
      call InputReadWord(string,word,PETSC_TRUE,ierr)
      call StringToUpper(word)
      select case(word)
        case('TIME_UNITS')
          call InputReadWord(string,time_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'TIME_UNITS','CONDITION FILE')
          call StringToLower(time_units) 
          cycle
        case('DATA_UNITS')
          call InputReadWord(string,data_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'DATA_UNITS','CONDITION FILE')
          call StringToLower(data_units) 
          cycle
      end select
    endif
    count = count + 1
    call InputReadDouble(input,option,temp_times(count))
    call InputErrorMsg(input,option,'time','CONDITION FILE')   
    call InputReadDouble(input,option,temp_array1(count))
    call InputErrorMsg(input,option,'array1','CONDITION FILE')
    if (flow_dataset%time_series%rank > 1) then
      call InputReadDouble(input,option,temp_array2(count))
      call InputErrorMsg(input,option,'array2','CONDITION FILE') 
    endif
    if (flow_dataset%time_series%rank > 2) then
      call InputReadDouble(input,option,temp_array3(count))
      call InputErrorMsg(input,option,'array3','CONDITION FILE') 
    endif
    if (count+1 > max_size) then
      temp_max_size = max_size
      call reallocateRealArray(temp_times,max_size) 
      ! careful.  reallocateRealArray double max_size every time.
      i = temp_max_size
      call reallocateRealArray(temp_array1,i) 
      if (flow_dataset%time_series%rank > 1) then
        i = temp_max_size
        call reallocateRealArray(temp_array2,i)
      endif
      if (flow_dataset%time_series%rank > 2) then
        i = temp_max_size
        call reallocateRealArray(temp_array3,i)
      endif
    endif  
  enddo
  
  if (associated(flow_dataset%time_series%times)) then
    if (count /= size(flow_dataset%time_series%times,1) .and. &
        OptionPrintToScreen(option)) then
      print *, 'Number of times (', count, ') in ', trim(input%filename), &
               ' does not match previous allocation: ', &
               size(flow_dataset%time_series%times,1)
      stop
    endif
    do i=1,count
      if (dabs(flow_dataset%time_series%times(i)-temp_times(i)) > 1.d-8 .and. &
          OptionPrintToScreen(option)) then
        print *, 'Time (', temp_times(i), ') in ', trim(input%filename), &
                 ' does not match previous allocation time: ', &
                 flow_dataset%time_series%times(i), i
        stop
      endif
    enddo
  else
    allocate(flow_dataset%time_series%times(count))
  endif

  if (associated(flow_dataset%time_series%values)) &
    deallocate(flow_dataset%time_series%values)
  allocate(flow_dataset%time_series%values(flow_dataset%time_series%rank,count))

  flow_dataset%time_series%times(1:count) = temp_times(1:count)
  flow_dataset%time_series%values(1,1:count) = temp_array1(1:count)
  if (flow_dataset%time_series%rank > 1) &
    flow_dataset%time_series%values(2,1:count) = temp_array2(1:count)
  if (flow_dataset%time_series%rank > 2) &
    flow_dataset%time_series%values(3,1:count) = temp_array3(1:count)
  
  deallocate(temp_times)
  deallocate(temp_array1)
  if (flow_dataset%time_series%rank > 1) deallocate(temp_array2)
  if (flow_dataset%time_series%rank > 2) deallocate(temp_array3)

  if (len_trim(time_units) > 0) then
    ! Times
    conversion = UnitsConvertToInternal(time_units,option)
    flow_dataset%time_series%times(1:count) = conversion * &
                                   flow_dataset%time_series%times(1:count)
  endif
  if (len_trim(data_units) > 0) then
    ! Data
    conversion = UnitsConvertToInternal(data_units,option)
    flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1:count) = &
      conversion * &
      flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1:count)
  endif
  
end subroutine FlowConditionReadValuesFromFile
! ************************************************************************** !
! 
! FlowConditionReadValuesFromFile: Read values from a external file with 4 more columns
! author: Chuan Lu
! date: 5/31/11
!
! ************************************************************************** !
subroutine FlowConditionReadValuesFromFile2(input,flow_dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module

  implicit none
  
  type(input_type) :: input
  type(flow_condition_dataset_type) :: flow_dataset
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: time_units, data_units
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_times(:), temp_array(:,:)
  PetscReal :: temp_time
  PetscReal :: conversion
  PetscInt :: max_size
  PetscInt :: temp_max_size
  PetscInt :: count, i, status
  PetscErrorCode :: ierr
  
  time_units = ''
  data_units = ''
  max_size = 1000
  allocate(temp_times(max_size))
  allocate(temp_array(flow_dataset%time_series%rank, max_size)) 

  count = 0
  ierr = 0
  do
    call InputReadFlotranString(input,option)
    ! reach the end of file or close out block
    if (InputError(input)) exit  ! check for end of file
    if (InputCheckExit(input,option)) exit  ! check for end of list
    ! check for units on first line
    if (count == 0) then
      string = input%buf
      call InputReadWord(string,word,PETSC_TRUE,ierr)
      call StringToUpper(word)
      select case(word)
        case('TIME_UNITS')
          call InputReadWord(string,time_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'TIME_UNITS','CONDITION FILE')
          call StringToLower(time_units) 
          cycle
        case('DATA_UNITS')
          call InputReadWord(string,data_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'DATA_UNITS','CONDITION FILE')
          call StringToLower(data_units) 
          cycle
      end select
    endif
    count = count + 1
    call InputReadDouble(input,option,temp_times(count))
    call InputErrorMsg(input,option,'time','CONDITION FILE')   


   do i =1, flow_dataset%time_series%rank
      call InputReadDouble(input,option,temp_array(i,count))
      call InputErrorMsg(input,option,'array:',' CONDITION FILE') 
   enddo

 enddo
 if (count+1 > max_size) then
   print *, 'Number of times (', count, ') in ', trim(input%filename), &
          ' exceed 1000 '
   stop
 endif 
  
 if (associated(flow_dataset%time_series%times)) then
   if (count /= size(flow_dataset%time_series%times,1) .and. &
        OptionPrintToScreen(option)) then
      print *, 'Number of times (', count, ') in ', trim(input%filename), &
               ' does not match previous allocation: ', &
               size(flow_dataset%time_series%times,1)
      stop
    endif
    do i=1,count
      if (dabs(flow_dataset%time_series%times(i)-temp_times(i)) > 1.d-8 .and. &
          OptionPrintToScreen(option)) then
        print *, 'Time (', temp_times(i), ') in ', trim(input%filename), &
                 ' does not match previous allocation time: ', &
                 flow_dataset%time_series%times(i), i
        stop
      endif
    enddo
  else
    allocate(flow_dataset%time_series%times(count))
  endif

  if (associated(flow_dataset%time_series%values)) &
    deallocate(flow_dataset%time_series%values)
  allocate(flow_dataset%time_series%values(flow_dataset%time_series%rank,count))

  flow_dataset%time_series%times(1:count) = temp_times(1:count)
  flow_dataset%time_series%values(:,1:count) = temp_array(:,1:count)
  
    
  deallocate(temp_times)
  deallocate(temp_array)

  if (len_trim(time_units) > 0) then
    ! Times
    conversion = UnitsConvertToInternal(time_units,option)
    flow_dataset%time_series%times(1:count) = conversion * &
      flow_dataset%time_series%times(1:count)
  endif
  if (len_trim(data_units) > 0) then
    ! Data
    conversion = UnitsConvertToInternal(data_units,option)
    flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1:count) = &
      conversion * &
      flow_dataset%time_series%values(1:flow_dataset%time_series%rank,1:count)
  endif
  
end subroutine FlowConditionReadValuesFromFile2

! ************************************************************************** !
!
! FlowConditionPrint: Prints flow condition info
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine FlowConditionPrint(condition,option)

  use Option_module

  implicit none
  
  type(flow_condition_type) :: condition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Flow Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') trim(condition%length_units)
  
  do i=1, condition%num_sub_conditions
    call FlowConditionPrintSubCondition(condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)
  
end subroutine FlowConditionPrint

! ************************************************************************** !
!
! FlowConditionPrintSubCondition: Prints flow subcondition info
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine FlowConditionPrintSubCondition(subcondition,option)

  use Option_module

  implicit none
  
  type(flow_sub_condition_type) :: subcondition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  select case(subcondition%itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(DIRICHLET_ZERO_GRADIENT_BC)
      string = 'dirichlet-zero gradient'
    case(MASS_RATE_SS)
      string = 'mass_rate'
    case(WELL_SS)
      string = 'well'
    case(HYDROSTATIC_BC)
      string = 'hydrostatic'
    case(CONDUCTANCE_BC)
      string = 'conductance'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
    case(SEEPAGE_BC)
      string = 'seepage'
    case(VOLUMETRIC_RATE_SS)
      string = 'volumetric rate'
    case(EQUILIBRIUM_SS)
      string = 'equilibrium'
    case(UNIT_GRADIENT_BC)
      string = 'unit gradient'
    case(SCALED_MASS_RATE_SS)
      string = 'scaled mass rate'
    case(SCALED_VOLUMETRIC_RATE_SS)
      string = 'scaled volumetric rate'
    case(DISTRIBUTED_VOLUMETRIC_RATE_SS)
      string = 'distributed volumetric rate'
    case(DISTRIBUTED_MASS_RATE_SS)
      string = 'distributed mass rate'
  end select
  100 format(6x,'Type: ',a)  
  write(option%fid_out,100) trim(string)
  
  110 format(6x,a)  

  write(option%fid_out,110) 'Datum:'
  if (associated(subcondition%datum%time_series)) then
    call TimeSeriesPrint(subcondition%datum%time_series,option)
  endif
  if (associated(subcondition%datum%dataset)) then
    call DatasetPrint(subcondition%datum%dataset,option)
  endif
  
  write(option%fid_out,110) 'Gradient:'
  if (associated(subcondition%gradient%time_series)) then
    call TimeSeriesPrint(subcondition%gradient%time_series,option)
  endif
  if (associated(subcondition%gradient%dataset)) then
    call DatasetPrint(subcondition%gradient%dataset,option)
  endif

  write(option%fid_out,110) 'Dataset:'
  if (associated(subcondition%flow_dataset%time_series)) then
    call TimeSeriesPrint(subcondition%flow_dataset%time_series,option)
  endif
  if (associated(subcondition%flow_dataset%dataset)) then
    call DatasetPrint(subcondition%flow_dataset%dataset,option)
  endif
            
end subroutine FlowConditionPrintSubCondition
 
! ************************************************************************** !
!
! FlowConditionDatasetPrint: Prints flow condition dataset info
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine FlowConditionDatasetPrint(flow_dataset,option)

  use Option_module

  implicit none
  
  type(flow_condition_dataset_type) :: flow_dataset
  type(option_type) :: option
  
  if(associated(flow_dataset%time_series)) then
    call TimeSeriesPrint(flow_dataset%time_series,option)
  endif
  if(associated(flow_dataset%dataset)) then
    !TODO(geh): setup
  endif
  
end subroutine FlowConditionDatasetPrint

! ************************************************************************** !
!
! FlowConditionUpdate: Updates a transient condition
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine FlowConditionUpdate(condition_list,option,time)

  use Option_module
  
  implicit none
  
  type(condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  
  type(flow_condition_type), pointer :: condition
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition  
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr
      
      if (associated(sub_condition)) then
        call FlowSubConditionUpdateDataset(option,time, &
                                           sub_condition%flow_dataset)
        call FlowSubConditionUpdateDataset(option,time, &
                                           sub_condition%datum)
        call FlowSubConditionUpdateDataset(option,time, &
                                           sub_condition%gradient)
      endif
      
    enddo
      
    condition => condition%next
    
  enddo
  
end subroutine FlowConditionUpdate

! ************************************************************************** !
!
! FlowSubConditionUpdateDataset: Updates a transient condition dataset
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine FlowSubConditionUpdateDataset(option,time,flow_condition_dataset)

  use Option_module
  use Dataset_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: time
  type(flow_condition_dataset_type) :: flow_condition_dataset

  if (associated(flow_condition_dataset%time_series)) then
    if (time < 1.d-40 .or. &
        flow_condition_dataset%time_series%is_transient) then
      call TimeSeriesUpdate(option,time,flow_condition_dataset%time_series)
    endif
  endif
  
  if (associated(flow_condition_dataset%dataset)) then
    if ((time < 1.d-40 .or. &
         flow_condition_dataset%dataset%is_transient) .and. &
        .not.flow_condition_dataset%dataset%is_cell_indexed) then
      call DatasetLoad(flow_condition_dataset%dataset,option)
    endif
  endif
  
end subroutine FlowSubConditionUpdateDataset

! ************************************************************************** !
!
! TranConditionUpdate: Updates a transient transport condition
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine TranConditionUpdate(condition_list,option,time)

  use Option_module
  
  implicit none
  
  type(tran_condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  
  type(tran_condition_type), pointer :: condition
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    do
      if (associated(condition%cur_constraint_coupler%next)) then
        if (time >= condition%cur_constraint_coupler%next%time) then
          condition%cur_constraint_coupler => &
            condition%cur_constraint_coupler%next
        else
          exit
        endif
      else 
        exit
      endif
    enddo
    condition => condition%next
    
  enddo
  
end subroutine TranConditionUpdate

! ************************************************************************** !
!
! FlowConditionInitList: Initializes a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine FlowConditionInitList(list)

  implicit none

  type(condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine FlowConditionInitList

! ************************************************************************** !
!
! FlowConditionAddToList: Adds a new condition to a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine FlowConditionAddToList(new_condition,list)

  implicit none
  
  type(flow_condition_type), pointer :: new_condition
  type(condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine FlowConditionAddToList

! ************************************************************************** !
!
! FlowConditionGetPtrFromList: Returns a pointer to the condition matching &
!                          condition_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function FlowConditionGetPtrFromList(condition_name,condition_list)

  use String_module
  
  implicit none
  
  type(flow_condition_type), pointer :: FlowConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list
 
  PetscInt :: length
  type(flow_condition_type), pointer :: condition
    
  nullify(FlowConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      FlowConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function FlowConditionGetPtrFromList

! ************************************************************************** !
!
! TranConditionInitList: Initializes a transport condition list
! author: Glenn Hammond
! date: 10/13/08
!
! ************************************************************************** !
subroutine TranConditionInitList(list)

  implicit none

  type(tran_condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine TranConditionInitList

! ************************************************************************** !
!
! TranConstraintInitList: Initializes a transport constraint list
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintInitList(list)

  implicit none

  type(tran_constraint_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_constraints = 0

end subroutine TranConstraintInitList

! ************************************************************************** !
!
! TranConditionAddToList: Adds a new condition to a transport condition list
! author: Glenn Hammond
! date: 10/13/08
!
! ************************************************************************** !
subroutine TranConditionAddToList(new_condition,list)

  implicit none
  
  type(tran_condition_type), pointer :: new_condition
  type(tran_condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine TranConditionAddToList

! ************************************************************************** !
!
! TranConstraintAddToList: Adds a new constraint to a transport constraint
!                          list
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintAddToList(new_constraint,list)

  implicit none
  
  type(tran_constraint_type), pointer :: new_constraint
  type(tran_constraint_list_type) :: list
  
  list%num_constraints = list%num_constraints + 1
  new_constraint%id = list%num_constraints
  if (.not.associated(list%first)) list%first => new_constraint
  if (associated(list%last)) list%last%next => new_constraint
  list%last => new_constraint
  
end subroutine TranConstraintAddToList

! ************************************************************************** !
!
! TranConditionGetPtrFromList: Returns a pointer to the condition matching
!                              condition_name
! author: Glenn Hammond
! date: 10/13/08
!
! ************************************************************************** !
function TranConditionGetPtrFromList(condition_name,condition_list)

  use String_module

  implicit none
  
  type(tran_condition_type), pointer :: TranConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(tran_condition_list_type) :: condition_list
 
  PetscInt :: length
  type(tran_condition_type), pointer :: condition
    
  nullify(TranConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                        length)) then
      TranConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function TranConditionGetPtrFromList

! ************************************************************************** !
!
! TranConstraintGetPtrFromList: Returns a pointer to the constraint matching
!                               constraint_name
! author: Glenn Hammond
! date: 10/13/08
!
! ************************************************************************** !
function TranConstraintGetPtrFromList(constraint_name,constraint_list)

  use String_module

  implicit none
  
  type(tran_constraint_type), pointer :: TranConstraintGetPtrFromList
  character(len=MAXWORDLENGTH) :: constraint_name
  type(tran_constraint_list_type) :: constraint_list
 
  PetscInt :: length
  type(tran_constraint_type), pointer :: constraint
    
  nullify(TranConstraintGetPtrFromList)
  constraint => constraint_list%first
  
  do 
    if (.not.associated(constraint)) exit
    length = len_trim(constraint_name)
    if (length == len_trim(constraint%name) .and. &
        StringCompare(constraint%name,constraint_name, &
                        length)) then
      TranConstraintGetPtrFromList => constraint
      return
    endif
    constraint => constraint%next
  enddo
  
end function TranConstraintGetPtrFromList

! ************************************************************************** !
!
! FlowConditionIsTransient: Returns PETSC_TRUE
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
function FlowConditionIsTransient(condition)

  implicit none
  
  type(flow_condition_type) :: condition
  
  PetscBool :: FlowConditionIsTransient
  
  FlowConditionIsTransient = PETSC_FALSE

  !TODO(geh): add check for general condition
  
  ! pressure
  if (FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%concentration) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%well) .or. &
      FlowSubConditionIsTransient(condition%enthalpy)) then
    FlowConditionIsTransient = PETSC_TRUE
  endif
  
end function FlowConditionIsTransient

! ************************************************************************** !
!
! FlowSubConditionIsTransient: Returns PETSC_TRUE
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
function FlowSubConditionIsTransient(sub_condition)

  implicit none
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  PetscBool :: FlowSubConditionIsTransient
  
  FlowSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (FlowDatasetIsTransient(sub_condition%flow_dataset) .or. &
        FlowDatasetIsTransient(sub_condition%gradient) .or. &
        FlowDatasetIsTransient(sub_condition%datum)) then
      FlowSubConditionIsTransient = PETSC_TRUE
    endif
  endif  
  
end function FlowSubConditionIsTransient

! ************************************************************************** !
!
! FlowDatasetIsTransient: Returns PETSC_TRUE
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
function FlowDatasetIsTransient(flow_dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: flow_dataset
  
  PetscBool :: FlowDatasetIsTransient
  
  FlowDatasetIsTransient = PETSC_FALSE

  if (associated(flow_dataset%time_series)) then
    if (flow_dataset%time_series%is_transient) then
      FlowDatasetIsTransient = PETSC_TRUE
    endif
  endif
  if (associated(flow_dataset%dataset)) then
    if (flow_dataset%dataset%is_transient) then
      FlowDatasetIsTransient = PETSC_TRUE
    endif
  endif
  
end function FlowDatasetIsTransient

! ************************************************************************** !
!
! FlowConditionDestroyList: Deallocates a list of conditions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine FlowConditionDestroyList(condition_list)

  implicit none
  
  type(condition_list_type), pointer :: condition_list
  
  type(flow_condition_type), pointer :: condition, prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call FlowConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine FlowConditionDestroyList

! ************************************************************************** !
!
! FlowConditionDestroy: Deallocates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine FlowConditionDestroy(condition)

  implicit none
  
  type(flow_condition_type), pointer :: condition
  
  PetscInt :: i
  
  if (.not.associated(condition)) return
  
  if (associated(condition%sub_condition_ptr)) then
    do i=1,condition%num_sub_conditions
      call FlowSubConditionDestroy(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  
  nullify(condition%pressure)
  nullify(condition%rate)
  nullify(condition%well)
  nullify(condition%saturation)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)

  call FlowGeneralConditionDestroy(condition%general)
  
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine FlowConditionDestroy

! ************************************************************************** !
!
! FlowGeneralConditionDestroy:Destroys a general mode condition
! author: Glenn Hammond
! date: 05/26/11
!
! ************************************************************************** !
subroutine FlowGeneralConditionDestroy(general_condition)

  use Option_module
  
  implicit none
  
  type(flow_general_condition_type), pointer :: general_condition
  
  if (.not.associated(general_condition)) return

  call FlowSubConditionDestroy(general_condition%liquid_pressure)
  call FlowSubConditionDestroy(general_condition%gas_pressure)
  call FlowSubConditionDestroy(general_condition%gas_saturation)
  call FlowSubConditionDestroy(general_condition%mole_fraction)
  call FlowSubConditionDestroy(general_condition%temperature)
  call FlowSubConditionDestroy(general_condition%flux)

  deallocate(general_condition)
  nullify(general_condition)

end subroutine FlowGeneralConditionDestroy

! ************************************************************************** !
!
! FlowSubConditionDestroy: Destroys a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine FlowSubConditionDestroy(sub_condition)

  implicit none
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  if (.not.associated(sub_condition)) return
  
  call FlowConditionDatasetDestroy(sub_condition%flow_dataset)
  call FlowConditionDatasetDestroy(sub_condition%datum)
  call FlowConditionDatasetDestroy(sub_condition%gradient)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine FlowSubConditionDestroy

! ************************************************************************** !
!
! FlowConditionDatasetDestroy: Destroys a dataset associated with a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine FlowConditionDatasetDestroy(flow_condition_dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: flow_condition_dataset
  
  if (associated(flow_condition_dataset%time_series)) then
    call TimeSeriesDestroy(flow_condition_dataset%time_series)
  endif
  nullify(flow_condition_dataset%time_series)
  ! since datasets reside in the dataset list, we simple unlink the pointer
  ! without destroying the object.
  if (associated(flow_condition_dataset%dataset)) &
    nullify(flow_condition_dataset%dataset)

end subroutine FlowConditionDatasetDestroy

! ************************************************************************** !
!
! TranConditionDestroyList: Deallocates a list of conditions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine TranConditionDestroyList(condition_list)

  implicit none
  
  type(tran_condition_list_type), pointer :: condition_list
  
  type(tran_condition_type), pointer :: condition, prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call TranConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine TranConditionDestroyList

! ************************************************************************** !
!
! TranConstraintDestroy: Deallocates a constraint
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintDestroy(constraint)

  implicit none
  
  type(tran_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return

  if (associated(constraint%aqueous_species)) &
    call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
  nullify(constraint%aqueous_species)
  if (associated(constraint%minerals)) &
    call MineralConstraintDestroy(constraint%minerals)
  nullify(constraint%minerals)
  if (associated(constraint%surface_complexes)) &
    call SurfaceComplexConstraintDestroy(constraint%surface_complexes)
  nullify(constraint%surface_complexes)
  if (associated(constraint%colloids)) &
    call ColloidConstraintDestroy(constraint%colloids)
  nullify(constraint%colloids)

  deallocate(constraint)
  nullify(constraint)

end subroutine TranConstraintDestroy

! ************************************************************************** !
!
! TranConstraintDestroyList: Deallocates a list of constraints
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintDestroyList(constraint_list)

  implicit none
  
  type(tran_constraint_list_type), pointer :: constraint_list
  
  type(tran_constraint_type), pointer :: constraint, prev_constraint
  
  if (.not.associated(constraint_list)) return
  
  constraint => constraint_list%first
  do 
    if (.not.associated(constraint)) exit
    prev_constraint => constraint
    constraint => constraint%next
    call TranConstraintDestroy(prev_constraint)
  enddo
  
  constraint_list%num_constraints = 0
  nullify(constraint_list%first)
  nullify(constraint_list%last)
  if (associated(constraint_list%array)) deallocate(constraint_list%array)
  nullify(constraint_list%array)
  
  deallocate(constraint_list)
  nullify(constraint_list)

end subroutine TranConstraintDestroyList

! ************************************************************************** !
!
! TranConditionDestroy: Deallocates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine TranConditionDestroy(condition)

  implicit none
  
  type(tran_condition_type), pointer :: condition
  
  if (.not.associated(condition)) return
  
  if (associated(condition%constraint_coupler_list)) &
    call TranConstraintCouplerDestroy(condition%constraint_coupler_list)

  deallocate(condition)
  nullify(condition)

end subroutine TranConditionDestroy

! ************************************************************************** !
!
! TranConstraintCouplerDestroy: Destroys a constraint coupler linked list
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintCouplerDestroy(coupler_list)

  use Option_module
  
  implicit none
  
  type(tran_constraint_coupler_type), pointer :: coupler_list
  
  type(tran_constraint_coupler_type), pointer :: cur_coupler, prev_coupler
  
  cur_coupler => coupler_list
  
  do
    if (.not.associated(cur_coupler)) exit
    prev_coupler => cur_coupler
    cur_coupler => cur_coupler%next
    if (associated(prev_coupler%rt_auxvar)) then
      call RTAuxVarDestroy(prev_coupler%rt_auxvar)
    endif
    nullify(prev_coupler%rt_auxvar)
    if (associated(prev_coupler%global_auxvar)) then
      call GlobalAuxVarDestroy(prev_coupler%global_auxvar)
    endif
    nullify(prev_coupler%global_auxvar)
    nullify(prev_coupler%aqueous_species)
    nullify(prev_coupler%minerals)
    nullify(prev_coupler%surface_complexes)
    nullify(prev_coupler%next)
    deallocate(prev_coupler)
    nullify(prev_coupler)
  enddo
  
  nullify(coupler_list)
  
end subroutine TranConstraintCouplerDestroy

end module Condition_module
