module Condition_module
 
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2

  type, public :: flow_condition_dataset_type
    PetscInt :: rank
    PetscBool :: is_transient
    PetscBool :: is_cyclic
    PetscInt :: interpolation_method
    PetscReal, pointer :: times(:)
    PetscReal, pointer :: values(:,:)
    PetscReal, pointer :: cur_value(:)
    PetscInt :: cur_time_index
    PetscInt :: max_time_index
    PetscReal :: time_shift
    PetscReal :: lame_aux_variable_remove_me
  end type flow_condition_dataset_type
  
  type, public :: flow_condition_type
    PetscInt :: id                                 ! id from which condition can be referenced
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_general_condition_type), pointer :: general
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type flow_condition_type

  ! data structure for general phase
  type, public :: flow_general_condition_type
    type(flow_sub_condition_type), pointer :: liquid_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: mole_fraction
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: flux
  end type flow_general_condition_type
    
  type, public :: flow_sub_condition_type
    PetscInt :: itype                  ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    type(flow_condition_dataset_type) :: datum
    type(flow_condition_dataset_type) :: gradient
    type(flow_condition_dataset_type) :: dataset
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
    PetscInt :: id                                ! id from which condition can be referenced
    PetscInt :: itype                  ! integer describing type of condition
    PetscBool :: is_transient
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
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
    PetscInt :: iflag
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
  nullify(condition%rate)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%sub_condition_ptr)
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
  coupler%iflag = 0
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
  sub_condition%ctype = ''
  sub_condition%name = ''

  call FlowConditionDatasetInit(sub_condition%dataset)
  sub_condition%dataset%rank = ndof
  call FlowConditionDatasetInit(sub_condition%gradient)
  sub_condition%gradient%rank = 3
  call FlowConditionDatasetInit(sub_condition%datum)
  sub_condition%datum%rank = 3

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
subroutine FlowConditionDatasetInit(dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: dataset

  nullify(dataset%times)
  nullify(dataset%values)
  nullify(dataset%cur_value)
  dataset%cur_time_index = 0
  dataset%max_time_index = 0
  dataset%rank = 0
  dataset%is_cyclic = PETSC_FALSE
  dataset%is_transient = PETSC_FALSE
  dataset%interpolation_method = NULL
  dataset%time_shift = 0.d0
  dataset%lame_aux_variable_remove_me = 0.d0
    
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
                                  default_dataset, &
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
  type(flow_condition_dataset_type) :: default_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  PetscBool :: destroy_if_null

  PetscInt :: array_size

  if (.not.associated(sub_condition)) return
  
  if (.not.associated(sub_condition%dataset%values)) then
    if (destroy_if_null) call FlowSubConditionDestroy(sub_condition)
    return
  endif
  
  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call printErrMsg(option)
  endif
  
  call FlowConditionDatasetVerify(option,condition%name,sub_condition_name, &
                                  default_time,sub_condition%dataset, &
                                  default_dataset)
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
subroutine FlowConditionDatasetVerify(option, condition_name, sub_condition_name, &
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
  
  PetscInt :: array_size
  
  if (default_dataset%is_cyclic) dataset%is_cyclic = PETSC_TRUE
  if (dataset%interpolation_method == NULL) &
    dataset%interpolation_method = default_dataset%interpolation_method
  
  dataset%max_time_index = 1
  if (.not.associated(dataset%times)) then
    if (.not.associated(default_dataset%times)) then
      array_size = 1
      allocate(dataset%times(array_size))
      dataset%times = default_time
    else
      array_size = size(default_dataset%times,1)
      allocate(dataset%times(array_size))
      dataset%times(1:array_size) = default_dataset%times(1:array_size)
    endif
  endif
  if (.not.associated(dataset%values)) then
    if (.not.associated(default_dataset%values)) then
      ! allocate to size 1
      array_size = 1
      allocate(dataset%values(1:dataset%rank,array_size))
      dataset%values(1:dataset%rank,1:array_size) = 0.d0
    else
      ! copy default
      array_size = size(default_dataset%values,2)
      allocate(dataset%values(1:dataset%rank,array_size))
      dataset%values(1:dataset%rank,1:array_size) = &
        default_dataset%values(1:dataset%rank,1:array_size)
    endif
  endif
  dataset%max_time_index = size(dataset%times,1) 
  if (dataset%max_time_index > 1) then
    if (size(dataset%values,2) /= & ! size of 2nd rank
        dataset%max_time_index) then
      write(size1,*) size(dataset%times,1)
      write(size2,*) size(dataset%values,2)
      option%io_buffer = 'times/values ('//trim(size1)//'/'//trim(size2) // &
                         ') array size mismatch in ' // &
                         'condition: ' // trim(condition_name) // &
                         'subcondition: ' // trim(sub_condition_name)
      call printErrMsg(option) 
    endif
    dataset%is_transient = PETSC_TRUE
  else
    dataset%is_transient = PETSC_FALSE
  endif  
  dataset%cur_time_index = 1
  
  allocate(dataset%cur_value(dataset%rank))
  dataset%cur_value(1:dataset%rank) = dataset%values(1:dataset%rank,1)

  dataset%time_shift = dataset%times(dataset%max_time_index)
! print *,'condition: ',dataset%max_time_index,dataset%time_shift

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
  
  type(flow_condition_dataset_type), pointer :: dataset
  PetscInt :: num_times
  PetscInt :: itime
  PetscReal :: time_shift
  PetscReal, allocatable :: temp_times(:)

  dataset => sub_condition%dataset

  if (.not.dataset%is_cyclic .or. dataset%max_time_index == 1) then
    allocate(times(dataset%max_time_index))
    times =  dataset%times
  else ! cyclic
    num_times = (int(max_sim_time/dataset%times(dataset%max_time_index))+1)* &
                dataset%max_time_index
    allocate(temp_times(num_times))
    temp_times = 0.d0

    num_times = 0
    itime = 0
    time_shift = 0.d0
    do
      num_times = num_times + 1
      itime = itime + 1
      ! exist for non-cyclic
      if (itime > dataset%max_time_index) exit
      temp_times(num_times) = dataset%times(itime) + time_shift
      if (mod(itime,dataset%max_time_index) == 0) then
        itime = 0
        time_shift = time_shift + dataset%times(dataset%max_time_index) 
      endif 
      ! exit for cyclic
      if (temp_times(num_times) >= max_sim_time) exit
    enddo

    allocate(times(num_times))
    times(:) = temp_times(1:num_times)
    deallocate(temp_times)
  endif
 
end subroutine FlowConditionDatasetGetTimes

#ifndef GLENN
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
                                       concentration, enthalpy, rate, &
                                       sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, idof
  PetscBool :: found
  PetscBool :: destroy_if_null
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read,ierr)

  default_time = 0.d0
  default_iphase = 0
  call FlowConditionDatasetInit(default_dataset)
  default_dataset%rank = 1
  default_dataset%interpolation_method = STEP
  default_dataset%is_cyclic = PETSC_FALSE
  call FlowConditionDatasetInit(default_datum)
  default_datum%rank = 3
  call FlowConditionDatasetInit(default_gradient)
  default_gradient%rank = 3
  
  pressure => FlowSubConditionCreate(option%nphase)
  pressure%name = 'pressure'
  flux => pressure
  rate => FlowSubConditionCreate(option%nflowspec)
  rate%name = 'rate'
  temperature => FlowSubConditionCreate(ONE_INTEGER)
  temperature%name = 'temperature'
  concentration => FlowSubConditionCreate(ONE_INTEGER)
  concentration%name = 'concentration'
  enthalpy => FlowSubConditionCreate(option%nphase)
  enthalpy%name = 'enthalpy'

  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  rate%units = 'kg/s'
  temperature%units = 'C'
  concentration%units = 'M'
  enthalpy%units = 'KJ/mol'
  
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
        default_dataset%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_dataset%interpolation_method = STEP
          case('linear') 
            default_dataset%interpolation_method = LINEAR
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
            case('FLUX')
              sub_condition_ptr => flux
            case('TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONCENTRATION','SATURATION')
              sub_condition_ptr => concentration
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
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
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('production_well')
              sub_condition_ptr%itype = PRODUCTION_WELL
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
            case('scaled_volumetric_rate')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('unit_gradient')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may ' // &
                  'only be associated with a PRESSURE flow condition.'
                call printErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC
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
            case('FLUX')
              sub_condition_ptr => flux
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONC','CONCENTRATION','SATURATION')
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
                                     temperature%dataset, &
                                     temperature%units)
      case('ENTHALPY','H')
        call FlowConditionReadValues(input,option,word,string, &
                                     enthalpy%dataset, &
                                     enthalpy%units)
      case('PRESSURE','PRES','PRESS')
        call FlowConditionReadValues(input,option,word,string, &
                                     pressure%dataset, &
                                     pressure%units)
      case('RATE')
        call FlowConditionReadValues(input,option,word,string, &
                                     rate%dataset, &
                                     rate%units)
      case('FLUX','VELOCITY','VEL')
        call FlowConditionReadValues(input,option,word,string, &
                                     pressure%dataset, &
                                     pressure%units)
      case('CONC','CONCENTRATION','SATURATION')
        call FlowConditionReadValues(input,option,word,string, &
                                     concentration%dataset, &
                                     concentration%units)
      case('CONDUCTANCE')
        call InputReadDouble(input,option,pressure%dataset%lame_aux_variable_remove_me)
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
  
  ! update datum and gradient defaults, if null, based on dataset default
  if (default_dataset%is_cyclic) then
    default_datum%is_cyclic = PETSC_TRUE
    default_gradient%is_cyclic = PETSC_TRUE
  endif
  if (default_datum%interpolation_method == NULL) &
    default_datum%interpolation_method = default_dataset%interpolation_method
  if (default_gradient%interpolation_method == NULL) &
    default_gradient%interpolation_method = default_dataset%interpolation_method

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
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS,PRODUCTION_WELL)
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
                              default_dataset, &
                              default_datum, default_gradient, PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'concentration'
  call FlowSubConditionVerify(option,condition,word,concentration,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)

  select case(option%iflowmode)
    case(THC_MODE, MPH_MODE, IMS_MODE, FLASH2_MODE, G_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)) then
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

      if (option%iflowmode == G_MODE) then
        if (associated(pressure)) condition%sub_condition_ptr(GENERAL_LIQUID_PRESSURE_DOF)%ptr => pressure
        if (associated(rate)) condition%sub_condition_ptr(GENERAL_LIQUID_FLUX_DOF)%ptr => rate
        condition%sub_condition_ptr(GENERAL_TEMPERATURE_DOF)%ptr => temperature
        condition%sub_condition_ptr(GENERAL_CONCENTRATION_DOF)%ptr => concentration
        if (associated(enthalpy)) condition%sub_condition_ptr(GENERAL_ENTHALPY_DOF)%ptr => enthalpy
        
        allocate(condition%itype(FOUR_INTEGER))
        condition%itype = 0
        if (associated(pressure)) condition%itype(GENERAL_LIQUID_PRESSURE_DOF) = pressure%itype
        if (associated(rate)) condition%itype(GENERAL_LIQUID_FLUX_DOF) = rate%itype
        condition%itype(GENERAL_TEMPERATURE_DOF) = temperature%itype
        condition%itype(GENERAL_CONCENTRATION_DOF) = concentration%itype
        if (associated(enthalpy)) condition%itype(GENERAL_ENTHALPY_DOF) = concentration%itype
      else
        ! must be in this order, which matches the dofs i problem
        if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
        if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
        condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
        if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
        
        allocate(condition%itype(FIVE_INTEGER))
        condition%itype = 0
        if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
        if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
        condition%itype(TWO_INTEGER) = temperature%itype
        condition%itype(THREE_INTEGER) = concentration%itype
        if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype
      endif
    
    case(RICHARDS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate) .and. &
          .not.associated(concentration)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)      
      endif                         
      if (associated(concentration)) then
        condition%concentration => concentration
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
      elseif (associated(concentration)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => concentration
      elseif (associated(rate)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      endif                         

      allocate(condition%itype(ONE_INTEGER))
      if (associated(pressure)) then 
        condition%itype(ONE_INTEGER) = pressure%itype
      else if (associated(concentration)) then
        condition%itype(ONE_INTEGER) = concentration%itype
      else if (associated(rate)) then
        condition%itype(ONE_INTEGER) = rate%itype
      endif
      
      ! these are not used with richards
      if (associated(temperature)) call FlowSubConditionDestroy(temperature)
      if (associated(enthalpy)) call FlowSubConditionDestroy(enthalpy)
      
  end select
  
  call FlowConditionDatasetDestroy(default_dataset)
  call FlowConditionDatasetDestroy(default_datum)
  call FlowConditionDatasetDestroy(default_gradient)
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr)

end subroutine FlowConditionRead

#else
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
  type(flow_general_condition_type), pointer :: general
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_dataset
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
  call FlowConditionDatasetInit(default_dataset)
  default_dataset%rank = 1
  default_dataset%interpolation_method = STEP
  default_dataset%is_cyclic = PETSC_FALSE
  call FlowConditionDatasetInit(default_datum)
  default_datum%rank = 3
  call FlowConditionDatasetInit(default_gradient)
  default_gradient%rank = 3
  
  select case(option%iflowmode)
    case(G_MODE)
      general => FlowGeneralConditionCreate(option)
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
        default_dataset%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_dataset%interpolation_method = STEP
          case('linear') 
            default_dataset%interpolation_method = LINEAR
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
                                       sub_condition_ptr%dataset, &
                                       sub_condition_ptr%units)
        select case(word)
          case('LIQUID_SATURATION') ! convert to gas saturation
            do i = 1, size(sub_condition_ptr%dataset%values)
              sub_condition_ptr%dataset%values(i,:) = &
                1.d0 - sub_condition_ptr%dataset%values(i,:)
            enddo
        end select
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized in flow condition'
        call printErrMsg(option)                                 
    end select 
  
  enddo  
  
  ! update datum and gradient defaults, if null, based on dataset default
  if (default_dataset%is_cyclic) then
    default_datum%is_cyclic = PETSC_TRUE
    default_gradient%is_cyclic = PETSC_TRUE
  endif
  if (default_datum%interpolation_method == NULL) &
    default_datum%interpolation_method = default_dataset%interpolation_method
  if (default_gradient%interpolation_method == NULL) &
    default_gradient%interpolation_method = default_dataset%interpolation_method

  select case(option%iflowmode)
    case(G_MODE)
      ! need mole fraction and some sort of saturation
      if (associated(general%flux)) then
        ! neumann or mass/volumetric flux
        ! need temperature
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
        if (associated(general%gas_saturation)) then
          ! two phase condition

        else if (associated(general%mole_fraction)) then

        else

        endif
      else
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
        if (associated(general%gas_saturation)) then
          ! two phase condition
        else if (associated(general%mole_fraction)) then

        else

        endif
      endif
  end select 
    


#if 0
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
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS,PRODUCTION_WELL)
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
                              default_dataset, &
                              default_datum, default_gradient, PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'concentration'
  call FlowSubConditionVerify(option,condition,word,concentration,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient,PETSC_TRUE)

  select case(option%iflowmode)
    case(THC_MODE, MPH_MODE, IMS_MODE, FLASH2_MODE, G_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)) then
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

      if (option%iflowmode == G_MODE) then
        if (associated(pressure)) condition%sub_condition_ptr(GENERAL_LIQUID_PRESSURE_DOF)%ptr => pressure
        if (associated(rate)) condition%sub_condition_ptr(GENERAL_LIQUID_FLUX_DOF)%ptr => rate
        condition%sub_condition_ptr(GENERAL_TEMPERATURE_DOF)%ptr => temperature
        condition%sub_condition_ptr(GENERAL_CONCENTRATION_DOF)%ptr => concentration
        if (associated(enthalpy)) condition%sub_condition_ptr(GENERAL_ENTHALPY_DOF)%ptr => enthalpy
        
        allocate(condition%itype(FOUR_INTEGER))
        condition%itype = 0
        if (associated(pressure)) condition%itype(GENERAL_LIQUID_PRESSURE_DOF) = pressure%itype
        if (associated(rate)) condition%itype(GENERAL_LIQUID_FLUX_DOF) = rate%itype
        condition%itype(GENERAL_TEMPERATURE_DOF) = temperature%itype
        condition%itype(GENERAL_CONCENTRATION_DOF) = concentration%itype
        if (associated(enthalpy)) condition%itype(GENERAL_ENTHALPY_DOF) = concentration%itype
      else
        ! must be in this order, which matches the dofs i problem
        if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
        if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
        condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
        if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
        
        allocate(condition%itype(FIVE_INTEGER))
        condition%itype = 0
        if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
        if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
        condition%itype(TWO_INTEGER) = temperature%itype
        condition%itype(THREE_INTEGER) = concentration%itype
        if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype
      endif
    
    case(RICHARDS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate) .and. &
          .not.associated(concentration)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)      
      endif                         
      if (associated(concentration)) then
        condition%concentration => concentration
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
      elseif (associated(concentration)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => concentration
      elseif (associated(rate)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      endif                         

      allocate(condition%itype(ONE_INTEGER))
      if (associated(pressure)) then 
        condition%itype(ONE_INTEGER) = pressure%itype
      else if (associated(concentration)) then
        condition%itype(ONE_INTEGER) = concentration%itype
      else if (associated(rate)) then
        condition%itype(ONE_INTEGER) = rate%itype
      endif
      
      ! these are not used with richards
      if (associated(temperature)) call FlowSubConditionDestroy(temperature)
      if (associated(enthalpy)) call FlowSubConditionDestroy(enthalpy)
      
  end select
  
  call FlowConditionDatasetDestroy(default_dataset)
  call FlowConditionDatasetDestroy(default_datum)
  call FlowConditionDatasetDestroy(default_gradient)
#endif
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr)

end subroutine FlowConditionRead
#endif

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
  PetscInt :: icomp
  PetscInt :: isrfcplx
  PetscInt :: length
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  type(srfcplx_constraint_type), pointer :: srfcplx_constraint
  type(colloid_constraint_type), pointer :: colloid_constraint
  PetscErrorCode :: ierr

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

        aq_species_constraint => AqueousSpeciesConstraintCreate(reaction,option)

        icomp = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT, CONCENTRATIONS')
          
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
          call InputReadDouble(input,option,aq_species_constraint%constraint_conc(icomp))
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
              case('S','TOTAL_SORB')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL_SORB
              case('P','PH')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_PH
              case('L','LOG')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_LOG
              case('M','MINERAL','MNRL') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_MINERAL
              case('G','GAS') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_GAS
              case('SC','CONSTRAINT_SUPERCRIT_CO2') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_SUPERCRIT_CO2
              case('Z','CHG') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_CHARGE_BAL
              case default
                option%io_buffer = 'Keyword: ' // trim(word) // &
                         ' not recognized in constraint,concentration'
                call printErrMsg(option)
            end select 
            if (aq_species_constraint%constraint_type(icomp) == CONSTRAINT_MINERAL .or. &
                aq_species_constraint%constraint_type(icomp) == CONSTRAINT_GAS .or.&
                aq_species_constraint%constraint_type(icomp) == CONSTRAINT_SUPERCRIT_CO2) then
              call InputReadWord(input,option,aq_species_constraint%constraint_spec_name(icomp), &
                             PETSC_FALSE)
              call InputErrorMsg(input,option,'constraint name', &
                              'CONSTRAINT, CONCENTRATIONS') 
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
        
        if (associated(constraint%aqueous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint 
        
      case('MNRL','MINERALS')

        mineral_constraint => MineralConstraintCreate(reaction,option)

        icomp = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT, MINERALS')
          
          if (InputCheckExit(input,option)) exit          
          
          icomp = icomp + 1

          if (icomp > reaction%nkinmnrl) then
            option%io_buffer = &
                     'Number of mineral constraints exceeds number of ' // &
                     'kinetic minerals in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,mineral_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral name', &
                          'CONSTRAINT, MINERALS')  
          option%io_buffer = 'Constraint Minerals: ' // &
                             trim(mineral_constraint%names(icomp))
          call printMsg(option)
          call InputReadDouble(input,option,mineral_constraint%constraint_vol_frac(icomp))
          call InputErrorMsg(input,option,'volume fraction', &
                          'CONSTRAINT, MINERALS')          
          call InputReadDouble(input,option,mineral_constraint%constraint_area(icomp))
          call InputErrorMsg(input,option,'area', &
                          'CONSTRAINT, MINERALS')          
        
        enddo  
        
        if (icomp < reaction%nkinmnrl) then
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

        srfcplx_constraint => SurfaceComplexConstraintCreate(reaction,option)

        isrfcplx = 0
        do
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT, SURFACE_COMPLEXES')
          
          if (InputCheckExit(input,option)) exit          
          
          isrfcplx = isrfcplx + 1

          if (isrfcplx > reaction%nkinsrfcplx) then
            option%io_buffer = &
                     'Number of surface complex constraints exceeds number of ' // &
                     'kinetic surface complexes in constraint: ' // &
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
          call InputReadDouble(input,option,srfcplx_constraint%constraint_conc(isrfcplx))
          call InputErrorMsg(input,option,'concentration', &
                          'CONSTRAINT, SURFACE COMPLEX')          
        enddo  
        
        if (isrfcplx < reaction%nkinsrfcplx) then
          option%io_buffer = &
                   'Number of surface complex constraints is less than ' // &
                   'number of kinetic surface complexes in surface complex ' // &
                   'constraint.'
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
          call InputReadDouble(input,option,colloid_constraint%constraint_conc_mob(icomp))
          call InputErrorMsg(input,option,'mobile concentration', &
                          'CONSTRAINT, COLLOIDS')          
          call InputReadDouble(input,option,colloid_constraint%constraint_conc_imb(icomp))
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
subroutine FlowConditionReadValues(input,option,keyword,string,dataset,units)

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
  type(flow_condition_dataset_type) :: dataset
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
  if (length==FOUR_INTEGER .and. StringCompare(word,'file',length)) then  !sp 
    input%err_buf2 = trim(keyword) // ', FILE'
    input%err_buf = 'keyword'
    call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
    if (input%ierr == 0) then
      filename = string2
    else
      do
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,'FLOW_CONDITION')
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword','OUTPUT') 
        call StringToUpper(word)
        select case(trim(word))
          case('FILENAME')
            call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_TRUE)
            call InputErrorMsg(input,option,'filename','CONDITION')
          case('REALIZATION_DEPENDENT')
            ! we only want realization dependent if a realization id exists
            if (option%id > 0) then
              write(word,*) option%id
              realization_word = adjustl(word)
            else
              realization_word = ''
            endif
          case('HDF5_PATH')
            ! we assume that the remainder of the string is the path
            hdf5_path = adjustl(input%buf)
          case('UNITS')
        end select          
      enddo      
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
      if (dims(1)-1 == dataset%rank) then
        ! alright, the 2d data is layed out in C-style.  now place it in
        ! the appropriate arrays
        allocate(dataset%times(dims(2)))
        dataset%times = -999.d0
        allocate(dataset%values(dataset%rank,dims(2))) 
        dataset%values = -999.d0
        icount = 1
        do i = 1, dims(2)
          dataset%times(i) = real_buffer(icount)
          icount = icount + 1
          do irank = 1, dataset%rank
            dataset%values(irank,i) = real_buffer(icount)
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
      input2 => InputCreate(IUNIT_TEMP,filename)
      call FlowConditionReadValuesFromFile(input2,dataset,option)
      call InputDestroy(input2)
    endif
  else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then  !sp 
    call FlowConditionReadValuesFromFile(input,dataset,option)
  else
    input%buf = trim(string2)
    allocate(dataset%values(dataset%rank,1))
    do irank=1,dataset%rank
      call InputReadDouble(input,option,dataset%values(irank,1))
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
      dataset%values(1:dataset%rank,1) = &
        UnitsConvertToInternal(units,option) * &
        dataset%values(1:dataset%rank,1)
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
subroutine FlowConditionReadValuesFromFile(input,dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module
  use Units_module

  implicit none
  
  type(input_type) :: input
  type(flow_condition_dataset_type) :: dataset
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

  if (dataset%rank > 1) then
    allocate(temp_array2(max_size))
    temp_array2 = 0.d0
  endif
  
  if (dataset%rank > 2) then
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
    if (dataset%rank > 1) then
      call InputReadDouble(input,option,temp_array2(count))
      call InputErrorMsg(input,option,'array2','CONDITION FILE') 
    endif
    if (dataset%rank > 2) then
      call InputReadDouble(input,option,temp_array3(count))
      call InputErrorMsg(input,option,'array3','CONDITION FILE') 
    endif
    if (count+1 > max_size) then
      temp_max_size = max_size
      call reallocateRealArray(temp_times,max_size) 
      ! careful.  reallocateRealArray double max_size every time.
      i = temp_max_size
      call reallocateRealArray(temp_array1,i) 
      if (dataset%rank > 1) then
        i = temp_max_size
        call reallocateRealArray(temp_array2,i)
      endif
      if (dataset%rank > 2) then
        i = temp_max_size
        call reallocateRealArray(temp_array3,i)
      endif
    endif  
  enddo
  
  if (associated(dataset%times)) then
    if (count /= size(dataset%times,1) .and. &
        OptionPrintToScreen(option)) then
      print *, 'Number of times (', count, ') in ', trim(input%filename), &
               ' does not match previous allocation: ', size(dataset%times,1)
      stop
    endif
    do i=1,count
      if (dabs(dataset%times(i)-temp_times(i)) > 1.d-8 .and. &
          OptionPrintToScreen(option)) then
        print *, 'Time (', temp_times(i), ') in ', trim(input%filename), &
                 ' does not match previous allocation time: ', &
                 dataset%times(i), i
        stop
      endif
    enddo
  else
    allocate(dataset%times(count))
  endif

  if (associated(dataset%values)) deallocate(dataset%values)
  allocate(dataset%values(dataset%rank,count))

  dataset%times(1:count) = temp_times(1:count)
  dataset%values(1,1:count) = temp_array1(1:count)
  if (dataset%rank > 1) dataset%values(2,1:count) = temp_array2(1:count)
  if (dataset%rank > 2) dataset%values(3,1:count) = temp_array3(1:count)
  
  deallocate(temp_times)
  deallocate(temp_array1)
  if (dataset%rank > 1) deallocate(temp_array2)
  if (dataset%rank > 2) deallocate(temp_array3)

  if (len_trim(time_units) > 0) then
    ! Times
    conversion = UnitsConvertToInternal(time_units,option)
    dataset%times(1:count) = conversion * dataset%times(1:count)
  endif
  if (len_trim(data_units) > 0) then
    ! Data
    conversion = UnitsConvertToInternal(data_units,option)
    dataset%values(1:dataset%rank,1:count) = conversion * dataset%values(1:dataset%rank,1:count)
  endif
  
end subroutine FlowConditionReadValuesFromFile

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
    case(MASS_RATE_SS)
      string = 'mass_rate'
    case(HYDROSTATIC_BC)
      string = 'hydrostatic'
    case(CONDUCTANCE_BC)
      string = 'conductance'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
    case(PRODUCTION_WELL)
      string = 'production well'
    case(SEEPAGE_BC)
      string = 'seepage'
    case(VOLUMETRIC_RATE_SS)
      string = 'volumetric rate'
    case(EQUILIBRIUM_SS)
      string = 'equilibrium'
    case(UNIT_GRADIENT_BC)
      string = 'unit gradient'
  end select
  100 format(6x,'Type: ',a)  
  write(option%fid_out,100) trim(string)
  
  110 format(6x,a)  

  write(option%fid_out,110) 'Datum:'
  call FlowConditionPrintDataset(subcondition%datum,option)
  write(option%fid_out,110) 'Gradient:'
  call FlowConditionPrintDataset(subcondition%gradient,option)
  write(option%fid_out,110) 'Dataset:'
  call FlowConditionPrintDataset(subcondition%dataset,option)
            
end subroutine FlowConditionPrintSubCondition
 
! ************************************************************************** !
!
! FlowConditionPrintDataset: Prints flow condition dataset info
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine FlowConditionPrintDataset(dataset,option)

  use Option_module

  implicit none
  
  type(flow_condition_dataset_type) :: dataset
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  if (dataset%is_cyclic) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(8x,''Is cyclic: '',a)') trim(string)

  100 format(8x,'Is transient: ', a)
  if (dataset%is_transient) then
    write(option%fid_out,100) 'yes'
    write(option%fid_out,'(8x,''  Number of values: '', i7)') &
      dataset%max_time_index
    select case(dataset%interpolation_method)
      case(STEP)
        string = 'step'
      case(LINEAR)
        string = 'linear'
    end select
    write(option%fid_out,'(8x,''  Interpolation method: '', a)') trim(string)
    write(option%fid_out,'(8x,''Start value:'',es16.8)') dataset%cur_value(1)
  else
    write(option%fid_out,100) 'no'
    write(option%fid_out,'(8x,''Value:'',es16.8)') dataset%cur_value(1)
  endif

            
end subroutine FlowConditionPrintDataset

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
        call FlowSubConditionUpdateDataset(option,time,sub_condition%dataset)
        call FlowSubConditionUpdateDataset(option,time,sub_condition%datum)
        call FlowSubConditionUpdateDataset(option,time,sub_condition%gradient)
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
subroutine FlowSubConditionUpdateDataset(option,time,dataset)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: time
  PetscBool :: is_cyclic
  PetscInt :: interpolation_method
  type(flow_condition_dataset_type) :: dataset
  
  PetscInt :: irank
  PetscInt :: cur_time_index
  PetscInt :: next_time_index
  PetscReal :: time_fraction

 ! potentially for initial condition
  if (time < 1.d-40 .or. .not.dataset%is_transient) return  

  ! cycle times if at max_time_index and cyclic
  if (dataset%cur_time_index == dataset%max_time_index .and. &
      dataset%is_cyclic .and. dataset%max_time_index > 1) then
      
    do cur_time_index = 1, dataset%max_time_index
      dataset%times(cur_time_index) = dataset%times(cur_time_index) + &     
                               dataset%time_shift
    enddo
    dataset%cur_time_index = 1
  endif
 
  cur_time_index = dataset%cur_time_index
  next_time_index = min(dataset%cur_time_index+1, &
                        dataset%max_time_index)

  ! ensure that condition has started
  if (time >= dataset%times(cur_time_index) .or. &
      dabs(time-dataset%times(cur_time_index)) < 1.d-40) then

    ! find appropriate time interval
    do
      if (time < dataset%times(next_time_index) .or. &
          cur_time_index == next_time_index) &
        exit
      cur_time_index = next_time_index
      ! ensure that time index does not go beyond end of array
      if (next_time_index < dataset%max_time_index) &
        next_time_index = next_time_index + 1
    enddo
    
    dataset%cur_time_index = cur_time_index
    
    select case(dataset%interpolation_method)
      case(STEP) ! just use the current value
        do irank=1,dataset%rank
          dataset%cur_value(irank) =  &
                             dataset%values(irank,dataset%cur_time_index) 
        enddo 
      case(LINEAR) ! interpolate the value between times
        do irank=1,dataset%rank
          ! interpolate value based on time
          dataset%cur_time_index = cur_time_index
          if (cur_time_index < dataset%max_time_index) then
            time_fraction = (time-dataset%times(cur_time_index)) / &
                              (dataset%times(next_time_index) - &
                               dataset%times(cur_time_index))
            dataset%cur_value(irank) = dataset%values(irank,cur_time_index) + &
                                  time_fraction * &
                                  (dataset%values(irank,next_time_index) - &
                                   dataset%values(irank,cur_time_index))
          else
            dataset%cur_value(irank) =  &
                               dataset%values(irank,dataset%max_time_index) 
          endif
        enddo
    end select 
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

  ! pressure
  if (associated(condition%pressure)) then
    if (condition%pressure%dataset%is_transient .or. &
        condition%pressure%gradient%is_transient .or. &
        condition%pressure%datum%is_transient) &
      FlowConditionIsTransient = PETSC_TRUE
  endif
  
  ! temperature
  if (associated(condition%temperature)) then
    if (condition%temperature%dataset%is_transient .or. &
        condition%temperature%gradient%is_transient .or. &
        condition%temperature%datum%is_transient) &
      FlowConditionIsTransient = PETSC_TRUE
  endif

  ! concentration
  if (associated(condition%concentration)) then
    if (condition%concentration%dataset%is_transient .or. &
        condition%concentration%gradient%is_transient .or. &
        condition%concentration%datum%is_transient) &
      FlowConditionIsTransient = PETSC_TRUE
  endif

  ! rate
  if (associated(condition%rate)) then
    if (condition%rate%dataset%is_transient .or. &
        condition%rate%gradient%is_transient .or. &
        condition%rate%datum%is_transient) &
      FlowConditionIsTransient = PETSC_TRUE
  endif
  
  ! enthalpy
  if (associated(condition%enthalpy)) then
    if (condition%enthalpy%dataset%is_transient .or. &
        condition%enthalpy%gradient%is_transient .or. &
        condition%enthalpy%datum%is_transient) &
      FlowConditionIsTransient = PETSC_TRUE
  endif
   
end function FlowConditionIsTransient

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
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)

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
  
  call FlowConditionDatasetDestroy(sub_condition%dataset)
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
subroutine FlowConditionDatasetDestroy(dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: dataset
  
  if (associated(dataset%times)) deallocate(dataset%times)
  nullify(dataset%times)
  if (associated(dataset%values)) deallocate(dataset%values)
  nullify(dataset%values)  
  if (associated(dataset%cur_value)) deallocate(dataset%cur_value)
  nullify(dataset%cur_value)  

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
