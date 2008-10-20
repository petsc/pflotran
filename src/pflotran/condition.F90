module Condition_module
 
  use Reaction_Aux_module
  
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2

  type, public :: flow_condition_dataset_type
    PetscInt :: rank
    PetscTruth :: is_transient
    PetscTruth :: is_cyclic
    PetscInt :: interpolation_method
    PetscReal, pointer :: times(:)
    PetscReal, pointer :: values(:,:)
    PetscReal, pointer :: cur_value(:)
    PetscInt :: cur_time_index
    PetscInt :: max_time_index
  end type flow_condition_dataset_type
  
  type, public :: flow_condition_type
    PetscInt :: id                                 ! id from which condition can be referenced
    PetscTruth :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: mass_rate
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type flow_condition_type
  
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
    type(condition_ptr_type), pointer :: array(:)    
  end type condition_list_type
  
  type, public :: tran_condition_type
    PetscInt :: id                                ! id from which condition can be referenced
    PetscInt :: itype                  ! integer describing type of condition
    PetscTruth :: sync_time_with_update
    PetscTruth :: is_transient
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
    character(len=MAXWORDLENGTH) :: time_units
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(mineral_constraint_type), pointer :: minerals
    type(tran_constraint_coupler_type), pointer :: next   
  end type tran_constraint_coupler_type
      
  public :: ConditionCreate, ConditionDestroy, ConditionRead, &
            ConditionAddToList, ConditionInitList, ConditionDestroyList, &
            ConditionGetPtrFromList, ConditionUpdate, &
            TranConditionCreate, TranConstraintCreate, &
            TranConditionAddToList, TranConditionInitList, &
            TranConditionDestroyList, TranConditionGetPtrFromList, &
            TranConstraintAddToList, TranConstraintInitList, &
            TranConstraintDestroyList, TranConstraintGetPtrFromList, &
            TranConditionRead, TranConstraintRead, &
            TranConditionUpdate 
    
contains

! ************************************************************************** !
!
! ConditionCreate: Creates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function ConditionCreate(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_condition_type), pointer :: ConditionCreate
  
  type(flow_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%mass_rate)
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
  
  ConditionCreate => condition

end function ConditionCreate

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
  condition%sync_time_with_update = PETSC_FALSE
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
  nullify(coupler%next)
  coupler%constraint_name = ''
  coupler%time = 0.d0
  coupler%time_units = ''
  
  TranConstraintCouplerCreate => coupler

end function TranConstraintCouplerCreate

! ************************************************************************** !
!
! SubConditionCreate: Creates a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
function SubConditionCreate(ndof)

  use Option_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: SubConditionCreate
  
  PetscInt :: ndof
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''

  call ConditionDatasetInit(sub_condition%dataset)
  sub_condition%dataset%rank = ndof
  call ConditionDatasetInit(sub_condition%gradient)
  sub_condition%gradient%rank = 3
  call ConditionDatasetInit(sub_condition%datum)
  sub_condition%datum%rank = 3

  SubConditionCreate => sub_condition

end function SubConditionCreate

! ************************************************************************** !
!
! GetSubConditionFromArrayByName: returns a pointer to a subcondition with
!                                 matching name
! author: Glenn Hammond
! date: 06/02/08
!
! ************************************************************************** !
function GetSubConditionFromArrayByName(sub_condition_ptr_list,name)

  use Fileio_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: GetSubConditionFromArrayByName
  type(sub_condition_ptr_type), pointer :: sub_condition_ptr_list(:)
  character(len=MAXWORDLENGTH) :: name
  
  PetscInt :: idof
  
  nullify(GetSubConditionFromArrayByName)
  do idof = 1, size(sub_condition_ptr_list)
    if (len_trim(name) == len_trim(sub_condition_ptr_list(idof)%ptr%name) .and. &
        fiStringCompare(name,sub_condition_ptr_list(idof)%ptr%name,len_trim(name))) then
      GetSubConditionFromArrayByName => sub_condition_ptr_list(idof)%ptr
      return
    endif
  enddo
  
end function GetSubConditionFromArrayByName
          
! ************************************************************************** !
!
! ConditionDatasetInit: Initializes a dataset
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine ConditionDatasetInit(dataset)

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
    
end subroutine ConditionDatasetInit

! ************************************************************************** !
!
! SubConditionVerify: Verifies the data in a subcondition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine SubConditionVerify(option, condition, sub_condition_name, &
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
  PetscTruth :: default_cyclic
  PetscInt :: default_interpolation
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(flow_condition_dataset_type) :: default_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  PetscTruth :: destroy_if_null

  PetscInt :: array_size

  if (.not.associated(sub_condition)) return
  
  if (.not.associated(sub_condition%dataset%values)) then
    if (destroy_if_null) call SubConditionDestroy(sub_condition)
    return
  endif
  
  if (len_trim(sub_condition%ctype) < 1) then
    call printWrnMsg(option,'type of ' // trim(condition%name) // ':' // &
                     trim(sub_condition_name) // ' set to default')
    sub_condition%ctype = default_ctype
    sub_condition%itype = default_itype
  endif
  
  call ConditionDatasetVerify(option,condition%name,sub_condition_name, &
                              default_time,sub_condition%dataset,default_dataset)
  call ConditionDatasetVerify(option,condition%name,sub_condition_name, &
                              default_time,sub_condition%datum,default_datum)
  call ConditionDatasetVerify(option,condition%name,sub_condition_name, &
                              default_time,sub_condition%gradient,default_gradient)

end subroutine SubConditionVerify

! ************************************************************************** !
!
! ConditionDatasetVerify: Verifies the data in a dataset
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine ConditionDatasetVerify(option, condition_name, sub_condition_name, &
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
      dataset%values(1:dataset%rank,1:array_size) = default_dataset%values(1:dataset%rank,1:array_size)
    endif
  endif
  dataset%max_time_index = size(dataset%times,1) 
  if (dataset%max_time_index > 1) then
    if (size(dataset%values,2) /= & ! size of 2nd rank
        dataset%max_time_index) then
      write(size1,*) size(dataset%times,1)
      write(size2,*) size(dataset%values,2)
      call printErrMsg(option,'times/values ('//trim(size1)//'/'//trim(size2) // &
                       ') array size mismatch in ' // &
                       'condition: ' // trim(condition_name) // &
                       'subcondition: ' // trim(sub_condition_name)) 
    endif
    dataset%is_transient = PETSC_TRUE
  else
    dataset%is_transient = PETSC_FALSE
  endif  
  dataset%cur_time_index = 1
  
  allocate(dataset%cur_value(dataset%rank))
  dataset%cur_value(1:dataset%rank) = dataset%values(1:dataset%rank,1)


end subroutine ConditionDatasetVerify

! ************************************************************************** !
!
! ConditionRead: Reads a condition from the input file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine ConditionRead(condition,option,fid)

  use Option_module
  use Fileio_module
  use Logging_module  
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(flow_sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, mass_rate, &
                                       sub_condition_ptr
  PetscReal :: default_time = 0.d0
  PetscInt :: default_iphase = 0
  type(flow_condition_dataset_type) :: default_dataset
  type(flow_condition_dataset_type) :: default_datum
  type(flow_condition_dataset_type) :: default_gradient
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, length, idof
  logical :: found
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  call ConditionDatasetInit(default_dataset)
  default_dataset%rank = 1
  default_dataset%interpolation_method = STEP
  default_dataset%is_cyclic = PETSC_FALSE
  call ConditionDatasetInit(default_datum)
  default_datum%rank = 3
  call ConditionDatasetInit(default_gradient)
  default_gradient%rank = 3
  
  pressure => SubConditionCreate(option%nphase)
  flux => pressure
  mass_rate => SubConditionCreate(option%nflowspec)
  temperature => SubConditionCreate(ONE_INTEGER)
  concentration => SubConditionCreate(ONE_INTEGER)
  enthalpy => SubConditionCreate(option%nphase)

  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  mass_rate%units = 'kg/s'
  temperature%units = 'C'
  concentration%units = 'M'
  enthalpy%units = 'KJ/mol'
  
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
    if (fiCheckExit(string)) exit  

    call fiReadWord(string,word,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CONDITION', ierr)   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          if (ierr /= 0) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
            case('Pa','KPa')
              pressure%units = trim(word)
            case('kg/s','kg/yr')
              mass_rate%units = trim(word)
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
      case('INTERPOLATION')
        call fiReadWord(string,word,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'INTERPOLATION','CONDITION', ierr)   
        length = len_trim(word)
        call fiCharsToLower(word,length)
        select case(word)
          case('step')
            default_dataset%interpolation_method = STEP
          case('linear') 
            default_dataset%interpolation_method = LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call fiReadFlotranString(fid,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
          if (fiCheckExit(string)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'keyword','CONDITION,TYPE', ierr)   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
            case('MASS','MASS_RATE')
              sub_condition_ptr => mass_rate
            case('FLUX')
              sub_condition_ptr => flux
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONC','CONCENTRATION')
              sub_condition_ptr => concentration
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
            case default
              call printErrMsg(option,'keyword not recognized in condition,type')
          end select
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'TYPE','CONDITION', ierr)   
          length = len_trim(word)
          call fiCharsToLower(word,length)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('mass','mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
            case('hydrostatic','hydro','hydrostat','static')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('prod','production_well')
              sub_condition_ptr%itype = PRODUCTION_WELL
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('volume','volumetric','volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case default
              string = 'bc type "' // trim(word) // '" not recognized in condition,type'
              call printErrMsg(option,string)
          end select
        enddo
      case('TIME','TIMES')
        call fiReadDouble(string,default_time,ierr)
        call fiErrorMsg(option%myrank,'TIME','CONDITION', ierr)   
      case('IPHASE')
        call fiReadInt(string,default_iphase,ierr)
        call fiErrorMsg(option%myrank,'IPHASE','CONDITION', ierr)   
      case('DATUM','DATM')
        call ConditionReadValues(option,word,string,default_datum,word)
      case('GRADIENT','GRAD')
        do
          call fiReadFlotranString(fid,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
          if (fiCheckExit(string)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'keyword','CONDITION,TYPE', ierr)   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
            case('MASS','MASS_RATE')
              sub_condition_ptr => mass_rate
            case('FLUX')
              sub_condition_ptr => flux
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONC','CONCENTRATION')
              sub_condition_ptr => concentration
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
            case default
              call printErrMsg(option,'keyword not recognized in condition,type')
          end select
          call ConditionReadValues(option,word,string,sub_condition_ptr%gradient,word)
          nullify(sub_condition_ptr)
        enddo
      case('TEMPERATURE','TEMP')
        call ConditionReadValues(option,word,string,temperature%dataset, &
                                 temperature%units)
      case('ENTHALPY','H')
        call ConditionReadValues(option,word,string,enthalpy%dataset, &
                                 enthalpy%units)
      case('PRESSURE','PRES','PRESS')
        call ConditionReadValues(option,word,string,pressure%dataset, &
                                 pressure%units)
      case('MASS','MASS_RATE')
        call ConditionReadValues(option,word,string,mass_rate%dataset, &
                                 mass_rate%units)
      case('FLUX','VELOCITY','VEL')
        call ConditionReadValues(option,word,string,pressure%dataset, &
                                 pressure%units)
      case('CONC','CONCENTRATION')
        call ConditionReadValues(option,word,string,concentration%dataset, &
                                 sub_condition_ptr%units)
      case default
        string = 'Keyword: ' // trim(word) // &
                 ' not recognized in flow condition'
        call printErrMsg(option,string)                                 
    end select 
  
  enddo  
  
  ! check whether
  if (default_iphase == 0) then
    call printWrnMsg(option,'"iphase" not set in condition; set to 1')
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

  ! verify the datasets
  word = 'pressure'
  call SubConditionVerify(option,condition,word,pressure,default_time, &
                          default_ctype, default_itype, &
                          default_dataset, &
                          default_datum, default_gradient,PETSC_TRUE)
  word = 'mass_rate'
  call SubConditionVerify(option,condition,word,mass_rate,default_time, &
                          default_ctype, default_itype, &
                          default_dataset, &
                          default_datum, default_gradient,PETSC_TRUE)
  word = 'temperature'
  call SubConditionVerify(option,condition,word,temperature,default_time, &
                          default_ctype, default_itype, &
                          default_dataset, &
                          default_datum, default_gradient,PETSC_TRUE)
  word = 'concentration'
  call SubConditionVerify(option,condition,word,concentration,default_time, &
                          default_ctype, default_itype, &
                          default_dataset, &
                          default_datum, default_gradient,PETSC_TRUE)
  word = 'enthalpy'
  call SubConditionVerify(option,condition,word,enthalpy,default_time, &
                          default_ctype, default_itype, &
                          default_dataset, &
                          default_datum, default_gradient,PETSC_TRUE)
    
  select case(option%iflowmode)
    case(THC_MODE,MPH_MODE)
      if (.not.associated(pressure) .and. .not.associated(mass_rate)) then
        call printErrMsg(option,'pressure and mass_rate condition null in condition: ' // &
                         condition%name)
      endif                         
      if (associated(pressure)) then
        condition%pressure => pressure
      endif                         
      if (associated(mass_rate)) then
        condition%mass_rate => mass_rate
      endif                         
      if (.not.associated(temperature)) then
        call printErrMsg(option,'temperature condition null in condition: ' // &
                         condition%name)
      endif                         
      condition%temperature => temperature
      if (.not.associated(concentration)) then
        call printErrMsg(option,'concentration condition null in condition: ' // &
                         condition%name)
      endif                         
      condition%concentration => concentration
      if (.not.associated(enthalpy)) then
        call printWrnMsg(option,'enthalpy condition null in condition: ' // &
                         condition%name)
      endif                         
      condition%enthalpy => enthalpy
      condition%num_sub_conditions = 4
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 4
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo
      ! must be in this order, which matches the dofs i problem
      if (associated(mass_rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => mass_rate
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
      if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
      
      allocate(condition%itype(FIVE_INTEGER))
      condition%itype = 0
      if (associated(mass_rate)) condition%itype(ONE_INTEGER) = mass_rate%itype
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      condition%itype(TWO_INTEGER) = temperature%itype
      condition%itype(THREE_INTEGER) = concentration%itype
      if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype
      
    case(RICHARDS_MODE)
      if (.not.associated(pressure) .and. .not.associated(mass_rate)) then
        call printErrMsg(option,'pressure and mass_rate condition null in condition: ' // &
                         condition%name)
      endif                         
      if (associated(pressure)) then
        condition%pressure => pressure
      endif                         
      if (associated(mass_rate)) then
        condition%mass_rate => mass_rate
      endif                         
      condition%num_sub_conditions = 1
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure

      allocate(condition%itype(ONE_INTEGER))
      if (associated(mass_rate)) condition%itype(ONE_INTEGER) = mass_rate%itype
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      
      ! these are not used with richards
      if (associated(temperature)) call SubConditionDestroy(temperature)
      if (associated(concentration)) call SubConditionDestroy(concentration)
      if (associated(enthalpy)) call SubConditionDestroy(enthalpy)
      
  end select
  
  call ConditionDatasetDestroy(default_dataset)
  call ConditionDatasetDestroy(default_datum)
  call ConditionDatasetDestroy(default_gradient)
    
  call PetscLogEventEnd(logging%event_flow_condition_read, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine ConditionRead

! ************************************************************************** !
!
! TranConditionRead: Reads a transport condition from the input file
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConditionRead(condition,constraint_list,option,fid)

  use Option_module
  use Fileio_module
  use Logging_module  
  
  implicit none
  
  type(tran_condition_type) :: condition
  type(tran_constraint_list_type) :: constraint_list
  type(option_type) :: option
  PetscInt :: fid
  
  type(tran_constraint_type), pointer :: constraint
  type(tran_constraint_coupler_type), pointer :: constraint_coupler, cur_coupler
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: default_time = 0.d0
  character(len=MAXWORDLENGTH) :: default_time_units
  PetscInt :: default_iphase = 0
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  logical :: found
  PetscInt :: icomp
  PetscInt :: length
  PetscTruth :: minerals_exist
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_tran_condition_read, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC
  default_time_units = ''

  ! read the condition
  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
    if (fiCheckExit(string)) exit  

    call fiReadWord(string,word,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CONDITION', ierr)   
      
    select case(trim(word))
    
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call fiReadWord(string,word,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'INTERPOLATION','CONDITION', ierr)   
        length = len_trim(word)
        call fiCharsToLower(word,length)
        select case(word)
            case('dirichlet')
              condition%itype = DIRICHLET_BC
            case('neumann')
              condition%itype = NEUMANN_BC
            case('mole','mole_rate')
              condition%itype = MASS_RATE_SS
            case('zero_gradient')
              condition%itype = ZERO_GRADIENT_BC
            case default
              call printErrMsg(option,'keyword not recognized in condition,type')
        end select
      case('TIME')
        call fiReadDouble(string,default_time,ierr)
        call fiErrorMsg(option%myrank,'TIME','CONDITION', ierr) 
      case('UNITS') 
        call fiReadWord(string,word,PETSC_TRUE,ierr) 
        call fiErrorMsg(option%myrank,'UNITS','CONDITION', ierr)   
        call fiWordToLower(word)
        select case(trim(word))     
          case('s','sec','min','hr','d','day','y','yr')
            default_time_units = trim(word)         
        end select          
      case('CONSTRAINT_LIST')
        do
          call fiReadFlotranString(fid,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONSTRAINT',ierr)
              
          if (fiCheckExit(string)) exit  

          constraint_coupler => TranConstraintCouplerCreate(option)
          call fiReadDouble(string,constraint_coupler%time,ierr)
          call fiErrorMsg(option%myrank,'time','CONSTRAINT_LIST', ierr) 
          ! time units are optional  
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'constraint name','CONSTRAINT_LIST', ierr) 
          ! read constraint name
          call fiReadWord(string,constraint_coupler%constraint_name,PETSC_TRUE,ierr)
          if (ierr /= 0) then
            constraint_coupler%time_units = default_time_units
            constraint_coupler%constraint_name = trim(word)
          else
            constraint_coupler%time_units = word
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
        call fiReadWord(string,constraint%name,PETSC_TRUE,ierr)
        call fiErrorMsg(option%myrank,'constraint','name',ierr) 
        call printMsg(option,constraint%name)
        call TranConstraintRead(constraint,option,fid)
        call TranConstraintAddToList(constraint,constraint_list)
        constraint_coupler%aqueous_species => constraint%aqueous_species
        constraint_coupler%minerals => constraint%minerals
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
        string = 'Keyword: ' // trim(word) // &
                 ' not recognized in transport condition'
        call printErrMsg(option,string)
    end select 
  
  enddo  
  
  call PetscLogEventEnd(logging%event_tran_condition_read, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine TranConditionRead

! ************************************************************************** !
!
! TranConstraintRead: Reads a transport constraint from the input file
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine TranConstraintRead(constraint,option,fid)

  use Option_module
  use Fileio_module
  use Logging_module  
  
  implicit none
  
  type(tran_constraint_type) :: constraint
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: length
  PetscInt :: icomp
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_tran_constraint_read, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  ! read the constraint
  ierr = 0
  do
  
    call fiReadFlotranString(fid,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'CONSTRAINT',ierr)
        
    if (fiCheckExit(string)) exit  

    call fiReadWord(string,word,PETSC_TRUE,ierr)
    call fiErrorMsg(option%myrank,'keyword','CONSTRAINT', ierr)   
      
    select case(trim(word))

      case('CONC','CONCENTRATIONS')

        aq_species_constraint => AqueousSpeciesConstraintCreate(option)

        icomp = 0
        do
          call fiReadFlotranString(fid,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONSTRAINT,CONCENTRATIONS',ierr)
          
          if (fiCheckExit(string)) exit  
          
          icomp = icomp + 1        
          
          call fiReadWord(string,aq_species_constraint%names(icomp), &
                          PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'aqueous species name', &
                          'CONSTRAINT,CONCENTRATIONS', ierr) 
          call printMsg(option,trim(aq_species_constraint%names(icomp)))
          call fiReadDouble(string,aq_species_constraint%conc(icomp),ierr)
          call fiErrorMsg(option%myrank,'concentration', &
                          'CONSTRAINT,CONCENTRATIONS', ierr)          
          call fiReadWord(string,word,PETSC_TRUE,ierr)
          call fiDefaultMsg(option%myrank, &
                            'CONSTRAINT,CONCENTRATION,constraint_type', ierr)
          length = len_trim(word)
          if (length > 0) then
            call fiCharsToUpper(word,length)
            select case(word)
              case('F','FREE')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_FREE
              case('T','Total')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
              case('P')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_P
              case('MINERAL','MNRL') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_MINERAL
              case('GAS') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_GAS
              case default
                string = 'Keyword: ' // trim(word) // &
                         ' not recognized in constraint,concentration'
                call printErrMsg(option,string)
            end select 
            if (aq_species_constraint%constraint_type(icomp) == CONSTRAINT_MINERAL .or. &
                aq_species_constraint%constraint_type(icomp) == CONSTRAINT_GAS) then
             call fiReadWord(string,aq_species_constraint%constraint_spec_name(icomp), &
                             PETSC_FALSE,ierr)
              call fiErrorMsg(option%myrank,'constraint name', &
                              'CONSTRAINT,CONCENTRATIONS', ierr) 
            endif
          else
            aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
          endif  
        
        enddo  
        
        if (associated(constraint%aqueous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint 
        
      case('MNRL','MINERALS')

        mineral_constraint => MineralConstraintCreate(option)

        icomp = 0
        do
          call fiReadFlotranString(fid,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONSTRAINT,MINERALS',ierr)
          
          if (fiCheckExit(string)) exit          
          
          icomp = icomp + 1
          
          call fiReadWord(string,mineral_constraint%names(icomp), &
                          PETSC_TRUE,ierr)
          call fiErrorMsg(option%myrank,'mineral name', &
                          'CONSTRAINT,CONCENTRATIONS', ierr)  
          call printMsg(option,trim(mineral_constraint%names(icomp)))
          call fiReadDouble(string,mineral_constraint%conc(icomp),ierr)
          call fiErrorMsg(option%myrank,'concentration', &
                          'CONSTRAINT,CONCENTRATIONS', ierr)          
        
        enddo  
        
        if (associated(constraint%minerals)) then
          call MineralConstraintDestroy(constraint%minerals)
        endif
        constraint%minerals => mineral_constraint 
                            
      case default
        string = 'Keyword: ' // trim(word) // &
                 ' not recognized in transport constraint'
        call printErrMsg(option,string)
    end select 
  
  enddo  
  
  call PetscLogEventEnd(logging%event_tran_constraint_read, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine TranConstraintRead

! ************************************************************************** !
!
! ConditionReadValues: Read the value(s) of a condition variable
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine ConditionReadValues(option,keyword,string,dataset,units)

  use Fileio_module
  use Option_module
  use Logging_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword
  type(flow_condition_dataset_type) :: dataset
  character(len=MAXWORDLENGTH) :: units
  
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, irank
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read_values, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
  ierr = 0
  string2 = trim(string)
  call fiReadWord(string,word,PETSC_TRUE,ierr)
  call fiErrorMsg(option%myrank,'file or value','CONDITION', ierr)
  length = len_trim(word)
  call fiCharsToLower(word,length)
  if (fiStringCompare(word,'file',FOUR_INTEGER)) then
    call fiReadWord(string,word,PETSC_TRUE,ierr)
    error_string = keyword // ' FILE'
    call fiErrorMsg(option%myrank,error_string,'CONDITION', ierr)
    call ConditionReadValuesFromFile(word,dataset,option)
  else
    string = trim(string2)
    allocate(dataset%values(dataset%rank,1))
    do irank=1,dataset%rank
      call fiReadDouble(string,dataset%values(irank,1),ierr)
      write(error_string,*) trim(keyword) // ' dataset_values, irank = ', irank
      call fiErrorMsg(option%myrank,error_string,'CONDITION', ierr) 
    enddo
  endif
  call fiReadWord(string,word,PETSC_TRUE,ierr)
  if (ierr /= 0) then
    word = trim(keyword) // ' UNITS'
    call fiDefaultMsg(option%myrank,word, ierr)
  else
    units = trim(word)
  endif

  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

end subroutine ConditionReadValues

! ************************************************************************** !
!
! ConditionReadValuesFromFile: Read values from a external file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine ConditionReadValuesFromFile(filename,dataset,option)

  use Fileio_module
  use Utility_module
  use Option_module

  implicit none
  
  type(option_type) :: option
  type(flow_condition_dataset_type) :: dataset
  character(len=MAXWORDLENGTH) :: filename
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: temp_times(:), temp_array1(:), temp_array2(:), &
                        temp_array3(:)
  PetscReal :: temp_time
  PetscInt :: max_size = 1000
  PetscInt :: temp_max_size
  PetscInt :: fid
  PetscInt :: count, i, status
  PetscErrorCode :: ierr

  fid = 86
  open(unit=fid,file=filename,status="old",iostat=status)
  if (status /= 0) then
    print *, 'file: ', trim(filename), ' not found'
    stop
  endif
  
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
  do
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit
    count = count + 1
    call fiReadDouble(string,temp_times(count),ierr)
    call fiErrorMsg(option%myrank,'time','CONDITION FILE', ierr)   
    call fiReadDouble(string,temp_array1(count),ierr)
    call fiErrorMsg(option%myrank,'array1','CONDITION FILE', ierr)
!    print *, 'RCF:', temp_times(count),  temp_array1(count)
!   I have commented out the above line because it creates an impossible 
!   amount of console output for even moderately-sized parallel runs!
!   If someone needs this data, we need to add a -print_debug flag or 
!   some such thing.  --RTM
    if (dataset%rank > 1) then
      call fiReadDouble(string,temp_array2(count),ierr)
      call fiErrorMsg(option%myrank,'array2','CONDITION FILE', ierr) 
    endif
    if (dataset%rank > 2) then
      call fiReadDouble(string,temp_array3(count),ierr)
      call fiErrorMsg(option%myrank,'array3','CONDITION FILE', ierr) 
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
    if (count /= size(dataset%times,1)) then
      print *, 'Number of times (', count, ') in ', trim(filename), &
               ' does not match previous allocation: ', size(dataset%times,1)
      stop
    endif
    do i=1,count
      if (dabs(dataset%times(i)-temp_times(i)) > 1.d-8) then
        print *, 'Time (', temp_times(i), ') in ', trim(filename), &
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
  
  close(fid)

end subroutine ConditionReadValuesFromFile

! ************************************************************************** !
!
! ConditionUpdate: Updates a transient condition
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine ConditionUpdate(condition_list,option,time)

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
        call SubConditionUpdateDataset(option,time,sub_condition%dataset)
        call SubConditionUpdateDataset(option,time,sub_condition%datum)
        call SubConditionUpdateDataset(option,time,sub_condition%gradient)
      endif
      
    enddo
      
    condition => condition%next
    
  enddo
  
end subroutine ConditionUpdate

! ************************************************************************** !
!
! SubConditionUpdateDataset: Updates a transient condition dataset
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine SubConditionUpdateDataset(option,time,dataset)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: time
  logical :: is_cyclic
  PetscInt :: interpolation_method
  type(flow_condition_dataset_type) :: dataset
  
  PetscInt :: irank
  PetscInt :: cur_time_index
  PetscInt :: next_time_index
  PetscReal :: time_fraction

 ! potentially for initial condition
  if (time < 1.d-40 .and. .not.dataset%is_transient) return  

  ! cycle times if at max_time_index and cyclic
  if (dataset%cur_time_index == dataset%max_time_index .and. &
      dataset%is_cyclic .and. dataset%max_time_index > 1) then
      
    do cur_time_index = 1, dataset%max_time_index
      dataset%times(cur_time_index) = dataset%times(cur_time_index) + &     
                               dataset%times(dataset%max_time_index)
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
  
end subroutine SubConditionUpdateDataset

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
! ConditionInitList: Initializes a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine ConditionInitList(list)

  implicit none

  type(condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine ConditionInitList

! ************************************************************************** !
!
! ConditionAddToList: Adds a new condition to a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine ConditionAddToList(new_condition,list)

  implicit none
  
  type(flow_condition_type), pointer :: new_condition
  type(condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine ConditionAddToList

! ************************************************************************** !
!
! ConditionGetPtrFromList: Returns a pointer to the condition matching &
!                          condition_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function ConditionGetPtrFromList(condition_name,condition_list)

  use Fileio_module

  implicit none
  
  type(flow_condition_type), pointer :: ConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list
 
  PetscInt :: length
  type(flow_condition_type), pointer :: condition
    
  nullify(ConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        fiStringCompare(condition%name,condition_name, &
                        length)) then
      ConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function ConditionGetPtrFromList

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

  use Fileio_module

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
        fiStringCompare(condition%name,condition_name, &
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

  use Fileio_module

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
        fiStringCompare(constraint%name,constraint_name, &
                        length)) then
      TranConstraintGetPtrFromList => constraint
      return
    endif
    constraint => constraint%next
  enddo
  
end function TranConstraintGetPtrFromList

! ************************************************************************** !
!
! ConditionDestroyList: Deallocates a list of conditions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine ConditionDestroyList(condition_list)

  implicit none
  
  type(condition_list_type), pointer :: condition_list
  
  type(flow_condition_type), pointer :: condition, prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call ConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine ConditionDestroyList

! ************************************************************************** !
!
! ConditionDestroy: Deallocates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine ConditionDestroy(condition)

  implicit none
  
  type(flow_condition_type), pointer :: condition
  
  PetscInt :: i
  
  if (.not.associated(condition)) return
  
  if (associated(condition%sub_condition_ptr)) then
    do i=1,condition%num_sub_conditions
      call SubConditionDestroy(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  
  nullify(condition%pressure)
  nullify(condition%mass_rate)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine ConditionDestroy

! ************************************************************************** !
!
! SubConditionDestroy: Destroys a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine SubConditionDestroy(sub_condition)

  implicit none
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  if (.not.associated(sub_condition)) return
  
  call ConditionDatasetDestroy(sub_condition%dataset)
  call ConditionDatasetDestroy(sub_condition%datum)
  call ConditionDatasetDestroy(sub_condition%gradient)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine SubConditionDestroy

! ************************************************************************** !
!
! ConditionDatasetDestroy: Destroys a dataset associated with a sub_condition
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine ConditionDatasetDestroy(dataset)

  implicit none
  
  type(flow_condition_dataset_type) :: dataset
  
  if (associated(dataset%times)) deallocate(dataset%times)
  nullify(dataset%times)
  if (associated(dataset%values)) deallocate(dataset%values)
  nullify(dataset%values)  

end subroutine ConditionDatasetDestroy

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
    nullify(prev_coupler%aqueous_species)
    nullify(prev_coupler%minerals)
    nullify(prev_coupler%next)
    deallocate(prev_coupler)
    nullify(prev_coupler)
  enddo
  
  nullify(coupler_list)
  
end subroutine TranConstraintCouplerDestroy

end module Condition_module
