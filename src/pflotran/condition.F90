module Condition_module
 
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2

  type, public :: condition_dataset_type
    PetscInt :: rank
    logical :: is_transient
    logical :: is_cyclic
    PetscInt :: interpolation_method
    PetscReal, pointer :: times(:)
    PetscReal, pointer :: values(:,:)
    PetscReal, pointer :: cur_value(:)
    PetscInt :: cur_time_index
    PetscInt :: max_time_index
  end type condition_dataset_type
  
  type, public :: condition_type
    PetscInt :: id                                 ! id from which condition can be referenced
    character(len=MAXWORDLENGTH) :: class         ! character string describing class of condition
    PetscInt :: iclass                            ! integer id for class
    logical :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(sub_condition_type), pointer :: pressure
    type(sub_condition_type), pointer :: mass_rate
    type(sub_condition_type), pointer :: temperature
    type(sub_condition_type), pointer :: concentration
    type(sub_condition_ptr_type), pointer :: transport_concentrations(:)
    type(sub_condition_type), pointer :: enthalpy
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type condition_type
  
  type, public :: sub_condition_type
    PetscInt :: itype                  ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    type(condition_dataset_type) :: datum
    type(condition_dataset_type) :: gradient
    type(condition_dataset_type) :: dataset

  end type sub_condition_type
  
  type, public :: sub_condition_ptr_type
    type(sub_condition_type), pointer :: ptr
  end type sub_condition_ptr_type
    
  type, public :: condition_ptr_type
    type(condition_type), pointer :: ptr
  end type condition_ptr_type
  
  type, public :: condition_list_type
    PetscInt :: num_conditions
    type(condition_type), pointer :: first
    type(condition_type), pointer :: last
    type(condition_ptr_type), pointer :: array(:)    
  end type condition_list_type
  
  PetscInt, save :: condition_count = 0
  
  public :: ConditionCreate, ConditionDestroy, ConditionRead, &
            ConditionAddToList, ConditionInitList, ConditionDestroyList, &
            ConditionGetPtrFromList, ConditionUpdate
    
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
  type(condition_type), pointer :: ConditionCreate
  
  type(condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%mass_rate)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%transport_concentrations)
  nullify(condition%sub_condition_ptr)
  nullify(condition%itype)
  nullify(condition%next)
  condition%sync_time_with_update = .false.
  condition%time_units = ""
  condition%length_units = ""
  condition%id = 0
  condition%iphase = 0
  condition%num_sub_conditions = 0
  condition%class = ""
  condition%iclass = NULL_CLASS
  condition%name = ""

  condition_count = condition_count + 1
  condition%id = condition_count
  
  ConditionCreate => condition

end function ConditionCreate

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
  
  type(sub_condition_type), pointer :: SubConditionCreate
  
  PetscInt :: ndof
  
  type(sub_condition_type), pointer :: sub_condition
  
  allocate(sub_condition)
  sub_condition%units = ""
  sub_condition%itype = 0
  sub_condition%ctype = ""
  sub_condition%name = ""

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
  
  type(sub_condition_type), pointer :: GetSubConditionFromArrayByName
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
  
  type(condition_dataset_type) :: dataset

  nullify(dataset%times)
  nullify(dataset%values)
  nullify(dataset%cur_value)
  dataset%cur_time_index = 0
  dataset%max_time_index = 0
  dataset%rank = 0
  dataset%is_cyclic = .false.
  dataset%is_transient = .false.
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
                              default_datum, default_gradient)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(sub_condition_type), pointer :: sub_condition
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscTruth :: default_cyclic
  PetscInt :: default_interpolation
  PetscReal :: default_time
  PetscInt :: default_iphase
  type(condition_dataset_type) :: default_dataset
  type(condition_dataset_type) :: default_datum
  type(condition_dataset_type) :: default_gradient

  PetscInt :: array_size

  if (.not.associated(sub_condition)) return
  
  if (.not.associated(sub_condition%dataset%values)) then
    call SubConditionDestroy(sub_condition)
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
  type(condition_dataset_type) :: dataset
  type(condition_dataset_type) :: default_dataset
  
  PetscInt :: array_size
  
  if (default_dataset%is_cyclic) dataset%is_cyclic = .true.
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
    dataset%is_transient = .true.
  else
    dataset%is_transient = .false.
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
  
  type(condition_type) :: condition
  type(option_type) :: option
  PetscInt :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  type(sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, mass_rate, &
                                       sub_condition_ptr
  type(sub_condition_ptr_type), pointer :: transport_concentrations(:)
  PetscReal :: default_time = 0.d0
  PetscInt :: default_iphase = 0
  type(condition_dataset_type) :: default_dataset
  type(condition_dataset_type) :: default_datum
  type(condition_dataset_type) :: default_gradient
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, length, idof
  logical :: found
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_condition_read, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  call ConditionDatasetInit(default_dataset)
  default_dataset%rank = 1
  default_dataset%interpolation_method = STEP
  default_dataset%is_cyclic = .false.
  call ConditionDatasetInit(default_datum)
  default_datum%rank = 3
  call ConditionDatasetInit(default_gradient)
  default_gradient%rank = 3
  
  pressure => SubConditionCreate(option%nphase)
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
  
  if (option%ntrandof > 0) then
    allocate(transport_concentrations(option%ntrandof))
    do idof = 1, option%ntrandof
      transport_concentrations(idof)%ptr => SubConditionCreate(ONE_INTEGER)
      transport_concentrations(idof)%ptr%name = option%comp_names(idof)
      transport_concentrations(idof)%ptr%units = 'M'
    enddo
  else
    nullify(transport_concentrations)
  endif
  
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','CONDITION', ierr)   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call fiReadWord(string,word,.true.,ierr)
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
              if (condition%iclass == TRANSPORT_CLASS) then
                do idof = 1, option%ntrandof
                  transport_concentrations(idof)%ptr%units = trim(word)
                enddo
              else
                concentration%units = trim(word)
              endif
            case('KJ/mol')
              enthalpy%units = trim(word)
          end select
        enddo
      case('CLASS') ! read condition class (flow vs. transport)
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(option%myrank,'CLASS','CONDITION', ierr)   
        length = len_trim(word)
        call fiCharsToLower(word,length)
        condition%class = word
        select case(word)
          case('flow')
            condition%iclass = FLOW_CLASS
          case('tran','transport')
            condition%iclass = TRANSPORT_CLASS
          case default
            call printErrMsg(option,'class: '//word//' not recognized in condition')
        end select
      case('CYCLIC')
        default_dataset%is_cyclic = .true.
      case('INTERPOLATION')
        call fiReadWord(string,word,.true.,ierr)
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
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,.true.,ierr)
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
              if (condition%iclass == TRANSPORT_CLASS) then
                call fiReadWord(string,word,.true.,ierr)
                call fiErrorMsg(option%myrank,'name','CONDITION,CONCENTRATION', ierr)
                nullify(sub_condition_ptr)
                if (option%ntrandof > 0) then
                  sub_condition_ptr => &
                    GetSubConditionFromArrayByName(transport_concentrations,word)
                  if (.not.associated(sub_condition_ptr)) then
                    string = 'solute name "' // trim(word) // &
                             '" not recognized in condition'
                    call printErrMsg(option,string)
                  endif
                else
                  sub_condition_ptr => concentration
                endif
              else            
                sub_condition_ptr => concentration
              endif
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
            case default
              call printErrMsg(option,'keyword not recognized in condition,type')
          end select
          call fiReadWord(string,word,.true.,ierr)
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
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
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
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg(option%myrank,'CONDITION',ierr)
          
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',THREE_INTEGER)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,.true.,ierr)
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
              if (condition%iclass == TRANSPORT_CLASS) then
                call fiReadWord(string,word,.true.,ierr)
                call fiErrorMsg(option%myrank,'name','CONDITION,CONCENTRATION', ierr)
                nullify(sub_condition_ptr)
                if (option%ntrandof > 0) then                
                  sub_condition_ptr => &
                    GetSubConditionFromArrayByName(transport_concentrations,word)
                  if (.not.associated(sub_condition_ptr)) then
                    string = 'solute name "' // trim(word) // &
                             '" not recognized in condition'
                    call printErrMsg(option,string)
                  endif
                else
                  sub_condition_ptr => concentration
                endif
              else            
                sub_condition_ptr => concentration
              endif
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
        if (condition%iclass == NULL_CLASS) condition%iclass = FLOW_CLASS
        call ConditionReadValues(option,word,string,pressure%dataset, &
                                 pressure%units)
      case('MASS','MASS_RATE')
        if (condition%iclass == NULL_CLASS) condition%iclass = FLOW_CLASS
        call ConditionReadValues(option,word,string,mass_rate%dataset, &
                                 mass_rate%units)
      case('FLUX','VELOCITY','VEL')
        call ConditionReadValues(option,word,string,pressure%dataset, &
                                 pressure%units)
      case('CONC','CONCENTRATION')
        if (condition%iclass == TRANSPORT_CLASS) then
          call fiReadWord(string,word,.true.,ierr)
          call fiErrorMsg(option%myrank,'name','CONDITION,CONCENTRATION', ierr)
          nullify(sub_condition_ptr)
          if (option%ntrandof > 0) then           
            sub_condition_ptr => &
              GetSubConditionFromArrayByName(transport_concentrations,word)
            if (.not.associated(sub_condition_ptr)) then
              string = 'solute name "' // trim(word) // &
                       '" not recognized in condition'
              call printErrMsg(option,string)
            endif
          else
            sub_condition_ptr => concentration
          endif
        else
          sub_condition_ptr => concentration
        endif
        call ConditionReadValues(option,word,string,sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units)
    end select 
  
  enddo  
  
  ! check to ensure that class and type have been set
  if (len_trim(condition%class) < 1) then
    call printErrMsg(option,'"class" not set in condition')
  endif

  ! check whether
  if (default_iphase == 0) then
    call printWrnMsg(option,'"iphase" not set in condition; set to 1')
    condition%iphase = 1
  else
    condition%iphase = default_iphase    
  endif
  
  ! update datum and gradient defaults, if null, based on dataset default
  if (default_dataset%is_cyclic) then
    default_datum%is_cyclic = .true.
    default_gradient%is_cyclic = .true.
  endif
  if (default_datum%interpolation_method == NULL) &
    default_datum%interpolation_method = default_dataset%interpolation_method
  if (default_gradient%interpolation_method == NULL) &
    default_gradient%interpolation_method = default_dataset%interpolation_method

  ! verify the datasets
  if (condition%iclass == FLOW_CLASS) then
    word = 'pressure'
    call SubConditionVerify(option,condition,word,pressure,default_time, &
                            default_ctype, default_itype, &
                            default_dataset, &
                            default_datum, default_gradient)
    word = 'mass_rate'
    call SubConditionVerify(option,condition,word,mass_rate,default_time, &
                            default_ctype, default_itype, &
                            default_dataset, &
                            default_datum, default_gradient)
    word = 'temperature'
    call SubConditionVerify(option,condition,word,temperature,default_time, &
                            default_ctype, default_itype, &
                            default_dataset, &
                            default_datum, default_gradient)
    word = 'concentration'
    call SubConditionVerify(option,condition,word,concentration,default_time, &
                            default_ctype, default_itype, &
                            default_dataset, &
                            default_datum, default_gradient)
    word = 'enthalpy'
    call SubConditionVerify(option,condition,word,enthalpy,default_time, &
                            default_ctype, default_itype, &
                            default_dataset, &
                            default_datum, default_gradient)
    ! these are not used with transport
    do idof = 1, option%ntrandof
      if (associated(transport_concentrations(idof)%ptr)) &
        call SubConditionDestroy(transport_concentrations(idof)%ptr)
    enddo
    if (associated(transport_concentrations)) deallocate(transport_concentrations)
    nullify(transport_concentrations)
  else
    do idof = 1, option%ntrandof
      word = 'solute concentration: ' // trim(transport_concentrations(idof)%ptr%name)
      call SubConditionVerify(option,condition,word,transport_concentrations(idof)%ptr,default_time, &
                              default_ctype, default_itype, &
                              default_dataset, &
                              default_datum, default_gradient)
    enddo
    ! these are not used with transport
    if (associated(pressure)) call SubConditionDestroy(pressure)
    if (associated(mass_rate)) call SubConditionDestroy(mass_rate)
    if (associated(temperature)) call SubConditionDestroy(temperature)
    if (associated(concentration)) call SubConditionDestroy(concentration)
    if (associated(enthalpy)) call SubConditionDestroy(enthalpy)    
  endif
    
  if (condition%iclass == FLOW_CLASS) then
    select case(option%iflowmode)
      case(RICHARDS_MODE,MPH_MODE)
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
        
      case(RICHARDS_LITE_MODE)
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
        
        ! these are not used with richards_lite
        if (associated(temperature)) call SubConditionDestroy(temperature)
        if (associated(concentration)) call SubConditionDestroy(concentration)
        if (associated(enthalpy)) call SubConditionDestroy(enthalpy)
        
    end select
  else
    condition%num_sub_conditions = option%ntrandof
    allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
    allocate(condition%itype(condition%num_sub_conditions))
    do idof = 1, option%ntrandof
      condition%sub_condition_ptr(idof)%ptr => transport_concentrations(idof)%ptr
      condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
    enddo
  endif
  
  call ConditionDatasetDestroy(default_dataset)
  call ConditionDatasetDestroy(default_datum)
  call ConditionDatasetDestroy(default_gradient)
    
  call PetscLogEventEnd(logging%event_condition_read, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine ConditionRead

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
  type(condition_dataset_type) :: dataset
  character(len=MAXWORDLENGTH) :: units
  
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, irank
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_condition_read_values, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    
  ierr = 0
  string2 = trim(string)
  call fiReadWord(string,word,.true.,ierr)
  call fiErrorMsg(option%myrank,'file or value','CONDITION', ierr)
  length = len_trim(word)
  call fiCharsToLower(word,length)
  if (fiStringCompare(word,'file',FOUR_INTEGER)) then
    call fiReadWord(string,word,.true.,ierr)
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
  call fiReadWord(string,word,.true.,ierr)
  if (ierr /= 0) then
    word = trim(keyword) // ' UNITS'
    call fiDefaultMsg(option%myrank,word, ierr)
  else
    units = trim(word)
  endif

  call PetscLogEventEnd(logging%event_condition_read_values, &
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
  type(condition_dataset_type) :: dataset
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
  print *, 'Read condition from file:',filename
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
subroutine ConditionUpdate(condition_list,option,time,iclass)

  use Option_module
  
  implicit none
  
  type(condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  PetscInt :: iclass
  
  type(condition_type), pointer :: condition
  type(sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition  
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    if (iclass == NULL_CLASS .or. &
        iclass == condition%iclass) then
    
      do isub_condition = 1, condition%num_sub_conditions

        sub_condition => condition%sub_condition_ptr(isub_condition)%ptr
        
        if (associated(sub_condition)) then
          call SubConditionUpdateDataset(option,time,sub_condition%dataset)
          call SubConditionUpdateDataset(option,time,sub_condition%datum)
          call SubConditionUpdateDataset(option,time,sub_condition%gradient)
        endif
        
      enddo
      
    endif
        
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
  type(condition_dataset_type) :: dataset
  
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
  
  type(condition_type), pointer :: new_condition
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
  
  type(condition_type), pointer :: ConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list
 
  PetscInt :: length
  type(condition_type), pointer :: condition
    
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
! ConditionDestroyList: Deallocates a list of conditions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine ConditionDestroyList(condition_list)

  implicit none
  
  type(condition_list_type), pointer :: condition_list
  
  type(condition_type), pointer :: condition, prev_condition
  
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
  
  type(condition_type), pointer :: condition
  
  PetscInt :: i
  
  if (.not.associated(condition)) return
  
  do i=1,condition%num_sub_conditions
    call SubConditionDestroy(condition%sub_condition_ptr(i)%ptr)
  enddo
  deallocate(condition%sub_condition_ptr)
  
  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  
  nullify(condition%pressure)
  nullify(condition%mass_rate)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  
  nullify(condition%sub_condition_ptr)
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
  
  type(sub_condition_type), pointer :: sub_condition
  
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
  
  type(condition_dataset_type) :: dataset
  
  if (associated(dataset%times)) deallocate(dataset%times)
  nullify(dataset%times)
  if (associated(dataset%values)) deallocate(dataset%values)
  nullify(dataset%values)  

end subroutine ConditionDatasetDestroy

end module Condition_module
