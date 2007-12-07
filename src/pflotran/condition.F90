module Condition_module
 
  implicit none

  private
  
#include "definitions.h"

!#define X_DIRECTION 1 ! now in definitions.h
!#define Y_DIRECTION 2
!#define Z_DIRECTION 3

#define STEP 0
#define LINEAR 1
 
  type, public :: condition_type
    integer :: id                                 ! id from which condition can be referenced
    integer, pointer :: itype(:)                  ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: class         ! character string describing class of condition
    character(len=MAXWORDLENGTH), pointer :: ctype(:) ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    character(len=MAXWORDLENGTH), pointer :: units(:)      ! units
      ! units(1) = time
      ! units(2) = length
      ! units(3:ndof) = dofs in problem
    integer :: num_values                         ! number of entries in the arrays of values
    integer :: num_dof                            ! number of degrees of freedom in condtion
    logical :: is_cyclic                          ! cycles after last time
    integer :: interpolation_method               ! method of interplating condition based on time
    integer :: iphase
    logical :: is_transient                       ! tells whether condition will change over time
    real*8, pointer :: times(:)                   ! array of times between which linear interpolation of values occurs
    real*8, pointer :: values(:,:)                ! array of condition values, size(ndof,max_time_index)
    real*8, pointer :: cur_value(:)               ! current value of condition a time t, size(ndof)
    real*8 :: datum(3)                            ! location of reference value(s) in domain
    real*8, pointer :: gradient(:,:)              ! rate at which reference value(s) change(s) over 3D space
    integer :: cur_time_index, max_time_index     ! current and maximum time index in arrays
    type(condition_type), pointer :: next         ! pointer to next condition_type for linked-lists
  end type condition_type
  
  type, public :: condition_ptr_type
    type(condition_type), pointer :: ptr
  end type condition_ptr_type
  
  type, public :: condition_list_type
    integer :: num_conditions
    type(condition_type), pointer :: first
    type(condition_type), pointer :: last
    type(condition_ptr_type), pointer :: array(:)    
  end type condition_list_type
  
  integer, save :: condition_count = 0
  
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
  nullify(condition%units)
  nullify(condition%times)
  nullify(condition%values)
  nullify(condition%cur_value)
  nullify(condition%next)
  condition%id = 0
  condition%iphase = 0
  condition%num_values = 0
  condition%num_dof = 0
  condition%is_cyclic = .false.
  condition%is_transient = .false.
!  condition%interpolation_method = LINEAR
  condition%interpolation_method = STEP ! default to step for src/sinks
  nullify(condition%itype)
  condition%class = ""
  nullify(condition%ctype)
  condition%name = ""
  condition%datum = 0.d0
  condition%cur_time_index = 0
  condition%max_time_index = 0
  
  condition_count = condition_count + 1
  condition%id = condition_count
  
  allocate(condition%units(option%ndof+2))
  select case(option%imode)
    case(RICHARDS_MODE,MPH_MODE)
      condition%units(1) = 'yr'
      condition%units(2) = 'm'
      condition%units(3) = 'Pa'
      condition%units(4) = 'C'
      condition%units(5) = 'M'
      allocate(condition%gradient(3,option%ndof))
  end select
  
  condition%gradient = 0.d0

  ConditionCreate => condition

end function ConditionCreate

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
  
  implicit none
  
  type(condition_type) :: condition
  type(option_type) :: option
  integer :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  real*8, pointer :: pressure(:), flux(:), temperature(:), &
                     concentration(:), times(:), enthalpy(:)
  character(len=MAXWORDLENGTH), pointer :: units(:), ctype(:)
  integer, pointer :: itype(:)
  integer :: time_index, length_index, pres_index, temp_index, conc_index, iphase
  integer :: enthalpy_index
  integer :: max_size, index, idof, ierr, max_index, max_dof, max_required_dof

  nullify(times)
  nullify(pressure)
  nullify(temperature)
  nullify(flux)
  nullify(concentration)
  nullify(enthalpy)
  
  ! need to set up indices for potential data contained in condition
  select case(option%imode)
    case(RICHARDS_MODE,MPH_MODE)
      pres_index = 1
      temp_index = 2
      conc_index = 3
    max_required_dof = 3
      enthalpy_index = 4
    max_dof = 4
      time_index = 5
      length_index = 6
      max_index = 6
    case default
      pres_index = 1
      temp_index = 2
      conc_index = 3
    max_required_dof = 3
    max_dof = 3
      max_dof = 4
      max_index = 6  
  end select
  allocate(units(max_index))
  allocate(itype(max_dof))
  allocate(ctype(max_dof))  
  units = ""
  itype = 0                
  ctype = ""
  iphase = 0

  ! read the condition
  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    call fiReadStringErrorMsg('CONDITION',ierr)
          
    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',3)) exit  

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg('keyword','CONDITION', ierr)   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%units(time_index) = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%units(length_index) = trim(word)
            case('Pa','KPa')
              condition%units(pres_index) = trim(word)
            case('C','K')
              condition%units(temp_index) = trim(word)
            case('M','mol/L')
              condition%units(conc_index) = trim(word)
            case('KJ/mol')
              condition%units(enthalpy_index) = trim(word)
          end select
        enddo
      case('CLASS') ! read condition class (flow vs. transport)
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg('CLASS','CONDITION', ierr)   
        call fiCharsToLower(word,len_trim(word))
        condition%class = word
      case('CYCLIC') ! read condition class (flow vs. transport)
        condition%is_cyclic = .true.
      case('INTERPOLATION') ! read condition class (flow vs. transport)
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg('INTERPOLATION','CONDITION', ierr)   
        call fiCharsToLower(word,len_trim(word))
        select case(word)
          case('step')
            condition%interpolation_method = STEP
          case('linear') 
            condition%interpolation_method = LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('CONDITION',ierr)
          
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',3)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,.true.,ierr)
          call fiErrorMsg('keyword','CONDITION,TYPE', ierr)   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              index = pres_index
            case('TEMP','TEMPERATURE')
              index = temp_index
            case('CONC','CONCENTRATION')
              index = conc_index
            case('H','ENTHALPY')
              index = enthalpy_index
            case default
              call printErrMsg(option,'index not recognized in condition,type')
          end select
          call fiReadWord(string,word,.true.,ierr)
          call fiErrorMsg('TYPE','CONDITION', ierr)   
          call fiCharsToLower(word,len_trim(word))
          ctype(index) = word
          select case(word)
            case('dirichlet')
              itype(index) = DIRICHLET_BC
            case('neumann')
              itype(index) = NEUMANN_BC
            case('mass')
              itype(index) = MASS_RATE
            case('hydrostatic','hydro','hydrostat','static')
              itype(index) = HYDROSTATIC_BC
            case('zero_gradient')
              itype(index) = ZERO_GRADIENT_BC
            case default
              call printErrMsg(option,'bc type not recognized in condition,type')
          end select
        enddo
      case('TIME','TIMES')
        if (.not.associated(times)) then
          allocate(times(1))
        endif
        call fiReadDouble(string,times(1),ierr)
        call fiErrorMsg('TIME','CONDITION', ierr)   
      case('IPHASE')
        call fiReadInt(string,iphase,ierr)
        call fiErrorMsg('IPHASE','CONDITION', ierr)   
      case('DATUM','DATM')
        call fiReadDouble(string,condition%datum(X_DIRECTION),ierr)
        call fiErrorMsg('X Datum','CONDITION', ierr)   
        call fiReadDouble(string,condition%datum(Y_DIRECTION),ierr)
        call fiErrorMsg('Y Datum','CONDITION', ierr)   
        call fiReadDouble(string,condition%datum(Z_DIRECTION),ierr)
        call fiErrorMsg('Z Datum','CONDITION', ierr)   
      case('GRADIENT','GRAD')
        do
          call fiReadFlotranString(IUNIT1,string,ierr)
          call fiReadStringErrorMsg('CONDITION',ierr)
          
          if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
              fiStringCompare(string,'END',3)) exit          
          
          if (ierr /= 0) exit
          call fiReadWord(string,word,.true.,ierr)
          call fiErrorMsg('keyword','CONDITION,TYPE', ierr)   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              index = pres_index
            case('TEMP','TEMPERATURE')
              index = temp_index
            case('CONC','CONCENTRATION')
              index = conc_index
            case('H','ENTHALPY')
              index = enthalpy_index
            case default
              call printErrMsg(option,'index not recognized in condition,gradient')
          end select
          call fiReadDouble(string,condition%gradient(index,X_DIRECTION),ierr)
          call fiErrorMsg('X Gradient','CONDITION', ierr)   
          call fiReadDouble(string,condition%gradient(index,Y_DIRECTION),ierr)
          call fiErrorMsg('Y Gradient','CONDITION', ierr)   
          call fiReadDouble(string,condition%gradient(index,Z_DIRECTION),ierr)
          call fiErrorMsg('Z Gradient','CONDITION', ierr)   
        enddo
      case('TEMPERATURE','TEMP')
        call ConditionReadValues(option,word,string,times,temperature,units(temp_index))
      case('ENTHALPY','H')
        call ConditionReadValues(option,word,string,times,enthalpy,units(enthalpy_index))
      case('PRESSURE','PRES','PRESS')
        call ConditionReadValues(option,word,string,times,pressure,units(pres_index))
      case('FLUX','VELOCITY','VEL')
        call ConditionReadValues(option,word,string,times,flux,units(pres_index))
      case('CONC','CONCENTRATION')
        call ConditionReadValues(option,word,string,times,concentration,units(conc_index))
    end select 
  
  enddo  
  
  ! check to ensure that class and type have been set
  if (len_trim(condition%class) < 1) then
    call printErrMsg(option,'"class" not set in condition')
  endif
  do idof=1,max_required_dof
    if (len_trim(ctype(idof)) < 1) then
      call printWrnMsg(option,'"type" not set in condition; set to dirichlet')
      ctype(idof) = 'dirichlet'
      itype(idof) = DIRICHLET_BC
    endif
  enddo
  ! check whether
  if (iphase == 0) then
    call printWrnMsg(option,'"iphase" not set in condition; set to 1')
    condition%iphase = 1
  endif
  
  allocate(condition%units(max_index))
  condition%units(1:max_index) = units(1:max_index)
  if (associated(times)) then
    max_size = size(times)
  else
    max_size = 1
  endif
  
  ! determine whether transient
  if (max_size > 1) condition%is_transient = .true.
  if (associated(times)) then
    if (times(1) > 1.d-40) condition%is_transient = .true.
  endif

  ! initialize time indices
  condition%cur_time_index = 1
  condition%max_time_index = max_size
  
  select case(option%imode)
    case(RICHARDS_MODE,MPH_MODE)
      if (.not.associated(enthalpy)) then
        max_index = max_index-1
        max_dof = max_dof-1
      endif
  end select
  condition%num_dof = max_dof
  
  allocate(condition%times(max_size),condition%values(max_dof,max_size))
  allocate(condition%cur_value(max_dof))
  allocate(condition%ctype(max_dof))
  allocate(condition%itype(max_dof))
  condition%times = -1.d20
  condition%values = -1.d20
  condition%cur_value = -1.d20
  condition%ctype = ""
  condition%itype = DIRICHLET_BC
  
  if (associated(times)) then
    condition%times(1:max_size) = times(1:max_size)
  else
    condition%times = 0.d0
  endif
  
  if (associated(pressure)) then
    if (size(pressure) < max_size) then
      condition%values(pres_index,1:max_size) = pressure(1)
    else
      condition%values(pres_index,1:max_size) = pressure(1:max_size)
    endif
  else if (associated(flux)) then
    if (size(flux) < max_size) then
      condition%values(pres_index,1:max_size) = flux(1)
    else
      condition%values(pres_index,1:max_size) = flux(1:max_size)
    endif
  else
    call printErrMsg(option,'Pressure conditon not set')
  endif
  
  if (associated(temperature)) then
    if (size(temperature) < max_size) then
      condition%values(temp_index,1:max_size) = temperature(1)
    else
      condition%values(temp_index,1:max_size) = temperature(1:max_size)
    endif
  else
    call printErrMsg(option,'Temperature conditon not set')
  endif
  
  condition%itype(1:max_dof) = itype(1:max_dof)
  condition%ctype(1:max_dof) = ctype(1:max_dof)
  
  if (associated(concentration)) then
    if (size(concentration) < max_size) then
      condition%values(conc_index,1:max_size) = concentration(1)
    else
      condition%values(conc_index,1:max_size) = concentration(1:max_size)
    endif
  else
    call printErrMsg(option,'Concentration conditon not set') 
  endif
  
  if (associated(enthalpy)) then
    if (size(enthalpy) < max_size) then
      condition%values(enthalpy_index,1:max_size) = enthalpy(1)
    else
      condition%values(enthalpy_index,1:max_size) = enthalpy(1:max_size)
    endif
  else
    call printWrnMsg(option,'Enthalpy conditon not set')   
  endif
  
  if (minval(condition%values) < -1.d19) then
    call printErrMsg(option,'Something is wrong in conditions....')
  endif
  
  condition%cur_value(1:max_dof) = condition%values(1:max_dof,1)
  
  if (associated(times)) deallocate(times)
  nullify(times)
  if (associated(pressure)) deallocate(pressure)
  nullify(pressure)
  if (associated(temperature)) deallocate(temperature)
  nullify(temperature)
  if (associated(flux)) deallocate(flux)
  nullify(flux)
  if (associated(concentration)) deallocate(concentration)
  nullify(concentration)
  
  if (associated(units)) deallocate(units)
  nullify(units)
  if (associated(itype)) deallocate(itype)
  nullify(itype)
  if (associated(ctype)) deallocate(ctype)
  nullify(ctype)

end subroutine ConditionRead

! ************************************************************************** !
!
! ConditionReadValues: Read the value(s) of a condition variable
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine ConditionReadValues(option,keyword,string,times,values,units)

  use Fileio_module
  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword
  real*8, pointer :: times(:), values(:)
  character(len=MAXWORDLENGTH) :: units
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  integer :: ierr
  
  ierr = 0
  call fiReadWord(string,word,.true.,ierr)
  call fiErrorMsg('file or value','CONDITION', ierr)
  call fiCharsToLower(word,len_trim(word))
  if (fiStringCompare(word,'file',4)) then
    call fiReadWord(string,word,.true.,ierr)
    error_string = keyword // ' FILE'
    call fiErrorMsg(error_string,'CONDITION', ierr)
    call ConditionReadValuesFromFile(word,times,values)
  else
    allocate(values(1))
    call fiReadDouble(word,values(1),ierr)
    call fiErrorMsg('value','CONDITION', ierr) 
  endif
  call fiReadWord(string,word,.true.,ierr)
  if (ierr /= 0) then
    word = trim(keyword) // ' UNITS'
    call fiDefaultMsg(word, ierr)
  else
    units = trim(word)
  endif

end subroutine ConditionReadValues

! ************************************************************************** !
!
! ConditionReadValuesFromFile: Read values from a external file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine ConditionReadValuesFromFile(filename,times,values)

  use Fileio_module
  use Utility_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: filename
  real*8, pointer :: times(:), values(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  real*8, pointer :: temp_times(:), temp_values(:)
  real*8 :: temp_time
  integer :: max_size = 1000
  integer :: fid
  integer :: count, i, status, ierr
  
  fid = 86
  open(unit=fid,file=filename,status="old",iostat=status)
  if (status /= 0) then
    print *, 'file: ', trim(filename), ' not found'
    stop
  endif
  
  allocate(temp_times(max_size))
  allocate(temp_values(max_size))
  
  temp_times = 0.d0
  temp_values = 0.d0
  
  count = 0
  do
    call fiReadFlotranString(fid,string,ierr)
    if (ierr /= 0) exit
    count = count + 1
    call fiReadDouble(string,temp_times(count),ierr)
    call fiErrorMsg('time','CONDITION FILE', ierr)   
    call fiReadDouble(string,temp_values(count),ierr)
    call fiErrorMsg('value','CONDITION FILE', ierr) 
    if (count+1 > max_size) then
      i = max_size
      call reallocateRealArray(temp_times,max_size) 
      ! careful.  max_size is already double, thus use i for now
      call reallocateRealArray(temp_values,i) 
    endif  
  enddo
  
  if (associated(times)) then
    if (count /= size(times)) then
      print *, 'Number of times (', count, ') in ', trim(filename), &
               ' does not match previous allocation: ', size(times)
      stop
    endif
    do i=1,count
      if (abs(times(i)-temp_times(i)) > 1.d-8) then
        print *, 'Time (', temp_times(i), ') in ', trim(filename), &
                 ' does not match previous allocation time: ', &
                 times(i), i
        stop
      endif
    enddo
  else
    allocate(times(count))
  endif

  if (associated(values)) then
    deallocate(values)
    allocate(values(count))
  else
    allocate(values(count))
  endif
  
  times(1:count) = temp_times(1:count)
  values(1:count) = temp_values(1:count)
  
  deallocate(temp_times)
  deallocate(temp_values)
  
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
  real*8 :: time
  
  type(condition_type), pointer :: condition
  integer :: cur_time_index, next_time_index, idof
  real*8 :: time_fraction
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
        ! potentially for initial condition
    if (time < 1.d-40 .or. condition%is_transient) then
    
      ! cycle times if at max_time_index and cyclic
      if (condition%cur_time_index == condition%max_time_index .and. &
          condition%is_cyclic .and. condition%max_time_index > 1) then
          
        do cur_time_index = 1, condition%max_time_index
          condition%times(cur_time_index) = condition%times(cur_time_index) + &     
                                   condition%times(condition%max_time_index)
        enddo
        condition%cur_time_index = 1
      endif
     
      cur_time_index = condition%cur_time_index
      next_time_index = min(condition%cur_time_index+1, &
                            condition%max_time_index)

      ! ensure that condition has started
      if (time >= condition%times(cur_time_index) .or. &
          abs(time-condition%times(cur_time_index)) < 1.d-40) then

        ! find appropriate time interval
        do
          if (time < condition%times(next_time_index) .or. &
              cur_time_index == next_time_index) &
            exit
          cur_time_index = next_time_index
          ! ensure that time index does not go beyond end of array
          if (next_time_index < condition%max_time_index) &
            next_time_index = next_time_index + 1
        enddo
        
        condition%cur_time_index = cur_time_index
        
        select case(condition%interpolation_method)
          case(STEP) ! just use the current value
            do idof=1,condition%num_dof
              condition%cur_value(idof) =  &
                                 condition%values(idof,condition%cur_time_index) 
            enddo 
          case(LINEAR) ! interpolate the value between times
            do idof=1,condition%num_dof
              ! interpolate value based on time
              condition%cur_time_index = cur_time_index
              if (cur_time_index < condition%max_time_index) then
                time_fraction = (time-condition%times(cur_time_index)) / &
                                  (condition%times(next_time_index) - &
                                   condition%times(cur_time_index))
                condition%cur_value(idof) = condition%values(idof,cur_time_index) + &
                                      time_fraction * &
                                      (condition%values(idof,next_time_index) - &
                                       condition%values(idof,cur_time_index))
              else
                condition%cur_value(idof) =  &
                                   condition%values(idof,condition%max_time_index) 
              endif
            enddo
        end select 
      endif 
    endif
    
    condition => condition%next
    
  enddo
  
end subroutine ConditionUpdate

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
  character(len=MAXNAMELENGTH) :: condition_name
  type(condition_list_type) :: condition_list

  type(condition_type), pointer :: condition
    
  nullify(ConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    if (len_trim(condition_name) == len_trim(condition%name) .and. &
        fiStringCompare(condition%name,condition_name, &
                        len_trim(condition_name))) then
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
  
  if (.not.associated(condition)) return
  
  if (associated(condition%times)) deallocate(condition%times)
  nullify(condition%times)
  if (associated(condition%values)) deallocate(condition%values)
  nullify(condition%values)
  if (associated(condition%cur_value)) deallocate(condition%cur_value)
  nullify(condition%cur_value)
  if (associated(condition%units)) deallocate(condition%units)
  nullify(condition%units)
  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  if (associated(condition%ctype)) deallocate(condition%ctype)
  nullify(condition%ctype)
  if (associated(condition%gradient)) deallocate(condition%gradient)
  nullify(condition%gradient)
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine ConditionDestroy
  
end module Condition_module
