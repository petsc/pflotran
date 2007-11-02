module Condition_module
 
  implicit none

  private
  
#include "definitions.h"

#define X_DIRECTION 1
#define Y_DIRECTION 2
#define Z_DIRECTION 3
 
  type, public :: condition_type
    integer :: id                                 ! id from which condition can be referenced
    integer :: itype                              ! integer describing type of condition
    character(len=MAXWORDLENGTH) :: class         ! character string describing class of condition
    character(len=MAXWORDLENGTH) :: ctype      ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: name          ! name of condition (e.g. initial, recharge)
    character(len=MAXWORDLENGTH), pointer :: units(:)      ! units
      ! units(1) = time
      ! units(2) = length
      ! units(3:ndof) = dofs in problem
    integer :: num_values                         ! number of entries in the arrays of values
    real*8, pointer :: times(:)                   ! array of times between which linear interpolation of values occurs
    real*8, pointer :: values(:,:)                ! array of condition values, size(ndof,max_time_index)
    real*8, pointer :: cur_value(:)               ! current value of condition a time t, size(ndof)
    real*8 :: datum(3)                            ! location of reference value(s) in domain
    real*8 :: gradient(3)                         ! rate at which reference value(s) change(s) over 3D space
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
  
  public :: createCondition, destroyCondition, readCondition, &
            addConditionToList, initConditionList, destroyConditionList, &
            getConditionPtrFromList
    
contains

! ************************************************************************** !
!
! createCondition: Creates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createCondition(option)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(condition_type), pointer :: createCondition
  
  type(condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%units)
  nullify(condition%times)
  nullify(condition%values)
  nullify(condition%cur_value)
  nullify(condition%next)
  condition%id = 0
  condition%num_values = 0
  condition%itype = 0
  condition%class = ""
  condition%ctype = ""
  condition%name = ""
  condition%datum = 0.d0
  condition%gradient = 0.d0
  condition%cur_time_index = 0
  condition%max_time_index = 0
  
  condition_count = condition_count + 1
  condition%id = condition_count
  
  allocate(condition%units(option%ndof+2))
  select case(option%imode)
    case(RICHARDS_MODE)
      condition%units(1) = 'yr'
      condition%units(2) = 'm'
      condition%units(3) = 'Pa'
      condition%units(4) = 'C'
      condition%units(5) = 'M'
  end select
  
  createCondition => condition

end function createCondition

! ************************************************************************** !
!
! readCondition: Reads a condition from the input file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine readCondition(condition,option,fid)

  use Option_module
  use Fileio_module
  
  implicit none
  
  type(condition_type) :: condition
  type(option_type) :: option
  integer :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  real*8, pointer :: pressure(:), flux(:), temperature(:), &
                     concentration(:), times(:)
  character(len=MAXWORDLENGTH), pointer :: units(:)
  integer :: time_dof, length_dof, pres_dof, temp_dof, conc_dof
  integer :: max_size, ierr

  nullify(times)
  nullify(pressure)
  nullify(temperature)
  nullify(flux)
  nullify(concentration)
  
  select case(option%imode)
    case(RICHARDS_MODE)
      allocate(units(option%ndof+2))
      pres_dof = 1
      temp_dof = 2
      conc_dof = 3
      time_dof = 4
      length_dof = 5
  end select
  units = ""

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg('keyword','CONDITION', ierr)   
      
    select case(trim(word))
    
      case('UNITS')
        do
          call fiReadWord(string,word,.true.,ierr)
          if (ierr /= 0) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%units(time_dof) = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%units(length_dof) = trim(word)
            case('Pa','KPa')
              condition%units(pres_dof) = trim(word)
            case('C','K')
              condition%units(temp_dof) = trim(word)
            case('M','mol/L')
              condition%units(conc_dof) = trim(word)
          end select
        enddo
      case('CLASS')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg('CLASS','CONDITION', ierr)   
        call fiCharsToLower(word,len_trim(word))
        condition%class = word
        if (fiStringCompare(word,'flow',4)) then
          condition%ctype = 'dirichlet'
        else
          condition%ctype = 'dirichlet'
        endif
        call fiReadWord(string,word,.true.,ierr)
        if (ierr /= 0) then
          call fiDefaultMsg('condition type',ierr)  
        else
          call fiCharsToLower(word,len_trim(word))
          condition%ctype = word
        endif
      case('TIME','TIMES')
        if (.not.associated(times)) then
          allocate(times(1))
        endif
        call fiReadDouble(string,times(1),ierr)
        call fiErrorMsg('TIME','CONDITION', ierr)   
      case('DATUM','DATM')
        call fiReadDouble(string,condition%datum(X_DIRECTION),ierr)
        call fiErrorMsg('X Datum','CONDITION', ierr)   
        call fiReadDouble(string,condition%datum(Y_DIRECTION),ierr)
        call fiErrorMsg('Y Datum','CONDITION', ierr)   
        call fiReadDouble(string,condition%datum(Z_DIRECTION),ierr)
        call fiErrorMsg('Z Datum','CONDITION', ierr)   
      case('GRADIENT','GRAD')
        call fiReadDouble(string,condition%gradient(X_DIRECTION),ierr)
        call fiErrorMsg('X Gradient','CONDITION', ierr)   
        call fiReadDouble(string,condition%gradient(Y_DIRECTION),ierr)
        call fiErrorMsg('Y Gradient','CONDITION', ierr)   
        call fiReadDouble(string,condition%gradient(Z_DIRECTION),ierr)
        call fiErrorMsg('Z Gradient','CONDITION', ierr)   
      case('TEMPERATURE','TEMP')
        call readValues(option,word,string,times,temperature,units(temp_dof))
      case('PRESSURE','PRES','PRESS')
        call readValues(option,word,string,times,pressure,units(pres_dof))
      case('FLUX','VELOCITY','VEL')
        call readValues(option,word,string,times,flux,units(pres_dof))
      case('CONC','CONCENTRATION')
        call readValues(option,word,string,times,concentration,units(conc_dof))
      case('END')
        exit
    end select 
  
  enddo  
  
  ! check to ensure that class and type have been set
  if (len_trim(condition%class) < 1) then
    call printErrMsg(option,'"class" not set in condition')
  endif
  if (len_trim(condition%ctype) < 1) then
    call printErrMsg(option,'"type" not set in condition')
  endif
  
  allocate(condition%units(2+option%ndof))
  condition%units(1:2+option%ndof) = units(1:2+option%ndof)
  if (associated(times)) then
    max_size = size(times)
  else
    max_size = 1
  endif
  allocate(condition%times(max_size),condition%values(option%ndof,max_size))
  condition%times = -999.d0
  condition%values = -999.d0
  
  if (associated(times)) then
    condition%times(1:max_size) = times(1:max_size)
  else
    condition%times = 0.d0
  endif
  
  if (associated(pressure)) then
    if (size(pressure) < max_size) then
      condition%values(pres_dof,1:max_size) = pressure(1)
    else
      condition%values(pres_dof,1:max_size) = pressure(1:max_size)
    endif
  else if (associated(flux)) then
    if (size(flux) < max_size) then
      condition%values(pres_dof,1:max_size) = flux(1)
    else
      condition%values(pres_dof,1:max_size) = flux(1:max_size)
    endif
  else
    print *, 'Error: Pressure conditon not set'
    stop
  endif
  
  if (associated(temperature)) then
    if (size(temperature) < max_size) then
      condition%values(temp_dof,1:max_size) = temperature(1)
    else
      condition%values(temp_dof,1:max_size) = temperature(1:max_size)
    endif
  else
    print *, 'Error: Temperature conditon not set'
    stop  
  endif
  
  if (associated(concentration)) then
    if (size(concentration) < max_size) then
      condition%values(conc_dof,1:max_size) = concentration(1)
    else
      condition%values(conc_dof,1:max_size) = concentration(1:max_size)
    endif
  else
    print *, 'Error: Concentration conditon not set'
    stop  
  endif
  
  if (minval(condition%values) < -998.d0) then
    print *, 'Error: Something is wrong in conditions....'
    stop  
  endif

end subroutine readCondition

! ************************************************************************** !
!
! readValues: Read the value(s) of a condition variable
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine readValues(option,keyword,string,times,values,units)

  use Fileio_module
  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword
  real*8, pointer :: times(:), values(:)
  character(len=MAXWORDLENGTH) :: units
  
  character(len=MAXWORDLENGTH) :: word
  integer :: ierr
  
  ierr = 0
  call fiReadWord(string,word,.true.,ierr)
  call fiErrorMsg('file or value','CONDITION', ierr)
  call fiCharsToLower(word,len_trim(word))
  if (fiStringCompare(word,'file',4)) then
    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg('PRESSURE FILE','CONDITION', ierr)
    call readValuesFromFile(word,times,values)
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

end subroutine readValues

! ************************************************************************** !
!
! readValuesFromFile: Read values from a external file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine readValuesFromFile(filename,times,values)

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
  open(unit=fid,file=filename,iostat=status)
  if (status /= 0) then
    print *, 'file: ', trim(filename), ' not found'.  
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

end subroutine

! ************************************************************************** !
!
! initConditionList: Initializes a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine initConditionList(list)

  implicit none

  type(condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine initConditionList

! ************************************************************************** !
!
! addConditionToList: Adds a new condition to a condition list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine addConditionToList(new_condition,list)

  implicit none
  
  type(condition_type), pointer :: new_condition
  type(condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine addConditionToList

! ************************************************************************** !
!
! getConditionPtrFromList: Returns a pointer to the condition matching &
!                          condition_name
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
function getConditionPtrFromList(condition_name,condition_list)

  use Fileio_module

  implicit none
  
  type(condition_type), pointer :: getConditionPtrFromList
  character(len=MAXNAMELENGTH) :: condition_name
  type(condition_list_type) :: condition_list

  type(condition_type), pointer :: condition
    
  nullify(getConditionPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    if (fiStringCompare(condition%name,condition_name, &
                        len_trim(condition_name))) then
      getConditionPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function getConditionPtrFromList

! ************************************************************************** !
!
! destroyConditionList: Deallocates a list of conditions
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine destroyConditionList(condition_list)

  implicit none
  
  type(condition_list_type), pointer :: condition_list
  
  type(condition_type), pointer :: condition, prev_condition
  
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call destroyCondition(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine destroyConditionList

! ************************************************************************** !
!
! destroyCondition: Deallocates a condition
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine destroyCondition(condition)

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
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine destroyCondition
  
end module Condition_module
