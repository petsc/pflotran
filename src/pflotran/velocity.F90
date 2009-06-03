module Velocity_module
 
  implicit none

  private
  
#include "definitions.h"

  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2
  
  type, public :: velocity_dataset_type
    PetscInt :: rank
    PetscTruth :: is_transient
    PetscTruth :: is_cyclic
    PetscInt :: interpolation_method
    PetscReal, pointer :: times(:)
    PetscReal, pointer :: values(:,:)
    PetscReal, pointer :: cur_value(:)
    PetscInt :: cur_time_index
    PetscInt :: max_time_index
    PetscReal :: time_shift
  end type velocity_dataset_type
  
      
  public :: VelocityDatasetCreate, VelocityDatasetDestroy, VelocityDatasetRead, &
            VelocityDatasetUpdate, VelocityDatasetVerify
    
contains

! ************************************************************************** !
!
! VelocityDatasetCreate: Creates a velocity data set
! author: Glenn Hammond
! date: 06/02/09
!
! ************************************************************************** !
function VelocityDatasetCreate()

  implicit none
  
  type(velocity_dataset_type), pointer :: VelocityDatasetCreate
  
  type(velocity_dataset_type), pointer :: dataset

  allocate(dataset)
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
   
  VelocityDatasetCreate => dataset
    
end function VelocityDatasetCreate

! ************************************************************************** !
!
! VelocityDatasetRead: Reads a velocity data set from the input file
! author: Glenn Hammond
! date: 06/02/09
!
! ************************************************************************** !
subroutine VelocityDatasetRead(dataset,input,option)

  use Option_module
  use Input_module
  use String_module
  use Logging_module
  use Units_module 
  
  implicit none
  
  type(velocity_dataset_type) :: dataset
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, units
  PetscReal :: units_conversion

  PetscErrorCode :: ierr

!  call PetscLogEventBegin(logging%event_flow_condition_read, &
!                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
!                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

  dataset%rank = 3
  dataset%interpolation_method = STEP
  dataset%is_cyclic = PETSC_FALSE
  
  units = ''
  
  ! read the velocity data set
  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,'VELOCITY_DATASET')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','VELOCITY_DATASET')   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        call InputReadWord(input,option,units,PETSC_TRUE)
        call InputErrorMsg(input,option,'UNITS','VELOCITY_DATASET')       
      case('CYCLIC')
        dataset%is_cyclic = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','VELOCITY_DATASET')   
        call StringToLower(word)
        select case(word)
          case('step')
            dataset%interpolation_method = STEP
          case('linear') 
            dataset%interpolation_method = LINEAR
        end select
      case('VELOCITY')
        call VelocityDatasetReadValues(input,option,word,string, &
                                       dataset,units)
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized in VELOCITY_DATASET'
        call printErrMsg(option)                                 
    end select 
  
  enddo

  if (len_trim(units) > 1) then
    units_conversion = UnitsConvertToInternal(units,option)
    dataset%values = dataset%values * units_conversion
    word = units(index(units,'/')+1:)
    units_conversion = UnitsConvertToInternal(word,option)
    dataset%times = dataset%times * units_conversion
  endif
  
  call VelocityDatasetVerify(option, dataset)
  
!  call PetscLogEventEnd(logging%event_flow_condition_read, &
!                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
!                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine VelocityDatasetRead

! ************************************************************************** !
!
! VelocityDatasetReadValues: Read the value(s) of a velocity data set
! author: Glenn Hammond
! date: 06/02/09
!
! ************************************************************************** !
subroutine VelocityDatasetReadValues(input,option,keyword,string,dataset,units)

  use Input_module
  use String_module
  use Option_module
  use Logging_module

  implicit none
  
  type(input_type), target :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword
  type(velocity_dataset_type) :: dataset
  character(len=MAXWORDLENGTH) :: units
  
  type(input_type), pointer :: input2
  character(len=MAXSTRINGLENGTH) :: string2, filename
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscTruth :: read_from_file
  PetscTruth :: read_multiple_values
  PetscInt :: irank
  PetscErrorCode :: ierr

!  call PetscLogEventBegin(logging%event_flow_condition_read_values, &
!                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
!                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

  read_from_file = PETSC_FALSE
  read_multiple_values = PETSC_TRUE
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  if (input%ierr == 0) then
    call StringToLower(word)
    if (StringCompare(word,'file',FOUR_INTEGER)) then
      call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_TRUE)
      input%err_buf = trim(keyword) // ' FILE'
      input%err_buf2 = 'VELOCITY_DATASET'
      call InputErrorMsg(input,option)
      read_from_file = PETSC_TRUE
    else
      read_multiple_values = PETSC_FALSE
    endif
  endif
  
  if (read_from_file .or. read_multiple_values) then
    if (read_from_file) then
      input2 => InputCreate(IUNIT_TEMP,filename)
    else
      input2 => input
    endif
    call VelocityDatasetReadFromFile(input2,dataset,option)
    if (read_from_file) call InputDestroy(input2)
  else
    input%buf = trim(string2)
    allocate(dataset%times(1))
    dataset%times = 0.d0
    allocate(dataset%values(dataset%rank,1))
    dataset%values = 0.d0
    do irank=1,dataset%rank
      call InputReadDouble(input,option,dataset%values(irank,1))
      write(input%err_buf,'(a,i2)') trim(keyword) // ' dataset_values, irank = ', irank
      input%err_buf2 = 'VELCITY_DATASET'
      call InputErrorMsg(input,option) 
    enddo
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) then
      word = trim(keyword) // ' UNITS'
      call InputDefaultMsg(input,option,word)
    else
      units = trim(word)
    endif
  endif

!  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
!                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
!                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)    

end subroutine VelocityDatasetReadValues

! ************************************************************************** !
!
! VelocityDatasetReadFromFile: Read values from a external file
! author: Glenn Hammond
! date: 10/31/07
!
! ************************************************************************** !
subroutine VelocityDatasetReadFromFile(input,dataset,option)

  use Input_module
  use String_module
  use Utility_module
  use Option_module

  implicit none
  
  type(option_type) :: option
  type(velocity_dataset_type) :: dataset
  character(len=MAXSTRINGLENGTH) :: filename
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: temp_times(:), temp_array1(:), temp_array2(:), &
                        temp_array3(:)
  PetscReal :: temp_time
  PetscInt :: max_size
  PetscInt :: temp_max_size
  PetscInt :: count, i, status
  type(input_type), pointer :: input
  
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
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    count = count + 1
    call InputReadDouble(input,option,temp_times(count))
    call InputErrorMsg(input,option,'time','VELOCITY_DATASET FILE')   
    call InputReadDouble(input,option,temp_array1(count))
    call InputErrorMsg(input,option,'array1','VELOCITY_DATASET FILE')
    if (dataset%rank > 1) then
      call InputReadDouble(input,option,temp_array2(count))
      call InputErrorMsg(input,option,'array2','VELOCITY_DATASET FILE') 
    endif
    if (dataset%rank > 2) then
      call InputReadDouble(input,option,temp_array3(count))
      call InputErrorMsg(input,option,'array3','VELOCITY_DATASET FILE') 
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
  
  if (associated(dataset%times)) deallocate(dataset%times)
  nullify(dataset%times)
  if (associated(dataset%values)) deallocate(dataset%values)
  nullify(dataset%values)

  allocate(dataset%times(count))
  allocate(dataset%values(dataset%rank,count))

  dataset%times(1:count) = temp_times(1:count)
  dataset%values(1,1:count) = temp_array1(1:count)
  if (dataset%rank > 1) dataset%values(2,1:count) = temp_array2(1:count)
  if (dataset%rank > 2) dataset%values(3,1:count) = temp_array3(1:count)
  
  deallocate(temp_times)
  deallocate(temp_array1)
  if (dataset%rank > 1) deallocate(temp_array2)
  if (dataset%rank > 2) deallocate(temp_array3)
  
end subroutine VelocityDatasetReadFromFile

! ************************************************************************** !
!
! VelocityDatasetVerify: Verifies the data in a dataset
! author: Glenn Hammond
! date: 02/04/08
!
! ************************************************************************** !
subroutine VelocityDatasetVerify(option, dataset)
  use Option_module

  implicit none
  
  type(option_type) :: option
  type(velocity_dataset_type) :: dataset
  
  PetscInt :: array_size
  character(len=MAXWORDLENGTH) :: size1, size2
  
  dataset%max_time_index = 1
  if (.not.associated(dataset%times)) then
    array_size = 1
    allocate(dataset%times(array_size))
    dataset%times = 0.d0
  endif
  if (.not.associated(dataset%values)) then
    ! allocate to size 1
    array_size = 1
    allocate(dataset%values(1:dataset%rank,array_size))
    dataset%values(1:dataset%rank,1:array_size) = 0.d0
  endif
  dataset%max_time_index = size(dataset%times,1) 
  if (dataset%max_time_index > 1) then
    if (size(dataset%values,2) /= & ! size of 2nd rank
        dataset%max_time_index) then
      write(size1,*) size(dataset%times,1)
      write(size2,*) size(dataset%values,2)
      option%io_buffer = 'times/values ('//trim(size1)//'/'//trim(size2) // &
                         ') array size mismatch in velocity dataaset'
      call printErrMsg(option) 
    endif
    dataset%is_transient = PETSC_TRUE
  else
    dataset%is_transient = PETSC_FALSE
  endif  
  dataset%cur_time_index = 1
  
  if (associated(dataset%cur_value)) deallocate(dataset%cur_value)
  allocate(dataset%cur_value(dataset%rank))
  dataset%cur_value(1:dataset%rank) = dataset%values(1:dataset%rank,1)

  dataset%time_shift = dataset%times(dataset%max_time_index)

end subroutine VelocityDatasetVerify

! ************************************************************************** !
!
! VelocityDatasetUpdate: Updates a velocity dataset
! author: Glenn Hammond
! date: 06/02/09
!
! ************************************************************************** !
subroutine VelocityDatasetUpdate(option,time,dataset)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: time
  PetscTruth :: is_cyclic
  PetscInt :: interpolation_method
  type(velocity_dataset_type) :: dataset
  
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
  
end subroutine VelocityDatasetUpdate

! ************************************************************************** !
!
! VelocityDatasetDestroy: Destroys a velocity dataset
! author: Glenn Hammond
! date: 06/02/09
!
! ************************************************************************** !
subroutine VelocityDatasetDestroy(dataset)

  implicit none
  
  type(velocity_dataset_type) :: dataset
  
  if (associated(dataset%times)) deallocate(dataset%times)
  nullify(dataset%times)
  if (associated(dataset%values)) deallocate(dataset%values)
  nullify(dataset%values)  
  if (associated(dataset%cur_value)) deallocate(dataset%cur_value)
  nullify(dataset%cur_value)  

end subroutine VelocityDatasetDestroy

end module Velocity_module
