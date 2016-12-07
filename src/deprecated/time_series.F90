module Time_Series_module
 
#include "finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: time_series_type
    character(len=MAXWORDLENGTH) :: name
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
    PetscReal :: lame_auxvariable_remove_me
  end type time_series_type
  
  public :: TimeSeriesCreate, &
            TimeSeriesGetTimes, &
            TimeSeriesVerify, &
            TimeSeriesUpdate, &
            TimeSeriesPrint, &
            TimeSeriesDestroy

contains

! ************************************************************************** !

function TimeSeriesCreate()
  ! 
  ! Initializes a time series
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  implicit none
  
  type(time_series_type), pointer :: time_series
  type(time_series_type), pointer :: TimeSeriesCreate

  allocate(time_series)
  time_series%name = ''
  nullify(time_series%times)
  nullify(time_series%values)
  nullify(time_series%cur_value)
  time_series%cur_time_index = 0
  time_series%max_time_index = 0
  time_series%rank = 0
  time_series%is_cyclic = PETSC_FALSE
  time_series%is_transient = PETSC_FALSE
  time_series%interpolation_method = INTERPOLATION_NULL
  time_series%time_shift = 0.d0
  time_series%lame_auxvariable_remove_me = 0.d0
  
  TimeSeriesCreate => time_series
    
end function TimeSeriesCreate

! ************************************************************************** !

subroutine TimeSeriesVerify(option, default_time, time_series, &
                            default_time_series)
  ! 
  ! Verifies the data in a time_series
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 
  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: condition_name
  character(len=MAXWORDLENGTH) :: sub_condition_name
  character(len=MAXWORDLENGTH) :: size1, size2
  PetscReal :: default_time
  type(time_series_type) :: time_series
  type(time_series_type) :: default_time_series
  
  PetscInt :: array_size
  
  if (default_time_series%is_cyclic) time_series%is_cyclic = PETSC_TRUE
  if (time_series%interpolation_method == INTERPOLATION_NULL) &
    time_series%interpolation_method = default_time_series%interpolation_method
  
  time_series%max_time_index = 1
  if (.not.associated(time_series%times)) then
    if (.not.associated(default_time_series%times)) then
      array_size = 1
      allocate(time_series%times(array_size))
      time_series%times = default_time
    else
      array_size = size(default_time_series%times,1)
      allocate(time_series%times(array_size))
      time_series%times(1:array_size) = default_time_series%times(1:array_size)
    endif
  endif
  if (.not.associated(time_series%values)) then
    if (.not.associated(default_time_series%values)) then
      ! allocate to size 1
      array_size = 1
      allocate(time_series%values(1:time_series%rank,array_size))
      time_series%values(1:time_series%rank,1:array_size) = 0.d0
    else
      ! copy default
      array_size = size(default_time_series%values,2)
      allocate(time_series%values(1:time_series%rank,array_size))
      time_series%values(1:time_series%rank,1:array_size) = &
        default_time_series%values(1:time_series%rank,1:array_size)
    endif
  endif
  time_series%max_time_index = size(time_series%times,1) 
  if (time_series%max_time_index > 1) then
    if (size(time_series%values,2) /= & ! size of 2nd rank
        time_series%max_time_index) then
      write(size1,*) size(time_series%times,1)
      write(size2,*) size(time_series%values,2)
      option%io_buffer = 'times/values ('//trim(size1)//'/'//trim(size2) // &
                         ') array size mismatch in time series "' // &
                         trim(time_series%name) // '".'
      call printErrMsg(option) 
    endif
    time_series%is_transient = PETSC_TRUE
  else
    time_series%is_transient = PETSC_FALSE
  endif  
  time_series%cur_time_index = 1
  
  allocate(time_series%cur_value(time_series%rank))
  time_series%cur_value(1:time_series%rank) = &
    time_series%values(1:time_series%rank,1)

  time_series%time_shift = time_series%times(time_series%max_time_index)
! print *,'condition: ',time_series%max_time_index,time_series%time_shift

end subroutine TimeSeriesVerify

! ************************************************************************** !

subroutine TimeSeriesGetTimes(option, time_series, max_sim_time, times)
  ! 
  ! Fills an array of times based on time_series
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(time_series_type) :: time_series
  PetscReal :: max_sim_time
  PetscReal, pointer :: times(:)
  
  PetscInt :: num_times
  PetscInt :: itime
  PetscReal :: time_shift
  PetscReal, allocatable :: temp_times(:)

  if (.not.time_series%is_cyclic .or. time_series%max_time_index == 1) then
    allocate(times(time_series%max_time_index))
    times =  time_series%times
  else ! cyclic
    num_times = (int(max_sim_time / &
                     time_series%times(time_series%max_time_index))+1)* &
                time_series%max_time_index
    allocate(temp_times(num_times))
    temp_times = 0.d0

    num_times = 0
    itime = 0
    time_shift = 0.d0
    do
      num_times = num_times + 1
      itime = itime + 1
      ! exist for non-cyclic
      if (itime > time_series%max_time_index) exit
      temp_times(num_times) = time_series%times(itime) + time_shift
      if (mod(itime,time_series%max_time_index) == 0) then
        itime = 0
        time_shift = time_shift + time_series%times(time_series%max_time_index) 
      endif 
      ! exit for cyclic
      if (temp_times(num_times) >= max_sim_time) exit
    enddo

    allocate(times(num_times))
    times(:) = temp_times(1:num_times)
    deallocate(temp_times)
  endif
 
end subroutine TimeSeriesGetTimes

! ************************************************************************** !

subroutine TimeSeriesPrint(time_series,option)
  ! 
  ! Prints flow condition time_series info
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module

  implicit none
  
  type(time_series_type) :: time_series
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string

  if (time_series%is_cyclic) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(8x,''Is cyclic: '',a)') trim(string)

  100 format(8x,'Is transient: ', a)
  if (time_series%is_transient) then
    write(option%fid_out,100) 'yes'
    write(option%fid_out,'(8x,''  Number of values: '', i7)') &
      time_series%max_time_index
    select case(time_series%interpolation_method)
      case(INTERPOLATION_STEP)
        string = 'step'
      case(INTERPOLATION_LINEAR)
        string = 'linear'
    end select
    write(option%fid_out,'(8x,''  Interpolation method: '', a)') &
      trim(string)
    write(option%fid_out,'(8x,''Start value:'',es16.8)') &
      time_series%cur_value(1)
  else
    write(option%fid_out,100) 'no'
    write(option%fid_out,'(8x,''Value:'',es16.8)') time_series%cur_value(1)
  endif

            
end subroutine TimeSeriesPrint

! ************************************************************************** !

subroutine TimeSeriesUpdate(option,time,time_series)
  ! 
  ! Updates a transient condition time_series
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: time
  PetscBool :: is_cyclic
  PetscInt :: interpolation_method
  type(time_series_type) :: time_series
  
  PetscInt :: irank
  PetscInt :: cur_time_index
  PetscInt :: next_time_index
  PetscReal :: time_fraction

  ! cycle times if at max_time_index and cyclic
  if (time_series%cur_time_index == time_series%max_time_index .and. &
      time_series%is_cyclic .and. time_series%max_time_index > 1) then
      
    do cur_time_index = 1, time_series%max_time_index
      time_series%times(cur_time_index) = time_series%times(cur_time_index) + &   
                               time_series%time_shift
    enddo
    time_series%cur_time_index = 1
  endif
 
  cur_time_index = time_series%cur_time_index
  next_time_index = min(time_series%cur_time_index+1, &
                        time_series%max_time_index)

  ! ensure that condition has started
  if (time >= time_series%times(cur_time_index) .or. &
      dabs(time-time_series%times(cur_time_index)) < 1.d-40) then

    ! find appropriate time interval
    do
      if (time < time_series%times(next_time_index) .or. &
          cur_time_index == next_time_index) &
        exit
      cur_time_index = next_time_index
      ! ensure that time index does not go beyond end of array
      if (next_time_index < time_series%max_time_index) &
        next_time_index = next_time_index + 1
    enddo
    
    time_series%cur_time_index = cur_time_index
    
    select case(time_series%interpolation_method)
      case(INTERPOLATION_STEP) ! just use the current value
        do irank=1,time_series%rank
          time_series%cur_value(irank) =  &
                             time_series%values(irank,time_series%cur_time_index) 
        enddo 
      case(INTERPOLATION_LINEAR) ! interpolate the value between times
        do irank=1,time_series%rank
          ! interpolate value based on time
          time_series%cur_time_index = cur_time_index
          if (cur_time_index < time_series%max_time_index) then
            time_fraction = (time-time_series%times(cur_time_index)) / &
                              (time_series%times(next_time_index) - &
                               time_series%times(cur_time_index))
            time_series%cur_value(irank) = &
              time_series%values(irank,cur_time_index) + &
              time_fraction * &
                (time_series%values(irank,next_time_index) - &
                 time_series%values(irank,cur_time_index))
          else
            time_series%cur_value(irank) =  &
              time_series%values(irank,time_series%max_time_index) 
          endif
        enddo
    end select 
  endif    
  
end subroutine TimeSeriesUpdate

! ************************************************************************** !

subroutine TimeSeriesDestroy(time_series)
  ! 
  ! Destroys a time_series associated with a sub_condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  ! 

  implicit none
  
  type(time_series_type), pointer :: time_series
  
  if (.not.associated(time_series)) return
  
  if (associated(time_series%times)) deallocate(time_series%times)
  nullify(time_series%times)
  if (associated(time_series%values)) deallocate(time_series%values)
  nullify(time_series%values)  
  if (associated(time_series%cur_value)) deallocate(time_series%cur_value)
  nullify(time_series%cur_value)  
  
  deallocate(time_series)
  nullify(time_series)

end subroutine TimeSeriesDestroy

end module Time_Series_module
