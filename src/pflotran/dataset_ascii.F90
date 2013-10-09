module Dataset_Ascii_class
 
  use Dataset_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public, extends(dataset_base_type) :: dataset_ascii_type
    PetscInt :: array_rank
  end type dataset_ascii_type
  
  interface DatasetAsciiRead
    module procedure DatasetAsciiOpenAndLoad
    module procedure DatasetAsciiLoad
  end interface
  
  public :: DatasetAsciiCreate, &
            DatasetAsciiInit, &
            DatasetAsciiVerify, &
            DatasetAsciiCast, &
            DatasetAsciiRead, &
            DatasetAsciiDestroy
  
contains

! ************************************************************************** !
!
! DatasetAsciiCreate: Creates ascii dataset class
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
function DatasetAsciiCreate()
  
  implicit none
  
  class(dataset_ascii_type), pointer :: dataset

  class(dataset_ascii_type), pointer :: DatasetAsciiCreate
  
  allocate(dataset)
  call DatasetAsciiInit(dataset)

  DatasetAsciiCreate => dataset
    
end function DatasetAsciiCreate

! ************************************************************************** !
!
! DatasetAsciiCast: Casts a dataset_base_type to database_ascii_type
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
function DatasetAsciiCast(this)

  use Dataset_Base_class
  
  implicit none

  class(dataset_base_type), pointer :: this

  class(dataset_ascii_type), pointer :: DatasetAsciiCast
  
  nullify(DatasetAsciiCast)
  select type (this)
    class is (dataset_ascii_type)
      DatasetAsciiCast => this
  end select
    
end function DatasetAsciiCast

! ************************************************************************** !
!
! DatasetAsciiInit: Initializes members of ascii dataset class
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
subroutine DatasetAsciiInit(this)
  
  implicit none
  
  class(dataset_ascii_type) :: this
  
  call DatasetBaseInit(this)
  this%array_rank = 0
    
end subroutine DatasetAsciiInit

! ************************************************************************** !
!
! DatasetAsciiOpenandLoad: Opens a file and calls the load routine.
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
subroutine DatasetAsciiOpenandLoad(this,filename,option)

  use Input_module
  use Option_module
  
  implicit none
  
  class(dataset_ascii_type) :: this
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  
  input => InputCreate(IUNIT_TEMP,filename,option)
  call DatasetAsciiLoad(this,input,option)
  call InputDestroy(input)

end subroutine DatasetAsciiOpenandLoad

! ************************************************************************** !
!
! DatasetAsciiLoad: Reads a text-based dataset from an ASCII file.
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
subroutine DatasetAsciiLoad(this,input,option)

  use Input_module
  use String_module
  use Utility_module, only : reallocateRealArray
  use Option_module
  use Units_module, only : UnitsConvertToInternal
  use Time_Storage_module  

  implicit none
  
  class(dataset_ascii_type) :: this
  type(input_type) :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: time_units
  character(len=MAXSTRINGLENGTH) :: string, data_units
  character(len=MAXWORDLENGTH) :: word
  PetscReal, pointer :: temp_array(:,:)
  PetscReal :: temp_time
  PetscReal :: conversion
  PetscInt :: max_size, offset
  PetscInt :: row_count, column_count, data_count, i
  PetscBool :: force_units_for_all_data
  PetscErrorCode :: ierr
  
  time_units = ''
  data_units = ''
  max_size = 1000

  row_count = 0
  ierr = 0
  do
    call InputReadFlotranString(input,option)
    ! reach the end of file or close out block
    if (InputError(input)) exit  ! check for end of file
    if (InputCheckExit(input,option)) exit  ! check for end of list
    ! check for units on first line
    if (row_count == 0) then
      string = input%buf
      ierr = 0
      call InputReadWord(string,word,PETSC_TRUE,ierr)
      call InputErrorMsg(input,option,'KEYWORD','CONDITION (LIST or FILE)')
      call StringToUpper(word)
      select case(word)
        case('TIME_UNITS')
          call InputReadWord(string,time_units,PETSC_TRUE,ierr)
          input%ierr = ierr
          call InputErrorMsg(input,option,'TIME_UNITS', &
                             'CONDITION (LIST or FILE)')
          call StringToLower(time_units) 
          cycle
        case('DATA_UNITS')
          ! it is possible to have more than one data unit. therefore, read the
          ! entire string
          data_units = adjustl(string)
          if (len_trim(data_units) < 1) then
            call InputErrorMsg(input,option,'DATA_UNITS', &
                               'CONDITION (LIST or FILE)')
          endif
          call StringToLower(data_units) 
          cycle
        case default
          ! copy the first row of actual data and count up the number of 
          ! columns.
          string = input%buf
          column_count = 0
          do
            ierr = 0
            call InputReadWord(string,word,PETSC_TRUE,ierr)
            if (ierr /= 0) exit   
            column_count = column_count + 1
          enddo
          ! allocate the 2d array to max_size rows and col_count columns.
          allocate(temp_array(column_count,max_size))
          temp_array = 0.d0
          ! do not cycle, as we now need to proceed.
      end select
    endif

    row_count = row_count + 1
    
    ! read columns of data, including the time in the first column
    do i = 1, column_count
      call InputReadDouble(input,option,temp_array(i,row_count))
      call InputErrorMsg(input,option,'column data','ascii dataset file') 
    enddo
    
    ! enlarge the array as needed.
    if (row_count+1 > max_size) then
      call reallocateRealArray(temp_array,max_size) 
    endif  
  enddo
  
  this%data_type = DATASET_REAL
  this%rank = 2
  allocate(this%dims(this%rank))
  data_count = column_count - 1 ! subtract 1 for time column
  this%dims(1) = data_count
  this%dims(2) = row_count
  this%time_storage => TimeStorageCreate()
  this%time_storage%max_time_index = row_count
  allocate(this%time_storage%times(row_count))
  this%time_storage%times = temp_array(1,1:row_count)
  allocate(this%rbuffer(data_count*row_count))
  this%rbuffer = 0.d0 ! we copy after units conversion for efficiency sake
  
  ! time units conversion
  if (len_trim(time_units) > 0) then
    conversion = UnitsConvertToInternal(time_units,option)
    this%time_storage%times(:) = conversion * &
                                 this%time_storage%times(:)
  endif
  ! data units conversion
  if (len_trim(data_units) > 0) then
    ! set flag to determine whether we check for data units for each
    ! data column.  if only one data unit is provided, it is applied
    ! to all columns by default.  otherwise, data units must be defined
    ! for each column - geh.
    force_units_for_all_data = PETSC_FALSE
    do i = 1, data_count ! number of data columns
      if (len_trim(data_units) > 0 .or. force_units_for_all_data) then
        ! the conditional immediately below will force 'conversion' to be
        ! calculated for each column. if a unit does not exist, the input
        ! error below will be spawned.
        if (i > 1) force_units_for_all_data = PETSC_TRUE 
        ierr = 0
        call InputReadWord(data_units,word,PETSC_TRUE,ierr)
        input%ierr = ierr
        call InputErrorMsg(input,option,'DATA_UNITS','CONDITION FILE')
        conversion = UnitsConvertToInternal(word,option)
      endif
      temp_array(i+1,:) = conversion * temp_array(i+1,:)
    enddo
  endif

  ! now that the data units conversion has taken place with temp_array, copy
  ! over to rbuffer.
  offset = 0
  do i = 1, row_count
    this%rbuffer(offset + 1:offset + data_count) = &
      temp_array(2:column_count,i)
    offset = offset + data_count
  enddo
  
  deallocate(temp_array)
  nullify(temp_array)

  if (this%array_rank > 0) then
    if (this%array_rank /= data_count) then
      write(word,*) this%array_rank
      option%io_buffer = 'Inconsistency between dataset prescribed rank (' // &
        trim(word) // ') and rank in file ('
      write(word,*) data_count
      option%io_buffer = trim(option%io_buffer) // ').'
      call printErrMsg(option)
    endif
  else
    this%array_rank = data_count
  endif
  
end subroutine DatasetAsciiLoad

! ************************************************************************** !
!
! DatasetAsciiVerify: Verifies that data structure is properly set up.
! author: Glenn Hammond
! date: 10/08/13
!
! ************************************************************************** !
subroutine DatasetAsciiVerify(this,option)
  
  use Option_module
  
  implicit none
  
  class(dataset_ascii_type) :: this
  type(option_type) :: option
  
  if (len_trim(this%name) < 1) then
    this%name = 'Unknown Ascii Dataset'
  endif
  call DatasetBaseVerify(this,option)
  if (associated(this%rbuffer)) then
    if (this%array_rank /= this%dims(1)) then
      option%io_buffer = '"array_rank" is not equal to "dims(1)" in dataset: ' &
                         // trim(this%name)
      call printErrMsg(option)
    endif
    ! set initial values
    this%rarray(:) = this%rbuffer(1:this%array_rank)
  endif
    
end subroutine DatasetAsciiVerify
  
! ************************************************************************** !
!
! DatasetAsciiStrip: Strips allocated objects within Ascii dataset object
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
subroutine DatasetAsciiStrip(this)

  implicit none
  
  class(dataset_ascii_type)  :: this
  
  call DatasetBaseStrip(this)
  
end subroutine DatasetAsciiStrip

! ************************************************************************** !
!
! DatasetAsciiDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 10/03/13
!
! ************************************************************************** !
subroutine DatasetAsciiDestroy(this)

  implicit none
  
  class(dataset_ascii_type), pointer :: this
  
  if (.not.associated(this)) return
  
  call DatasetAsciiStrip(this)
  
  deallocate(this)
  nullify(this)
  
end subroutine DatasetAsciiDestroy

end module Dataset_Ascii_class
