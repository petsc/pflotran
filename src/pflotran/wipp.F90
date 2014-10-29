module Creep_Closure_module
  
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none
  
  private

#include "finclude/petscsys.h"

  type, public :: creep_closure_type
    character(len=MAXWORDLENGTH) :: material_name
    PetscInt :: imat
    PetscInt :: num_times
    PetscInt :: num_values_per_time
    class(lookup_table_general_type), pointer :: lookup_table
  contains
    procedure, public :: Read => CreepClosureRead
    procedure, public :: Evaluate => CreepClosureEvaluate
    procedure, public :: Test => CreepClosureTest
  end type creep_closure_type
  
  class(creep_closure_type), pointer, public :: creep_closure
  
  interface CreepClosureDestroy
    module procedure CreepClosureDestroy1
    module procedure CreepClosureDestroy2
  end interface
  
  public :: CreepClosureInit, &
            CreepClosureCreate, &
            CreepClosureDestroy
  
contains

! ************************************************************************** !

subroutine CreepClosureInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: CreepClosureCreate
  
  if (associated(creep_closure)) then
    call CreepClosureDestroy(creep_closure)
  endif
  nullify(creep_closure)
  
end subroutine CreepClosureInit

! ************************************************************************** !

function CreepClosureCreate()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: CreepClosureCreate
  
  allocate(CreepClosureCreate)
  CreepClosureCreate%material_name = ''
  CreepClosureCreate%imat = UNINITIALIZED_INTEGER
  CreepClosureCreate%num_times = UNINITIALIZED_INTEGER
  CreepClosureCreate%num_values_per_time = UNINITIALIZED_INTEGER
  nullify(CreepClosureCreate%lookup_table)
  
end function CreepClosureCreate

! ************************************************************************** !

subroutine CreepClosureRead(this,input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  
  implicit none
  
  class(creep_closure_type) :: this
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string = 'CREEP_CLOSURE'
  type(input_type), pointer :: input2
  PetscInt :: temp_int
  PetscReal :: time_units_conversion

  time_units_conversion = 1.d0
  filename = ''
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('FILENAME') 
        call InputReadWord(input,option,filename,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename',error_string)
      case('MATERIAL') 
        call InputReadWord(input,option,this%material_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'material',error_string)
     case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in creep closure'    
        call printErrMsg(option)
    end select
  enddo
  
  if (len_trim(filename) < 1) then
    option%io_buffer = 'FILENAME must be specified for CREEP_CLOSURE.'
    call printErrMsg(option)
  endif
  
  this%lookup_table => LookupTableCreateGeneral(2)
  error_string = 'CREEP_CLOSURE file'
  input2 => InputCreate(IUNIT_TEMP,filename,option)
  input2%ierr = 0
  do
  
    call InputReadPflotranString(input2,option)

    if (InputError(input2)) exit

    call InputReadWord(input2,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input2,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('NUM_TIMES') 
        call InputReadInt(input2,option,this%num_times)
        call InputErrorMsg(input2,option,'number of times',error_string)
      case('NUM_VALUES_PER_TIME') 
        call InputReadInt(input2,option,this%num_values_per_time)
        call InputErrorMsg(input2,option,'number of pressure',error_string)
      case('TIME_UNITS') 
        call InputReadWord(input2,option,word,PETSC_TRUE) 
        call InputErrorMsg(input2,option,'UNITS','CONDITION')   
        call StringToLower(word)
        time_units_conversion = UnitsConvertToInternal(word,option)
      case('TIME')
        if (Uninitialized(this%num_times) .or. &
            Uninitialized(this%num_values_per_time)) then
          option%io_buffer = 'NUM_TIMES and NUM_VALUES_PER_TIME must be ' // &
            'specified prior to reading the corresponding arrays.'
          call printErrMsg(option)
        endif
        this%lookup_table%dims(1) = this%num_times
        this%lookup_table%dims(2) = this%num_values_per_time
        temp_int = this%num_times*this%num_values_per_time
        allocate(this%lookup_table%axis1%values(this%num_times))
        allocate(this%lookup_table%axis2%values(temp_int))
        allocate(this%lookup_table%data(temp_int))
        string = 'TIME in CREEP_CLOSURE'
        call UtilityReadRealArray(this%lookup_table%axis1%values, &
                                  -1,string, &
                                  input2,option)
        this%lookup_table%axis1%values = this%lookup_table%axis1%values * &
          time_units_conversion
      case('PRESSURE') 
        string = 'PRESSURE in CREEP_CLOSURE'
        call UtilityReadRealArray(this%lookup_table%axis2%values, &
                                  -1, &
                                  string,input2,option)
      case('POROSITY') 
        string = 'POROSITY in CREEP_CLOSURE'
        call UtilityReadRealArray(this%lookup_table%data, &
                                  -1, &
                                  string,input2,option)
     case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in creep closure file.'    
        call printErrMsg(option)
    end select
  enddo
  call InputDestroy(input2)
  
  if (size(this%lookup_table%axis1%values) /= this%num_times) then
    option%io_buffer = 'Number of times does not match NUM_TIMES.'
    call printErrMsg(option)
  endif  
  if (size(this%lookup_table%axis2%values) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of pressures does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif
  if (size(this%lookup_table%data) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of porosities does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif
  
end subroutine CreepClosureRead

! ************************************************************************** !

function CreepClosureEvaluate(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  PetscReal :: CreepClosureEvaluate
  
  CreepClosureEvaluate = this%lookup_table%Sample(time,pressure)
  
end function CreepClosureEvaluate


! ************************************************************************** !

subroutine CreepClosureTest(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  PetscReal :: CreepClosureEvaluate
  
  print *, time, pressure, this%Evaluate(time,pressure)
  
end subroutine CreepClosuretest

! ************************************************************************** !

subroutine CreepClosureDestroy1()
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  call CreepClosureDestroy(creep_closure)

end subroutine CreepClosureDestroy1

! ************************************************************************** !

subroutine CreepClosureDestroy2(creep_closure)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: creep_closure
  
  if (.not.associated(creep_closure)) return

  call LookupTableDestroy(creep_closure%lookup_table)
  deallocate(creep_closure)
  nullify(creep_closure)

end subroutine CreepClosureDestroy2

end module Creep_Closure_module
