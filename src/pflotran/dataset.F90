module Dataset_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: dataset_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: itype
    PetscBool :: realization_dependent
    PetscReal :: scalar
    PetscReal :: vector(3)
    PetscReal :: tensor(3:3)
    type(dataset_type), pointer :: next
  end type dataset_type

  type, public :: dataset_list_type
    PetscInt :: num_datasets
    type(dataset_type), pointer :: first
    type(dataset_type), pointer :: last
    type(dataset_type), pointer :: array(:)    
  end type dataset_list_type
    
  public :: DatasetCreate, &
            DatasetRead, &
            DatasetAddToList, &
            DatasetGetPointer, &
            DatasetDestroy

contains

! ************************************************************************** !
!
! DatasetCreate: Creates a dataset object
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
function DatasetCreate()
  
  implicit none

  type(dataset_type), pointer :: DatasetCreate
  
  type(dataset_type), pointer :: dataset
  
  allocate(dataset)
  dataset%name = ''
  dataset%filename = ''
  dataset%itype = DATASET_HETEROGENEOUS
  dataset%realization_dependent = PETSC_FALSE
  dataset%scalar = 0.d0
  dataset%vector = 0.d0
  dataset%tensor = 0.d0
  nullify(dataset%next)
  DatasetCreate => dataset

end function DatasetCreate

! ************************************************************************** !
!
! DatasetRead: Reads in contents of a dataset card
! author: Glenn Hammond
! date: 01/12/11
! 
! ************************************************************************** !
subroutine DatasetRead(dataset,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(dataset_type) :: dataset
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DATASET')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('NAME') 
        call InputReadWord(input,option,dataset%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','DATASET')
      case('FILENAME') 
        call InputReadWord(input,option,dataset%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','DATASET')
      case('TYPE') 
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'type','DATASET')
        call StringToUpper(word)
        select case (trim(word))
          case('HETEROGENEOUS')
            dataset%itype = DATASET_HETEROGENEOUS
          case default
            option%io_buffer = 'Dataset type: ' // trim(word) // &
                               ' not recognized in dataset'    
            call printErrMsg(option)
        end select
      case('REALIZATION_DEPENDENT')
        dataset%realization_dependent = PETSC_TRUE
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in dataset'    
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine DatasetRead

! ************************************************************************** !
!
! DatasetAddToList: Adds a dataset to linked list
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine DatasetAddToList(dataset,list)

  implicit none
  
  type(dataset_type), pointer :: dataset
  type(dataset_type), pointer :: list

  type(dataset_type), pointer :: cur_dataset
  
  if (associated(list)) then
    cur_dataset => list
    ! loop to end of list
    do
      if (.not.associated(cur_dataset%next)) exit
      cur_dataset => cur_dataset%next
    enddo
    cur_dataset%next => dataset
  else
    list => dataset
  endif
  
end subroutine DatasetAddToList

! ************************************************************************** !
!
! DatasetGetPointer: Returns the pointer to the dataset named "dataset_name"
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
function DatasetGetPointer(dataset_list, dataset_name, debug_string, option)

  use Option_module
  use String_module
  
  type(dataset_type), pointer :: dataset_list
  character(len=MAXWORDLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: debug_string
  type(option_type) :: option

  type(dataset_type), pointer :: DatasetGetPointer
  PetscBool :: found
  type(dataset_type), pointer :: cur_dataset

  found = PETSC_FALSE
  cur_dataset => dataset_list
  do 
    if (.not.associated(cur_dataset)) exit
    if (StringCompare(dataset_name, &
                      cur_dataset%name,MAXWORDLENGTH)) then
      found = PETSC_TRUE
      DatasetGetPointer => cur_dataset
      return
    endif
    cur_dataset => cur_dataset%next
  enddo
  if (.not.found) then
    option%io_buffer = 'Dataset "' // trim(dataset_name) // '" in "' // &
             trim(debug_string) // '" not found among available dataset.'
    call printErrMsg(option)    
  endif

end function DatasetGetPointer

! ************************************************************************** !
!
! DatasetDestroy: Destroys a dataset
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
recursive subroutine DatasetDestroy(dataset)

  implicit none
  
  type(dataset_type), pointer :: dataset
  
  if (.not.associated(dataset)) return
  
  deallocate(dataset)
  nullify(dataset)
  
end subroutine DatasetDestroy

end module Dataset_module
