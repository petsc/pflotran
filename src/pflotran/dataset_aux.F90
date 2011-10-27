module Dataset_Aux_module
 
  implicit none

  private

#include "definitions.h"
 
  type, public :: dataset_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: h5_dataset_name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: data_type
    PetscInt :: data_dim ! dimensions of data: XY, X, YXZ, etc.
    PetscBool :: realization_dependent
    ! all data stored internally as a 1D array
    PetscReal, pointer :: time_array(:)
    PetscInt, pointer :: iarray(:)
    PetscReal, pointer :: rarray(:)
    PetscReal :: rmax ! maximum rarray value in dataset
    PetscReal :: rmin ! maximum rarray value in dataset
    PetscInt :: ndims
    PetscInt, pointer :: dims(:)
    PetscReal, pointer :: origin(:)
    PetscReal, pointer :: discretization(:)
    type(dataset_type), pointer :: next
  end type dataset_type

  ! dataset types
  PetscInt, parameter :: DATASET_INTEGER = 1
  PetscInt, parameter :: DATASET_REAL = 2
  
  PetscInt, parameter :: DIM_NULL = 0
  PetscInt, parameter :: DIM_X = 1
  PetscInt, parameter :: DIM_Y = 2
  PetscInt, parameter :: DIM_XY = 3
    
  public :: DatasetCreate, &
            DatasetRead, &
            DatasetAddToList, &
            DatasetGetPointer, &
            DatasetInterpolateReal, &
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
  dataset%h5_dataset_name = ''
  dataset%filename = ''
  dataset%realization_dependent = PETSC_FALSE
  dataset%data_type = 0
  dataset%ndims = 0
  dataset%rmax = -1.d20
  dataset%rmin = 1.d20
  dataset%data_dim = DIM_NULL
  nullify(dataset%time_array)
  nullify(dataset%iarray)
  nullify(dataset%rarray)
  nullify(dataset%dims)
  nullify(dataset%discretization)
  nullify(dataset%origin)
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
      case('HDF5_DATASET_NAME') 
        call InputReadWord(input,option,dataset%h5_dataset_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'hdf5_dataset_name','DATASET')
      case('FILENAME') 
        call InputReadNChars(input,option,dataset%filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','DATASET')
!TODO(geh): remove is here after 10/30/11
#if 0        
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
#endif
      case('REALIZATION_DEPENDENT')
        dataset%realization_dependent = PETSC_TRUE
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in dataset'    
        call printErrMsg(option)
    end select 
  
  enddo
  
  if (len_trim(dataset%h5_dataset_name) < 1) then
    dataset%h5_dataset_name = dataset%name
  endif

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
! DatasetGetPointer: Returns the pointer to the dataset named "name"
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
             trim(debug_string) // '" not found among available datasets.'
    call printErrMsg(option)    
  endif

end function DatasetGetPointer

! ************************************************************************** !
!
! DatasetInterpolate: Interpolates data from the dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetInterpolateReal(dataset,x,y,z,time,real_value,option)

  use Utility_module, only : InterpolateBilinear
  use Option_module
  
  implicit none
  
  type(dataset_type) :: dataset
  PetscReal :: x, y, z, time
  PetscReal :: real_value
  type(option_type) :: option
  
  PetscInt :: interpolation_method
  PetscInt :: i, j, k
  PetscInt :: direction
  PetscReal :: x1,x2,y1,y2,v1,v2,v3,v4
  PetscInt :: index, nx
  PetscReal :: dx, dy
  
  call DatasetGetIJK(dataset,x,y,z,i,j,k)
  
  interpolation_method = INTERPOLATION_LINEAR
  !TODO(geh): fix
  dataset%data_dim = DIM_XY
  
  select case(interpolation_method)
    case(INTERPOLATION_STEP)
    case(INTERPOLATION_LINEAR)
      select case(dataset%data_dim)
        case(DIM_X,DIM_Y)
          select case(dataset%data_dim)
            case(DIM_X)
              direction = X_DIRECTION
            case(DIM_Y)
              direction = Y_DIRECTION
              i = j
          end select              
          ! why not have x return the weight?
          if (i <= 0 .or. i+1 > dataset%dims(direction)) then 
            select case(dataset%data_dim)
              case(X_DIRECTION)
                option%io_buffer = 'Out of x bounds'
              case(Y_DIRECTION)
                option%io_buffer = 'Out of y bounds'
            end select
            call printErrMsg(option)
          endif
          dx = dataset%discretization(direction)
          x1 = dataset%origin(direction) + (i-1)*dx
          v1 = dataset%rarray(i)
          v2 = dataset%rarray(i+1)
          real_value = v1 + (x-x1)/dx*(v2-v1)
        case(DIM_XY)
          if (i <= 0 .or. i+1 > dataset%dims(X_DIRECTION)) then
            option%io_buffer = 'Out of x bounds'
            call printErrMsg(option)
          endif
          if (j <= 0 .or. j+1 > dataset%dims(Y_DIRECTION)) then
            option%io_buffer = 'Out of y bounds'
            call printErrMsg(option)
          endif

          nx = dataset%dims(X_DIRECTION)
          dx = dataset%discretization(X_DIRECTION)
          dy = dataset%discretization(Y_DIRECTION)

          x1 = dataset%origin(X_DIRECTION) + (i-1)*dx
          x2 = x1 + dx
          
          index = i + (j-1)*nx
          v1 = dataset%rarray(index)
          v2 = dataset%rarray(index+1)
          
          y1 = dataset%origin(Y_DIRECTION) + (j-1)*dy
          y2 = y1 + dy
          
           ! really (j-1+1)
          index = i + j*nx
          v3 = dataset%rarray(index)
          v4 = dataset%rarray(index+1)
          
          real_value = InterpolateBilinear(x,y,x1,x2,y1,y2,v1,v2,v3,v4)
      end select
  end select
  
end subroutine DatasetInterpolateReal

! ************************************************************************** !
!
! DatasetGetIJK: Returns bounding i, j, k indices for point in dataset
! author: Glenn Hammond
! date: 10/26/11
!
! ************************************************************************** !
subroutine DatasetGetIJK(dataset,x,y,z,i,j,k)

  implicit none

  type(dataset_type) :: dataset
  PetscReal :: x, y, z
  PetscInt :: i, j, k
  
  ! get i, j, k indices in dataset
  if (dataset%dims(X_DIRECTION) > 1) then 
    i = int((x - dataset%origin(X_DIRECTION))/ &
            dataset%discretization(X_DIRECTION) + 1.d0)
  else
    i = 1
  endif
  if (dataset%dims(Y_DIRECTION) > 1) then 
    j = int((y - dataset%origin(Y_DIRECTION))/ &
            dataset%discretization(Y_DIRECTION) + 1.d0)
  else
    j = 1
  endif  
  !TODO(geh): set up k
  
end subroutine DatasetGetIJK

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

  if (associated(dataset%time_array)) deallocate(dataset%time_array)
  nullify(dataset%time_array)
  if (associated(dataset%iarray)) deallocate(dataset%iarray)
  nullify(dataset%iarray)
  if (associated(dataset%rarray)) deallocate(dataset%rarray)
  nullify(dataset%rarray)
  if (associated(dataset%dims)) deallocate(dataset%dims)
  nullify(dataset%dims)
  if (associated(dataset%discretization)) deallocate(dataset%discretization)
  nullify(dataset%discretization)
  if (associated(dataset%discretization)) deallocate(dataset%discretization)
  nullify(dataset%discretization)
  
  deallocate(dataset)
  nullify(dataset)
  
end subroutine DatasetDestroy

end module Dataset_Aux_module
