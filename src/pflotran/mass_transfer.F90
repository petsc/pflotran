module Mass_Transfer_module
 
  use Dataset_Global_class
  
  implicit none

  private

#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
 
  type, public :: mass_transfer_type
    PetscInt :: idof
    character(len=MAXSTRINGLENGTH) :: filename
    character(len=MAXWORDLENGTH) :: dataset_name
    class(dataset_global_type), pointer :: dataset
    Vec :: vec
    type(mass_transfer_type), pointer :: next
  end type mass_transfer_type
  
  public :: MassTransferCreate, MassTransferDestroy, &
            MassTransferRead, MassTransferAddToList, &
            MassTransferUpdate

contains

! ************************************************************************** !
!
! MassTransferCreate: Creates a mass transfer object
! author: Glenn Hammond
! date: 05/01/13
!
! ************************************************************************** !
function MassTransferCreate()
  
  implicit none

  type(mass_transfer_type), pointer :: MassTransferCreate
  
  type(mass_transfer_type), pointer :: mass_transfer
  
  allocate(mass_transfer)
  mass_transfer%idof = 0
  mass_transfer%filename = ''
  mass_transfer%dataset_name = ''
  nullify(mass_transfer%dataset)
  nullify(mass_transfer%next)
  mass_transfer%vec = 0
  MassTransferCreate => mass_transfer

end function MassTransferCreate

! ************************************************************************** !
!
! MassTransferRead: Reads in contents of a mass transfer card
! author: Glenn Hammond
! date: 05/01/13
! 
! ************************************************************************** !
subroutine MassTransferRead(mass_transfer,input,option)

  use Option_module
  use Input_module
  use String_module

  implicit none
  
  type(mass_transfer_type) :: mass_transfer
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','MASS_TRANSFER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('IDOF') 
        call InputReadInt(input,option,mass_transfer%idof)
        call InputErrorMsg(input,option,'idof','MASS_TRANSFER')
      case('DATASET')
        call InputReadNChars(input,option, &
                             mass_transfer%dataset_name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'DATASET,NAME','MASS_TRANSFER')
      case('FILENAME') 
        call InputReadNChars(input,option,mass_transfer%filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'FILENAME','MASS_TRANSFER')        
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
                           ' not recognized in mass transfer'    
        call printErrMsg(option)
    end select
    
  enddo  

end subroutine MassTransferRead

! ************************************************************************** !
!
! MassTransferAddToList: Adds a mass transfer object to linked list
! author: Glenn Hammond
! date: 05/01/13
!
! ************************************************************************** !
subroutine MassTransferAddToList(mass_transfer,list)

  implicit none
  
  type(mass_transfer_type), pointer :: mass_transfer
  type(mass_transfer_type), pointer :: list

  type(mass_transfer_type), pointer :: cur_mass_transfer
  
  if (associated(list)) then
    cur_mass_transfer => list
    ! loop to end of list
    do
      if (.not.associated(cur_mass_transfer%next)) exit
      cur_mass_transfer => cur_mass_transfer%next
    enddo
    cur_mass_transfer%next => mass_transfer
  else
    list => mass_transfer
  endif
  
end subroutine MassTransferAddToList

! ************************************************************************** !
!
! MassTransferUpdate: Updates a mass transfer object transfering data from
!                     the buffer into the PETSc Vec
! author: Glenn Hammond
! date: 05/01/13
!
! ************************************************************************** !
recursive subroutine MassTransferUpdate(mass_transfer, discretization, &
                                        grid, option)
  use Discretization_module
  use Grid_module
  use Option_module
  
  implicit none
  
  type(mass_transfer_type), pointer :: mass_transfer
  type(discretization_type) :: discretization
  type(grid_type) :: grid
  type(option_type) :: option  
  PetscReal :: time
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return
  
  if (.not.associated(mass_transfer%dataset)) then
    mass_transfer%dataset => DatasetGlobalCreate()
    mass_transfer%dataset%filename = mass_transfer%filename
    mass_transfer%dataset%dataset_name = mass_transfer%dataset_name
  endif
!  call mass_transfer%dataset%Load(discretization,grid,option)
  call DatasetGlobalLoad(mass_transfer%dataset,discretization,grid,option)

  call VecGetArrayF90(mass_transfer%vec,vec_ptr,ierr)
  vec_ptr(:) = mass_transfer%dataset%rarray(:)
  call VecRestoreArrayF90(mass_transfer%vec,vec_ptr,ierr)
  
end subroutine MassTransferUpdate

! ************************************************************************** !
!
! MassTransferDestroy: Destroys a mass transfer object
! author: Glenn Hammond
! date: 05/01/13
!
! ************************************************************************** !
recursive subroutine MassTransferDestroy(mass_transfer)

  implicit none
  
  type(mass_transfer_type), pointer :: mass_transfer
  
  PetscErrorCode :: ierr
  
  if (.not.associated(mass_transfer)) return
  
  ! Simply nullify the pointer as the dataset resides in a list to be
  ! destroyed separately.
!  nullify(mass_transfer%dataset)
  call DatasetGlobalDestroy(mass_transfer%dataset)
  if (mass_transfer%vec /= 0) &
    call VecDestroy(mass_transfer%vec ,ierr)
  call MassTransferDestroy(mass_transfer%next)

  deallocate(mass_transfer)
  nullify(mass_transfer)
  
end subroutine MassTransferDestroy

end module Mass_Transfer_module
