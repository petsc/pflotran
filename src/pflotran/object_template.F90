module Object_module

  ! add accessory modules here

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
  
  private

  type, public :: object_type
  
    ! add object members here
!   type(member_object_type), pointer :: member_object

    ! add variables/arrays here
    PetscInt :: variable
    PetscReal, pointer :: array(:)

  end type object_type
  
  ! add interface definitions here
! interface ObjectCreate
!   module procedure ObjectCreate1
!   module procedure ObjectCreate2
! end interface
  
  ! add public definitions here
  public :: ObjectCreate, &
            ObjectDestroy
  
contains

! ************************************************************************** !
!
! ObjectCreate: Allocates and initializes a new XXX object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
function ObjectCreate()

  implicit none
  
  type(object_type), pointer :: ObjectCreate
  
  type(object_type), pointer :: object
  
  allocate(object)
  ! create object member objects here
! object%member_object => MemberObjectCreate()
  ! nullify variables/arrays here
  object%variable = 0
  nullify(object%array)
  
  ObjectCreate => object  
  
end function ObjectCreate

! ************************************************************************** !
!
! ObjectDestroy: Deallocates an XXX object
! author: Glenn Hammond
! date: 03/14/13
!
! ************************************************************************** !
subroutine ObjectDestroy(object)

  use Utility_module, only: DeallocateArray 

  implicit none
  
  type(object_type), pointer :: object
  
  if (.not.associated(object)) return

  ! destroy member objects
! call MemberObjectDestroy(object%member_object)

  ! deallocate arrays
  call DeallocateArray(object%array)
  
  ! deallocate and nullify the object itself
  deallocate(object)
  nullify(object)
  
end subroutine ObjectDestroy
  
end module Object_module
