module Breakthrough_module

  use Region_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: breakthrough_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: region_name
    type(region_type), pointer :: region
    type(breakthrough_type), pointer :: next
  end type breakthrough_type
  
  type, public :: breakthrough_list_type
    PetscInt :: num_breakthroughs
    type(breakthrough_type), pointer :: first
    type(breakthrough_type), pointer :: last
    type(breakthrough_type), pointer :: array(:)
  end type breakthrough_list_type

  public :: BreakthroughCreate, BreakthroughDestroy, BreakthroughRead, &
            BreakthroughAddToList, BreakthroughInitList, BreakthroughDestroyList, &
            BreakthroughGetPtrFromList
  
contains

! ************************************************************************** !
!
! BreakthroughCreate: Create object that stores breakthrough regions
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function BreakthroughCreate()

  implicit none
  
  type(breakthrough_type), pointer :: BreakthroughCreate
  
  type(breakthrough_type), pointer :: breakthrough
  
  allocate(breakthrough)
  
  breakthrough%name = ""
  breakthrough%region_name = ""
  breakthrough%id = 0
  nullify(breakthrough%region)
  nullify(breakthrough%next)
  
  BreakthroughCreate => breakthrough

end function BreakthroughCreate

! ************************************************************************** !
!
! BreakthroughRead: Reads breakthrough data from the input file
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine BreakthroughRead(breakthrough,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(breakthrough_type) :: breakthrough
  PetscInt :: fid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    
    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','BREAKTHROUGH', ierr)   
      
    select case(trim(keyword))
    
      case('REGION')
        call fiReadWord(string,breakthrough%region_name,.true.,ierr)
        call fiErrorMsg(option%myrank,'region name','BREAKTHROUGH', ierr)
    end select 
  
  enddo  

end subroutine BreakthroughRead

! ************************************************************************** !
!
! BreakthroughInitList: Initializes a breakthrough list
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine BreakthroughInitList(list)

  implicit none

  type(breakthrough_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_breakthroughs = 0

end subroutine BreakthroughInitList

! ************************************************************************** !
!
! BreakthroughAddToList: Adds a new breakthrough to a breakthrough list
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine BreakthroughAddToList(new_breakthrough,list)

  implicit none
  
  type(breakthrough_type), pointer :: new_breakthrough
  type(breakthrough_list_type) :: list
  
  list%num_breakthroughs = list%num_breakthroughs + 1
  new_breakthrough%id = list%num_breakthroughs
  if (.not.associated(list%first)) list%first => new_breakthrough
  if (associated(list%last)) list%last%next => new_breakthrough
  list%last => new_breakthrough
  
end subroutine BreakthroughAddToList

! ************************************************************************** !
!
! BreakthroughGetPtrFromList: Returns a pointer to the breaktrough matching &
!                             breakthrough_name
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function BreakthroughGetPtrFromList(breakthrough_name,breakthrough_list)

  use Fileio_module

  implicit none
  
  type(breakthrough_type), pointer :: BreakthroughGetPtrFromList
  character(len=MAXNAMELENGTH) :: breakthrough_name
  type(breakthrough_list_type) :: breakthrough_list
 
  PetscInt :: length
  type(breakthrough_type), pointer :: breakthrough
    
  nullify(BreakthroughGetPtrFromList)
  breakthrough => breakthrough_list%first
  
  do 
    if (.not.associated(breakthrough)) exit
    length = len_trim(breakthrough_name)
    if (length == len_trim(breakthrough%name) .and. &
        fiStringCompare(breakthrough%name,breakthrough_name, &
                        length)) then
      BreakthroughGetPtrFromList => breakthrough
      return
    endif
    breakthrough => breakthrough%next
  enddo
  
end function BreakthroughGetPtrFromList

! ************************************************************************** !
!
! BreakthroughDestroyList: Deallocates a list of breakthroughs
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
subroutine BreakthroughDestroyList(breakthrough_list)

  implicit none
  
  type(breakthrough_list_type), pointer :: breakthrough_list
  
  type(breakthrough_type), pointer :: breakthrough, prev_breakthrough
  
  if (.not.associated(breakthrough_list)) return
  
  breakthrough => breakthrough_list%first
  do 
    if (.not.associated(breakthrough)) exit
    prev_breakthrough => breakthrough
    breakthrough => breakthrough%next
    call BreakthroughDestroy(prev_breakthrough)
  enddo
  
  breakthrough_list%num_breakthroughs = 0
  nullify(breakthrough_list%first)
  nullify(breakthrough_list%last)
  if (associated(breakthrough_list%array)) deallocate(breakthrough_list%array)
  nullify(breakthrough_list%array)
  
  deallocate(breakthrough_list)
  nullify(breakthrough_list)

end subroutine BreakthroughDestroyList

! ************************************************************************** !
!
! BreakthroughDestroy: Deallocates a breakthrough
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine BreakthroughDestroy(breakthrough)

  implicit none
  
  type(breakthrough_type), pointer :: breakthrough
  
  PetscInt :: i
  
  if (.not.associated(breakthrough)) return
  
  nullify(breakthrough%region)
  deallocate(breakthrough)
  nullify(breakthrough)

end subroutine BreakthroughDestroy

end module Breakthrough_module
