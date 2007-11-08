module Strata_module

  use Region_module
  use Material_module
 
  implicit none

  private
 
#include "definitions.h"
 
  type, public :: strata_type
    integer :: id                                       ! id of strata
    character(len=MAXWORDLENGTH) :: material_name       ! character string defining name of material to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    integer :: imaterial                                ! id of material in material array/list
    integer :: iregion                                  ! id of region in region array/list
    type(material_type), pointer :: material            ! pointer to material in material array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(strata_type), pointer :: next            ! pointer to next strata
  end type strata_type
  
  type, public :: strata_ptr_type
    type(strata_type), pointer :: ptr
  end type strata_ptr_type
    
  type, public :: strata_list_type
    integer :: num_strata
    type(strata_type), pointer :: first
    type(strata_type), pointer :: last
    type(strata_ptr_type), pointer :: array(:)    
  end type strata_list_type
  
  integer, save :: num_strata = 0
  
  public :: StrataCreate, StrataDestroy, StrataInitList, &
            StrataAddToList, StrataRead, StrataDestroyList
  
contains

! ************************************************************************** !
!
! StrataCreate: Creates a strata
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function StrataCreate()

  implicit none

  type(strata_type), pointer :: StrataCreate
  
  type(strata_type), pointer :: strata
  
  allocate(strata)
  strata%id = 0
  strata%material_name = ""
  strata%region_name = ""
  strata%iregion = 0

  nullify(strata%region)
  nullify(strata%next)
  
  num_strata = num_strata + 1
  strata%id = num_strata
  
  StrataCreate => strata

end function StrataCreate

! ************************************************************************** !
!
! StrataInitList: Initializes a strata list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StrataInitList(list)

  implicit none

  type(strata_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_strata = 0

end subroutine StrataInitList

! ************************************************************************** !
!
! StrataRead: Reads a strata from the input file
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StrataRead(strata,fid)

  use Fileio_module
  
  implicit none
  
  type(strata_type) :: strata
  integer :: fid
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  integer :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    if (ierr /= 0) exit

    call fiReadWord(string,word,.true.,ierr)
    call fiErrorMsg('keyword','UNIT', ierr)   
      
    select case(trim(word))
    
      case('REGION')
        call fiReadWord(string,strata%region_name,.true.,ierr)
      case('MATERIAL')
        call fiReadWord(string,strata%material_name,.true.,ierr)
      case('END')
        exit
    end select 
  
  enddo  

end subroutine StrataRead

! ************************************************************************** !
!
! StrataAddToList: Adds a new strata to a strata list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StrataAddToList(new_strata,list)

  implicit none
  
  type(strata_type), pointer :: new_strata
  type(strata_list_type) :: list
  
  list%num_strata = list%num_strata + 1
  new_strata%id = list%num_strata
  if (.not.associated(list%first)) list%first => new_strata
  if (associated(list%last)) list%last%next => new_strata
  list%last => new_strata
  
end subroutine StrataAddToList

! ************************************************************************** !
!
! StrataDestroyList: Deallocates a list of stratas
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine StrataDestroyList(strata_list)

  implicit none
  
  type(strata_list_type), pointer :: strata_list
  
  type(strata_type), pointer :: strata, prev_strata
  
  
  strata => strata_list%first
  do 
    if (.not.associated(strata)) exit
    prev_strata => strata
    strata => strata%next
    call StrataDestroy(prev_strata)
  enddo
  
  strata_list%num_strata = 0
  nullify(strata_list%first)
  nullify(strata_list%last)
  if (associated(strata_list%array)) deallocate(strata_list%array)
  nullify(strata_list%array)
  
  deallocate(strata_list)
  nullify(strata_list)

end subroutine StrataDestroyList

! ************************************************************************** !
!
! StrataDestroy: Destroys a strata
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine StrataDestroy(strata)

  implicit none
  
  type(strata_type), pointer :: strata
  
  if (.not.associated(strata)) return
  
  ! since strata%region is a pointer to a region in a list, nullify instead
  ! of destroying since the list will be destroyed separately
  nullify(strata%region)
  
  deallocate(strata)
  nullify(strata)

end subroutine StrataDestroy

end module Strata_module
