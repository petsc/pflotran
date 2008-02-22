module Strata_module

  use Region_module
  use Material_module
 
  implicit none

  private
 
#include "definitions.h"

 
  type, public :: strata_type
    PetscInt :: id                                       ! id of strata
    logical :: active
    character(len=MAXWORDLENGTH) :: material_name       ! character string defining name of material to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    PetscInt :: imaterial                                ! id of material in material array/list
    PetscInt :: iregion                                  ! id of region in region array/list
    type(material_type), pointer :: material            ! pointer to material in material array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(strata_type), pointer :: next            ! pointer to next strata
  end type strata_type
  
  type, public :: strata_ptr_type
    type(strata_type), pointer :: ptr
  end type strata_ptr_type
    
  type, public :: strata_list_type
    PetscInt :: num_strata
    type(strata_type), pointer :: first
    type(strata_type), pointer :: last
    type(strata_ptr_type), pointer :: array(:)    
  end type strata_list_type
  
  interface StrataCreate
    module procedure StrataCreate1
    module procedure StrataCreateFromStrata
  end interface
  
  public :: StrataCreate, StrataDestroy, StrataInitList, &
            StrataAddToList, StrataRead, StrataDestroyList
  
contains

! ************************************************************************** !
!
! StrataCreate1: Creates a strata
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function StrataCreate1()

  implicit none

  type(strata_type), pointer :: StrataCreate1
  
  type(strata_type), pointer :: strata
  
  allocate(strata)
  strata%id = 0
  strata%active = .true.
  strata%material_name = ""
  strata%region_name = ""
  strata%iregion = 0

  nullify(strata%region)
  nullify(strata%material)
  nullify(strata%next)
  
  StrataCreate1 => strata

end function StrataCreate1

! ************************************************************************** !
!
! StrataCreateFromStrata: Creates a strata
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function StrataCreateFromStrata(strata)

  implicit none

  type(strata_type), pointer :: StrataCreateFromStrata
  type(strata_type), pointer :: strata

  type(strata_type), pointer :: new_strata
  
  new_strata => StrataCreate1()
  
  new_strata%id = strata%id
  new_strata%active = strata%active
  new_strata%material_name = strata%material_name
  new_strata%region_name = strata%region_name
  new_strata%iregion = strata%iregion

  ! keep these null
  nullify(new_strata%region)
  nullify(new_strata%material)
  nullify(new_strata%next)
  
  StrataCreateFromStrata => new_strata

end function StrataCreateFromStrata

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
subroutine StrataRead(strata,fid,option)

  use Fileio_module
  use Option_module
  
  implicit none
  
  type(strata_type) :: strata
  PetscInt :: fid
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string, error_string
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  PetscReal :: value
  PetscInt :: count1, material_file_id = 86
  PetscErrorCode :: ierr

  ierr = 0
  do
  
    call fiReadFlotranString(IUNIT1,string,ierr)
    
    if (string(1:1) == '.' .or. string(1:1) == '/' .or. &
        fiStringCompare(string,'END',THREE_INTEGER)) exit  

    call fiReadWord(string,keyword,.true.,ierr)
    call fiErrorMsg(option%myrank,'keyword','STRATA', ierr)   
      
    select case(trim(keyword))
    
      case('REGION')
        call fiReadWord(string,strata%region_name,.true.,ierr)
        call fiErrorMsg(option%myrank,'region name','STRATA', ierr)
      case('MATERIAL')
        call fiReadWord(string,word,.true.,ierr)
        call fiErrorMsg(option%myrank,'material name','STRATA', ierr)
        strata%material_name = word
      case('INACTIVE')
        strata%active = .false.
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
