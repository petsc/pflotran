module Strata_module

  use Region_module
  use Material_module
 
  implicit none

  private
 
#include "definitions.h"

 
  type, public :: strata_type
    PetscInt :: id                                       ! id of strata
    PetscTruth :: active
    character(len=MAXWORDLENGTH) :: material_property_name  ! character string defining name of material to be applied
    character(len=MAXSTRINGLENGTH) :: material_property_filename  ! character string defining name of file containing materia ids
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    PetscInt :: imaterial_property                       ! id of material in material array/list
    PetscInt :: iregion                                  ! id of region in region array/list
    type(material_property_type), pointer :: material_property ! pointer to material in material array/list
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
  strata%active = PETSC_TRUE
  strata%material_property_name = ""
  strata%material_property_filename = ""
  strata%region_name = ""
  strata%iregion = 0
  strata%imaterial_property = 0

  nullify(strata%region)
  nullify(strata%material_property)
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
  new_strata%material_property_name = strata%material_property_name
  new_strata%material_property_filename = strata%material_property_filename
  new_strata%region_name = strata%region_name
  new_strata%iregion = strata%iregion

  ! keep these null
  nullify(new_strata%region)
  nullify(new_strata%material_property)
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
subroutine StrataRead(strata,input,option)

  use Input_module
  use Option_module
  
  implicit none
  
  type(strata_type) :: strata
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadFlotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','STRATA')   
      
    select case(trim(keyword))
    
      case('REGION')
        call InputReadWord(input,option,strata%region_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'region name','STRATA')
      case('MATERIAL')
        call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'material property name','STRATA')
        strata%material_property_name = string
        strata%material_property_filename = string
      case('INACTIVE')
        strata%active = PETSC_FALSE
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
