module Level_module

#include "finclude/petscsys.h"
  use petscsys
  use Patch_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: level_type 
    
    PetscInt :: id
    type(patch_list_type), pointer :: patch_list
    type(level_type), pointer :: next

  end type level_type

  ! pointer data structure required for making an array of level pointers in F90
  type, public :: level_ptr_type
    type(level_type), pointer :: ptr           ! pointer to the level_type
  end type level_ptr_type 

  type, public :: level_list_type
    PetscInt :: num_level_objects
    type(level_type), pointer :: first
    type(level_type), pointer :: last
    type(level_ptr_type), pointer :: array(:)
  end type level_list_type
    
  public :: LevelCreate, LevelDestroy, LevelCreateList, LevelDestroyList, &
            LevelAddToList, LevelConvertListToArray

contains

! ************************************************************************** !

function LevelCreate()
  ! 
  ! Allocates and initializes a new Level object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none
  
  type(level_type), pointer :: LevelCreate
  
  type(level_type), pointer :: level
  
  allocate(level)
  level%id = 0
  level%patch_list => PatchCreateList()
  
  nullify(level%next)
  
  LevelCreate => level
  
end function LevelCreate

! ************************************************************************** !

function LevelCreateList()
  ! 
  ! LevelListCreate: Creates a level list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none

  type(level_list_type), pointer :: LevelCreateList

  type(level_list_type), pointer :: level_list
  
  allocate(level_list)
  nullify(level_list%first)
  nullify(level_list%last)
  nullify(level_list%array)
  level_list%num_level_objects = 0

  LevelCreateList => level_list

end function LevelCreateList

! ************************************************************************** !

subroutine LevelAddPatch(new_patch,level)
  ! 
  ! Adds a new level to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none
  
  type(patch_type), pointer :: new_patch
  type(level_type) :: level

  call PatchAddToList(new_patch,level%patch_list)
  
end subroutine LevelAddPatch

! ************************************************************************** !

subroutine LevelAddToList(new_level,level_list)
  ! 
  ! Adds a new level to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none
  
  type(level_type), pointer :: new_level
  type(level_list_type) :: level_list
  
  level_list%num_level_objects = level_list%num_level_objects + 1
  new_level%id = level_list%num_level_objects
  if (.not.associated(level_list%first)) level_list%first => new_level
  if (associated(level_list%last)) level_list%last%next => new_level
  level_list%last => new_level
  
end subroutine LevelAddToList

! ************************************************************************** !

subroutine LevelConvertListToArray(level_list)
  ! 
  ! Creates an array of pointers to the
  ! levels in the level list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none
  
  type(level_list_type) :: level_list
    
  PetscInt :: count
  type(level_type), pointer :: cur_level
  
  
  allocate(level_list%array(level_list%num_level_objects))
  
  cur_level => level_list%first
  do 
    if (.not.associated(cur_level)) exit
    level_list%array(cur_level%id)%ptr => cur_level
    cur_level => cur_level%next
  enddo

end subroutine LevelConvertListToArray

! ************************************************************************** !

subroutine LevelDestroyList(level_list)
  ! 
  ! Deallocates a level list and array of leveles
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  ! 

  implicit none
  
  type(level_list_type), pointer :: level_list
    
  type(level_type), pointer :: cur_level, prev_level
  
  if (.not.associated(level_list)) return
  
  if (associated(level_list%array)) deallocate(level_list%array)
  nullify(level_list%array)
  
  cur_level => level_list%first
  do 
    if (.not.associated(cur_level)) exit
    prev_level => cur_level
    cur_level => cur_level%next
    call LevelDestroy(prev_level)
  enddo
  
  nullify(level_list%first)
  nullify(level_list%last)
  level_list%num_level_objects = 0
  
  deallocate(level_list)
  nullify(level_list)

end subroutine LevelDestroyList

! ************************************************************************** !

subroutine LevelDestroy(level)
  ! 
  ! Deallocates a level object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  implicit none
  
  type(level_type), pointer :: level
  
  call PatchDestroyList(level%patch_list)
      
  deallocate(level)
  nullify(level)
  
end subroutine LevelDestroy

end module Level_module
