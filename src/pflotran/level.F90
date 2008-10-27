module Level_module

  use Patch_module
 
  implicit none

  private

#include "definitions.h"

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
!
! LevelCreate: Allocates and initializes a new Level object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function LevelCreate()

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
!
! LevelListCreate: Creates a level list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function LevelCreateList()

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
!
! LevelAddPatch: Adds a new level to list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelAddPatch(new_patch,level)

  implicit none
  
  type(patch_type), pointer :: new_patch
  type(level_type) :: level

  call PatchAddToList(new_patch,level%patch_list)
  
end subroutine LevelAddPatch

! ************************************************************************** !
!
! LevelAddToList: Adds a new level to list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelAddToList(new_level,level_list)

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
!
! LevelConvertListToArray: Creates an array of pointers to the 
!                               levels in the level list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelConvertListToArray(level_list)

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
#if 0
! ************************************************************************** !
!
! LevelProcessCouplers: Deallocates a realization
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelProcessCouplers(level,flow_conditions,transport_conditions, &
                                material_array,option)

  use Option_module
  use Material_module
  use Condition_module
  
  implicit none
  
  type(level_type) :: level
  type(material_ptr_type), pointer :: material(:)
  type(condition_list_type) :: flow_conditions
  type(condition_list_type) :: transport_conditions
  type(option_type) :: option  
  
  type(patch_type), pointer :: cur_patch

  cur_patch => level%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    call PatchProcessCouplers(cur_patch,flow_conditions,transport_conditions, &
                              material_array,option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine LevelProcessCouplers

! ************************************************************************** !
!
! LevelInitAllCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelInitAllCouplerAuxVars(level,option)

  use Option_module
  
  implicit none
  
  type(level_type) :: level
  type(option_type) :: option
  
  type(patch_type), pointer :: cur_patch

  cur_patch => level%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    call PatchInitAllCouplerAuxVars(cur_patch,option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine LevelInitAllCouplerAuxVars

! ************************************************************************** !
!
! LevelInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelInitCouplerAuxVars(level,option)

  use Option_module
  
  implicit none
  
  type(level_type) :: level
  type(option_type) :: option
  
  type(patch_type), pointer :: cur_patch

  cur_patch => level%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    call PatchInitCouplerAuxVars(cur_patch,option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine LevelInitCouplerAuxVars

! ************************************************************************** !
!
! LevelUpdateAllCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in lis
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelUpdateAllCouplerAuxVars(level,force_update_flag,option)

  use Option_module

  implicit none
  
  type(level_type) :: level
  logical :: force_update_flag
  type(option_type) :: option
    
  type(patch_type), pointer :: cur_patch

  cur_patch => level%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    call PatchUpdateAllCouplerAuxVars(cur_patch,force_update_flag,option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine LevelUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! LevelUpdateCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in lis
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelUpdateCouplerAuxVars(level,force_update_flag,option)

  use Option_module

  implicit none
  
  type(level_type) :: level
  logical :: force_update_flag
  type(option_type) :: option
    
  type(patch_type), pointer :: cur_patch

  cur_patch => level%patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    call PatchUpdateCouplerAuxVars(cur_patch,force_update_flag,option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine LevelUpdateCouplerAuxVars
#endif
! ************************************************************************** !
!
! LevelDestroyList: Deallocates a level list and array of leveles
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine LevelDestroyList(level_list)

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
!
! LevelDestroy: Deallocates a level object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine LevelDestroy(level)

  implicit none
  
  type(level_type), pointer :: level
  
  call PatchDestroyList(level%patch_list)
      
  deallocate(level)
  nullify(level)
  
end subroutine LevelDestroy

end module Level_module
