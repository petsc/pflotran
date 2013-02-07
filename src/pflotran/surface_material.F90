#ifdef SURFACE_FLOW

module Surface_Material_module

  implicit none
  
  private
  
#include "definitions.h"

  type, public :: surface_material_property_type
    
    character(len=MAXWORDLENGTH) :: name
    PetscInt                     :: id
    PetscReal                    :: mannings
    
    type(surface_material_property_type), pointer :: next
  end type surface_material_property_type
  
  type, public :: surface_material_property_ptr_type
    type(surface_material_property_type), pointer :: ptr
  end type surface_material_property_ptr_type

  public :: SurfaceMaterialPropertyCreate, &
            SurfaceMaterialPropertyDestroy, &
            SurfaceMaterialPropertyAddToList, &
            SurfaceMaterialPropertyRead, &
            SurfaceMaterialPropConvertListToArray, &
            SurfaceMaterialPropGetPtrFromArray

  contains
  
! ************************************************************************** !
!> This routine creates a surface material property
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
function SurfaceMaterialPropertyCreate()

  implicit none
  
  type(surface_material_property_type), pointer :: SurfaceMaterialPropertyCreate
  type(surface_material_property_type), pointer :: surf_material_property
  
  allocate(surf_material_property)

  surf_material_property%name      = ''
  surf_material_property%id        = 0
  surf_material_property%mannings  = 0.d0
  
  nullify(surf_material_property%next)
  
  SurfaceMaterialPropertyCreate => surf_material_property

end function SurfaceMaterialPropertyCreate

! ************************************************************************** !
!> This routine reads in contents of a surface material property
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceMaterialPropertyRead(surf_material_property,input,option)

  use Option_module
  use Input_module
  use String_module
  
  implicit none
  
  type(surface_material_property_type) :: surf_material_property
  type(input_type)                     :: input
  type(option_type)                    :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string

  do
    call InputReadFlotranString(input,option)
    
    if(InputCheckExit(input,option)) exit
  
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SURFACE_MATERIAL_PROPERTY')
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('ID')
        call InputReadInt(input,option,surf_material_property%id)
        call InputErrorMsg(input,option,'id','SURFACE_MATERIAL_PROPERTY')
      case('MANNINGS')
        call InputReadDouble(input,option,surf_material_property%mannings)
        call InputErrorMsg(input,option,'MANNINGS','SURFACE_MATERIAL_PROPERTY')
      case default
        option%io_buffer = 'Keyword: ' // trim(keyword) // &
          ' not recognized in surface_material_property'
        call printErrMsg(option)
      end select
  enddo
  
end subroutine SurfaceMaterialPropertyRead

! ************************************************************************** !
!> This routine adds a surface material property to a linked list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
subroutine SurfaceMaterialPropertyAddToList(surf_material_property,list)

  implicit none
  
  type(surface_material_property_type), pointer :: surf_material_property
  type(surface_material_property_type), pointer :: list
  type(surface_material_property_type), pointer :: cur_surf_material_property
  
  if (associated(list)) then
    cur_surf_material_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_surf_material_property%next)) exit
      cur_surf_material_property => cur_surf_material_property%next
    enddo
    cur_surf_material_property%next => surf_material_property
  else
    list => surf_material_property
  endif
  
end subroutine SurfaceMaterialPropertyAddToList

! ************************************************************************** !
!> This routine destroys a surface material property
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/09/12
! ************************************************************************** !
recursive subroutine SurfaceMaterialPropertyDestroy(surf_material_property)

  implicit none
  
  type(surface_material_property_type), pointer :: surf_material_property
  
  if(.not.associated(surf_material_property)) return
  
  call SurfaceMaterialPropertyDestroy(surf_material_property%next)
  
  deallocate(surf_material_property)
  nullify(surf_material_property)
  
end subroutine SurfaceMaterialPropertyDestroy

! ************************************************************************** !
!> This routine creates an array of pointers to the surface_material_properties
!! in the list (similar to subroutine MaterialPropConvertListToArray)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/11/12
! ************************************************************************** !
subroutine SurfaceMaterialPropConvertListToArray(list,array,option)

  use Option_module
  use String_module

  implicit none

  type(surface_material_property_type), pointer :: list
  type(surface_material_property_ptr_type), pointer :: array(:)
  type(option_type) :: option

  type(surface_material_property_type), pointer :: cur_material_property
  type(surface_material_property_type), pointer :: prev_material_property
  type(surface_material_property_type), pointer :: next_material_property
  PetscInt :: i, j, length1,length2, max_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

  max_id = 0
  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    max_id = max(max_id,cur_material_property%id)
    cur_material_property => cur_material_property%next
  enddo

  allocate(array(max_id))
  do i = 1, max_id
    nullify(array(i)%ptr)
  enddo

  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_id))
  id_count = 0

  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    id_count(cur_material_property%id) = &
      id_count(cur_material_property%id) + 1
    array(cur_material_property%id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo

  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Material ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call printMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Material IDs.'
    call printErrMsg(option)
  endif

  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, max_id
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Material name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call printMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Material names.'
    call printErrMsg(option)
  endif

end subroutine SurfaceMaterialPropConvertListToArray

! ************************************************************************** !
!> This routine returns a pointer to the surface material property matching
!! surface_material_propertry_name (similar to subroutine
!! MaterialPropGetPtrFromArray)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/11/12
! ************************************************************************** !
function SurfaceMaterialPropGetPtrFromArray(surf_material_property_name, &
                                            surf_material_property_array)

  use String_module

  implicit none

  type(surface_material_property_type), pointer     :: SurfaceMaterialPropGetPtrFromArray
  type(surface_material_property_ptr_type), pointer :: surf_material_property_array(:)
  character(len=MAXWORDLENGTH)                      :: surf_material_property_name
  PetscInt :: length
  PetscInt :: isurf_material_property

  nullify(SurfaceMaterialPropGetPtrFromArray)

  do isurf_material_property = 1, size(surf_material_property_array)
    length = len_trim(surf_material_property_name)
    if (.not.associated(surf_material_property_array(isurf_material_property)%ptr)) cycle
    if (length == &
        len_trim(surf_material_property_array(isurf_material_property)%ptr%name) .and. &
        StringCompare(surf_material_property_array(isurf_material_property)%ptr%name, &
                        surf_material_property_name,length)) then
      SurfaceMaterialPropGetPtrFromArray => &
        surf_material_property_array(isurf_material_property)%ptr
      return
    endif
  enddo

end function SurfaceMaterialPropGetPtrFromArray

end module Surface_Material_module

#endif
! SURFACE_FLOW