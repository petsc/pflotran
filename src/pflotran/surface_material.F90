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
            SurfaceMaterialPropertyRead

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

end module Surface_Material_module

#endif ! SURFACE_FLOW