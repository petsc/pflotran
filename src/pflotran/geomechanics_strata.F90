#ifdef GEOMECH

module Geomechanics_Strata_module

  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use PFLOTRAN_Constants_module
 
  implicit none

  private
 
#include "finclude/petscsys.h"
 
  type, public :: geomech_strata_type
    PetscInt :: id                                                        ! id of strata
    PetscBool :: active
    character(len=MAXWORDLENGTH) :: material_property_name                ! character string defining name of material to be applied
    character(len=MAXSTRINGLENGTH) :: material_property_filename          ! character string defining name of file containing materia ids
    PetscBool :: realization_dependent
    character(len=MAXWORDLENGTH) :: region_name                           ! character string defining name of region to be applied
    PetscInt :: imaterial_property                                        ! id of material in material array/list
    PetscInt :: iregion                                                   ! id of region in region array/list
    type(geomech_material_property_type), pointer :: material_property    ! pointer to material in material array/list
    type(gm_region_type), pointer :: region                               ! pointer to region in region array/list
    type(geomech_strata_type), pointer :: next                            ! pointer to next strata
  end type geomech_strata_type
  
  type, public :: geomech_strata_ptr_type
    type(geomech_strata_type), pointer :: ptr
  end type geomech_strata_ptr_type
    
  type, public :: geomech_strata_list_type
    PetscInt :: num_strata
    type(geomech_strata_type), pointer :: first
    type(geomech_strata_type), pointer :: last
    type(geomech_strata_ptr_type), pointer :: array(:)    
  end type geomech_strata_list_type
  
  interface GeomechStrataCreate
    module procedure GeomechStrataCreate1
    module procedure GeomechStrataCreateFromGeomechStrata
  end interface
  
  public :: GeomechStrataCreate, &
            GeomechStrataDestroy, &
            GeomechStrataInitList, &
            GeomechStrataAddToList, &
            GeomechStrataRead, &
            GeomechStrataDestroyList
  
contains

! ************************************************************************** !
!
! GeomechStrataCreate1: Creates a geomechanics strata
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
function GeomechStrataCreate1()

  implicit none

  type(geomech_strata_type), pointer :: GeomechStrataCreate1
  
  type(geomech_strata_type), pointer :: strata
  
  allocate(strata)
  strata%id = 0
  strata%active = PETSC_TRUE
  strata%material_property_name = ""
  strata%material_property_filename = ""
  strata%realization_dependent = PETSC_FALSE
  strata%region_name = ""
  strata%iregion = 0
  strata%imaterial_property = 0

  nullify(strata%region)
  nullify(strata%material_property)
  nullify(strata%next)
  
  GeomechStrataCreate1 => strata

end function GeomechStrataCreate1

! ************************************************************************** !
!
! GeomechStrataCreateFromGeomechStrata: Creates a geomechanics strata
!                                       from another geomechanics strata
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
function GeomechStrataCreateFromGeomechStrata(strata)

  implicit none

  type(geomech_strata_type), pointer :: GeomechStrataCreateFromGeomechStrata
  type(geomech_strata_type), pointer :: strata

  type(geomech_strata_type), pointer :: new_strata
  
  new_strata => GeomechStrataCreate1()
  
  new_strata%id = strata%id
  new_strata%active = strata%active
  new_strata%material_property_name = strata%material_property_name
  new_strata%material_property_filename = strata%material_property_filename
  new_strata%realization_dependent = strata%realization_dependent
  new_strata%region_name = strata%region_name
  new_strata%iregion = strata%iregion
  ! keep these null
  nullify(new_strata%region)
  nullify(new_strata%material_property)
  nullify(new_strata%next)
  
  GeomechStrataCreateFromGeomechStrata => new_strata

end function GeomechStrataCreateFromGeomechStrata

! ************************************************************************** !
!
! GeomechStrataInitList: Initializes a geomechanics strata list
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
subroutine GeomechStrataInitList(list)

  implicit none

  type(geomech_strata_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_strata = 0

end subroutine GeomechStrataInitList

! ************************************************************************** !
!
! GeomechStrataRead: Reads a geomechanics strata from the input file
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
subroutine GeomechStrataRead(strata,input,option)

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(geomech_strata_type) :: strata
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS STRATA')   
      
    select case(trim(keyword))
    
      case('GEOMECHANICS_REGION')
        call InputReadWord(input,option,strata%region_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'geomechanics region name', &
                           'GEOMECHANICS STRATA')
      case('GEOMECHANICS_MATERIAL')
        call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option, &
                   'geomechancis material property name','GEOMECHANICS STRATA')
        if (StringCompareIgnoreCase(string,'realization_dependent')) then
          strata%realization_dependent = PETSC_TRUE
          call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_TRUE)
          call InputErrorMsg(input,option, &
                   'geomechancis material property name','GEOMECHANICS STRATA')
        endif
        strata%material_property_name = trim(string)
        strata%material_property_filename = string
      case('INACTIVE')
        strata%active = PETSC_FALSE
      case default
        option%io_buffer = 'Keyword "' // trim(keyword) // &
          '" in GEOMECHANICS STRATA block not recognized.'
        call printErrMsg(option)
    end select 
  
  enddo  

end subroutine GeomechStrataRead

! ************************************************************************** !
!
! GeomechStrataAddToList: Adds a new geomechanics strata to a geomechanics
!                         strata list
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
subroutine GeomechStrataAddToList(new_strata,list)

  implicit none
  
  type(geomech_strata_type), pointer :: new_strata
  type(geomech_strata_list_type) :: list
  
  list%num_strata = list%num_strata + 1
  new_strata%id = list%num_strata
  if (.not.associated(list%first)) list%first => new_strata
  if (associated(list%last)) list%last%next => new_strata
  list%last => new_strata
  
end subroutine GeomechStrataAddToList

! ************************************************************************** !
!
! GeomechStrataDestroyList: Deallocates a list of geomechanics stratas
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
subroutine GeomechStrataDestroyList(strata_list)

  implicit none
  
  type(geomech_strata_list_type), pointer :: strata_list
  
  type(geomech_strata_type), pointer :: strata, prev_strata
  
  
  strata => strata_list%first
  do 
    if (.not.associated(strata)) exit
    prev_strata => strata
    strata => strata%next
    call GeomechStrataDestroy(prev_strata)
  enddo
  
  strata_list%num_strata = 0
  nullify(strata_list%first)
  nullify(strata_list%last)
  if (associated(strata_list%array)) deallocate(strata_list%array)
  nullify(strata_list%array)
  
  deallocate(strata_list)
  nullify(strata_list)

end subroutine GeomechStrataDestroyList

! ************************************************************************** !
!
! GeomechStrataDestroy: Destroys a geomechanics strata
! author: Satish Karra, LANL
! date: 06/07/2013
!
! ************************************************************************** !
subroutine GeomechStrataDestroy(strata)

  implicit none
  
  type(geomech_strata_type), pointer :: strata
  
  if (.not.associated(strata)) return
  
  ! since strata%region is a pointer to a region in a list, nullify instead
  ! of destroying since the list will be destroyed separately
  nullify(strata%region)
  
  deallocate(strata)
  nullify(strata)

end subroutine GeomechStrataDestroy

end module Geomechanics_Strata_module

#endif
