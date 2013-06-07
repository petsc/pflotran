#ifdef GEOMECH
module Geomechanics_Patch_module

  use Option_module
  use Geomechanics_Grid_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Region_module
  use Geomechanics_Strata_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: geomech_patch_type
    PetscInt :: id
    PetscInt, pointer :: imat(:)
    type(geomech_grid_type), pointer :: geomech_grid
    type(geomech_material_property_type), pointer :: geomech_material_properties
    type(geomech_material_property_ptr_type),&
       pointer :: geomech_material_property_array(:)
    type(geomech_strata_list_type), pointer :: geomech_strata
    type(gm_region_list_type), pointer :: geomech_regions
  end type geomech_patch_type


  public :: GeomechanicsPatchCreate, &
            GeomechanicsPatchDestroy

contains

! ************************************************************************** !
!
! GeomechanicsPatchCreate: Allocates and initializes a new geomechanics 
! patch  object
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
function GeomechanicsPatchCreate()

  implicit none
  
  type(geomech_patch_type), pointer :: GeomechanicsPatchCreate
  type(geomech_patch_type), pointer :: patch
  
  allocate(patch)
  
  patch%id = 0
  nullify(patch%imat)
  nullify(patch%geomech_grid)
  
  nullify(patch%geomech_material_properties)
  nullify(patch%geomech_material_property_array)
  
  allocate(patch%geomech_strata)
  call GeomechStrataInitList(patch%geomech_strata)

  allocate(patch%geomech_regions)
  call GeomechRegionInitList(patch%geomech_regions)
  
  GeomechanicsPatchCreate => patch
  
end function GeomechanicsPatchCreate

! ************************************************************************** !
!
! GeomechanicsPatchDestroy: Destroys a new geomechanics patch  object
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechanicsPatchDestroy(geomech_patch)

  implicit none
  
  type(geomech_patch_type), pointer :: geomech_patch
  if (associated(geomech_patch%imat)) deallocate(geomech_patch%imat)
  nullify(geomech_patch%imat)
  
  if (associated(geomech_patch%geomech_material_property_array)) &
    deallocate(geomech_patch%geomech_material_property_array)
  nullify(geomech_patch%geomech_material_property_array)
  if (associated(geomech_patch%geomech_material_properties)) &
    deallocate(geomech_patch%geomech_material_properties)
  nullify(geomech_patch%geomech_material_properties)

  call GeomechStrataDestroyList(geomech_patch%geomech_strata)
  call GeomechRegionDestroyList(geomech_patch%geomech_regions)

  deallocate(geomech_patch)
  nullify(geomech_patch)

end subroutine GeomechanicsPatchDestroy

end module Geomechanics_Patch_module
#endif