#ifdef GEOMECH
module Geomechanics_Patch_module

  use Option_module
  use Geomechanics_Grid_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Region_module
  use Geomechanics_Strata_module
  use Geomechanics_Coupler_module
  use Geomechanics_Field_module
  
  implicit none
  
  private
  
#include "definitions.h"

  type, public :: geomech_patch_type
    PetscInt                                      :: id
    PetscInt, pointer                             :: imat(:)
    type(geomech_grid_type), pointer              :: geomech_grid
    type(geomech_material_property_type), pointer :: geomech_material_properties
    type(geomech_material_property_ptr_type),&
       pointer :: geomech_material_property_array(:)
    type(geomech_strata_list_type), pointer       :: geomech_strata
    type(gm_region_list_type), pointer            :: geomech_regions
    type(geomech_coupler_list_type), pointer      :: geomech_boundary_conditions
    type(geomech_coupler_list_type), pointer      :: geomech_source_sinks
    type(geomech_field_type), pointer             :: geomech_field
  end type geomech_patch_type


  public :: GeomechanicsPatchCreate, &
            GeomechPatchLocalizeRegions, &
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
  
  allocate(patch%geomech_boundary_conditions)
  call GeomechCouplerInitList(patch%geomech_boundary_conditions)
  allocate(patch%geomech_source_sinks)
  call GeomechCouplerInitList(patch%geomech_source_sinks)  
  
  nullify(patch%geomech_material_properties)
  nullify(patch%geomech_material_property_array)
  
  allocate(patch%geomech_strata)
  call GeomechStrataInitList(patch%geomech_strata)

  allocate(patch%geomech_regions)
  call GeomechRegionInitList(patch%geomech_regions)
  
  nullify(patch%geomech_field)
  
  GeomechanicsPatchCreate => patch
  
end function GeomechanicsPatchCreate

! ************************************************************************** !
!
! GeomechPatchLocalizeRegions: Localizes regions within each patch
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine GeomechPatchLocalizeRegions(geomech_patch,regions,option)

  use Option_module
  use Geomechanics_Region_module

  implicit none
  
  type(geomech_patch_type)           :: geomech_patch
  type(gm_region_list_type)          :: regions
  type(option_type)                  :: option
  
  type(gm_region_type), pointer      :: cur_region
  type(gm_region_type), pointer      :: patch_region
  
  cur_region => regions%first
  do
    if (.not.associated(cur_region)) exit
    patch_region => GeomechRegionCreate(cur_region)
    call GeomechRegionAddToList(patch_region,geomech_patch%geomech_regions)
    cur_region => cur_region%next
  enddo
  
 ! Need a call to a subroutine similar to GridlocalizeRegions 
 ! call GridLocalizeRegions(patch%grid,patch%regions,option)
  call GeomechGridLocalizeRegions(geomech_patch%geomech_grid, &
                                  geomech_patch%geomech_regions, &
                                  option)
 
end subroutine GeomechPatchLocalizeRegions


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
  
  call GeomechCouplerDestroyList(geomech_patch%geomech_boundary_conditions)
  call GeomechCouplerDestroyList(geomech_patch%geomech_source_sinks)
  
  nullify(geomech_patch%geomech_field)

  deallocate(geomech_patch)
  nullify(geomech_patch)

end subroutine GeomechanicsPatchDestroy

end module Geomechanics_Patch_module
#endif