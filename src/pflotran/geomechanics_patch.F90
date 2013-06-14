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
            GeomechPatchProcessGeomechCouplers, &
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
! GeomechPatchProcessGeomechCouplers: Assigns conditions and regions to couplers
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine GeomechPatchProcessGeomechCouplers(patch,conditions,option)

  use Option_module
  use Geomechanics_Material_module
  use Geomechanics_Condition_module
  
  implicit none
  
  type(geomech_patch_type)                         :: patch
  type(geomech_condition_list_type)                :: conditions
  type(option_type)                                :: option
  
  type(geomech_coupler_type), pointer              :: coupler
  type(geomech_coupler_list_type), pointer         :: coupler_list 
  type(geomech_strata_type), pointer               :: strata
 ! type(geomech_observation_type), pointer          :: observation, &
 !                                                     next_observation
  
  PetscInt                                         :: temp_int, isub
  
  ! boundary conditions
  coupler => patch%geomech_boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => GeomechRegionGetPtrFromList(coupler%region_name, &
                                                  patch%geomech_regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Geomech Region "' // trim(coupler%region_name) // &
                 '" in Geomech boundary condition "' // &
                 trim(coupler%name) // &
                 '" not found in Geomech region list'
      call printErrMsg(option)
    endif

    ! pointer to flow condition
    if (option%ngeomechdof > 0) then
      if (len_trim(coupler%geomech_condition_name) > 0) then
        coupler%geomech_condition => &
          GeomechConditionGetPtrFromList(coupler%geomech_condition_name, &
                                         conditions)
        if (.not.associated(coupler%geomech_condition)) then
          option%io_buffer = 'Geomech condition "' // &
                   trim(coupler%geomech_condition_name) // &
                   '" in Geomech boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in geomech condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A GEOMECHANICS_CONDITION must be specified in ' // &
                           'GEOMECHANICS_BOUNDARY_CONDITION: ' // &
                            trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo


  ! SK: There are no initial conditions (at this point)

  ! source/sinks
  coupler => patch%geomech_source_sinks%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => GeomechRegionGetPtrFromList(coupler%region_name, &
                                                  patch%geomech_regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Geomech Region "' // trim(coupler%region_name) // &
                 '" in geomech source/sink "' // &
                 trim(coupler%name) // &
                 '" not found in geomech region list'
      call printErrMsg(option)
    endif
    ! pointer to geomech condition
    if (option%ngeomechdof > 0) then    
      if (len_trim(coupler%geomech_condition_name) > 0) then
        coupler%geomech_condition => &
          GeomechConditionGetPtrFromList(coupler%geomech_condition_name, &
                                         conditions)
        if (.not.associated(coupler%geomech_condition)) then
          option%io_buffer = 'Geomech condition "' // &
                   trim(coupler%geomech_condition_name) // &
                   '" in geomech source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in geomech condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A GEOMECHANICS_CONDITION must be specified in ' // &
                           'GEOMECHANICS_SOURCE_SINK: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

!----------------------------  
! AUX  
    
  ! strata
  ! connect pointers from strata to regions
  strata => patch%geomech_strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 1) then
      strata%region => GeomechRegionGetPtrFromList(strata%region_name, &
                                                   patch%geomech_regions)
      if (.not.associated(strata%region)) then
        option%io_buffer = 'Geomech Region "' // trim(strata%region_name) // &
                 '" in geomech strata not found in geomech region list'
        call printErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material_property => &
            GeomechanicsMaterialPropGetPtrFromArray( &
                                        strata%material_property_name, &
                                        patch%geomech_material_property_array)
        if (.not.associated(strata%material_property)) then
          option%io_buffer = 'Geomech Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in geomech material list'
          call printErrMsg(option)
        endif
      endif
    else
      nullify(strata%region)
      nullify(strata%material_property)
    endif
    strata => strata%next
  enddo

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
#if 0
  observation => patch%observation%first
  do
    if (.not.associated(observation)) exit
    next_observation => observation%next
    select case(observation%itype)
      case(OBSERVATION_SCALAR)
        ! pointer to region
        observation%region => RegionGetPtrFromList(observation%linkage_name, &
                                                    patch%regions)
        if (.not.associated(observation%region)) then
          option%io_buffer = 'Region "' // &
                   trim(observation%linkage_name) // &
                 '" in observation point "' // &
                 trim(observation%name) // &
                 '" not found in region list'                   
          call printErrMsg(option)
        endif
        if (observation%region%num_cells == 0) then
          ! remove the observation object
          call ObservationRemoveFromList(observation,patch%observation)
        endif
      case(OBSERVATION_FLUX)
        coupler => CouplerGetPtrFromList(observation%linkage_name, &
                                         patch%boundary_conditions)
        if (associated(coupler)) then
          observation%connection_set => coupler%connection_set
        else
          option%io_buffer = 'Boundary Condition "' // &
                   trim(observation%linkage_name) // &
                   '" not found in Boundary Condition list'
          call printErrMsg(option)
        endif
        if (observation%connection_set%num_connections == 0) then
          ! cannot remove from list, since there must be a global reduction
          ! across all procs
          ! therefore, just nullify connection set
          nullify(observation%connection_set)
        endif                                      
    end select
    observation => next_observation
  enddo
#endif
 
end subroutine GeomechPatchProcessGeomechCouplers

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