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
  use Geomechanics_Auxiliary_module
  use Dataset_Base_class
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
#include "finclude/petscsys.h"

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
    class(dataset_base_type), pointer             :: geomech_datasets
    type(geomech_auxiliary_type)                  :: geomech_aux
  end type geomech_patch_type


  public :: GeomechanicsPatchCreate, &
            GeomechPatchLocalizeRegions, &
            GeomechPatchProcessGeomechCouplers, &
            GeomechPatchInitAllCouplerAuxVars, &
            GeomechPatchGetDataset, &
            GeomechanicsPatchDestroy

contains

! ************************************************************************** !
!
! GeomechanicsPatchCreate: Allocates and initializes a new geomechanics 
!                          patch object
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
  
  call GeomechAuxInit(patch%geomech_aux)
  
  nullify(patch%geomech_field)
  nullify(patch%geomech_datasets)
  
  GeomechanicsPatchCreate => patch
  
end function GeomechanicsPatchCreate

! ************************************************************************** !
!
! GeomechPatchLocalizeRegions: Localizes regions within each patch
! author: Satish Karra, LANL
! date: 05/23/13
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
! author: Satish Karra, LANL
! date: 05/23/13
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

    ! pointer to geomech condition
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
! GeomechPatchInitAllCouplerAuxVars: Initializes coupler auxillary variables 
!                                    within list
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechPatchInitAllCouplerAuxVars(patch,option)

  use Option_module
  
  implicit none
  
  type(geomech_patch_type), pointer :: patch
  type(option_type) :: option
  
  PetscBool :: force_update_flag = PETSC_TRUE
  
  call GeomechPatchInitCouplerAuxVars(patch%geomech_boundary_conditions,patch, &
                                      option)
  call GeomechPatchInitCouplerAuxVars(patch%geomech_source_sinks,patch, &
                                      option)

  call GeomechPatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

end subroutine GeomechPatchInitAllCouplerAuxVars

! ************************************************************************** !
!
! GeomechPatchInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                 within list
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechPatchInitCouplerAuxVars(coupler_list,patch,option)

  use Option_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Condition_module
  
  implicit none
  
  type(geomech_coupler_list_type), pointer :: coupler_list
  type(geomech_patch_type), pointer :: patch
  type(option_type) :: option
  
  PetscInt :: num_verts
  PetscBool :: force_update_flag
  
  type(geomech_coupler_type), pointer :: coupler
  PetscInt :: idof
  character(len=MAXSTRINGLENGTH) :: string
  
  if (.not.associated(coupler_list)) return
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%region)) then
      num_verts = coupler%region%num_verts
      if (associated(coupler%geomech_condition) .and. &
          coupler%itype == GM_BOUNDARY_COUPLER_TYPE) then
        if (associated(coupler%geomech_condition%displacement_x) .or. &
            associated(coupler%geomech_condition%displacement_y) .or. &
            associated(coupler%geomech_condition%displacement_z) .or. &
            associated(coupler%geomech_condition%force_x) .or. &
            associated(coupler%geomech_condition%force_y) .or. &
            associated(coupler%geomech_condition%force_z)) then
          ! allocate arrays that match the number of boundary vertices
          allocate(coupler%geomech_aux_real_var(option%ngeomechdof,num_verts))
          allocate(coupler%geomech_aux_int_var(1,num_verts))
          coupler%geomech_aux_real_var = 0.d0
          coupler%geomech_aux_int_var = 0
        endif ! associated(coupler%geomech_condition%displacement_x)
      else if (coupler%itype == GM_SRC_SINK_COUPLER_TYPE) then
        option%io_buffer='Source/Sink not implemented for geomechanics.'
        call printErrMsg(option)
      endif ! coupler%itype == GM_SRC_SINK_COUPLER_TYPE
    endif ! associated(coupler%region)
      
    coupler => coupler%next
    
  enddo
  
end subroutine GeomechPatchInitCouplerAuxVars

! ************************************************************************** !
!
! GeomechPatchUpdateAllCouplerAuxVars: Updates auxiliary variables associated 
!                                      with couplers in list
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechPatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

  use Option_module
  
  implicit none
  
  type(geomech_patch_type) :: patch
  PetscBool :: force_update_flag
  type(option_type) :: option

  !geh: no need to update initial conditions as they only need updating
  !     once as performed in PatchInitCouplerAuxVars()
  call GeomechPatchUpdateCouplerAuxVars(patch, &
                                        patch%geomech_boundary_conditions, &
                                        force_update_flag,option)
  call GeomechPatchUpdateCouplerAuxVars(patch,patch%geomech_source_sinks, &
                                        force_update_flag,option)

end subroutine GeomechPatchUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! GeomechPatchUpdateCouplerAuxVars: Updates auxiliary variables associated 
!                                   with couplers in list
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechPatchUpdateCouplerAuxVars(patch,coupler_list, &
                                            force_update_flag,option)
                                     
  use Option_module
  use Geomechanics_Condition_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module

  implicit none
  
  type(geomech_patch_type) :: patch
  type(geomech_coupler_list_type), pointer :: coupler_list
  PetscBool :: force_update_flag
  type(option_type) :: option
  
  type(geomech_coupler_type), pointer :: coupler
  type(geomech_condition_type), pointer :: geomech_condition
  PetscBool :: update
  character(len=MAXSTRINGLENGTH) :: string,string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof,num_verts
  PetscInt :: ivertex,local_id,ghosted_id
  

  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first
    
  do
    if (.not.associated(coupler)) exit    
    if (associated(coupler%geomech_aux_real_var)) then
      num_verts = coupler%region%num_verts
      geomech_condition => coupler%geomech_condition
      if (force_update_flag .or. &
          GeomechConditionIsTransient(geomech_condition)) then
        if (associated(geomech_condition%displacement_x)) then
          select case(geomech_condition%displacement_x%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_X_DOF, &
                                           1:num_verts) = &
              geomech_condition%displacement_x%geomech_dataset% &
                time_series%cur_value(1)
          end select
        endif
        if (associated(geomech_condition%displacement_y)) then
          select case(geomech_condition%displacement_y%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_Y_DOF, &
                                           1:num_verts) = &
              geomech_condition%displacement_y%geomech_dataset% &
                time_series%cur_value(1)
          end select
        endif
        if (associated(geomech_condition%displacement_z)) then
          select case(geomech_condition%displacement_z%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                                           1:num_verts) = &
              geomech_condition%displacement_z%geomech_dataset% &
                time_series%cur_value(1)
           end select
        endif
        if (associated(geomech_condition%force_x)) then
          select case(geomech_condition%force_x%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                                           1:num_verts) = &
              geomech_condition%force_x%geomech_dataset% &
                time_series%cur_value(1)
            end select
        endif
        if (associated(geomech_condition%force_y)) then
          select case(geomech_condition%force_y%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                                           1:num_verts) = &
              geomech_condition%force_y%geomech_dataset% &
                time_series%cur_value(1)
             end select
        endif
        if (associated(geomech_condition%force_z)) then
          select case(geomech_condition%force_z%itype)
            case(DIRICHLET_BC)
              coupler%geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                                           1:num_verts) = &
              geomech_condition%force_z%geomech_dataset% &
                time_series%cur_value(1)
          end select
        endif        
      endif
    endif

    coupler => coupler%next
  enddo

end subroutine GeomechPatchUpdateCouplerAuxVars

! ************************************************************************** !
!
! GeomechPatchGetDataset: Extracts variables indexed by ivar and isubvar 
!                         from a geomechanics patch
! author: Satish Karra, LANL
! date: 07/02/13
!
! ************************************************************************** !
subroutine GeomechPatchGetDataset(patch,geomech_field,option,output_option, &
                                  vec,ivar,isubvar,isubvar1)

  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Output_Aux_module
  use Geomechanics_Field_module
  use Variables_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  !type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: iphase

  PetscInt :: local_id
  type(geomech_grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  grid => patch%geomech_grid

  call GeomechGridVecGetArrayF90(grid,vec,vec_ptr,ierr)
  
  iphase = 1
  
  select case(ivar)
    case(GEOMECH_DISP_X)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%disp_vector(1)
      enddo
    case(GEOMECH_DISP_Y)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%disp_vector(2)
      enddo
    case(GEOMECH_DISP_Z)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%disp_vector(3)
      enddo
    case(STRAIN_XX)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(1)
      enddo
    case(STRAIN_YY)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(2)
      enddo
    case(STRAIN_ZZ)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(3)
      enddo
    case(STRAIN_XY)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(4)
      enddo
    case(STRAIN_YZ)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(5)
      enddo
    case(STRAIN_ZX)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%strain(6)
      enddo
    case(STRESS_XX)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(1)
      enddo
    case(STRESS_YY)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(2)
      enddo
    case(STRESS_ZZ)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(3)
      enddo
    case(STRESS_XY)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(4)
      enddo
    case(STRESS_YZ)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(5)
      enddo
    case(STRESS_ZX)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%geomech_aux%GeomechGlobal%aux_vars(grid%nL2G(local_id))%stress(6)
      enddo
    case(GEOMECH_MATERIAL_ID)
      do local_id=1,grid%nlmax_node
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in GeomechPatchGetDataset'')') ivar
      call printErrMsg(option)
  end select

end subroutine GeomechPatchGetDataset

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
  nullify(geomech_patch%geomech_datasets)
  
  call GeomechAuxDestroy(geomech_patch%geomech_aux)

  deallocate(geomech_patch)
  nullify(geomech_patch)

end subroutine GeomechanicsPatchDestroy

end module Geomechanics_Patch_module
#endif
