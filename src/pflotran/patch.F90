module Patch_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Observation_module
  use Strata_module
  use Region_module
  use Reaction_Aux_module
  use Dataset_Base_class
  use Material_module
  use Field_module
  use Saturation_Function_module
#ifdef SURFACE_FLOW
  use Surface_Field_module
  use Surface_Material_module
  use Surface_Auxiliary_module
#endif
  
  use Auxiliary_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: patch_type 
    
    PetscInt :: id
    
    ! These arrays will be used by all modes, mode-specific arrays should
    ! go in the auxiliary data stucture for that mode
    PetscInt, pointer :: imat(:)
    PetscInt, pointer :: sat_func_id(:)

    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    PetscReal, pointer :: internal_fluxes(:,:,:)    
    PetscReal, pointer :: boundary_fluxes(:,:,:)  
    PetscReal, pointer :: internal_tran_coefs(:,:)
    PetscReal, pointer :: boundary_tran_coefs(:,:)
    PetscReal, pointer :: ss_fluid_fluxes(:,:)

    type(grid_type), pointer :: grid

    type(region_list_type), pointer :: regions

    type(coupler_list_type), pointer :: boundary_conditions
    type(coupler_list_type), pointer :: initial_conditions
    type(coupler_list_type), pointer :: source_sinks

    type(material_property_type), pointer :: material_properties
    type(material_property_ptr_type), pointer :: material_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)

    type(strata_list_type), pointer :: strata
    type(observation_list_type), pointer :: observation

    ! Pointers to objects in mother realization object
    type(field_type), pointer :: field 
    type(reaction_type), pointer :: reaction
    class(dataset_base_type), pointer :: datasets
    
    type(auxiliary_type) :: aux
    
    type(patch_type), pointer :: next

    PetscInt :: surf_or_subsurf_flag  ! Flag to identify if the current patch
                                      ! is a surface or subsurface (default)
#ifdef SURFACE_FLOW
    type(surface_material_property_type), pointer     :: surf_material_properties
    type(surface_material_property_ptr_type), pointer :: surf_material_property_array(:)
    type(surface_field_type),pointer                  :: surf_field
    type(surface_auxiliary_type) :: surf_aux
    
    PetscReal,pointer :: surf_internal_fluxes(:,:)
    PetscReal,pointer :: surf_boundary_fluxes(:,:)
#endif

  end type patch_type

  ! pointer data structure required for making an array of patch pointers in F90
  type, public :: patch_ptr_type
    type(patch_type), pointer :: ptr           ! pointer to the patch_type
  end type patch_ptr_type 

  type, public :: patch_list_type
    PetscInt :: num_patch_objects
    type(patch_type), pointer :: first
    type(patch_type), pointer :: last
    type(patch_ptr_type), pointer :: array(:)
  end type patch_list_type

  PetscInt, parameter, public :: INT_VAR = 0
  PetscInt, parameter, public :: REAL_VAR = 1
    
  interface PatchGetVariable
    module procedure PatchGetVariable1
#ifdef SURFACE_FLOW
    module procedure PatchGetVariable2
#endif
  end interface

  public :: PatchCreate, PatchDestroy, PatchCreateList, PatchDestroyList, &
            PatchAddToList, PatchConvertListToArray, PatchProcessCouplers, &
            PatchUpdateAllCouplerAuxVars, PatchInitAllCouplerAuxVars, &
            PatchLocalizeRegions, PatchUpdateUniformVelocity, &
            PatchGetVariable, PatchGetVariableValueAtCell, &
            PatchSetVariable, &
            PatchInitConstraints, &
            PatchCountCells, PatchGetIvarsFromKeyword, &
            PatchGetVarNameFromKeyword, &
            PatchCalculateCFL1Timestep

contains

! ************************************************************************** !
!
! PatchCreate: Allocates and initializes a new Patch object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function PatchCreate()

  implicit none
  
  type(patch_type), pointer :: PatchCreate
  
  type(patch_type), pointer :: patch
  
  allocate(patch)

  patch%id = 0
  patch%surf_or_subsurf_flag = SUBSURFACE
  nullify(patch%imat)
  nullify(patch%sat_func_id)
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)
  nullify(patch%internal_fluxes)
  nullify(patch%boundary_fluxes)
  nullify(patch%internal_tran_coefs)
  nullify(patch%boundary_tran_coefs)
  nullify(patch%ss_fluid_fluxes)

  nullify(patch%grid)

  allocate(patch%regions)
  call RegionInitList(patch%regions)
  
  allocate(patch%boundary_conditions)
  call CouplerInitList(patch%boundary_conditions)
  allocate(patch%initial_conditions)
  call CouplerInitList(patch%initial_conditions)
  allocate(patch%source_sinks)
  call CouplerInitList(patch%source_sinks)

  nullify(patch%material_properties)
  nullify(patch%material_property_array)
  nullify(patch%saturation_functions)
  nullify(patch%saturation_function_array)

  allocate(patch%observation)
  call ObservationInitList(patch%observation)

  allocate(patch%strata)
  call StrataInitList(patch%strata)
  
  call AuxInit(patch%aux)
  
  nullify(patch%field)
  nullify(patch%reaction)
  nullify(patch%datasets)
  
  nullify(patch%next)
  
#ifdef SURFACE_FLOW
    nullify(patch%surf_material_properties)
    nullify(patch%surf_material_property_array)
    nullify(patch%surf_field)
    nullify(patch%surf_internal_fluxes)
    nullify(patch%surf_boundary_fluxes)
    call SurfaceAuxInit(patch%surf_aux)
#endif

  PatchCreate => patch
  
end function PatchCreate

! ************************************************************************** !
!
! PatchListCreate: Creates a patch list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
function PatchCreateList()

  implicit none

  type(patch_list_type), pointer :: PatchCreateList

  type(patch_list_type), pointer :: patch_list
  
  allocate(patch_list)
  nullify(patch_list%first)
  nullify(patch_list%last)
  nullify(patch_list%array)
  patch_list%num_patch_objects = 0

  PatchCreateList => patch_list

end function PatchCreateList

! ************************************************************************** !
!
! PatchAddToList: Adds a new patch to list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchAddToList(new_patch,patch_list)

  implicit none
  
  type(patch_type), pointer :: new_patch
  type(patch_list_type) :: patch_list
  
  if(associated(new_patch)) then
     patch_list%num_patch_objects = patch_list%num_patch_objects + 1
     new_patch%id = patch_list%num_patch_objects
     if (.not.associated(patch_list%first)) patch_list%first => new_patch
     if (associated(patch_list%last)) patch_list%last%next => new_patch
     patch_list%last => new_patch
  end if
end subroutine PatchAddToList

! ************************************************************************** !
!
! PatchConvertListToArray: Creates an array of pointers to the 
!                               patchs in the patch list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchConvertListToArray(patch_list)

  implicit none
  
  type(patch_list_type) :: patch_list
    
  PetscInt :: count
  type(patch_type), pointer :: cur_patch
  
  
  allocate(patch_list%array(patch_list%num_patch_objects))
  
  cur_patch => patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    patch_list%array(cur_patch%id)%ptr => cur_patch
    cur_patch => cur_patch%next
  enddo

end subroutine PatchConvertListToArray

! ************************************************************************** !
!
! PatchLocalizeRegions: Localizes regions within each patch
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchLocalizeRegions(patch,regions,option)

  use Option_module
  use Region_module

  implicit none
  
  type(patch_type) :: patch
  type(region_list_type) :: regions
  type(option_type) :: option
  
  type(region_type), pointer :: cur_region
  type(region_type), pointer :: patch_region
  
  cur_region => regions%first
  do
    if (.not.associated(cur_region)) exit
    patch_region => RegionCreate(cur_region)
    call RegionAddToList(patch_region,patch%regions)
    cur_region => cur_region%next
  enddo
  
  !geh: All grids must be localized through GridLocalizeRegions.  Patch
  !     should not differentiate between structured/unstructured, etc.
  call GridLocalizeRegions(patch%grid,patch%regions,option)
 
end subroutine PatchLocalizeRegions

! ************************************************************************** !
!
! PatchProcessCouplers: Assigns conditions and regions to couplers
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
                                option)

  use Option_module
  use Material_module
  use Condition_module
  use Constraint_module
  use Connection_module

  implicit none
  
  type(patch_type) :: patch
  type(condition_list_type) :: flow_conditions
  type(tran_condition_list_type) :: transport_conditions
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list 
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation, next_observation
  
  PetscInt :: temp_int, isub
  
  ! boundary conditions
  coupler => patch%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in boundary condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
    if (associated(patch%grid%structured_grid)) then
      if (coupler%region%num_cells > 0 .and. &
          (coupler%region%iface == 0 .and. &
           .not.associated(coupler%region%faces))) then
        option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '", which is tied to a boundary condition, has not ' // &
                 'been assigned a face in the structured grid. '
        call printErrMsg(option)
      endif
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
           option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in ' // &
                           'BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => patch%initial_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in initial condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name,transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in ' // &
                           'INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => patch%source_sinks%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in source/sink "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then    
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name,flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'SOURCE_SINK: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then    
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in ' // &
                           'SOURCE_SINK: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

!----------------------------  
! AUX  
    
  ! strata
  ! connect pointers from strata to regions
  strata => patch%strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 1) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  patch%regions)
      if (.not.associated(strata%region)) then
        option%io_buffer = 'Region "' // trim(strata%region_name) // &
                 '" in strata not found in region list'
        call printErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        ! gb: Depending on a surface/subsurface patch, use corresponding
        !     material properties
        if (patch%surf_or_subsurf_flag == SUBSURFACE) then
          strata%material_property => &
            MaterialPropGetPtrFromArray(strata%material_property_name, &
                                        patch%material_property_array)
          if (.not.associated(strata%material_property)) then
            option%io_buffer = 'Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in material list'
            call printErrMsg(option)
          endif
        endif

#ifdef SURFACE_FLOW
        if(patch%surf_or_subsurf_flag == SURFACE) then
          strata%surf_material_property => &
            SurfaceMaterialPropGetPtrFromArray(strata%material_property_name, &
                                            patch%surf_material_property_array)
          if (.not.associated(strata%surf_material_property)) then
            option%io_buffer = 'Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in material list'
            call printErrMsg(option)
          endif
        endif
#endif

      endif
    else
      nullify(strata%region)
      nullify(strata%material_property)
    endif
    strata => strata%next
  enddo

  ! connectivity between initial conditions, boundary conditions, srcs/sinks, etc and grid
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%initial_conditions)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%boundary_conditions)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%source_sinks)

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
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
 
  temp_int = ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
  allocate(patch%internal_velocities(option%nphase,temp_int))
  patch%internal_velocities = 0.d0
  allocate(patch%internal_tran_coefs(option%nphase,temp_int))
  patch%internal_tran_coefs = 0.d0
  if (option%store_solute_fluxes) then
    allocate(patch%internal_fluxes(option%nphase,option%ntrandof,temp_int))
    patch%internal_fluxes = 0.d0
  endif
  if (option%store_flowrate) then
    if(option%store_solute_fluxes) then
      option%io_buffer='Model does not support store_solute_fluxes and flowrate ' // &
      ' options together. If you run into this message, complain on pflotran-dev@googlegroups.com'
      call printErrMsg(option)
    endif
    allocate(patch%internal_fluxes(option%nflowdof,1,temp_int))
    allocate(patch%boundary_fluxes(option%nflowdof,1,temp_int))
    patch%internal_fluxes = 0.d0
    patch%boundary_fluxes = 0.d0
  endif
#ifdef SURFACE_FLOW
  if (patch%surf_or_subsurf_flag == SURFACE) then
    !if (option%store_flowrate) then
      allocate(patch%surf_internal_fluxes(option%nflowdof,temp_int))
      patch%surf_internal_fluxes = 0.d0
    !endif
  endif
  ! Always allocate the array to store boundary fluxes as they are needed
  ! to store data for hydrograph output
  allocate(patch%surf_boundary_fluxes(option%nflowdof,temp_int))
  patch%surf_boundary_fluxes = 0.d0
#endif
 
  if (patch%grid%itype == STRUCTURED_GRID_MIMETIC.or. &
      patch%grid%discretization_itype == UNSTRUCTURED_GRID_MIMETIC ) then
    temp_int = CouplerGetNumBoundConnectionsInListMFD(patch%grid, &
                                                 patch%boundary_conditions, &
                                                 option)
  else  
    temp_int = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  end if

  if (temp_int > 0) then
    allocate(patch%boundary_velocities(option%nphase,temp_int)) 
    patch%boundary_velocities = 0.d0
    allocate(patch%boundary_tran_coefs(option%nphase,temp_int))
    patch%boundary_tran_coefs = 0.d0
    if (option%store_solute_fluxes) then
      allocate(patch%boundary_fluxes(option%nphase,option%ntrandof,temp_int))
      patch%boundary_fluxes = 0.d0
    endif
  endif

  temp_int = CouplerGetNumConnectionsInList(patch%source_sinks)
  if (temp_int > 0) then
    allocate(patch%ss_fluid_fluxes(option%nphase,temp_int))
    patch%ss_fluid_fluxes = 0.d0
  endif

end subroutine PatchProcessCouplers

! ************************************************************************** !
!
! PatchInitAllCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchInitAllCouplerAuxVars(patch,option)

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  type(patch_type), pointer :: patch
  type(option_type) :: option
  
  PetscBool :: force_update_flag = PETSC_TRUE
  
  call PatchInitCouplerAuxVars(patch%initial_conditions,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%boundary_conditions,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%source_sinks,patch, &
                               option)

  !geh: This should not be included in PatchUpdateAllCouplerAuxVars
  ! as it will result in excessive updates to initial conditions
  ! that are not necessary after the simulation has started time stepping.
  call PatchUpdateCouplerAuxVars(patch,patch%initial_conditions, &
                                 force_update_flag,option)
  call PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

end subroutine PatchInitAllCouplerAuxVars

! ************************************************************************** !
!
! PatchInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchInitCouplerAuxVars(coupler_list,patch,option)

  use Option_module
  use Connection_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Condition_module
  use Constraint_module
  
  implicit none
  
  type(coupler_list_type), pointer :: coupler_list
  type(patch_type), pointer :: patch
  type(option_type) :: option
  
  PetscInt :: num_connections
  PetscBool :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  PetscInt :: idof
  character(len=MAXSTRINGLENGTH) :: string
  
  if (.not.associated(coupler_list)) return
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%connection_set)) then
      num_connections = coupler%connection_set%num_connections
      
      ! FLOW
      if (associated(coupler%flow_condition) .and. &
          (coupler%itype == INITIAL_COUPLER_TYPE .or. &
           coupler%itype == BOUNDARY_COUPLER_TYPE)) then

        if (associated(coupler%flow_condition%pressure) .or. &
            associated(coupler%flow_condition%concentration) .or. &
            associated(coupler%flow_condition%saturation) .or. &
            associated(coupler%flow_condition%rate) .or. &
            associated(coupler%flow_condition%temperature) .or. &
            associated(coupler%flow_condition%general)) then

          ! allocate arrays that match the number of connections
          select case(option%iflowmode)

            case(RICHARDS_MODE)
!geh              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              if (option%mimetic) then
                 if (coupler%itype == INITIAL_COUPLER_TYPE) then 
                     num_connections = coupler%numfaces_set + coupler%region%num_cells
                 else 
                     num_connections = coupler%numfaces_set
                 end if
              end if
              allocate(coupler%flow_aux_real_var(2,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0

            case(TH_MODE)
              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0

            case(THC_MODE)
              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
              
            case(MPH_MODE, IMS_MODE, FLASH2_MODE, MIS_MODE)
!geh              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_real_var(option%nflowdof,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
                
            case(G_MODE)
              allocate(coupler%flow_aux_real_var(FOUR_INTEGER,num_connections))
              allocate(coupler%flow_aux_int_var(ONE_INTEGER,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
                
            case default
          end select
      
        endif ! associated(coupler%flow_condition%pressure)
      
      else if (coupler%itype == SRC_SINK_COUPLER_TYPE) then

        if (associated(coupler%flow_condition%rate)) then

          select case(coupler%flow_condition%rate%itype)
            case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS, &
                 HET_VOL_RATE_SS,HET_MASS_RATE_SS)
              select case(option%iflowmode)
                case(RICHARDS_MODE)
                  allocate(coupler%flow_aux_real_var(1,num_connections))
                  coupler%flow_aux_real_var = 0.d0
                case(TH_MODE)
                  allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
                  coupler%flow_aux_real_var = 0.d0
                case default
                  string = GetSubConditionName(coupler%flow_condition%rate%itype)
                  option%io_buffer='Source/Sink of rate%itype = "' // &
                    trim(adjustl(string)) // '", not implemented in this mode.'
                  call printErrMsg(option)
              end select
          end select
        endif ! associated(coupler%flow_condition%rate)
      endif ! coupler%itype == SRC_SINK_COUPLER_TYPE
    endif ! associated(coupler%connection_set)

    ! TRANSPORT   
    if (associated(coupler%tran_condition)) then
      cur_constraint_coupler => &
        coupler%tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        ! Setting option%iflag = 0 ensures that the "mass_balance" array
        ! is not allocated.
        option%iflag = 0
        ! Only allocate the XXX_auxvar objects if they have not been allocated.
        ! Since coupler%tran_condition is a pointer to a separate list of
        ! tran conditions, the XXX_auxvar object may already be allocated.
        if (.not.associated(cur_constraint_coupler%global_auxvar)) then
          allocate(cur_constraint_coupler%global_auxvar)
          call GlobalAuxVarInit(cur_constraint_coupler%global_auxvar,option)
        endif
        if (.not.associated(cur_constraint_coupler%rt_auxvar)) then
          allocate(cur_constraint_coupler%rt_auxvar)
          call RTAuxVarInit(cur_constraint_coupler%rt_auxvar,patch%reaction, &
                            option)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
      
    coupler => coupler%next
  enddo
  
end subroutine PatchInitCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateAllCouplerAuxVars: Updates auxiliary variables associated 
!                                  with couplers in list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

  use Option_module
  
  implicit none
  
  type(patch_type) :: patch
  PetscBool :: force_update_flag
  type(option_type) :: option

  PetscInt :: iconn
  
  !geh: no need to update initial conditions as they only need updating
  !     once as performed in PatchInitCouplerAuxVars()
  call PatchUpdateCouplerAuxVars(patch,patch%boundary_conditions, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%source_sinks, &
                                 force_update_flag,option)

!  stop
end subroutine PatchUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVars: Updates auxiliary variables associated 
!                                  with couplers in list
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVars(patch,coupler_list,force_update_flag, &
                                     option)
  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  PetscBool :: force_update_flag
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: flow_condition

  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      flow_condition => coupler%flow_condition
      if (force_update_flag .or. FlowConditionIsTransient(flow_condition)) then
        select case(option%iflowmode)
          case(G_MODE)
            call PatchUpdateCouplerAuxVarsG(patch,coupler,option)
          case(MPH_MODE)
            call PatchUpdateCouplerAuxVarsMPH(patch,coupler,option)
          case(IMS_MODE)
            call PatchUpdateCouplerAuxVarsIMS(patch,coupler,option)
          case(FLASH2_MODE)
            call PatchUpdateCouplerAuxVarsFLASH2(patch,coupler,option)
          case(THC_MODE)
            call PatchUpdateCouplerAuxVarsTHC(patch,coupler,option)
          case(TH_MODE)
            call PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
          case(MIS_MODE)
            call PatchUpdateCouplerAuxVarsMIS(patch,coupler,option)
          case(RICHARDS_MODE)
            call PatchUpdateCouplerAuxVarsRich(patch,coupler,option)
        end select
      endif
    endif
      
    ! TRANSPORT
    ! nothing for transport at this point in time
    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for G_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsG(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  
  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  general => flow_condition%general
  dof1 = PETSC_FALSE
  dof2 = PETSC_FALSE
  dof3 = PETSC_FALSE
  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
    flow_condition%iphase            
  select case(flow_condition%iphase)
    case(TWO_PHASE_STATE)
      select case(general%gas_pressure%itype)
        case(DIRICHLET_BC)
          coupler%flow_aux_real_var(GENERAL_GAS_PRESSURE_DOF,1:num_connections) = &
            general%gas_pressure%dataset%rarray(1)
          dof1 = PETSC_TRUE
        case default
      end select
      select case(general%gas_saturation%itype)
        case(DIRICHLET_BC)
          coupler%flow_aux_real_var(GENERAL_GAS_SATURATION_DOF,1:num_connections) = &
            general%gas_saturation%dataset%rarray(1)
          dof2 = PETSC_TRUE
      end select
      select case(general%temperature%itype)
        case(DIRICHLET_BC)
          temperature = general%temperature%dataset%rarray(1)
          call psat(temperature,p_sat,ierr)
          coupler%flow_aux_real_var(GENERAL_AIR_PRESSURE_DOF,1:num_connections) = &
            general%gas_pressure%dataset%rarray(1) - p_sat
          dof3 = PETSC_TRUE
      end select
    case(LIQUID_STATE)
      if (general%liquid_pressure%itype == HYDROSTATIC_BC) then
        if (general%mole_fraction%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure bc for flow condition "' // &
            trim(flow_condition%name) // '" requires a mole fraction bc of type dirichlet'
          call printErrMsg(option)
        endif
        if (general%temperature%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure bc for flow condition "' // &
            trim(flow_condition%name) // '" requires a temperature bc of type dirichlet'
          call printErrMsg(option)
        endif
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
#if 0                  
        do iconn=1,coupler%connection_set%num_connections
          local_id = coupler%connection_set%id_dn(iconn)
          ghosted_id = patch%grid%nL2G(local_id)
          x(1:option%nflowdof) = &
            coupler%flow_aux_real_var(1:option%nflowdof,iconn)
          patch%aux%Global%aux_vars(ghosted_id)%istate = LIQUID_STATE
          call GeneralAuxVarUpdateState(x, &
                 patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id), &
                 patch%aux%Global%aux_vars(ghosted_id), &
                 patch%saturation_function_array( &
                   patch%sat_func_id(ghosted_id))%ptr, &
                 0.d0,0.d0,ghosted_id,option)
          coupler%flow_aux_real_var(1:option%nflowdof,iconn) = &
            x(1:option%nflowdof)
          coupler%flow_aux_int_var(ONE_INTEGER,iconn) = &
            patch%aux%Global%aux_vars(ghosted_id)%istate  
          patch%aux%Global%aux_vars(ghosted_id)%istate = &
            patch%aux%Global%aux_vars(ghosted_id)%istate    
        enddo
#endif                  
      else
        select case(general%liquid_pressure%itype)
          case(DIRICHLET_BC)
            coupler%flow_aux_real_var(GENERAL_LIQUID_PRESSURE_DOF,1:num_connections) = &
              general%liquid_pressure%dataset%rarray(1)
            dof1 = PETSC_TRUE
        end select
        select case(general%mole_fraction%itype)
          case(DIRICHLET_BC)
            coupler%flow_aux_real_var(GENERAL_LIQUID_STATE_MOLE_FRACTION_DOF,1:num_connections) = &
              general%mole_fraction%dataset%rarray(1)
            dof2 = PETSC_TRUE
        end select
        select case(general%temperature%itype)
          case(DIRICHLET_BC)
            coupler%flow_aux_real_var(GENERAL_LIQUID_STATE_TEMPERATURE_DOF,1:num_connections) = &
              general%temperature%dataset%rarray(1)
            dof3 = PETSC_TRUE
        end select
      endif
    case(GAS_STATE)
      select case(general%gas_pressure%itype)
        case(DIRICHLET_BC)
          coupler%flow_aux_real_var(GENERAL_GAS_PRESSURE_DOF,1:num_connections) = &
            general%gas_pressure%dataset%rarray(1)
          dof1 = PETSC_TRUE
      end select
      select case(general%mole_fraction%itype)
        case(DIRICHLET_BC)
          coupler%flow_aux_real_var(GENERAL_AIR_PRESSURE_DOF,1:num_connections) = &
            general%mole_fraction%dataset%rarray(1) * &
            general%gas_pressure%dataset%rarray(1)
          dof2 = PETSC_TRUE
      end select                
      select case(general%temperature%itype)
        case(DIRICHLET_BC)
          coupler%flow_aux_real_var(GENERAL_GAS_STATE_TEMPERATURE_DOF,1:num_connections) = &
            general%temperature%dataset%rarray(1)
          dof3 = PETSC_TRUE
      end select
    case(ANY_STATE)
      do iconn = 1, num_connections
        coupler%flow_aux_real_var(1:3,iconn) = &
            general%temperature%dataset%rarray(1:3)
      enddo
      dof1 = PETSC_TRUE
      dof2 = PETSC_TRUE
      dof3 = PETSC_TRUE
  end select  
  !geh: is this really correct, or should it be .or.
  if (.not.dof1 .or. .not.dof2 .or. .not.dof3) then
    option%io_buffer = 'Error with general phase boundary condition'
  endif

end subroutine PatchUpdateCouplerAuxVarsG

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for MPH_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsMPH(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
 !  case(SATURATION_BC)
    end select
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%temperature%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
        endif
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%concentration%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
        endif
    end select
  else
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
         coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsMPH

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for IMS_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsIMS(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  
  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
 !  case(SATURATION_BC)
    end select
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%temperature%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
        endif
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%concentration%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
        endif
    end select
  else
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
         coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsIMS

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for FLASH2_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsFLASH2(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
 !  case(SATURATION_BC)
    end select
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%temperature%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
        endif
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%concentration%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
        endif
    end select
  else
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
         coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsFLASH2

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for THC_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsTHC(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
 !  case(SATURATION_BC)
    end select
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%temperature%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
        endif
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%concentration%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
        endif
    end select
  else
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
         coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsTHC

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for TH_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsTH(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  
  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(TH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
      case(HET_DIRICHLET)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%temperature%dataset, &
                num_connections,TH_PRESSURE_DOF,option)
      case(HET_SURF_SEEPAGE_BC)
        ! Do nothing, since this BC type is only used for coupling of
        ! surface-subsurface model
      case default
        string = GetSubConditionName(flow_condition%pressure%itype)
        option%io_buffer='For TH mode: flow_condition%pressure%itype = "' // &
          trim(adjustl(string)) // '", not implemented.'
          write(*,*), trim(string)
        call printErrMsg(option)
    end select
    if(associated(flow_condition%temperature)) then
      select case(flow_condition%temperature%itype)
        case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
          if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
             (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
             flow_condition%temperature%itype /= DIRICHLET_BC)) then
            coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
                    flow_condition%temperature%dataset%rarray(1)
          endif
        case (HET_DIRICHLET)
          call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                  flow_condition%temperature%dataset, &
                  num_connections,TH_TEMPERATURE_DOF,option)
        case default
          string = GetSubConditionName(flow_condition%temperature%itype)
          option%io_buffer='For TH mode: flow_condition%temperature%itype = "' // &
            trim(adjustl(string)) // '", not implemented.'
          call printErrMsg(option)
      end select
    endif
  else
    if(associated(flow_condition%temperature)) then
      select case(flow_condition%temperature%itype)
        case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
          coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
                    flow_condition%temperature%dataset%rarray(1)
        case (HET_DIRICHLET)
          call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                  flow_condition%temperature%dataset, &
                  num_connections,TH_TEMPERATURE_DOF,option)
        case default
          write(string,*),flow_condition%temperature%itype
          string = GetSubConditionName(flow_condition%temperature%itype)
          option%io_buffer='For TH mode: flow_condition%temperature%itype = "' // &
            trim(adjustl(string)) // '", not implemented.'
          call printErrMsg(option)
      end select
    endif
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case (HET_MASS_RATE_SS,HET_VOL_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                  flow_condition%rate%dataset, &
                  num_connections,TH_PRESSURE_DOF,option)
      case (MASS_RATE_SS)
          coupler%flow_aux_real_var(TH_PRESSURE_DOF,1:num_connections) = &
                  flow_condition%rate%dataset%rarray(1)
      case default
        write(string,*),flow_condition%rate%itype
        string = GetSubConditionName(flow_condition%rate%itype)
        option%io_buffer='For TH mode: flow_condition%rate%itype = "' // &
          trim(adjustl(string)) // '", not implemented.'
        call printErrMsg(option)
    end select
  endif
  if(associated(flow_condition%energy_rate)) then
    select case (flow_condition%energy_rate%itype)
      case (ENERGY_RATE_SS)
        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
      case (HET_ENERGY_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%energy_rate%dataset, &
                num_connections,TH_TEMPERATURE_DOF,option)
      case default
        string = GetSubConditionName(flow_condition%energy_rate%itype)
        option%io_buffer='For TH mode: flow_condition%energy_rate%itype = "' // &
          trim(adjustl(string)) // '", not implemented.'
          write(*,*), trim(string)
        call printErrMsg(option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsTH

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for MIS_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsMIS(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  
  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif
  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MIS_PRESSURE_DOF, &
                                  1:num_connections) = &
          flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
    end select
  endif
  if (associated(flow_condition%concentration)) then
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (associated(flow_condition%concentration%dataset)) then
          coupler%flow_aux_real_var(MIS_CONCENTRATION_DOF, &
                                    1:num_connections) = &
            flow_condition%concentration%dataset%rarray(1)
        endif
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
    end select
  endif  

end subroutine PatchUpdateCouplerAuxVarsMIS

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVarsG: Updates flow auxiliary variables associated
!                             with a coupler for RICHARDS_MODE
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerAuxVarsRich(patch,coupler,option)

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module
  use Water_EOS_module
  
  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_class

  implicit none
  
  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  
  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
  if (option%mimetic) then
    num_connections = coupler%numfaces_set
  end if
#endif
  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (associated(flow_condition%pressure%dataset)) then
          coupler%flow_aux_real_var(RICHARDS_PRESSURE_DOF, &
                                    1:num_connections) = &
            flow_condition%pressure%dataset%rarray(1)
        else
          select type(dataset => &
                      flow_condition%pressure%dataset)
            class is(dataset_gridded_type)
              call PatchUpdateCouplerFromDataset(coupler,option, &
                                              patch%grid,dataset, &
                                              RICHARDS_PRESSURE_DOF)
            class default
          end select
        endif
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%saturation_function_array, &
                                 patch%sat_func_id)
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,option)
      case (HET_VOL_RATE_SS,HET_MASS_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%rate%dataset, &
                num_connections,RICHARDS_PRESSURE_DOF,option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsRich

! ************************************************************************** !
!
! PatchUpdateCouplerFromDataset: Updates auxiliary variables from dataset.
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine PatchUpdateCouplerFromDataset(coupler,option,grid,dataset,dof)

  use Option_module
  use Grid_module
  use Coupler_module
  use Dataset_Gridded_class
  
  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  class(dataset_gridded_type) :: dataset
  PetscInt :: dof
  
  PetscReal :: temp_real
  PetscInt :: iconn
  PetscInt :: local_id
  PetscInt :: ghosted_id
  
  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    call DatasetGriddedInterpolateReal(dataset, &
                                       grid%x(ghosted_id), &
                                       grid%y(ghosted_id), &
                                       grid%z(ghosted_id), &
                                       0.d0,temp_real,option)
    coupler%flow_aux_real_var(dof,iconn) = temp_real
  enddo
  
end subroutine PatchUpdateCouplerFromDataset

! ************************************************************************** !
!
! PatchScaleSourceSink: Scales select source/sinks based on perms*volume
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine PatchScaleSourceSink(patch,source_sink,option)

  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdmda.h"
  
  type(patch_type) :: patch
  type(coupler_type) :: source_sink
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(connection_set_type), pointer :: cur_connection_set
  type(field_type), pointer :: field
  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: perm_loc_ptr(:)
  PetscReal, pointer :: vol_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscInt :: iscale_type
  PetscReal :: scale, sum
  PetscInt :: icount, x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  PetscInt :: ghosted_neighbors(27)
    
  field => patch%field
  grid => patch%grid

  call VecGetArrayF90(field%perm_xx_loc,perm_loc_ptr,ierr)
  call VecGetArrayF90(field%volume,vol_ptr,ierr)

  grid => patch%grid

  call VecZeroEntries(field%work,ierr)
  call VecGetArrayF90(field%work,vec_ptr,ierr)

  cur_connection_set => source_sink%connection_set

  iscale_type = source_sink%flow_condition%rate%isubtype
  
  select case(iscale_type)
    case(SCALE_BY_VOLUME)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        vec_ptr(local_id) = vec_ptr(local_id) + vol_ptr(local_id)
      enddo
    case(SCALE_BY_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = vec_ptr(local_id) + perm_loc_ptr(ghosted_id) * &
                                                vol_ptr(local_id)
      enddo
    case(SCALE_BY_NEIGHBOR_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                    x_width,y_width,z_width, &
                                    x_count,y_count,z_count, &
                                    ghosted_neighbors,option)
        ! ghosted neighbors is ordered first in x, then, y, then z
        icount = 0
        sum = 0.d0
        ! x-direction
        do while (icount < x_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                      grid%structured_grid%dy(neighbor_ghosted_id)* &
                      grid%structured_grid%dz(neighbor_ghosted_id)
             
        enddo
        ! y-direction
        do while (icount < x_count + y_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)                 
          sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                      grid%structured_grid%dx(neighbor_ghosted_id)* &
                      grid%structured_grid%dz(neighbor_ghosted_id)
             
        enddo
        ! z-direction
        do while (icount < x_count + y_count + z_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)                 
          sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                      grid%structured_grid%dx(neighbor_ghosted_id)* &
                      grid%structured_grid%dy(neighbor_ghosted_id)
        enddo
        vec_ptr(local_id) = vec_ptr(local_id) + sum
      enddo
  end select

  call VecRestoreArrayF90(field%work,vec_ptr,ierr)
  call VecNorm(field%work,NORM_1,scale,ierr)
  scale = 1.d0/scale
  call VecScale(field%work,scale,ierr)

  call VecGetArrayF90(field%work,vec_ptr, ierr)
  do iconn = 1, cur_connection_set%num_connections      
    local_id = cur_connection_set%id_dn(iconn)
    select case(option%iflowmode)
      case(RICHARDS_MODE,G_MODE)
        source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
          vec_ptr(local_id)
      case(TH_MODE)
      case(THC_MODE)
      case(MPH_MODE)
      case(IMS_MODE)
      case(MIS_MODE)
      case(FLASH2_MODE)
    end select 

  enddo
  call VecRestoreArrayF90(field%work,vec_ptr,ierr)

  call VecRestoreArrayF90(field%perm_xx_loc,perm_loc_ptr, ierr)
  call VecRestoreArrayF90(field%volume,vol_ptr, ierr)
   
end subroutine PatchScaleSourceSink

! ************************************************************************** !
!> This subroutine updates aux vars for distributed copuler_type
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 10/03/2012
! ************************************************************************** !
subroutine PatchUpdateHetroCouplerAuxVars(patch,coupler,dataset_base, &
                                          sum_connection,isub_condition,option)

  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Dataset_module
  use Dataset_Map_class
  use Dataset_Base_class
  use Dataset_Ascii_class

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdmda.h"

  type(patch_type) :: patch
  type(coupler_type) :: coupler
  class(dataset_base_type), pointer :: dataset_base
  PetscInt :: isub_condition
  type(option_type) :: option

  type(connection_set_type), pointer :: cur_connection_set
  type(grid_type),pointer :: grid
  PetscErrorCode :: ierr
  PetscInt       :: iconn,sum_connection
  PetscInt :: ghosted_id,local_id
  PetscInt,pointer::cell_ids_nat(:)
  type(flow_sub_condition_type) :: flow_sub_condition

  class(dataset_map_type), pointer :: dataset_map
  class(dataset_ascii_type), pointer :: dataset_ascii

  grid => patch%grid
  
  if (isub_condition>option%nflowdof*option%nphase) then
    option%io_buffer='ERROR: PatchUpdateHetroCouplerAuxVars  '// &
      'isub_condition > option%nflowdof*option%nphase.'
    call printErrMsg(option)
  endif
  
  if(option%iflowmode/=RICHARDS_MODE.and.option%iflowmode/=TH_MODE) then
    option%io_buffer='PatchUpdateHetroCouplerAuxVars only implemented '// &
      ' for RICHARDS or TH mode.'
    call printErrMsg(option)
  endif

  cur_connection_set => coupler%connection_set

  select type(selector=>dataset_base)
    class is(dataset_map_type)
      dataset_map => selector

      ! If called for the first time, create the map
      if (dataset_map%first_time) then
        allocate(cell_ids_nat(cur_connection_set%num_connections))
        do iconn=1,cur_connection_set%num_connections
          sum_connection = sum_connection + 1
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          cell_ids_nat(iconn)=grid%nG2A(ghosted_id)
        enddo

        call PatchCreateFlowConditionDatasetMap(patch%grid,dataset_map,&
                cell_ids_nat,cur_connection_set%num_connections,option)

        dataset_map%first_time = PETSC_FALSE
        deallocate(cell_ids_nat)

      endif
    
      ! Save the data in the array
      do iconn=1,cur_connection_set%num_connections
        coupler%flow_aux_real_var(isub_condition,iconn) = &
          dataset_map%rarray(dataset_map%datatocell_ids(iconn))
      enddo

    class is(dataset_ascii_type)
      dataset_ascii => selector

      do iconn=1,cur_connection_set%num_connections
          coupler%flow_aux_real_var(isub_condition,iconn) = &
            dataset_ascii%rarray(1)
      enddo

    class default
      option%io_buffer = 'Incorrect dataset class (' // &
        trim(DatasetGetClass(dataset_base)) // &
        ') for coupler "' // trim(coupler%name) // &
        '" in PatchUpdateHetroCouplerAuxVars.'
      call printErrMsg(option)
  end select
  
end subroutine PatchUpdateHetroCouplerAuxVars

! ************************************************************************** !
!> This routine creates dataset-map for flow condition
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 10/26/12
! ************************************************************************** !
subroutine PatchCreateFlowConditionDatasetMap(grid,dataset_map,cell_ids,ncells,option)

  use Grid_module
  use Dataset_Map_class
  use Option_module
  
  implicit none
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(grid_type) :: grid
  class(dataset_map_type) :: dataset_map
  type(option_type):: option
  PetscInt,pointer :: cell_ids(:)
  PetscInt :: ncells
  
  PetscInt, allocatable :: int_array(:)
  PetscInt :: ghosted_id,local_id
  PetscInt :: ii,count
  PetscReal, pointer :: vec_ptr(:)  
  PetscErrorCode :: ierr
  PetscInt :: nloc,nglo
  PetscInt :: istart
  
  IS :: is_from, is_to
  Vec:: map_ids_1, map_ids_2,map_ids_3
  VecScatter::vec_scatter
  PetscViewer :: viewer
  
  ! Step-1: Rearrange map dataset
  nloc = maxval(dataset_map%mapping(2,:))
  call MPI_Allreduce(nloc,nglo,ONE_INTEGER,MPIU_INTEGER,MPI_Max,option%mycomm,ierr)
  call VecCreateMPI(option%mycomm,dataset_map%map_dims_local(2),&
                    PETSC_DETERMINE,map_ids_1,ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE,nglo,map_ids_2,ierr)
  call VecSet(map_ids_2,0,ierr)

  istart = 0
  call MPI_Exscan(dataset_map%map_dims_local(2), istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  allocate(int_array(dataset_map%map_dims_local(2)))
  do ii=1,dataset_map%map_dims_local(2)
    int_array(ii)=ii+istart
  enddo
  int_array=int_array-1
  
  call ISCreateBlock(option%mycomm,1,dataset_map%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_from,ierr)
  deallocate(int_array)
  
  allocate(int_array(dataset_map%map_dims_local(2)))
  do ii=1,dataset_map%map_dims_local(2)
    int_array(ii)=dataset_map%mapping(2,ii)
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,dataset_map%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_to,ierr)
  deallocate(int_array)

  !call VecCreateSeq(PETSC_COMM_SELF,dataset_map%map_dims_global(2),map_ids_1,ierr)
  !call VecCreateSeq(PETSC_COMM_SELF,maxval(dataset_map%map(2,:)),map_ids_2,ierr)
  !call VecSet(map_ids_2,0,ierr)

  call VecScatterCreate(map_ids_1,is_from,map_ids_2,is_to,vec_scatter,ierr)
  call ISDestroy(is_from,ierr)
  call ISDestroy(is_to,ierr)

  call VecGetArrayF90(map_ids_1,vec_ptr,ierr)
  do ii=1,dataset_map%map_dims_local(2)
    vec_ptr(ii)=dataset_map%mapping(1,ii)
  enddo
  call VecRestoreArrayF90(map_ids_1,vec_ptr,ierr)

  call VecScatterBegin(vec_scatter,map_ids_1,map_ids_2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,map_ids_1,map_ids_2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

  ! Step-2: Get ids in map dataset for cells
  allocate(int_array(ncells))
  allocate(dataset_map%cell_ids_local(ncells))
  int_array=cell_ids-1

  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_from,ierr)
    
  istart = 0
  call MPI_Exscan(ncells, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  do local_id=1,ncells
    int_array(local_id)=local_id+istart
  enddo
  int_array=int_array-1
  
  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_to,ierr)
  deallocate(int_array)
  
  !call VecCreateSeq(PETSC_COMM_SELF,ncells,map_ids_3,ierr)
  call VecCreateMPI(option%mycomm,ncells,PETSC_DETERMINE,map_ids_3,ierr)
  
  call VecScatterCreate(map_ids_2,is_from,map_ids_3,is_to,vec_scatter,ierr)
  call ISDestroy(is_from,ierr)
  call ISDestroy(is_to,ierr)

  call VecScatterBegin(vec_scatter,map_ids_2,map_ids_3, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,map_ids_2,map_ids_3, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

  ! Step-3: Save the datatocell_ids
  allocate(dataset_map%datatocell_ids(ncells))
  call VecGetArrayF90(map_ids_3,vec_ptr,ierr)
  do local_id=1,ncells
    dataset_map%datatocell_ids(local_id) = int(vec_ptr(local_id))
  enddo
  call VecRestoreArrayF90(map_ids_3,vec_ptr,ierr)
  
  call VecDestroy(map_ids_1,ierr)
  call VecDestroy(map_ids_2,ierr)
  call VecDestroy(map_ids_3,ierr)

end subroutine PatchCreateFlowConditionDatasetMap

! ************************************************************************** !
!
! PatchInitConstraints: Initializes constraint concentrations
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine PatchInitConstraints(patch,reaction,option)

  use Reaction_Aux_module
    
  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  type(reaction_type), pointer :: reaction
  
  call PatchInitCouplerConstraints(patch%initial_conditions, &
                                   reaction,option)
  
  call PatchInitCouplerConstraints(patch%boundary_conditions, &
                                   reaction,option)
  
  call PatchInitCouplerConstraints(patch%source_sinks, &
                                   reaction,option)

end subroutine PatchInitConstraints

! ************************************************************************** !
!
! PatchInitCouplerConstraints: Initializes constraint concentrations
!                              for a given coupler
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine PatchInitCouplerConstraints(coupler_list,reaction,option)

  use Reaction_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Constraint_module
  use Water_EOS_module
    
  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(option_type) :: option
  type(reaction_type), pointer :: reaction

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(coupler_type), pointer :: cur_coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  PetscReal :: r1, r2, r3, r4, r5, r6
  PetscErrorCode :: ierr
  
  cur_coupler => coupler_list%first
  do
    if (.not.associated(cur_coupler)) exit

    if (.not.associated(cur_coupler%tran_condition)) then
      option%io_buffer = 'Null transport condition found in coupler'
      if (len_trim(cur_coupler%name) > 1) then
        option%io_buffer = trim(option%io_buffer) // &
                           ' "' // trim(cur_coupler%name) // '"'
      endif
      call printErrMsg(option)
    endif

    cur_constraint_coupler => &
      cur_coupler%tran_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      global_auxvar => cur_constraint_coupler%global_auxvar
      rt_auxvar => cur_constraint_coupler%rt_auxvar
      if (associated(cur_coupler%flow_condition)) then
        if (associated(cur_coupler%flow_condition%pressure)) then
          if (associated(cur_coupler%flow_condition%pressure%dataset)) then
            global_auxvar%pres = &
              cur_coupler%flow_condition%pressure%dataset%rarray(1)
          else
            global_auxvar%pres = option%reference_pressure
          endif
        else
          global_auxvar%pres = option%reference_pressure
        endif
        if (associated(cur_coupler%flow_condition%temperature)) then
          if (associated(cur_coupler%flow_condition%temperature%dataset)) then
            global_auxvar%temp = &
              cur_coupler%flow_condition%temperature%dataset%rarray(1)
          else
            global_auxvar%temp = option%reference_temperature
          endif
        else
          global_auxvar%temp = option%reference_temperature
        endif

#ifndef DONT_USE_WATEOS
        call wateos(global_auxvar%temp(1),global_auxvar%pres(1), &
                    global_auxvar%den_kg(1),r1,r2,r3,r4,r5,r6, &
                    option%scale,ierr)
#else
        call density(global_auxvar%temp(1),global_auxvar%pres(1), &
                     global_auxvar%den_kg(1))
#endif                     
      else
        global_auxvar%pres = option%reference_pressure
        global_auxvar%temp = option%reference_temperature
        global_auxvar%den_kg = option%reference_water_density
      endif     
      global_auxvar%sat = option%reference_saturation  
  
      call ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                            reaction,cur_constraint_coupler%constraint_name, &
                            cur_constraint_coupler%aqueous_species, &
                            cur_constraint_coupler%minerals, &
                            cur_constraint_coupler%surface_complexes, &
                            cur_constraint_coupler%colloids, &
                            cur_constraint_coupler%immobile_species, &
                            option%reference_porosity, &
                            cur_constraint_coupler%num_iterations, &
                            PETSC_FALSE,option)
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    cur_coupler => cur_coupler%next
  enddo

end subroutine PatchInitCouplerConstraints

! ************************************************************************** !
!
! PatchUpdateUniformVelocity: Assigns uniform velocity in connection list
!                        darcy velocities
! author: Glenn Hammond
! date: 02/20/08
!
! ************************************************************************** !
subroutine PatchUpdateUniformVelocity(patch,velocity,option)

  use Option_module
  use Coupler_module
  use Condition_module
  use Connection_module
  
  implicit none
  
  type(patch_type), pointer :: patch   
  PetscReal :: velocity(3)
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: vdarcy

  grid => patch%grid
    
  ! Internal Flux Terms -----------------------------------
  cur_connection_set => grid%internal_connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = dot_product(velocity, &
                           cur_connection_set%dist(1:3,iconn))
      patch%internal_velocities(1,sum_connection) = vdarcy
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = dot_product(velocity, &
                           cur_connection_set%dist(1:3,iconn))
      patch%boundary_velocities(1,sum_connection) = vdarcy
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PatchUpdateUniformVelocity

! ************************************************************************** !
!
! PatchAuxVarsUpToDate: Checks to see if aux vars are up to date
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function PatchAuxVarsUpToDate(patch)

  use Grid_module
  use Option_module
  use Field_module
  
  use Mphase_Aux_module
  use TH_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module  
  
  type(patch_type) :: patch
  
  PetscBool :: PatchAuxVarsUpToDate
  PetscBool :: flow_up_to_date
  PetscBool :: transport_up_to_date
  PetscInt :: dummy
  dummy = 1

  if (associated(patch%aux%THC)) then
    flow_up_to_date = patch%aux%THC%aux_vars_up_to_date
  else if (associated(patch%aux%TH)) then
    flow_up_to_date = patch%aux%TH%aux_vars_up_to_date
  else if (associated(patch%aux%Richards)) then
    flow_up_to_date = patch%aux%Richards%aux_vars_up_to_date
  else if (associated(patch%aux%Mphase)) then
    flow_up_to_date = patch%aux%Mphase%aux_vars_up_to_date
  else if (associated(patch%aux%Flash2)) then
    flow_up_to_date = patch%aux%Flash2%aux_vars_up_to_date
  else if (associated(patch%aux%Immis)) then
    flow_up_to_date = patch%aux%Immis%aux_vars_up_to_date
  endif

  if (associated(patch%aux%RT)) then
    transport_up_to_date = patch%aux%RT%aux_vars_up_to_date
  endif
  
  PatchAuxVarsUpToDate = flow_up_to_date .and. transport_up_to_date
  
end function PatchAuxVarsUpToDate

! ************************************************************************** !
!
! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchGetVariable1(patch,field,reaction,option,output_option,vec,ivar, &
                           isubvar,isubvar1)

  use Grid_module
  use Option_module
  use Field_module
  
  use Immis_Aux_module
  use Miscible_Aux_module
  use Mphase_Aux_module
  use TH_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Mineral_module
  use Reaction_module
  use Reactive_Transport_Aux_module  
  use Surface_Complexation_Aux_module
  use General_Aux_module
  use Output_Aux_module
  use Variables_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass, lnQKgas, ehfac, eh0, pe0, ph0, tk
  PetscReal :: tempreal
  PetscInt :: tempint
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  PetscErrorCode :: ierr

  grid => patch%grid

  call VecGetArrayF90(vec,vec_ptr,ierr)

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY,GAS_VISCOSITY, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,STATE,ICE_DENSITY)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY,GAS_VISCOSITY) ! still needs implementation
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = 0.d0
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
#ifdef ICE
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat_gas
#else
              vec_ptr(local_id) = 0.d0
#endif 
            enddo
#ifdef ICE
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat_ice
            enddo
          case(ICE_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_ice*FMWH2O
            enddo
#endif
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%vis
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%kvr
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%u
            enddo
        end select

      else if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY,GAS_VISCOSITY) ! still needs implementation
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = 0.d0
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
#ifdef ICE
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%sat_gas
#else
              vec_ptr(local_id) = 0.d0
#endif 
            enddo
#ifdef ICE
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%sat_ice
            enddo
          case(ICE_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%den_ice*FMWH2O
            enddo
#endif
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%vis
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%kvr
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%aux_vars(grid%nL2G(local_id))%u
            enddo
        end select
        
      else if (associated(patch%aux%Richards)) then
      
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(ICE_SATURATION)
            call printErrMsg(option,'ICE_SATURATION not supported by Richards')
          case(ICE_DENSITY)
            call printErrMsg(option,'ICE_DENSITY not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(LIQUID_VISCOSITY)
            call printErrMsg(option,'LIQUID_VISCOSITY not supported by Richards')
          case(GAS_VISCOSITY)
            call printErrMsg(option,'GAS_VISCOSITY not supported by Richards')
          case(LIQUID_MOBILITY)
            call printErrMsg(option,'LIQUID_MOBILITY not supported by Richards')
          case(GAS_MOBILITY)
            call printErrMsg(option,'GAS_MOBILITY not supported by Richards')
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
        end select
        
      else if (associated(patch%aux%Flash2)) then
      
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(2)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(2)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
            enddo
          case(GAS_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(2)
            enddo
          case(GAS_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(2)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(2)
            enddo
          case(GAS_DENSITY_MOL) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den(2)
            enddo
          case(SC_FUGA_COEFF)
            if (.not.associated(patch%aux%Global%aux_vars(1)%fugacoeff) .and. &
                OptionPrintToScreen(option))then
               print *,'ERRor, fugacoeff not allocated for ', option%iflowmode, 1
            endif
            do local_id=1,grid%nlmax
             vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%fugacoeff(1)
            enddo 
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(1)
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(1)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1)
            enddo
        end select
        
      else if (associated(patch%aux%Mphase)) then
        
        select case(ivar)
        
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(2)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(1)
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(1)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(2)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
            enddo
          case(GAS_VISCOSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(2)
            enddo
          case(GAS_MOBILITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(2)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(2)
            enddo
          case(GAS_DENSITY_MOL) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den(2)
            enddo
          case(SC_FUGA_COEFF)
            if (.not.associated(patch%aux%Global%aux_vars(1)%fugacoeff) .and. &
                OptionPrintToScreen(option))then
               print *,'ERRor, fugacoeff not allocated for ', option%iflowmode, 1
            endif
            do local_id=1,grid%nlmax
             vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%fugacoeff(1)
            enddo 
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1)
            enddo
        end select
        
      else if (associated(patch%aux%Miscible)) then
        
        select case(ivar)
        
!         case(TEMPERATURE)
!           do local_id=1,grid%nlmax
!             vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
!           enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Miscible%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(1)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Miscible%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar)
            enddo
        end select
        
      else if (associated(patch%aux%immis)) then
      
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(2)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(2)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(2)
            enddo
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%temp
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%pres(option%liquid_phase)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%pres(option%gas_phase)
            enddo
          case(STATE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%istate
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%liquid_phase)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%U(option%liquid_phase)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%liquid_phase)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%gas_phase)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%U(option%gas_phase)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%gas_phase)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%gas_phase)
            enddo
        end select         
      endif
      
    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY,TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_RATE,MINERAL_VOLUME_FRACTION,MINERAL_SATURATION_INDEX, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK, &
         IMMOBILE_SPECIES)

      select case(ivar)

        case(PH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              vec_ptr(local_id) = &
               -log10(patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)* &
                      patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo

        case(EH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then

              ph0 = -log10(patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)* &
                      patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))

              ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
              lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
              if (reaction%eqgash2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%eqgasspecid(0,ifo2)
                comp_id = reaction%eqgasspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
              enddo

              tk = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)+273.15d0
              ehfac = IDEAL_GAS_CONST*tk*LOG_TO_LN/faraday
              eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
              pe0 = eh0/ehfac
              vec_ptr(local_id) = eh0

            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo

        case(PE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then

              ph0 = -log10(patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)* &
                      patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))

              ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
              lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
              if (reaction%eqgash2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%eqgasspecid(0,ifo2)
                comp_id = reaction%eqgasspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
              enddo

              tk = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)+273.15d0
              ehfac = IDEAL_GAS_CONST*tk*LOG_TO_LN/faraday
              eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
              pe0 = eh0/ehfac
              vec_ptr(local_id) = pe0
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo

        case(O2)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then


              ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
              lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
              if (reaction%eqgash2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%eqgasspecid(0,ifo2)
                comp_id = reaction%eqgasspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
              enddo

              vec_ptr(local_id) = lnQKgas * LN_TO_LOG
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo

        case(PRIMARY_MOLALITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%pri_molal(isubvar)
          enddo
        case(PRIMARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%aux_vars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%aux_vars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) * xmass * &
              (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(SECONDARY_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)
          enddo
        case(SECONDARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%aux_vars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%aux_vars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar) * xmass * &
              (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(TOTAL_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id =grid%nL2G(local_id)
            if (associated(patch%aux%Global%aux_vars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%aux_vars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            if (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase) > 0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase) / xmass / &
                (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(TOTAL_MOLARITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%total(isubvar,iphase)
          enddo
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase) * &
              vec_ptr2(ghosted_id) * &
              patch%aux%Global%aux_vars(ghosted_id)%sat(iphase) * 1.d-3 ! mol/L -> mol/m^3
          enddo
          call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
          ! add in total sorbed.  already in mol/m^3 bulk
          if (patch%reaction%nsorb > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (patch%reaction%surface_complexation%neqsrfcplxrxn > 0) then
                vec_ptr(local_id) = vec_ptr(local_id) + &
                  patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
              endif
              if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
                do irxn = 1, &
                   patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                  do irate = 1, &
                     patch%reaction%surface_complexation%kinmr_nrate(irxn)
                    vec_ptr(local_id) = vec_ptr(local_id) + &
                      patch%aux%RT%aux_vars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                  enddo 
                enddo
              endif
            enddo
          endif
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_volfrac(isubvar)
          enddo
        case(MINERAL_RATE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_rate(isubvar)
          enddo
        case(MINERAL_SATURATION_INDEX)
          do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = RMineralSaturationIndex(isubvar, &
                                                  patch%aux%RT%aux_vars(ghosted_id), &
                                                  patch%aux%Global%aux_vars(ghosted_id), &
                                                  reaction,option)
          enddo
        case(IMMOBILE_SPECIES)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%immobile(isubvar)
          enddo          
        case(SURFACE_CMPLX)
          if (associated(patch%aux%RT%aux_vars(1)%eqsrfcplx_conc)) then
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))% &
                                    eqsrfcplx_conc(isubvar)
            enddo
          else
            vec_ptr = -999.d0
          endif
        case(SURFACE_SITE_DENSITY)
          tempreal = reaction%surface_complexation%srfcplxrxn_site_density(isubvar)
          select case(reaction%surface_complexation%srfcplxrxn_surf_type(isubvar))
            case(MINERAL_SURFACE)
              tempint = reaction%surface_complexation%srfcplxrxn_to_surf(isubvar)
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = tempreal* &
                                    patch%aux%RT%aux_vars(grid%nL2G(local_id))% &
                                      mnrl_volfrac(tempint)
              enddo
            case(COLLOID_SURFACE)
                option%io_buffer = 'Printing of surface site density for colloidal surfaces ' // &
                  'not implemented.'
                call printErrMsg(option)
            case(NULL_SURFACE)
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = tempreal
              enddo
          end select          
        case(SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))% &
              srfcplxrxn_free_site_conc(isubvar)
          enddo
        case(KIN_SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))% &
              kinsrfcplx_conc(isubvar,1)
          enddo
        case(KIN_SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))% &
              kinsrfcplx_free_site_conc(isubvar)
          enddo
        case(PRIMARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)
          enddo
        case(SECONDARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%aux_vars(ghosted_id)%sec_act_coef(isubvar)
          enddo
        case(PRIMARY_KD)
          call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call ReactionComputeKd(isubvar,vec_ptr(local_id), &
                                   patch%aux%RT%aux_vars(ghosted_id), &
                                   patch%aux%Global%aux_vars(ghosted_id), &
                                   vec_ptr2(ghosted_id),patch%reaction,option)
          enddo
          call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
              enddo
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = 0.d0
                do irxn = 1, &
                  patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                  do irate = 1, &
                    patch%reaction%surface_complexation%kinmr_nrate(irxn)
                    vec_ptr(local_id) = vec_ptr(local_id) + &
                      patch%aux%RT%aux_vars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                  enddo            
                enddo            
              enddo
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = patch%aux%RT%aux_vars(ghosted_id)%colloid% &
                total_eq_mob(isubvar)
            enddo
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase) > 0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%aux_vars(grid%nL2G(local_id))%colloid%conc_mob(isubvar) / &
                  (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%RT%aux_vars(grid%nL2G(local_id))%colloid%conc_mob(isubvar)
            enddo
          endif      
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase) > 0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%aux_vars(grid%nL2G(local_id))%colloid%conc_imb(isubvar) / &
                  (patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%RT%aux_vars(grid%nL2G(local_id))%colloid%conc_imb(isubvar)
            enddo
          endif
        case(AGE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) / &
                patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar1) / &
                output_option%tconv
            endif
          enddo        
      end select
    case(POROSITY)
      call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
    case(PERMEABILITY,PERMEABILITY_X)
      call VecGetArrayF90(field%perm_xx_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_xx_loc,vec_ptr2,ierr)
    case(PERMEABILITY_Y)
      call VecGetArrayF90(field%perm_yy_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_yy_loc,vec_ptr2,ierr)
    case(PERMEABILITY_Z)
      call VecGetArrayF90(field%perm_zz_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_zz_loc,vec_ptr2,ierr)
    case(PERMEABILITY_XY)
      call VecGetArrayF90(field%perm_xy_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_xy_loc,vec_ptr2,ierr)
    case(PERMEABILITY_XZ)
      call VecGetArrayF90(field%perm_xz_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_xz_loc,vec_ptr2,ierr)
    case(PERMEABILITY_YZ)
      call VecGetArrayF90(field%perm_yz_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%perm_yz_loc,vec_ptr2,ierr)
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
    case(PROCESSOR_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariable'')') ivar
      call printErrMsg(option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine PatchGetVariable1

! ************************************************************************** !
!
! PatchGetVariableValueAtCell: Returns variables indexed by ivar,
!                             isubvar, local id from Reactive Transport type
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function PatchGetVariableValueAtCell(patch,field,reaction,option, &
                                    output_option, &
                                    ivar,isubvar,ghosted_id,isubvar1)

  use Grid_module
  use Option_module
  use Field_module

  use Mphase_Aux_module
  use TH_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Miscible_Aux_module
  use Reactive_Transport_Aux_module  
  use Mineral_module
  use Reaction_module
  use Mineral_Aux_module
  use Surface_Complexation_Aux_module
  use Output_Aux_module
  use Variables_module
  use General_Aux_module  

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscReal :: PatchGetVariableValueAtCell
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: iphase
  PetscInt :: ghosted_id, local_id

  PetscReal :: value, xmass, lnQKgas, tk, ehfac, eh0, pe0, ph0
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr2(:)  
  PetscErrorCode :: ierr

  grid => patch%grid
  
  value = -999.99d0

  ! inactive grid cell
  if (patch%imat(ghosted_id) <= 0) then
    PatchGetVariableValueAtCell = 0.d0
    return
  endif

  iphase = 1
  xmass = 1.d0
  if (associated(patch%aux%Global%aux_vars(ghosted_id)%xmass)) &
    xmass = patch%aux%Global%aux_vars(ghosted_id)%xmass(iphase)
             
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY,GAS_VISCOSITY, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,STATE,ICE_DENSITY, &
         SECONDARY_TEMPERATURE)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%THC%aux_vars(ghosted_id)%vis
          case(LIQUID_MOBILITY)
            value = patch%aux%THC%aux_vars(ghosted_id)%kvr
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            value = 0.d0
          case(GAS_SATURATION)
#ifdef ICE
            value = patch%aux%THC%aux_vars(ghosted_id)%sat_gas
#else
            value = 0.d0
#endif
#ifdef ICE
          case(ICE_SATURATION)
            value = patch%aux%THC%aux_vars(ghosted_id)%sat_ice
          case(ICE_DENSITY)
            value = patch%aux%THC%aux_vars(ghosted_id)%den_ice*FMWH2O
#endif
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%THC%aux_vars(ghosted_id)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%THC%aux_vars(ghosted_id)%u
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
        end select
     else if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%TH%aux_vars(ghosted_id)%vis
          case(LIQUID_MOBILITY)
            value = patch%aux%TH%aux_vars(ghosted_id)%kvr
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            value = 0.d0
          case(GAS_SATURATION)
#ifdef ICE
            value = patch%aux%TH%aux_vars(ghosted_id)%sat_gas
#else
            value = 0.d0
#endif
#ifdef ICE
          case(ICE_SATURATION)
            value = patch%aux%TH%aux_vars(ghosted_id)%sat_ice
          case(ICE_DENSITY)
            value = patch%aux%TH%aux_vars(ghosted_id)%den_ice*FMWH2O
#endif
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%TH%aux_vars(ghosted_id)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%TH%aux_vars(ghosted_id)%u
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
        end select
      else if (associated(patch%aux%Flash2)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1)
          case(LIQUID_MOBILITY)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(2)
          case(GAS_MOLE_FRACTION)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar)
          case(GAS_ENERGY)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(2)
          case(GAS_DENSITY_MOL) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den(2)
          case(GAS_VISCOSITY) 
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%vis(2)
          case(GAS_MOBILITY) 
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(2)
          case(SC_FUGA_COEFF)
            value = patch%aux%Global%aux_vars(ghosted_id)%fugacoeff(1)   
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
          case(GAS_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1)
          case(LIQUID_MOBILITY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(2)
          case(GAS_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar)
          case(GAS_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(2)
          case(GAS_VISCOSITY) 
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%vis(2)
          case(GAS_MOBILITY) 
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(2)
          case(GAS_DENSITY_MOL) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den(2)
          case(SC_FUGA_COEFF)
            value = patch%aux%Global%aux_vars(ghosted_id)%fugacoeff(1)   
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
        end select
      else if (associated(patch%aux%Immis)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(LIQUID_ENERGY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1)
          case(LIQUID_MOBILITY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(2)
          case(GAS_ENERGY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(2)
          case(GAS_DENSITY_MOL) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den(2)
          case(GAS_VISCOSITY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%vis(2)
          case(GAS_MOBILITY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(2)
        end select
      else if (associated(patch%aux%Miscible)) then
        select case(ivar)
!         case(TEMPERATURE)
!           value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(1)
!         case(LIQUID_SATURATION)
!           value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
!         case(LIQUID_ENERGY)
!           value = patch%aux%Miscible%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Miscible%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1)
!         case(LIQUID_MOBILITY)
!           value = patch%aux%Miscible%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1)
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%Miscible%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%temp
          case(LIQUID_PRESSURE)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%liquid_phase)
          case(GAS_PRESSURE)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%gas_phase)
          case(STATE)
            value = patch%aux%Global%aux_vars(ghosted_id)%istate
          case(LIQUID_SATURATION)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%sat(option%liquid_phase)
          case(LIQUID_DENSITY)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%den_kg(option%liquid_phase)
          case(LIQUID_ENERGY)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%U(option%liquid_phase)
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%xmol(isubvar,option%liquid_phase)
          case(GAS_SATURATION)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%sat(option%gas_phase)
          case(GAS_DENSITY) 
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%den_kg(option%gas_phase)
          case(GAS_ENERGY)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%U(option%gas_phase)
          case(GAS_MOLE_FRACTION)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%xmol(isubvar,option%gas_phase)
        end select        
      endif
      
    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY, &
         TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_VOLUME_FRACTION,MINERAL_RATE,MINERAL_SATURATION_INDEX, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK, &
         IMMOBILE_SPECIES)
         
      select case(ivar)
        case(PH)
          value = -log10(patch%aux%RT%aux_vars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))
        case(EH)
          ph0 = -log10(patch%aux%RT%aux_vars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))

          ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
          lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
          if (reaction%eqgash2oid(ifo2) > 0) then
            lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
          endif
          do jcomp = 1, reaction%eqgasspecid(0,ifo2)
            comp_id = reaction%eqgasspecid(jcomp,ifo2)
            lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
          enddo

          tk = patch%aux%Global%aux_vars(grid%nL2G(ghosted_id))%temp(1)+273.15d0
          ehfac = IDEAL_GAS_CONST*tk*LOG_TO_LN/faraday
          eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0

          value = eh0

        case(PE)
          ph0 = -log10(patch%aux%RT%aux_vars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))

          ifo2 = reaction%species_idx%o2_gas_id
      
      ! compute gas partial pressure
          lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
          if (reaction%eqgash2oid(ifo2) > 0) then
            lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
          endif
          do jcomp = 1, reaction%eqgasspecid(0,ifo2)
            comp_id = reaction%eqgasspecid(jcomp,ifo2)
            lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
          enddo

          tk = patch%aux%Global%aux_vars(grid%nL2G(ghosted_id))%temp(1)+273.15d0
          ehfac = IDEAL_GAS_CONST*tk*LOG_TO_LN/faraday
          eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
          pe0 = eh0/ehfac
          value = pe0

        case(O2)
      
      ! compute gas partial pressure
              ifo2 = reaction%species_idx%o2_gas_id
              lnQKgas = -reaction%eqgas_logK(ifo2)*LOG_TO_LN
      
      ! activity of water
              if (reaction%eqgash2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%eqgash2ostoich(ifo2) * &
                    patch%aux%RT%aux_vars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%eqgasspecid(0,ifo2)
                comp_id = reaction%eqgasspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%eqgasstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%aux_vars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(comp_id))
              enddo
           value = lnQKgas * LN_TO_LOG
        case(PRIMARY_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar)
        case(PRIMARY_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) * xmass * &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(SECONDARY_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar)
        case(SECONDARY_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_molal(isubvar) * xmass * &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(TOTAL_MOLALITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase) / xmass / &
                  patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)*1000.d0
        case(TOTAL_MOLARITY)
          value = patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase)
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
          value = &
              patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase) * &
              vec_ptr2(ghosted_id) * &
              patch%aux%Global%aux_vars(ghosted_id)%sat(iphase) * 1.d-3 ! mol/L -> mol/m^3
          call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
          ! add in total sorbed.  already in mol/m^3 bulk
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%surface_complexation%neqsrfcplxrxn > 0) then
              value = value + &
                patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              do irxn = 1, patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                do irate = 1, &
                  patch%reaction%surface_complexation%kinmr_nrate(irxn)
                  value = value + &
                      patch%aux%RT%aux_vars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                enddo
              enddo            
            endif
          endif          
        case(MINERAL_VOLUME_FRACTION)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(isubvar)
        case(MINERAL_RATE)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_rate(isubvar)
        case(MINERAL_SATURATION_INDEX)
          value = RMineralSaturationIndex(isubvar,patch%aux%RT%aux_vars(ghosted_id), &
                                          patch%aux%Global%aux_vars(ghosted_id), &
                                          reaction,option)
        case(IMMOBILE_SPECIES)
          value = patch%aux%RT%aux_vars(ghosted_id)%immobile(isubvar)
        case(SURFACE_CMPLX)
          if (associated(patch%aux%RT%aux_vars(ghosted_id)%eqsrfcplx_conc)) then
            value = patch%aux%RT%aux_vars(ghosted_id)%eqsrfcplx_conc(isubvar)
          else
            value = -999.d0
          endif
        case(SURFACE_CMPLX_FREE)
          value = patch%aux%RT%aux_vars(ghosted_id)%srfcplxrxn_free_site_conc(isubvar)
        case(SURFACE_SITE_DENSITY)
          select case(reaction%surface_complexation%srfcplxrxn_surf_type(isubvar))
            case(MINERAL_SURFACE)
              value = reaction%surface_complexation%srfcplxrxn_site_density(isubvar)* &
                      patch%aux%RT%aux_vars(ghosted_id)% &
                        mnrl_volfrac(reaction%surface_complexation%srfcplxrxn_to_surf(isubvar))
            case(COLLOID_SURFACE)
                option%io_buffer = 'Printing of surface site density for colloidal surfaces ' // &
                  'not implemented.'
                call printErrMsg(option)
            case(NULL_SURFACE)
              value = reaction%surface_complexation%srfcplxrxn_site_density(isubvar)
          end select
        case(KIN_SURFACE_CMPLX)
          value = patch%aux%RT%aux_vars(ghosted_id)%kinsrfcplx_conc(isubvar,1)
        case(KIN_SURFACE_CMPLX_FREE)
          value = patch%aux%RT%aux_vars(ghosted_id)%kinsrfcplx_free_site_conc(isubvar)
        case(PRIMARY_ACTIVITY_COEF)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)
        case(SECONDARY_ACTIVITY_COEF)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_act_coef(isubvar)
        case(PRIMARY_KD)
          call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
          call ReactionComputeKd(isubvar,value, &
                                 patch%aux%RT%aux_vars(ghosted_id), &
                                 patch%aux%Global%aux_vars(ghosted_id), &
                                 vec_ptr2(ghosted_id),patch%reaction,option)
          call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              value = patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              value = 0.d0
              do irxn = 1, patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                do irate = 1, &
                  patch%reaction%surface_complexation%kinmr_nrate(irxn)
                  value = value + &
                    patch%aux%RT%aux_vars(ghosted_id)% &
                      kinmr_total_sorb(isubvar,irate,irxn)
                enddo
              enddo            
            endif
          endif          
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            value = &
              patch%aux%RT%aux_vars(ghosted_id)%colloid%total_eq_mob(isubvar)
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_mob(isubvar) / &
                    patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_mob(isubvar)
          endif
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_imb(isubvar) / &
                    patch%aux%Global%aux_vars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%aux_vars(ghosted_id)%colloid%conc_imb(isubvar)
          endif
        case(AGE)
          if (patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) > &
              0.d0) then
            value = patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) / &
            patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar1) / &
            output_option%tconv
          endif
      end select
    case(POROSITY)
      call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
    case(PERMEABILITY,PERMEABILITY_X)
      call VecGetArrayF90(field%perm_xx_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%perm_xx_loc,vec_ptr2,ierr)
    case(PERMEABILITY_Y)
      call VecGetArrayF90(field%perm_yy_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%perm_yy_loc,vec_ptr2,ierr)
    case(PERMEABILITY_Z)
      call VecGetArrayF90(field%perm_zz_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%perm_zz_loc,vec_ptr2,ierr)
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      value = patch%imat(ghosted_id)
    case(PROCESSOR_ID)
      value = option%myrank
    ! Need to fix the below two cases (they assume only one component) -- SK 02/06/13  
    case(SECONDARY_CONCENTRATION)
      ! Note that the units are in mol/kg
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%pri_molal(isubvar1)
    case(SEC_MIN_VOLFRAC)
      local_id = grid%nG2L(ghosted_id)        
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%mnrl_volfrac(isubvar1)
     case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariableValueAtCell'')') &
            ivar
      call printErrMsg(option)
  end select

  PatchGetVariableValueAtCell = value
 
end function PatchGetVariableValueAtCell

! ************************************************************************** !
!
! PatchSetVariable: Sets variables indexed by ivar and isubvar within a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchSetVariable(patch,field,option,vec,vec_format,ivar,isubvar)

  use Grid_module
  use Option_module
  use Field_module
  use Variables_module
  use General_Aux_module  

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase, istate

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid

  call VecGetArrayF90(vec,vec_ptr,ierr)

  if (vec_format == NATURAL) then
    call printErrMsg(option,&
                     'NATURAL vector format not supported by PatchSetVariable')
  endif

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY,GAS_VISCOSITY, &
         LIQUID_MOBILITY,GAS_MOBILITY,STATE)
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%sat = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%den_kg = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
          case(GAS_SATURATION)
#ifdef ICE
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat_gas = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%sat_gas = vec_ptr(ghosted_id)
              enddo
            endif
#else
#endif
#ifdef ICE
          case(ICE_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat_ice = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%sat_ice = vec_ptr(ghosted_id)
              enddo
            endif
          case(ICE_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_ice = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%den_ice = vec_ptr(ghosted_id)
              enddo
            endif
#endif
          case(LIQUID_VISCOSITY)
          case(GAS_VISCOSITY)
          case(LIQUID_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%xmol(isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%xmol(isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%u = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%u = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%sat = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%den_kg = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
          case(GAS_SATURATION)
#ifdef ICE
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%aux_vars(grid%nL2G(local_id))%sat_gas = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%aux_vars(ghosted_id)%sat_gas = vec_ptr(ghosted_id)
              enddo
            endif
#else
#endif
#ifdef ICE
          case(ICE_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%aux_vars(grid%nL2G(local_id))%sat_ice = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%aux_vars(ghosted_id)%sat_ice = vec_ptr(ghosted_id)
              enddo
            endif
          case(ICE_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%aux_vars(grid%nL2G(local_id))%den_ice = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%aux_vars(ghosted_id)%den_ice = vec_ptr(ghosted_id)
              enddo
            endif
#endif
          case(LIQUID_VISCOSITY)
          case(GAS_VISCOSITY)
          case(LIQUID_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%aux_vars(grid%nL2G(local_id))%xmol(isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%aux_vars(ghosted_id)%xmol(isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%aux_vars(grid%nL2G(local_id))%u = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%aux_vars(ghosted_id)%u = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call printErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call printErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call printErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_VISCOSITY)
            call printErrMsg(option,'LIQUID_VISCOSITY not supported by Richards')
          case(GAS_VISCOSITY)
            call printErrMsg(option,'GAS_VISCOSITY not supported by Richards')
          case(LIQUID_MOBILITY)
            call printErrMsg(option,'LIQUID_MOBILITY not supported by Richards')
          case(GAS_MOBILITY)
            call printErrMsg(option,'GAS_MOBILITY not supported by Richards')
          case(LIQUID_ENERGY)
            call printErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call printErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%pres(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%sat(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%sat(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(ghosted_id)%den_kg(1) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Flash2)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%sat(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%den(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_VISCOSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOBILITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%sat(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%u(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_DENSITY, GAS_DENSITY_MOL) 
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%den(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_VISCOSITY) 
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%vis(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOBILITY) 
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flash2%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flash2%aux_vars(ghosted_id)%aux_var_elem(0)%u(1) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(grid%nL2G(ghosted_id))%pres(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%aux_vars(grid%nL2G(local_id))%pres(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%aux_vars(grid%nL2G(ghosted_id))%pres(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%sat(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_VISCOSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%vis(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%vis(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOBILITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%kvr(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%kvr(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%sat(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_DENSITY, GAS_DENSITY_MOL) 
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%den(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(1) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Immis)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%sat(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%den(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%den(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%sat(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%sat(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(1) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%temp = vec_ptr(local_id)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%pres(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%pres(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(STATE)
            do local_id=1,grid%nlmax
              patch%aux%Global%aux_vars(grid%nL2G(local_id))%istate = int(vec_ptr(local_id)+1.d-10)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%sat(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%den_kg(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%U(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%xmol(isubvar,option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%sat(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%den_kg(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%U(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%xmol(isubvar,option%gas_phase) = vec_ptr(local_id)
            enddo
        end select         
      endif
    case(PRIMARY_MOLALITY,TOTAL_MOLARITY,MINERAL_VOLUME_FRACTION, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,IMMOBILE_SPECIES)
      select case(ivar)
        case(PRIMARY_MOLALITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%pri_molal(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(TOTAL_MOLARITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%total(isubvar,iphase) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%total(isubvar,iphase) = vec_ptr(ghosted_id)
            enddo
          endif
        case(MINERAL_VOLUME_FRACTION)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_volfrac(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(IMMOBILE_SPECIES)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%immobile(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%immobile(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(PRIMARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%pri_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(SECONDARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%aux_vars(grid%nL2G(local_id))%sec_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%aux_vars(ghosted_id)%sec_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
      end select
    case(PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY,TOTAL_MOLALITY, &
         COLLOID_MOBILE,COLLOID_IMMOBILE)
      select case(ivar)
        case(PRIMARY_MOLARITY)
          call printErrMsg(option,'Setting of primary molarity at grid cell not supported.')
        case(SECONDARY_MOLALITY)
          call printErrMsg(option,'Setting of secondary molality at grid cell not supported.')
        case(SECONDARY_MOLARITY)
          call printErrMsg(option,'Setting of secondary molarity at grid cell not supported.')
        case(TOTAL_MOLALITY)
          call printErrMsg(option,'Setting of total molality at grid cell not supported.')
        case(COLLOID_MOBILE)
          call printErrMsg(option,'Setting of mobile colloid concentration at grid cell not supported.')
        case(COLLOID_IMMOBILE)
          call printErrMsg(option,'Setting of immobile colloid concentration at grid cell not supported.')
      end select
    case(POROSITY)
      if (vec_format == GLOBAL) then
        call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
        do local_id=1,grid%nlmax
          vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
        call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
      else if (vec_format == LOCAL) then
        call VecGetArrayF90(field%porosity_loc,vec_ptr2,ierr)
        vec_ptr2(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
        call VecRestoreArrayF90(field%porosity_loc,vec_ptr2,ierr)
      endif
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y,PERMEABILITY_Z)
      option%io_buffer = 'Setting of permeability in "PatchSetVariable"' // &
        ' not supported.'
      call printErrMsg(option)
    case(PHASE)
      if (vec_format == GLOBAL) then
        call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr)
        do local_id=1,grid%nlmax
          vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
        call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr)
      else if (vec_format == LOCAL) then
        call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr)
        vec_ptr2(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
        call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr)
      endif
    case(MATERIAL_ID)
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          patch%imat(grid%nL2G(local_id)) = int(vec_ptr(local_id))
        enddo
      else if (vec_format == LOCAL) then
        patch%imat(1:grid%ngmax) = int(vec_ptr(1:grid%ngmax))
      endif
    case(PROCESSOR_ID)
      call printErrMsg(option, &
                       'Cannot set PROCESSOR_ID through PatchSetVariable()')
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr)
  
end subroutine PatchSetVariable

! ************************************************************************** !
!
! PatchCountCells: Counts # of active and inactive grid cells 
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine PatchCountCells(patch,total_count,active_count)

  use Option_module

  implicit none
  
  type(patch_type) :: patch
  PetscInt :: total_count
  PetscInt :: active_count
  
  type(grid_type), pointer :: grid
  PetscInt :: local_id
  
  grid => patch%grid
  
  total_count = grid%nlmax
  
  active_count = 0
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    active_count = active_count + 1
  enddo

end subroutine PatchCountCells

! ************************************************************************** !
!
! PatchCalculateCFL1Timestep: Calculates largest time step to preserves a 
!                                CFL # of 1 in a patch
! author: Glenn Hammond
! date: 10/06/11
!
! ************************************************************************** !
subroutine PatchCalculateCFL1Timestep(patch,option,max_dt_cfl_1)

  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Global_Aux_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(patch_type) :: patch
  type(option_type) :: option
  PetscReal :: max_dt_cfl_1
  
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(global_auxvar_type), pointer :: global_aux_vars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: por_sat_ave, por_sat_min, v_darcy, v_pore_ave, v_pore_max
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase

  PetscReal, pointer :: porosity_loc_p(:)
  PetscReal :: dt_cfl_1
  PetscErrorCode :: ierr

  field => patch%field
  global_aux_vars => patch%aux%Global%aux_vars
  grid => patch%grid

  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)

  max_dt_cfl_1 = 1.d20
  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
      distance = cur_connection_set%dist(0,iconn)
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      do iphase = 1, option%nphase
        por_sat_min = min(porosity_loc_p(ghosted_id_up)* &
                          global_aux_vars(ghosted_id_up)%sat(iphase), &
                          porosity_loc_p(ghosted_id_dn)* &
                          global_aux_vars(ghosted_id_dn)%sat(iphase))
        por_sat_ave = (fraction_upwind*porosity_loc_p(ghosted_id_up)* &
                       global_aux_vars(ghosted_id_up)%sat(iphase) + &
                      (1.d0-fraction_upwind)*porosity_loc_p(ghosted_id_dn)* &
                      global_aux_vars(ghosted_id_dn)%sat(iphase))
        v_darcy = patch%internal_velocities(iphase,sum_connection)
        v_pore_max = v_darcy / por_sat_min
        v_pore_ave = v_darcy / por_sat_ave
        !geh: I use v_por_max to ensure that we limit the cfl based on the
        !     highest velocity through the face.  If porosity*saturation
        !     varies, the pore water velocity will be highest on the side
        !     of the face with the smalled value of porosity*saturation.
        dt_cfl_1 = distance / dabs(v_pore_max)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
      enddo
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
      if (patch%imat(ghosted_id_dn) <= 0) cycle
      !geh: since on boundary, dist must be scaled by 2.d0
      distance = 2.d0*cur_connection_set%dist(0,iconn)
      do iphase = 1, option%nphase
        por_sat_ave = porosity_loc_p(ghosted_id_dn)* &
                      global_aux_vars(ghosted_id_dn)%sat(iphase)
        v_darcy = patch%boundary_velocities(iphase,sum_connection)
        v_pore_ave = v_darcy / por_sat_ave
        dt_cfl_1 = distance / dabs(v_pore_ave)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)

end subroutine PatchCalculateCFL1Timestep

! ************************************************************************** !
!
! PatchGetVarNameFromKeyword: Returns the name of variable defined by keyword
! author: Glenn Hammond
! date: 07/28/11
!
! ************************************************************************** !
function PatchGetVarNameFromKeyword(keyword,option)
 
  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: PatchGetVarNameFromKeyword
  character(len=MAXSTRINGLENGTH) :: var_name

  select case(keyword)
    case('PROCESSOR_ID')
      var_name = 'Processor ID'
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call printErrMsg(option)
  end select

  PatchGetVarNameFromKeyword = var_name

end function PatchGetVarNameFromKeyword

! ************************************************************************** !
!
! PatchGetIvarsFromKeyword: Returns the ivar and isubvars for extracting
!                            datasets using PatchGet/PatchSet routines
! author: Glenn Hammond
! date: 07/28/11
!
! ************************************************************************** !
subroutine PatchGetIvarsFromKeyword(keyword,ivar,isubvar,var_type,option)
 
  use Option_module
  use Variables_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: var_type
  type(option_type) :: option

  select case(keyword)
    case('PROCESSOR_ID')
      ivar = PROCESSOR_ID
      isubvar = ZERO_INTEGER
      var_type = INT_VAR
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call printErrMsg(option)
  end select

end subroutine

! ************************************************************************** !
!
! PatchDestroyList: Deallocates a patch list and array of patches
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine PatchDestroyList(patch_list)

  implicit none
  
  type(patch_list_type), pointer :: patch_list
    
  type(patch_type), pointer :: cur_patch, prev_patch
  
  if (.not.associated(patch_list)) return
  
  if (associated(patch_list%array)) deallocate(patch_list%array)
  nullify(patch_list%array)
  
  cur_patch => patch_list%first
  do 
    if (.not.associated(cur_patch)) exit
    prev_patch => cur_patch
    cur_patch => cur_patch%next
    call PatchDestroy(prev_patch)
  enddo
  
  nullify(patch_list%first)
  nullify(patch_list%last)
  patch_list%num_patch_objects = 0
  
  deallocate(patch_list)
  nullify(patch_list)

end subroutine PatchDestroyList

! ************************************************************************** !
!
! PatchDestroy: Deallocates a patch object
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine PatchDestroy(patch)

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(patch_type), pointer :: patch
  
  call DeallocateArray(patch%imat)
  call DeallocateArray(patch%sat_func_id)
  call DeallocateArray(patch%internal_velocities)
  call DeallocateArray(patch%boundary_velocities)
  call DeallocateArray(patch%internal_fluxes)
  call DeallocateArray(patch%boundary_fluxes)
  call DeallocateArray(patch%internal_tran_coefs)
  call DeallocateArray(patch%boundary_tran_coefs)
  call DeallocateArray(patch%ss_fluid_fluxes)

  if (associated(patch%material_property_array)) &
    deallocate(patch%material_property_array)
  nullify(patch%material_property_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%material_properties)
  if (associated(patch%saturation_function_array)) &
    deallocate(patch%saturation_function_array)
  nullify(patch%saturation_function_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%saturation_functions)

#ifdef SURFACE_FLOW
  nullify(patch%surf_field)
  if (associated(patch%surf_material_property_array)) &
    deallocate(patch%surf_material_property_array)
  nullify(patch%surf_material_property_array)
  nullify(patch%surf_material_properties)
  if (associated(patch%surf_internal_fluxes)) deallocate(patch%surf_internal_fluxes)
  if (associated(patch%surf_boundary_fluxes)) deallocate(patch%surf_boundary_fluxes)
  nullify(patch%surf_internal_fluxes)
  nullify(patch%surf_boundary_fluxes)
#endif

  ! solely nullify grid since destroyed in discretization
  nullify(patch%grid)
  call RegionDestroyList(patch%regions)
  call CouplerDestroyList(patch%boundary_conditions)
  call CouplerDestroyList(patch%initial_conditions)
  call CouplerDestroyList(patch%source_sinks)
  
  
  call ObservationDestroyList(patch%observation)
  call StrataDestroyList(patch%strata)
  
  call AuxDestroy(patch%aux)
  
  call ObservationDestroyList(patch%observation)
  
  ! these are solely pointers, must not destroy.
  nullify(patch%reaction)
  nullify(patch%datasets)
  nullify(patch%field)
  
  deallocate(patch)
  nullify(patch)
  
end subroutine PatchDestroy

#ifdef SURFACE_FLOW
! ************************************************************************** !
!
! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchGetVariable2(patch,surf_field,option,output_option,vec,ivar, &
                           isubvar,isubvar1)

  use Grid_module
  use Option_module
  use Output_Aux_module
  use Surface_Field_module
  use Variables_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
  !type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(surface_field_type), pointer :: surf_field
  type(patch_type), pointer :: patch  
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass
  PetscReal :: tempreal
  PetscInt :: tempint
  PetscInt :: irate, istate, irxn
  PetscErrorCode :: ierr

  grid => patch%grid

  call VecGetArrayF90(vec,vec_ptr,ierr)
  
  iphase = 1
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%surf_aux%SurfaceGlobal%aux_vars(grid%nL2G(local_id))%head(1)
      enddo
    case(SURFACE_LIQUID_TEMPERATURE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%surf_aux%SurfaceGlobal%aux_vars(grid%nL2G(local_id))%temp(1)
      enddo
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
    case(PROCESSOR_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariable'')') ivar
      call printErrMsg(option)
  end select

end subroutine PatchGetVariable2

#endif
! SURFACE_FLOW

end module Patch_module
