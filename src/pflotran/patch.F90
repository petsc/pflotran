module Patch_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Observation_module
  use Strata_module
  use Region_module
  use Reaction_Aux_module
  use Material_module
  use Field_module
  
  use Auxilliary_module

  implicit none

  private

#include "definitions.h"

  type, public :: patch_type 
    
    PetscInt :: id
    
    ! These arrays will be used by all modes, mode-specific arrays should
    ! go in the auxilliary data stucture for that mode
    PetscInt, pointer :: imat(:)
      ! Integer array of material ids of size ngmax.
    type(material_property_ptr_type), pointer :: material_property_array(:)

#ifdef SUBCONTINUUM_MODEL
    !These arrays will hold the no. of subcontinuum types at each cell 
    PetscInt, pointer :: num_subcontinuum_type(:,:)
    
    !These arrays will hold the subcontinuum types ids
    PetscInt, pointer :: subcontinuum_type_ids(:)

    type(subcontinuum_property_ptr_type), pointer ::  &
                          subcontinuum_property_array(:)
    type(subcontinuum_field_typep), pointer :: subcontinuum_field_patch
#endif

    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    PetscReal, pointer :: internal_fluxes(:,:,:)    
    PetscReal, pointer :: boundary_fluxes(:,:,:)  
    PetscReal, pointer :: internal_tran_coefs(:,:)
    PetscReal, pointer :: boundary_tran_coefs(:,:)

    type(grid_type), pointer :: grid

    type(region_list_type), pointer :: regions

    type(coupler_list_type), pointer :: boundary_conditions
    type(coupler_list_type), pointer :: initial_conditions
    type(coupler_list_type), pointer :: source_sinks

    ! pointer to field object in mother realization object
    type(field_type), pointer :: field 
    type(strata_list_type), pointer :: strata
    type(observation_list_type), pointer :: observation

    type(reaction_type), pointer :: reaction
    type(auxilliary_type) :: aux
    
    type(patch_type), pointer :: next

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
    
  public :: PatchCreate, PatchDestroy, PatchCreateList, PatchDestroyList, &
            PatchAddToList, PatchConvertListToArray, PatchProcessCouplers, &
            PatchUpdateAllCouplerAuxVars, PatchInitAllCouplerAuxVars, &
            PatchLocalizeRegions, PatchUpdateUniformVelocity, &
            PatchGetDataset, PatchGetDatasetValueAtCell, &
            PatchSetDataset, &
            PatchInitConstraints, &
            PatchCountCells

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
  nullify(patch%imat)
  nullify(patch%material_property_array)
#ifdef SUBCONTINUUM_MODEL
  nullify(patch%subcontinuum_count)
  nullify(patch%subcontinnuum_ids)
  nullify(patch%subcontinuum_property_array)  
  nullify(patch%subcontinuum_field_patch)  
#endif
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)
  nullify(patch%internal_fluxes)
  nullify(patch%boundary_fluxes)
  nullify(patch%internal_tran_coefs)
  nullify(patch%boundary_tran_coefs)

  nullify(patch%grid)

  allocate(patch%regions)
  call RegionInitList(patch%regions)
  
  allocate(patch%boundary_conditions)
  call CouplerInitList(patch%boundary_conditions)
  allocate(patch%initial_conditions)
  call CouplerInitList(patch%initial_conditions)
  allocate(patch%source_sinks)
  call CouplerInitList(patch%source_sinks)

  nullify(patch%field)
  allocate(patch%observation)
  call ObservationInitList(patch%observation)

  allocate(patch%strata)
  call StrataInitList(patch%strata)
  
  call AuxInit(patch%aux)
  
  nullify(patch%reaction)
  
  nullify(patch%next)
  
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
  
  call GridLocalizeRegions(patch%grid,patch%regions,option)
 
end subroutine PatchLocalizeRegions

! ************************************************************************** !
!
! PatchProcessCouplers: Assigns conditions and regions to couplers
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
!subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
!                                material_properties, subcontinuum_properties, option)
subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
                                material_properties,option)

  use Option_module
  use Material_module
  use Condition_module
  use Connection_module

#ifdef SUBCONTINUUM_MODEL
  use Subcontinuum_module
#endif

  implicit none
  
  type(patch_type) :: patch
  type(material_property_type), pointer :: material_properties
#ifdef SUBCONTINUUM_MODEL
  type(subcontinuum_property_type), pointer :: subcontinuum_properties
#endif
  type(condition_list_type) :: flow_conditions
  type(tran_condition_list_type) :: transport_conditions
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list 
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation, next_observation
  
  PetscInt :: temp_int, isub
  
  call MaterialPropConvertListToArray(material_properties,patch%material_property_array,option)
#ifdef SUBCONTINUUM_MODEL
  call SubcontinuumPropConvertListToArray(subcontinuum_properties,patch%subcontinuum_property_array,option)
#endif
  
  ! boundary conditions
  coupler => patch%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%regions)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region ' // trim(coupler%region_name) // &
                 '" in boundary condition ' // &
                 trim(coupler%name) // &
                 ' not found in region list'
      call printErrMsg(option)
    endif
    if (associated(patch%grid%structured_grid)) then
      if (coupler%region%num_cells > 0 .and. &
          (coupler%region%iface == 0 .and. &
           .not.associated(coupler%region%faces))) then
        option%io_buffer = 'Region ' // trim(coupler%region_name) // &
                 ', which is tied to a boundary condition, has not ' // &
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
                   '" in boundary condition ' // &
                   trim(coupler%name) // &
                   ' not found in flow condition list'
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
                   trim(coupler%flow_condition_name) // &
                   '" in boundary condition ' // &
                   trim(coupler%name) // &
                   ' not found in transport condition list'
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
      option%io_buffer = 'Region ' // trim(coupler%region_name) // &
                 '" in initial condition ' // &
                 trim(coupler%name) // &
                 ' not found in region list'
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
                   '" in initial condition ' // &
                   trim(coupler%name) // &
                   ' not found in flow condition list'
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
                   '" in initial condition ' // &
                   trim(coupler%name) // &
                   ' not found in transport condition list'
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
      option%io_buffer = 'Region ' // trim(coupler%region_name) // &
                 '" in source/sink ' // &
                 trim(coupler%name) // &
                 ' not found in region list'
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
                   '" in source/sink ' // &
                   trim(coupler%name) // &
                   ' not found in flow condition list'
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
                   '" in source/sink ' // &
                   trim(coupler%name) // &
                   ' not found in transport condition list'
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
      option%io_buffer = 'Region ' // trim(coupler%region_name) // &
                 '" in strata not found in region list'
        call printErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material_property => &
          MaterialPropGetPtrFromArray(strata%material_property_name, &
                                      patch%material_property_array)
        if (.not.associated(strata%material_property)) then
          option%io_buffer = 'Material ' // &
                             trim(strata%material_property_name) // &
                             ' not found in material list'
          call printErrMsg(option)
        endif

#ifdef SUBCONTINUUM_PROPERTY
        ! connect subcontinuum properties pointers
        ! allocate storage to hold subcontinuum pointers
        if (strata%material_property%num_subcontinuum_type > 0) then
          allocate(strata%subcontinuum_property( & 
                      strata%material_property%subcontinuum_type_count))
          ! loop over each subcontinuum
          do isub=1,strata%material_property%num_subcontinuum_type
            strata%subcontinuum_property(isub) => &
              SubcontinuumPropGetPtrFromArray( & 
               strata%material_property%subcontinuum_property_name(isub), &
               patch%subcontinuum_property_array)
            if (.not.associated(strata%subcontinuum_property(isub))) then
              option%io_buffer = 'Subcontinuum ' // &
                trim(strata%material_property%subcontinuum_property_name(isub)) // &
                             ' not found in subcontinuum list'
              call printErrMsg(option)
            endif
          enddo
        endif
#endif

      endif
    else
      nullify(strata%region)
      nullify(strata%material_property)
#ifdef SUBCONTINUUM_MODEL
      nullify(strata%subcontinuum_property)
#endif
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
                 '" in observation point ' // &
                 trim(observation%name) // &
                 ' not found in region list'                   
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
          option%io_buffer = 'Boundary Condition ' // &
                   trim(observation%linkage_name) // &
                   ' not found in Boundary Condition list'
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
 
  if (patch%grid%itype == STRUCTURED_GRID_MIMETIC) then
!#ifdef DASVYAT
!    temp_int = ConnectionGetNumberInList(patch%grid%boundary_connection_set_list)
    temp_int = CouplerGetNumBoundConnectionsInListMFD(patch%grid, &
                            patch%boundary_conditions, &
                           option)
!#endif
  else  
    temp_int = CouplerGetNumConnectionsInList(patch%boundary_conditions)
  end if
  allocate(patch%boundary_velocities(option%nphase,temp_int)) 
  patch%boundary_velocities = 0.d0
  allocate(patch%boundary_tran_coefs(option%nphase,temp_int))
  patch%boundary_tran_coefs = 0.d0
  if (option%store_solute_fluxes) then
    allocate(patch%boundary_fluxes(option%nphase,option%ntrandof,temp_int))
    patch%boundary_fluxes = 0.d0
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
subroutine PatchInitAllCouplerAuxVars(patch,reaction,option)

  use Option_module
  use Reaction_Aux_module
  
  implicit none
  
  type(patch_type) :: patch
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  
  PetscBool :: force_update_flag = PETSC_TRUE
  
  call PatchInitCouplerAuxVars(patch%initial_conditions,reaction, &
                               option)
  call PatchInitCouplerAuxVars(patch%boundary_conditions,reaction, &
                               option)
  call PatchInitCouplerAuxVars(patch%source_sinks,reaction, &
                               option)
  
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
subroutine PatchInitCouplerAuxVars(coupler_list,reaction,option)

  use Option_module
  use Connection_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Condition_module
  
  implicit none
  
  type(coupler_list_type), pointer :: coupler_list
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  
  PetscInt :: num_connections
  PetscBool :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  
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

        if (associated(coupler%flow_condition%pressure)) then

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

            case(THC_MODE)
              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0

            case(MPH_MODE, IMS_MODE, FLASH2_MODE)
!geh              allocate(coupler%flow_aux_real_var(option%nflowdof*option%nphase,num_connections))
              allocate(coupler%flow_aux_real_var(option%nflowdof,num_connections))
              allocate(coupler%flow_aux_int_var(1,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
                
            case(G_MODE)
              allocate(coupler%flow_aux_real_var(FOUR_INTEGER,num_connections))
              allocate(coupler%flow_aux_int_var(0,num_connections))
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_int_var = 0
                
            case default
          end select
      
        endif ! associated(coupler%flow_condition%pressure)
      
      else if (coupler%itype == SRC_SINK_COUPLER_TYPE) then

        if (associated(coupler%flow_condition%rate)) then

          select case(coupler%flow_condition%rate%itype)
            case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)

              select case(option%iflowmode)
                case(RICHARDS_MODE)
                  allocate(coupler%flow_aux_real_var(1,num_connections))
                  coupler%flow_aux_real_var = 0.d0
                  
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
        allocate(cur_constraint_coupler%global_auxvar)
        allocate(cur_constraint_coupler%rt_auxvar)
        option%iflag = 0 ! be sure not to allocate mass_balance array
        call GlobalAuxVarInit(cur_constraint_coupler%global_auxvar,option)
        call RTAuxVarInit(cur_constraint_coupler%rt_auxvar,reaction,option)
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
      
    coupler => coupler%next
  enddo
  
end subroutine PatchInitCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateAllCouplerAuxVars: Updates auxilliary variables associated 
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
  
  call PatchUpdateCouplerAuxVars(patch,patch%initial_conditions, &
                                 force_update_flag,option)

  call PatchUpdateCouplerAuxVars(patch,patch%boundary_conditions, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%source_sinks, &
                                 force_update_flag,option)

!  stop
end subroutine PatchUpdateAllCouplerAuxVars

! ************************************************************************** !
!
! PatchUpdateCouplerAuxVars: Updates auxilliary variables associated 
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

  implicit none
  
  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  PetscBool :: force_update_flag
  type(option_type) :: option
  
  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  PetscBool :: update
  
  PetscInt :: idof, num_connections
  
  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first
  
  do
    if (.not.associated(coupler)) exit
    
    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then


      num_connections = coupler%connection_set%num_connections
#ifdef DASVYAT      
      if (option%mimetic) then
        num_connections = coupler%numfaces_set
      end if
#endif

      flow_condition => coupler%flow_condition

      if (force_update_flag .or. FlowConditionIsTransient(flow_condition)) then
        if (associated(flow_condition%pressure)) then
          select case(flow_condition%pressure%itype)
            case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
!              do idof = 1, condition%num_sub_conditions
!                if (associated(condition%sub_condition_ptr(idof)%ptr)) then
!                  coupler%flow_aux_real_var(idof,1:num_connections) = &
!                    condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
!                endif
!              enddo
              select case(option%iflowmode)
                case(FLASH2_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)  ! <-- Chuan Fix
                   coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(MPH_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)  ! <-- Chuan Fix
                   coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(IMS_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)  ! <-- Chuan Fix
                   coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)! <-- Chuan Fix
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(THC_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)
                  coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
                    flow_condition%temperature%dataset%cur_value(1)
                  coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                    flow_condition%concentration%dataset%cur_value(1)
                  coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                    flow_condition%iphase
                case(RICHARDS_MODE)
                  coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                    flow_condition%pressure%dataset%cur_value(1)
                case(G_MODE)
                  do idof = 1, option%nflowdof
                    coupler%flow_aux_real_var(idof,1:num_connections) = &
                      flow_condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
                  enddo
              end select
            case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
              call HydrostaticUpdateCoupler(coupler,option,patch%grid)
          end select
        endif
        if (associated(flow_condition%rate)) then
          select case(flow_condition%rate%itype)
            case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
              call PatchScaleSourceSink(patch,option)
          end select
        endif
      endif
     
    endif
      
    ! TRANSPORT
    ! nothing for transport at this point in time

    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !
!
! PatchScaleSourceSink: Scales select source/sinks based on perms*volume
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine PatchScaleSourceSink(patch,option)

  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(patch_type) :: patch
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(field_type), pointer :: field
  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: perm_loc_ptr(:)
  PetscReal, pointer :: vol_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount, x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  PetscInt :: ghosted_neighbors(27)
    
  field => patch%field
  grid => patch%grid

  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_loc_ptr,ierr)
  call GridVecGetArrayF90(grid,field%volume,vol_ptr,ierr)

  grid => patch%grid

  cur_source_sink => patch%source_sinks%first
  do
    if (.not.associated(cur_source_sink)) exit

    call VecZeroEntries(field%work,ierr)
    call GridVecGetArrayF90(grid,field%work,vec_ptr,ierr)

    cur_connection_set => cur_source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      select case(option%iflowmode)
        case(RICHARDS_MODE,G_MODE)
           call GridGetGhostedNeighbors(grid,ghosted_id,STAR_STENCIL, &
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
        case(THC_MODE)
        case(MPH_MODE)
        case(IMS_MODE)
        case(FLASH2_MODE)
      end select 

    enddo
    
    call GridVecRestoreArrayF90(grid,field%work,vec_ptr,ierr)
    call VecNorm(field%work,NORM_1,scale,ierr)
    scale = 1.d0/scale
    call VecScale(field%work,scale,ierr)

    call GridVecGetArrayF90(grid,field%work,vec_ptr, ierr)
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      select case(option%iflowmode)
        case(RICHARDS_MODE,G_MODE)
          cur_source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
            vec_ptr(local_id)
        case(THC_MODE)
        case(MPH_MODE)
        case(IMS_MODE)
        case(FLASH2_MODE)
      end select 

    enddo
    call GridVecRestoreArrayF90(grid,field%work,vec_ptr,ierr)
    
    cur_source_sink => cur_source_sink%next
  enddo

  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_loc_ptr, ierr)
  call GridVecRestoreArrayF90(grid,field%volume,vol_ptr, ierr)
   
end subroutine PatchScaleSourceSink

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
  use Condition_module
  use water_eos_module
    
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
          global_auxvar%pres = &
            cur_coupler%flow_condition%pressure%dataset%cur_value(1)
        else
          global_auxvar%pres = option%reference_pressure
        endif
        if (associated(cur_coupler%flow_condition%temperature)) then
          global_auxvar%temp = &
            cur_coupler%flow_condition%temperature%dataset%cur_value(1)
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
                            cur_constraint_coupler%surface_complexes, &
                            cur_constraint_coupler%colloids, &
                            cur_constraint_coupler%num_iterations, &
                            PETSC_TRUE,option)
      ! turn on flag indicating constraint has not yet been used

      cur_constraint_coupler%iflag = ONE_INTEGER
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
! PatchGetDataset: Extracts variables indexed by ivar and isubvar from a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchGetDataset(patch,field,option,output_option,vec,ivar, &
  isubvar,isubvar1)

  use Grid_module
  use Option_module
  use Field_module
  
  use Mphase_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module  
  use Reaction_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(option_type), pointer :: option
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
  PetscReal :: xmass
  PetscInt :: irate, istate
  PetscErrorCode :: ierr

  grid => patch%grid

  call GridVecGetArrayF90(grid,vec,vec_ptr,ierr)

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,SC_FUGA_COEFF,STATE)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%temp
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%pres
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_kg
            enddo
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = 0.d0
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
          case(PRESSURE)
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
          case(PRESSURE)
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
          case(PRESSURE)
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
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%xmol(2+isubvar)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
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
      else if (associated(patch%aux%immis)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%temp(1)
            enddo
          case(PRESSURE)
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
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(2)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%den_kg(2)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Immis%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%u(1)
            enddo
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%temp
            enddo
          case(PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              istate = patch%aux%Global%aux_vars(ghosted_id)%istate
              select case(istate)
                case(LIQUID_STATE,TWO_PHASE_STATE)
                  vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%liquid_phase)
                case(GAS_STATE)
                  vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%gas_phase)
              end select
            enddo
          case(STATE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%aux_vars(grid%nL2G(local_id))%istate
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%den_kg(option%liquid_phase)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%U(option%liquid_phase)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%xmol(isubvar,option%liquid_phase)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%sat(option%gas_phase)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%U(option%gas_phase)
            enddo
          case(GAS_DENSITY) 
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%den_kg(option%gas_phase)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%aux_vars(ZERO_INTEGER,grid%nL2G(local_id))%xmol(isubvar,option%gas_phase)
            enddo
        end select         
      endif
      
    case(PH,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY,TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_RATE,MINERAL_VOLUME_FRACTION,SURFACE_CMPLX,SURFACE_CMPLX_FREE, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE)
         
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
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%total(isubvar,iphase)
          enddo
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_volfrac(isubvar)
          enddo
        case(MINERAL_RATE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%mnrl_rate(isubvar)
          enddo
        case(SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%eqsrfcplx_conc(isubvar)
          enddo
        case(SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%eqsrfcplx_free_site_conc(isubvar)
          enddo
        case(KIN_SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%kinsrfcplx_conc(isubvar)
          enddo
        case(KIN_SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%aux_vars(grid%nL2G(local_id))%kinsrfcplx_free_site_conc(isubvar)
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
          call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call ReactionComputeKd(isubvar,vec_ptr(local_id), &
                                   patch%aux%RT%aux_vars(ghosted_id), &
                                   patch%aux%Global%aux_vars(ghosted_id), &
                                   vec_ptr2(ghosted_id),patch%reaction,option)
          enddo
          call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
        case(TOTAL_SORBED)
          if (patch%reaction%neqsorb > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (patch%reaction%kinmr_nrate > 0) then
                vec_ptr(local_id) = 0.d0
                do irate = 1, patch%reaction%kinmr_nrate
                  vec_ptr(local_id) = vec_ptr(local_id) + &
                    patch%aux%RT%aux_vars(ghosted_id)%kinmr_total_sorb(isubvar,irate)
                enddo            
              else
                vec_ptr(local_id) = patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
              endif
            enddo
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%neqsorb > 0 .and. patch%reaction%ncollcomp > 0) then
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
      call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
    case(PHASE)
      call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%imat(grid%nL2G(local_id))
      enddo
    case default
      call printErrMsg(option,'IVAR not found in OutputGetVarFromArray')
  end select

  call GridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)
  
end subroutine PatchGetDataset

! ************************************************************************** !
!
! PatchGetDatasetValueAtCell: Returns variables indexed by ivar,
!                             isubvar, local id from Reactive Transport type
! author: Glenn Hammond
! date: 02/11/08
!
! ************************************************************************** !
function PatchGetDatasetValueAtCell(patch,field,option,output_option, &
                                    ivar,isubvar,ghosted_id,isubvar1)

  use Grid_module
  use Option_module
  use Field_module

  use Mphase_Aux_module
  use THC_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module  
  use Reaction_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscReal :: PatchGetDatasetValueAtCell
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  PetscInt :: iphase
  PetscInt :: ghosted_id

  PetscReal :: value, xmass
  PetscInt :: irate, istate
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr2(:)  
  PetscErrorCode :: ierr

  grid => patch%grid
  
  value = -999.99d0

  ! inactive grid cell
  if (patch%imat(ghosted_id) <= 0) then
    PatchGetDatasetValueAtCell = 0.d0
    return
  endif

  iphase = 1
  xmass = 1.d0
  if (associated(patch%aux%Global%aux_vars(ghosted_id)%xmass)) &
    xmass = patch%aux%Global%aux_vars(ghosted_id)%xmass(iphase)
             
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,SC_FUGA_COEFF,STATE)
         
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%THC%aux_vars(ghosted_id)%temp
          case(PRESSURE)
            value = patch%aux%THC%aux_vars(ghosted_id)%pres
          case(LIQUID_SATURATION)
            value = patch%aux%THC%aux_vars(ghosted_id)%sat
          case(LIQUID_DENSITY)
            value = patch%aux%THC%aux_vars(ghosted_id)%den_kg
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            value = 0.d0
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%THC%aux_vars(ghosted_id)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%THC%aux_vars(ghosted_id)%u
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
          case(PRESSURE)
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
          case(PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
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
          case(PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(2)
          case(GAS_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(2+isubvar)
          case(GAS_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(2)
          case(GAS_DENSITY_MOL) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den(2)
          case(SC_FUGA_COEFF)
            value = patch%aux%Global%aux_vars(ghosted_id)%fugacoeff(1)   
          case(LIQUID_MOLE_FRACTION)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%xmol(isubvar)
          case(LIQUID_ENERGY)
            value = patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
        end select
      else if (associated(patch%aux%Immis)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%temp(1)
          case(PRESSURE)
            value = patch%aux%Global%aux_vars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%aux_vars(ghosted_id)%sat(2)
          case(GAS_ENERGY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(2)
          case(GAS_DENSITY) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den_kg(2)
          case(GAS_DENSITY_MOL) 
            value = patch%aux%Global%aux_vars(ghosted_id)%den(2)
          case(LIQUID_ENERGY)
            value = patch%aux%Immis%aux_vars(ghosted_id)%aux_var_elem(0)%u(1)
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%temp
          case(PRESSURE)
            istate = patch%aux%Global%aux_vars(ghosted_id)%istate
            select case(istate)
              case(LIQUID_STATE,TWO_PHASE_STATE)
                value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%liquid_phase)
              case(GAS_STATE)
                value = patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%gas_phase)
            end select
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
      
    case(PH,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY, &
         TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_VOLUME_FRACTION,MINERAL_RATE, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE)
         
      select case(ivar)
        case(PH)
          value = -log10(patch%aux%RT%aux_vars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%aux_vars(ghosted_id)%pri_molal(isubvar))
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
        case(MINERAL_VOLUME_FRACTION)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_volfrac(isubvar)
        case(MINERAL_RATE)
          value = patch%aux%RT%aux_vars(ghosted_id)%mnrl_rate(isubvar)
        case(SURFACE_CMPLX)
          value = patch%aux%RT%aux_vars(ghosted_id)%eqsrfcplx_conc(isubvar)
        case(SURFACE_CMPLX_FREE)
          value = patch%aux%RT%aux_vars(ghosted_id)%eqsrfcplx_free_site_conc(isubvar)
        case(KIN_SURFACE_CMPLX)
          value = patch%aux%RT%aux_vars(ghosted_id)%kinsrfcplx_conc(isubvar)
        case(KIN_SURFACE_CMPLX_FREE)
          value = patch%aux%RT%aux_vars(ghosted_id)%kinsrfcplx_free_site_conc(isubvar)
        case(PRIMARY_ACTIVITY_COEF)
          value = patch%aux%RT%aux_vars(ghosted_id)%pri_act_coef(isubvar)
        case(SECONDARY_ACTIVITY_COEF)
          value = patch%aux%RT%aux_vars(ghosted_id)%sec_act_coef(isubvar)
        case(PRIMARY_KD)
          call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
          call ReactionComputeKd(isubvar,value, &
                                 patch%aux%RT%aux_vars(ghosted_id), &
                                 patch%aux%Global%aux_vars(ghosted_id), &
                                 vec_ptr2(ghosted_id),patch%reaction,option)
          call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
        case(TOTAL_SORBED)
          if (patch%reaction%neqsorb > 0) then
            if (patch%reaction%kinmr_nrate > 0) then
              value = 0.d0
              do irate = 1, patch%reaction%kinmr_nrate
                value = value + &
                  patch%aux%RT%aux_vars(ghosted_id)%kinmr_total_sorb(isubvar,irate)
              enddo            
            else
              value = patch%aux%RT%aux_vars(ghosted_id)%total_sorb_eq(isubvar)
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%neqsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            value = patch%aux%RT%aux_vars(ghosted_id)%colloid%total_eq_mob(isubvar)
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
      call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
    case(PHASE)
      call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      value = vec_ptr2(ghosted_id)
      call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
    case(MATERIAL_ID)
      value = patch%imat(ghosted_id)
  end select

  PatchGetDatasetValueAtCell = value
 
end function PatchGetDatasetValueAtCell

! ************************************************************************** !
!
! PatchSetDataset: Sets variables indexed by ivar and isubvar within a patch
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine PatchSetDataset(patch,field,option,vec,vec_format,ivar,isubvar)

  use Grid_module
  use Option_module
  use Field_module

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

  call GridVecGetArrayF90(grid,vec,vec_ptr,ierr)

  if (vec_format == NATURAL) then
    call printErrMsg(option,&
                     'NATURAL vector format not supported by PatchSetDataset')
  endif

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,PRESSURE,LIQUID_SATURATION,GAS_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,STATE)
      if (associated(patch%aux%THC)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%sat = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%sat = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%THC%aux_vars(grid%nL2G(local_id))%den_kg = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%THC%aux_vars(ghosted_id)%den_kg = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_SATURATION,GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
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
          case(PRESSURE)
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
          case(PRESSURE)
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
          case(PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%aux_vars(grid%nL2G(local_id))%aux_var_elem(0)%pres = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%aux_vars(ghosted_id)%aux_var_elem(0)%pres = vec_ptr(ghosted_id)
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
          case(PRESSURE)
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
          case(PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              istate = patch%aux%Global%aux_vars(ghosted_id)%istate
              select case(istate)
                case(LIQUID_STATE,TWO_PHASE_STATE)
                  patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%liquid_phase) = vec_ptr(local_id)
                case(GAS_STATE)
                  patch%aux%General%aux_vars(ZERO_INTEGER,ghosted_id)%pres(option%gas_phase) = vec_ptr(local_id)
              end select
            enddo
          case(STATE)
            do local_id=1,grid%nlmax
              patch%aux%Global%aux_vars(ghosted_id)%istate = int(vec_ptr(local_id)+1.d-10)
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
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF)
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
        call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
        do local_id=1,grid%nlmax
          vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
      else if (vec_format == LOCAL) then
        call GridVecGetArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
        vec_ptr2(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
        call GridVecRestoreArrayF90(grid,field%porosity_loc,vec_ptr2,ierr)
      endif
    case(PHASE)
      if (vec_format == GLOBAL) then
        call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
        do local_id=1,grid%nlmax
          vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
        call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      else if (vec_format == LOCAL) then
        call GridVecGetArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
        vec_ptr2(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
        call GridVecRestoreArrayF90(grid,field%iphas_loc,vec_ptr2,ierr)
      endif
    case(MATERIAL_ID)
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          patch%imat(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
      else if (vec_format == LOCAL) then
        patch%imat(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
      endif
  end select

  call GridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)
  
end subroutine PatchSetDataset

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

  implicit none
  
  type(patch_type), pointer :: patch
  
  if (associated(patch%imat)) deallocate(patch%imat)
  nullify(patch%imat)
  if (associated(patch%internal_velocities)) deallocate(patch%internal_velocities)
  nullify(patch%internal_velocities)
  if (associated(patch%boundary_velocities)) deallocate(patch%boundary_velocities)
  nullify(patch%boundary_velocities)
  if (associated(patch%internal_fluxes)) deallocate(patch%internal_fluxes)
  nullify(patch%internal_fluxes)
  if (associated(patch%boundary_fluxes)) deallocate(patch%boundary_fluxes)
  nullify(patch%boundary_fluxes)
  if (associated(patch%internal_tran_coefs)) deallocate(patch%internal_tran_coefs)
  nullify(patch%internal_tran_coefs)
  if (associated(patch%boundary_tran_coefs)) deallocate(patch%boundary_tran_coefs)
  nullify(patch%boundary_tran_coefs)

  call GridDestroy(patch%grid)
  call RegionDestroyList(patch%regions)
  call CouplerDestroyList(patch%boundary_conditions)
  call CouplerDestroyList(patch%initial_conditions)
  call CouplerDestroyList(patch%source_sinks)
  
  nullify(patch%field)
  
  call ObservationDestroyList(patch%observation)
  call StrataDestroyList(patch%strata)
  
  call AuxDestroy(patch%aux)
  
  call ObservationDestroyList(patch%observation)
  
  nullify(patch%reaction)
  
  deallocate(patch)
  nullify(patch)
  
end subroutine PatchDestroy

end module Patch_module
