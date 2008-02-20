module Realization_module

  use Grid_module
  use Option_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Material_module
  use Strata_module
  use Breakthrough_module
  use Field_module
  use Debug_module
  use Waypoint_module
  
  use Reactive_Transport_Aux_module
  
  implicit none

private

#include "definitions.h"

  type, public :: realization_type

    type(grid_type), pointer :: grid
    type(option_type), pointer :: option
    type(field_type), pointer :: field
    type(pflow_debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    type(region_list_type), pointer :: regions
    
    type(condition_list_type), pointer :: flow_conditions
    type(coupler_list_type), pointer :: flow_boundary_conditions
    type(coupler_list_type), pointer :: flow_initial_conditions
    type(coupler_list_type), pointer :: flow_source_sinks
    
    type(condition_list_type), pointer :: transport_conditions
    type(coupler_list_type), pointer :: transport_boundary_conditions
    type(coupler_list_type), pointer :: transport_initial_conditions
    type(coupler_list_type), pointer :: transport_source_sinks
    
    type(breakthrough_list_type), pointer :: breakthrough
    type(strata_list_type), pointer :: strata
    
    type(material_type), pointer :: materials
    type(material_ptr_type), pointer :: material_array(:)
    type(thermal_property_type), pointer :: thermal_properties
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    
    type(reactive_transport_aux_type), pointer :: RTaux
    
    type(waypoint_list_type), pointer :: waypoints
    
  end type realization_type

  public :: RealizationCreate, RealizationDestroy, &
            RealizationProcessCouplers, &
            RealizationInitCouplerAuxVars, &
            RealizationUpdate, RealizationAddWaypointsToList
  
contains
  
! ************************************************************************** !
!
! RealizationCreate: Allocates and initializes a new Realization object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function RealizationCreate()

  implicit none
  
  type(realization_type), pointer :: RealizationCreate
  
  type(realization_type), pointer :: realization
  
  allocate(realization)
  realization%option => OptionCreate()
  realization%field => FieldCreate()
  realization%debug => DebugCreatePflow()
  realization%output_option => OutputOptionCreate()
  nullify(realization%grid)
  allocate(realization%regions)
  call RegionInitList(realization%regions)

  allocate(realization%flow_conditions)
  call ConditionInitList(realization%flow_conditions)
  allocate(realization%flow_boundary_conditions)
  call CouplerInitList(realization%flow_boundary_conditions)
  allocate(realization%flow_initial_conditions)
  call CouplerInitList(realization%flow_initial_conditions)
  allocate(realization%flow_source_sinks)
  call CouplerInitList(realization%flow_source_sinks)

  allocate(realization%transport_conditions)
  call ConditionInitList(realization%transport_conditions)
  allocate(realization%transport_boundary_conditions)
  call CouplerInitList(realization%transport_boundary_conditions)
  allocate(realization%transport_initial_conditions)
  call CouplerInitList(realization%transport_initial_conditions)
  allocate(realization%transport_source_sinks)
  call CouplerInitList(realization%transport_source_sinks)

  allocate(realization%strata)
  call StrataInitList(realization%strata)
  allocate(realization%breakthrough)
  call BreakthroughInitList(realization%breakthrough)
  
  nullify(realization%materials)
  nullify(realization%material_array)
  nullify(realization%thermal_properties)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  nullify(realization%RTaux)
  
  RealizationCreate => realization
  
end function RealizationCreate  

! ************************************************************************** !
!
! RealizationProcessCouplers: Deallocates a realization
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RealizationProcessCouplers(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list 
  type(strata_type), pointer :: strata
  type(breakthrough_type), pointer :: breakthrough
  
  ! boundary conditions
  coupler => realization%flow_boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           realization%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in boundary condition list'
      call printErrMsg(realization%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 realization%flow_conditions)
    if (.not.associated(coupler%condition)) then
      coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                   realization%transport_conditions)
    endif
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in boundary condition list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => realization%flow_initial_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           realization%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in initial condition list'
      call printErrMsg(realization%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 realization%flow_conditions)
    if (.not.associated(coupler%condition)) then
      coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                   realization%transport_conditions)
    endif
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in initial condition list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => realization%flow_source_sinks%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           realization%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in source/sink list'
      call printErrMsg(realization%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 realization%flow_conditions)
    if (.not.associated(coupler%condition)) then
      coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                   realization%transport_conditions)
    endif
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in source/sink list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo
  
! Initially, all couplers are in flow lists.  Need to separate
! them into flow and transport lists. 
  call CouplerListSplitFlowAndTran(realization%flow_boundary_conditions, &
                                   realization%transport_boundary_conditions)
  call CouplerListSplitFlowAndTran(realization%flow_initial_conditions, &
                                   realization%transport_initial_conditions)
  call CouplerListSplitFlowAndTran(realization%flow_source_sinks, &
                                   realization%transport_source_sinks)

  if (realization%option%nflowdof == 0) then
    call CouplerDestroyList(realization%flow_boundary_conditions)
    call CouplerDestroyList(realization%flow_initial_conditions)
    call CouplerDestroyList(realization%flow_source_sinks)
  endif
  if (realization%option%ntrandof == 0) then
    call CouplerDestroyList(realization%transport_boundary_conditions)
    call CouplerDestroyList(realization%transport_initial_conditions)
    call CouplerDestroyList(realization%transport_source_sinks)
  endif

!----------------------------  
! AUX  
    
  ! strata
  ! connect pointers from strata to regions
  strata => realization%strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 1) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  realization%regions)
      if (.not.associated(strata%region)) then
        string = 'Region ' // trim(strata%region_name) // &
                 ' not found in region list'
        call printErrMsg(realization%option,string)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material => MaterialGetPtrFromList(strata%material_name, &
                                                  realization%materials)
        if (.not.associated(strata%material)) then
          string = 'Material ' // trim(strata%material_name) // &
                   ' not found in material list'
          call printErrMsg(realization%option,string)
        endif
      endif
    else
      nullify(strata%region)
      nullify(strata%material)
    endif
    strata => strata%next
  enddo 

  ! breakthrough
  breakthrough => realization%breakthrough%first
  do
    if (.not.associated(breakthrough)) exit
    ! pointer to region
    breakthrough%region => RegionGetPtrFromList(breakthrough%region_name, &
                                                realization%regions)
    if (.not.associated(breakthrough%region)) then
      string = 'Region ' // trim(breakthrough%region_name) // &
               ' not found in region list'
      call printErrMsg(realization%option,string)
    endif
    breakthrough => breakthrough%next
  enddo
 
end subroutine RealizationProcessCouplers

! ************************************************************************** !
!
! RealizationInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine RealizationInitCouplerAuxVars(realization,coupler_list)

  use Connection_module
  
  implicit none
  
  type(realization_type) :: realization
  type(coupler_list_type), pointer :: coupler_list
  
  PetscInt :: num_connections
  logical :: force_update_flag
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: coupler
    
  option => realization%option
  
  if (.not.associated(coupler_list)) return
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%connection)) then
      num_connections = coupler%connection%num_connections

      if (coupler%condition%iclass == FLOW_CLASS) then

        ! allocate arrays that match the number of connections
        select case(option%iflowmode)

          case(RICHARDS_MODE,RICHARDS_LITE_MODE)
         
            allocate(coupler%aux_real_var(option%nflowdof*option%nphase,num_connections))
            allocate(coupler%aux_int_var(1,num_connections))
            coupler%aux_real_var = 0.d0
            coupler%aux_int_var = 0

          case(MPH_MODE)

            allocate(coupler%aux_real_var(option%nflowdof*option%nphase,num_connections))
            allocate(coupler%aux_int_var(1,num_connections))
            coupler%aux_real_var = 0.d0
            coupler%aux_int_var = 0
              
          case default
        end select
         
      else ! TRANSPORT_CLASS

        allocate(coupler%aux_real_var(option%ntrandof,num_connections))
        coupler%aux_real_var = 0.d0

      endif
      
    endif
    coupler => coupler%next
  enddo
  
  force_update_flag = .true.
  call RealizationUpdateCouplerAuxVars(realization,coupler_list, &
                                       force_update_flag)

end subroutine RealizationInitCouplerAuxVars

! ************************************************************************** !
!
! RealizationUpdateCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in list
! author: Glenn Hammond
! date: 11/26/07
!
! ************************************************************************** !
subroutine RealizationUpdateCouplerAuxVars(realization,coupler_list, &
                                           force_update_flag)

  use Hydrostatic_module

  implicit none
  
  type(realization_type) :: realization
  type(coupler_list_type), pointer :: coupler_list
  logical :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  type(condition_type), pointer :: condition
  logical :: update
  
  PetscInt :: idof, num_connections
  
  if (.not.associated(coupler_list)) return
 
  coupler => coupler_list%first
  
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%aux_real_var)) then
        
      num_connections = coupler%connection%num_connections

      condition => coupler%condition
      
      if (condition%iclass == FLOW_CLASS) then

        update = .false.
        select case(realization%option%iflowmode)
          case(RICHARDS_MODE,MPH_MODE)
            if (force_update_flag .or. &
                condition%pressure%dataset%is_transient .or. &
                condition%pressure%gradient%is_transient .or. &
                condition%pressure%datum%is_transient .or. &
                condition%temperature%dataset%is_transient .or. &
                condition%temperature%gradient%is_transient .or. &
                condition%temperature%datum%is_transient .or. &
                condition%concentration%dataset%is_transient .or. &
                condition%concentration%gradient%is_transient .or. &
                condition%concentration%datum%is_transient) then
              update = .true.
            endif
          case(RICHARDS_LITE_MODE)
            if (force_update_flag .or. &
                condition%pressure%dataset%is_transient .or. &
                condition%pressure%gradient%is_transient .or. &
                condition%pressure%datum%is_transient) then
              update = .true.
            endif
        end select
        
        if (update) then
          select case(condition%pressure%itype)
            case(DIRICHLET_BC,NEUMANN_BC,MASS_RATE,ZERO_GRADIENT_BC)
              do idof = 1, condition%num_sub_conditions
                coupler%aux_real_var(idof,1:num_connections) = &
                  condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
              enddo
              coupler%aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                condition%iphase
            case(HYDROSTATIC_BC,SEEPAGE_BC)
    !          call HydrostaticUpdateCoupler(coupler,realization%option,realization%grid)
              call HydrostaticUpdateCouplerBetter(coupler,realization%option,realization%grid)
          end select
        endif
        
      else ! TRANSPORT_CLASS

        update = .false.
        if (force_update_flag .or. &
            condition%concentration%dataset%is_transient .or. &
            condition%concentration%gradient%is_transient .or. &
            condition%concentration%datum%is_transient) then
          update = .true.
        endif
        
        if (update) then ! for now, everything transport is dirichlet-type
          do idof = 1, condition%num_sub_conditions
            coupler%aux_real_var(idof,1:num_connections) = &
              condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
          enddo
        endif

      endif
      
    endif

    coupler => coupler%next
  enddo

end subroutine RealizationUpdateCouplerAuxVars

! ************************************************************************** !
!
! RealizationUpdate: Update parameters in realization (e.g. conditions, bcs, srcs)
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine RealizationUpdate(realization)

  implicit none
  
  type(realization_type) :: realization
  
  logical :: force_update_flag = .false.
  
  ! must update conditions first
  call ConditionUpdate(realization%flow_conditions,realization%option, &
                       realization%option%time)
  call ConditionUpdate(realization%transport_conditions,realization%option, &
                       realization%option%time)
  call RealizationUpdateCouplerAuxVars(realization, &
                                       realization%flow_boundary_conditions, &
                                       force_update_flag)
  call RealizationUpdateCouplerAuxVars(realization, &
                                       realization%transport_boundary_conditions, &
                                       force_update_flag)
! currently don't use aux_vars, just condition for src/sinks
!  call RealizationUpdateSrcSinks(realization)

end subroutine RealizationUpdate

! ************************************************************************** !
!
! RealizationAddWaypointsToList: Creates waypoints assoiciated with source/sinks
!                             boundary conditions, etc. and add to list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RealizationAddWaypointsToList(realization)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string
  type(waypoint_list_type), pointer :: waypoint_list
  type(coupler_type), pointer :: coupler
  type(sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint
  PetscInt :: itime, isub_condition

  waypoint_list => realization%waypoints

  if (realization%option%nflowdof > 0) then
  ! FLOW ----------------
    ! boundary conditions
    coupler => realization%flow_boundary_conditions%first
    do
      if (.not.associated(coupler)) exit
      ! the only way a boundary condition is included as a waypoint is if
      ! its size is 1 and time > 0.
      do isub_condition = 1, coupler%condition%num_sub_conditions
        sub_condition => coupler%condition%sub_condition_ptr(isub_condition)%ptr
        itime = 1
        if (sub_condition%dataset%max_time_index == 1 .and. &
            sub_condition%dataset%times(itime) > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = coupler%condition%pressure%dataset%times(itime)
          waypoint%update_bcs = .true.
          call WaypointInsertInList(waypoint,waypoint_list)
          exit
        endif
      enddo
      coupler => coupler%next
    enddo

    ! source/sinks
    coupler => realization%flow_source_sinks%first
    do
      if (.not.associated(coupler)) exit
      do isub_condition = 1, coupler%condition%num_sub_conditions
        sub_condition => coupler%condition%sub_condition_ptr(isub_condition)%ptr
        do itime=1,sub_condition%dataset%max_time_index
          if (sub_condition%dataset%times(itime) > 1.d-40) then
            waypoint => WaypointCreate()
            waypoint%time = coupler%condition%pressure%dataset%times(itime)
            waypoint%update_srcs = .true.
            call WaypointInsertInList(waypoint,waypoint_list)
          endif
        enddo
      enddo
      coupler => coupler%next
    enddo
  endif
    
  if (realization%option%ntrandof > 0) then
  ! TRANSPORT ----------------
    ! boundary conditions
    coupler => realization%transport_boundary_conditions%first
    do
      if (.not.associated(coupler)) exit
      ! the only way a boundary condition is included as a waypoint is if
      ! its size is 1 and time > 0.
      do isub_condition = 1, coupler%condition%num_sub_conditions
        sub_condition => coupler%condition%sub_condition_ptr(isub_condition)%ptr
        itime = 1
        if (sub_condition%dataset%max_time_index == 1 .and. &
            sub_condition%dataset%times(itime) > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = coupler%condition%pressure%dataset%times(itime)
          waypoint%update_bcs = .true.
          call WaypointInsertInList(waypoint,waypoint_list)
          exit
        endif
      enddo
      coupler => coupler%next
    enddo

    ! source/sinks
    coupler => realization%transport_source_sinks%first
    do
      if (.not.associated(coupler)) exit
      do isub_condition = 1, coupler%condition%num_sub_conditions
        sub_condition => coupler%condition%sub_condition_ptr(isub_condition)%ptr
        do itime=1,sub_condition%dataset%max_time_index
          if (sub_condition%dataset%times(itime) > 1.d-40) then
            waypoint => WaypointCreate()
            waypoint%time = coupler%condition%pressure%dataset%times(itime)
            waypoint%update_srcs = .true.
            call WaypointInsertInList(waypoint,waypoint_list)
          endif
        enddo
      enddo
      coupler => coupler%next
    enddo
  endif
  
end subroutine RealizationAddWaypointsToList

! ************************************************************************** !
!
! RealizationDestroy: Deallocates a realization
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine RealizationDestroy(realization)

  implicit none
  
  type(realization_type), pointer :: realization
  
  if (.not.associated(realization)) return
    
  call GridDestroy(realization%grid)
  call FieldDestroy(realization%field)
  call OptionDestroy(realization%option)
  call RegionDestroyList(realization%regions)
  
  call ConditionDestroyList(realization%flow_conditions)
  call CouplerDestroyList(realization%flow_boundary_conditions)
  call CouplerDestroyList(realization%flow_initial_conditions)
  call CouplerDestroyList(realization%flow_source_sinks)
  
  call ConditionDestroyList(realization%transport_conditions)
  call CouplerDestroyList(realization%transport_boundary_conditions)
  call CouplerDestroyList(realization%transport_initial_conditions)
  call CouplerDestroyList(realization%transport_source_sinks)
  
  call StrataDestroyList(realization%strata)
  call BreakthroughDestroyList(realization%breakthrough)
  
  if (associated(realization%debug)) deallocate(realization%debug)
  nullify(realization%debug)
  
  if (associated(realization%material_array)) &
    deallocate(realization%material_array)
  nullify(realization%material_array)
  call MaterialDestroy(realization%materials)
  
  if (associated(realization%saturation_function_array)) &
    deallocate(realization%saturation_function_array)
  nullify(realization%saturation_function_array)
  call SaturationFunctionDestroy(realization%saturation_functions)
  
  call ReactiveTransportAuxDestroy(realization%RTaux)
    
end subroutine RealizationDestroy
  
end module Realization_module
