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
  
  use Level_module
  use Patch_module
  
  use Reactive_Transport_Aux_module
  
  implicit none

private

#include "definitions.h"

  type, public :: realization_type

    type(level_list_type), pointer :: level_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(field_type), pointer :: field
    type(pflow_debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    type(region_list_type), pointer :: regions
    
    type(condition_list_type), pointer :: flow_conditions
    type(condition_list_type), pointer :: transport_conditions
    
    type(material_type), pointer :: materials
    type(material_ptr_type), pointer :: material_array(:)
    type(thermal_property_type), pointer :: thermal_properties
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    
    type(waypoint_list_type), pointer :: waypoints
    
  end type realization_type

  public :: RealizationCreate, RealizationDestroy, &
            RealizationProcessCouplers, &
            RealizationInitAllCouplerAuxVars, &
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

  allocate(realization%regions)
  call RegionInitList(realization%regions)

  allocate(realization%flow_conditions)
  call ConditionInitList(realization%flow_conditions)

  allocate(realization%transport_conditions)
  call ConditionInitList(realization%transport_conditions)

  nullify(realization%materials)
  nullify(realization%material_array)
  nullify(realization%thermal_properties)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  
  RealizationCreate => realization
  
end function RealizationCreate  

! ************************************************************************** !
!
! RealizationProcessCouplers: Deallocates a realization
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RealizationProcessCouplers(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchProcessCouplers(cur_patch,realization%flow_conditions, &
                                realization%transport_conditions, &
                                realization%materials,realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
 
end subroutine RealizationProcessCouplers

! ************************************************************************** !
!
! RealizationInitCouplerAuxVars: Initializes coupler auxillary variables 
!                                within list
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RealizationInitAllCouplerAuxVars(realization)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

 cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchInitAllCouplerAuxVars(cur_patch,realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
end subroutine RealizationInitAllCouplerAuxVars

! ************************************************************************** !
!
! RealizUpdateAllCouplerAuxVars: Updates auxilliary variables associated 
!                                  with couplers in lis
! author: Glenn Hammond
! date: 00/00/00
!
! ************************************************************************** !
subroutine RealizUpdateAllCouplerAuxVars(realization,force_update_flag)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  logical :: force_update_flag

  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

 cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchUpdateAllCouplerAuxVars(cur_patch,force_update_flag, &
                                     realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
end subroutine RealizUpdateAllCouplerAuxVars

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
  call RealizUpdateAllCouplerAuxVars(realization,force_update_flag)
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
  type(condition_type), pointer :: cur_condition
  type(sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint
  PetscInt :: itime, isub_condition

  waypoint_list => realization%waypoints

  if (realization%option%nflowdof > 0) then
  ! FLOW ----------------
    ! boundary conditions
    cur_condition => realization%flow_conditions%first
    do
      if (.not.associated(cur_condition)) exit
      if (cur_condition%sync_time_with_update) then
        do isub_condition = 1, cur_condition%num_sub_conditions
          sub_condition => cur_condition%sub_condition_ptr(isub_condition)%ptr
          itime = 1
          if (sub_condition%dataset%max_time_index == 1 .and. &
              sub_condition%dataset%times(itime) > 1.d-40) then
            waypoint => WaypointCreate()
            waypoint%time = cur_condition%pressure%dataset%times(itime)
            waypoint%update_bcs = .true.
            call WaypointInsertInList(waypoint,waypoint_list)
            exit
          endif
        enddo
      endif
      cur_condition => cur_condition%next
    enddo
  endif
    
  if (realization%option%ntrandof > 0) then
  ! TRANSPORT ----------------
    ! boundary conditions
    cur_condition => realization%transport_conditions%first
    do
      if (.not.associated(cur_condition)) exit
      if (cur_condition%sync_time_with_update) then
        do isub_condition = 1, cur_condition%num_sub_conditions
          sub_condition => cur_condition%sub_condition_ptr(isub_condition)%ptr
          itime = 1
          if (sub_condition%dataset%max_time_index == 1 .and. &
              sub_condition%dataset%times(itime) > 1.d-40) then
            waypoint => WaypointCreate()
            waypoint%time = cur_condition%pressure%dataset%times(itime)
            waypoint%update_bcs = .true.
            call WaypointInsertInList(waypoint,waypoint_list)
            exit
          endif
        enddo
      endif
      cur_condition => cur_condition%next
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
    
  call FieldDestroy(realization%field)
  call OptionDestroy(realization%option)
  call RegionDestroyList(realization%regions)
  
  call ConditionDestroyList(realization%flow_conditions)
  call ConditionDestroyList(realization%transport_conditions)

  call LevelDestroyList(realization%level_list)

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
  
end subroutine RealizationDestroy
  
end module Realization_module
