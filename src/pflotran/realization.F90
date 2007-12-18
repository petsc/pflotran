module Realization_module

  use Grid_module
  use Option_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Material_module
  use Strata_module
  use Field_module

  implicit none

private

#include "definitions.h"

  type, public :: realization_type

    type(grid_type), pointer :: grid
    type(option_type), pointer :: option
    type(field_type), pointer :: field
    type(output_option_type), pointer :: output_option
    type(region_list_type), pointer :: regions
    type(condition_list_type), pointer :: conditions
    type(coupler_list_type), pointer :: boundary_conditions
    type(coupler_list_type), pointer :: initial_conditions
    type(coupler_list_type), pointer :: source_sinks
    type(strata_list_type), pointer :: strata
    
    type(material_type), pointer :: materials
    type(material_ptr_type), pointer :: material_array(:)
    type(thermal_property_type), pointer :: thermal_properties
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)

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
  realization%output_option => OutputOptionCreate()
  nullify(realization%grid)
  allocate(realization%regions)
  call RegionInitList(realization%regions)
  allocate(realization%conditions)
  call ConditionInitList(realization%conditions)
  allocate(realization%boundary_conditions)
  call CouplerInitList(realization%boundary_conditions)
  allocate(realization%initial_conditions)
  call CouplerInitList(realization%initial_conditions)
  allocate(realization%source_sinks)
  call CouplerInitList(realization%source_sinks)
  allocate(realization%strata)
  call StrataInitList(realization%strata)
  
  nullify(realization%materials)
  nullify(realization%thermal_properties)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  
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
  type(strata_type), pointer :: strata


  ! boundary conditions
  coupler => realization%boundary_conditions%first
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
                                                 realization%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in boundary condition list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => realization%initial_conditions%first
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
                                                 realization%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in initial condition list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => realization%source_sinks%first
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
                                                 realization%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in source/sink list'
      call printErrMsg(realization%option,string)
    endif
    coupler => coupler%next
  enddo
  
    
  ! strata
  ! connect pointers from strata to regions
  strata => realization%strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    strata%region => RegionGetPtrFromList(strata%region_name, &
                                                realization%regions)
    if (.not.associated(strata%region)) then
      string = 'Region ' // trim(strata%region_name) // &
               ' not found in strata list'
      call printErrMsg(realization%option,string)
    endif
    if (strata%active) then
      ! pointer to material
      strata%material => &
                            MaterialGetPtrFromList(strata%material_name, &
                                                   realization%materials)
      if (.not.associated(strata%material)) then
        string = 'Material ' // trim(strata%material_name) // &
                 ' not found in unit list'
        call printErrMsg(realization%option,string)
      endif
    else
      nullify(strata%material)
    endif
    strata => strata%next
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
  type(coupler_list_type) :: coupler_list
  
  integer :: num_connections
  logical :: force_update_flag
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: coupler
    
  option => realization%option
    
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    
    if (associated(coupler%connection)) then
      num_connections = coupler%connection%num_connections

      ! allocate arrays that match the number of connections
      select case(option%imode)

        case(FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
       
          allocate(coupler%aux_real_var(option%ndof*option%nphase,num_connections))
          allocate(coupler%aux_int_var(1,num_connections))
          coupler%aux_real_var = 0.d0
          coupler%aux_int_var = 0

        case(MPH_MODE)

          allocate(coupler%aux_real_var(option%ndof*option%nphase,num_connections))
          allocate(coupler%aux_int_var(1,num_connections))
          coupler%aux_real_var = 0.d0
          coupler%aux_int_var = 0
            
        case default
      end select 
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
  type(coupler_list_type) :: coupler_list
  logical :: force_update_flag
  
  type(coupler_type), pointer :: coupler
  
  integer :: idof, num_connections
 
  coupler => coupler_list%first
  
  do
    if (.not.associated(coupler)) exit
    
    if ((coupler%condition%is_transient .or. &
         force_update_flag) .and. &
        associated(coupler%aux_real_var)) then
        
      select case(realization%option%imode)
        case(RICHARDS_MODE,MPH_MODE)
          select case(coupler%condition%itype(RICHARDS_PRESSURE_DOF))
            case(DIRICHLET_BC,NEUMANN_BC,MASS_RATE,ZERO_GRADIENT_BC)
              num_connections = coupler%connection%num_connections
              do idof = 1,realization%option%ndof
                coupler%aux_real_var(idof,1:num_connections) = &
                  coupler%condition%cur_value(idof)
              enddo
              coupler%aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                coupler%condition%iphase
            case(HYDROSTATIC_BC)
!              call HydrostaticUpdateCoupler(coupler,realization%option,realization%grid)
              call HydrostaticUpdateCouplerBetter(coupler,realization%option,realization%grid)
          end select
      end select
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
  call ConditionUpdate(realization%conditions,realization%option,realization%option%time)
  call RealizationUpdateCouplerAuxVars(realization,realization%boundary_conditions, &
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
subroutine RealizationAddWaypointsToList(realization,waypoint_list)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  type(realization_type) :: realization
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(waypoint_type), pointer :: waypoint
  integer :: itime

#if 0  
  ! Ignore boundary conditions for now
  ! boundary conditions
  coupler => realization%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    endif
    coupler => coupler%next
  enddo
#endif  

  ! source/sinks
  coupler => realization%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    do itime=1,coupler%condition%num_values
      if (coupler%condition%times(itime) > 1.d-40) then
        waypoint => WaypointCreate()
        waypoint%time = coupler%condition%times(itime)
        waypoint%update_srcs = .true.
        call WaypointInsertInList(waypoint,waypoint_list)
      endif
    enddo
    coupler => coupler%next
  enddo
  
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
  call ConditionDestroyList(realization%conditions)
  call CouplerDestroyList(realization%boundary_conditions)
  call CouplerDestroyList(realization%initial_conditions)
  call CouplerDestroyList(realization%source_sinks)
  call StrataDestroyList(realization%strata)
  
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
