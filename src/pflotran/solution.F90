module Solution_module

  use Grid_module
  use Option_module
  use Region_module
  use Condition_module
  use Coupler_module
  use Material_module
  use Strata_module

  implicit none

private

#include "definitions.h"

  type, public :: solution_type

    type(grid_type), pointer :: grid
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    type(region_list_type), pointer :: regions
    type(condition_list_type), pointer :: conditions
    type(coupler_list_type), pointer :: boundary_conditions
    type(coupler_list_type), pointer :: initial_conditions
    type(coupler_list_type), pointer :: source_sinks
    type(strata_list_type), pointer :: strata
    
    type(material_type), pointer :: materials
    type(thermal_property_type), pointer :: thermal_properties
    type(saturation_function_type), pointer :: saturation_functions

  end type solution_type

  public :: SolutionCreate, SolutionDestroy, &
            SolutionProcessCouplers, &
            SolutionInitBoundConditions, &
            SolutionUpdate, SolutionAddWaypointsToList
  
contains
  
! ************************************************************************** !
!
! SolutionCreate: Allocates and initializes a new Solution object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function SolutionCreate()

  implicit none
  
  type(solution_type), pointer :: SolutionCreate
  
  type(solution_type), pointer :: solution
  
  allocate(solution)
  solution%option => OptionCreate()
  solution%output_option => OutputOptionCreate()
  nullify(solution%grid)
  allocate(solution%regions)
  call RegionInitList(solution%regions)
  allocate(solution%conditions)
  call ConditionInitList(solution%conditions)
  allocate(solution%boundary_conditions)
  call CouplerInitList(solution%boundary_conditions)
  allocate(solution%initial_conditions)
  call CouplerInitList(solution%initial_conditions)
  allocate(solution%source_sinks)
  call CouplerInitList(solution%source_sinks)
  allocate(solution%strata)
  call StrataInitList(solution%strata)
  
  nullify(solution%materials)
  nullify(solution%thermal_properties)
  nullify(solution%saturation_functions)
  
  SolutionCreate => solution
  
end function SolutionCreate  

! ************************************************************************** !
!
! SolutionProcessCouplers: Deallocates a solution
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SolutionProcessCouplers(solution)

  use Option_module

  implicit none
  
  type(solution_type) :: solution
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata


  ! boundary conditions
  coupler => solution%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           solution%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in boundary condition list'
      call printErrMsg(solution%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 solution%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in boundary condition list'
      call printErrMsg(solution%option,string)
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => solution%initial_conditions%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           solution%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in initial condition list'
      call printErrMsg(solution%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 solution%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in initial condition list'
      call printErrMsg(solution%option,string)
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => solution%source_sinks%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           solution%regions)
    if (.not.associated(coupler%region)) then
      string = 'Region ' // trim(coupler%region_name) // &
               ' not found in source/sink list'
      call printErrMsg(solution%option,string)
    endif
    ! pointer to flow condition
    coupler%condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                 solution%conditions)
    if (.not.associated(coupler%condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in source/sink list'
      call printErrMsg(solution%option,string)
    endif
    coupler => coupler%next
  enddo
  
    
  ! strata
  ! connect pointers from strata to regions
  strata => solution%strata%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    strata%region => RegionGetPtrFromList(strata%region_name, &
                                                solution%regions)
    if (.not.associated(strata%region)) then
      string = 'Region ' // trim(strata%region_name) // &
               ' not found in strata list'
      call printErrMsg(solution%option,string)
    endif
    ! pointer to material
    strata%material => &
                          MaterialGetPtrFromList(strata%material_name, &
                                                 solution%materials)
    if (.not.associated(strata%material)) then
      string = 'Material ' // trim(strata%material_name) // &
               ' not found in unit list'
      call printErrMsg(solution%option,string)
    endif
    strata => strata%next
  enddo 
    
end subroutine SolutionProcessCouplers

! ************************************************************************** !
!
! SolutionInitBoundConditions: Initializes boundary conditions within model
! author: Glenn Hammond
! date: 11/06/07
!
! ************************************************************************** !
subroutine SolutionInitBoundConditions(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: num_connections
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
    
  option => solution%option
    
  ! sum the number of connections among all boundary conditions
  num_connections = 0
  boundary_condition => solution%boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    num_connections = num_connections + &
                      boundary_condition%connection%num_connections
    boundary_condition => boundary_condition%next
  enddo
  
  ! allocate arrays that match the number of connections
  select case(option%imode)

    case(MPH_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
  
      allocate(option%xxbc(option%ndof,num_connections))
      allocate(option%iphasebc(num_connections))
      allocate(option%velocitybc(option%nphase,num_connections))
      option%xxbc = 0.d0
      option%iphasebc = 0
      option%velocitybc = 0.d0
  
    case default
    
      allocate(option%pressurebc(option%nphase,num_connections))
      allocate(option%tempbc(num_connections))
      allocate(option%sgbc(num_connections))
      allocate(option%concbc(num_connections))
      allocate(option%velocitybc(option%nphase,num_connections))
      allocate(option%iphasebc(num_connections))
      option%pressurebc = 0.d0
      option%tempbc = 0.d0
      option%concbc = 0.d0
      option%sgbc = 0.d0
      option%velocitybc = 0.d0
      option%iphasebc = 0

  end select 
  
  call SolutionUpdateBoundConditions(solution)

end subroutine SolutionInitBoundConditions

! ************************************************************************** !
!
! SolutionUpdateBoundConditions: Updates boundary conditions within model
! author: Glenn Hammond
! date: 11/06/07
!
! ************************************************************************** !
subroutine SolutionUpdateBoundConditions(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: icell, idof, count
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
    
  option => solution%option
 
  boundary_condition => solution%boundary_conditions%first
 
  count = 0
  do
  
    if (.not.associated(boundary_condition)) exit
  
    select case(option%imode)

      case(MPH_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
  
        do icell=1,boundary_condition%region%num_cells
          count = count + 1
          option%iphasebc(count) = boundary_condition%condition%iphase
          do idof=1,option%ndof
            select case(boundary_condition%condition%itype(idof))
              case(DIRICHLET_BC)
                option%xxbc(idof,count) = &
                  boundary_condition%condition%cur_value(idof)
              case(NEUMANN_BC)
                option%velocitybc(1:option%nphase,count) = &
                  boundary_condition%condition%cur_value(idof)
            end select
          enddo
        enddo
      
      case default

        do icell=1,boundary_condition%region%num_cells
          count = count + 1
          option%iphasebc(count) = boundary_condition%condition%iphase
          if (boundary_condition%condition%itype(1) == DIRICHLET_BC) then
            option%pressurebc(:,count) = boundary_condition%condition%cur_value(1)
          else
            option%velocitybc(:,count) = boundary_condition%condition%cur_value(1)
          endif
          option%tempbc(icell) = boundary_condition%condition%cur_value(2)
          option%concbc(icell) = boundary_condition%condition%cur_value(3)
          option%sgbc(icell) = 1.d0-boundary_condition%condition%cur_value(4) ! read in as sl
        enddo

    end select 
    boundary_condition => boundary_condition%next
  enddo

end subroutine SolutionUpdateBoundConditions

! ************************************************************************** !
!
! SolutionInitSrcSinks: Initializes source/sinks within model
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine SolutionInitSrcSinks(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: num_connections
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: source_sink
    
  option => solution%option
    
  ! sum the number of connections among all boundary conditions
  num_connections = 0
  source_sink => solution%source_sinks%first
  do
    if (.not.associated(source_sink)) exit
    num_connections = num_connections + &
                      source_sink%connection%num_connections
    source_sink => source_sink%next
  enddo
  
  ! allocate arrays that match the number of connections
  select case(option%imode)

    case(MPH_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
  
      allocate(option%xxbc(option%ndof,num_connections))
      allocate(option%iphasebc(num_connections))
      allocate(option%velocitybc(option%nphase,num_connections))
      option%xxbc = 0.d0
      option%iphasebc = 0
      option%velocitybc = 0.d0
  
    case default
    
      allocate(option%pressurebc(option%nphase,num_connections))
      allocate(option%tempbc(num_connections))
      allocate(option%sgbc(num_connections))
      allocate(option%concbc(num_connections))
      allocate(option%velocitybc(option%nphase,num_connections))
      allocate(option%iphasebc(num_connections))
      option%pressurebc = 0.d0
      option%tempbc = 0.d0
      option%concbc = 0.d0
      option%sgbc = 0.d0
      option%velocitybc = 0.d0
      option%iphasebc = 0

  end select 
  
  call SolutionUpdateSrcSinks(solution)

end subroutine SolutionInitSrcSinks

! ************************************************************************** !
!
! SolutionUpdateSrcSinks: Updates source/sinks within model
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine SolutionUpdateSrcSinks(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: icell, idof, count
  
  type(option_type), pointer :: option
  type(coupler_type), pointer :: source_sink
    
  option => solution%option
 
  source_sink => solution%source_sinks%first
 
  count = 0
  do
  
    if (.not.associated(source_sink)) exit
  
    select case(option%imode)

      case(MPH_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
  
        do icell=1,source_sink%region%num_cells
          count = count + 1
          option%iphasebc(count) = source_sink%condition%iphase
          do idof=1,option%ndof
            select case(source_sink%condition%itype(idof))
              case(DIRICHLET_BC)
                option%xxbc(idof,count) = &
                  source_sink%condition%cur_value(idof)
              case(NEUMANN_BC)
                option%velocitybc(1:option%nphase,count) = &
                  source_sink%condition%cur_value(idof)
            end select
          enddo
        enddo
      
      case default

        do icell=1,source_sink%region%num_cells
          count = count + 1
          option%iphasebc(count) = source_sink%condition%iphase
          if (source_sink%condition%itype(1) == DIRICHLET_BC) then
            option%pressurebc(:,count) = source_sink%condition%cur_value(1)
          else
            option%velocitybc(:,count) = source_sink%condition%cur_value(1)
          endif
          option%tempbc(icell) = source_sink%condition%cur_value(2)
          option%concbc(icell) = source_sink%condition%cur_value(3)
          option%sgbc(icell) = 1.d0-source_sink%condition%cur_value(4) ! read in as sl
        enddo

    end select 
    source_sink => source_sink%next
  enddo

end subroutine SolutionUpdateSrcSinks

#if 0
! NO LONGER NEEDED
! ************************************************************************** !
!
! SolutionSetIBNDTYPE: Sets values in ibndtyp array
! ibndtyp needs to be done away with!!!!
! author: Glenn Hammond
! date: 11/06/07
!
! ************************************************************************** !
subroutine SolutionSetIBNDTYPE(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: count, num_conditions
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
    
  option => solution%option
  grid => solution%grid

  num_conditions = 0
  boundary_condition => solution%boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    num_conditions = num_conditions + 1
    boundary_condition => boundary_condition%next
  enddo
    
  allocate(option%ibndtyp(num_conditions))
  option%ibndtyp = 0
    
  count = 0
  boundary_condition => solution%boundary_conditions%first
  do
    if (.not.associated(boundary_condition)) exit
    count = count + 1
    option%ibndtyp(count) = boundary_condition%condition%itype(1)  ! hardwired to dof=1 (pressure)
    boundary_condition => boundary_condition%next
  enddo

end subroutine SolutionSetIBNDTYPE
#endif

! ************************************************************************** !
!
! SolutionUpdate: Update parameters in solution (e.g. conditions, bcs, srcs)
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine SolutionUpdate(solution)

  implicit none
  
  type(solution_type) :: solution
  
  ! must update conditions first
  call ConditionUpdate(solution%conditions,solution%option,solution%option%time)
  call SolutionUpdateBoundConditions(solution)
  call SolutionUpdateSrcSinks(solution)

end subroutine SolutionUpdate

! ************************************************************************** !
!
! SolutionAddWaypointsToList: Creates waypoints assoiciated with source/sinks
!                             boundary conditions, etc. and add to list
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SolutionAddWaypointsToList(solution,waypoint_list)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  type(solution_type) :: solution
  
  character(len=MAXSTRINGLENGTH) :: string
  type(coupler_type), pointer :: coupler
  type(waypoint_type), pointer :: waypoint
  integer :: itime

  ! Ignore boundary conditions for now
  ! boundary conditions
#if 0  
  coupler => solution%boundary_conditions%first
  do
    if (.not.associated(coupler)) exit
    endif
    coupler => coupler%next
  enddo
#endif  

  ! source/sinks
  coupler => solution%boundary_conditions%first
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
  
end subroutine SolutionAddWaypointsToList

! ************************************************************************** !
!
! SolutionDestroy: Deallocates a solution
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine SolutionDestroy(solution)

  implicit none
  
  type(solution_type), pointer :: solution
  
  if (.not.associated(solution)) return
    
  call GridDestroy(solution%grid)
  call OptionDestroy(solution%option)
  call RegionDestroyList(solution%regions)
  call ConditionDestroyList(solution%conditions)
  call CouplerDestroyList(solution%boundary_conditions)
  call CouplerDestroyList(solution%initial_conditions)
  call CouplerDestroyList(solution%source_sinks)
  call StrataDestroyList(solution%strata)
    
end subroutine SolutionDestroy
  
end module Solution_module
