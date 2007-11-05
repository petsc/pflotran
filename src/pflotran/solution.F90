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

  public :: SolutionCreate, SolutionDestroy, SolutionProcessCouplers
  
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
    coupler%flow_condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                      solution%conditions)
    if (.not.associated(coupler%flow_condition)) then
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
    coupler%flow_condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                      solution%conditions)
    if (.not.associated(coupler%flow_condition)) then
      string = 'Condition ' // trim(coupler%condition_name) // &
               ' not found in initial condition list'
      call printErrMsg(solution%option,string)
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => solution%boundary_conditions%first
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
    coupler%flow_condition => ConditionGetPtrFromList(coupler%condition_name, &
                                                      solution%conditions)
    if (.not.associated(coupler%flow_condition)) then
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
  
    
  call GridLocalizeRegions(solution%regions,solution%grid,solution%option)
  call GridComputeBoundaryConnect2(solution%grid,solution%option, &
                                    solution%boundary_conditions)
    
end subroutine SolutionProcessCouplers

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
  call CouplerCreateList(solution%boundary_conditions)
  call CouplerCreateList(solution%initial_conditions)
  call CouplerCreateList(solution%source_sinks)
  call StrataDestroyList(solution%strata)
    
end subroutine SolutionDestroy
  
end module Solution_module
