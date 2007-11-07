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

  public :: SolutionCreate, SolutionDestroy, &
            SolutionProcessCouplers, &
            SolutionUpdateBoundaryConditions, &
            SolutionSetIBNDTYPE
  
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
    
end subroutine SolutionProcessCouplers

! ************************************************************************** !
!
! SolutionUpdateBoundaryConditions: Updates boundary conditions within model
! author: Glenn Hammond
! date: 11/06/07
!
! ************************************************************************** !
subroutine SolutionUpdateBoundaryConditions(solution)

  implicit none
  
  type(solution_type) :: solution
  
  integer :: icell, idof, count, num_connections
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
    
  option => solution%option
  grid => solution%grid
    
  num_connections = grid%boundary_connection_list%first%num_connections
  select case(option%imode)

    case(MPH_MODE,FLASH_MODE,RICHARDS_MODE,OWG_MODE,VADOSE_MODE)
  
      allocate(option%xxbc(option%ndof,num_connections))
      allocate(option%iphasebc(num_connections))
      allocate(option%velocitybc(option%nphase,num_connections))
      option%xxbc = 0.d0
      option%iphasebc = 0
      option%velocitybc = 0.d0
        
      count = 0
      boundary_condition => solution%boundary_conditions%first
      do
        if (.not.associated(boundary_condition)) exit
        do icell=1,boundary_condition%region%num_cells
          count = count + 1
          option%iphasebc(count) = boundary_condition%flow_condition%iphase
          do idof=1,option%ndof
            select case(boundary_condition%flow_condition%itype(idof))
              case(DIRICHLET_BC)
                option%xxbc(idof,count) = &
                  boundary_condition%flow_condition%cur_value(idof)
              case(NEUMANN_BC)
                option%velocitybc(1:option%nphase,count) = &
                  boundary_condition%flow_condition%cur_value(idof)
            end select
          enddo
        enddo
        boundary_condition => boundary_condition%next
      enddo
      
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

      count = 0
      boundary_condition => solution%boundary_conditions%first
      do
        if (.not.associated(boundary_condition)) exit
        do icell=1,boundary_condition%region%num_cells
          count = count + 1
          option%iphasebc(count) = boundary_condition%flow_condition%iphase
          if (boundary_condition%flow_condition%itype(1) == DIRICHLET_BC) then
            option%pressurebc(:,count) = boundary_condition%flow_condition%cur_value(1)
          else
            option%velocitybc(:,count) = boundary_condition%flow_condition%cur_value(1)
          endif
          option%tempbc(icell) = boundary_condition%flow_condition%cur_value(2)
          option%concbc(icell) = boundary_condition%flow_condition%cur_value(3)
          option%sgbc(icell) = 1.d0-boundary_condition%flow_condition%cur_value(4) ! read in as sl
        enddo
        boundary_condition => boundary_condition%next
      enddo

  end select 

end subroutine SolutionUpdateBoundaryConditions

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
    option%ibndtyp(count) = boundary_condition%flow_condition%itype(1)  ! hardwired to dof=1 (pressure)
    boundary_condition => boundary_condition%next
  enddo

end subroutine SolutionSetIBNDTYPE

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
