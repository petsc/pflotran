module Realization_module

  use Option_module
  use Region_module
  use Condition_module
  use Material_module
  use Discretization_module
  use Field_module
  use Debug_module
  use Waypoint_module
  
  use Level_module
  use Patch_module
  
  implicit none

private

#include "definitions.h"

  type, public :: realization_type

    type(discretization_type), pointer :: discretization
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
            RealizationUpdate, RealizationAddWaypointsToList, &
            RealizationCreateDiscretization, &
            RealizationLocalizeRegions, &
            RealizationAddCoupler, RealizationAddStrata, &
            RealizationAddBreakthrough, RealizAssignInitialConditions, &
            RealizAssignUniformVelocity, RealizAssignTransportInitCond
  
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
  realization%discretization => DiscretizationCreate()
  realization%option => OptionCreate()
  realization%field => FieldCreate()
  realization%debug => DebugCreatePflow()
  realization%output_option => OutputOptionCreate()

  realization%level_list => LevelCreateList()

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
! RealizationCreateDiscretization: Creates grid
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationCreateDiscretization(realization)

  use Grid_module

  implicit none
  
  type(realization_type) :: realization
  
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  option => realization%option
  field => realization%field
  
  discretization => realization%discretization

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
    
      grid => discretization%grid
      
      call DiscretizationCreateDMs(discretization,option)
      
      ! 1 degree of freedom, global
      call DiscretizationCreateVector(discretization,ONEDOF,field%porosity0,GLOBAL)
      call DiscretizationDuplicateVector(discretization,field%porosity0, field%volume)
      
      ! 1 degree of freedom, local
      call DiscretizationCreateVector(discretization,ONEDOF,field%porosity_loc,LOCAL)
      call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%tor_loc)
      
      if (option%nflowdof > 0) then

        ! 1-dof global  
        call DiscretizationDuplicateVector(discretization,field%porosity0, field%perm0_xx)
        call DiscretizationDuplicateVector(discretization,field%porosity0, field%perm0_yy)
        call DiscretizationDuplicateVector(discretization,field%porosity0, field%perm0_zz)
        call DiscretizationDuplicateVector(discretization,field%porosity0, field%perm_pow)

        ! 1-dof local
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%ithrm_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%icap_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%iphas_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%iphas_old_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%perm_xx_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%perm_yy_loc)
        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%perm_zz_loc)

        ! ndof degrees of freedom, global
        call DiscretizationCreateVector(discretization,NFLOWDOF, field%flow_xx, GLOBAL)
        call DiscretizationDuplicateVector(discretization,field%flow_xx, field%flow_yy)
        call DiscretizationDuplicateVector(discretization,field%flow_xx, field%flow_dxx)
        call DiscretizationDuplicateVector(discretization,field%flow_xx, field%flow_r)
        call DiscretizationDuplicateVector(discretization,field%flow_xx, field%flow_accum)

        ! ndof degrees of freedom, local
        call DiscretizationCreateVector(discretization,NFLOWDOF, field%flow_xx_loc, LOCAL)
      endif

      if (option%ntrandof > 0) then
        ! ndof degrees of freedom, global
        call DiscretizationCreateVector(discretization,NTRANDOF, field%tran_xx, GLOBAL)
        call DiscretizationDuplicateVector(discretization,field%tran_xx, field%tran_yy)
        call DiscretizationDuplicateVector(discretization,field%tran_xx, field%tran_dxx)
        call DiscretizationDuplicateVector(discretization,field%tran_xx, field%tran_r)
        call DiscretizationDuplicateVector(discretization,field%tran_xx, field%tran_accum)

        call DiscretizationDuplicateVector(discretization,field%porosity_loc, field%saturation_loc)
        
        ! ndof degrees of freedom, local
        call DiscretizationCreateVector(discretization,NTRANDOF, field%tran_xx_loc, LOCAL)
      endif

      ! set up nG2L, NL2G, etc.
      call GridMapIndices(grid)
      call GridComputeSpacing(grid)
      call GridComputeCoordinates(grid,option)
      call GridComputeVolumes(grid,field%volume,option)
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option)

    case(AMR_GRID)
    
  end select      

end subroutine RealizationCreateDiscretization

! ************************************************************************** !
!
! RealizationLocalizeRegions: Localizes regions within each patch
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationLocalizeRegions(realization)

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
      call PatchLocalizeRegions(cur_patch,realization%regions, &
                                realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
 
end subroutine RealizationLocalizeRegions

! ************************************************************************** !
!
! RealizationAddCoupler: Adds a copy of a coupler to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddCoupler(realization,coupler)

  use Coupler_module

  implicit none
  
  type(realization_type) :: realization
  type(coupler_type), pointer :: coupler
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      ! only add to flow list for now, since they will be split out later
      select case(coupler%itype)
        case(BOUNDARY_COUPLER_TYPE)
          call CouplerAddToList(CouplerCreate(coupler),cur_patch%flow_boundary_conditions)
        case(INITIAL_COUPLER_TYPE)
          call CouplerAddToList(CouplerCreate(coupler),cur_patch%flow_initial_conditions)
        case(SRC_SINK_COUPLER_TYPE)
          call CouplerAddToList(CouplerCreate(coupler),cur_patch%flow_source_sinks)
      end select
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  call CouplerDestroy(coupler)
 
end subroutine RealizationAddCoupler

! ************************************************************************** !
!
! RealizationAddStrata: Adds a copy of a strata to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddStrata(realization,strata)

  use Strata_module

  implicit none
  
  type(realization_type) :: realization
  type(strata_type), pointer :: strata
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call StrataAddToList(StrataCreate(strata),cur_patch%strata)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  call StrataDestroy(strata)
 
end subroutine RealizationAddStrata

! ************************************************************************** !
!
! RealizationAddBreakthrough: Adds a copy of a breakthrough object to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddBreakthrough(realization,breakthrough)

  use Breakthrough_module

  implicit none
  
  type(realization_type) :: realization
  type(breakthrough_type), pointer :: breakthrough
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call BreakthroughAddToList(BreakthroughCreate(breakthrough), &
                                 cur_patch%breakthrough)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  call BreakthroughDestroy(breakthrough)
 
end subroutine RealizationAddBreakthrough

! ************************************************************************** !
!
! RealizationProcessCouplers: Sets connectivity and pointers for couplers
! author: Glenn Hammond
! date: 02/22/08
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
! date: 02/22/08
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
! date: 02/22/08
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
! RealizAssignInitialConditions: Assigns initial conditions to model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine RealizAssignInitialConditions(realization)
  
  implicit none

  type(realization_type) :: realization

  call RealizAssignFlowInitCond(realization)
  call RealizAssignTransportInitCond(realization)

end subroutine RealizAssignInitialConditions

! ************************************************************************** !
!
! RealizAssignFlowInitCond: Assigns flow initial conditions to model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine RealizAssignFlowInitCond(realization)

  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  use Patch_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      if (option%nflowdof > 0) then
      
        select case(option%iflowmode)
          case(RICHARDS_LITE_MODE)
          case(RICHARDS_MODE)
          case(MPH_MODE)
!            call pflow_mphase_setupini(realization)
        end select 

        ! assign initial conditions values to domain
        call VecGetArrayF90(field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
        call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr)
        
        xx_p = -999.d0
        
        initial_condition => cur_patch%flow_initial_conditions%first
        do
        
          if (.not.associated(initial_condition)) exit

          if (.not.associated(initial_condition%connection)) then
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = patch%grid%nL2G(local_id)
              iend = local_id*option%nflowdof
              ibegin = iend-option%nflowdof+1
              if (associated(patch%imat)) then
                if (patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  iphase_loc_p(ghosted_id) = 0
                  cycle
                endif
              endif
              do idof = 1, option%nflowdof
                xx_p(ibegin+idof-1) = &
                  initial_condition%condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
              enddo
              iphase_loc_p(ghosted_id)=initial_condition%condition%iphase
            enddo
          else
            do iconn=1,initial_condition%connection%num_connections
              local_id = initial_condition%connection%id_dn(iconn)
              ghosted_id = patch%grid%nL2G(local_id)
              iend = local_id*option%nflowdof
              ibegin = iend-option%nflowdof+1
              if (associated(patch%imat)) then
                if (patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  iphase_loc_p(ghosted_id) = 0
                  cycle
                endif
              endif
              xx_p(ibegin:iend) = &
                initial_condition%aux_real_var(1:option%nflowdof,iconn)
              iphase_loc_p(ghosted_id)=initial_condition%aux_int_var(1,iconn)
            enddo
          endif
          initial_condition => initial_condition%next
        enddo
        
        call VecRestoreArrayF90(field%flow_xx,xx_p, ierr)
        call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr)
        
        ! update dependent vectors
        call DiscretizationGlobalToLocal(discretization,field%flow_xx,field%flow_xx_loc,NFLOWDOF)  
        call VecCopy(field%flow_xx, field%flow_yy, ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)  
        call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_old_loc,ONEDOF)
        
      endif
      
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
end subroutine RealizAssignFlowInitCond

! ************************************************************************** !
!
! RealizAssignTransportInitCond: Assigns transport initial conditions to model
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
subroutine RealizAssignTransportInitCond(realization)

  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  use Patch_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      if (option%ntrandof > 0) then

        ! assign initial conditions values to domain
        call VecGetArrayF90(field%tran_xx,xx_p, ierr); CHKERRQ(ierr)
        
        xx_p = -999.d0
        
        initial_condition => cur_patch%transport_initial_conditions%first
        do
        
          if (.not.associated(initial_condition)) exit

          if (.not.associated(initial_condition%connection)) then
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = patch%grid%nL2G(local_id)
              iend = local_id*option%ntrandof
              ibegin = iend-option%ntrandof+1
              if (associated(patch%imat)) then
                if (patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  cycle
                endif
              endif
              do idof = 1, option%ntrandof
                xx_p(ibegin+idof-1) = &
                  initial_condition%condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
              enddo
            enddo
          else
            do iconn=1,initial_condition%connection%num_connections
              local_id = initial_condition%connection%id_dn(iconn)
              ghosted_id = patch%grid%nL2G(local_id)
              iend = local_id*option%ntrandof
              ibegin = iend-option%ntrandof+1
              if (associated(patch%imat)) then
                if (patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  cycle
                endif
              endif
              xx_p(ibegin:iend) = &
                initial_condition%aux_real_var(1:option%ntrandof,iconn)
            enddo
          endif
          initial_condition => initial_condition%next
        enddo
        
        call VecRestoreArrayF90(field%tran_xx,xx_p, ierr)
        
        ! update dependent vectors
        call DiscretizationGlobalToLocal(discretization,field%tran_xx,field%tran_xx_loc,NTRANDOF)  
        call VecCopy(field%tran_xx, field%tran_yy, ierr)

      endif

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
end subroutine RealizAssignTransportInitCond

! ************************************************************************** !
!
! RealizAssignUniformVelocity: Assigns uniform velocity for transport
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizAssignUniformVelocity(realization)

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
      call PatchAssignUniformVelocity(cur_patch,realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
 
end subroutine RealizAssignUniformVelocity

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
