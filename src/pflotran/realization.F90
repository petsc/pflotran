module Realization_module

  use Option_module
  use Region_module
  use Condition_module
  use Material_module
  use Discretization_module
  use Field_module
  use Debug_module
  use Waypoint_module
  
  use Reaction_Aux_module
  
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
    type(tran_condition_list_type), pointer :: transport_conditions
    type(tran_constraint_list_type), pointer :: transport_constraints
    
    type(reaction_type), pointer :: reaction
    
    type(material_type), pointer :: materials
    type(material_ptr_type), pointer :: material_array(:)
    type(thermal_property_type), pointer :: thermal_properties
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    type(fluid_property_type), pointer :: fluid_properties
    
    type(waypoint_list_type), pointer :: waypoints
    
  end type realization_type

  public :: RealizationCreate, &
            RealizationDestroy, &
            RealizationProcessCouplers, &
            RealizationInitAllCouplerAuxVars, &
            RealizationProcessConditions, &
            RealizationUpdate, &
            RealizationAddWaypointsToList, &
            RealizationCreateDiscretization, &
            RealizationLocalizeRegions, &
            RealizationAddCoupler, &
            RealizationAddStrata, &
            RealizationAddBreakthrough, &
            RealizAssignInitialConditions, &
            RealizAssignFlowInitCond, &
            RealizAssignTransportInitCond, &
            RealizAssignUniformVelocity, &
            RealizationRevertFlowParameters, &
            RealizBridgeFlowAndTransport, &
            RealizationGetDataset, &
            RealizGetDatasetValueAtCell, &
            RealizationSetDataset, &
            RealizationPrintCouplers
            
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
  call TranConditionInitList(realization%transport_conditions)
  allocate(realization%transport_constraints)
  call TranConstraintInitList(realization%transport_constraints)

  nullify(realization%materials)
  nullify(realization%material_array)
  nullify(realization%thermal_properties)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  nullify(realization%fluid_properties)
  
  nullify(realization%reaction)
  
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
  use AMR_Grid_module
  implicit none
  
  type(realization_type) :: realization
  
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  option => realization%option
  field => realization%field
  
  discretization => realization%discretization
  
  call DiscretizationCreateDMs(discretization,option)
  
  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity0, &
                                  GLOBAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%volume)
  
  ! 1 degree of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity_loc, &
                                  LOCAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%tor_loc)
  
  if (option%nflowdof > 0) then

    ! 1-dof global  
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_xx)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_yy)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_zz)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm_pow)

    ! 1-dof local
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%ithrm_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%icap_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%iphas_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%iphas_old_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_xx_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_yy_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_zz_loc)

    ! ndof degrees of freedom, global
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx, &
                                    GLOBAL,option)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_yy)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_dxx)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_r)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_accum)

    ! ndof degrees of freedom, local
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx_loc, &
                                    LOCAL,option)
  endif

  if (option%ntrandof > 0) then
    ! ndof degrees of freedom, global
    call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                    GLOBAL,option)
    call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                       field%tran_yy)
    call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                       field%tran_dxx)
    call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                       field%tran_r)
    call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                       field%tran_accum)

    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%saturation_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%saturation0_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%density_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%density0_loc)
    
    ! ndof degrees of freedom, local
    call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                    LOCAL,option)
                                    
    if (realization%reaction%use_log_formulation) then
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_log_xx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx_loc, &
                                         field%tran_work_loc)
    endif
    
  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
    
      grid => discretization%grid

      ! set up nG2L, NL2G, etc.
      call GridMapIndices(grid)
      call GridComputeSpacing(grid)
      call GridComputeCoordinates(grid,discretization%origin,option)
      call GridComputeVolumes(grid,field%volume,option)
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option)

    case(AMR_GRID)
       call AMRGridComputeGeometryInformation(discretization%amrgrid, &
                                              discretization%origin, &
                                              field,option)
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
  type(coupler_type), pointer :: new_coupler
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      ! only add to flow list for now, since they will be split out later
      new_coupler => CouplerCreate(coupler)
      select case(coupler%itype)
        case(BOUNDARY_COUPLER_TYPE)
          call CouplerAddToList(new_coupler,cur_patch%boundary_conditions)
        case(INITIAL_COUPLER_TYPE)
          call CouplerAddToList(new_coupler,cur_patch%initial_conditions)
        case(SRC_SINK_COUPLER_TYPE)
          call CouplerAddToList(new_coupler,cur_patch%source_sinks)
      end select
      nullify(new_coupler)
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
  type(strata_type), pointer :: new_strata
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      new_strata => StrataCreate(strata)
      call StrataAddToList(new_strata,cur_patch%strata)
      nullify(new_strata)
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
  type(breakthrough_type), pointer :: new_breakthrough
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      new_breakthrough => BreakthroughCreate(breakthrough)
      call BreakthroughAddToList(new_breakthrough, &
                                 cur_patch%breakthrough)
      nullify(new_breakthrough)
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
                                realization%material_array, &
                                realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
 
end subroutine RealizationProcessCouplers

! ************************************************************************** !
!
! RealizationProcessConditions: Sets up auxilliary data associated with 
!                               conditions
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine RealizationProcessConditions(realization)

  implicit none
  
  type(realization_type) :: realization
  
  call RealProcessTranConditions(realization)
 
end subroutine RealizationProcessConditions


! ************************************************************************** !
!
! RealProcessTranConditions: Sets up auxilliary data associated with 
!                            transport conditions
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine RealProcessTranConditions(realization)

  use Fileio_module
  use Reaction_module
  
  implicit none
  
  type(realization_type) :: realization
  
  
  PetscTruth :: found
  character(len=MAXSTRINGLENGTH) :: string
  type(tran_condition_type), pointer :: cur_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(tran_constraint_type), pointer :: cur_constraint, another_constraint
  
  ! check for duplicate constraint names
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
      another_constraint => cur_constraint%next
      ! now compare names
      found = PETSC_FALSE
      do
        if (.not.associated(another_constraint)) exit
        if (fiStringCompare(cur_constraint%name,another_constraint%name, &
            MAXWORDLENGTH)) then
          found = PETSC_TRUE
        endif
        another_constraint => another_constraint%next
      enddo
      if (found) then
        string = 'Duplicate transport constraints named "' // &
                 trim(cur_constraint%name) // '"'
        call printErrMsg(realization%option,string)
      endif
    cur_constraint => cur_constraint%next
  enddo
  
  ! initialize constraints
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
    call ReactionInitializeConstraint(realization%reaction, &
                                      cur_constraint%name, &
                                      cur_constraint%aqueous_species, &
                                      cur_constraint%minerals, &
                                      realization%option)
    cur_constraint => cur_constraint%next
  enddo
  
  ! tie constraints to couplers, if not already associated
  cur_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_condition)) exit
    cur_constraint_coupler => cur_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      if (.not.associated(cur_constraint_coupler%aqueous_species)) then
        cur_constraint => realization%transport_constraints%first
        do
          if (.not.associated(cur_constraint)) exit
          if (fiStringCompare(cur_constraint%name, &
                              cur_constraint_coupler%constraint_name, &
                              MAXWORDLENGTH)) then
            cur_constraint_coupler%aqueous_species => cur_constraint%aqueous_species
            cur_constraint_coupler%minerals => cur_constraint%minerals
            exit
          endif
          cur_constraint => cur_constraint%next
        enddo
        if (.not.associated(cur_constraint_coupler%aqueous_species)) then
          string = 'Transport constraint "' // &
                   trim(cur_constraint_coupler%constraint_name) // &
                   '" not found in input file constraints.'
          call printErrMsg(realization%option,string)
        endif
      endif
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    if (associated(cur_condition%constraint_coupler_list%next)) then ! more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    cur_condition => cur_condition%next
  enddo
 
  ! final details for setup
  cur_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_condition)) exit
    ! is the condition transient?
    if (associated(cur_condition%constraint_coupler_list%next)) then ! more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    ! set pointer to first constraint coupler
    cur_condition%cur_constraint_coupler => cur_condition%constraint_coupler_list
    
    cur_condition => cur_condition%next
  enddo

end subroutine RealProcessTranConditions

! ************************************************************************** !
!
! RealizationPrintCouplers: Print boundary and initial condition data
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine RealizationPrintCouplers(realization)

  use Coupler_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(coupler_type), pointer :: cur_coupler
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
 
  option => realization%option
  reaction => realization%reaction
 
  if (option%myrank /= 0) return
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      cur_coupler => cur_patch%initial_conditions%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo
     
      cur_coupler => cur_patch%boundary_conditions%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo
     
      cur_coupler => cur_patch%source_sinks%first
      do
        if (.not.associated(cur_coupler)) exit
        call RealizationPrintCoupler(cur_coupler,reaction,option)    
        cur_coupler => cur_coupler%next
      enddo

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo            
 
end subroutine RealizationPrintCouplers

! ************************************************************************** !
!
! RealizationPrintCoupler: Prints boundary and initial condition coupler 
! author: Glenn Hammond
! date: 10/28/08
!
! ************************************************************************** !
subroutine RealizationPrintCoupler(coupler,reaction,option)

  use Coupler_module
  use Reaction_module
  
  implicit none
  
  type(coupler_type) :: coupler
  type(option_type) :: option
  type(reaction_type), pointer :: reaction
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(region_type), pointer :: region
  type(tran_constraint_coupler_type), pointer :: constraint_coupler

   
98 format(40('=+'))
99 format(80('-'))
100 format(a)
  
  flow_condition => coupler%flow_condition
  tran_condition => coupler%tran_condition
  region => coupler%region

  write(option%fid_out,*)
  write(option%fid_out,98)
101 format(5x,'     Flow Condition: ',2x,a)
  if (associated(flow_condition)) write(option%fid_out,101) trim(flow_condition%name)
102 format(5x,'Transport Condition: ',2x,a)
  if (associated(tran_condition)) write(option%fid_out,102) trim(tran_condition%name)
103 format(5x,'             Region: ',2x,a)
  if (associated(region)) write(option%fid_out,103) trim(region%name)
  write(option%fid_out,98)
  
  if (associated(flow_condition)) then
!    write(option%fid_out,99)
    write(option%fid_out,100) '  Flow Parameters:'
104 format(a,f10.1)
    write(option%fid_out,104) '    pressure = ', &
      flow_condition%pressure%dataset%cur_value(1)
105 format(a,f5.1)
    if (associated(flow_condition%temperature)) then
      write(option%fid_out,105) '    temperature = ', &
        flow_condition%temperature%dataset%cur_value(1)
    endif
106 format(a,es12.4)
    if (associated(flow_condition%concentration)) then
      write(option%fid_out,106) '    concentration = ', &
        flow_condition%concentration%dataset%cur_value(1)
    endif
  endif
  if (associated(tran_condition)) then
    constraint_coupler => tran_condition%cur_constraint_coupler
    write(option%fid_out,99) 
    write(option%fid_out,100) '  Transport Parameters:'
    if (associated(reaction)) then
      call RPrintConstraint(constraint_coupler,reaction, &
                            option)
    endif
  endif
 
end subroutine RealizationPrintCoupler

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
  
  logical :: force_update_flag = PETSC_FALSE
  
  ! must update conditions first
  call ConditionUpdate(realization%flow_conditions,realization%option, &
                       realization%option%time)
  call TranConditionUpdate(realization%transport_conditions, &
                           realization%option, &
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
  
  if (realization%option%nflowdof > 0) &
    call RealizAssignFlowInitCond(realization)
  if (realization%option%ntrandof > 0) &
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

      grid => cur_patch%grid

      select case(option%iflowmode)
        case(RICHARDS_MODE)
        case(THC_MODE)
        case(MPH_MODE)
        case(IMS_MODE)
!            call pflow_mphase_setupini(realization)
      end select 

      ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
      call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
      
      xx_p = -999.d0
      
      initial_condition => cur_patch%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit

        if (.not.associated(initial_condition%flow_aux_real_var)) then
          do icell=1,initial_condition%region%num_cells
            local_id = initial_condition%region%cell_ids(icell)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%nflowdof
            ibegin = iend-option%nflowdof+1
            if (associated(cur_patch%imat)) then
              if (cur_patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 0.d0
                iphase_loc_p(ghosted_id) = 0
                cycle
              endif
            endif
            do idof = 1, option%nflowdof
              xx_p(ibegin+idof-1) = &
                initial_condition%flow_condition%sub_condition_ptr(idof)%ptr%dataset%cur_value(1)
            enddo
            iphase_loc_p(ghosted_id)=initial_condition%flow_condition%iphase
          enddo
        else
          do iconn=1,initial_condition%connection_set%num_connections
            local_id = initial_condition%connection_set%id_dn(iconn)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%nflowdof
            ibegin = iend-option%nflowdof+1
            if (associated(cur_patch%imat)) then
              if (cur_patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 0.d0
                iphase_loc_p(ghosted_id) = 0
                cycle
              endif
            endif
            xx_p(ibegin:iend) = &
              initial_condition%flow_aux_real_var(1:option%nflowdof,iconn)
            iphase_loc_p(ghosted_id)=initial_condition%flow_aux_int_var(1,iconn)
          enddo
        endif
        initial_condition => initial_condition%next
      enddo
      
      call GridVecRestoreArrayF90(grid,field%flow_xx,xx_p, ierr)
      call GridVecRestoreArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
        
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx,field%flow_xx_loc,NFLOWDOF)  
  call VecCopy(field%flow_xx, field%flow_yy, ierr)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)  
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_old_loc,ONEDOF)
  
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
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof, isub_condition
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
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: aux_vars(:)
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      aux_vars => cur_patch%aux%RT%aux_vars

      ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid, field%tran_xx,xx_p, ierr); CHKERRQ(ierr)
      
      xx_p = -999.d0
      
      initial_condition => cur_patch%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit

        if (.not.associated(initial_condition%tran_aux_real_var)) then
          do icell=1,initial_condition%region%num_cells
            local_id = initial_condition%region%cell_ids(icell)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%ntrandof
            ibegin = iend-option%ntrandof+1
            if (associated(cur_patch%imat)) then
              if (cur_patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 1.d-200
                cycle
              endif
            endif
            do idof = 1, option%ntrandof ! primary aqueous concentrations
              xx_p(ibegin+idof-1) = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                  aqueous_species%basis_molarity(idof) / &
                  aux_vars(ghosted_id)%den(1)*1000 ! convert molarity -> molality
            enddo
            if (associated(initial_condition%tran_condition%cur_constraint_coupler%minerals)) then
              do idof = 1, reaction%nmnrl
                aux_vars(ghosted_id)%mnrl_volfrac(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_mol_frac(idof)
              enddo
            endif
          enddo
        else
          do iconn=1,initial_condition%connection_set%num_connections
            local_id = initial_condition%connection_set%id_dn(iconn)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%ntrandof
            ibegin = iend-option%ntrandof+1
            if (associated(cur_patch%imat)) then
              if (cur_patch%imat(ghosted_id) <= 0) then
                xx_p(ibegin:iend) = 1.d-200
                cycle
              endif
            endif
            xx_p(ibegin:iend) = &
              initial_condition%tran_aux_real_var(1:option%ntrandof,iconn) / &
              aux_vars(ghosted_id)%den(1)*1000 ! convert molarity -> molality
              ! minerals 
            if (associated(initial_condition%tran_condition%cur_constraint_coupler%minerals)) then
              do idof = 1, reaction%nmnrl
                aux_vars(ghosted_id)%mnrl_volfrac(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_mol_frac(idof)
              enddo
            endif
          enddo
        endif
        initial_condition => initial_condition%next
      enddo
      
      call GridVecRestoreArrayF90(grid,field%tran_xx,xx_p, ierr)

      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx,field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr)

end subroutine RealizAssignTransportInitCond

! ************************************************************************** !
!
! RealizBridgeFlowAndTransport: Maps auxilliary data (e.g. density) from flow
!                               to transport
! author: Glenn Hammond
! date: 09/03/08
!
! ************************************************************************** !
subroutine RealizBridgeFlowAndTransport(realization)
  
  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  option => realization%option

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchBridgeFlowAndTransport(cur_patch,option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RealizBridgeFlowAndTransport

! ************************************************************************** !
!
! RealizationRevertFlowParameters: Assigns initial porosity/perms to vecs
! author: Glenn Hammond
! date: 05/09/08
!
! ************************************************************************** !
subroutine RealizationRevertFlowParameters(realization)

  use Option_module
  use Field_module
  use Discretization_module

  implicit none
  
  type(realization_type) :: realization
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization

  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                           field%perm_xx_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                           field%perm_yy_loc,ONEDOF)  
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                           field%perm_zz_loc,ONEDOF)   
  endif   
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                         field%porosity_loc,ONEDOF)
                           
end subroutine RealizationRevertFlowParameters

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
  type(flow_condition_type), pointer :: cur_flow_condition
  type(tran_condition_type), pointer :: cur_tran_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(waypoint_type), pointer :: waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition

  option => realization%option
  waypoint_list => realization%waypoints

  cur_flow_condition => realization%flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    if (cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        itime = 1
        if (sub_condition%dataset%max_time_index == 1 .and. &
            sub_condition%dataset%times(itime) > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = sub_condition%dataset%times(itime)
          waypoint%update_bcs = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
          exit
        endif
      enddo
    endif
    cur_flow_condition => cur_flow_condition%next
  enddo
      
  cur_tran_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_tran_condition)) exit
    if (cur_tran_condition%sync_time_with_update .and. &
        cur_tran_condition%is_transient) then
      cur_constraint_coupler => cur_tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        if (cur_constraint_coupler%time > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = cur_constraint_coupler%time
          waypoint%update_bcs = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
    cur_tran_condition => cur_tran_condition%next
  enddo
      
end subroutine RealizationAddWaypointsToList

! ************************************************************************** !
!
! RealizationGetDataset: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationGetDataset(realization,vec,ivar,isubvar)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchGetDataset(cur_patch,realization%field,realization%option, &
                           vec,ivar,isubvar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RealizationGetDataset

! ************************************************************************** !
!
! RealizGetDatasetValueAtCell: Extracts variables indexed by ivar and isubvar from a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
function RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetDatasetValueAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: ghosted_id
  
  PetscReal :: value
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      value = PatchGetDatasetValueAtCell(cur_patch,realization%field, &
                                         realization%option, &
                                         ivar,isubvar,ghosted_id)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  RealizGetDatasetValueAtCell = value

end function RealizGetDatasetValueAtCell

! ************************************************************************** !
!
! RealizationSetDataset: Sets variables indexed by ivar and isubvar in a 
!                        realization
! author: Glenn Hammond
! date: 09/12/08
!
! ************************************************************************** !
subroutine RealizationSetDataset(realization,vec,ivar,isubvar)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchSetDataset(cur_patch,realization%field,realization%option, &
                           vec,ivar,isubvar)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine RealizationSetDataset

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
  call TranConditionDestroyList(realization%transport_conditions)
  call TranConstraintDestroyList(realization%transport_constraints)

  call LevelDestroyList(realization%level_list)

  if (associated(realization%debug)) deallocate(realization%debug)
  nullify(realization%debug)
  
  call FluidPropertyDestroy(realization%fluid_properties)
  
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
