module Realization_module

  use Option_module
  use Input_module
  use Region_module
  use Condition_module
  use Material_module
  use Saturation_Function_module
  use Fluid_module
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

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    type(level_list_type), pointer :: level_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(input_type), pointer :: input
    type(field_type), pointer :: field
    type(pflow_debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option

    type(region_list_type), pointer :: regions
    type(condition_list_type), pointer :: flow_conditions
    type(tran_condition_list_type), pointer :: transport_conditions
    type(tran_constraint_list_type), pointer :: transport_constraints
    
    type(reaction_type), pointer :: reaction
    
    type(material_property_type), pointer :: material_properties
    type(material_property_ptr_type), pointer :: material_property_array(:)
    type(fluid_property_type), pointer :: fluid_properties
    type(fluid_property_type), pointer :: fluid_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    
    type(waypoint_list_type), pointer :: waypoints
    
  end type realization_type

  interface RealizationCreate
    module procedure RealizationCreate1
    module procedure RealizationCreate2
  end interface
  
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
            RealizationAddObservation, &
            RealizAssignInitialConditions, &
            RealizAssignFlowInitCond, &
            RealizAssignTransportInitCond, &
            RealizAssignUniformVelocity, &
            RealizationRevertFlowParameters, &
            RealizationGetDataset, &
            RealizGetDatasetValueAtCell, &
            RealizationSetDataset, &
            RealizationPrintCouplers, &
            RealizationInitConstraints, &
            RealProcessMatPropAndSatFunc, &
            RealProcessFluidProperties
            
            
contains
  
! ************************************************************************** !
!
! RealizationCreate1: Allocates and initializes a new Realization object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function RealizationCreate1()

  implicit none
  
  type(realization_type), pointer :: RealizationCreate1
  
  type(realization_type), pointer :: realization
  type(option_type), pointer :: option
  
  nullify(option)
  RealizationCreate1 => RealizationCreate2(option)
  
end function RealizationCreate1  

! ************************************************************************** !
!
! RealizationCreate2: Allocates and initializes a new Realization object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function RealizationCreate2(option)

  implicit none
  
  type(option_type), pointer :: option
  
  type(realization_type), pointer :: RealizationCreate2
  
  type(realization_type), pointer :: realization
  
  allocate(realization)
  realization%discretization => DiscretizationCreate()
  if (associated(option)) then
    realization%option => option
  else
    realization%option => OptionCreate()
  endif
  nullify(realization%input)
  realization%field => FieldCreate()
  realization%debug => DebugCreatePflow()
  realization%output_option => OutputOptionCreate()

  realization%level_list => LevelCreateList()

  allocate(realization%regions)
  call RegionInitList(realization%regions)

  allocate(realization%flow_conditions)
  call FlowConditionInitList(realization%flow_conditions)
  allocate(realization%transport_conditions)
  call TranConditionInitList(realization%transport_conditions)
  allocate(realization%transport_constraints)
  call TranConstraintInitList(realization%transport_constraints)

  nullify(realization%material_properties)
  nullify(realization%material_property_array)
  nullify(realization%fluid_properties)
  nullify(realization%fluid_property_array)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  
  nullify(realization%reaction)
  
  RealizationCreate2 => realization
  
end function RealizationCreate2 

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
  
  option%ivar_centering = CELL_CENTERED
  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity0, &
                                  GLOBAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%volume)

  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%work)
  
  ! 1 degree of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity_loc, &
                                  LOCAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%tor_loc)

  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%work_loc)
  
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

    if(option%use_samr) then
       option%ivar_centering = SIDE_CENTERED
       call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_face_fluxes, &
                                    GLOBAL,option)
      option%ivar_centering = CELL_CENTERED
    endif
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

    ! ndof degrees of freedom, local
    call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                    LOCAL,option)
                                    
    if (realization%reaction%use_log_formulation) then
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_log_xx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx_loc, &
                                         field%tran_work_loc)
    endif
    
    if(option%use_samr) then
       option%ivar_centering = SIDE_CENTERED
        call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_face_fluxes, &
                                    GLOBAL,option)
      option%ivar_centering = CELL_CENTERED
    endif

  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
    
      grid => discretization%grid

      ! set up nG2L, NL2G, etc.
      call GridMapIndices(grid)
      call GridComputeSpacing(grid,option)
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
! RealizationAddObservation: Adds a copy of a observation object to a list
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizationAddObservation(realization,observation)

  use Observation_module

  implicit none
  
  type(realization_type) :: realization
  type(observation_type), pointer :: observation
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type(observation_type), pointer :: new_observation
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      new_observation => ObservationCreate(observation)
      call ObservationAddToList(new_observation, &
                                 cur_patch%observation)
      nullify(new_observation)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  call ObservationDestroy(observation)
 
end subroutine RealizationAddObservation

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
                                realization%material_properties, &
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
  
  if (realization%option%ntrandof > 0) then
    call RealProcessTranConditions(realization)
  endif
 
end subroutine RealizationProcessConditions

! ************************************************************************** !
!
! RealProcessMatPropAndSatFunc: Sets up linkeage between material properties
!                               and saturation function and auxilliary arrays
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
subroutine RealProcessMatPropAndSatFunc(realization)

  use String_module
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscTruth :: found
  PetscInt :: i
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  type(saturation_function_type), pointer :: cur_saturation_function
  
  option => realization%option
  
  ! organize lists
  call MaterialPropConvertListToArray(realization%material_properties, &
                                      realization%material_property_array)
  call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                     realization%saturation_function_array, &
                                     option) 
    
  cur_material_property => realization%material_properties                            
  do                                      
    if (.not.associated(cur_material_property)) exit
    found = PETSC_FALSE
    cur_saturation_function => realization%saturation_functions
    do 
      if (.not.associated(cur_saturation_function)) exit
      if (StringCompare(cur_material_property%saturation_function_name, &
                        cur_saturation_function%name,MAXWORDLENGTH)) then
        found = PETSC_TRUE
        cur_material_property%saturation_function_id = &
          cur_saturation_function%id
        exit
      endif
      cur_saturation_function => cur_saturation_function%next
    enddo
    if (.not.found) then
      option%io_buffer = 'Saturation function "' // &
               trim(cur_material_property%saturation_function_name) // &
               '" in material property "' // &
               trim(cur_material_property%name) // &
               '" not found among available saturation functions.'
      call printErrMsg(realization%option)    
    endif
    cur_material_property => cur_material_property%next
  enddo
  
end subroutine RealProcessMatPropAndSatFunc

! ************************************************************************** !
!
! RealProcessFluidProperties: Sets up linkeage with fluid properties
! author: Glenn Hammond
! date: 01/21/09
!
! ************************************************************************** !
subroutine RealProcessFluidProperties(realization)
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscTruth :: found
  type(option_type), pointer :: option
  type(fluid_property_type), pointer :: cur_fluid_property
  
  option => realization%option
  
  found = PETSC_FALSE
  cur_fluid_property => realization%fluid_properties                            
  do                                      
    if (.not.associated(cur_fluid_property)) exit
    found = PETSC_TRUE
    select case(trim(cur_fluid_property%phase_name))
      case('LIQUID_PHASE')
        cur_fluid_property%phase_id = LIQUID_PHASE
      case('GAS_PHASE')
        cur_fluid_property%phase_id = GAS_PHASE
      case default
        cur_fluid_property%phase_id = LIQUID_PHASE
    end select
    cur_fluid_property => cur_fluid_property%next
  enddo
  
  if (option%ntrandof > 0 .and. .not.found) then
    option%io_buffer = 'A fluid property must be present in input file' // &
                       ' for solute transport'
  endif
  
end subroutine RealProcessFluidProperties

! ************************************************************************** !
!
! RealProcessTranConditions: Sets up auxilliary data associated with 
!                            transport conditions
! author: Glenn Hammond
! date: 10/14/08
!
! ************************************************************************** !
subroutine RealProcessTranConditions(realization)

  use String_module
  use Reaction_module
  
  implicit none
  
  type(realization_type) :: realization
  
  
  PetscTruth :: found
  type(option_type), pointer :: option
  type(tran_condition_type), pointer :: cur_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(tran_constraint_type), pointer :: cur_constraint, another_constraint
  
  option => realization%option
  
  ! check for duplicate constraint names
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
      another_constraint => cur_constraint%next
      ! now compare names
      found = PETSC_FALSE
      do
        if (.not.associated(another_constraint)) exit
        if (StringCompare(cur_constraint%name,another_constraint%name, &
            MAXWORDLENGTH)) then
          found = PETSC_TRUE
        endif
        another_constraint => another_constraint%next
      enddo
      if (found) then
        option%io_buffer = 'Duplicate transport constraints named "' // &
                 trim(cur_constraint%name) // '"'
        call printErrMsg(realization%option)
      endif
    cur_constraint => cur_constraint%next
  enddo
  
  ! initialize constraints
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
    call ReactionProcessConstraint(realization%reaction, &
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
          if (StringCompare(cur_constraint%name, &
                             cur_constraint_coupler%constraint_name, &
                             MAXWORDLENGTH)) then
            cur_constraint_coupler%aqueous_species => cur_constraint%aqueous_species
            cur_constraint_coupler%minerals => cur_constraint%minerals
            exit
          endif
          cur_constraint => cur_constraint%next
        enddo
        if (.not.associated(cur_constraint_coupler%aqueous_species)) then
          option%io_buffer = 'Transport constraint "' // &
                   trim(cur_constraint_coupler%constraint_name) // &
                   '" not found in input file constraints.'
          call printErrMsg(realization%option)
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
! RealizationInitConstraints: Initializes constraint concentrations
! author: Glenn Hammond
! date: 12/04/08
!
! ************************************************************************** !
subroutine RealizationInitConstraints(realization)

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
      call PatchInitConstraints(cur_patch,realization%reaction, &
                                realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo            
 
end subroutine RealizationInitConstraints

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
 
  if (.not.OptionPrintToFile(option)) return
  
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


  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      string = 'Initial Condition'
    case(BOUNDARY_COUPLER_TYPE)
      string = 'Boundary Condition'
    case(SRC_SINK_COUPLER_TYPE)
      string = 'Source Sink'
  end select
  write(option%fid_out,'(/,2x,a,/)') trim(string)

  write(option%fid_out,99)
101 format(5x,'     Flow Condition: ',2x,a)
  if (associated(flow_condition)) write(option%fid_out,101) trim(flow_condition%name)
102 format(5x,'Transport Condition: ',2x,a)
  if (associated(tran_condition)) write(option%fid_out,102) trim(tran_condition%name)
103 format(5x,'             Region: ',2x,a)
  if (associated(region)) write(option%fid_out,103) trim(region%name)
  write(option%fid_out,99)
  
  if (associated(flow_condition)) then
    call FlowConditionPrint(flow_condition,option)
  endif
  if (associated(tran_condition)) then
    constraint_coupler => tran_condition%cur_constraint_coupler
    write(option%fid_out,'(/,2x,''Transport Condition: '',a)') &
      trim(tran_condition%name)
    if (associated(reaction)) then
      call ReactionPrintConstraint(constraint_coupler,reaction,option)
      write(option%fid_out,'(/)')
      write(option%fid_out,99)
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
      call PatchInitAllCouplerAuxVars(cur_patch,realization%reaction, &
                                      realization%option)
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
  PetscTruth :: force_update_flag

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
  
  PetscTruth :: force_update_flag = PETSC_FALSE
  
  ! must update conditions first
  call FlowConditionUpdate(realization%flow_conditions,realization%option, &
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

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
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
  use Global_Aux_module
  use Reaction_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
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
  type(reactive_transport_auxvar_type), pointer :: rt_aux_vars(:)
  type(global_auxvar_type), pointer :: global_aux_vars(:)
 
  PetscInt :: iphase
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  
  iphase = 1
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit

      grid => cur_patch%grid
      rt_aux_vars => cur_patch%aux%RT%aux_vars
      global_aux_vars => cur_patch%aux%Global%aux_vars

      ! assign initial conditions values to domain
      call GridVecGetArrayF90(grid, field%tran_xx,xx_p, ierr); CHKERRQ(ierr)
      
      xx_p = -999.d0
      
      initial_condition => cur_patch%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit

!        if (.not.associated(initial_condition%tran_aux_real_var)) then
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
            if (.not.option%use_isothermal) then
              if (icell == 1) then
                call ReactionEquilibrateConstraint(rt_aux_vars(ghosted_id), &
                                                   global_aux_vars(ghosted_id), &
                                                   reaction, &
                    initial_condition%tran_condition%cur_constraint_coupler%constraint_name, &
                    initial_condition%tran_condition%cur_constraint_coupler%aqueous_species, &
                    initial_condition%tran_condition%cur_constraint_coupler%num_iterations, &
                                                   PETSC_TRUE,option)
              else
                call RTAuxVarCopy(rt_aux_vars(ghosted_id), &
                  rt_aux_vars(grid%nL2G(initial_condition%region%cell_ids(icell-1))), &
                  option)
                call ReactionEquilibrateConstraint(rt_aux_vars(ghosted_id), &
                                                   global_aux_vars(ghosted_id), &
                                                   reaction, &
                    initial_condition%tran_condition%cur_constraint_coupler%constraint_name, &
                    initial_condition%tran_condition%cur_constraint_coupler%aqueous_species, &
                    initial_condition%tran_condition%cur_constraint_coupler%num_iterations, &
                                                   PETSC_FALSE,option)
              endif
            endif
            do idof = 1, option%ntrandof ! primary aqueous concentrations
              xx_p(ibegin+idof-1) = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                  aqueous_species%basis_molarity(idof) / &
                  global_aux_vars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
            enddo
            if (associated(initial_condition%tran_condition%cur_constraint_coupler%minerals)) then
              do idof = 1, reaction%nkinmnrl
                rt_aux_vars(ghosted_id)%mnrl_volfrac0(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_vol_frac(idof)
                rt_aux_vars(ghosted_id)%mnrl_volfrac(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_vol_frac(idof)
                rt_aux_vars(ghosted_id)%mnrl_area0(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_area(idof)
                rt_aux_vars(ghosted_id)%mnrl_area(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                    minerals%basis_area(idof)
              enddo
            endif
            ! this is for the multi-rate surface complexation model
            if (reaction%kinmr_nrate > 0) then
              rt_aux_vars(ghosted_id)%kinmr_total_sorb = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                rt_auxvar%kinmr_total_sorb
            endif
          enddo
!        endif
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
! RealizationAddWaypointsToList: Creates waypoints associated with source/sinks
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
  
  type(waypoint_list_type), pointer :: waypoint_list
  type(flow_condition_type), pointer :: cur_flow_condition
  type(tran_condition_type), pointer :: cur_tran_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time

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
          pause
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
          pause
          call WaypointInsertInList(waypoint,waypoint_list)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
    cur_tran_condition => cur_tran_condition%next
  enddo

  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_output = realization%output_option%print_final
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) final_time = cur_waypoint%time

  ! add waypoints for periodic output
  if (realization%output_option%periodic_output_time_incr > 0.d0 .or. &
      realization%output_option%periodic_tr_output_time_incr > 0.d0) then

    if (realization%output_option%periodic_output_time_incr > 0.d0) then
      ! standard output
      temp_real = 0.d0
      do
        temp_real = temp_real + realization%output_option%periodic_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,realization%waypoints)
      enddo
    endif
    
    if (realization%output_option%periodic_tr_output_time_incr > 0.d0) then
      ! transient observation output
      temp_real = 0.d0
      do
        temp_real = temp_real + realization%output_option%periodic_tr_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_tr_output = PETSC_TRUE 
        call WaypointInsertInList(waypoint,realization%waypoints)
      enddo
    endif

  endif

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
subroutine RealizationSetDataset(realization,vec,vec_format,ivar,isubvar)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: vec_format
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
                           vec,vec_format,ivar,isubvar)
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
!  call OptionDestroy(realization%option) !geh it will be destroy externally
  call RegionDestroyList(realization%regions)
  
  call FlowConditionDestroyList(realization%flow_conditions)
  call TranConditionDestroyList(realization%transport_conditions)
  call TranConstraintDestroyList(realization%transport_constraints)

  call LevelDestroyList(realization%level_list)

  if (associated(realization%debug)) deallocate(realization%debug)
  nullify(realization%debug)
  
  if (associated(realization%fluid_property_array)) &
    deallocate(realization%fluid_property_array)
  nullify(realization%fluid_property_array)
  call FluidPropertyDestroy(realization%fluid_properties)
  
  if (associated(realization%material_property_array)) &
    deallocate(realization%material_property_array)
  nullify(realization%material_property_array)
  call MaterialPropertyDestroy(realization%material_properties)
  
  if (associated(realization%saturation_function_array)) &
    deallocate(realization%saturation_function_array)
  nullify(realization%saturation_function_array)
  call SaturationFunctionDestroy(realization%saturation_functions)
  
end subroutine RealizationDestroy
  
end module Realization_module
