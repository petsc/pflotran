module Realization_module

  use Option_module
  use Input_module
  use Region_module
  use Condition_module
  use Material_module
#ifdef SUBCONTINUUM_MODEL
  use Subcontinuum_module
#endif
  use Saturation_Function_module
  use Dataset_module
  use Fluid_module
  use Discretization_module
  use Field_module
  use Debug_module
  use Velocity_module
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
#ifdef SUBCONTINUUM_MODEL
    type(subcontinuum_property_type), pointer :: subcontinuum_properties
    type(subcontinuum_property_ptr_type), pointer ::  &
                               subcontinuum_property_array(:)
#endif
    type(fluid_property_type), pointer :: fluid_properties
    type(fluid_property_type), pointer :: fluid_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    type(dataset_type), pointer :: datasets
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    
    type(velocity_dataset_type), pointer :: velocity_dataset
    
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
            RealizUpdateUniformVelocity, &
            RealizationRevertFlowParameters, &
            RealizationGetDataset, &
            RealizGetDatasetValueAtCell, &
            RealizationSetDataset, &
            RealizationPrintCouplers, &
            RealizationInitConstraints, &
            RealProcessMatPropAndSatFunc, &
            RealProcessSubcontinuumProp, &
            RealProcessFluidProperties, &
            RealizationUpdateProperties, &
            RealizationCountCells, &
            RealizationPrintGridStatistics, &
            RealizationSetUpBC4Faces, &
            RealizatonPassFieldPtrToPatches
            
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
#ifdef SUBCONTINUUM_MODEL
  nullify(realization%subcontinuum_properties)
  nullify(realization%subcontinuum_property_array)
#endif
  nullify(realization%fluid_properties)
  nullify(realization%fluid_property_array)
  nullify(realization%saturation_functions)
  nullify(realization%saturation_function_array)
  nullify(realization%datasets)
  nullify(realization%velocity_dataset)
  
  nullify(realization%reaction)

  nullify(realization%patch)    
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
  use Unstructured_Grid_module, only : UGridMapIndices
  use AMR_Grid_module
  use MFD_module
  use Coupler_module
  
  implicit none
  
  type(realization_type) :: realization
  
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  PetscErrorCode :: ierr
  PetscInt, allocatable :: int_tmp(:)
  PetscInt :: test,j
  PetscOffset :: i_da
  PetscReal, pointer :: real_tmp(:)
  
  option => realization%option
  field => realization%field
 
 
  discretization => realization%discretization
  
  call DiscretizationCreateDMs(discretization,option)

  option%ivar_centering = CELL_CENTERED
  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity0, &
                                  GLOBAL,option)
 
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%tortuosity0)
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%volume)

  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%work)
  ! temporary for samr testing
  call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                     field%work_samr)
  
  ! 1 degree of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field%porosity_loc, &
                                  LOCAL,option)
  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%tortuosity_loc)

  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%work_loc)
  
  ! temporary for samr testing
  call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                     field%work_samr_loc)
  
  if (option%nflowdof > 0) then

    ! 1-dof global  
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_xx)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_yy)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_zz)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_xz)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_xy)
    call DiscretizationDuplicateVector(discretization,field%porosity0, &
                                       field%perm0_yz)
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
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_xz_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_xy_loc)
    call DiscretizationDuplicateVector(discretization,field%porosity_loc, &
                                       field%perm_yz_loc)

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
    if (option%reactive_transport_coupling == GLOBAL_IMPLICIT) then
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
    else ! operator splitting
      ! ndof degrees of freedom, global
      ! create the 1 dof vector for solving the individual linear systems
      call DiscretizationCreateVector(discretization,ONEDOF,field%tran_rhs_coef, &
                                      GLOBAL,option)
      ! create the ntran dof vector for storage of the solution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_rhs)

      ! ndof degrees of freedom, local
      ! again, just for storage of the current colution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                      LOCAL,option)
      if(option%use_samr) then
         option%ivar_centering = SIDE_CENTERED
         call DiscretizationCreateVector(discretization,ONEDOF,field%tran_face_fluxes, &
                                      GLOBAL,option)
         option%ivar_centering = CELL_CENTERED
      endif 

    endif
    
  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID, STRUCTURED_GRID_MIMETIC)
      grid => discretization%grid
      ! set up nG2L, nL2G, etc.
      call GridMapIndices(grid, discretization%dm_1dof%sgdm)
      call GridComputeSpacing(grid,option)
      call GridComputeCoordinates(grid,discretization%origin,option)
      call GridComputeVolumes(grid,field%volume,option)
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option)
      if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
          call GridComputeCell2FaceConnectivity(grid, discretization%MFD, option)
      end if
    case(UNSTRUCTURED_GRID)
      grid => discretization%grid
      ! set up nG2L, NL2G, etc.
      call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
                           grid%nG2L,grid%nL2G,grid%nL2A,grid%nG2A)
      call GridComputeCoordinates(grid,discretization%origin,option, & 
                                   discretization%dm_1dof%ugdm)  !sp 
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option, discretization%dm_1dof%ugdm) !sp 
      call GridComputeVolumes(grid,field%volume,option)
    case(AMR_GRID)
       call AMRGridComputeGeometryInformation(discretization%amrgrid, &
                                              discretization%origin, &
                                              field,option)
  end select 
 
  ! Vectors with face degrees of freedom
#ifdef DASVYAT
   if (discretization%itype==STRUCTURED_GRID_MIMETIC) then

     if (option%nflowdof > 0) then
   
       call VecCreateMPI(option%mycomm,grid%nlmax_faces*option%nflowdof, &
                    PETSC_DETERMINE,field%flow_xx_faces,ierr)



       call VecSetBlockSize(field%flow_xx_faces,option%nflowdof,ierr)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_r_faces)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_dxx_faces)

       call DiscretizationDuplicateVector(discretization, field%flow_xx_faces, &
                                                        field%flow_yy_faces)


       call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*option%nflowdof, field%flow_xx_loc_faces, ierr)
       call VecSetBlockSize(field%flow_xx_loc_faces,option%nflowdof,ierr)

!       call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*option%nflowdof, field%flow_r_loc_faces, ierr)
!       call VecSetBlockSize(field%flow_r_loc_faces,NFLOWDOF,ierr)

!       call VecCreateSeq(PETSC_COMM_SELF, grid%ngmax_faces*option%nflowdof, field%flow_bc_loc_faces, ierr)
!       call VecSetBlockSize(field%flow_bc_loc_faces,NFLOWDOF,ierr)

        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%flow_r_loc_faces) 

        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%flow_bc_loc_faces)
   
        call DiscretizationDuplicateVector(discretization, field%flow_xx_loc_faces, field%work_loc_faces)

!       call VecGetArrayF90(field%volume, real_tmp, ierr)
!       call VecRestoreArrayF90(field%volume, real_tmp, ierr)
 


     end if

     call GridComputeGlobalCell2FaceConnectivity(grid, discretization%MFD, NFLOWDOF, option)
  
   end if
#endif
 
  ! initialize to -999.d0 for check later that verifies all values 
  ! have been set
  call VecSet(field%porosity0,-999.d0,ierr)
       

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
! RealizatonPassFieldPtrToPatches: Sets patch%field => realization%field
! author: Glenn Hammond
! date: 01/12/11
!
! ************************************************************************** !
subroutine RealizatonPassFieldPtrToPatches(realization)

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
      cur_patch%field => realization%field
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine RealizatonPassFieldPtrToPatches

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
!      call PatchProcessCouplers(cur_patch,realization%flow_conditions, &
!                                realization%transport_conditions, &
!                                realization%material_properties, &
!                                realization%subcontinuum_properties, &
!                                realization%option)
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
!                               and saturation function, auxilliary arrays
!                               and datasets
! author: Glenn Hammond
! date: 01/21/09, 01/12/11
!
! ************************************************************************** !
subroutine RealProcessMatPropAndSatFunc(realization)

  use String_module
  
  implicit none
  
  type(realization_type) :: realization
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  character(len=MAXSTRINGLENGTH) :: string
  
  option => realization%option
  
  ! organize lists
  call MaterialPropConvertListToArray(realization%material_properties, &
                                      realization%material_property_array, &
                                      option)
  call SaturatFuncConvertListToArray(realization%saturation_functions, &
                                     realization%saturation_function_array, &
                                     option) 
    
  cur_material_property => realization%material_properties                            
  do                                      
    if (.not.associated(cur_material_property)) exit

    ! obtain saturation function id
    cur_material_property%saturation_function_id = &
      SaturationFunctionGetID(realization%saturation_functions, &
                              cur_material_property%saturation_function_name, &
                              cur_material_property%name,option)
    
    ! if named, link dataset to property
    if (.not.StringNull(cur_material_property%porosity_dataset_name)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),POROSITY'
      cur_material_property%porosity_dataset => &
        DatasetGetPointer(realization%datasets, &
                          cur_material_property%porosity_dataset_name, &
                          string,option)
    endif
    if (.not.StringNull(cur_material_property%permeability_dataset_name)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY'
      cur_material_property%permeability_dataset => &
        DatasetGetPointer(realization%datasets, &
                          cur_material_property%permeability_dataset_name, &
                          string,option)
    endif
    
    cur_material_property => cur_material_property%next
  enddo
  
  
end subroutine RealProcessMatPropAndSatFunc


! *********************************************************************** !
!
! RealProcessSubcontinuumProp: Sets up linkeage between material properties
!                         and subcontinuum properties and auxilliary arrays
! author: Jitendra Kumar 
! date: 10/07/2010 
!
! *********************************************************************** !
subroutine RealProcessSubcontinuumProp(realization)

  use String_module
  
  implicit none
  type(realization_type) :: realization
#ifdef SUBCONTINUUM_MODEL 
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  type(subcontinuum_property_type), pointer :: cur_subcontinuum_property
  
  option => realization%option
  
  ! organize lists
  call MaterialPropConvertListToArray(realization%material_properties, &
                                      realization%material_property_array)
  call SubcontinuumPropConvertListToArray(realization%subcontinuum_properties, &
                                     realization%subcontinuum_properties_array, &
                                     option) 
    
  cur_material_property => realization%material_properties                            
  do                                      
    if (.not.associated(cur_material_property)) exit
    found = PETSC_FALSE
    cur_subcontinuum_property => realization%subcontinuum_properties
    do 
      if (.not.associated(cur_subcontinuum_property)) exit
      do i=1,cur_material_property%subcontinuum_type_count
        if (StringCompare(cur_material_property%subcontinuum_type_name(i), &
                        cur_subcontinuum_property%name,MAXWORDLENGTH)) then
          found = PETSC_TRUE
          cur_material_property%subcontinuum_property_id(i) = &
            cur_subcontinuum_property%id
          exit
        endif
        ! START TODO: check if this error check block is at right place. might have
        ! to push it to upper do loop: Jitu 10/07/2010 
        if (.not.found) then
          option%io_buffer = 'Saturation function "' // &
               trim(cur_material_property%subcontinuum_type_name(i)) // &
               '" in material property "' // &
               trim(cur_material_property%name) // &
               '" not found among available subcontinuum types.'
          call printErrMsg(realization%option)    
        endif
        ! END TODO
      enddo  
      cur_subcontinuum_property => cur_subcontinuum_property%next
    enddo
    cur_material_property => cur_material_property%next
  enddo
#endif 
end subroutine RealProcessSubcontinuumProp

! ************************************************************************** !
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
  
  PetscBool :: found
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
  
  
  PetscBool :: found
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
                                   cur_constraint%surface_complexes, &
                                   cur_constraint%colloids, &
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
            cur_constraint_coupler%surface_complexes => cur_constraint%surface_complexes
            cur_constraint_coupler%colloids => cur_constraint%colloids
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
  PetscBool :: force_update_flag

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
  
  PetscBool :: force_update_flag = PETSC_FALSE
  
  ! must update conditions first
  call FlowConditionUpdate(realization%flow_conditions,realization%option, &
                           realization%option%time)
  call TranConditionUpdate(realization%transport_conditions, &
                           realization%option, &
                           realization%option%time)
  call RealizUpdateAllCouplerAuxVars(realization,force_update_flag)
  if (associated(realization%velocity_dataset)) then
    call VelocityDatasetUpdate(realization%option,realization%option%time, &
                               realization%velocity_dataset)
    call RealizUpdateUniformVelocity(realization)
  endif
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
#ifdef DASVYAT
  use MFD_module, only :MFDInitializeMassMatrices
#endif
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscInt :: icell, iconn, idof, iface
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
        case(FLASH2_MODE)
!            call pflow_mphase_setupini(realization)
      end select 

      ! assign initial conditions values to domain
      if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
        call GridVecGetArrayF90(grid,field%flow_xx_faces, xx_p, ierr); CHKERRQ(ierr)
      else
        call GridVecGetArrayF90(grid,field%flow_xx,xx_p, ierr); CHKERRQ(ierr)
      end if
      call GridVecGetArrayF90(grid,field%iphas_loc,iphase_loc_p,ierr)
      
      xx_p = -999.d0
      
      initial_condition => cur_patch%initial_conditions%first
      do
      
        if (.not.associated(initial_condition)) exit

        if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
#ifdef DASVYAT
           if (.not.associated(initial_condition%flow_aux_real_var)) then
             do icell=1,initial_condition%region%num_cells
               local_id = initial_condition%region%cell_ids(icell)
               ghosted_id = grid%nL2G(local_id)
               if (associated(cur_patch%imat)) then
                 if (cur_patch%imat(ghosted_id) <= 0) then
                   iphase_loc_p(ghosted_id) = 0
                   cycle
                 endif
               endif
               iphase_loc_p(ghosted_id)=initial_condition%flow_condition%iphase
             enddo
           else
             do icell=1,initial_condition%region%num_cells
               local_id = initial_condition%region%cell_ids(icell)
               ghosted_id = grid%nL2G(local_id)
               if (associated(cur_patch%imat)) then
                 if (cur_patch%imat(ghosted_id) <= 0) then
                   iphase_loc_p(ghosted_id) = 0
                   cycle
                 endif
               endif
               iphase_loc_p(ghosted_id)=initial_condition%flow_aux_int_var(1,icell)
             enddo
             do iface=1,initial_condition%numfaces_set
                ghosted_id = initial_condition%faces_set(iface)
                local_id = grid%fG2L(ghosted_id)
                if (local_id > 0) then
                  iend = local_id*option%nflowdof
                  ibegin = iend-option%nflowdof+1
                  xx_p(ibegin:iend) = &
                 initial_condition%flow_aux_real_var(1:option%nflowdof,iface)
                end if
             end do
           
           end if
#endif
        else 
           if (.not.associated(initial_condition%flow_aux_real_var)) then
          if (.not.associated(initial_condition%flow_condition)) then
            option%io_buffer = 'Flow condition is NULL in initial condition'
            call printErrMsg(option)
          endif
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
        end if
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

#ifdef DASVYAT
  if (discretization%itype == STRUCTURED_GRID_MIMETIC) then
   call DiscretizationGlobalToLocalFaces(discretization, field%flow_xx_faces, field%flow_xx_loc_faces, NFLOWDOF)

   call VecCopy(field%flow_xx_faces, field%flow_yy_faces, ierr)
   call MFDInitializeMassMatrices(realization%discretization%grid,&
                                  realization%field, &
                                  realization%discretization%MFD, realization%option)
  end if
#endif
!  stop 
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
  PetscInt :: irxn, isite
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
  PetscInt :: offset
  
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
                  global_aux_vars(ghosted_id),reaction, &
                  initial_condition%tran_condition%cur_constraint_coupler%constraint_name, &
                  initial_condition%tran_condition%cur_constraint_coupler%aqueous_species, &
                  initial_condition%tran_condition%cur_constraint_coupler%surface_complexes, &
                  initial_condition%tran_condition%cur_constraint_coupler%colloids, &
                  initial_condition%tran_condition%cur_constraint_coupler%num_iterations, &
                  PETSC_TRUE,option)
              else
                call RTAuxVarCopy(rt_aux_vars(ghosted_id), &
                  rt_aux_vars(grid%nL2G(initial_condition%region%cell_ids(icell-1))), &
                  option)
                call ReactionEquilibrateConstraint(rt_aux_vars(ghosted_id), &
                  global_aux_vars(ghosted_id),reaction, &
                  initial_condition%tran_condition%cur_constraint_coupler%constraint_name, &
                  initial_condition%tran_condition%cur_constraint_coupler%aqueous_species, &
                  initial_condition%tran_condition%cur_constraint_coupler%surface_complexes, &
                  initial_condition%tran_condition%cur_constraint_coupler%colloids, &
                  initial_condition%tran_condition%cur_constraint_coupler%num_iterations, &
                  PETSC_FALSE,option)
              endif
            endif
            offset = ibegin + reaction%offset_aq - 1
            do idof = 1, reaction%naqcomp ! primary aqueous concentrations
              xx_p(offset+idof) = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                aqueous_species%basis_molarity(idof) / &
                global_aux_vars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
            enddo
            ! colloids fractions
            if (associated(initial_condition%tran_condition%cur_constraint_coupler%colloids)) then
              offset = ibegin + reaction%offset_coll - 1
              do idof = 1, reaction%ncoll ! primary aqueous concentrations
                xx_p(offset+idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                  colloids%basis_conc_mob(idof) / &
                  global_aux_vars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
                rt_aux_vars(ghosted_id)%colloid%conc_imb(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                  colloids%basis_conc_imb(idof)
              enddo
            endif
            ! mineral volume fractions
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
            ! kinetic surface complexes
            if (associated(initial_condition%tran_condition%cur_constraint_coupler%surface_complexes)) then
              do idof = 1, reaction%nkinsrfcplx
                rt_aux_vars(ghosted_id)%kinsrfcplx_conc(idof) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                  surface_complexes%basis_conc(idof)
              enddo
              do irxn = 1, reaction%nkinsrfcplxrxn
                isite = reaction%kinsrfcplx_rxn_to_site(irxn)
                rt_aux_vars(ghosted_id)%kinsrfcplx_free_site_conc(isite) = &
                  initial_condition%tran_condition%cur_constraint_coupler% &
                  surface_complexes%basis_free_site_conc(isite)
              enddo
            endif
            ! this is for the multi-rate surface complexation model
            if (reaction%kinmr_nrate > 0) then
              ! copy over total sorbed concentration
              rt_aux_vars(ghosted_id)%kinmr_total_sorb = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                rt_auxvar%kinmr_total_sorb
              ! copy over free site concentration
              rt_aux_vars(ghosted_id)%eqsrfcplx_free_site_conc = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                rt_auxvar%eqsrfcplx_free_site_conc
              ! copy over surface complex concentrations
              rt_aux_vars(ghosted_id)%eqsrfcplx_conc = &
                initial_condition%tran_condition%cur_constraint_coupler% &
                rt_auxvar%eqsrfcplx_conc
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
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr)


end subroutine RealizAssignTransportInitCond

! ************************************************************************** !
!
! RealizationScaleSourceSink: Scales select source/sinks based on perms
! author: Glenn Hammond
! date: 09/03/08
!
! ************************************************************************** !
subroutine RealizationScaleSourceSink(realization)

  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Patch_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(realization_type) :: realization
  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: perm_loc_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  
  PetscInt :: ghosted_neighbors(0:27)
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_loc_ptr,ierr)

  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      ! BIG-TIME warning here.  I assume that all source/sink cells are within 
      ! a single patch - geh

      grid => cur_patch%grid

      cur_source_sink => cur_patch%source_sinks%first
      do
        if (.not.associated(cur_source_sink)) exit

        call VecZeroEntries(field%work,ierr)
        call GridVecGetArrayF90(grid,field%work,vec_ptr,ierr)

        cur_connection_set => cur_source_sink%connection_set
    
        do iconn = 1, cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)

          select case(option%iflowmode)
            case(RICHARDS_MODE)
               call GridGetGhostedNeighbors(grid,ghosted_id,STAR_STENCIL, &
                                            x_width,y_width,z_width, &
                                            ghosted_neighbors,option)
               ! ghosted neighbors is ordered first in x, then, y, then z
               icount = 0
               sum = 0.d0
               ! x-direction
               do while (icount < 2*x_width)
                 icount = icount + 1
                 neighbor_ghosted_id = ghosted_neighbors(icount)
                 sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                             grid%structured_grid%dy(neighbor_ghosted_id)* &
                             grid%structured_grid%dz(neighbor_ghosted_id)
                 
               enddo
               ! y-direction
               do while (icount < 2*(x_width+y_width))
                 icount = icount + 1
                 neighbor_ghosted_id = ghosted_neighbors(icount)                 
                 sum = sum + perm_loc_ptr(neighbor_ghosted_id)* &
                             grid%structured_grid%dx(neighbor_ghosted_id)* &
                             grid%structured_grid%dz(neighbor_ghosted_id)
                 
               enddo
               ! z-direction
               do while (icount < 2*(x_width+y_width+z_width))
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
            case(RICHARDS_MODE)
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
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_loc_ptr, ierr)
   
end subroutine RealizationScaleSourceSink

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
  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
                                   field%tortuosity_loc,ONEDOF)
                           
end subroutine RealizationRevertFlowParameters

! ************************************************************************** !
!
! RealizUpdateUniformVelocity: Assigns uniform velocity for transport
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RealizUpdateUniformVelocity(realization)

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
      call PatchUpdateUniformVelocity(cur_patch, &
                                      realization%velocity_dataset%cur_value, &
                                      realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
 
end subroutine RealizUpdateUniformVelocity

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

  if (associated(realization%velocity_dataset)) then
    if (realization%velocity_dataset%times(1) > 1.d-40 .or. &
        size(realization%velocity_dataset%times) > 1) then
      do itime = 1, size(realization%velocity_dataset%times)
        waypoint => WaypointCreate()
        waypoint%time = realization%velocity_dataset%times(itime)
        waypoint%update_srcs = PETSC_TRUE
        call WaypointInsertInList(waypoint,waypoint_list)
      enddo
    endif
  endif
  
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
subroutine RealizationGetDataset(realization,vec,ivar,isubvar,isubvar1)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchGetDataset(cur_patch,realization%field,realization%option, &
         realization%output_option,vec,ivar,isubvar,isubvar1)
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
function RealizGetDatasetValueAtCell(realization,ivar,isubvar,ghosted_id,isubvar1)

  use Option_module

  implicit none
  
  PetscReal :: RealizGetDatasetValueAtCell
  type(realization_type) :: realization
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
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
                realization%option,realization%output_option, &
                ivar,isubvar,ghosted_id,isubvar1)
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
! RealizationUpdateProperties: Updates coupled properties at each grid cell
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RealizationUpdateProperties(realization)

  type(realization_type) :: realization
  
  type(option_type), pointer :: option  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  PetscReal :: min_value  
  PetscInt :: ivalue
  PetscErrorCode :: ierr
  
  option => realization%option
    
  cur_level => realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      realization%patch => cur_patch
      call RealizationUpdatePropertiesPatch(realization)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  ! perform check to ensure that porosity is bounded between 0 and 1
  ! since it is calculated as 1.d-sum_volfrac, it cannot be > 1
  call VecMin(realization%field%porosity_loc,ivalue,min_value,ierr)
  if (min_value < 0.d0) then
    write(option%io_buffer,*) 'Sum of mineral volume fractions has ' // &
      'exceeded 1.d0 at cell (note PETSc numbering): ', ivalue
    call printErrMsg(option)
  endif
   
end subroutine RealizationUpdateProperties

! ************************************************************************** !
!
! RealizationUpdatePropertiesPatch: Updates coupled properties at each grid cell 
! author: Glenn Hammond
! date: 08/05/09
!
! ************************************************************************** !
subroutine RealizationUpdatePropertiesPatch(realization)

  use Grid_module
  use Reactive_Transport_Aux_module
 
  implicit none
  
  type(realization_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:) 
  type(discretization_type), pointer :: discretization

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl
  PetscReal :: sum_volfrac
  PetscReal :: scale, porosity_scale, volfrac_scale
  PetscBool :: porosity_updated
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: porosity_loc_p(:), porosity0_p(:)
  PetscReal, pointer :: tortuosity_loc_p(:), tortuosity0_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  reaction => realization%reaction
  grid => patch%grid
  material_property_array => realization%material_property_array
  rt_auxvars => patch%aux%RT%aux_vars

  if (.not.associated(patch%imat)) then
    option%io_buffer = 'Materials IDs not present in run.  Material ' // &
      ' properties cannot be updated without material ids ask Glenn'
    call printErrMsg(option)
  endif

  porosity_updated = PETSC_FALSE
  if (realization%option%update_porosity) then
    porosity_updated = PETSC_TRUE
  
    call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)

    if (reaction%nkinmnrl > 0) then
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)

        ! Go ahead and compute for inactive cells since their porosity does
        ! not matter (avoid check on active/inactive)
        sum_volfrac = 0.d0
        do imnrl = 1, reaction%nkinmnrl
          sum_volfrac = sum_volfrac + &
                        rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)
        enddo 
        porosity_loc_p(ghosted_id) = 1.d0-sum_volfrac
      enddo
    endif

    call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    
  endif
  
  if ((porosity_updated .and. &
       (option%update_tortuosity .or. &
        option%update_permeability)) .or. &
      ! if porosity ratio is used in mineral surface area update, we must
      ! recalculate it every time.
      (option%update_mineral_surface_area .and. &
       option%update_mnrl_surf_with_porosity)) then
    call GridVecGetArrayF90(grid,field%porosity0,porosity0_p,ierr)
    call GridVecGetArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = porosity_loc_p(ghosted_id)/porosity0_p(local_id)
    enddo
    call GridVecRestoreArrayF90(grid,field%porosity0,porosity0_p,ierr)
    call GridVecRestoreArrayF90(grid,field%porosity_loc,porosity_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
  endif      

  if (option%update_mineral_surface_area) then
    porosity_scale = 1.d0
    if (option%update_mnrl_surf_with_porosity) then
      ! placing the get/restore array calls within the condition will 
      ! avoid improper access.
      call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    endif
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (option%update_mnrl_surf_with_porosity) then
        porosity_scale = vec_p(local_id)** &
          material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_porosity_pwr
      endif
      do imnrl = 1, reaction%nkinmnrl
        if (rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) > 0.d0) then
          volfrac_scale = (rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)/ &
                         rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl))** &
            material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_volfrac_pwr
          rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
            rt_auxvars(ghosted_id)%mnrl_area0(imnrl)*porosity_scale*volfrac_scale
        else
          rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
            rt_auxvars(ghosted_id)%mnrl_area0(imnrl)
        endif
      enddo
    enddo
    if (option%update_mnrl_surf_with_porosity) then
      call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)
    endif

    call DiscretizationGlobalToLocal(discretization,field%tortuosity_loc, &
                                     field%tortuosity_loc,ONEDOF)
  endif
      
  if (option%update_tortuosity) then
    call GridVecGetArrayF90(grid,field%tortuosity_loc,tortuosity_loc_p,ierr)  
    call GridVecGetArrayF90(grid,field%tortuosity0,tortuosity0_p,ierr)  
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = vec_p(local_id)** &
        material_property_array(patch%imat(ghosted_id))%ptr%tortuosity_pwr
      tortuosity_loc_p(ghosted_id) = tortuosity0_p(local_id)*scale
    enddo
    call GridVecRestoreArrayF90(grid,field%tortuosity_loc,tortuosity_loc_p,ierr)  
    call GridVecRestoreArrayF90(grid,field%tortuosity0,tortuosity0_p,ierr)  
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)

    call DiscretizationGlobalToLocal(discretization,field%tortuosity_loc, &
                                     field%tortuosity_loc,ONEDOF)
  endif
      
  if (option%update_permeability) then
    call GridVecGetArrayF90(grid,field%perm0_xx,perm0_xx_p,ierr)
    call GridVecGetArrayF90(grid,field%perm0_zz,perm0_zz_p,ierr)
    call GridVecGetArrayF90(grid,field%perm0_yy,perm0_yy_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_zz_loc,perm_zz_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%perm_yy_loc,perm_yy_loc_p,ierr)
    call GridVecGetArrayF90(grid,field%work,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = vec_p(local_id)** &
              material_property_array(patch%imat(ghosted_id))%ptr%permeability_pwr
      perm_xx_loc_p(ghosted_id) = perm0_xx_p(local_id)*scale
      perm_yy_loc_p(ghosted_id) = perm0_yy_p(local_id)*scale
      perm_zz_loc_p(ghosted_id) = perm0_zz_p(local_id)*scale
    enddo
    call GridVecRestoreArrayF90(grid,field%perm0_xx,perm0_xx_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm0_zz,perm0_zz_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm0_yy,perm0_yy_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_xx_loc,perm_xx_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_zz_loc,perm_zz_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%perm_yy_loc,perm_yy_loc_p,ierr)
    call GridVecRestoreArrayF90(grid,field%work,vec_p,ierr)

    call DiscretizationGlobalToLocal(discretization,field%perm_xx_loc, &
                                     field%perm_xx_loc,ONEDOF)
    call DiscretizationGlobalToLocal(discretization,field%perm_yy_loc, &
                                     field%perm_yy_loc,ONEDOF)
    call DiscretizationGlobalToLocal(discretization,field%perm_zz_loc, &
                                     field%perm_zz_loc,ONEDOF)
  endif  
  
end subroutine RealizationUpdatePropertiesPatch
 
! ************************************************************************** !
!
! RealizationCountCells: Counts # of active and inactive grid cells 
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine RealizationCountCells(realization,global_total_count, &
                                 global_active_count,total_count,active_count)

  use Option_module

  implicit none
  
  type(realization_type) :: realization
  PetscInt :: global_total_count
  PetscInt :: global_active_count
  PetscInt :: total_count
  PetscInt :: active_count
  
  PetscInt :: patch_total_count
  PetscInt :: patch_active_count
  PetscInt :: temp_int_in(2), temp_int_out(2)
  PetscErrorCode :: ierr
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  total_count = 0
  active_count = 0
    
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchCountCells(cur_patch,patch_total_count,patch_active_count)
      total_count = total_count + patch_total_count
      active_count = active_count + patch_active_count
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
  temp_int_in(1) = total_count
  temp_int_in(2) = active_count
  call MPI_Allreduce(temp_int_in,temp_int_out,TWO_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,realization%option%mycomm,ierr)
  global_total_count = temp_int_out(1)
  global_active_count = temp_int_out(2)

end subroutine RealizationCountCells



subroutine RealizationSetUpBC4Faces(realization)




  use Connection_module
  use Coupler_module
  use Patch_module
  use Grid_module
  use Field_module
  use MFD_Aux_module
  

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(realization_type) :: realization

#ifdef DASVYAT

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  type(mfd_auxvar_type), pointer :: aux_var
  type(connection_set_type), pointer :: conn
  type(coupler_type), pointer ::  boundary_condition

  PetscReal, pointer :: bc_faces_p(:)
  PetscInt :: iconn, sum_connection, bc_type
  PetscInt :: local_id, ghosted_id, ghost_face_id, j, jface
  PetscErrorCode :: ierr


  patch => realization%patch
  grid => patch%grid
  field => realization%field



  call VecGetArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)

  boundary_condition => patch%boundary_conditions%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    bc_type = boundary_condition%flow_condition%itype(RICHARDS_PRESSURE_DOF)

    do iconn = 1, boundary_condition%numfaces_set
      sum_connection = sum_connection + 1

      local_id = boundary_condition%region%cell_ids(iconn)
      ghosted_id = grid%nL2G(local_id)

      aux_var => grid%MFD%aux_vars(local_id)
      do j = 1, aux_var%numfaces
        ghost_face_id = aux_var%face_id_gh(j)
        conn => grid%faces(ghost_face_id)%conn_set_ptr
        jface = grid%faces(ghost_face_id)%id
        if (boundary_condition%faces_set(iconn) == ghost_face_id) then
           if ((bc_type == DIRICHLET_BC).or.(bc_type == HYDROSTATIC_BC)  &
              .or.(bc_type == SEEPAGE_BC).or.(bc_type == CONDUCTANCE_BC) ) then
                    bc_faces_p(ghost_face_id) = boundary_condition%flow_aux_real_var(1,iconn)*conn%area(jface)
           else if ((bc_type == NEUMANN_BC)) then
                    bc_faces_p(ghost_face_id) = boundary_condition%flow_aux_real_var(1,iconn)*conn%area(jface)
           end if 
  !         write(*,*) ghost_face_id, boundary_condition%flow_aux_real_var(1,iconn), conn%cntr(3,jface)     
  !            bc_faces_p(ghost_face_id) = conn%cntr(3,jface)*conn%area(jface) 
        end if
      end do
    end do
    boundary_condition => boundary_condition%next
  end do


  call VecRestoreArrayF90(field%flow_bc_loc_faces, bc_faces_p, ierr)

!  write(*,*) "RealizationSetUpBC4Faces Finished"
!  read(*,*)
#endif

end subroutine RealizationSetUpBC4Faces

! ************************************************************************** !
!
! RealizationPrintGridStatistics: Prints statistics regarding the numerical
!                                 discretization 
! author: Glenn Hammond
! date: 06/01/10
!
! ************************************************************************** !
subroutine RealizationPrintGridStatistics(realization)

  use Grid_module

  implicit none

  type(realization_type) :: realization
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid

  PetscInt :: i1, i2, i3
  PetscReal :: r1, r2, r3
  PetscInt :: global_total_count, global_active_count
  PetscInt :: total_count, active_count
  PetscReal :: total_min, total_max, total_mean, total_variance
  PetscReal :: active_min, active_max, active_mean, active_variance
  PetscInt :: inactive_histogram(12), temp_int_out(12)
  PetscReal :: inactive_percentages(12)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%patch%grid

  ! print # of active and inactive grid cells
  call RealizationCountCells(realization,global_total_count, &
                             global_active_count,total_count,active_count)
  r1 = dble(total_count)
  call OptionMaxMinMeanVariance(r1,total_max, &
                                total_min,total_mean, &
                                total_variance,PETSC_TRUE,option)
  r1 = dble(active_count)
  call OptionMaxMinMeanVariance(r1,active_max, &
                                active_min,active_mean, &
                                active_variance,PETSC_TRUE,option)
                  
  r1 = dble(active_count) / dble(total_count)    
  inactive_histogram = 0                          
  if (r1 >= (1.d0-1.d-8)) then
    inactive_histogram(12) = 1
  else if (r1 >= .9d0 .and. r1 < (1.d0-1.d-8)) then
    inactive_histogram(11) = 1
  else if (r1 >= .8d0 .and. r1 < .9d0) then
    inactive_histogram(10) = 1
  else if (r1 >= .7d0 .and. r1 < .8d0) then
    inactive_histogram(9) = 1
  else if (r1 >= .6d0 .and. r1 < .7d0) then
    inactive_histogram(8) = 1
  else if (r1 >= .5d0 .and. r1 < .6d0) then
    inactive_histogram(7) = 1
  else if (r1 >= .4d0 .and. r1 < .5d0) then
    inactive_histogram(6) = 1
  else if (r1 >= .3d0 .and. r1 < .4d0) then
    inactive_histogram(5) = 1
  else if (r1 >= .2d0 .and. r1 < .3d0) then
    inactive_histogram(4) = 1
  else if (r1 >= .1d0 .and. r1 < .2d0) then
    inactive_histogram(3) = 1
  else if (r1 > 1.d-20 .and. r1 < .1d0) then
    inactive_histogram(2) = 1
  else if (r1 < 1.d-20) then
    inactive_histogram(1) = 1
  endif
  
  call MPI_Allreduce(inactive_histogram,temp_int_out,TWELVE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! why I cannot use *100, I do not know....geh
  inactive_percentages = dble(temp_int_out)/dble(option%mycommsize)*10.d0
  inactive_percentages = inactive_percentages+1.d-8

  r1 = 0.d0
  do i1 = 1, 12
    r1 = r1 + inactive_percentages(i1)
  enddo
                                
  i1 = -999
  i2 = -999
  i3 = -999
  if (associated(grid%structured_grid)) then
    i1 = grid%structured_grid%npx_final
    i2 = grid%structured_grid%npy_final
    i3 = grid%structured_grid%npz_final
  endif
  if (OptionPrintToScreen(option)) then
    write(*,'(/," Grid Stats:",/, &
              & "                       Global # cells: ",i12,/, &
              & "                Global # active cells: ",i12,/, &
              & "                              # cores: ",i12,/, &
              & "         Processor core decomposition: ",3i6,/, &
              & "               Maximum # cells / core: ",i12,/, &
              & "               Minimum # cells / core: ",i12,/, &
              & "               Average # cells / core: ",1pe12.4,/, &
              & "               Std Dev # cells / core: ",1pe12.4,/, &
              & "        Maximum # active cells / core: ",i12,/, &
              & "        Minimum # active cells / core: ",i12,/, &
              & "        Average # active cells / core: ",1pe12.4,/, &
              & "        Std Dev # active cells / core: ",1pe12.4,/,/, &
              & "        % cores with % active cells =       0%: ",1f7.2,/, &
              & "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
              & "        % cores with % active cells =   10-20%: ",1f7.2,/, &
              & "        % cores with % active cells =   20-30%: ",1f7.2,/, &
              & "        % cores with % active cells =   30-40%: ",1f7.2,/, &
              & "        % cores with % active cells =   40-50%: ",1f7.2,/, &
              & "        % cores with % active cells =   50-60%: ",1f7.2,/, &
              & "        % cores with % active cells =   60-70%: ",1f7.2,/, &
              & "        % cores with % active cells =   70-80%: ",1f7.2,/, &
              & "        % cores with % active cells =   80-90%: ",1f7.2,/, &
              & "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
              & "        % cores with % active cells =     100%: ",1f7.2,/, &
              & "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'(/," Grid Stats:",/, &
                "                       Global # cells: ",i12,/, &
                "                Global # active cells: ",i12,/, &
                "                              # cores: ",i12,/, &
                "         Processor core decomposition: ",3i6,/, &
                "               Maximum # cells / core: ",i12,/, &
                "               Minimum # cells / core: ",i12,/, &
                "               Average # cells / core: ",1pe12.4,/, &
                "               Std Dev # cells / core: ",1pe12.4,/, &
                "        Maximum # active cells / core: ",i12,/, &
                "        Minimum # active cells / core: ",i12,/, &
                "        Average # active cells / core: ",1pe12.4,/, &
                "        Std Dev # active cells / core: ",1pe12.4,/,/, &
                "        % cores with % active cells =       0%: ",1f7.2,/, &
                "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
                "        % cores with % active cells =   10-20%: ",1f7.2,/, &
                "        % cores with % active cells =   20-30%: ",1f7.2,/, &
                "        % cores with % active cells =   30-40%: ",1f7.2,/, &
                "        % cores with % active cells =   40-50%: ",1f7.2,/, &
                "        % cores with % active cells =   50-60%: ",1f7.2,/, &
                "        % cores with % active cells =   60-70%: ",1f7.2,/, &
                "        % cores with % active cells =   70-80%: ",1f7.2,/, &
                "        % cores with % active cells =   80-90%: ",1f7.2,/, &
                "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
                "        % cores with % active cells =     100%: ",1f7.2,/, &
                "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif

end subroutine RealizationPrintGridStatistics

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

  call DatasetDestroy(realization%datasets)
  
  call VelocityDatasetDestroy(realization%velocity_dataset)
  
  call DiscretizationDestroy(realization%discretization)
  
end subroutine RealizationDestroy
  
end module Realization_module
