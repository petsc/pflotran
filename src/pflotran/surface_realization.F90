#ifdef SURFACE_FLOW

module Surface_Realization_class

  use Realization_Base_class
  
  use Condition_module
  use Debug_module
  use Discretization_module
  use Input_module
  use Level_module
  use Option_module
  use Patch_module
  use Region_module
  use Surface_Field_module
  use Surface_Material_module
  use Waypoint_module
  use Dataset_Base_class
  use Reaction_Aux_module
  use Output_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
private


#include "finclude/petscsys.h"

  PetscReal, parameter :: eps       = 1.D-8

  type, public, extends(realization_base_type) :: surface_realization_type

    type(waypoint_list_type), pointer  :: waypoints
    
    type(surface_field_type), pointer                 :: surf_field
    type(region_list_type), pointer                   :: surf_regions
    type(condition_list_type),pointer                 :: surf_flow_conditions
    type(tran_condition_list_type),pointer            :: surf_transport_conditions
    type(surface_material_property_type), pointer     :: surf_material_properties
    type(surface_material_property_ptr_type), pointer :: surf_material_property_array(:)
    type(reaction_type), pointer                      :: surf_reaction
    character(len=MAXSTRINGLENGTH)                    :: surf_filename
    character(len=MAXSTRINGLENGTH)                    :: subsurf_filename

    class(dataset_base_type), pointer :: datasets
    
    PetscReal :: dt_max
    PetscReal :: dt_min
    PetscReal :: dt_coupling
    
    PetscInt :: iter_count
    PetscBool :: first_time

  end type surface_realization_type

  !         123456789+123456789+123456789+1
  public :: SurfRealizCreate, &
            SurfRealizDestroy, &
            SurfRealizAddWaypointsToList, &
            SurfRealizCreateDiscretization, &
            SurfRealizAddCoupler, &
            SurfRealizAddStrata, &
            SurfRealizLocalizeRegions, &
            SurfRealizPassFieldPtrToPatches, &
            SurfRealizProcessCouplers, &
            SurfRealizLocalToLocalWithArray, &
            SurfRealizProcessConditions, &
            SurfRealizProcessFlowConditions, &
            SurfRealizMapSurfSubsurfGrids, &
            SurfRealizInitAllCouplerAuxVars, &
            SurfRealizAllCouplerAuxVars, &
            SurfRealizProcessMatProp, &
            SurfRealizUpdate, &
!            SurfRealizCreateSurfSubsurfVec, &
!            SurfRealizUpdateSubsurfBC, &
!            SurfRealizUpdateSurfBC, &
!            SurfRealizSurf2SubsurfFlux, &
            SurfRealizGetVariable

  !TODO(intel)
!  public :: SurfaceRealizationGetVariable     
  
!  interface SurfaceRealizationGetVariable
!    module procedure :: RealizationGetVariable ! from Realization_Base_class
!  end interface
  
contains

! ************************************************************************** !
!> This routine allocates and initializes a new SurfaceRealization object
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/16/12
! ************************************************************************** !
function SurfRealizCreate(option)

  implicit none

  type(option_type), pointer             :: option
  type(surface_realization_type),pointer :: SurfRealizCreate
  type(surface_realization_type),pointer :: surf_realization
  
  allocate(surf_realization)
  call RealizationBaseInit(surf_realization,option)
  surf_realization%option => option
  nullify(surf_realization%input)

  surf_realization%surf_field => SurfaceFieldCreate()
  surf_realization%debug => DebugCreateFlow()
  surf_realization%output_option => OutputOptionCreate()
  surf_realization%level_list => LevelCreateList()

  nullify(surf_realization%surf_material_properties)

  allocate(surf_realization%surf_regions)
  call RegionInitList(surf_realization%surf_regions)
  
  allocate(surf_realization%surf_flow_conditions)
  call FlowConditionInitList(surf_realization%surf_flow_conditions)
  allocate(surf_realization%surf_transport_conditions)
  call TranConditionInitList(surf_realization%surf_transport_conditions)
  
  nullify(surf_realization%surf_reaction)
  nullify(surf_realization%datasets)

  surf_realization%iter_count = 0
  surf_realization%dt_min = 1.d0
  surf_realization%dt_max = 1.d0
  surf_realization%dt_coupling = 0.d0
  
  surf_realization%first_time = PETSC_TRUE
  SurfRealizCreate => surf_realization

end function SurfRealizCreate

! ************************************************************************** !
!> This routine adds a copy of a coupler to a list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/10/12
! ************************************************************************** !
subroutine SurfRealizAddCoupler(surf_relation,coupler)

  use Coupler_module

  implicit none
  
  type(surface_realization_type) :: surf_relation
  type(coupler_type), pointer    :: coupler
  
  type(level_type), pointer      :: cur_level
  type(patch_type), pointer      :: cur_patch
  type(coupler_type), pointer    :: new_coupler
  
  cur_level => surf_relation%level_list%first
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
 
end subroutine SurfRealizAddCoupler

! ************************************************************************** !
!> This routine sets connectivity and pointers for couplers related to 
!! surface flow.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfRealizProcessCouplers(surf_realization)

  use Option_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  type(level_type), pointer      :: cur_level
  type(patch_type), pointer      :: cur_patch
  
  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchProcessCouplers(cur_patch,surf_realization%surf_flow_conditions, &
                                surf_realization%surf_transport_conditions, &
                                surf_realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfRealizProcessCouplers

! ************************************************************************** !
!> This routine sets up linkeage between surface material properties
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfRealizProcessMatProp(surf_realization)

  use String_module
  
  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string

  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  
  option => surf_realization%option
  
  ! organize lists
  call SurfaceMaterialPropConvertListToArray( &
                                surf_realization%surf_material_properties, &
                                surf_realization%surf_material_property_array, &
                                option)
  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      cur_patch%surf_material_properties => surf_realization%surf_material_properties
      call SurfaceMaterialPropConvertListToArray( &
                                      cur_patch%surf_material_properties, &
                                      cur_patch%surf_material_property_array, &
                                      option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo 
  
end subroutine SurfRealizProcessMatProp

! ************************************************************************** !
!
!> This routine localizes surface regions within each patch
!! (similar to RealizationLocalizeRegions)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
!
! ************************************************************************** !
subroutine SurfRealizLocalizeRegions(surf_realization)

  use Option_module
  use String_module
  use Grid_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch
  type (region_type), pointer :: cur_region
  type(option_type), pointer :: option
  type(region_type), pointer :: patch_region

  option => surf_realization%option

  ! localize the regions on each patch
  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchLocalizeRegions(cur_patch,surf_realization%surf_regions, &
                                surf_realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfRealizLocalizeRegions


! ************************************************************************** !
!> This routine adds a copy of a strata to a list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfRealizAddStrata(surf_realization,strata)

  use Strata_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  type(strata_type), pointer     :: strata
  
  type(level_type), pointer      :: cur_level
  type(patch_type), pointer      :: cur_patch
  type(strata_type), pointer     :: new_strata
  
  cur_level => surf_realization%level_list%first
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
 
end subroutine SurfRealizAddStrata

! ************************************************************************** !
!> This routine creates grid
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfRealizCreateDiscretization(surf_realization)

  use Grid_module
  use Unstructured_Grid_Aux_module, only : UGridMapIndices
  use Unstructured_Grid_module, only     : UGridEnsureRightHandRule
  use Coupler_module
  use Discretization_module
  use Unstructured_Cell_module
  use DM_Kludge_module
  
  implicit none

  type(surface_realization_type)      :: surf_realization
  type(discretization_type), pointer  :: discretization
  type(grid_type), pointer            :: grid
  type(surface_field_type), pointer   :: surf_field
  type(option_type), pointer          :: option
  type(dm_ptr_type), pointer          :: dm_ptr

  PetscErrorCode :: ierr

  option => surf_realization%option
  surf_field => surf_realization%surf_field
  discretization => surf_realization%discretization

  call DiscretizationCreateDMs(discretization,option)

  ! n degree of freedom, global
  call DiscretizationCreateVector(discretization,NFLOWDOF,surf_field%flow_xx, &
                                  GLOBAL,option)
  call VecSet(surf_field%flow_xx,0.d0,ierr)

  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_yy)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_dxx)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_r)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_accum)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%work)

  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%exchange_subsurf_2_surf)

  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%mannings0, &
                                  GLOBAL,option)
  call VecSet(surf_field%mannings0,0.d0,ierr)

  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%area)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%Dq)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%perm_xx)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%perm_yy)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%perm_zz)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%por)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%icap_loc)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%ithrm_loc)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%subsurf_xx)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%subsurf_yy)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%subsurf_zz)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%surf2subsurf_dist)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%surf2subsurf_dist_gravity)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%press_subsurf)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%temp_subsurf)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%sat_ice)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%ckwet)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%ckdry)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%ckice)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%th_alpha)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%th_alpha_fr)

  ! n degrees of freedom, local
  call DiscretizationCreateVector(discretization,NFLOWDOF,surf_field%flow_xx_loc, &
                                  LOCAL,option)
  call VecSet(surf_field%flow_xx_loc,0.d0,ierr)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx_loc, &
                                     surf_field%work_loc)

  ! 1-dof degrees of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%mannings_loc, &
                                  LOCAL,option)
  call VecSet(surf_field%mannings_loc,0.d0,ierr)

  grid => discretization%grid

  ! set up nG2L, NL2G, etc.
  call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
                        grid%nG2L,grid%nL2G,grid%nG2A,grid%nG2P,option)
  call GridComputeCoordinates(grid,discretization%origin,option, & 
                              discretization%dm_1dof%ugdm) 
  call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                grid%y,grid%z,grid%nG2A,grid%nL2G,option)

  ! set up internal connectivity, distance, etc.
  call GridComputeInternalConnect(grid,option,discretization%dm_1dof%ugdm) 
  call GridComputeAreas(grid,surf_field%area,option)

  ! Allocate vectors to hold flowrate quantities
  if(surf_realization%output_option%print_hdf5_mass_flowrate.or. &
     surf_realization%output_option%print_hdf5_energy_flowrate.or. &
     surf_realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
     surf_realization%output_option%print_hdf5_aveg_energy_flowrate) then

    call VecCreateMPI(option%mycomm, &
         (option%nflowdof*MAX_FACE_PER_CELL_SURF+1)*surf_realization%patch%grid%nlmax, &
          PETSC_DETERMINE,surf_field%flowrate_inst,ierr)
    call VecSet(surf_field%flowrate_inst,0.d0,ierr)

    ! If average flowrate has to be saved, create a vector for it
    if(surf_realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
       surf_realization%output_option%print_hdf5_aveg_energy_flowrate) then
      call VecCreateMPI(option%mycomm, &
          (option%nflowdof*MAX_FACE_PER_CELL_SURF+1)*surf_realization%patch%grid%nlmax, &
          PETSC_DETERMINE,surf_field%flowrate_aveg,ierr)
    call VecSet(surf_field%flowrate_aveg,0.d0,ierr)
    endif
  endif

end subroutine SurfRealizCreateDiscretization

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/19/12
! ************************************************************************** !
subroutine SurfRealizPassFieldPtrToPatches(surf_realization)

  use Option_module

  implicit none
  
  type(surface_realization_type) :: surf_realization

  type(level_type), pointer      :: cur_level
  type(patch_type), pointer      :: cur_patch

  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      cur_patch%surf_field => surf_realization%surf_field
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  
end subroutine SurfRealizPassFieldPtrToPatches

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/19/12
! ************************************************************************** !
subroutine SurfRealizProcessConditions(surf_realization)

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  if (surf_realization%option%nflowdof > 0) then
    call SurfRealizProcessFlowConditions(surf_realization)
  endif
  if (surf_realization%option%ntrandof > 0) then
    !call SurfaceRealProcessTranConditions(surf_realization)
  endif
 
end subroutine SurfRealizProcessConditions

! ************************************************************************** !
!> This routine takes an F90 array that is ghosted and updates the ghosted
!! values (similar to RealLocalToLocalWithArray)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/13/12
! ************************************************************************** !
subroutine SurfRealizLocalToLocalWithArray(surf_realization,array_id)

  use Grid_module
  use Surface_Field_module

  implicit none

  type(surface_realization_type) :: surf_realization
  PetscInt                       :: array_id
  
  type(level_type), pointer         :: cur_level
  type(patch_type), pointer         :: cur_patch
  type(grid_type), pointer          :: grid
  type(surface_field_type), pointer :: surf_field

  surf_field => surf_realization%surf_field

  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      grid => cur_patch%grid
      select case(array_id)
        case(MATERIAL_ID_ARRAY)
          call GridCopyIntegerArrayToVec(grid, cur_patch%imat,surf_field%work_loc, &
                                         grid%ngmax)
      end select
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo
  call DiscretizationLocalToLocal(surf_realization%discretization, &
                                  surf_field%work_loc, &
                                  surf_field%work_loc,ONEDOF)
  cur_level => surf_realization%level_list%first
  do
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      grid => cur_patch%grid

      select case(array_id)
        case(MATERIAL_ID_ARRAY)
          call GridCopyVecToIntegerArray(grid, cur_patch%imat,surf_field%work_loc, &
                                         grid%ngmax)
      end select
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfRealizLocalToLocalWithArray

! ************************************************************************** !
!> This routine sets linkage of flow conditions to dataset
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/20/12
! ************************************************************************** !
subroutine SurfRealizProcessFlowConditions(surf_realization)

  use Dataset_Base_class
  use Dataset_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(flow_condition_type), pointer     :: cur_surf_flow_condition
  type(flow_sub_condition_type), pointer :: cur_surf_flow_sub_condition
  type(option_type), pointer             :: option
  character(len=MAXSTRINGLENGTH)         :: string
  character(len=MAXWORDLENGTH)           :: dataset_name
  class(dataset_base_type), pointer      :: dataset
  PetscInt                               :: i
  
  option => surf_realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_surf_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_surf_flow_condition)) exit
    !TODO(geh): could destroy the time_series here if dataset allocated
    select case(option%iflowmode)
      case(RICHARDS_MODE,TH_MODE)
        do i = 1, size(cur_surf_flow_condition%sub_condition_ptr)
          ! check for dataset in flow_dataset
          if (associated(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          flow_dataset%dataset)) then
            dataset_name = cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          flow_dataset%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                                flow_dataset%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_surf_flow_condition%name)
            dataset => &
              DatasetBaseGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%flow_dataset%dataset => &
              dataset
            !call DatasetLoad(dataset,surf_realization%discretization%dm_1dof, &
            !                 option)
          endif
          if (associated(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          datum%dataset)) then
            dataset_name = cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          datum%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                                datum%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_surf_flow_condition%name)
            dataset => &
              DatasetBaseGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%datum%dataset => &
              dataset
            !call DatasetXYZLoad(dataset,option)
          endif
          if (associated(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          gradient%dataset)) then
            dataset_name = cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                          gradient%dataset%name
            ! delete the dataset since it is solely a placeholder
            call DatasetDestroy(cur_surf_flow_condition%sub_condition_ptr(i)%ptr% &
                                gradient%dataset)
            ! get dataset from list
            string = 'flow_condition ' // trim(cur_surf_flow_condition%name)
            dataset => &
              DatasetBaseGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%gradient%dataset => &
              dataset
            !call DatasetLoad(dataset,surf_realization%discretization%dm_1dof, &
            !                 option)
          endif
        enddo
      case default
        option%io_buffer='SurfRealizProcessFlowConditions not implemented in this mode'
        call printErrMsg(option)
    end select
    cur_surf_flow_condition => cur_surf_flow_condition%next
  enddo

end subroutine SurfRealizProcessFlowConditions

! ************************************************************************** !
!> This routine initializez coupler auxillary variables within list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/21/12
! ************************************************************************** !
subroutine SurfRealizInitAllCouplerAuxVars(surf_realization)

  use Option_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

  call FlowConditionUpdate(surf_realization%surf_flow_conditions, &
                           surf_realization%option, &
                           surf_realization%option%time)

  cur_level => surf_realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      call PatchInitAllCouplerAuxVars(cur_patch, &
                                      surf_realization%option)
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

end subroutine SurfRealizInitAllCouplerAuxVars

! ************************************************************************** !
!> This routine updates auxiliary variables associated with couplers in the
!! list.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 04/18/13
! ************************************************************************** !
subroutine SurfRealizAllCouplerAuxVars(surf_realization,force_update_flag)

  use Option_module

  implicit none

  type(surface_realization_type) :: surf_realization
  PetscBool :: force_update_flag

  call PatchUpdateAllCouplerAuxVars(surf_realization%patch,force_update_flag, &
                                    surf_realization%option)

end subroutine SurfRealizAllCouplerAuxVars

! ************************************************************************** !
!> This routine creates vector scatter contexts between surface and subsurface 
!! grids.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! Algorithm:
!!  - It uses a similar logic of Matrix-Vector multiplication used in 
!!    UGridMapSideSet() subroutine. The algorithm here is extended to use 
!!    Matrix-Matrix mulitplication
!!
!! date: 01/17/12
! ************************************************************************** !
subroutine SurfRealizMapSurfSubsurfGrids(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Level_module
  use Patch_module
  use Region_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type), pointer         :: realization
  type(surface_realization_type), pointer :: surf_realization

  type(option_type), pointer           :: option
  type(unstructured_grid_type),pointer :: subsurf_grid
  type(unstructured_grid_type),pointer :: surf_grid
  type(level_type), pointer            :: cur_level
  type(patch_type), pointer            :: cur_patch 
  type(region_type), pointer           :: cur_region, top_region
  type(region_type), pointer           :: patch_region

  Mat :: Mat_vert_to_face_subsurf
  Mat :: Mat_vert_to_face_subsurf_transp
  Mat :: Mat_vert_to_face_surf
  Mat :: Mat_vert_to_face_surf_transp
  Mat :: prod
  Vec :: subsurf_petsc_ids, surf_petsc_ids

  PetscViewer :: viewer


  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4,1)
  PetscInt :: nvertices
  PetscInt :: iface
  PetscInt :: local_id, ii, jj
  PetscInt :: cell_type
  PetscInt :: ivertex, vertex_id_local
  PetscReal :: real_array4(4)
  PetscReal, pointer :: vec_ptr(:)

  PetscErrorCode :: ierr
  PetscBool :: found
  
  found = PETSC_FALSE
  
  if (.not.associated(realization)) return
  
  option => realization%option
  subsurf_grid => realization%discretization%grid%unstructured_grid
  surf_grid    => surf_realization%discretization%grid%unstructured_grid

  ! localize the regions on each patch
  cur_level => realization%level_list%first
  do 
    if (.not.associated(cur_level)) exit
    cur_patch => cur_level%patch_list%first
    do
      if (.not.associated(cur_patch)) exit
      cur_region => cur_patch%regions%first
        do
          if (.not.associated(cur_region)) exit
          if (StringCompare(cur_region%name,'top')) then
            found = PETSC_TRUE
            top_region => cur_region
            exit
          endif
          cur_region => cur_region%next
        enddo
      cur_patch => cur_patch%next
    enddo
    cur_level => cur_level%next
  enddo

  if(found.eqv.PETSC_FALSE) then
    option%io_buffer = 'When running with -DSURFACE_FLOW need to specify ' // &
      ' in the inputfile explicitly region: top '
    call printErrMsg(option)
  endif

  call MatCreateAIJ(option%mycomm, &
                       top_region%num_cells, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_subsurf, &
                       ierr)

   call MatCreateAIJ(option%mycomm, &
                     PETSC_DECIDE, &
                     top_region%num_cells, &
                     subsurf_grid%num_vertices_global, &
                     PETSC_DETERMINE, &
                     12, &
                     PETSC_NULL_INTEGER, &
                     12, &
                     PETSC_NULL_INTEGER, &
                     Mat_vert_to_face_subsurf_transp, &
                     ierr)

  call VecCreateMPI(option%mycomm,top_region%num_cells,PETSC_DETERMINE, &
                    subsurf_petsc_ids,ierr)  
  call MatZeroEntries(Mat_vert_to_face_subsurf,ierr)
  real_array4 = 1.d0

  call VecGetArrayF90(subsurf_petsc_ids,vec_ptr,ierr)

  offset=0
  call MPI_Exscan(top_region%num_cells,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii = 1, top_region%num_cells
    local_id = top_region%cell_ids(ii)
    vec_ptr(ii) = subsurf_grid%cell_ids_petsc(local_id)
    iface    = top_region%faces(ii)
    cell_type = subsurf_grid%cell_type(local_id)
    !nfaces = UCellGetNFaces(cell_type,option)

    call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                    int_array4)
    ! For this matrix:
    !   irow = local face id
    !   icol = natural (global) vertex id
    do ivertex = 1, nvertices
      vertex_id_local = &
        subsurf_grid%cell_vertices(int_array4(ivertex),local_id)
      int_array4_0(ivertex,1) = &
        subsurf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo
    call MatSetValues(Mat_vert_to_face_subsurf, &
                      1,ii-1+offset, &
                      nvertices,int_array4_0, &
                      real_array4, &
                      INSERT_VALUES,ierr)
    call MatSetValues(Mat_vert_to_face_subsurf_transp, &
                      nvertices,int_array4_0, &
                      1,ii-1+offset, &
                      real_array4, &
                      INSERT_VALUES,ierr)
  enddo

  call VecRestoreArrayF90(subsurf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyBegin(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_subsurf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_subsurf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_subsurf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_subsurf_transp,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  string = 'subsurf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(subsurf_petsc_ids,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call MatCreateAIJ(option%mycomm, &
                       surf_grid%nlmax, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_surf, &
                       ierr)

  call MatCreateAIJ(option%mycomm, &
                    PETSC_DECIDE, &
                    surf_grid%nlmax, &
                    subsurf_grid%num_vertices_global, &
                    PETSC_DETERMINE, &
                    12, &
                    PETSC_NULL_INTEGER, &
                    12, &
                    PETSC_NULL_INTEGER, &
                    Mat_vert_to_face_surf_transp, &
                    ierr)

  call VecCreateMPI(option%mycomm,surf_grid%nlmax,PETSC_DETERMINE, &
                    surf_petsc_ids,ierr)  
  offset=0
  call MPI_Exscan(surf_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  call VecGetArrayF90(surf_petsc_ids,vec_ptr,ierr)

  do local_id = 1, surf_grid%nlmax
    cell_type = surf_grid%cell_type(local_id)
    vec_ptr(local_id) = surf_grid%cell_ids_petsc(local_id)
    
    int_array4_0 = 0
    nvertices = surf_grid%cell_vertices(0,local_id)
    do ivertex = 1, nvertices
      vertex_id_local = surf_grid%cell_vertices(ivertex,local_id)
      int_array4_0(ivertex,1) = &
        surf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo    
   call MatSetValues(Mat_vert_to_face_surf,1,local_id-1+offset, &
                     nvertices,int_array4_0,real_array4, &
                     INSERT_VALUES,ierr)
   call MatSetValues(Mat_vert_to_face_surf_transp, &
                     nvertices,int_array4_0, &
                     1,local_id-1+offset, &
                     real_array4, &
                     INSERT_VALUES,ierr)
  enddo

  call VecRestoreArrayF90(surf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_surf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  string = 'surf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(surf_petsc_ids,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_surf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_surf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_surf_transp,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call MatMatMult(Mat_vert_to_face_subsurf,Mat_vert_to_face_surf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,prod,ierr)

#if UGRID_DEBUG
  string = 'prod.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(prod,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  call SurfRealizMapSurfSubsurfGrid(realization, surf_realization, prod, TWO_DIM_GRID, &
                                     surf_petsc_ids)
  call MatDestroy(prod,ierr)

  call MatMatMult(Mat_vert_to_face_surf,Mat_vert_to_face_subsurf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,prod,ierr)

#if UGRID_DEBUG
  string = 'prod_2.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(prod,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  call SurfRealizMapSurfSubsurfGrid(realization, surf_realization, prod, THREE_DIM_GRID, &
                                        subsurf_petsc_ids)

  call MatDestroy(prod,ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr)

  call VecDestroy(subsurf_petsc_ids,ierr)
  call VecDestroy(surf_petsc_ids,ierr)
  
end subroutine SurfRealizMapSurfSubsurfGrids

! ************************************************************************** !
!
!> This subroutine creates a single vector scatter context
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! Algorithm:
!!  - 
!!
!! date: 01/18/12
! ************************************************************************** !
subroutine SurfRealizMapSurfSubsurfGrid( &
              realization, &       !<
              surf_realization, &  !<
              prod_mat, &          !< Mat-Mat-Mult matrix
              source_grid_flag, &  !< To identify a surface or subsurface grid
              source_petsc_ids &   !< MPI-Vector containing cell ids in PETSc order
              )

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Cell_module
  use Realization_class
  use Option_module
  use Field_module
  use Surface_Field_module
  use Unstructured_Grid_module
  use Discretization_module
  use Unstructured_Grid_Aux_module
  use DM_Kludge_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type), pointer         :: realization
  type(surface_realization_type), pointer :: surf_realization
  Mat :: prod_mat
  PetscInt  :: source_grid_flag
  Vec :: source_petsc_ids

  Mat :: prod_loc_mat
  Vec :: source_loc_vec
  Vec :: corr_dest_ids_vec
  Vec :: corr_dest_ids_vec_ndof
  Vec :: source_petsc_ids_ndof
  IS  :: is_tmp1, is_tmp2
  IS  :: is_tmp3, is_tmp4
  PetscInt,pointer :: corr_v2_ids(:)
  VecScatter :: scatter
  VecScatter :: scatter_ndof

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(surface_field_type), pointer :: surf_field

  type(dm_ptr_type), pointer :: dm_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4,1)
  PetscReal :: real_array4(4)
  PetscInt :: ii, jj
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ivertex, cell_id, vertex_id_local
  PetscReal :: max_value

  PetscInt, pointer               :: ia_p(:), ja_p(:)
  PetscInt                        :: nrow,rstart,rend,icol(1)
  PetscInt                        :: index
  PetscInt                        :: vertex_id
  PetscOffset                     :: iia,jja,aaa,iicol
  PetscBool                       :: done
  PetscScalar                     :: aa(1)

  PetscErrorCode :: ierr
  PetscBool :: found
  PetscInt :: nlocal

  option     => realization%option
  field      => realization%field
  surf_field => surf_realization%surf_field
  
  if(option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(prod_mat,MAT_INITIAL_MATRIX,prod_loc_mat,ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_loc_mat, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(prod_loc_mat, aa, aaa, ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_mat, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArray(prod_mat, aa, aaa, ierr)
  endif

  ! For each row of the local-matrix, find the column with the largest value
  allocate(corr_v2_ids(nrow))
  do ii = 1, nrow
    max_value = 0.d0
    do jj = ia_p(ii), ia_p(ii + 1) - 1
      if(aa(aaa+ jj) > max_value) then
        corr_v2_ids(ii) = ja_p(jj)
        max_value = aa(aaa+ jj)
      endif
    enddo
    if(max_value<3) then
      option%io_buffer = 'Atleast three vertices need to form a face'
      call printErrMsg(option)
    endif
  enddo

  offset = 0
  call MPI_Exscan(nrow,offset,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(int_array(nrow))
  do ii = 1, nrow
    int_array(ii) = (ii-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp3,ierr)

  do ii = 1, nrow
    int_array(ii) = corr_v2_ids(ii)-1
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp4,ierr)
  deallocate(int_array)
  
  call VecCreateMPI(option%mycomm,nrow,PETSC_DETERMINE, &
                    corr_dest_ids_vec,ierr)
  call VecScatterCreate(source_petsc_ids,is_tmp2,corr_dest_ids_vec,is_tmp1, &
                        scatter,ierr)
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)
  
  call VecScatterBegin(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  select case(source_grid_flag)
    case(TWO_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids,ierr)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids_1dof,ierr)
    case(THREE_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids,ierr)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids_1dof,ierr)
  end select
  call VecScatterDestroy(scatter,ierr)

#if UGRID_DEBUG
  if(source_grid_flag==TWO_DIM_GRID) write(string,*) 'surf'
  if(source_grid_flag==THREE_DIM_GRID) write(string,*) 'subsurf'
  string = adjustl(string)
  string = 'corr_dest_ids_vec_' // trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(corr_dest_ids_vec,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecDestroy(corr_dest_ids_vec,ierr)
  if(option%mycommsize>1) call MatDestroy(prod_loc_mat,ierr)

#if UGRID_DEBUG
  if(source_grid_flag==TWO_DIM_GRID) write(string,*) 'surf'
  if(source_grid_flag==THREE_DIM_GRID) write(string,*) 'subsurf'
  string = adjustl(string)
  string = 'scatter_bet_grids_' // trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecScatterView(dm_ptr%ugdm%scatter_bet_grids,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! Create stridded vectors
  call VecCreate(option%mycomm,corr_dest_ids_vec_ndof,ierr)
  call VecSetSizes(corr_dest_ids_vec_ndof,nrow*option%nflowdof, &
                  PETSC_DECIDE,ierr)  
  call VecSetBlockSize(corr_dest_ids_vec_ndof,option%nflowdof,ierr)
  call VecSetFromOptions(corr_dest_ids_vec_ndof,ierr)

  call VecGetLocalSize(source_petsc_ids,nlocal,ierr)
  call VecCreate(option%mycomm,source_petsc_ids_ndof,ierr)
  call VecSetSizes(source_petsc_ids_ndof,nlocal*option%nflowdof, &
                  PETSC_DECIDE,ierr)  
  call VecSetBlockSize(source_petsc_ids_ndof,option%nflowdof,ierr)
  call VecSetFromOptions(source_petsc_ids_ndof,ierr)

  ! Create stridded vectors-scatter context
  call VecScatterCreate(source_petsc_ids_ndof,is_tmp4, &
                        corr_dest_ids_vec_ndof,is_tmp3, &
                        scatter_ndof,ierr)

  ! Save the stridded vectors-scatter context
  select case(source_grid_flag)
    case(TWO_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,NFLOWDOF)
    case(THREE_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,NFLOWDOF)
  end select
  call VecScatterCopy(scatter_ndof,dm_ptr%ugdm%scatter_bet_grids_ndof,ierr)

  ! Cleanup
  call VecScatterDestroy(scatter_ndof,ierr)
  call VecDestroy(source_petsc_ids_ndof,ierr)
  call VecDestroy(corr_dest_ids_vec_ndof,ierr)

end subroutine SurfRealizMapSurfSubsurfGrid

! ************************************************************************** !
!> This routine destroys SurfRealiz object
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/16/12
! ************************************************************************** !
subroutine SurfRealizDestroy(surf_realization)

  implicit none
  
  type(surface_realization_type), pointer :: surf_realization
  
  if(.not.associated(surf_realization)) return
  
  call SurfaceFieldDestroy(surf_realization%surf_field)
  
  call RegionDestroyList(surf_realization%surf_regions)
  
  call FlowConditionDestroyList(surf_realization%surf_flow_conditions)
  
  call LevelDestroyList(surf_realization%level_list)
  
  if(associated(surf_realization%debug)) deallocate(surf_realization%debug)
  nullify(surf_realization%debug)
  
  if (associated(surf_realization%surf_material_property_array)) &
    deallocate(surf_realization%surf_material_property_array)
  nullify(surf_realization%surf_material_property_array)
  call SurfaceMaterialPropertyDestroy(surf_realization%surf_material_properties)
  
  call DiscretizationDestroy(surf_realization%discretization)

end subroutine SurfRealizDestroy


! ************************************************************************** !
!> This routine updates parameters in realization (eg. conditions, bcs, srcs)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/22/12
! ************************************************************************** !
subroutine SurfRealizUpdate(surf_realization)

  implicit none
  
  type(surface_realization_type) :: surf_realization

  PetscBool :: force_update_flag = PETSC_FALSE

  ! must update conditions first
  call FlowConditionUpdate(surf_realization%surf_flow_conditions, &
                           surf_realization%option, &
                           surf_realization%option%time)

  call SurfRealizAllCouplerAuxVars(surf_realization,force_update_flag)

end subroutine SurfRealizUpdate

! ************************************************************************** !
!> This routine extracts variables indexed by ivar and isubvar from surface
!! realization.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 05/22/12
! ************************************************************************** !
subroutine SurfRealizGetVariable(surf_realization,vec,ivar,isubvar,isubvar1)

  use Option_module
  use Surface_Field_module

  implicit none

  type(surface_realization_type) :: surf_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call PatchGetVariable(surf_realization%patch, &
                       surf_realization%surf_field, &
                       !surf_realization%reaction, &
                       surf_realization%option, &
                       surf_realization%output_option, &
                       vec,ivar,isubvar,isubvar1)

end subroutine SurfRealizGetVariable

! ************************************************************************** !
!> This routine creates waypoints assocated with source/sink, boundary 
!! condition, etc. and adds to a list
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/15/13
! ************************************************************************** !
subroutine SurfRealizAddWaypointsToList(surf_realization)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  type(waypoint_list_type), pointer :: waypoint_list
  type(flow_condition_type), pointer :: cur_flow_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => surf_realization%option
  waypoint_list => surf_realization%waypoints
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_output = surf_realization%output_option%print_final
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in SurfRealizAddWaypointsToList'
    call printErrMsg(option)
  endif

  ! add update of flow conditions
  cur_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    if (cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        !TODO(geh): check if this updated more than simply the flow_dataset (i.e. datum and gradient)
        !geh: followup - no, datum/gradient are not considered.  Should they be considered?
        call FlowConditionDatasetGetTimes(option,sub_condition,final_time, &
                                          times)
        if (size(times) > 1000) then
          option%io_buffer = 'For flow condition "' // &
            trim(cur_flow_condition%name) // &
            '" dataset "' // trim(sub_condition%name) // &
            '", the number of times is excessive for synchronization ' // &
            'with waypoints.'
          call printErrMsg(option)
        endif
        do itime = 1, size(times)
          waypoint => WaypointCreate()
          waypoint%time = times(itime)
          waypoint%update_conditions = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        enddo
        deallocate(times)
        nullify(times)
      enddo
    endif
    cur_flow_condition => cur_flow_condition%next
  enddo
      
  ! add waypoints for periodic output
  if (surf_realization%output_option%periodic_output_time_incr > 0.d0 .or. &
      surf_realization%output_option%periodic_tr_output_time_incr > 0.d0) then

    if (surf_realization%output_option%periodic_output_time_incr > 0.d0) then
      ! standard output
      temp_real = 0.d0
      do
        temp_real = temp_real + surf_realization%output_option%periodic_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,surf_realization%waypoints)
      enddo
    endif
    
    if (surf_realization%output_option%periodic_tr_output_time_incr > 0.d0) then
      ! transient observation output
      temp_real = 0.d0
      do
        temp_real = temp_real + surf_realization%output_option%periodic_tr_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_tr_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,surf_realization%waypoints)
      enddo
    endif

  endif

end subroutine SurfRealizAddWaypointsToList

end module Surface_Realization_class

#endif
