#ifdef SURFACE_FLOW

module Surface_Realization_module

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
  use Dataset_Aux_module
  use Reaction_Aux_module
  
  implicit none
  
private

#include "definitions.h"

  type, public :: surface_realization_type

    PetscInt :: id
    
    type(discretization_type), pointer :: discretization
    type(level_list_type),pointer      :: level_list
    type(patch_type), pointer          :: patch

    type(option_type), pointer         :: option
    type(input_type), pointer          :: input
    type(flow_debug_type), pointer     :: debug
    type(output_option_type), pointer  :: output_option
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

    type(dataset_type), pointer                       :: datasets

  end type surface_realization_type

  public :: SurfaceRealizationCreate, &
            SurfaceRealizationDestroy, &
            SurfaceRealizationCreateDiscretization, &
            SurfaceRealizationAddCoupler, &
            SurfaceRealizationAddStrata, &
            SurfaceRealizationLocalizeRegions, &
            SurfaceRealizatonPassFieldPtrToPatches, &
            SurfaceRealizationProcessCouplers, &
            SurfaceRealizationLocalToLocalWithArray, &
            SurfaceRealizationProcessConditions, &
            SurfaceRealizationProcessFlowConditions, &
            SurfaceRealizationMapSurfSubsurfaceGrids, &
            SurfaceRealizationInitAllCouplerAuxVars, &
            SurfaceRealizationProcessMatProp

contains

! ************************************************************************** !
!> This routine allocates and initializes a new SurfaceRealization object
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/16/12
! ************************************************************************** !
function SurfaceRealizationCreate(option)

  implicit none

  type(option_type), pointer             :: option
  type(surface_realization_type),pointer :: SurfaceRealizationCreate
  type(surface_realization_type),pointer :: surf_realization
  
  allocate(surf_realization)
  surf_realization%discretization => DiscretizationCreate()
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

  SurfaceRealizationCreate => surf_realization

end function SurfaceRealizationCreate

! ************************************************************************** !
!> This routine adds a copy of a coupler to a list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/10/12
! ************************************************************************** !
subroutine SurfaceRealizationAddCoupler(surf_relation,coupler)

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
 
end subroutine SurfaceRealizationAddCoupler

! ************************************************************************** !
!> This routine sets connectivity and pointers for couplers related to 
!! surface flow.
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfaceRealizationProcessCouplers(surf_realization)

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

end subroutine SurfaceRealizationProcessCouplers

! ************************************************************************** !
!> This routine sets up linkeage between surface material properties
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfaceRealizationProcessMatProp(surf_realization)

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
  
end subroutine SurfaceRealizationProcessMatProp

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
subroutine SurfaceRealizationLocalizeRegions(surf_realization)

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

end subroutine SurfaceRealizationLocalizeRegions


! ************************************************************************** !
!> This routine adds a copy of a strata to a list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfaceRealizationAddStrata(surf_realization,strata)

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
 
end subroutine SurfaceRealizationAddStrata

! ************************************************************************** !
!> This routine creates grid
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/17/12
! ************************************************************************** !
subroutine SurfaceRealizationCreateDiscretization(surf_realization)

  use Grid_module
  use Unstructured_Grid_Aux_module, only : UGridMapIndices
  use Unstructured_Grid_module, only     : UGridEnsureRightHandRule
  use Coupler_module
  use Discretization_module
  
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

  ! Create Global vectors
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%flow_xx, &
                                  GLOBAL,option)
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
                                     surf_field%mannings0)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%area)

  ! Create ghosted-local vectors
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%flow_xx_loc, &
                                  LOCAL,option)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx_loc, &
                                     surf_field%work_loc)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx_loc, &
                                     surf_field%mannings_loc)

  grid => discretization%grid

  ! set up nG2L, NL2G, etc.
  call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
                        grid%nG2L,grid%nL2G,grid%nG2A)
  call GridComputeCoordinates(grid,discretization%origin,option, & 
                              discretization%dm_1dof%ugdm) 
  call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                grid%y,grid%z,grid%nG2A,grid%nL2G,option)

  ! set up internal connectivity, distance, etc.
  call GridComputeInternalConnect(grid,option,discretization%dm_1dof%ugdm) 
  call GridComputeAreas(grid,surf_field%area,option)

end subroutine SurfaceRealizationCreateDiscretization

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/19/12
! ************************************************************************** !
subroutine SurfaceRealizatonPassFieldPtrToPatches(surf_realization)

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
  
end subroutine SurfaceRealizatonPassFieldPtrToPatches

! ************************************************************************** !
!> This routine 
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/19/12
! ************************************************************************** !
subroutine SurfaceRealizationProcessConditions(surf_realization)

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  if (surf_realization%option%nflowdof > 0) then
    call SurfaceRealizationProcessFlowConditions(surf_realization)
  endif
  if (surf_realization%option%ntrandof > 0) then
    write(*,*), '>>>>>>>> Add SurfaceRealizationProcessTranConditions>>>>>>'
    !call SurfaceRealProcessTranConditions(surf_realization)
  endif
 
end subroutine SurfaceRealizationProcessConditions

! ************************************************************************** !
!> This routine takes an F90 array that is ghosted and updates the ghosted
!! values (similar to RealLocalToLocalWithArray)
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/13/12
! ************************************************************************** !
subroutine SurfaceRealizationLocalToLocalWithArray(surf_realization,array_id)

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

end subroutine SurfaceRealizationLocalToLocalWithArray

! ************************************************************************** !
!> This routine sets linkage of flow conditions to dataset
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/20/12
! ************************************************************************** !
subroutine SurfaceRealizationProcessFlowConditions(surf_realization)

  use Dataset_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(flow_condition_type), pointer     :: cur_surf_flow_condition
  type(flow_sub_condition_type), pointer :: cur_surf_flow_sub_condition
  type(option_type), pointer             :: option
  character(len=MAXSTRINGLENGTH)         :: string
  character(len=MAXWORDLENGTH)           :: dataset_name
  type(dataset_type), pointer            :: dataset
  PetscInt                               :: i
  
  option => surf_realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_surf_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_surf_flow_condition)) exit
    !TODO(geh): could destroy the time_series here if dataset allocated
    select case(option%iflowmode)
      case(G_MODE)
      case(RICHARDS_MODE)
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
              DatasetGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%flow_dataset%dataset => &
              dataset
            call DatasetLoad(dataset,option)
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
              DatasetGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%datum%dataset => &
              dataset
            call DatasetLoad(dataset,option)
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
              DatasetGetPointer(surf_realization%datasets,dataset_name,string,option)
            cur_surf_flow_condition%sub_condition_ptr(i)%ptr%gradient%dataset => &
              dataset
            call DatasetLoad(dataset,option)
          endif
        enddo
    end select
    cur_surf_flow_condition => cur_surf_flow_condition%next
  enddo

end subroutine SurfaceRealizationProcessFlowConditions

! ************************************************************************** !
!> This routine initializez coupler auxillary variables within list
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/21/12
! ************************************************************************** !
subroutine SurfaceRealizationInitAllCouplerAuxVars(surf_realization)

  use Option_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  type(level_type), pointer :: cur_level
  type(patch_type), pointer :: cur_patch

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
   
end subroutine SurfaceRealizationInitAllCouplerAuxVars

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
subroutine SurfaceRealizationMapSurfSubsurfaceGrids(realization,surf_realization)

  use Grid_module
  use String_module
  use Unstructured_Grid_module
  use Unstructured_Grid_Aux_module
  use Unstructured_Cell_module
  use Realization_module
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
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_subsurf, &
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
    call MatSetValues(Mat_vert_to_face_subsurf,1,ii-1+offset, &
                      nvertices,int_array4_0,real_array4, &
                      INSERT_VALUES,ierr)
  enddo

  call VecRestoreArrayF90(subsurf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face_subsurf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  call MatTranspose(Mat_vert_to_face_subsurf,MAT_INITIAL_MATRIX, &
                    Mat_vert_to_face_subsurf_transp,ierr)

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
                       PETSC_DETERMINE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_surf, &
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
  enddo

  call VecRestoreArrayF90(surf_petsc_ids,vec_ptr,ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY,ierr)

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

  call MatTranspose(Mat_vert_to_face_surf,MAT_INITIAL_MATRIX, &
                    Mat_vert_to_face_surf_transp,ierr)

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

  call SurfaceRealizationMapSurfSubsurfaceGrid(realization, surf_realization, prod, TWO_DIM_GRID, &
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
  call SurfaceRealizationMapSurfSubsurfaceGrid(realization, surf_realization, prod, THREE_DIM_GRID, &
                                        subsurf_petsc_ids)

  call MatDestroy(prod,ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr)

  call VecDestroy(subsurf_petsc_ids,ierr)
  call VecDestroy(surf_petsc_ids,ierr)

end subroutine SurfaceRealizationMapSurfSubsurfaceGrids

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
subroutine SurfaceRealizationMapSurfSubsurfaceGrid( &
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
  use Realization_module
  use Option_module
  use Field_module
  use Surface_Field_module
  use Unstructured_Grid_module
  use Discretization_module
  use Unstructured_Grid_Aux_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(realization_type), pointer         :: realization
  type(surface_realization_type), pointer :: surf_realization


  Mat :: prod_mat, prod_loc_mat
  Vec :: source_petsc_ids
  Vec :: source_loc_vec
  Vec :: corr_dest_ids_vec
  IS  :: is_tmp1, is_tmp2
  PetscInt,pointer :: corr_v2_ids(:)
  VecScatter :: scatter

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
  PetscInt  :: source_grid_flag

  PetscErrorCode :: ierr
  PetscBool :: found

  option     => realization%option
  field      => realization%field
  surf_field => surf_realization%surf_field
  
  if(option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(prod_mat,MAT_INITIAL_MATRIX,prod_loc_mat,ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_loc_mat, 1, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr)
    ! Get values stored in the local-matrix
    call MatGetArray(prod_loc_mat, aa, aaa, ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_mat, 1, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr)
    ! Get values stored in the local-matrix
    call MatGetArray(prod_mat, aa, aaa, ierr)
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

  do ii = 1, nrow
    int_array(ii) = corr_v2_ids(ii)-1
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)
  deallocate(int_array)
  
  call VecCreateMPI(option%mycomm,nrow,PETSC_DETERMINE, &
                    corr_dest_ids_vec,ierr)
  call VecScatterCreate(source_petsc_ids,is_tmp1,corr_dest_ids_vec,is_tmp2, &
                        scatter,ierr)
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)
  
  call VecScatterBegin(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
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

  call VecCreateSeq(PETSC_COMM_SELF,nrow,source_loc_vec,ierr)
  
  allocate(int_array(nrow))
  call VecGetArrayF90(corr_dest_ids_vec,vec_ptr,ierr)
  do ii = 1,nrow
    int_array(ii) = vec_ptr(ii)-1
  enddo
  call VecRestoreArrayF90(corr_dest_ids_vec,vec_ptr,ierr)
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr)
  call VecDestroy(corr_dest_ids_vec,ierr)
  
  do ii = 1,nrow
    int_array(ii) = ii-1
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr)
  
  select case(source_grid_flag)
    case(TWO_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)
      call VecScatterCreate(surf_field%flow_r,is_tmp1,source_loc_vec,is_tmp2, &
                            dm_ptr%ugdm%scatter_bet_grids,ierr)
    case(THREE_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
      call VecScatterCreate(field%flow_r,is_tmp1,source_loc_vec,is_tmp2, &
                            dm_ptr%ugdm%scatter_bet_grids,ierr)    
  end select
  
  call ISDestroy(is_tmp1,ierr)
  call ISDestroy(is_tmp2,ierr)
  call VecDestroy(source_loc_vec,ierr)
  
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

end subroutine SurfaceRealizationMapSurfSubsurfaceGrid

! ************************************************************************** !
!> This routine destroys SurfaceRealization object
!!
!> @author
!! Gautam Bisht, ORNL
!!
!! date: 02/16/12
! ************************************************************************** !
subroutine SurfaceRealizationDestroy(surf_realization)

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

end subroutine SurfaceRealizationDestroy


end module Surface_Realization_module
#endif