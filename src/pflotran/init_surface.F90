module Init_Surface_module

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"

  public :: InitSurfaceSetupRealization, &
            InitSurfaceSetupSolvers
contains

! ************************************************************************** !

subroutine InitSurfaceSetupRealization(simulation)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Simulation_module
  
  use Surface_Flow_module
  use Surface_Realization_class
  use Surface_TH_module
  use Surface_Global_module
  
  use Option_module
  use Waypoint_module
  use Condition_Control_module
  
  implicit none
  
  type(simulation_type) :: simulation
  
  type(option_type), pointer :: option
  
  option => simulation%realization%option
  
  if (option%surf_flow_on) then
    ! Check if surface-flow is compatible with the given flowmode
    select case(option%iflowmode)
      case(RICHARDS_MODE,TH_MODE)
      case default
        option%io_buffer = 'For surface-flow only RICHARDS and TH mode implemented'
        call printErrMsgByRank(option)
    end select

    call SurfaceInitReadRegionFiles(simulation%surf_realization)
    call SurfRealizMapSurfSubsurfGrids(simulation%realization, &
                                       simulation%surf_realization)
    call SurfRealizLocalizeRegions(simulation%surf_realization)
    call SurfRealizPassFieldPtrToPatches(simulation%surf_realization)
    call SurfRealizProcessMatProp(simulation%surf_realization)
    call SurfRealizProcessCouplers(simulation%surf_realization)
    call SurfRealizProcessConditions(simulation%surf_realization)
    !call RealProcessFluidProperties(simulation%surf_realization)
    call SurfaceInitMatPropToRegions(simulation%surf_realization)
    call SurfRealizInitAllCouplerAuxVars(simulation%surf_realization)
    !call SurfaceRealizationPrintCouplers(simulation%surf_realization)

    ! add waypoints associated with boundary conditions, source/sinks etc. to list
    call SurfRealizAddWaypointsToList(simulation%surf_realization)
    call WaypointListFillIn(option,simulation%surf_realization%waypoint_list)
    call WaypointListRemoveExtraWaypnts(option,simulation%surf_realization%waypoint_list)
    if (associated(simulation%flow_timestepper)) then
      simulation%surf_flow_timestepper%cur_waypoint => simulation%surf_realization%waypoint_list%first
    endif

    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call SurfaceFlowSetup(simulation%surf_realization)
      case default
      case(TH_MODE)
        call SurfaceTHSetup(simulation%surf_realization)
    end select

    call SurfaceGlobalSetup(simulation%surf_realization)
    ! initialize FLOW
    ! set up auxillary variable arrays

    ! assign initial conditionsRealizAssignFlowInitCond
    call CondControlAssignFlowInitCondSurface(simulation%surf_realization)

    ! override initial conditions if they are to be read from a file
    if (len_trim(option%surf_initialize_flow_filename) > 1) then
      option%io_buffer = 'For surface-flow initial conditions cannot be read from file'
      call printErrMsgByRank(option)
    endif
  
    select case(option%iflowmode)
      case(RICHARDS_MODE)
        call SurfaceFlowUpdateAuxVars(simulation%surf_realization)
      case(TH_MODE)
        call SurfaceTHUpdateAuxVars(simulation%surf_realization)
      case default
        option%io_buffer = 'For surface-flow only RICHARDS and TH mode implemented'
        call printErrMsgByRank(option)
    end select
  endif ! option%surf_flow_on
  
end subroutine InitSurfaceSetupRealization

! ************************************************************************** !

subroutine InitSurfaceSetupSolvers(surf_realization,solver)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Surface_Realization_class
  use Option_module
  
  use Solver_module
  use Convergence_module
  use Discretization_module
  use Surface_Flow_module
  use Surface_TH_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"
  
  type(surface_realization_type) :: surf_realization
  type(solver_type), pointer :: solver
  
  type(option_type), pointer :: option
  type(convergence_context_type), pointer :: convergence_context
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  option => surf_realization%option
  
  call printMsg(option,"  Beginning setup of FLOW SNES ")

  ! Setup PETSc TS for explicit surface flow solution
  call printMsg(option,"  Beginning setup of SURF FLOW TS ")

  call SolverCreateTS(solver,option%mycomm)
  call TSSetProblemType(solver%ts,TS_NONLINEAR, &
                        ierr);CHKERRQ(ierr)
  select case(option%iflowmode)
    case (RICHARDS_MODE)
      call TSSetRHSFunction(solver%ts,PETSC_NULL_OBJECT, &
                            SurfaceFlowRHSFunction, &
                            surf_realization, &
                            ierr);CHKERRQ(ierr)
    case (TH_MODE)
      call TSSetRHSFunction(solver%ts,PETSC_NULL_OBJECT, &
                            SurfaceTHRHSFunction, &
                            surf_realization, &
                            ierr);CHKERRQ(ierr)
      call TSSetIFunction(solver%ts,PETSC_NULL_OBJECT, &
                          SurfaceTHIFunction, &
                          surf_realization, &
                          ierr);CHKERRQ(ierr)
  end select
  call TSSetDuration(solver%ts,ONE_INTEGER, &
                      surf_realization%waypoint_list%last%time, &
                      ierr);CHKERRQ(ierr)
  
end subroutine InitSurfaceSetupSolvers

! ************************************************************************** !

subroutine SurfaceInitMatPropToRegions(surf_realization)
  ! 
  ! This routine assigns surface material properties to associated regions in
  ! the model (similar to assignMaterialPropToRegions)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/13/12
  ! 

  use Surface_Realization_class
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Surface_Field_module
  use Surface_Material_module
  
  use HDF5_module

  implicit none
  
  type(surface_realization_type) :: surf_realization
  
  PetscReal, pointer :: man0_p(:)
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, surf_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: patch  
  type(patch_type), pointer :: cur_patch

  type(surface_material_property_type), pointer :: surf_material_property
  type(surface_material_property_type), pointer :: null_surf_material_property
  type(region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  
  option => surf_realization%option
  discretization => surf_realization%discretization
  surf_field => surf_realization%surf_field

  ! loop over all patches and allocation material id arrays
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    if (.not.associated(cur_patch%imat)) then
      allocate(cur_patch%imat(cur_patch%grid%ngmax))
      ! initialize to "unset"
      cur_patch%imat = UNINITIALIZED_INTEGER
      ! also allocate saturation function id
      allocate(cur_patch%sat_func_id(cur_patch%grid%ngmax))
      cur_patch%sat_func_id = UNINITIALIZED_INTEGER
    endif
    cur_patch => cur_patch%next
  enddo

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    strata => cur_patch%strata_list%first
    do
      if (.not.associated(strata)) exit
      ! Read in cell by cell material ids if they exist
      if (.not.associated(strata%region) .and. strata%active) then
        option%io_buffer = 'Reading of material prop from file for' // &
          ' surface flow is not implemented.'
        call printErrMsgByRank(option)
        !call readMaterialsFromFile(realization,strata%realization_dependent, &
        !                           strata%material_property_filename)
      ! Otherwise, set based on region
      else if (strata%active) then
        update_ghosted_material_ids = PETSC_TRUE
        region => strata%region
        surf_material_property => strata%surf_material_property
        if (associated(region)) then
          istart = 1
          iend = region%num_cells
        else
          istart = 1
          iend = grid%nlmax
        endif
        do icell=istart, iend
          if (associated(region)) then
            local_id = region%cell_ids(icell)
          else
            local_id = icell
          endif
          ghosted_id = grid%nL2G(local_id)
          cur_patch%imat(ghosted_id) = surf_material_property%internal_id
        enddo
      endif
      strata => strata%next
    enddo
    cur_patch => cur_patch%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call SurfRealizLocalToLocalWithArray(surf_realization,MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_surf_material_property => SurfaceMaterialPropertyCreate()
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    call VecGetArrayF90(surf_field%mannings0,man0_p,ierr);CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      surf_material_id = cur_patch%imat(ghosted_id)
      if (surf_material_id == 0) then ! accomodate inactive cells
        surf_material_property = null_surf_material_property
      else if ( surf_material_id > 0 .and. &
                surf_material_id <= &
                size(surf_realization%surf_material_property_array)) then
        surf_material_property => &
          surf_realization%surf_material_property_array(surf_material_id)%ptr
        if (.not.associated(surf_material_property)) then
          write(dataset_name,*) surf_material_id
          option%io_buffer = 'No material property for surface material id ' // &
                              trim(adjustl(dataset_name)) &
                              //  ' defined in input file.'
          call printErrMsgByRank(option)
        endif
      else if (Uninitialized(surf_material_id)) then 
        write(dataset_name,*) grid%nG2A(ghosted_id)
        option%io_buffer = 'Uninitialized surface material id in patch at cell ' // &
                            trim(adjustl(dataset_name))
        call printErrMsgByRank(option)
      else if (surf_material_id > size(surf_realization%surf_material_property_array)) then
        write(option%io_buffer,*) surf_material_id
        option%io_buffer = 'Unmatched surface material id in patch:' // &
          adjustl(trim(option%io_buffer))
        call printErrMsgByRank(option)
      else
        option%io_buffer = 'Something messed up with surface material ids. ' // &
          ' Possibly material ids not assigned to all grid cells. ' // &
          ' Contact Glenn!'
        call printErrMsgByRank(option)
      endif
      man0_p(local_id) = surf_material_property%mannings
    enddo ! local_id - loop

    call VecRestoreArrayF90(surf_field%mannings0,man0_p,ierr);CHKERRQ(ierr)
      
    cur_patch => cur_patch%next
  enddo ! looping over patches
  
  call SurfaceMaterialPropertyDestroy(null_surf_material_property)
  nullify(null_surf_material_property)

  call DiscretizationGlobalToLocal(discretization,surf_field%mannings0, &
                                   surf_field%mannings_loc,ONEDOF)

end subroutine SurfaceInitMatPropToRegions

! ************************************************************************** !

subroutine SurfaceInitReadRegionFiles(surf_realization)
  ! 
  ! This routine reads surface region files
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/20/12
  ! 

  use Surface_Realization_class
  use Region_module
  use HDF5_module
  use Grid_module

  implicit none

  type(surface_realization_type) :: surf_realization
  
  type(region_type), pointer :: surf_region
  
  surf_region => surf_realization%surf_regions%first
  do 
    if (.not.associated(surf_region)) exit
    if (len_trim(surf_region%filename) > 1) then
      if (index(surf_region%filename,'.h5') > 0) then
        if (surf_region%grid_type == STRUCTURED_GRID) then
          !call HDF5ReadRegionFromFile(surf_realization,surf_region,surf_region%filename)
        else
#if defined(PETSC_HAVE_HDF5)
          call HDF5ReadUnstructuredGridRegionFromFile(surf_realization%option, &
                                                      surf_region, &
                                                      surf_region%filename)
#endif      
        endif
      else if (index(surf_region%filename,'.ss') > 0) then
        surf_region%sideset => RegionCreateSideset()
        call RegionReadFromFile(surf_region%sideset,surf_region%filename, &
                                surf_realization%option)
      else
        call RegionReadFromFile(surf_region,surf_realization%option, &
                                surf_region%filename)
      endif
    endif
    surf_region => surf_region%next
  enddo

end subroutine SurfaceInitReadRegionFiles

end module Init_Surface_module
