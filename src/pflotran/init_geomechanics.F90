module Init_Geomechanics_module

  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "finclude/petscsys.h"

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
#include "finclude/petscts.h"


  public :: InitGeomechSetupRealization, &
            InitGeomechSetupSolvers

contains

! ************************************************************************** !

subroutine InitGeomechSetupRealization(simulation)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Simulation_module
  
  use Geomechanics_Realization_class
  use Geomechanics_Global_module
  use Geomechanics_Force_module
  
  use Option_module
  use Waypoint_module
  
  implicit none
  
  type(simulation_type) :: simulation
  
  type(option_type), pointer :: option
  
  option => simulation%realization%option
  
  if (option%ngeomechdof > 0) then
    if (option%geomech_subsurf_coupling /= 0) then
      call GeomechCreateGeomechSubsurfVec(simulation%realization, &
                                          simulation%geomech_realization)
      call GeomechCreateSubsurfStressStrainVec(simulation%realization, &
                                               simulation%geomech_realization)

      call GeomechRealizMapSubsurfGeomechGrid(simulation%realization, &
                                              simulation%geomech_realization, &
                                              option)
    endif
    call GeomechRealizLocalizeRegions(simulation%geomech_realization)
    call GeomechRealizPassFieldPtrToPatch(simulation%geomech_realization)
    call GeomechRealizProcessMatProp(simulation%geomech_realization)
    call GeomechRealizProcessGeomechCouplers(simulation%geomech_realization)
    call GeomechRealizProcessGeomechConditions(simulation%geomech_realization)
    call InitGeomechMatPropToGeomechRegions(simulation%geomech_realization)
    call GeomechRealizInitAllCouplerAuxVars(simulation%geomech_realization)  
    call GeomechRealizPrintCouplers(simulation%geomech_realization)  
    call GeomechRealizAddWaypointsToList(simulation%geomech_realization)
    call GeomechGridElemSharedByNodes(simulation%geomech_realization)
    call WaypointListFillIn(option,simulation%geomech_realization%waypoint_list)
    call WaypointListRemoveExtraWaypnts(option, &
                                    simulation%geomech_realization%waypoint_list)
    call GeomechForceSetup(simulation%geomech_realization)
    call GeomechGlobalSetup(simulation%geomech_realization)
    
    ! SK: We are solving quasi-steady state solution for geomechanics.
    ! Initial condition is not needed, hence CondControlAssignFlowInitCondGeomech
    ! is not needed, at this point.
    call GeomechForceUpdateAuxVars(simulation%geomech_realization)
  endif
  
end subroutine InitGeomechSetupRealization

! ************************************************************************** !

subroutine InitGeomechSetupSolvers(geomech_realization,realization,solver)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_class
  use Geomechanics_Realization_class
  use Option_module
  
  use Solver_module
  use Convergence_module
  use Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Discretization_module
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscpc.h"
  
  type(geomech_realization_type) :: geomech_realization
  type(realization_type) :: realization
  type(solver_type), pointer :: solver

  type(option_type), pointer :: option
  type(convergence_context_type), pointer :: convergence_context
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  option => realization%option
  
  call printMsg(option,"  Beginning setup of GEOMECH SNES ")
    
  if (solver%J_mat_type == MATAIJ) then
    option%io_buffer = 'AIJ matrix not supported for geomechanics.'
    call printErrMsg(option)
  endif

  call SolverCreateSNES(solver,option%mycomm)  
  call SNESSetOptionsPrefix(solver%snes, "geomech_", &
                            ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)
        
  if (solver%Jpre_mat_type == '') then
    solver%Jpre_mat_type = solver%J_mat_type
  endif
  call GeomechDiscretizationCreateJacobian(geomech_realization% &
                                            geomech_discretization,NGEODOF, &
                                            solver%Jpre_mat_type, &
                                            solver%Jpre,option)

  solver%J = solver%Jpre
  call MatSetOptionsPrefix(solver%Jpre,"geomech_", &
                            ierr);CHKERRQ(ierr)
    

  call SNESSetFunction(solver%snes,geomech_realization%geomech_field%disp_r, &
                       GeomechForceResidual, &
                       realization,ierr);CHKERRQ(ierr)

  call SNESSetJacobian(solver%snes,solver%J, &
                       solver%Jpre,GeomechForceJacobian, &
                       realization,ierr);CHKERRQ(ierr)
  ! by default turn off line search
  call SNESGetLineSearch(solver%snes,linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC, &
                              ierr);CHKERRQ(ierr)

  ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
  if (option%verbosity >= 1) then
    string = '-geomech_snes_view'
    call PetscOptionsInsertString(string, ierr);CHKERRQ(ierr)
  endif

  call SolverSetSNESOptions(solver)

  option%io_buffer = 'Solver: ' // trim(solver%ksp_type)
  call printMsg(option)
  option%io_buffer = 'Preconditioner: ' // trim(solver%pc_type)
  call printMsg(option)

  ! shell for custom convergence test.  The default SNES convergence test
  ! is call within this function.
  !TODO(geh): free this convergence context somewhere!
  option%io_buffer = 'DEALLOCATE GEOMECH CONVERGENCE CONTEXT somewhere!!!'
  convergence_context => ConvergenceContextCreate(solver,option, &
                                                  realization%patch%grid)
  call SNESSetConvergenceTest(solver%snes,ConvergenceTest, &
                              convergence_context, &
                              PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)                                                  

  call printMsg(option,"  Finished setting up GEOMECH SNES ")
    
end subroutine InitGeomechSetupSolvers

! ************************************************************************** !

subroutine InitGeomechMatPropToGeomechRegions(geomech_realization)
  ! 
  ! This routine assigns geomech material
  ! properties to associated regions
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Strata_module
  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Option_module

  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: ivertex, local_id, ghosted_id, natural_id, geomech_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(geomech_strata_type), pointer :: strata
  type(geomech_patch_type), pointer :: patch  

  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: null_geomech_material_property
  type(gm_region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  PetscReal, pointer :: imech_loc_p(:)
  
  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch

  ! loop over all patches and allocation material id arrays
  if (.not.associated(patch%imat)) then
    allocate(patch%imat(patch%geomech_grid%ngmax_node))
    ! initialize to "unset"
    patch%imat = UNINITIALIZED_INTEGER
  endif

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  grid => patch%geomech_grid
  strata => patch%geomech_strata_list%first
  do
    if (.not.associated(strata)) exit
    ! Read in cell by cell material ids if they exist
    if (.not.associated(strata%region) .and. strata%active) then
      option%io_buffer = 'Reading of material prop from file for' // &
        ' geomech is not implemented.'
      call printErrMsgByRank(option)
    ! Otherwise, set based on region
    else if (strata%active) then
      update_ghosted_material_ids = PETSC_TRUE
      region => strata%region
      geomech_material_property => strata%material_property
      if (associated(region)) then
        istart = 1
        iend = region%num_verts
      else
        istart = 1
        iend = grid%nlmax_node
      endif
      do ivertex = istart, iend
        if (associated(region)) then
          local_id = region%vertex_ids(ivertex)
        else
          local_id = ivertex
        endif
        ghosted_id = grid%nL2G(local_id)
        patch%imat(ghosted_id) = geomech_material_property%id
      enddo
    endif
    strata => strata%next
  enddo
    
  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call GeomechRealizLocalToLocalWithArray(geomech_realization, &
                                            MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_geomech_material_property => GeomechanicsMaterialPropertyCreate()
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    ghosted_id = grid%nL2G(local_id)
    geomech_material_id = patch%imat(ghosted_id)
    if (geomech_material_id == 0) then ! accomodate inactive cells
      geomech_material_property = null_geomech_material_property
    else if ( geomech_material_id > 0 .and. &
              geomech_material_id <= &
              size(geomech_realization%geomech_material_property_array)) then
      geomech_material_property => &
         geomech_realization% &
           geomech_material_property_array(geomech_material_id)%ptr
      if (.not.associated(geomech_material_property)) then
        write(dataset_name,*) geomech_material_id
        option%io_buffer = 'No material property for geomech material id ' // &
                            trim(adjustl(dataset_name)) &
                            //  ' defined in input file.'
        call printErrMsgByRank(option)
      endif
    else if (Uninitialized(geomech_material_id)) then 
      write(dataset_name,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized geomech material id in patch at cell ' // &
                         trim(adjustl(dataset_name))
      call printErrMsgByRank(option)
    else if (geomech_material_id > size(geomech_realization% &
      geomech_material_property_array)) then
      write(option%io_buffer,*) geomech_material_id
      option%io_buffer = 'Unmatched geomech material id in patch:' // &
        adjustl(trim(option%io_buffer))
      call printErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with geomech material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn/Satish!'
      call printErrMsgByRank(option)
    endif
    imech_loc_p(ghosted_id) = geomech_material_property%id
  enddo ! local_id - loop
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  
  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)
  
  call GeomechDiscretizationLocalToLocal(geomech_discretization,field%imech_loc, &
                                         field%imech_loc,ONEDOF)
  
end subroutine InitGeomechMatPropToGeomechRegions

end module Init_Geomechanics_module
