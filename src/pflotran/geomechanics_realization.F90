#ifdef GEOMECH

module Geomechanics_Realization_module

  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Material_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Region_module
  use Geomechanics_Condition_module
  use Input_module
  use Option_module
  use Output_Aux_module
  use Waypoint_module
  use Dataset_Base_class
  use PFLOTRAN_Constants_module

 
  implicit none
  
private


#include "finclude/petscsys.h"

  type, public :: geomech_realization_type
    PetscInt :: id
    type(geomech_discretization_type), pointer        :: discretization
    type(input_type), pointer                         :: input
    type(geomech_patch_type), pointer                 :: geomech_patch
    type(option_type), pointer                        :: option
    type(geomech_material_property_type), &
                           pointer       :: geomech_material_properties
    type(geomech_material_property_ptr_type), &
                           pointer       :: geomech_material_property_array(:)
    type(waypoint_list_type), pointer                 :: waypoints
    type(geomech_field_type), pointer                 :: geomech_field
    type(geomech_debug_type), pointer                 :: debug
    type(output_option_type), pointer                 :: output_option
    type(gm_region_list_type), pointer                :: geomech_regions
    type(geomech_condition_list_type),pointer         :: geomech_conditions
    class(dataset_base_type), pointer                 :: geomech_datasets
  end type geomech_realization_type

public :: GeomechRealizCreate, &
          GeomechRealizDestroy, &
          GeomechRealizAddStrata, &
          GeomechRealizAddGeomechCoupler, &
          GeomechRealizLocalizeRegions, &
          GeomechRealizPassFieldPtrToPatch, &
          GeomechRealizProcessMatProp, &
          GeomechRealizProcessGeomechCouplers, &
          GeomechRealizCreateDiscretization, &
          GeomechRealizProcessGeomechConditions, &
          GeomechRealizInitAllCouplerAuxVars, &
          GeomechRealizPrintCouplers, &
          GeomechRealizAddWaypointsToList, &
          GeomechRealizGetDataset, &
          GeomechRealizLocalToLocalWithArray, &
          GeomechRealizMapSubsurfGeomechGrid, &
          GeomechGridElemSharedByNodes

contains

! ************************************************************************** !
!
! GeomechRealizCreate: This subroutine creates realization for geomechanics
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
function GeomechRealizCreate(option)

  implicit none

  type(geomech_realization_type), pointer    :: GeomechRealizCreate
  type(geomech_realization_type), pointer    :: geomech_realization
  type(option_type), pointer                 :: option
  
  allocate(geomech_realization)
  geomech_realization%id = 0
  if (associated(option)) then
    geomech_realization%option => option
  else
    geomech_realization%option => OptionCreate()
  endif
  
  nullify(geomech_realization%input)
  geomech_realization%discretization => GeomechDiscretizationCreate()
  
  geomech_realization%geomech_field => GeomechFieldCreate()
  geomech_realization%output_option => OutputOptionCreate()
  geomech_realization%debug => GeomechDebugCreate()
  
  allocate(geomech_realization%geomech_regions)
  call GeomechRegionInitList(geomech_realization%geomech_regions)
  
  allocate(geomech_realization%geomech_conditions)
  call GeomechConditionInitList(geomech_realization%geomech_conditions)

  nullify(geomech_realization%geomech_material_properties)
  nullify(geomech_realization%geomech_material_property_array)
  
  nullify(geomech_realization%geomech_patch)

  GeomechRealizCreate => geomech_realization
  
end function GeomechRealizCreate

! ************************************************************************** !
!
! GeomechRealizAddStrata: Adds strata to a list
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechRealizAddStrata(geomech_realization,strata)

  use Geomechanics_Strata_module

  implicit none
  
  type(geomech_realization_type)          :: geomech_realization
  type(geomech_strata_type), pointer      :: strata
  
  type(geomech_patch_type), pointer       :: geomech_patch
  type(geomech_strata_type), pointer      :: new_strata
  
  geomech_patch => geomech_realization%geomech_patch

  if (.not.associated(geomech_patch)) return
 
  new_strata => GeomechStrataCreate(strata)
  call GeomechStrataAddToList(new_strata,geomech_patch%geomech_strata)
  nullify(new_strata)
  
  call GeomechStrataDestroy(strata)
 
end subroutine GeomechRealizAddStrata

! ************************************************************************** !
!
! GeomechRealizLocalizeRegions: This routine localizes geomechanics regions
!                               within each patch
! author: Satish Karra, LANL
! date: 06/07/13
!
! ************************************************************************** !
subroutine GeomechRealizLocalizeRegions(geomech_realization)

  use Option_module
  use String_module

  implicit none
  
  type(geomech_realization_type)               :: geomech_realization
  type(geomech_patch_type), pointer            :: patch
  type(option_type), pointer                   :: option

  option => geomech_realization%option

  ! localize the regions on each patch
  patch => geomech_realization%geomech_patch
  call GeomechPatchLocalizeRegions(patch, &
                                   geomech_realization%geomech_regions, &
                                   option)
                                   
end subroutine GeomechRealizLocalizeRegions

! ************************************************************************** !
!
! GeomechRealizProcessMatProp: Setup for material properties
! author: Satish Karra, LANL
! date: 06/13/13
!
! ************************************************************************** !
subroutine GeomechRealizProcessMatProp(geomech_realization)

  use String_module
  
  implicit none
  
  type(geomech_realization_type)               :: geomech_realization
  type(geomech_patch_type), pointer            :: patch  
  type(option_type), pointer                   :: option

  
  option => geomech_realization%option
  
  ! organize lists
  call GeomechanicsMaterialPropConvertListToArray( &
                          geomech_realization%geomech_material_properties, &
                          geomech_realization%geomech_material_property_array, &
                          option)
  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  patch => geomech_realization%geomech_patch
  patch%geomech_material_properties => geomech_realization% &
                                       geomech_material_properties
  call GeomechanicsMaterialPropConvertListToArray( &
                                    patch%geomech_material_properties, &
                                    patch%geomech_material_property_array, &
                                    option)
                                      
end subroutine GeomechRealizProcessMatProp

! ************************************************************************** !
!
! GeomechRealizCreateDiscretization: Creates grid
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechRealizCreateDiscretization(realization)

  use Geomechanics_Grid_Aux_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_realization_type)                 :: realization
  type(geomech_discretization_type), pointer     :: discretization
  type(geomech_grid_type), pointer               :: grid
  type(option_type), pointer                     :: option
  type(geomech_field_type), pointer              :: geomech_field
  type(gmdm_ptr_type), pointer                   :: dm_ptr
  PetscErrorCode                                 :: ierr

  discretization => realization%discretization
  grid => discretization%grid
  option => realization%option
  geomech_field => realization%geomech_field
  
  call GeomechDiscretizationCreateDMs(discretization,option)
  
  ! n degree of freedom, global
  call GeomechDiscretizationCreateVector(discretization,NGEODOF,geomech_field%disp_xx, &
                                         GLOBAL,option)
  call VecSet(geomech_field%disp_xx,0.d0,ierr)

  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%disp_xx, &
                                            geomech_field%disp_r)
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%disp_xx, &
                                            geomech_field%work)
  
  ! 1 degree of freedom, global                                                                                    
  call GeomechDiscretizationCreateVector(discretization,ONEDOF,geomech_field%press, &
                                         GLOBAL,option)
  call VecSet(geomech_field%press,0.d0,ierr)
  
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%press, &
                                            geomech_field%temp)
                                            
  ! n degrees of freedom, local
  call GeomechDiscretizationCreateVector(discretization,NGEODOF,geomech_field%disp_xx_loc, &
                                         LOCAL,option)
  call VecSet(geomech_field%disp_xx_loc,0.d0,ierr)
 
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%disp_xx_loc, &
                                            geomech_field%work_loc)

  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%disp_xx_loc, &
                                            geomech_field%disp_xx_init_loc)
                                            
  ! 1 degree of freedom, local
  call GeomechDiscretizationCreateVector(discretization,ONEDOF,geomech_field%press_loc, &
                                         LOCAL,option)

  call VecSet(geomech_field%press_loc,0.d0,ierr)
  
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%press_loc, &
                                            geomech_field%temp_loc)

  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%press_loc, &
                                            geomech_field%press_init_loc)

  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%press_loc, &
                                            geomech_field%temp_init_loc)

  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%press_loc, &
                                            geomech_field%imech_loc)

  ! 6 dof for strain and stress
  call GeomechDiscretizationCreateVector(discretization,SIX_INTEGER,geomech_field%strain_loc, &
                                         LOCAL,option)

  call VecSet(geomech_field%strain_loc,0.d0,ierr)
 
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%strain_loc, &
                                            geomech_field%stress_loc)

  call GeomechDiscretizationCreateVector(discretization,SIX_INTEGER,geomech_field%strain, &
                                         GLOBAL,option)

  call VecSet(geomech_field%strain,0.d0,ierr)
 
  call GeomechDiscretizationDuplicateVector(discretization,geomech_field%strain, &
                                            geomech_field%stress) 

  grid => discretization%grid
  
  ! set up nG2L, NL2G, etc.
  call GMGridMapIndices(grid,discretization%dm_1dof%gmdm, &
                        grid%nG2L,grid%nL2G,grid%nG2A,option)
                        
  ! SK, Need to add a subroutine to ensure right hand rule
  ! SK, Need to add a subroutine equivalent to UGridComputeCoord                      

  
end subroutine GeomechRealizCreateDiscretization

! ************************************************************************** !
!
! GeomechRealizMapSubsurfGeomechGrid: This routine creates scatter contexts
! betweeen subsurface and geomech grids
! author: Satish Karra, LANL
! date: 09/09/13
!
! ************************************************************************** !
subroutine GeomechRealizMapSubsurfGeomechGrid(realization,geomech_realization, &
                                              option)

  use Option_module
  use Geomechanics_Grid_Aux_module
  use Realization_class
  use Grid_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscdm.h"  
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(realization_type), pointer              :: realization
  type(geomech_realization_type), pointer      :: geomech_realization
  type(geomech_grid_type), pointer             :: geomech_grid
  type(option_type)                            :: option
  type(grid_type), pointer                     :: grid
  type(gmdm_type), pointer                     :: gmdm
  type(gmdm_ptr_type), pointer                 :: dm_ptr
  IS                                           :: is_geomech, is_subsurf
  IS                                           :: is_subsurf_natural
  IS                                           :: is_subsurf_petsc
  PetscViewer                                  :: viewer
  PetscErrorCode                               :: ierr
  AO                                           :: ao_geomech_to_subsurf_natural
  PetscInt, allocatable                        :: int_array(:)
  PetscInt                                     :: local_id
  VecScatter                                   :: scatter
  IS                                           :: is_geomech_petsc
  PetscInt, pointer                            :: int_ptr(:)
  IS                                           :: is_geomech_petsc_block
  IS                                           :: is_subsurf_petsc_block

  geomech_grid => geomech_realization%discretization%grid
  grid => realization%discretization%grid    
    
  ! Convert from 1-based to 0-based  
  call ISCreateGeneral(option%mycomm,geomech_grid%mapping_num_cells, &
                       geomech_grid%mapping_cell_ids_flow-1, &
                       PETSC_COPY_VALUES,is_subsurf,ierr)
                       
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_mapping_cell_ids_flow.out', &
                            viewer,ierr)
  call ISView(is_subsurf,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   
 
  call ISCreateGeneral(option%mycomm,geomech_grid%mapping_num_cells, &
                       geomech_grid%mapping_vertex_ids_geomech-1, &
                       PETSC_COPY_VALUES,is_geomech,ierr)   

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_mapping_vertex_ids_geomech.out', &
                            viewer,ierr)
  call ISView(is_geomech,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   

  call AOCreateMappingIS(is_geomech,is_subsurf,ao_geomech_to_subsurf_natural,ierr)
  call ISDestroy(is_geomech,ierr)
  call ISDestroy(is_subsurf,ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_ao_geomech_to_subsurf_natural.out', &
                            viewer,ierr)
  call AOView(ao_geomech_to_subsurf_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  allocate(int_array(grid%nlmax))
  do local_id = 1, grid%nlmax
    int_array(local_id) = grid%nG2A(grid%nL2G(local_id)) - 1
  enddo

  call ISCreateGeneral(option%mycomm,grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_subsurf_natural,ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_subsurf_natural.out', &
                            viewer,ierr)
  call ISView(is_subsurf_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  allocate(int_array(grid%nlmax))
  do local_id = 1, grid%nlmax
    int_array(local_id) = (local_id-1) + grid%global_offset
  enddo

  call ISCreateGeneral(option%mycomm,grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_subsurf_petsc,ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_subsurf_petsc.out', &
                            viewer,ierr)
  call ISView(is_subsurf_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call ISDuplicate(is_subsurf_natural,is_geomech_petsc,ierr)
  call ISCopy(is_subsurf_natural,is_geomech_petsc,ierr)
  
  call AOPetscToApplicationIS(ao_geomech_to_subsurf_natural,is_geomech_petsc,ierr)
 
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_subsurf_petsc_geomech_natural.out', &
                            viewer,ierr)
  call ISView(is_geomech_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   

  call AOApplicationToPetscIS(geomech_grid%ao_natural_to_petsc_nodes, &
                              is_geomech_petsc,ierr)
                              
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_subsurf_petsc_geomech_petsc.out', &
                            viewer,ierr)
  call ISView(is_geomech_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif                              

  call VecScatterCreate(realization%field%porosity0,is_subsurf_petsc, &
                        geomech_realization%geomech_field%press, &
                        is_geomech_petsc,scatter,ierr)
                        
  if (ierr /= 0) then
    option%io_buffer = 'The number of cells specified in ' // &
                       'input file might not be same as the ' // &
                       'SUBSURF->GEOMECH mapping used.'
    call printErrMsg(option)
  endif

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_scatter_subsurf_to_geomech.out', &
                            viewer,ierr)
  call VecScatterView(scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   discretization,ONEDOF)

  call VecScatterCopy(scatter,dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof,ierr)
  
  call VecScatterDestroy(scatter,ierr)
  
  ! Geomech to subsurf scatter
  
  allocate(int_array(grid%nlmax))
  call ISGetIndicesF90(is_geomech_petsc,int_ptr,ierr)
  do local_id = 1, grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo  
  call ISRestoreIndicesF90(is_geomech_petsc,int_ptr,ierr)
  call ISCreateBlock(option%mycomm,SIX_INTEGER,grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,is_geomech_petsc_block,ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_geomech_petsc_block.out', &
                            viewer,ierr)
  call ISView(is_geomech_petsc_block,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif 

  allocate(int_array(grid%nlmax))
  call ISGetIndicesF90(is_subsurf_petsc,int_ptr,ierr)
  do local_id = 1, grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo  
  call ISRestoreIndicesF90(is_subsurf_petsc,int_ptr,ierr)
  call ISCreateBlock(option%mycomm,SIX_INTEGER,grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,is_subsurf_petsc_block,ierr)
  deallocate(int_array)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_is_subsurf_petsc_block.out', &
                            viewer,ierr)
  call ISView(is_subsurf_petsc_block,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call VecScatterCreate(geomech_realization%geomech_field%strain,is_geomech_petsc_block, &
                        geomech_realization%geomech_field%strain_subsurf, &
                        is_subsurf_petsc_block,scatter,ierr)
                        
  if (ierr /= 0) then
    option%io_buffer = 'The number of cells specified in ' // &
                       'input file might not be same as the ' // &
                       'GEOMECH->SUBSURF mapping used.'
    call printErrMsg(option)
  endif

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'geomech_scatter_geomech_to_subsurf_block.out', &
                            viewer,ierr)
  call VecScatterView(scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   discretization,ONEDOF)

  call VecScatterCopy(scatter,dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof,ierr)
  
  call ISDestroy(is_geomech,ierr)
  call ISDestroy(is_subsurf,ierr)
  call ISDestroy(is_subsurf_natural,ierr)
  call ISDestroy(is_geomech_petsc,ierr)
  call ISDestroy(is_subsurf_petsc,ierr)
  call AODestroy(ao_geomech_to_subsurf_natural,ierr) 
  call ISDestroy(is_subsurf_petsc_block,ierr)
  call ISDestroy(is_geomech_petsc_block,ierr)

end subroutine GeomechRealizMapSubsurfGeomechGrid

! ************************************************************************** !
!
! GeomechGridElemsSharedByNodes: Calculates the number of elements common
! to a node (vertex)
! author: Satish Karra
! date: 09/17/13
!
! ************************************************************************** !
subroutine GeomechGridElemSharedByNodes(realization)
  
  use Geomechanics_Grid_Aux_module

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_realization_type) :: realization
  type(geomech_grid_type), pointer :: grid
  
  PetscInt :: ielem
  PetscInt :: ivertex
  PetscInt :: ghosted_id
  PetscInt, allocatable :: elenodes(:)
  PetscReal, pointer :: elem_sharing_node_loc_p(:)
  PetscErrorCode :: ierr
  
  grid => realization%discretization%grid
  
  call VecGetArrayF90(grid%no_elems_sharing_node_loc,elem_sharing_node_loc_p, &
                      ierr)
  
  do ielem = 1, grid%nlmax_elem
    elenodes = grid%elem_nodes(1:grid%elem_nodes(0,ielem),ielem)
    do ivertex = 1, grid%elem_nodes(0,ielem)
      ghosted_id = elenodes(ivertex) 
      elem_sharing_node_loc_p(ghosted_id) = &
        elem_sharing_node_loc_p(ghosted_id) + 1
    enddo
  enddo
  
  call VecRestoreArrayF90(grid%no_elems_sharing_node_loc, &
                         elem_sharing_node_loc_p,ierr)

  
  ! Local to global scatter
  call GeomechDiscretizationLocalToGlobalAdd(realization%discretization, &
                                             grid%no_elems_sharing_node_loc, &
                                             grid%no_elems_sharing_node, &
                                             ONEDOF)  
                                             
end subroutine GeomechGridElemSharedByNodes

! ************************************************************************** !
!
! GeomechRealizInitAllCouplerAuxVars: This routine initializez coupler 
!                                     auxillary variables
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizInitAllCouplerAuxVars(geomech_realization)

  use Option_module

  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  type(geomech_patch_type), pointer :: patch
  
  patch => geomech_realization%geomech_patch

  call GeomechPatchInitAllCouplerAuxVars(patch,geomech_realization%option)

end subroutine GeomechRealizInitAllCouplerAuxVars

! ************************************************************************** !
!
! GeomechRealizLocalToLocalWithArray: This routine takes an F90 array that is 
!                                    ghosted and updates the ghosted values
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizLocalToLocalWithArray(geomech_realization,array_id)

  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module

  implicit none

  type(geomech_realization_type)            :: geomech_realization
  PetscInt                                  :: array_id
  
  type(geomech_patch_type), pointer         :: patch
  type(geomech_grid_type), pointer          :: grid
  type(geomech_field_type), pointer         :: geomech_field

  geomech_field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GeomechGridCopyIntegerArrayToVec(grid,patch%imat, &
                                            geomech_field%work_loc, &
                                            grid%ngmax_node)
  end select
    
  call GeomechDiscretizationLocalToLocal(geomech_realization%discretization, &
                                         geomech_field%work_loc, &
                                         geomech_field%work_loc,ONEDOF)
                                  
  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GeomechGridCopyVecToIntegerArray(grid,patch%imat, &
                                            geomech_field%work_loc, &
                                            grid%ngmax_node)
  end select
  
end subroutine GeomechRealizLocalToLocalWithArray

! ************************************************************************** !
!
! GeomechRealizPrintCouplers: Print boundary data for geomechanics
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizPrintCouplers(realization)

  use Geomechanics_Coupler_module
  
  implicit none
  
  type(geomech_realization_type) :: realization
  
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: cur_coupler
  type(option_type), pointer :: option
 
  option => realization%option
 
  if (.not.OptionPrintToFile(option)) return
  
  patch => realization%geomech_patch
   
  cur_coupler => patch%geomech_boundary_conditions%first
  do
    if (.not.associated(cur_coupler)) exit
    call GeomechRealizPrintCoupler(cur_coupler,option)    
    cur_coupler => cur_coupler%next
  enddo
     
  cur_coupler => patch%geomech_source_sinks%first
  do
    if (.not.associated(cur_coupler)) exit
    call GeomechRealizPrintCoupler(cur_coupler,option)    
    cur_coupler => cur_coupler%next
  enddo
 
end subroutine GeomechRealizPrintCouplers

! ************************************************************************** !
!
! GeomechRealizPrintCoupler: Prints boundary condition coupler for geomechanics
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizPrintCoupler(coupler,option)

  use Geomechanics_Coupler_module
  
  implicit none
  
  type(geomech_coupler_type) :: coupler
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(geomech_condition_type), pointer :: geomech_condition
  type(gm_region_type), pointer :: region
   
98 format(40('=+'))
99 format(80('-'))
  
  geomech_condition => coupler%geomech_condition
  region => coupler%region

  write(option%fid_out,*)
  write(option%fid_out,98)


  select case(coupler%itype)
    case(GM_BOUNDARY_COUPLER_TYPE)
      string = 'Geomech Boundary Condition'
    case(GM_SRC_SINK_COUPLER_TYPE)
      string = 'Geomech Source Sink'
  end select
  write(option%fid_out,'(/,2x,a,/)') trim(string)

  write(option%fid_out,99)
101 format(5x,'     Geomech Condition: ',2x,a)
  if (associated(geomech_condition)) &
    write(option%fid_out,101) trim(geomech_condition%name)
102 format(5x,'             Region: ',2x,a)
  if (associated(region)) &
    write(option%fid_out,102) trim(region%name)
  write(option%fid_out,99)
  
  if (associated(geomech_condition)) then
    call GeomechConditionPrint(geomech_condition,option)
  endif
 
end subroutine GeomechRealizPrintCoupler

! ************************************************************************** !
!
! GeomechRealizPassFieldPtrToPatch: This subroutine passes field to patch
! author: Satish Karra, LANL
! date: 06/13/13
!
! ************************************************************************** !
subroutine GeomechRealizPassFieldPtrToPatch(geomech_realization)

  use Option_module

  implicit none
  
  type(geomech_realization_type)           :: geomech_realization

  type(geomech_patch_type), pointer        :: patch

  patch => geomech_realization%geomech_patch
   
  patch%geomech_field => geomech_realization%geomech_field
  
end subroutine GeomechRealizPassFieldPtrToPatch

! ************************************************************************** !
!
! GeomechRealizProcessGeomechCouplers: This subroutine sets up couplers in 
!                                      geomech realization
! author: Satish Karra, LANL
! date: 06/14/13
!
! ************************************************************************** !
subroutine GeomechRealizProcessGeomechCouplers(geomech_realization)

  implicit none
  
  type(geomech_realization_type)           :: geomech_realization

  type(geomech_patch_type), pointer        :: patch
  
  patch => geomech_realization%geomech_patch
  
  call GeomechPatchProcessGeomechCouplers(patch, &
                                   geomech_realization%geomech_conditions, &
                                   geomech_realization%option)
 
end subroutine GeomechRealizProcessGeomechCouplers

! ************************************************************************** !
!
! GeomechRealizProcessGeomechConditions: This subroutine sets up condition in 
!                                      geomech realization
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizProcessGeomechConditions(geomech_realization)

  use Dataset_Base_class
  use Dataset_module
  

  implicit none

  type(geomech_realization_type), pointer   :: geomech_realization
  
  type(geomech_condition_type), pointer     :: cur_geomech_condition
  type(geomech_sub_condition_type), pointer :: cur_geomech_sub_condition
  type(option_type), pointer                :: option
  character(len=MAXSTRINGLENGTH)            :: string
  character(len=MAXWORDLENGTH)              :: dataset_name
  class(dataset_base_type), pointer         :: dataset
  PetscInt                                  :: i
  
  option => geomech_realization%option
  
  ! loop over geomech conditions looking for linkage to datasets
  cur_geomech_condition => geomech_realization%geomech_conditions%first
  do
    if (.not.associated(cur_geomech_condition)) exit
      do i = 1, size(cur_geomech_condition%sub_condition_ptr)
        ! check for dataset in dataset
        if (associated(cur_geomech_condition%sub_condition_ptr(i)%ptr% &
                        geomech_dataset%dataset)) then
          dataset_name = cur_geomech_condition%sub_condition_ptr(i)%ptr% &
                        geomech_dataset%dataset%name
          ! delete the dataset since it is solely a placeholder
          call DatasetDestroy(cur_geomech_condition%sub_condition_ptr(i)%ptr% &
                              geomech_dataset%dataset)
          ! get dataset from list
          string = 'geomech_condition ' // trim(cur_geomech_condition%name)
          dataset => &
            DatasetBaseGetPointer(geomech_realization%geomech_datasets, &
                                  dataset_name,string,option)
          cur_geomech_condition%sub_condition_ptr(i)%ptr%geomech_dataset &
            %dataset => dataset
        endif
      enddo
     cur_geomech_condition => cur_geomech_condition%next
  enddo

end subroutine GeomechRealizProcessGeomechConditions


! ************************************************************************** !
!
! GeomechRealizGetDataset: This routine extracts variables indexed by 
!                          ivar and isubvar from geomechanics realization
! author: Satish Karra, LANL
! date: 07/03/13
!
! ************************************************************************** !
subroutine GeomechRealizGetDataset(geomech_realization,vec,ivar,isubvar,isubvar1)

  implicit none

  type(geomech_realization_type) :: geomech_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call GeomechPatchGetDataset(geomech_realization%geomech_patch, &
                              geomech_realization%geomech_field, &
                              geomech_realization%option, &
                              geomech_realization%output_option, &
                              vec,ivar,isubvar,isubvar1)

end subroutine GeomechRealizGetDataset

! ************************************************************************** !
!
! GeomechRealizAddGeomechCoupler: This subroutine addes a geomechanics 
!                                 coupler to a geomechanics realization
! author: Satish Karra, LANL
! date: 06/13/13
!
! ************************************************************************** !
subroutine GeomechRealizAddGeomechCoupler(realization,coupler)

  use Geomechanics_Coupler_module

  implicit none
  
  type(geomech_realization_type)                   :: realization
  type(geomech_coupler_type), pointer              :: coupler
  
  type(geomech_patch_type), pointer                :: patch
  type(geomech_coupler_type), pointer              :: new_coupler
  
  patch => realization%geomech_patch
  
 ! only add to geomech list for now, since they will be split out later
  new_coupler => GeomechCouplerCreate(coupler)
  select case(coupler%itype)
    case(GM_BOUNDARY_COUPLER_TYPE)
      call GeomechCouplerAddToList(new_coupler, &
                                   patch%geomech_boundary_conditions)
    case(GM_SRC_SINK_COUPLER_TYPE)
      call GeomechCouplerAddToList(new_coupler,patch%geomech_source_sinks)
  end select
  nullify(new_coupler)
  
  call GeomechCouplerDestroy(coupler)
 
end subroutine GeomechRealizAddGeomechCoupler

! ************************************************************************** !
!
! GeomechRealizAddWaypointsToList: Adds waypoints from BCs and source/sink
!                                  to waypoint list
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechRealizAddWaypointsToList(geomech_realization)

  use Option_module
  use Waypoint_module

  implicit none
  
  type(geomech_realization_type) :: geomech_realization
  
  type(waypoint_list_type), pointer :: waypoint_list
  type(geomech_condition_type), pointer :: cur_geomech_condition
  type(geomech_sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => geomech_realization%option
  waypoint_list => geomech_realization%waypoints
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_output = geomech_realization%output_option%print_final
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in GeomechRealizAddWaypointsToList'
    call printErrMsg(option)
  endif

  ! add update of geomech conditions
  cur_geomech_condition => geomech_realization%geomech_conditions%first
  do
    if (.not.associated(cur_geomech_condition)) exit
    if (cur_geomech_condition%sync_time_with_update) then
      do isub_condition = 1, cur_geomech_condition%num_sub_conditions
        sub_condition => cur_geomech_condition% &
                         sub_condition_ptr(isub_condition)%ptr
        call GeomechConditionDatasetGetTimes(option,sub_condition,final_time, &
                                          times)
        if (size(times) > 1000) then
          option%io_buffer = 'For geomech condition "' // &
            trim(cur_geomech_condition%name) // &
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
    cur_geomech_condition => cur_geomech_condition%next
  enddo
      
  ! add waypoints for periodic output
  if (geomech_realization%output_option%periodic_output_time_incr > 0.d0 .or. &
      geomech_realization%output_option%periodic_tr_output_time_incr > 0.d0) then

    if (geomech_realization%output_option%periodic_output_time_incr > 0.d0) then
      ! standard output
      temp_real = 0.d0
      do
        temp_real = temp_real + geomech_realization% &
                      output_option%periodic_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,geomech_realization%waypoints)
      enddo
    endif
    
    if (geomech_realization%output_option%periodic_tr_output_time_incr > 0.d0) then
      ! transient observation output
      temp_real = 0.d0
      do
        temp_real = temp_real + geomech_realization%output_option%periodic_tr_output_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_tr_output = PETSC_TRUE
        call WaypointInsertInList(waypoint,geomech_realization%waypoints)
      enddo
    endif

  endif

end subroutine GeomechRealizAddWaypointsToList

! ************************************************************************** !
!
! GeomechRealizDestroy: This subroutine deallocates geomechanics realization
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechRealizDestroy(geomech_realization)

  implicit none
  
  type(geomech_realization_type), pointer :: geomech_realization
  
  if(.not.associated(geomech_realization)) return
  
  call GeomechFieldDestroy(geomech_realization%geomech_field)
  
! sk: need to put this back. Crashes for now.  
!  call OutputOptionDestroy(geomech_realization%output_option)
  
  call GeomechRegionDestroyList(geomech_realization%geomech_regions)
  call GeomechConditionDestroyList(geomech_realization%geomech_conditions)
  
  if (associated(geomech_realization%debug)) &
    deallocate(geomech_realization%debug)
  nullify(geomech_realization%debug)  
  
  if (associated(geomech_realization%geomech_material_property_array)) &
    deallocate(geomech_realization%geomech_material_property_array)
  nullify(geomech_realization%geomech_material_property_array)
  if (associated(geomech_realization%geomech_patch)) &
    call GeomechanicsPatchDestroy(geomech_realization%geomech_patch)                                       
  call GeomechanicsMaterialPropertyDestroy(geomech_realization% &
                                           geomech_material_properties)
  call GeomechDiscretizationDestroy(geomech_realization%discretization)
  
end subroutine GeomechRealizDestroy

end module Geomechanics_Realization_module
#endif
! GEOMECH
