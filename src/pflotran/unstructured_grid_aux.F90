module Unstructured_Grid_Aux_module

!  use Connection_module
  use Unstructured_Cell_module
  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private 
  
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#if defined(SCORPIO)
  include "scorpiof.h"
#endif
  
  type, public :: unstructured_grid_type
    ! variables for all unstructured grids
    PetscInt :: num_ghost_cells   ! number of ghost cells (only) on processor
    PetscInt :: global_offset ! offset in petsc ordering for the first cell on a processor???
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash
    PetscInt, pointer :: cell_ids_natural(:) ! natural 1d right-hand i,j,k ordering
    PetscInt, pointer :: cell_ids_petsc(:) ! petsc ordering of cell ids
    PetscInt, pointer :: ghost_cell_ids_petsc(:) ! petsc ordering of ghost cells ids
    AO :: ao_natural_to_petsc ! mapping of natural to petsc ordering
    type(unstructured_explicit_type), pointer :: explicit_grid
    type(unstructured_polyhedra_type), pointer :: polyhedra_grid
    ! variables for implicit unstructured grids
    PetscInt :: grid_type         ! 3D subsurface (default) or 2D surface grid
    PetscInt :: num_vertices_global ! number of vertices in entire problem domain
    PetscInt :: num_vertices_local  ! number of vertices in local grid cells
    PetscInt :: num_vertices_natural ! number of vertices read initially
    PetscInt :: max_ndual_per_cell
    PetscInt :: max_nvert_per_cell
    PetscInt :: max_cells_sharing_a_vertex
    PetscInt, pointer :: cell_type(:)
    PetscInt, pointer :: cell_vertices(:,:) ! vertices for each grid cell (NO LONGER zero-based)
    PetscInt, pointer :: face_to_cell_ghosted(:,:) !
    PetscInt, pointer :: connection_to_face(:)
!geh: Should not need face_to_vertex_nindex() as one could use face_to_vertex() 
!     and vertex_ids_nindex() to get the same result.
!gb: face_to_vertex_natural is required in GridLocalizeRegionsForUGrid() and needs
!   to be saved because:
!   (i) face_to_vertex() - Removes the duplicate faces, and the algorithm in
!       GridLocalizeRegionsForUGrid() assumes presence of ALL faces.
!   (ii) More importantly, face_to_vertex() eventually has vertex indices in
!       local index (different from natural index, when using multiple processors).
!   Note: A region in the inputfile will be described in terms of natural index.
    PetscInt, pointer :: face_to_vertex_natural(:,:)
    PetscInt, pointer :: face_to_vertex(:,:)
    PetscInt, pointer :: cell_to_face_ghosted(:,:)
    PetscInt, pointer :: vertex_ids_natural(:)
    PetscInt, pointer :: cell_neighbors_local_ghosted(:,:) ! local neighbors
    type(point_type), pointer :: vertices(:)
    type(point_type), pointer :: face_centroid(:)
    PetscReal, pointer :: face_area(:)
  end type unstructured_grid_type
  
  type, public :: unstructured_explicit_type
    PetscInt, pointer :: cell_ids(:)
    PetscReal, pointer :: cell_volumes(:)
    type(point3d_type), pointer :: cell_centroids(:)
    PetscInt, pointer :: connections(:,:)
    PetscReal, pointer :: face_areas(:)
    type(point3d_type), pointer :: face_centroids(:)
    PetscInt :: num_cells_global  ! Number of cells in the entire domain
    PetscInt :: num_elems
    PetscInt :: num_elems_local   ! Number of elements locally
    PetscInt, pointer :: cell_connectivity(:,:)   
    type(point3d_type), pointer :: vertex_coordinates(:)
  end type unstructured_explicit_type

  type, public :: unstructured_polyhedra_type
    PetscInt, pointer :: cell_ids(:)
    PetscInt, pointer :: cell_nfaces(:)
    PetscInt, pointer :: cell_nverts(:)
    PetscInt, pointer :: cell_faceids(:,:)
    PetscInt, pointer :: cell_vertids(:,:)
    PetscReal, pointer :: cell_volumes(:)
    type(point3d_type), pointer :: cell_centroids(:)
    PetscInt, pointer :: face_ids(:)
    PetscInt, pointer :: face_cellids(:)
    PetscInt, pointer :: face_nverts(:)
    PetscInt, pointer :: face_vertids(:,:)
    PetscReal, pointer :: face_areas(:)
    type(point3d_type), pointer :: face_centroids(:)
    type(point3d_type), pointer :: vertex_coordinates(:)
    PetscInt :: num_cells_global
    PetscInt :: num_cells_local
    PetscInt :: num_faces_global
    PetscInt :: num_faces_local
    PetscInt :: num_vertices_global
    PetscInt :: num_vertices_local
    PetscInt :: max_nface_per_cell
    PetscInt :: max_nvert_per_face
    PetscInt :: max_nvert_per_cell
  end type unstructured_polyhedra_type

  type, public :: ugdm_type
    ! local: included both local (non-ghosted) and ghosted cells
    ! global: includes only local (non-ghosted) cells
    PetscInt :: ndof
    ! for the below
    ! ghosted = local (non-ghosted) and ghosted cells
    ! local = local (non-ghosted) cells
    IS :: is_ghosted_local ! IS for ghosted cells with local on-processor numbering
    IS :: is_local_local ! IS for local cells with local on-processor numbering
    IS :: is_ghosted_petsc ! IS for ghosted cells with petsc numbering
    IS :: is_local_petsc ! IS for local cells with petsc numbering
    IS :: is_ghosts_local ! IS for ghost cells with local on-processor numbering
    IS :: is_ghosts_petsc ! IS for ghost cells with petsc numbering
    IS :: is_local_natural ! IS for local cells with natural (global) numbering
    VecScatter :: scatter_ltog ! scatter context for local to global updates
    VecScatter :: scatter_gtol ! scatter context for global to local updates
    VecScatter :: scatter_ltol ! scatter context for local to local updates
    VecScatter :: scatter_gton ! scatter context for global to natural updates
    VecScatter :: scatter_ntog ! scatter context for natural to global updates
    ISLocalToGlobalMapping :: mapping_ltog  ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb ! block form of mapping_ltog
    Vec :: global_vec ! global vec (no ghost cells), petsc-ordering
    Vec :: local_vec ! local vec (includes local and ghosted cells), local ordering
#ifdef SURFACE_FLOW
    VecScatter :: scatter_bet_grids ! scatter context between surface and subsurface
                                    ! grids
    VecScatter :: scatter_bet_grids_1dof ! scatter context between surface and
                                         ! subsurface grids for 1-DOF
    VecScatter :: scatter_bet_grids_ndof ! scatter context between surface and
                                         ! subsurface grids for N-DOFs
#endif
  end type ugdm_type

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: UGridCreate, &
            UGridExplicitCreate, &
            UGridPolyhedraCreate, &
            UGridMapIndices, &
            UGridDMCreateJacobian, &
            UGridDMCreateVector, &
            UGridDestroy, &
            UGridCreateUGDM, &
            UGridDMDestroy, &
            UGridPartition, &
            UGridNaturalToPetsc, &
            UGridCreateOldVec, &
            UGridExplicitDestroy

contains

! ************************************************************************** !
!
! UGDMCreate: Creates an unstructured grid distributed mesh object
! author: Glenn Hammond
! date: 10/21/09
!
! ************************************************************************** !
function UGDMCreate()

  implicit none
  
  type(ugdm_type), pointer :: UGDMCreate

  type(ugdm_type), pointer :: ugdm

  allocate(ugdm)
  ugdm%is_ghosted_local = 0
  ugdm%is_local_local = 0
  ugdm%is_ghosted_petsc = 0
  ugdm%is_local_petsc = 0
  ugdm%is_ghosts_local = 0
  ugdm%is_ghosts_petsc = 0
  ugdm%is_local_natural = 0
  ugdm%scatter_ltog = 0
  ugdm%scatter_gtol = 0
  ugdm%scatter_ltol  = 0
  ugdm%scatter_gton = 0
  ugdm%scatter_ntog = 0
  ugdm%mapping_ltog = 0
  ugdm%mapping_ltogb = 0
  ugdm%global_vec = 0
  ugdm%local_vec = 0
#ifdef SURFACE_FLOW
  ugdm%scatter_bet_grids = 0
  ugdm%scatter_bet_grids_1dof = 0
  ugdm%scatter_bet_grids_ndof = 0
#endif
  UGDMCreate => ugdm

end function UGDMCreate

! ************************************************************************** !
!
! UGridCreate: Creates an unstructured grid object
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
function UGridCreate()

  implicit none
  
  type(unstructured_grid_type), pointer :: UGridCreate

  type(unstructured_grid_type), pointer :: unstructured_grid

  allocate(unstructured_grid)

  ! variables for all unstructured grids
  unstructured_grid%num_ghost_cells = 0
  unstructured_grid%global_offset = 0
  unstructured_grid%nmax = 0
  unstructured_grid%nlmax = 0
  unstructured_grid%ngmax = 0
  nullify(unstructured_grid%hash)
  unstructured_grid%num_hash = 100
  nullify(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  unstructured_grid%ao_natural_to_petsc = 0
  nullify(unstructured_grid%explicit_grid)
  nullify(unstructured_grid%polyhedra_grid)

  ! variables for implicit unstructured grids
  unstructured_grid%grid_type = THREE_DIM_GRID
  unstructured_grid%num_vertices_global = 0
  unstructured_grid%num_vertices_local = 0
  unstructured_grid%max_ndual_per_cell = 0
  unstructured_grid%max_nvert_per_cell = 0
  unstructured_grid%max_cells_sharing_a_vertex = 24
  nullify(unstructured_grid%cell_type)
  nullify(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%cell_to_face_ghosted)
  nullify(unstructured_grid%vertex_ids_natural)
  nullify(unstructured_grid%vertices)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
  nullify(unstructured_grid%connection_to_face)
  nullify(unstructured_grid%face_centroid)
  nullify(unstructured_grid%face_area)
  
  UGridCreate => unstructured_grid
  
end function UGridCreate

! ************************************************************************** !
!
! UGridExplicitCreate: Creates an explicit unstructured grid object
! author: Glenn Hammond
! date: 05/14/12
!
! ************************************************************************** !
function UGridExplicitCreate()

  implicit none
  
  type(unstructured_explicit_type), pointer :: UGridExplicitCreate

  type(unstructured_explicit_type), pointer :: explicit_grid

  allocate(explicit_grid)

  nullify(explicit_grid%cell_ids)
  nullify(explicit_grid%cell_volumes)
  nullify(explicit_grid%cell_centroids)
  nullify(explicit_grid%connections)
  nullify(explicit_grid%face_areas)
  nullify(explicit_grid%face_centroids)

  UGridExplicitCreate => explicit_grid
  
end function UGridExplicitCreate

! ************************************************************************** !
!> This routine creates a polyhedra unstructured grid object.
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 09/29/13
! ************************************************************************** !
function UGridPolyhedraCreate()

  implicit none

  type(unstructured_polyhedra_type), pointer :: UGridPolyhedraCreate

  type(unstructured_polyhedra_type), pointer :: polyhedra_grid

  allocate(polyhedra_grid)

  polyhedra_grid%num_cells_global = 0
  polyhedra_grid%num_cells_local = 0
  polyhedra_grid%num_faces_global = 0
  polyhedra_grid%num_faces_local = 0
  polyhedra_grid%num_vertices_global = 0
  polyhedra_grid%num_vertices_local = 0
  polyhedra_grid%max_nface_per_cell = 0
  polyhedra_grid%max_nvert_per_face = 0
  polyhedra_grid%max_nvert_per_cell = 0

  nullify(polyhedra_grid%cell_ids)
  nullify(polyhedra_grid%cell_nfaces)
  nullify(polyhedra_grid%cell_vertids)
  nullify(polyhedra_grid%cell_faceids)
  nullify(polyhedra_grid%cell_nverts)
  nullify(polyhedra_grid%cell_volumes)
  nullify(polyhedra_grid%cell_centroids)
  nullify(polyhedra_grid%face_ids)
  nullify(polyhedra_grid%face_cellids)
  nullify(polyhedra_grid%face_nverts)
  nullify(polyhedra_grid%face_vertids)
  nullify(polyhedra_grid%face_areas)
  nullify(polyhedra_grid%face_centroids)
  nullify(polyhedra_grid%vertex_coordinates)

  UGridPolyhedraCreate => polyhedra_grid

end function UGridPolyhedraCreate

! ************************************************************************** !
!
! UGridCreateUGDM: Constructs mappings / scatter contexts for PETSc DM 
!                                                  object
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UGridCreateUGDM(unstructured_grid,ugdm,ndof,option)
  
  use Option_module
  use Utility_module, only: reallocateIntArray
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscdm.h"  
#include "finclude/petscdm.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type), pointer :: ugdm
  PetscInt :: ndof
  type(option_type) :: option
  
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof
  IS :: is_tmp
  Vec :: vec_tmp
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: ndof_word
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscViewer :: viewer

  PetscInt, allocatable :: int_array(:)
  
  ugdm => UGDMCreate()
  ugdm%ndof = ndof

#if UGRID_DEBUG
  write(ndof_word,*) ndof
  ndof_word = adjustl(ndof_word)
  ndof_word = '_' // trim(ndof_word)
  string = 'Vectors' // ndof_word
  call printMsg(option,string)
#endif

  ! create global vec
  !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
  !                  PETSC_DETERMINE,ugdm%global_vec,ierr)
  call VecCreate(option%mycomm,ugdm%global_vec,ierr)
  call VecSetSizes(ugdm%global_vec,unstructured_grid%nlmax*ndof, &
                  PETSC_DECIDE,ierr)  
  call VecSetBlockSize(ugdm%global_vec,ndof,ierr)
  call VecSetFromOptions(ugdm%global_vec,ierr)

  ! create local vec
  !call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax*ndof, &
  !                  ugdm%local_vec,ierr)
  call VecCreate(PETSC_COMM_SELF,ugdm%local_vec,ierr)
  call VecSetSizes(ugdm%local_vec,unstructured_grid%ngmax*ndof,PETSC_DECIDE,ierr)
  call VecSetBlockSize(ugdm%local_vec,ndof,ierr)
  call VecSetFromOptions(ugdm%local_vec,ierr)
  
  ! IS for global numbering of local, non-ghosted cells
!geh  call VecGetOwnershipRange(ugdm%global_vec,istart,iend,ierr)
  ! ISCreateBlock requires block ids, not indices.  Therefore, istart should be
  ! the offset of the block from the beginning of the vector.
!geh  istart = istart / ndof
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax
 !geh   int_array(local_id) = (local_id-1)+istart
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo

  ! arguments for ISCreateBlock():
  ! option%mycomm  - the MPI communicator
  ! ndof  - number of elements in each block
  ! unstructured_grid%nlmax  - the length of the index set
  !                                      (the number of blocks
  ! int_array  - the list of integers, one for each block and count
  !              of block not indices
  ! PETSC_COPY_VALUES  - see PetscCopyMode, only PETSC_COPY_VALUES and
  !                      PETSC_OWN_POINTER are supported in this routine
  ! ugdm%is_local_petsc - the new index set
  ! ierr - PETScErrorCode
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_petsc,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_local_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do ghosted_id = 1, unstructured_grid%num_ghost_cells
    int_array(ghosted_id) = (ghosted_id+unstructured_grid%nlmax-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_local,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosts_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_ghosts_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
#if UGRID_DEBUG
  string = 'Index Sets' // ndof_word
  call printMsg(option,'Index Sets')
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do ghosted_id = 1, unstructured_grid%num_ghost_cells
    int_array(ghosted_id) = &
      (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_petsc,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosts_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_ghosts_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! IS for local numbering of local, non-ghosted cells
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax
    int_array(local_id) = (local_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_local,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_local_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! IS for ghosted numbering of local ghosted cells
  allocate(int_array(unstructured_grid%ngmax))
  do ghosted_id = 1, unstructured_grid%ngmax
    int_array(ghosted_id) = (ghosted_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_local,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosted_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_ghosted_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
             
  ! IS for petsc numbering of local ghosted cells
  allocate(int_array(unstructured_grid%ngmax))
  do local_id = 1, unstructured_grid%nlmax
!geh    int_array(local_id) = istart+(local_id-1)
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo
  do ghosted_id = 1,unstructured_grid%num_ghost_cells
    int_array(unstructured_grid%nlmax+ghosted_id) = &
      (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_petsc,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosted_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_ghosted_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif    
                 
  ! create a local to global mapping
#if UGRID_DEBUG
  string = 'ISLocalToGlobalMapping' // ndof_word
  call printMsg(option,string)
#endif

  call ISLocalToGlobalMappingCreateIS(ugdm%is_ghosted_petsc, &
                                      ugdm%mapping_ltog,ierr)

#if UGRID_DEBUG
  string = 'mapping_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISLocalToGlobalMappingView(ugdm%mapping_ltog,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
               
  ! create a block local to global mapping 
#if UGRID_DEBUG
  string = 'ISLocalToGlobalMappingBlock' // ndof_word
  call printMsg(option,string)
#endif

  call ISLocalToGlobalMappingBlock(ugdm%mapping_ltog,ndof, &
                                   ugdm%mapping_ltogb,ierr)
                                      
#if UGRID_DEBUG
  string = 'mapping_ltogb' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISLocalToGlobalMappingView(ugdm%mapping_ltogb,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  string = 'local to global' // ndof_word
  call printMsg(option,string)
#endif

  ! Create local to global scatter
  call VecScatterCreate(ugdm%local_vec,ugdm%is_local_local,ugdm%global_vec, &
                        ugdm%is_local_petsc,ugdm%scatter_ltog,ierr)
                        
#if UGRID_DEBUG
  string = 'scatter_ltog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_ltog,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  string = 'global to local' // ndof_word
  call printMsg(option,string)
#endif

  ! Create global to local scatter
  call VecScatterCreate(ugdm%global_vec,ugdm%is_ghosted_petsc,ugdm%local_vec, &
                        ugdm%is_ghosted_local,ugdm%scatter_gtol,ierr)
                        
#if UGRID_DEBUG
  string = 'scatter_gtol' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_gtol,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  string = 'local to local' // ndof_word
  call printMsg(option,string)
#endif
  
  ! Create local to local scatter.  Essentially remap the global to local as
  ! PETSc does in daltol.c
  call VecScatterCopy(ugdm%scatter_gtol,ugdm%scatter_ltol,ierr)
  call ISGetIndicesF90(ugdm%is_local_local,int_ptr,ierr)
  call VecScatterRemap(ugdm%scatter_ltol,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(ugdm%is_local_local,int_ptr,ierr)

#if UGRID_DEBUG
  string = 'scatter_ltol' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_ltol,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! Set up global to natural scatter
  ! Create index set of local non-ghosted Petsc ordering
  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax, &
                    PETSC_DETERMINE,vec_tmp,ierr)
!geh  call VecGetOwnershipRange(vec_tmp,istart,iend,ierr)
  call VecDestroy(vec_tmp,ierr)
  allocate(int_array(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax 
!geh    int_array(local_id) = (local_id-1)+istart
    int_array(local_id) = (local_id-1) + unstructured_grid%global_offset
  enddo
  call ISCreateGeneral(option%mycomm,unstructured_grid%nlmax, &
                       int_array,PETSC_COPY_VALUES,is_tmp,ierr) 
  deallocate(int_array)
  call AOPetscToApplicationIS(unstructured_grid%ao_natural_to_petsc, &
                              is_tmp,ierr)
  ! remap for ndof > 1  !geh: no longer need to accommodate ndof > 1, but leave
  ! alone for now.
  allocate(int_array(unstructured_grid%nlmax))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr)
  do local_id = 1, unstructured_grid%nlmax
    int_array(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr)
  call ISDestroy(is_tmp,ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%nlmax, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_natural,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  string = 'is_local_natural' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
  !                  PETSC_DETERMINE,vec_tmp,ierr)
  call VecCreate(option%mycomm,vec_tmp,ierr)
  call VecSetSizes(vec_tmp,unstructured_grid%nlmax*ndof,PETSC_DECIDE,ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr)
  call VecSetFromOptions(vec_tmp,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_petsc,vec_tmp, &
                        ugdm%is_local_natural,ugdm%scatter_gton,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_natural,vec_tmp, &
                        ugdm%is_local_petsc,ugdm%scatter_ntog,ierr)
  call VecDestroy(vec_tmp,ierr)

#if UGRID_DEBUG
  string = 'scatter_gton' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_gton,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  string = 'scatter_ntog' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_ntog,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine UGridCreateUGDM

! ************************************************************************** !
!
! UGridDMCreateJacobian: Creates a Jacobian matrix based on the unstructured
!                     grid dual
! author: Glenn Hammond
! date: 11/05/09
!
! ************************************************************************** !
subroutine UGridDMCreateJacobian(unstructured_grid,ugdm,mat_type,J,option)

  use Option_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  MatType :: mat_type
  Mat :: J
  type(option_type) :: option

  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt :: local_id, ineighbor, neighbor_id
  PetscInt :: iconn, id_up, id_dn
  PetscInt :: ndof_local
  PetscErrorCode :: ierr
  
  allocate(d_nnz(unstructured_grid%nlmax))
  allocate(o_nnz(unstructured_grid%nlmax))
  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0
  if (associated(unstructured_grid%explicit_grid)) then
    do iconn = 1, size(unstructured_grid%explicit_grid%connections,2)
      id_up = unstructured_grid%explicit_grid%connections(1,iconn)
      id_dn = unstructured_grid%explicit_grid%connections(2,iconn)
      if (id_up <= unstructured_grid%nlmax) then ! local
        if (id_dn <= unstructured_grid%nlmax) then
          d_nnz(id_up) = d_nnz(id_up) + 1
        else
          o_nnz(id_up) = o_nnz(id_up) + 1
        endif
      endif
      if (id_dn <= unstructured_grid%nlmax) then ! local
        if (id_up <= unstructured_grid%nlmax) then
          d_nnz(id_dn) = d_nnz(id_dn) + 1
        else
          o_nnz(id_dn) = o_nnz(id_dn) + 1
        endif
      endif
    enddo
  else
    do local_id = 1, unstructured_grid%nlmax
      do ineighbor = 1, unstructured_grid% &
                          cell_neighbors_local_ghosted(0,local_id)
        neighbor_id = unstructured_grid% &
                        cell_neighbors_local_ghosted(ineighbor,local_id)
        if (neighbor_id > 0) then
          d_nnz(local_id) = d_nnz(local_id) + 1
        else
          o_nnz(local_id) = o_nnz(local_id) + 1
        endif
      enddo
    enddo
  endif

  ndof_local = unstructured_grid%nlmax*ugdm%ndof
!  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*ugdm%ndof
        o_nnz = o_nnz*ugdm%ndof
        call MatCreateAIJ(option%mycomm,ndof_local,ndof_local, &
                          PETSC_DETERMINE,PETSC_DETERMINE, &
                          PETSC_NULL_INTEGER,d_nnz, &
                          PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog, &
                                        ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMappingBlock(J,ugdm%mapping_ltogb, &
                                             ugdm%mapping_ltogb,ierr)
      case(MATBAIJ)
        call MatCreateBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
                           PETSC_DETERMINE,PETSC_DETERMINE, &
                           PETSC_NULL_INTEGER,d_nnz, &
                           PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog, &
                                        ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMappingBlock(J,ugdm%mapping_ltogb, &
                                             ugdm%mapping_ltogb,ierr)
      case default
        option%io_buffer = 'MatType not recognized in UGridDMCreateJacobian'
        call printErrMsg(option)
    end select
!  else
!    select case(mat_type)
!      case(MATAIJ)
!        d_nnz = d_nnz*ugdm%ndof
!        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
!                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
!      case(MATBAIJ)
!        call MatCreateSeqBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
!                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
!      case default
!        option%io_buffer = 'MatType not recognized in UGridDMCreateJacobian'
!        call printErrMsg(option)
!    end select
!  endif

  deallocate(d_nnz)
  deallocate(o_nnz)
  
end subroutine UGridDMCreateJacobian

! ************************************************************************** !
!
! UGridDMCreateVector: Creates a global vector with PETSc ordering
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
subroutine UGridDMCreateVector(unstructured_grid,ugdm,vec,vec_type,option)

  use Option_module

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  Vec :: vec
  PetscInt :: vec_type
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  select case(vec_type)
    case(GLOBAL)
      !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
      !                  ugdm%ndof, &
      !                  PETSC_DETERMINE,vec,ierr)
      call VecCreate(option%mycomm,vec,ierr)
      call VecSetSizes(vec,unstructured_grid%nlmax*ugdm%ndof, &
                       PETSC_DECIDE,ierr)  
      call VecSetLocalToGlobalMapping(vec,ugdm%mapping_ltog,ierr)
      call VecSetLocalToGlobalMappingBlock(vec,ugdm%mapping_ltogb,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
      call VecSetFromOptions(vec,ierr)
    case(LOCAL)
      !call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax* &
      !                  ugdm%ndof, &
      !                  vec,ierr)
      call VecCreate(PETSC_COMM_SELF,vec,ierr)
      call VecSetSizes(vec,unstructured_grid%ngmax*ugdm%ndof, &
                  PETSC_DECIDE,ierr)  
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
      call VecSetFromOptions(vec,ierr)
    case(NATURAL)
      !call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
      !                  ugdm%ndof, &
      !                  PETSC_DETERMINE,vec,ierr)
      call VecCreate(option%mycomm,vec,ierr)
      call VecSetSizes(vec,unstructured_grid%nlmax*ugdm%ndof, &
                       PETSC_DECIDE,ierr)  
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
      call VecSetFromOptions(vec,ierr)
  end select
    
end subroutine UGridDMCreateVector

! ************************************************************************** !
!
! UGridMapIndices: maps global, local and natural indices of cells to each other
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
subroutine UGridMapIndices(unstructured_grid,ugdm,nG2L,nL2G,nG2A,nG2P,option)

  use Option_module

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nG2A(:)
  PetscInt, pointer :: nG2P(:)
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id

  allocate(nG2L(unstructured_grid%ngmax))
  allocate(nL2G(unstructured_grid%nlmax))
  allocate(nG2A(unstructured_grid%ngmax))
  
  ! initialize ghosted to 0
  !geh: any index beyond %nlmax will be 0 indicating that there is no local
  !     counterpart (i.e., it is a ghost cell)
  nG2L = 0

  !geh: Yes, it seems redundant that that we are setting both nL2G and nG2L to 
  !     the same index, but keep in mind that nG2L extends beyond %nlmax and
  !     we need these arrays to provide seemless integration for structured and
  !     unstructured
  do local_id = 1, unstructured_grid%nlmax
    nL2G(local_id) = local_id
    nG2L(local_id) = local_id
  enddo

  call ISGetIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    nG2A(ghosted_id) = int_ptr(ghosted_id)+1
  enddo
  call ISRestoreIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  nG2A = nG2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%ngmax, &
                            nG2A,ierr)
  nG2A = nG2A + 1 ! 1-based

#if MFD_UGRID
  allocate(nG2P(unstructured_grid%ngmax))
  do local_id = 1,unstructured_grid%nlmax
    nG2P(local_id) = local_id-1+unstructured_grid%global_offset
  enddo

  do ghosted_id = unstructured_grid%nlmax+1,unstructured_grid%ngmax
    nG2P(ghosted_id) = &
      unstructured_grid%ghost_cell_ids_petsc(ghosted_id-unstructured_grid%nlmax)-1
  enddo
#endif

end subroutine UGridMapIndices

! ************************************************************************** !
!
! UGridGet_Dual_Part_IS: Given an adjacency matrix, calculates the dual
!                        partitions, and provides a new IS with the ids 
!                        of the local cells on the processor
! author: Glenn Hammond
! date: 10/05/12
!
! ************************************************************************** !
subroutine UGridPartition(ugrid,option,Dual_mat,is_new, &
                          num_cells_local_new)

  use Option_module
  
  implicit none

#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"
  
  type(unstructured_grid_type) :: ugrid
  type(option_type) :: option
  Mat :: Dual_mat
  IS :: is_new
  PetscInt :: num_cells_local_new

  MatPartitioning :: Part
  PetscInt, allocatable :: cell_counts(:)
  PetscInt :: iflag
  PetscViewer :: viewer
  PetscInt :: local_vertex_offset
  PetscErrorCode :: ierr

#if UGRID_DEBUG
  call printMsg(option,'Partitioning')
#endif

  ! create the partitioning
  call MatPartitioningCreate(option%mycomm,Part,ierr)
  ! MatPartitioningSetAdjacency sets the adjacency graph (matrix) of the 
  ! thing to be partitioned.  - petsc
  call MatPartitioningSetAdjacency(Part,Dual_mat,ierr)
  call MatPartitioningSetFromOptions(Part,ierr)
  ! MatPartitioningApply gets a partitioning for a matrix. For each local cell 
  ! this tells the processor number that that cell is assigned to. - petsc
  ! is_new holds this information
  call MatPartitioningApply(Part,is_new,ierr)
  call MatPartitioningDestroy(Part,ierr)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_surf.out',viewer,ierr)
  endif
  call ISView(is_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  
#endif

  ! calculate the number of local grid cells on each processor
  allocate(cell_counts(option%mycommsize))
  ! ISPartitioningCount takes a ISPartitioning and determines the number of  
  ! resulting elements on each (partition) process - petsc
  call ISPartitioningCount(is_new,option%mycommsize,cell_counts,ierr)
  num_cells_local_new = cell_counts(option%myrank+1) 
  call MPI_Allreduce(num_cells_local_new,iflag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MIN,option%mycomm,ierr)
  deallocate(cell_counts)
  if (iflag < 1) then
    option%io_buffer = 'A processor core has been assigned zero cells.'
    call printErrMsg(option)
  endif
  
end subroutine UGridPartition  

! ************************************************************************** !
!
! UGridNaturalToPetsc: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGridCreateOldVec(ugrid,option,elements_old, &
                             num_cells_local_old, &
                             is_new,is_scatter,stride)
  use Option_module                  

  implicit none

#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(unstructured_grid_type) :: ugrid
  type(option_type) :: option
  Vec :: elements_old
  PetscInt :: num_cells_local_old
  IS :: is_new 
  IS :: is_scatter
  PetscInt :: stride

  PetscViewer :: viewer
  IS :: is_num  
  PetscInt, pointer :: index_ptr(:)  
  PetscErrorCode :: ierr  
  
  ! calculate the global offsets in the new vector for each grid cell
  
  ! ISPartitioningToNumbering takes an ISPartitioning and on each processor 
  ! generates an IS that contains a new global node number for each index 
  ! based on the partitioning. - petsc
  call ISPartitioningToNumbering(is_new,is_num,ierr)
  call ISDestroy(is_new,ierr)
  call ISGetIndicesF90(is_num,index_ptr,ierr)

  ! Create a mapping of local indices to global strided (use block ids, not
  ! indices)
  call ISCreateBlock(option%mycomm,stride, &
                     num_cells_local_old, &
                     index_ptr,PETSC_COPY_VALUES,is_scatter,ierr)
  call ISRestoreIndicesF90(is_num,index_ptr,ierr)
  call ISDestroy(is_num,ierr)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_old_to_new_subsurf.out', &
                              viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_old_to_new_surf.out', &
                              viewer,ierr)
  endif
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! create another strided vector with the old cell/element distribution
  call VecCreate(option%mycomm,elements_old,ierr)
  call VecSetSizes(elements_old,stride*num_cells_local_old,PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_old,ierr)  

end subroutine UGridCreateOldVec

! ************************************************************************** !
!
! UGridNaturalToPetsc: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGridNaturalToPetsc(ugrid,option,elements_old,elements_local, &
                               num_cells_local_new,stride,dual_offset, &
                               natural_id_offset,is_scatter)

  use Option_module
  use Utility_module, only: reallocateIntArray, DeallocateArray
  
  implicit none

#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  type(unstructured_grid_type) :: ugrid
  type(option_type) :: option
  Vec :: elements_old, elements_local
  PetscInt :: num_cells_local_new
  PetscInt :: stride
  PetscInt :: dual_offset
  PetscInt :: natural_id_offset
  IS :: is_scatter  
  
  Vec :: elements_petsc, elements_natural
  PetscViewer :: viewer
  VecScatter :: vec_scatter
  IS :: is_gather  
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: global_offset_new
  PetscInt :: local_id, idual, count
  PetscInt :: max_ghost_cell_count
  PetscInt :: temp_int
  PetscInt :: ghosted_id
  PetscInt :: ghost_cell_count
  PetscInt :: ghost_offset_new
  PetscInt :: dual_id
  PetscBool :: found
  PetscReal, pointer :: vec_ptr(:)  
  PetscReal, pointer :: vec_ptr2(:)  
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: int_array5(:)  
  PetscInt, pointer :: int_array_pointer(:)  
  PetscErrorCode :: ierr

  ! create a petsc vec to store all the information for each element
  ! based on the stride calculated above.  
  call VecCreate(option%mycomm,elements_natural,ierr)
  call VecSetSizes(elements_natural, &
                   stride*num_cells_local_new, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_natural,ierr)  
  
#if UGRID_DEBUG
  call printMsg(option,'Before element scatter')
#endif

  ! scatter all the cell data from the old decomposition (as read in in 
  ! parallel) to the more parmetis-calculated decomposition
  call VecScatterCreate(elements_old,PETSC_NULL_OBJECT,elements_natural,is_scatter, &
                        vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_natural, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_natural, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

#if UGRID_DEBUG
  call printMsg(option,'After element scatter')
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'elements_old_suburf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'elements_old_surf.out',viewer,ierr)
  endif
  call VecView(elements_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  call VecDestroy(elements_old,ierr)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'elements_natural_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'elements_natural_surf.out',viewer,ierr)
  endif
  call VecView(elements_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! update global offset based on new partitioning
  global_offset_new = 0
  call MPI_Exscan(num_cells_local_new,global_offset_new, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  ugrid%global_offset = global_offset_new

  allocate(ugrid%cell_ids_natural(num_cells_local_new))
  ugrid%cell_ids_natural = 0
  
  ! look at all connections and determine how many are non-local, and create
  !  a listof indices
  call VecDuplicate(elements_natural,elements_petsc,ierr)
  call VecCopy(elements_natural,elements_petsc,ierr)

#if UGRID_DEBUG
  call printMsg(option,'Lists of ids')
#endif

  ! now we unpack the decomposed cell data
  
  ! store the natural grid cell id for each local cell as read from the grid 
  ! file
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  do local_id = 1, num_cells_local_new
    ! Cell id read from explicit grid input file is second entry in block
    ! It may match the first entry (the calculated natural id based on the
    ! order that cells were read), but it need not.
    ugrid%cell_ids_natural(local_id) = &
      int(abs(vec_ptr((local_id-1)*stride+natural_id_offset)))
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'natural_ids' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id = 1, num_cells_local_new
    write(86,'(i5)') ugrid%cell_ids_natural(local_id)
  enddo
  close(86)
#endif     

  ! make a list of petsc ids for each local cell (you simply take the global 
  ! offset and add it to the local contiguous cell ids on each processor
  allocate(int_array(num_cells_local_new))
  do local_id = 1, num_cells_local_new
    int_array(local_id) = local_id+global_offset_new
  enddo
  
  ! make the arrays zero-based
  int_array = int_array - 1
  ugrid%cell_ids_natural = ugrid%cell_ids_natural - 1
  ! create an application ordering (mapping of natural to petsc ordering)
  call AOCreateBasic(option%mycomm,num_cells_local_new, &
                     ugrid%cell_ids_natural,int_array, &
                     ugrid%ao_natural_to_petsc,ierr)
  deallocate(int_array)
  ! make cell_ids_natural 1-based again
  ugrid%cell_ids_natural = ugrid%cell_ids_natural + 1

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'ao_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'ao_surf.out',viewer,ierr)
  endif
  call AOView(ugrid%ao_natural_to_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! The below creates a list of cells ids for the duals and converts them
  ! to petsc ordering   
  
  ! count the number of cells and their duals  
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    do idual = 1, ugrid%max_ndual_per_cell
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      if (dual_id < 1) exit ! here we hit the 0 at the end of last dual
      count = count + 1
    enddo
  enddo     
               
  ! allocate and fill an array with the natural cell and dual ids
  allocate(int_array(count))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    int_array(count) = ugrid%cell_ids_natural(local_id)
    do idual = 1, ugrid%max_ndual_per_cell
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      if (dual_id < 1) exit ! again we hit the 0 
      count = count + 1
      int_array(count) = dual_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

#if UGRID_DEBUG
  call printMsg(option,'Application ordering')
#endif

  ! convert the dual ids in int_array from natural to petsc numbering
  int_array = int_array - 1             
  call AOApplicationToPetsc(ugrid%ao_natural_to_petsc,count, &
                            int_array,ierr)
  int_array = int_array + 1                

#if UGRID_DEBUG
  call printMsg(option,'PETSc-ordered duals')
#endif

  ! load mapped petsc-ordered dual ids back into duplicated vector
  ! exactly the opposite operation of when we loaded the temporary int_array
  ! vector
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
!geh: do not believe that we need elements_natural here
!  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  allocate(ugrid%cell_ids_petsc(num_cells_local_new))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    ! extract the petsc id for the cell
    ugrid%cell_ids_petsc(local_id) = int_array(count)
    ! store it in the elements_petsc vector too
    vec_ptr((local_id-1)*stride+1) = int_array(count)
    do idual = 1, ugrid%max_ndual_per_cell
!geh      dual_id = vec_ptr2(idual + dual_offset + (local_id-1)*stride)
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      if (dual_id < 1) exit
      count = count + 1
      ! store the petsc numbered duals in the vector also
      vec_ptr(idual + dual_offset + (local_id-1)*stride) = int_array(count)
    enddo
  enddo                
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
!geh  call VecRestoreArrayF90(elements_natural,vec_ptr2,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_surf.out',viewer,ierr)
  endif
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! make a list of ghosted ids in petsc numbering
#if UGRID_DEBUG
  call printMsg(option,'Renumbering ghost ids to petsc numbering')
#endif
  
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ghost_cell_count = 0
  ! allocate a temporarily-sized array
  ! geh: this assumes that the number of ghost cells will not exceed the number
  !      of local and 100 is used to ensure that if this is not true, the array
  !       is still large enough
  max_ghost_cell_count = max(num_cells_local_new,100)
  allocate(int_array_pointer(max_ghost_cell_count))
  int_array_pointer = 0
  ! loop over all duals and find the off-processor cells on the other
  ! end of a dual
  do local_id=1, num_cells_local_new
    do idual = 1, ugrid%max_ndual_per_cell
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      found = PETSC_FALSE
      if (dual_id < 1) exit
      if (dual_id <= global_offset_new .or. &
          dual_id > global_offset_new + num_cells_local_new) then
        ghost_cell_count = ghost_cell_count + 1
        ! reallocate the ghost cell array if necessary
        if (ghost_cell_count > max_ghost_cell_count) then
          call reallocateIntArray(int_array_pointer,max_ghost_cell_count)
        endif
        int_array_pointer(ghost_cell_count) = dual_id
        vec_ptr(idual + dual_offset + (local_id-1)*stride) = &
          -ghost_cell_count
        ! temporarily store the index of the int_array_pointer
        ! flag negative
      else
        ! if dual_id is in petsc numbering, then the local ids is:
        vec_ptr(idual + dual_offset + (local_id-1)*stride) = &
          dual_id - global_offset_new
      endif
    enddo
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm, &
                              'elements_local_dual_unsorted_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm, &
                              'elements_local_dual_unsorted_surf.out',viewer,ierr)
  endif
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
 

#if UGRID_DEBUG
  call printMsg(option,'  Sorting local ghost ids')
#endif

  if (ghost_cell_count > 0) then
    ! sort ghost cell ids
    allocate(int_array2(ghost_cell_count))
    do ghosted_id = 1, ghost_cell_count
      int_array2(ghosted_id) = ghosted_id
    enddo

    ! sort with permutation
    ! int_array_pointer = ghost ids unsorted before and after
    ! int_array2 = permutation
    ! convert to 0-based
    int_array_pointer = int_array_pointer-1
    int_array2 = int_array2 - 1
    call PetscSortIntWithPermutation(ghost_cell_count,int_array_pointer, &
                                     int_array2,ierr)
    ! convert to 1-based
    int_array_pointer = int_array_pointer+1
    int_array2 = int_array2 + 1

    ! determine how many duplicates
    allocate(int_array3(ghost_cell_count))
    allocate(int_array4(ghost_cell_count))
    allocate(int_array5(ghost_cell_count))
    int_array3 = 0
    temp_int = 1
    int_array3(temp_int) = int_array_pointer(int_array2(1))
    ! do not change ghosted_id = 1 to ghosted_id = 2 as the first value in
    ! int_array2() will not be set correctly.
    do ghosted_id = 1, ghost_cell_count
      if (int_array3(temp_int) < &
            int_array_pointer(int_array2(ghosted_id))) then
        temp_int = temp_int + 1
        int_array3(temp_int) = int_array_pointer(int_array2(ghosted_id))
      endif
      int_array5(int_array2(ghosted_id)) = ghosted_id
      int_array4(ghosted_id) = temp_int
    enddo

    ghost_cell_count = temp_int
    allocate(ugrid%ghost_cell_ids_petsc(ghost_cell_count))
    ugrid%ghost_cell_ids_petsc = 0

    ugrid%ghost_cell_ids_petsc(1:ghost_cell_count) = &
      int_array3(1:ghost_cell_count)

#if UGRID_DEBUG
  call printMsg(option,'  Remappping ghost ids')
#endif

    ! remap of duals of ghost cells
    call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
    do local_id=1, num_cells_local_new
      do idual = 1, ugrid%max_ndual_per_cell
        ! dual_id is now the negative of the local unsorted ghost cell id
        dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
        ! dual_id = 0: not assigned
        ! dual_id > 0: assigned to local cell
        ! dual_id < 0: assigned to ghost cell
        if (dual_id < 0) then
          vec_ptr(idual + dual_offset + (local_id-1)*stride) = &
            int_array4(int_array5(-dual_id)) + num_cells_local_new
        endif
      enddo
    enddo
    call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

    deallocate(int_array2)
    deallocate(int_array3)
    deallocate(int_array4)
    deallocate(int_array5)
  endif
  call DeallocateArray(int_array_pointer)

  ugrid%nlmax = num_cells_local_new
  ugrid%num_ghost_cells = ghost_cell_count
  ugrid%ngmax = &
    num_cells_local_new + ghost_cell_count


#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'elements_local_dual_subsurf.out', &
                              viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'elements_local_dual_surf.out', &
                              viewer,ierr)
  endif
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  call printMsg(option,'Resizing natural cell id array')
#endif

  ! Resize cell_ids_natural to include ghosted cells
  allocate(int_array(ugrid%nlmax))
  int_array(:) = ugrid%cell_ids_natural(:)
  deallocate(ugrid%cell_ids_natural)
  allocate(ugrid%cell_ids_natural(ugrid%ngmax))
  ugrid%cell_ids_natural(:) = -999
  ugrid%cell_ids_natural(1:ugrid%nlmax) = int_array(:)
  deallocate(int_array)
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  do local_id=1, ugrid%nlmax
    do idual = 1, ugrid%max_ndual_per_cell
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      if (dual_id < 1) exit
      if (dual_id > ugrid%nlmax) then
        ugrid%cell_ids_natural(dual_id) = &
          int(vec_ptr2(idual + dual_offset + (local_id-1)*stride))
      endif       
    enddo
  enddo
  if (minval(ugrid%cell_ids_natural) < 1) then
    write(string,*) minval( ugrid%cell_ids_natural)
    option%io_buffer = 'Negative natural id: ' // trim(adjustl(string))
    call printErrMsgByRank(option)
  endif
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  call VecRestoreArrayF90(elements_natural,vec_ptr2,ierr)
  call VecDestroy(elements_natural,ierr)  
  
  ! NOW START ON GHOSTING CELLS

#if UGRID_DEBUG
  call printMsg(option,'Ghosting local element vectors')
#endif

  ! load cell neighbors into array
  ! start first index at zero to store # duals for a cell
  allocate(ugrid%cell_neighbors_local_ghosted( &
           0:ugrid%max_ndual_per_cell,ugrid%nlmax))
  ugrid%cell_neighbors_local_ghosted = 0
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do local_id=1, ugrid%nlmax
    count = 0
    do idual = 1, ugrid%max_ndual_per_cell
      dual_id = int(vec_ptr(idual + dual_offset + (local_id-1)*stride))
      if (dual_id < 1) exit
      count = count + 1
      ! flag ghosted cells in dual as negative
      !geh: these negative dual ids are used later in UGridDMCreateJacobian() 
      !     to specify off processor connectity in the Jacobian
      if (dual_id > ugrid%nlmax) dual_id = -dual_id
      ugrid%cell_neighbors_local_ghosted(idual,local_id) = dual_id
    enddo
    ! set the # of duals in for the cell
    ugrid%cell_neighbors_local_ghosted(0,local_id) = count
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! need to create a local ghosted vector in which we can collect element info
  ! including ghost cells
  !call VecCreateSeq(PETSC_COMM_SELF,stride*ugrid%ngmax, &
  !                  elements_local,ierr)
  call VecCreate(PETSC_COMM_SELF,elements_local,ierr)
  call VecSetSizes(elements_local,stride*ugrid%ngmax,PETSC_DECIDE,ierr)
  call VecSetBlockSize(elements_local,stride,ierr)
  call VecSetFromOptions(elements_local,ierr)
  allocate(int_array(ugrid%ngmax))
  int_array(1:ugrid%nlmax) = &
    ugrid%cell_ids_petsc(:)
  if (ugrid%num_ghost_cells > 0) then
    int_array(ugrid%nlmax+1:ugrid%ngmax) = &
      ugrid%ghost_cell_ids_petsc(:)
  endif
  int_array = int_array-1
  call ISCreateBlock(option%mycomm,stride,ugrid%ngmax, &
                     int_array,PETSC_COPY_VALUES,is_scatter,ierr)
  do ghosted_id = 1, ugrid%ngmax
    int_array(ghosted_id) = ghosted_id-1
  enddo
  call ISCreateBlock(option%mycomm,stride,ugrid%ngmax, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_local_to_ghost_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_local_to_ghost_surf.out',viewer,ierr)
  endif
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  if (ugrid%grid_type == THREE_DIM_GRID) then
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_elem_local_to_ghost_subsurf.out',viewer,ierr)
  else
    call PetscViewerASCIIOpen(option%mycomm,'is_gather_elem_local_to_ghost_surf.out',viewer,ierr)
  endif
  call ISView(is_gather,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecScatterCreate(elements_petsc,is_scatter,elements_local,is_gather, &
                        vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call ISDestroy(is_gather,ierr)
  call VecScatterBegin(vec_scatter,elements_petsc,elements_local, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_petsc,elements_local, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)
  call VecDestroy(elements_petsc,ierr)
    
#if UGRID_DEBUG
  write(string,*) option%myrank
  if (ugrid%grid_type == THREE_DIM_GRID) then
    string = 'elements_local' // trim(adjustl(string)) // '_subsurf.out'
  else
    string = 'elements_local' // trim(adjustl(string)) // '_surf.out'
  endif
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call VecView(elements_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

  call printMsg(option,'Scatter/gathering local ghosted vertices')
#endif

end subroutine UGridNaturalToPetsc  
  
! ************************************************************************** !
!
! UGridDestroy: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGridDestroy(unstructured_grid)

  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(unstructured_grid_type), pointer :: unstructured_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(unstructured_grid)) return

  ! variables for all unstructured grids
  call DeallocateArray(unstructured_grid%hash)
  call DeallocateArray(unstructured_grid%cell_ids_natural)
  call DeallocateArray(unstructured_grid%cell_ids_petsc)
  call DeallocateArray(unstructured_grid%ghost_cell_ids_petsc)
  call UGridExplicitDestroy(unstructured_grid%explicit_grid)
  call UGridPolyhedraDestroy(unstructured_grid%polyhedra_grid)
  if (unstructured_grid%ao_natural_to_petsc /= 0) &
    call AODestroy(unstructured_grid%ao_natural_to_petsc,ierr)
  
  ! variables for implicit unstructured grids
  call DeallocateArray(unstructured_grid%cell_type)
  call DeallocateArray(unstructured_grid%cell_vertices)
  call DeallocateArray(unstructured_grid%face_to_cell_ghosted)
  call DeallocateArray(unstructured_grid%connection_to_face)
  call DeallocateArray(unstructured_grid%face_to_vertex_natural)
  call DeallocateArray(unstructured_grid%face_to_vertex)
  call DeallocateArray(unstructured_grid%cell_to_face_ghosted)
  call DeallocateArray(unstructured_grid%vertex_ids_natural)
  call DeallocateArray(unstructured_grid%cell_neighbors_local_ghosted)
  if (associated(unstructured_grid%vertices)) &
    deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)  
  if (associated(unstructured_grid%face_centroid)) &
    deallocate(unstructured_grid%face_centroid)
  nullify(unstructured_grid%face_centroid)  
  call DeallocateArray(unstructured_grid%face_area)
  
  deallocate(unstructured_grid)
  nullify(unstructured_grid)

end subroutine UGridDestroy

! ************************************************************************** !
!
! UGridDMDestroy: Deallocates a unstructured grid distributed mesh
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGridDMDestroy(ugdm)

  implicit none

  type(ugdm_type), pointer :: ugdm
  
  PetscErrorCode :: ierr
    
  if (.not.associated(ugdm)) return
  
  call ISDestroy(ugdm%is_ghosted_local,ierr)
  call ISDestroy(ugdm%is_local_local,ierr)
  call ISDestroy(ugdm%is_ghosted_petsc,ierr)
  call ISDestroy(ugdm%is_local_petsc,ierr)
  call ISDestroy(ugdm%is_ghosts_local,ierr)
  call ISDestroy(ugdm%is_ghosts_petsc,ierr)
  call ISDestroy(ugdm%is_local_natural,ierr)
  call VecScatterDestroy(ugdm%scatter_ltog,ierr)
  call VecScatterDestroy(ugdm%scatter_gtol,ierr)
  call VecScatterDestroy(ugdm%scatter_ltol,ierr)
  call VecScatterDestroy(ugdm%scatter_gton,ierr)
  call VecScatterDestroy(ugdm%scatter_ntog,ierr)
  call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltog,ierr)
  if (ugdm%mapping_ltogb /= 0) &
    call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltogb,ierr)
  call VecDestroy(ugdm%global_vec,ierr)
  call VecDestroy(ugdm%local_vec,ierr)
#ifdef SURFACE_FLOW
  call VecScatterDestroy(ugdm%scatter_bet_grids,ierr)
  call VecScatterDestroy(ugdm%scatter_bet_grids_1dof,ierr)
  call VecScatterDestroy(ugdm%scatter_bet_grids_ndof,ierr)
#endif
  deallocate(ugdm)
  nullify(ugdm)

end subroutine UGridDMDestroy

! ************************************************************************** !
!
! UGridExplicitDestroy: Deallocates an explicit unstructured grid object
! author: Glenn Hammond
! date: 05/14/12
!
! ************************************************************************** !
subroutine UGridExplicitDestroy(explicit_grid)

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(unstructured_explicit_type), pointer :: explicit_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(explicit_grid)) return

  call DeallocateArray(explicit_grid%cell_ids)
  call DeallocateArray(explicit_grid%cell_volumes)
  if (associated(explicit_grid%cell_centroids)) &
    deallocate(explicit_grid%cell_centroids)
  nullify(explicit_grid%cell_centroids)
  call DeallocateArray(explicit_grid%connections)
  call DeallocateArray(explicit_grid%face_areas)
  if (associated(explicit_grid%face_centroids)) &
    deallocate(explicit_grid%face_centroids)
  nullify(explicit_grid%face_centroids)
  
  deallocate(explicit_grid)
  nullify(explicit_grid)

end subroutine UGridExplicitDestroy

! ************************************************************************** !
!> This routine deallocates a polyhedra unstructured grid object
!!
!> @author
!! Gautam Bisht, LBL
!!
!! date: 09/29/13
! ************************************************************************** !
subroutine UGridPolyhedraDestroy(polyhedra_grid)

  use Utility_module, only : DeallocateArray

  implicit none

  type(unstructured_polyhedra_type), pointer :: polyhedra_grid

  PetscErrorCode :: ierr

  if (.not.associated(polyhedra_grid)) return

  call DeallocateArray(polyhedra_grid%cell_ids)
  call DeallocateArray(polyhedra_grid%cell_nfaces)
  call DeallocateArray(polyhedra_grid%cell_nverts)
  call DeallocateArray(polyhedra_grid%cell_faceids)
  call DeallocateArray(polyhedra_grid%cell_vertids)
  call DeallocateArray(polyhedra_grid%cell_volumes)
  call DeallocateArray(polyhedra_grid%face_ids)
  call DeallocateArray(polyhedra_grid%face_cellids)
  call DeallocateArray(polyhedra_grid%face_nverts)
  call DeallocateArray(polyhedra_grid%face_vertids)
  call DeallocateArray(polyhedra_grid%face_areas)

  if (associated(polyhedra_grid%cell_centroids)) &
    deallocate(polyhedra_grid%cell_centroids)
  nullify(polyhedra_grid%cell_centroids)
  if (associated(polyhedra_grid%face_centroids)) &
    deallocate(polyhedra_grid%face_centroids)
  nullify(polyhedra_grid%face_centroids)
  if (associated(polyhedra_grid%vertex_coordinates)) &
    deallocate(polyhedra_grid%vertex_coordinates)
  nullify(polyhedra_grid%vertex_coordinates)

  deallocate(polyhedra_grid)
  nullify(polyhedra_grid)

end subroutine UGridPolyhedraDestroy

end module Unstructured_Grid_Aux_module
