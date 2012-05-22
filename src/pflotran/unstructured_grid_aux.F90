module Unstructured_Grid_Aux_module

!  use Connection_module
  use Unstructured_Cell_module
  use Unstructured_Explicit_module
  
  implicit none

  private 
  
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#if defined(PARALLELIO_LIB)
  include "piof.h"
#endif

  PetscInt, parameter, public :: TWO_DIM_GRID = 1
  PetscInt, parameter, public :: THREE_DIM_GRID = 2 
  
  !geh: for debugging purposes make these large
  !TODO(geh): change values to 1, 2 respectively.
  PetscInt, parameter, public :: IMPLICIT_UNSTRUCTURED_GRID = 98
  PetscInt, parameter, public :: EXPLICIT_UNSTRUCTURED_GRID = 99

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
    ! variables for implicit unstructured grids
    PetscInt :: grid_type         ! 3D subsurface (default) or 2D surface grid
    PetscInt :: num_vertices_global ! number of vertices in entire problem domain
    PetscInt :: num_vertices_local  ! number of vertices in local grid cells
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
            UGridMapIndices, &
            UGridDMCreateJacobian, &
            UGridDMCreateVector, &
            UGridDestroy, &
            UGridCreateUGDM, &
            UGridDMDestroy

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
      if (id_up > 0) then
        if (id_dn > 0) then
          d_nnz(id_up) = d_nnz(id_up) + 1
        else
          o_nnz(id_up) = o_nnz(id_up) + 1
        endif
      endif
      if (id_dn > 0) then
        if (id_up > 0) then
          d_nnz(id_dn) = d_nnz(id_dn) + 1
        else
          o_nnz(id_dn) = o_nnz(id_dn) + 1
        endif
      endif
    enddo
  else
    do local_id = 1, unstructured_grid%nlmax
      do ineighbor = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
        neighbor_id = unstructured_grid%cell_neighbors_local_ghosted(ineighbor,local_id)
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
#ifdef MATCREATE_OLD      
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
#else
        call MatCreateAIJ(option%mycomm,ndof_local,ndof_local, &
#endif
                          PETSC_DETERMINE,PETSC_DETERMINE, &
                          PETSC_NULL_INTEGER,d_nnz, &
                          PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog, &
                                        ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMappingBlock(J,ugdm%mapping_ltogb, &
                                             ugdm%mapping_ltogb,ierr)
      case(MATBAIJ)
#ifdef MATCREATE_OLD      
        call MatCreateMPIBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
#else
        call MatCreateBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
#endif
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
subroutine UGridMapIndices(unstructured_grid,ugdm,nG2L,nL2G,nG2A)

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nG2A(:)
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

end subroutine UGridMapIndices

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
  call ExplicitUGridDestroy(unstructured_grid%explicit_grid)
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
#endif
  deallocate(ugdm)
  nullify(ugdm)

end subroutine UGridDMDestroy

end module Unstructured_Grid_Aux_module
