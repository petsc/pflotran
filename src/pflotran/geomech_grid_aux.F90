#ifdef GEOMECH
module Geomech_Grid_Aux_module

  use Unstructured_Cell_module
 
  implicit none

  private 
  
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#if defined(SCORPIO)
  include "scorpiof.h"
#endif
 
  type, public :: geomech_grid_type
    ! variables for geomechanics grid (unstructured) for finite element formulation
    ! The dofs (displacements, here) are solved at the nodes 
    ! Notation corresponding to finite volume: node implies cell vertex and element means cell
    PetscInt :: global_offset                ! offset in petsc ordering for the first cell on a processor
    PetscInt :: nmax_elem                    ! Total number of elements in the global domain
    PetscInt :: nlmax_elem                   ! Total number of non-ghosted elements on a processor
    PetscInt :: nmax_node                    ! Total number of nodes in the global domain 
    PetscInt :: nlmax_node                   ! Total number of non-ghosted nodes on a processor
    PetscInt :: ngmax_node                   ! Total number of ghosted nodes on a processor
    PetscInt, pointer :: elem_ids_natural(:) ! Natural numbering of elements on a processor
    PetscInt, pointer :: elem_ids_petsc(:)   ! Petsc numbering of elements on a processor
    AO :: ao_natural_to_petsc                ! mapping of natural to Petsc ordering
    PetscInt :: max_ndual_per_elem           ! Max. number of dual elements connected to an element
    PetscInt :: max_nnode_per_elem           ! Max. number of nodes per element
    PetscInt :: max_elem_sharing_a_node      
    PetscInt, pointer :: elem_type(:)        ! Type of element
    PetscInt, pointer :: elem_nodes(:,:)     ! Node number on each element
    type(point_type), pointer :: nodes(:)    ! Coordinates of the nodes
    PetscInt, pointer :: node_ids_natural(:) ! Natural ids of the nodes
    PetscInt, pointer :: ghosted_node_ids_natural(:) ! Natural ids of the ghosted nodes
  end type geomech_grid_type
  

  type, public :: gmdm_type                  ! Geomech. DM type
    ! local: included both local (non-ghosted) and ghosted nodes
    ! global: includes only local (non-ghosted) nodes
    PetscInt :: ndof
    ! for the below
    ! ghosted = local (non-ghosted) and ghosted nodes
    ! local = local (non-ghosted) nodes
    IS :: is_ghosted_local                   ! IS for ghosted nodes with local on-processor numbering
    IS :: is_local_local                     ! IS for local nodes with local on-processor numbering
    IS :: is_ghosted_petsc                   ! IS for ghosted nodes with petsc numbering
    IS :: is_local_petsc                     ! IS for local nodes with petsc numbering
    IS :: is_ghosts_local                    ! IS for ghost nodes with local on-processor numbering
    IS :: is_ghosts_petsc                    ! IS for ghost nodes with petsc numbering
    IS :: is_local_natural                   ! IS for local nodes with natural (global) numbering
    VecScatter :: scatter_ltog               ! scatter context for local to global updates
    VecScatter :: scatter_gtol               ! scatter context for global to local updates
    VecScatter :: scatter_ltol               ! scatter context for local to local updates
    VecScatter :: scatter_gton               ! scatter context for global to natural updates
    VecScatter :: scatter_ntog               ! scatter context for natural to global updates
    ISLocalToGlobalMapping :: mapping_ltog   ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb  ! block form of mapping_ltog
    Vec :: global_vec                        ! global vec (no ghost nodes), petsc-ordering
    Vec :: local_vec                         ! local vec (includes local and ghosted nodes), local ordering
  end type gmdm_type

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: GMGridCreate, &
            GMGridDestroy, &
            GMDMCreate, &
            GMDMDestroy, &
            GMCreateGMDM
            
  
contains

! ************************************************************************** !
!
! GMDMCreate: Creates a geomech grid distributed mesh object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
function GMDMCreate()

  implicit none
  
  type(gmdm_type), pointer :: GMDMCreate
  type(gmdm_type), pointer :: gmdm

  allocate(gmdm)
  gmdm%is_ghosted_local = 0
  gmdm%is_local_local = 0
  gmdm%is_ghosted_petsc = 0
  gmdm%is_local_petsc = 0
  gmdm%is_ghosts_local = 0
  gmdm%is_ghosts_petsc = 0
  gmdm%is_local_natural = 0
  gmdm%scatter_ltog = 0
  gmdm%scatter_gtol = 0
  gmdm%scatter_ltol  = 0
  gmdm%scatter_gton = 0
  gmdm%scatter_ntog = 0
  gmdm%mapping_ltog = 0
  gmdm%mapping_ltogb = 0
  gmdm%global_vec = 0
  gmdm%local_vec = 0

  GMDMCreate => gmdm

end function GMDMCreate

! ************************************************************************** !
!
! GMGridCreate: Creates a geomechanics grid object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
function GMGridCreate()

  implicit none
  
  type(geomech_grid_type), pointer :: GMGridCreate
  type(geomech_grid_type), pointer :: geomech_grid

  allocate(geomech_grid)

  ! variables for all unstructured grids
  
  geomech_grid%global_offset = 0
  geomech_grid%nmax_elem = 0
  geomech_grid%nlmax_elem = 0
  geomech_grid%nmax_node = 0
  geomech_grid%nlmax_node = 0
  nullify(geomech_grid%elem_ids_natural)
  nullify(geomech_grid%elem_ids_petsc)
  geomech_grid%ao_natural_to_petsc = 0
  geomech_grid%max_ndual_per_elem = 0
  geomech_grid%max_nnode_per_elem = 0
  geomech_grid%max_elem_sharing_a_node = 0
  nullify(geomech_grid%elem_type)
  nullify(geomech_grid%elem_nodes)
  nullify(geomech_grid%nodes)
  nullify(geomech_grid%node_ids_natural)
  nullify(geomech_grid%ghosted_node_ids_natural)

  GMGridCreate => geomech_grid
  
end function GMGridCreate

! ************************************************************************** !
!
! GMCreateGMDM: Mapping/scatter contexts are created for PETSc DM object
! author: Satish Karra, LANL
! date: 05/30/13
!
! ************************************************************************** !
subroutine GMCreateGMDM(geomech_grid,gmdm,ndof,option)

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

  type(geomech_grid_type)             :: geomech_grid
  type(option_type)                   :: option
  type(gmdm_type), pointer            :: gmdm
  PetscInt                            :: ndof
  PetscInt, pointer                   :: int_ptr(:)
  PetscInt                            :: local_id, ghosted_id
  PetscInt                            :: idof
  IS                                  :: is_tmp
  Vec                                 :: vec_tmp
  PetscErrorCode                      :: ierr
  character(len=MAXWORDLENGTH)        :: ndof_word
  character(len=MAXSTRINGLENGTH)      :: string
  
  PetscViewer :: viewer

  PetscInt, allocatable :: int_array(:)
  
  gmdm => GMDMCreate()
  gmdm%ndof = ndof
  
#if GEOMECH_DEBUG
  write(ndof_word,*) ndof
  ndof_word = adjustl(ndof_word)
  ndof_word = '_' // trim(ndof_word)
  string = 'Vectors' // ndof_word
  call printMsg(option,string)
#endif

  ! create global vec
  call VecCreate(option%mycomm,gmdm%global_vec,ierr)
  call VecSetSizes(gmdm%global_vec,geomech_grid%nlmax_node*ndof,PETSC_DECIDE,&
                    ierr)  
  call VecSetBlockSize(gmdm%global_vec,ndof,ierr)
  call VecSetFromOptions(gmdm%global_vec,ierr)

  ! create local vec
  call VecCreate(PETSC_COMM_SELF,gmdm%local_vec,ierr)
  call VecSetSizes(gmdm%local_vec,geomech_grid%ngmax_node*ndof,PETSC_DECIDE,ierr)
  call VecSetBlockSize(gmdm%local_vec,ndof,ierr)
  call VecSetFromOptions(gmdm%local_vec,ierr)
  
  ! IS for global numbering of local, non-ghosted vertices
  ! ISCreateBlock requires block ids, not indices.  Therefore, istart should be
  ! the offset of the block from the beginning of the vector.
  allocate(int_array(geomech_grid%nlmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset
  enddo

  ! arguments for ISCreateBlock():
  ! option%mycomm  - the MPI communicator
  ! ndof  - number of elements in each block
  ! geomech_grid%nlmax  - the length of the index set
  !                                      (the number of blocks
  ! int_array  - the list of integers, one for each block and count
  !              of block not indices
  ! PETSC_COPY_VALUES  - see PetscCopyMode, only PETSC_COPY_VALUES and
  !                      PETSC_OWN_POINTER are supported in this routine
  ! ugdm%is_local_petsc - the new index set
  ! ierr - PETScErrorCode
  call ISCreateBlock(option%mycomm,ndof,geomech_grid%nlmax_node, &
                     int_array,PETSC_COPY_VALUES,gmdm%is_local_petsc,ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  string = 'geomech_is_local_petsc' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(gmdm%is_local_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  



end subroutine GMCreateGMDM

! ************************************************************************** !
!
! GMDMDestroy: Deallocates a geomechanics grid distributed mesh object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
subroutine GMGridDestroy(geomech_grid)

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(geomech_grid_type), pointer :: geomech_grid
  
  PetscErrorCode :: ierr
  
  if (.not.associated(geomech_grid)) return
  
  call DeallocateArray(geomech_grid%elem_ids_natural)
  call DeallocateArray(geomech_grid%elem_ids_petsc)
  call DeallocateArray(geomech_grid%elem_type)
  call DeallocateArray(geomech_grid%elem_nodes)
  call DeallocateArray(geomech_grid%node_ids_natural)
  call DeallocateArray(geomech_grid%ghosted_node_ids_natural)
  
  if (associated(geomech_grid%nodes)) &
    deallocate(geomech_grid%nodes)
  nullify(geomech_grid%nodes)
  
  deallocate(geomech_grid)
  nullify(geomech_grid)
  
end subroutine GMGridDestroy

! ************************************************************************** !
!
! GMGridDestroy: Deallocates a geomechanics grid object
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
subroutine GMDMDestroy(gmdm)

  implicit none
  
  type(gmdm_type), pointer :: gmdm
  
  PetscErrorCode :: ierr
  
  if (.not.associated(gmdm)) return
  
  call ISDestroy(gmdm%is_ghosted_local,ierr)
  call ISDestroy(gmdm%is_local_local,ierr)
  call ISDestroy(gmdm%is_ghosted_petsc,ierr)
  call ISDestroy(gmdm%is_local_petsc,ierr)
  call ISDestroy(gmdm%is_ghosts_local,ierr)
  call ISDestroy(gmdm%is_ghosts_petsc,ierr)
  call ISDestroy(gmdm%is_local_natural,ierr)
  call VecScatterDestroy(gmdm%scatter_ltog,ierr)
  call VecScatterDestroy(gmdm%scatter_gtol,ierr)
  call VecScatterDestroy(gmdm%scatter_ltol,ierr)
  call VecScatterDestroy(gmdm%scatter_gton,ierr)
  call VecScatterDestroy(gmdm%scatter_ntog,ierr)
  call ISLocalToGlobalMappingDestroy(gmdm%mapping_ltog,ierr)
  if (gmdm%mapping_ltogb /= 0) &
    call ISLocalToGlobalMappingDestroy(gmdm%mapping_ltogb,ierr)
  call VecDestroy(gmdm%global_vec,ierr)
  call VecDestroy(gmdm%local_vec,ierr)
  
  deallocate(gmdm)
  nullify(gmdm)

end subroutine GMDMDestroy

end module Geomech_Grid_Aux_module
#endif 
!GEOMECH

