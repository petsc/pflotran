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

  PetscInt, parameter, public :: TWO_DIM_GRID = 1
  PetscInt, parameter, public :: THREE_DIM_GRID = 2 
  
  type, public :: geomech_grid_type
    ! variables for geomechanics grid (unstructured) for finite element formulation
    ! The dofs (displacements, here) are solved at the nodes 
    ! Notation corresponding to finite volume: node implies cell vertex and element means cell
    PetscInt :: global_offset                ! offset in petsc ordering for the first cell on a processor
    PetscInt :: nmax_elem                    ! Total number of elements in the global domain
    PetscInt :: nlmax_elem                   ! Total number of non-ghosted elements on a processor
    PetscInt :: nmax_node                    ! Total number of nodes in the global domain 
    PetscInt :: nlmax_node                   ! Total number of non-ghosted nodes on a processor
    PetscInt, pointer :: elem_ids_natural(:) ! Natural numbering of elements on a processor
    PetscInt, pointer :: elem_ids_petsc(:)   ! Petsc numbering of elements on a processor
    AO :: ao_natural_to_petsc                ! mapping of natural to Petsc ordering
    PetscInt :: max_ndual_per_elem           ! Max. number of dual elements connected to an element
    PetscInt :: max_nnode_per_elem           ! Max. number of nodes per element
    PetscInt :: max_elem_sharing_a_node      
    PetscInt, pointer :: elem_type(:)        ! Type of element
    PetscInt, pointer :: elem_node(:,:)      ! Node number on each element
    type(point_type), pointer :: nodes(:)             ! Coordinates of the nodes
  end type geomech_grid_type
  

  type, public :: gmdm_type                  ! Geomech. DM type
    ! Ghosting of elements not required
    PetscInt :: ndof
    ! local = local (non-ghosted) cells
    IS :: is_local_local                     ! IS for local cells with local on-processor numbering
    IS :: is_local_petsc                     ! IS for local cells with petsc numbering
    IS :: is_local_natural                   ! IS for local cells with natural (global) numbering
    VecScatter :: scatter_nton               ! scatter context for natural to natural updates
    VecScatter :: scatter_gton               ! scatter context for global to natural updates
    VecScatter :: scatter_ntog               ! scatter context for natural to global updates
    Vec :: global_vec                        ! global vec (no ghost cells), petsc-ordering
    Vec :: local_vec                         ! local vec (includes local cells), local ordering
  end type gmdm_type

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: GMGridCreate
  
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
  gmdm%is_local_local = 0
  gmdm%is_local_petsc = 0
  gmdm%is_local_natural = 0
  gmdm%scatter_nton = 0
  gmdm%scatter_gton = 0
  gmdm%scatter_ntog = 0
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
  nullify(geomech_grid%elem_node)
  nullify(geomech_grid%nodes)

  GMGridCreate => geomech_grid
  
end function GMGridCreate


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
  call DeallocateArray(geomech_grid%elem_node)

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
  
  call ISDestroy(gmdm%is_local_local,ierr)
  call ISDestroy(gmdm%is_local_petsc,ierr)
  call ISDestroy(gmdm%is_local_natural,ierr)
  call VecScatterDestroy(gmdm%scatter_nton,ierr)
  call VecScatterDestroy(gmdm%scatter_gton,ierr)
  call VecScatterDestroy(gmdm%scatter_ntog,ierr)
  call VecDestroy(gmdm%global_vec,ierr)
  call VecDestroy(gmdm%local_vec,ierr)
  
  deallocate(gmdm)
  nullify(gmdm)

end subroutine GMDMDestroy

end module Geomech_Grid_Aux_module
#endif !GEOMECH

