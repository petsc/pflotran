module Unstructured_Grid_module

  use Connection_module
  
  implicit none

  private 
  
#include "definitions.h"
#include "include/finclude/petsc.h"

  type, public :: unstructured_grid_type
    PetscInt :: num_cells
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash = 100
  end type
  
  public :: UnstructuredGridInit, &
            UnstructuredGridCreateDMs, &
            UnstGridComputeInternConnect, &
            UnstGridComputeBoundConnect, &
            UnstructGridGetGhostIdFromHash, &
            UnstructuredGridDestroy

contains

! ************************************************************************** !
!
! UnstructuredGridInit: Initializes an unstructured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine UnstructuredGridInit(unstructured_grid)

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid

end subroutine UnstructuredGridInit

! ************************************************************************** !
!
! StructuredGridCreateDMs: creates unstructured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine UnstructuredGridCreateDMs(unstructured_grid,option)
      
  use Option_module
      
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
end subroutine UnstructuredGridCreateDMs
  
! ************************************************************************** !
!
! UnstGridComputeInternConnect: computes internal connectivity of an  
!                                 unstructured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function UnstGridComputeInternConnect(unstructured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none
  
  type(connection_type), pointer :: UnstGridComputeInternConnect
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option

  nullify(UnstGridComputeInternConnect)
  
end function UnstGridComputeInternConnect

! ************************************************************************** !
!
! UnstGridComputeBoundConnect: computes boundary connectivity of an 
!                                 unstructured grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function UnstGridComputeBoundConnect(unstructured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none

  type(connection_type), pointer :: UnstGridComputeBoundConnect  
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
  nullify(UnstGridComputeBoundConnect)
  
end function UnstGridComputeBoundConnect

! ************************************************************************** !
!
! CreateNaturalToLocalHash: Creates a hash table for looking up the local 
!                           ghosted id of a natural id, if it exists
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine UnstGridCreateNatToGhostedHash(unstructured_grid,nG2A)

  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  PetscInt :: nG2A(:)

  PetscInt :: local_ghosted_id, natural_id 
  PetscInt :: num_in_hash, num_ids_per_hash, hash_id, id
  PetscInt, pointer :: temp_hash(:,:,:)

  ! initial guess of 10% of ids per hash
  num_ids_per_hash = max(unstructured_grid%nlmax/ &
                           (unstructured_grid%num_hash/10), &
                         unstructured_grid%nlmax)

  allocate(unstructured_grid%hash(2,0:num_ids_per_hash, &
                                  unstructured_grid%num_hash))
  unstructured_grid%hash(:,:,:) = 0

  ! recall that natural ids are zero-based
  do local_ghosted_id = 1, unstructured_grid%ngmax
    natural_id = nG2A(local_ghosted_id)+1
    hash_id = mod(natural_id,unstructured_grid%num_hash)+1 
    num_in_hash = unstructured_grid%hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,0:unstructured_grid%num_hash))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,unstructured_grid%num_hash) = &
                   unstructured_grid%hash(1:2,0:num_ids_per_hash, &
                                          unstructured_grid%num_hash)
      deallocate(unstructured_grid%hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(unstructured_grid%hash(1:2,0:num_ids_per_hash, &
                                      unstructured_grid%num_hash))
      ! copy old to new
      do hash_id = 1, unstructured_grid%num_hash
        do id = 1, temp_hash(1,0,hash_id)
          unstructured_grid%hash(1:2,id,hash_id) = temp_hash(1:2,id,hash_id)
        enddo
        unstructured_grid%hash(1,0,hash_id) = temp_hash(1,0,hash_id)
      enddo
      deallocate(temp_hash)
    endif
    unstructured_grid%hash(1,0,hash_id) = num_in_hash
    unstructured_grid%hash(1,num_in_hash,hash_id) = natural_id
    unstructured_grid%hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

!  if (grid%myrank == 0) print *, 'num_ids_per_hash:', num_ids_per_hash

end subroutine UnstGridCreateNatToGhostedHash

! ************************************************************************** !
!
! getLocalIdFromHash: Returns the local ghosted id of a natural id, if it 
!                     exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function UnstructGridGetGhostIdFromHash(unstructured_grid,natural_id)

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid

  PetscInt :: natural_id
  PetscInt :: hash_id, id

  UnstructGridGetGhostIdFromHash = 0
  hash_id = mod(natural_id,unstructured_grid%num_hash)+1 
  do id = 1, unstructured_grid%hash(1,0,hash_id)
    if (unstructured_grid%hash(1,id,hash_id) == natural_id) then
      UnstructGridGetGhostIdFromHash = unstructured_grid%hash(2,id,hash_id)
      return
    endif
  enddo

end function UnstructGridGetGhostIdFromHash

! ************************************************************************** !
!
! UnstructGridPrintHashTable: Prints the hashtable for viewing
! author: Glenn Hammond
! date: 03/09/07
!
! ************************************************************************** !
subroutine UnstructGridPrintHashTable(unstructured_grid)

  implicit none

  type(unstructured_grid_type) :: unstructured_grid

  PetscInt :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,unstructured_grid%num_hash
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         unstructured_grid%hash(1,0,ihash), &
                         '):', &
                         (unstructured_grid%hash(1,id,ihash),id=1, &
                          unstructured_grid%hash(1,0,ihash))
  enddo
  close(fid)

end subroutine UnstructGridPrintHashTable

! ************************************************************************** !
!
! UnstructuredGridDestroy: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine UnstructuredGridDestroy(unstructured_grid)

  implicit none
  
  type(unstructured_grid_type), pointer :: unstructured_grid
    
  if (.not.associated(unstructured_grid)) return
  

  if (associated(unstructured_grid%hash)) deallocate(unstructured_grid%hash)
  nullify(unstructured_grid%hash)

  deallocate(unstructured_grid)
  nullify(unstructured_grid)

end subroutine UnstructuredGridDestroy

end module Unstructured_Grid_module
