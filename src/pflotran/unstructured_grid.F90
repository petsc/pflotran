module Unstructured_Grid_module

  use Connection_module
  
  implicit none

  private 
  
#include "definitions.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"

  type, public :: unstructured_grid_type
    ! num_cells_ghosted =
    ! num_cells_global =
    ! num_cells_local = 
    PetscInt :: num_cells_global  ! number of cells in entire problem domain
    PetscInt :: num_cells_local   ! number of local (non-ghosted) cells
    PetscInt :: num_cells_ghosted ! number of local and ghosted cells on the process
    PetscInt :: num_ghost_cells   ! number of ghost cells (only) on processor
    PetscInt :: num_vertices_global ! number of vertices in entire problem domain
    PetscInt :: num_vertices_local  ! number of vertices in local grid cells
    PetscInt :: global_offset ! offset in petsc ordering for the first cell on a processor???
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash
    PetscInt, pointer :: cell_vertices_0(:,:) ! vertices for each grid cell (zero-based)
    PetscInt, pointer :: cell_ids_natural(:) ! natural 1d right-hand i,j,k ordering
    PetscInt, pointer :: cell_ids_petsc(:) ! petsc ordering of cell ids
    PetscInt, pointer :: ghost_cell_ids_natural(:) ! natural ordering of ghost cell ids
    PetscInt, pointer :: ghost_cell_ids_petsc(:) ! petsc ordering of ghost cells ids
    PetscInt, pointer :: cell_neighbors_local_ghosted(:,:) ! local neighbors
    type(point_type), pointer :: vertices(:)
    AO :: ao_natural_to_petsc ! mappsing of natural to petsc ordering
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
    IS :: is_ghosts_local ! IS for ghosted cells with local on-processor numbering
    IS :: is_ghosts_petsc ! IS for ghosted cells with petsc numbering
    IS :: is_local_natural ! IS for local cells with natural (global) numbering
    VecScatter :: scatter_ltog ! scatter context for local to global updates
    VecScatter :: scatter_gtol ! scatter context for global to local updates
    VecScatter :: scatter_ltol ! scatter context for local to local updates
    VecScatter :: scatter_gton ! scatter context for global to natural updates
    ISLocalToGlobalMapping :: mapping_ltog  ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltogb ! block form of mapping_ltog
    Vec :: global_vec ! global vec (no ghost cells), petsc-ordering
    Vec :: local_vec ! local vec (includes local and ghosted cells), local ordering
  end type ugdm_type

  type, public :: point_type
    PetscInt :: id
    PetscReal :: x
    PetscReal :: y
    PetscReal :: z
  end type point_type
  
  type, public :: plane_type
    PetscReal :: A
    PetscReal :: B
    PetscReal :: C
    PetscReal :: D
  end type plane_type

  PetscInt, parameter :: HEX_TYPE = 1

  public :: UnstructuredGridCreate, &
            UnstructuredGridCreateDM, &
            UnstGridComputeInternConnect, &
            UnstGridComputeBoundConnect, &
            UnstructGridGetGhostIdFromHash, &
            UnstructuredGridRead, &
            UnstructuredGridDecompose, &
            UGridComputeInternConnect, &
            UGridComputeCoord, &
            UGridComputeVolumes, &
            UGridMapIndices, &
            UGDMCreateJacobian, &
            UGDMCreateVector, &
            UnstructuredGridDestroy, &
            UnstructuredGridCreateUGDM, &
            UGDMDestroy

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
  ugdm%mapping_ltog = 0
  ugdm%mapping_ltogb = 0
  ugdm%global_vec = 0
  ugdm%local_vec = 0

  UGDMCreate => ugdm

end function UGDMCreate

! ************************************************************************** !
!
! UnstructuredGridCreate: Creates an unstructured grid object
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
function UnstructuredGridCreate()

  implicit none
  
  type(unstructured_grid_type), pointer :: UnstructuredGridCreate

  type(unstructured_grid_type), pointer :: unstructured_grid

  allocate(unstructured_grid)

  unstructured_grid%num_cells_global = 0
  unstructured_grid%num_vertices_global = 0
  unstructured_grid%num_cells_local = 0
  unstructured_grid%num_vertices_local = 0
  unstructured_grid%num_cells_ghosted = 0
  unstructured_grid%num_ghost_cells = 0
  unstructured_grid%global_offset = 0
  unstructured_grid%nmax = 0
  unstructured_grid%nlmax = 0
  unstructured_grid%ngmax = 0
  nullify(unstructured_grid%cell_vertices_0)
  nullify(unstructured_grid%vertices)
  nullify(unstructured_grid%hash)
  nullify(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_natural)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
  unstructured_grid%num_hash = 100
  unstructured_grid%ao_natural_to_petsc = 0

  UnstructuredGridCreate => unstructured_grid
  
end function UnstructuredGridCreate

! ************************************************************************** !
!
! StructuredGridCreateDMs: creates unstructured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine UnstructuredGridCreateDM()
      
  use Option_module
      
  implicit none
  
  PetscInt :: ndof
  PetscInt :: stencil_width
  
end subroutine UnstructuredGridCreateDM
  
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
  
  type(connection_set_type), pointer :: UnstGridComputeInternConnect
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

  type(connection_set_type), pointer :: UnstGridComputeBoundConnect  
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
! UnstructuredGridRead: Reads an unstructured grid
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UnstructuredGridRead(unstructured_grid,filename,option)

  use Input_module
  use Option_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card
  PetscInt :: num_cells_local_save
  PetscInt :: num_vertices_local_save
  PetscInt :: num_to_read
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscReal, allocatable :: vertex_coordinates(:,:)

  PetscInt :: icell, ivertex, idir, irank
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename)

  card = 'Unstructured Grid'

  call InputReadFlotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,card)  

  call InputReadInt(input,option,unstructured_grid%num_cells_global)
  call InputErrorMsg(input,option,'number of cells',card)
  call InputReadInt(input,option,unstructured_grid%num_vertices_global)
  call InputErrorMsg(input,option,'number of vertices',card)

  ! divide elements and vertices across processors
  unstructured_grid%num_cells_local = unstructured_grid%num_cells_global/ &
                                      option%mycommsize 
  num_cells_local_save = unstructured_grid%num_cells_local
  remainder = unstructured_grid%num_cells_global - &
              unstructured_grid%num_cells_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_cells_local = &
                                  unstructured_grid%num_cells_local + 1

  allocate(unstructured_grid%cell_vertices_0(8,unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices_0 = 0

  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(8,num_cells_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do icell = 1, num_to_read
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        do ivertex = 1, 8
          call InputReadInt(input,option,temp_int_array(ivertex,icell))
          call InputErrorMsg(input,option,'vertex id',card)
        enddo
      enddo
      
      if (irank == option%io_rank) then
print *, '0: ', unstructured_grid%num_cells_local, ' cells'
        unstructured_grid%cell_vertices_0(:,1:unstructured_grid%num_cells_local) = &
          temp_int_array(:,1:unstructured_grid%num_cells_local)
      else
print *, '0: ', num_to_read, ' cells sent'
        int_mpi = num_to_read*8
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
print *, option%myrank,': ',unstructured_grid%num_cells_local, ' cells recv'
    int_mpi = unstructured_grid%num_cells_local*8
    call MPI_Recv(unstructured_grid%cell_vertices_0,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif


  unstructured_grid%num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = unstructured_grid%num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              unstructured_grid%num_vertices_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_vertices_local = &
                                 unstructured_grid%num_vertices_local + 1

  allocate(vertex_coordinates(3,unstructured_grid%num_vertices_local))
  vertex_coordinates = 0.d0

  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(3,num_vertices_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_vertices_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do ivertex = 1, num_to_read
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        do idir = 1, 3
          call InputReadDouble(input,option,temp_real_array(idir,ivertex))
          call InputErrorMsg(input,option,'vertex coordinate',card)
        enddo
      enddo
      
      if (irank == option%io_rank) then
        vertex_coordinates(:,1:unstructured_grid%num_vertices_local) = &
          temp_real_array(:,1:unstructured_grid%num_vertices_local)
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    int_mpi = unstructured_grid%num_vertices_local*3
    call MPI_Recv(vertex_coordinates, &
                  int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif
  
  allocate(unstructured_grid%vertices(unstructured_grid%num_vertices_local))
  do ivertex = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ivertex)%id = 0
    unstructured_grid%vertices(ivertex)%x = vertex_coordinates(1,ivertex)
    unstructured_grid%vertices(ivertex)%y = vertex_coordinates(2,ivertex)
    unstructured_grid%vertices(ivertex)%z = vertex_coordinates(3,ivertex)
  enddo
  deallocate(vertex_coordinates)

  call InputDestroy(input)

end subroutine UnstructuredGridRead

! ************************************************************************** !
!
! UnstructuredGridDecompose: Decomposes an unstructured grid
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UnstructuredGridDecompose(unstructured_grid,option)
  
  use Option_module
  use Utility_module, only: reallocateIntArray, SearchOrderedArray
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  interface
  end interface
  
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
#ifdef ENABLE_UNSTRUCTURED  
  PetscInt :: icell, icell2
  PetscInt :: ivertex
  PetscInt :: vertex_id
  PetscInt :: count, vertex_count
  PetscInt :: vertex_offset, global_vertex_offset
  PetscInt :: max_vertex_count
  PetscInt :: max_dual
  PetscInt :: stride
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscInt, allocatable :: cell_counts(:)
  PetscInt, pointer :: index_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscInt, allocatable :: strided_indices(:)
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscInt :: num_rows, num_cols, istart, iend, icol
  PetscTruth :: success
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  PetscViewer :: viewer
  Mat :: Adj_mat
  Mat :: Dual_mat
  MatPartitioning :: Part
  Vec :: elements_natural
  Vec :: elements_petsc
  Vec :: elements_old
  Vec :: vertices_old
  Vec :: vertices_new
  IS :: is_new
  IS :: is_num
  IS :: is_scatter
  IS :: is_gather

  VecScatter :: vec_scatter
  
  PetscInt :: vertex_ids_offset
  PetscInt :: dual_offset

  PetscInt :: ghost_cell_count
  PetscInt :: max_ghost_cell_count
  PetscInt :: max_int_count
  PetscInt :: temp_int
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, pointer :: int_array_pointer(:)
  
  PetscInt :: idual, dual_id
  PetscTruth :: found


  
  ! cell distribution across processors (size = num_cores + 1)
  ! core i owns cells cell_distribution(i):cell_distribution(i+1), note
  ! the zero-based ordering
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%num_cells_local,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)
  
  ! provide ids of vertices for local cells

  vertex_ids_offset = 1 + 1 ! +1 for -777
  max_vertex_count = 8
  dual_offset = vertex_ids_offset + max_vertex_count + 1 ! +1 for -888
  max_dual = 6
  stride = dual_offset+ max_dual + 1 ! +1 for -999

  allocate(local_vertices(max_vertex_count*unstructured_grid%num_cells_local))
  allocate(local_vertex_offset(unstructured_grid%num_cells_local+1))
  count = 0
  local_vertex_offset(1) = 0
  do icell = 1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices_0(ivertex,icell)
    enddo
    local_vertex_offset(icell+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
  num_common_vertices = 3 ! cells must share at least vertices

  unstructured_grid%global_offset = 0
  call MPI_Exscan(unstructured_grid%num_cells_local, &
                  unstructured_grid%global_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

#if GEH_DEBUG  
  call printMsg(option,'Adjacency matrix')
#endif

  call MatCreateMPIAdj(option%mycomm,unstructured_grid%num_cells_local, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat,ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Adj.out',viewer,ierr)
  call MatView(Adj_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG
  call printMsg(option,'Dual matrix')
#endif

  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)
  call MatDestroy(Adj_mat,ierr)
  deallocate(local_vertices)
  deallocate(local_vertex_offset)
  
#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Dual.out',viewer,ierr)
  call MatView(Dual_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG
  call printMsg(option,'Partitioning')
#endif

  call MatPartitioningCreate(option%mycomm,Part,ierr)
  call MatPartitioningSetAdjacency(Part,Dual_mat,ierr)
  call MatPartitioningSetFromOptions(Part,ierr)
  call MatPartitioningApply(Part,is_new,ierr)
  call MatPartitioningDestroy(Part,ierr)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is.out',viewer,ierr)
  call ISView(is_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  
#endif

  allocate(cell_counts(option%mycommsize))
  call ISPartitioningCount(is_new,option%mycommsize,cell_counts,ierr)
  unstructured_grid%num_cells_local = cell_counts(option%myrank+1)
  deallocate(cell_counts)
    
  call VecCreate(option%mycomm,elements_natural,ierr)
  call VecSetSizes(elements_natural, &
                   stride*unstructured_grid%num_cells_local, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_natural,ierr)
  
  call ISPartitioningToNumbering(is_new,is_num,ierr)
  call ISDestroy(is_new,ierr)
  call ISGetIndicesF90(is_num,index_ptr,ierr)
  allocate(strided_indices(unstructured_grid%num_cells_local))
  do icell=1, unstructured_grid%num_cells_local
    strided_indices(icell) = stride*index_ptr(icell)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,stride, &
                     unstructured_grid%num_cells_local, &
                     strided_indices,is_scatter,ierr)
  call ISRestoreIndicesF90(is_num,index_ptr,ierr)
  deallocate(strided_indices)
  call ISDestroy(is_num,ierr)
  
  call VecCreate(option%mycomm,elements_old,ierr)
  call VecSetSizes(elements_old, &
                   stride*unstructured_grid%num_cells_local,PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_old,ierr)

  ! 0 = 0-based indexing
  call MatGetRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                      ja_ptr,success,ierr)

  if (success == PETSC_FALSE .or. &
      num_rows /= unstructured_grid%num_cells_local) then
    print *, option%myrank, num_rows, success, unstructured_grid%num_cells_local
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  call VecGetArrayF90(elements_old,vec_ptr,ierr)
  count = 0
  vertex_count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(unstructured_grid%global_offset+icell)
    count = count + 1
    vec_ptr(count) = -777  ! help differentiate
    do ivertex = 1, max_vertex_count
      count = count + 1
      vertex_count = vertex_count + 1
        ! increment for 1-based ordering
!      vec_ptr(count) = local_vertices(vertex_count) + 1
      vec_ptr(count) = unstructured_grid%cell_vertices_0(ivertex,icell) + 1
    enddo

    count = count + 1 
    vec_ptr(count) = -888  ! help differentiate

    istart = ia_ptr(icell)
    iend = ia_ptr(icell+1)-1
    num_cols = iend-istart+1
    do icol = 1, max_dual
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo
    count = count + 1 
    vec_ptr(count) = -999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call MatRestoreRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                          ja_ptr,success,ierr)
  call MatDestroy(Dual_mat,ierr)
 
#if GEH_DEBUG
  call printMsg(option,'Before element scatter')
#endif

  call VecScatterCreate(elements_old,PETSC_NULL,elements_natural,is_scatter,vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_natural,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_natural,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

#if GEH_DEBUG
  call printMsg(option,'After element scatter')
  call PetscViewerASCIIOpen(option%mycomm,'elements_old.out',viewer,ierr)
  call VecView(elements_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  call VecDestroy(elements_old,ierr)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'elements_natural.out',viewer,ierr)
  call VecView(elements_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  unstructured_grid%cell_vertices_0 = 0
  allocate(unstructured_grid%cell_ids_natural(unstructured_grid%num_cells_local))
  unstructured_grid%cell_ids_natural = 0
  
  ! look at all connections and determine how many are non-local, and create
  !  a listof indices
  call VecDuplicate(elements_natural,elements_petsc,ierr)
  call VecCopy(elements_natural,elements_petsc,ierr)

#if GEH_DEBUG
  call printMsg(option,'Lists of ids')
#endif

  ! make a list of local ids
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    unstructured_grid%cell_ids_natural(icell) = abs(vec_ptr((icell-1)*stride+1))
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

  ! make a list of petsc ids
  allocate(int_array(unstructured_grid%num_cells_local))
  do icell=1, unstructured_grid%num_cells_local
    int_array(icell) = icell+unstructured_grid%global_offset
  enddo
  
  int_array = int_array - 1
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural - 1
  call AOCreateBasic(option%mycomm,unstructured_grid%num_cells_local, &
                     unstructured_grid%cell_ids_natural,int_array, &
                     unstructured_grid%ao_natural_to_petsc,ierr)
  deallocate(int_array)
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural + 1

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'ao.out',viewer,ierr)
  call AOView(unstructured_grid%ao_natural_to_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! The below creates a list of cells ids for the duals and converts them
  ! to petsc ordering   
  
  ! count the number of duals  
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
    enddo
  enddo                  
  ! fill an array with the dual ids
  allocate(int_array(count))
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    int_array(count) = unstructured_grid%cell_ids_natural(icell)
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      int_array(count) = dual_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

#if GEH_DEBUG
  call printMsg(option,'Application ordering')
#endif

  ! convert the dual ids from natural to petsc
  int_array = int_array - 1             
  call AOApplicationToPetsc(unstructured_grid%ao_natural_to_petsc,count, &
                            int_array,ierr)
  int_array = int_array + 1                

#if GEH_DEBUG
  call printMsg(option,'PETSc-ordered duals')
#endif

  ! load mapped petsc-ordered dual ids back into duplicated vector
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  allocate(unstructured_grid%cell_ids_petsc(unstructured_grid%num_cells_local))
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    unstructured_grid%cell_ids_petsc(icell) = int_array(count)
    vec_ptr((icell-1)*stride+1) = int_array(count)
    do idual = 1, max_dual
      dual_id = vec_ptr2(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      vec_ptr(idual + dual_offset + (icell-1)*stride) = int_array(count)
    enddo
  enddo                
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  call VecRestoreArrayF90(elements_natural,vec_ptr2,ierr)
  deallocate(int_array)
  call VecDestroy(elements_natural,ierr)

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! make a list of ghosted ids in petsc numbering
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ghost_cell_count = 0
  max_ghost_cell_count = max(unstructured_grid%num_cells_local,100)
  allocate(unstructured_grid%ghost_cell_ids_petsc(max_ghost_cell_count))
  do icell=1, unstructured_grid%num_cells_local
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      found = PETSC_FALSE
      if (dual_id < 1) exit
#if 0  
      do icell2 = 1, unstructured_grid%num_cells_local
        if (dual_id == unstructured_grid%cell_ids_petsc(icell2)) then
          vec_ptr(idual + dual_offset + (icell-1)*stride) = icell2
          found = PETSC_TRUE
          exit
        endif
      enddo
#endif
                                     ! global_offset is zero-based
    if (dual_id <= unstructured_grid%global_offset .or. &
        dual_id > unstructured_grid%global_offset + &
                  unstructured_grid%num_cells_local) then
      ! not located on-processor
        
      do icell2 = 1, ghost_cell_count
          if (dual_id == unstructured_grid%ghost_cell_ids_petsc(icell2)) then
            ! flag the id as negative for ghost cell and set to current count
            vec_ptr(idual + dual_offset + (icell-1)*stride) = -icell2
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (.not.found) then
          ghost_cell_count = ghost_cell_count + 1
          if (ghost_cell_count > max_ghost_cell_count) then
            call reallocateIntArray(unstructured_grid%ghost_cell_ids_petsc, &
                                    max_ghost_cell_count)
          endif
          unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count) = dual_id
          ! flag the id as negative for ghost cell and set to current count
          vec_ptr(idual + dual_offset + (icell-1)*stride) = -ghost_cell_count
        endif
      endif
    enddo
  enddo
  
  call printMsg(option,'Add code adds all ghost cells and removed duplicates')
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  unstructured_grid%num_ghost_cells = ghost_cell_count
  unstructured_grid%num_cells_ghosted = &
    unstructured_grid%num_cells_local + ghost_cell_count

  ! sort ghost cell ids
  allocate(int_array(ghost_cell_count))
  allocate(int_array2(ghost_cell_count))
  do icell = 1, ghost_cell_count
    int_array(icell) = unstructured_grid%ghost_cell_ids_petsc(icell)
    int_array2(icell) = icell
  enddo
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(ghost_cell_count,int_array,int_array2,ierr)
  int_array2 = int_array2+1
  ! resize ghost cell array down to ghost_cell_count
  deallocate(unstructured_grid%ghost_cell_ids_petsc)
  allocate(unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count))
  do icell = 1, ghost_cell_count
    unstructured_grid%ghost_cell_ids_petsc(int_array2(icell)) = int_array(icell)
  enddo
  ! now, rearrange the ghost cell ids of the dual accordingly
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ! add num_cells_local to get in local numbering
  int_array2 = int_array2 + unstructured_grid%num_cells_local
  do icell=1, unstructured_grid%num_cells_local
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 0) then
        vec_ptr(idual + dual_offset + (icell-1)*stride) = int_array2(-dual_id)
      endif       
    enddo
  enddo
  deallocate(int_array)
  deallocate(int_array2)
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_local.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! load cell neighbors into array
  ! start first index at zero to store # duals for a cell
  allocate(unstructured_grid%cell_neighbors_local_ghosted(0:6, &
                    unstructured_grid%num_cells_local))
  unstructured_grid%cell_neighbors_local_ghosted = 0

  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    count = 0
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      ! flag ghosted cells as negative
      if (dual_id > unstructured_grid%num_cells_local) dual_id = -dual_id
      unstructured_grid%cell_neighbors_local_ghosted(idual,icell) = dual_id
    enddo
    unstructured_grid%cell_neighbors_local_ghosted(0,icell) = count
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! make a list of local vertices
  max_int_count = unstructured_grid%num_cells_local
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride)
      if (vertex_id < 1) exit
      vertex_count = vertex_count + 1
      if (vertex_count > max_int_count) then
        call reallocateIntArray(int_array_pointer,max_int_count)
      endif
      vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride) = vertex_count
      int_array_pointer(vertex_count) = vertex_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! sort the vertex ids
  allocate(int_array(vertex_count))
  int_array(1:vertex_count) = int_array_pointer(1:vertex_count)
  allocate(int_array2(vertex_count))
  do ivertex = 1, vertex_count
    int_array2(ivertex) = ivertex 
  enddo
  deallocate(int_array_pointer)
  nullify(int_array_pointer)
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(vertex_count,int_array,int_array2,ierr)
  int_array2 = int_array2+1

  ! remove duplicates
  allocate(int_array3(vertex_count))
  allocate(int_array4(vertex_count))
  int_array3 = 0
  int_array4 = 0
  int_array3(1) = int_array(int_array2(1))
  int_array4(1) = 1
  count = 1
  do ivertex = 2, vertex_count
    vertex_id = int_array(int_array2(ivertex))
    if (vertex_id > int_array3(count)) then
      count = count + 1
      int_array3(count) = vertex_id
    endif
    int_array4(int_array2(ivertex)) = count
  enddo
  vertex_count = count
  deallocate(int_array)

  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  deallocate(unstructured_grid%cell_vertices_0)
  allocate(unstructured_grid%cell_vertices_0(0:max_vertex_count, &
                                             unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices_0 = -999
  unstructured_grid%cell_vertices_0(0,:) = 0
  
  ! permute the local ids set earlier
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride)
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices_0(0,icell)+1
      unstructured_grid%cell_vertices_0(count,icell) = int_array4(vertex_id)-1
      unstructured_grid%cell_vertices_0(0,icell) = count
      vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride) = int_array4(vertex_id)
    enddo
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  deallocate(int_array2)
  deallocate(int_array4)

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_vert_local.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  call VecDestroy(elements_petsc,ierr)

  ! IS for gather - need local numbering
  allocate(strided_indices(vertex_count))
  do ivertex = 1, vertex_count
    strided_indices(ivertex) = 3*(ivertex-1)
  enddo
  deallocate(int_array3)
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     strided_indices,is_gather,ierr)
  deallocate(strided_indices)

  call VecCreateMPI(option%mycomm,unstructured_grid%num_vertices_local*3, &
                    PETSC_DETERMINE,vertices_old,ierr)
  call VecSetBlockSize(vertices_old,3,ierr)

  call VecCreateSeq(PETSC_COMM_SELF,vertex_count*3,vertices_new,ierr)
  call VecSetBlockSize(vertices_new,3,ierr)

!  call VecCreate(option%mycomm,vertices_new,ierr)
!  call VecSetSizes(vertices_new, &
!                   vertex_count*3,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_new,ierr)
  
!  call VecCreate(option%mycomm,vertices_old,ierr)
!  call VecSetSizes(vertices_old, &
!                   3*unstructured_grid%num_vertices_local,PETSC_DECIDE,ierr)
!  call VecSetFromOptions(vertices_old,ierr)
  call VecGetArrayF90(vertices_old,vec_ptr,ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    vec_ptr((ivertex-1)*3+1) = unstructured_grid%vertices(ivertex)%x
    vec_ptr((ivertex-1)*3+2) = unstructured_grid%vertices(ivertex)%y
    vec_ptr((ivertex-1)*3+3) = unstructured_grid%vertices(ivertex)%z
  enddo
  call VecRestoreArrayF90(vertices_old,vec_ptr,ierr)
  deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)

#if 0
  global_vertex_offset = 0
  call MPI_Exscan(unstructured_grid%num_vertices_local, &
                  global_vertex_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
#endif

  ! IS for scatter - provide petsc gobal numbering
  allocate(strided_indices(vertex_count))
  do ivertex = 1, vertex_count
    strided_indices(ivertex) = 3*(needed_vertices_petsc(ivertex)-1)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     strided_indices,is_scatter,ierr)
  deallocate(strided_indices)

  ! resize vertex array to new size
  unstructured_grid%num_vertices_local = vertex_count
  allocate(unstructured_grid%vertices(vertex_count))
  do ivertex = 1, vertex_count
    unstructured_grid%vertices(ivertex)%id = 0
    unstructured_grid%vertices(ivertex)%x = 0.d0
    unstructured_grid%vertices(ivertex)%y = 0.d0
    unstructured_grid%vertices(ivertex)%z = 0.d0
  enddo

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter.out',viewer,ierr)
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'is_gather.out',viewer,ierr)
  call ISView(is_gather,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecScatterCreate(vertices_old,is_scatter,vertices_new,is_gather, &
                        vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call ISDestroy(is_gather,ierr)
  call VecScatterBegin(vec_scatter,vertices_old,vertices_new, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,vertices_old,vertices_new, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old.out',viewer,ierr)
  call VecView(vertices_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecDestroy(vertices_old,ierr)


  call VecGetArrayF90(vertices_new,vec_ptr,ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    unstructured_grid%vertices(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    unstructured_grid%vertices(ivertex)%z = vec_ptr((ivertex-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_new,vec_ptr,ierr)
  
#if GEH_DEBUG  
  write(string,*) option%myrank
  string = 'vertex_coord_new' // trim(adjustl(string)) // '.out'
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call VecView(vertices_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecDestroy(vertices_new,ierr)

#if 0
  ! calculate faces
  allocate(local_vertices(max_vertex_count*unstructured_grid%num_cells_local))
  allocate(local_vertex_offset(unstructured_grid%num_cells_local+1))
  count = 0
  local_vertex_offset(1) = 0
  do icell = 1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices_0(ivertex,icell)
    enddo
    local_vertex_offset(icell+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
  num_common_vertices = 3 ! cells must share at least vertices

#if GEH_DEBUG
  call printMsg(option,'Local Adjacency matrix')
#endif

  call MatCreateMPIAdj(PETSC_COMM_SELF,unstructured_grid%num_cells_local, &
                       unstructured_grid%num_vertices_local, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat,ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Adj_local.out',viewer,ierr)
  call MatView(Adj_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG
  call printMsg(option,'Dual matrix')
#endif

  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)
  call MatDestroy(Adj_mat,ierr)
  deallocate(local_vertices)
  deallocate(local_vertex_offset)

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Dual_local.out',viewer,ierr)
  call MatView(Dual_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! 0 = 0-based indexing
  call MatGetRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                      ja_ptr,success,ierr)

  call VecGetArrayF90(elements_old,vec_ptr,ierr)
  count = 0
  vertex_count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(unstructured_grid%global_offset+icell)
    count = count + 1
    vec_ptr(count) = -777  ! help differentiate
    do ivertex = 1, max_vertex_count
      count = count + 1
      vertex_count = vertex_count + 1
        ! increment for 1-based ordering
!      vec_ptr(count) = local_vertices(vertex_count) + 1
      vec_ptr(count) = unstructured_grid%cell_vertices_0(ivertex,icell) + 1
    enddo

    count = count + 1 
    vec_ptr(count) = -888  ! help differentiate

    istart = ia_ptr(icell)
    iend = ia_ptr(icell+1)-1
    num_cols = iend-istart+1
    do icol = 1, max_dual
      count = count + 1
      if (icol <= num_cols) then
        ! increment for 1-based ordering
        vec_ptr(count) = ja_ptr(icol+istart) + 1
      else
        vec_ptr(count) = 0
      endif
    enddo
    count = count + 1 
    vec_ptr(count) = -999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call MatRestoreRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                          ja_ptr,success,ierr)
  call MatDestroy(Dual_mat,ierr)

#endif

  unstructured_grid%nlmax = unstructured_grid%num_cells_local
  unstructured_grid%ngmax = unstructured_grid%num_cells_local + &
                            unstructured_grid%num_ghost_cells

#endif
  
end subroutine UnstructuredGridDecompose

! ************************************************************************** !
!
! UnstructuredGridDecompose: Decomposes an unstructured grid
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UnstructuredGridCreateUGDM(unstructured_grid,ugdm,ndof,option)
  
  use Option_module
  use Utility_module, only: reallocateIntArray
  
  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscda.h"
#include "finclude/petscda.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscviewer.h"

  interface
  end interface
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type), pointer :: ugdm
  PetscInt :: ndof
  type(option_type) :: option
  
#ifdef ENABLE_UNSTRUCTURED  
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: istart, iend
  PetscInt :: icell
  PetscInt :: idof
  IS :: is_tmp
  Vec :: vec_tmp
  PetscErrorCode :: ierr
  
  PetscViewer :: viewer

  PetscInt, allocatable :: int_array(:)
  
  ugdm => UGDMCreate()
  ugdm%ndof = ndof

#if GEH_DEBUG
  call printMsg(option,'Vectors')
#endif

  ! create global vec
  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local*ndof, &
                    PETSC_DETERMINE,ugdm%global_vec,ierr)
  call VecSetBlockSize(ugdm%global_vec,ndof,ierr)
  ! create local vec
  call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%num_cells_ghosted*ndof, &
                    ugdm%local_vec,ierr)
  call VecSetBlockSize(ugdm%local_vec,ndof,ierr)
  
  ! IS for global numbering of local, non-ghosted cells
  call VecGetOwnershipRange(ugdm%global_vec,istart,iend,ierr)
!  call ISCreateStride(option%mycomm,unstructured_grid%num_cells_local, &
!                      istart,ndof,ugdm%is_local_petsc,ierr)
  allocate(int_array(unstructured_grid%num_cells_local))
  do icell = 1, unstructured_grid%num_cells_local
    int_array(icell) = (icell-1)*ndof+istart
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,ugdm%is_local_petsc,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_local_petsc.out',viewer,ierr)
  call ISView(ugdm%is_local_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do icell = 1, unstructured_grid%num_ghost_cells
    int_array(icell) = (icell+unstructured_grid%num_cells_local-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,ugdm%is_ghosts_local,ierr)
!  call ISCreateGeneral(option%mycomm,unstructured_grid%num_ghost_cells, &
!                       int_array,ugdm%is_ghosts_local,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_ghosts_local.out',viewer,ierr)
  call ISView(ugdm%is_ghosts_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
#if GEH_DEBUG
  call printMsg(option,'Index Sets')
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do icell = 1, unstructured_grid%num_ghost_cells
    int_array(icell) = (unstructured_grid%ghost_cell_ids_petsc(icell)-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,ugdm%is_ghosts_petsc,ierr)
!  call ISCreateGeneral(option%mycomm,unstructured_grid%num_ghost_cells, &
!                       int_array,ugdm%is_ghosts_petsc,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_ghosts_petsc.out',viewer,ierr)
  call ISView(ugdm%is_ghosts_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! IS for local numbering of local, non-ghosted cells
  allocate(int_array(unstructured_grid%num_cells_local))
  do icell = 1, unstructured_grid%num_cells_local
    int_array(icell) = (icell-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,ugdm%is_local_local,ierr)
!  call ISCreateGeneral(option%mycomm,unstructured_grid%num_cells_local, &
!                       int_array,ugdm%is_local_local,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_local_local.out',viewer,ierr)
  call ISView(ugdm%is_local_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! IS for ghosted numbering of local ghosted cells
  allocate(int_array(unstructured_grid%num_cells_ghosted))
  do icell = 1, unstructured_grid%num_cells_ghosted
    int_array(icell) = (icell-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
                     int_array,ugdm%is_ghosted_local,ierr)
!  call ISCreateGeneral(option%mycomm,unstructured_grid%num_cells_ghosted, &
!                       int_array,ugdm%is_ghosted_local,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_local.out',viewer,ierr)
  call ISView(ugdm%is_ghosted_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
             
  ! IS for petsc numbering of local ghosted cells
  allocate(int_array(unstructured_grid%num_cells_ghosted))
  do icell = 1, unstructured_grid%num_cells_local
    int_array(icell) = istart+(icell-1)*ndof
  enddo
  do icell = 1,unstructured_grid%num_ghost_cells
    int_array(unstructured_grid%num_cells_local+icell) = &
      (unstructured_grid%ghost_cell_ids_petsc(icell)-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
                     int_array,ugdm%is_ghosted_petsc,ierr)
!  call ISCreateGeneral(option%mycomm,unstructured_grid%num_cells_ghosted, &
!                       int_array,ugdm%is_ghosted_petsc,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_ghosted_petsc.out',viewer,ierr)
  call ISView(ugdm%is_ghosted_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif    
                 
  ! create a local to global mapping
#if GEH_DEBUG
  call printMsg(option,'ISLocalToGlobalMapping')
#endif

  call ISLocalToGlobalMappingCreateIS(ugdm%is_ghosted_petsc, &
                                      ugdm%mapping_ltog,ierr)

#if GEH_DEBUG                                      
  call PetscViewerASCIIOpen(option%mycomm,'mapping_ltog.out',viewer,ierr)
  call ISLocalToGlobalMappingView(ugdm%mapping_ltog,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
               
  ! create a block local to global mapping 
#if GEH_DEBUG
  call printMsg(option,'ISLocalToGlobalMappingBlock')
#endif

  call ISLocalToGlobalMappingBlock(ugdm%mapping_ltog,ndof, &
                                   ugdm%mapping_ltogb,ierr)
                                      
#if GEH_DEBUG                                      
  call PetscViewerASCIIOpen(option%mycomm,'mapping_ltogb.out',viewer,ierr)
  call ISLocalToGlobalMappingView(ugdm%mapping_ltogb,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG                            
  call printMsg(option,'local to global')
#endif

  ! Create local to global scatter
  call VecScatterCreate(ugdm%local_vec,ugdm%is_local_local,ugdm%global_vec, &
                        ugdm%is_local_petsc,ugdm%scatter_ltog,ierr)
                        
#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'scatter_ltog.out',viewer,ierr)
  call VecScatterView(ugdm%scatter_ltog,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG
  call printMsg(option,'global to local')
#endif

  ! Create global to local scatter
  call VecScatterCreate(ugdm%global_vec,ugdm%is_ghosted_petsc,ugdm%local_vec, &
                        ugdm%is_ghosted_local,ugdm%scatter_gtol,ierr)
                        
#if GEH_DEBUG                        
  call PetscViewerASCIIOpen(option%mycomm,'scatter_gtol.out',viewer,ierr)
  call VecScatterView(ugdm%scatter_gtol,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if GEH_DEBUG
  call printMsg(option,'local to local')
#endif
  
  ! Create local to local scatter.  Essentially remap the global to local as
  ! PETSc does in daltol.c
  call VecScatterCopy(ugdm%scatter_gtol,ugdm%scatter_ltol,ierr)
  call ISGetIndicesF90(ugdm%is_local_local,int_ptr,ierr)
  call VecScatterRemap(ugdm%scatter_ltol,int_ptr,PETSC_NULL_INTEGER,ierr)
  call ISRestoreIndicesF90(ugdm%is_local_local,int_ptr,ierr)
!  call VecScatterCreate(local_vec,is_local_petsc,local_vec,is_ghosts_petsc,scatter_ltol,ierr)

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'scatter_ltol.out',viewer,ierr)
  call VecScatterView(ugdm%scatter_ltol,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! Set up global to natural scatter
  ! Create index set of local non-ghosted Petsc ordering
  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local, &
                    PETSC_DETERMINE,vec_tmp,ierr)
  call VecGetOwnershipRange(vec_tmp,istart,iend,ierr)
  call VecDestroy(vec_tmp,ierr)
  allocate(int_array(unstructured_grid%num_cells_local))
  do icell = 1, unstructured_grid%num_cells_local
!    int_array(icell) = (icell-1)*ndof+istart
    int_array(icell) = (icell-1)+istart
  enddo
  call ISCreateGeneral(option%mycomm,unstructured_grid%num_cells_local, &
                       int_array,is_tmp,ierr)
  deallocate(int_array)
  call AOPetscToApplicationIS(unstructured_grid%ao_natural_to_petsc, &
                              is_tmp,ierr)
  ! remap for ndof > 1
  allocate(int_array(unstructured_grid%num_cells_local))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr)
  do icell = 1, unstructured_grid%num_cells_local
    int_array(icell) = int_ptr(icell)*ndof
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr)
  call ISDestroy(is_tmp,ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,ugdm%is_local_natural,ierr)
  deallocate(int_array)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_local_natural.out',viewer,ierr)
  call ISView(ugdm%is_local_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local*ndof, &
                    PETSC_DETERMINE,vec_tmp,ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_petsc,vec_tmp, &
                        ugdm%is_local_natural,ugdm%scatter_gton,ierr)
  call VecDestroy(vec_tmp,ierr)

#if GEH_DEBUG                        
  call PetscViewerASCIIOpen(option%mycomm,'scatter_gton.out',viewer,ierr)
  call VecScatterView(ugdm%scatter_gton,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  
    
#endif
  
end subroutine UnstructuredGridCreateUGDM

! ************************************************************************** !
!
! UGridComputeInternConnect: computes internal connectivity of an
!                            unstructured grid
! author: Glenn Hammond
! date: 10/21/09
!
! ************************************************************************** !
function UGridComputeInternConnect(unstructured_grid,grid_x,grid_y,grid_z, &
                                   option)

  use Connection_module
  use Option_module
  use Utility_module, only : DotProduct, CrossProduct

  implicit none

  type(connection_set_type), pointer :: UGridComputeInternConnect
  type(option_type) :: option
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(unstructured_grid_type) :: unstructured_grid

#ifdef ENABLE_UNSTRUCTURED
  PetscInt :: nconn, iconn
  PetscInt :: idual, dual_id

  PetscInt, allocatable :: face_to_vertex(:,:)
  PetscInt, allocatable :: cell_to_face(:,:)
  PetscInt, allocatable :: face_to_cell(:,:)
  PetscInt, allocatable :: vertex_to_cell(:,:)
  PetscInt, allocatable :: dual_to_face(:)
  PetscInt, allocatable :: temp_int(:)
  PetscInt, allocatable :: temp_int_2d(:,:)
  PetscInt :: num_match
  PetscInt :: found_count
  PetscTruth :: found
  PetscTruth :: match_found
  PetscInt :: face_count
  PetscInt :: count
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: icell, icell2
  PetscInt :: cell_id, cell_id2
  PetscInt :: dual_icell
  PetscInt :: ivertex, ivertex2
  PetscInt :: vertex_id, vertex_id2
  PetscInt :: vertex_ids4(4)
  
  PetscReal :: v1(3), v2(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  
  type(plane_type) :: plane1, plane2
  type(point_type) :: point1, point2, point3, point4
  type(point_type) :: point_up, point_dn
  type(point_type) :: intercept1, intercept2
  
  
  character(len=MAXSTRINGLENGTH) :: string  

  type(connection_set_type), pointer :: connections

  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(4,unstructured_grid%num_cells_ghosted*6))
  face_to_vertex = -999
  allocate(cell_to_face(6,unstructured_grid%num_cells_ghosted))
  cell_to_face = -999
  allocate(face_to_cell(2,6*unstructured_grid%num_cells_ghosted))
  face_to_cell = -999
  allocate(vertex_to_cell(0:8,unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  face_count = 0
  do icell = 1, unstructured_grid%num_cells_ghosted
    do iface = 1, 6
      face_count = face_count + 1
      cell_to_face(iface,icell) = face_count
      face_to_cell(1,face_count) = icell
      call UGridGetCellFaceVertices(option,HEX_TYPE,iface,vertex_ids4)
      do ivertex = 1, 4
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices_0(vertex_ids4(ivertex),icell)+1
      enddo
    enddo
  enddo

  ! remove duplicate faces
  ! fill face ids
  do icell = 1, unstructured_grid%num_cells_ghosted
    cell_id = icell
    do icell2 = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      cell_id2 = unstructured_grid%cell_neighbors_local_ghosted(icell2,icell)
      if (cell_id2 <= cell_id) cycle
      num_match = 0
      do ivertex = 1, unstructured_grid%cell_vertices_0(0,cell_id)
        vertex_id = unstructured_grid%cell_vertices_0(ivertex,cell_id)+1
        do ivertex2 = 1, unstructured_grid%cell_vertices_0(0,cell_id2)
          vertex_id2 = unstructured_grid%cell_vertices_0(ivertex2,cell_id2)+1
          if (vertex_id == vertex_id2) then
            num_match = num_match + 1
            vertex_ids4(num_match) = vertex_id
            exit
          endif
        enddo
        if (num_match > 3) exit
      enddo
      if (num_match > 3) then
        ! now find the shared face
        match_found = PETSC_TRUE
        do iface = 1, 6
          face_id = cell_to_face(iface,cell_id)
          found_count = 0
          do ivertex = 1, 4
            do ivertex2 = 1, 4
              if (face_to_vertex(ivertex,face_id) == vertex_ids4(ivertex2)) then
                found_count = found_count + 1
                exit
              endif
            enddo
          enddo
          if (found_count > 3) exit
        enddo
        if (found_count > 3) then
          do iface2 = 1, 6
            face_id2 = cell_to_face(iface2,cell_id2)
            found_count = 0
            do ivertex = 1, 4
              do ivertex2 = 1, 4
                if (face_to_vertex(ivertex,face_id2) == vertex_ids4(ivertex2)) then
                  found_count = found_count + 1
                  exit
                endif
              enddo
            enddo
            if (found_count > 3) exit
          enddo
          if (found_count <= 3) then
            match_found = PETSC_FALSE
          else
            ! remove duplicate face
            if (face_id2 > face_id) then
              write(string,*) face_id2, ' -> ', face_id
              option%io_buffer = 'Duplicated face removed:' // string
              call printMsg(option)
              cell_to_face(iface2,cell_id2) = face_id
              ! flag as removed
              face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
              face_to_cell(2,face_id) = cell_id2
            else
              write(string,*) face_id, ' -> ', face_id2
              option%io_buffer = 'Duplicated face removed:' // string
              call printMsg(option)
              cell_to_face(iface,cell_id) = face_id2
              ! flag as removed
              face_to_cell(1,face_id) = -face_to_cell(1,face_id)
              face_to_cell(2,face_id2) = cell_id
            endif
          endif
        else
          match_found = PETSC_FALSE
        endif  
        if (match_found == PETSC_FALSE) then
          option%io_buffer = 'Matching faces not found'
          call printErrMsg(option)
        endif
      endif    
    enddo
  enddo
  
  ! count up the # of faces
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) &
      face_count = face_count + 1
  enddo
  print *, 'Face count:', face_count
  ! reallocate face_to_vertex
  allocate(temp_int_2d(4,face_count))
  count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      count = count + 1
      temp_int_2d(:,count) = face_to_vertex(:,iface)
    endif
  enddo
  deallocate(face_to_vertex)
  allocate(face_to_vertex(4,face_count))
  face_to_vertex = temp_int_2d
  deallocate(temp_int_2d)
  ! reallocate face_to_cell to proper size
  allocate(temp_int_2d(2,face_count))
  allocate(temp_int(size(face_to_cell,2)))
  temp_int = 0
  count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      count = count + 1
      temp_int_2d(:,count) = face_to_cell(:,iface)
      temp_int(iface) = count
    endif
  enddo
  deallocate(face_to_cell)
  allocate(face_to_cell(2,face_count))
  face_to_cell = temp_int_2d
  deallocate(temp_int_2d)
  
  ! remap faces in cells using temp_int from above
  do iface = 1, size(face_to_cell,2)
    face_id = iface
    do icell = 1,2
      cell_id = face_to_cell(icell,face_id)
      ! check for exterior face
      if (cell_id < 0) cycle
      found = PETSC_FALSE
      do iface2 = 1, 6
        face_id2 = cell_to_face(iface2,cell_id)
        if (face_id == temp_int(face_id2)) then
          found = PETSC_TRUE
          cell_to_face(iface2,cell_id) = face_id
          exit
        endif
      enddo
      if (found == PETSC_FALSE) then
        option%io_buffer = 'Remapping of cell face id unsuccessful'
        call printErrMsg(option)
      endif
    enddo
  enddo
  deallocate(temp_int)
  
  
  do icell = 1, unstructured_grid%num_cells_ghosted
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,icell)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,icell)+1
      count = vertex_to_cell(0,vertex_id) + 1
      vertex_to_cell(count,vertex_id) = icell
      vertex_to_cell(0,vertex_id) = count
    enddo
  enddo
  
  nconn = 0
  do icell = 1, unstructured_grid%num_cells_local
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      dual_id = unstructured_grid%cell_neighbors_local_ghosted(idual,icell)
      ! count all ghosted connections (dual_id < 0)
      ! only count connection with cells of larger ids to avoid double counts
      if (dual_id < 0 .or. icell < dual_id) then
        nconn = nconn + 1
      endif
    enddo
  enddo

  connections => ConnectionCreate(nconn,option%nphase,INTERNAL_CONNECTION_TYPE)
  
  allocate(dual_to_face(nconn))
  dual_to_face = -999

  ! loop over connection again
  iconn = 0
  do icell = 1, unstructured_grid%num_cells_local
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      dual_icell = unstructured_grid%cell_neighbors_local_ghosted(idual,icell)
      if (icell < dual_icell) then 
        iconn = iconn + 1
        ! find face
        found = PETSC_FALSE
        do iface = 1, 6
          face_id = cell_to_face(iface,icell)
          do icell2 = 1,2
            cell_id2 = face_to_cell(icell2,face_id)
            if (cell_id2 == abs(dual_icell)) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (found == PETSC_TRUE) exit
        enddo
        if (found == PETSC_TRUE) then
          dual_to_face(iconn) = face_id
        else
          option%io_buffer = 'face not found in connection loop'
          call printErrMsg(option)
        endif
        connections%id_up(iconn) = icell
        connections%id_dn(iconn) = abs(dual_icell)
        ! need to add the surface areas, distance, etc.
        point1 = unstructured_grid%vertices(face_to_vertex(1,face_id))
        point2 = unstructured_grid%vertices(face_to_vertex(2,face_id))
        point3 = unstructured_grid%vertices(face_to_vertex(3,face_id))
        point4 = unstructured_grid%vertices(face_to_vertex(4,face_id))
        
        call ComputePlane(plane1,point1,point2,point3)
        call ComputePlane(plane2,point3,point4,point1)
        
        point_up%x = grid_x(icell)
        point_up%y = grid_y(icell)
        point_up%z = grid_z(icell)
        point_dn%x = grid_x(abs(dual_icell))
        point_dn%y = grid_y(abs(dual_icell))
        point_dn%z = grid_z(abs(dual_icell))
        v1(1) = point_dn%x-point_up%x
        v1(2) = point_dn%y-point_up%y
        v1(3) = point_dn%z-point_up%z
        n_up_dn = v1 / DotProduct(v1,v1)
        call GetPlaneIntercept(plane1,point_up,point_dn,intercept1)
        call GetPlaneIntercept(plane2,point_up,point_dn,intercept2)
        
        v1(1) = point3%x-point2%x
        v1(2) = point3%y-point2%y
        v1(3) = point3%z-point2%z
        v2(1) = point1%x-point2%x
        v2(2) = point1%y-point2%y
        v2(3) = point1%z-point2%z
        n1 = CrossProduct(v1,v2)
        area1 = 0.5d0*DotProduct(n1,n1)
        area1 = area1*DotProduct(n1,n_up_dn)
        
        v1(1) = point1%x-point4%x
        v1(2) = point1%y-point4%y
        v1(3) = point1%z-point4%z
        v2(1) = point3%x-point4%x
        v2(2) = point3%y-point4%y
        v2(3) = point3%z-point4%z
        n2 = CrossProduct(v1,v2)
        area2 = 0.5d0*DotProduct(n2,n2)
        area2 = area2*DotProduct(n2,n_up_dn)
       !area1 = 0.5*|(p2-p1)X(p3-p1)|.
        !  http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#Triangles
        
        !GetPlaneNormalArea(
        ! http://mathworld.wolfram.com/Plane.html

        v1(1) = intercept1%x-point_up%x
        v1(2) = intercept1%y-point_up%y
        v1(3) = intercept1%z-point_up%z
        v2(1) = point_dn%x-intercept1%x
        v2(2) = point_dn%y-intercept1%y
        v2(3) = point_dn%z-intercept1%z
        dist_up = sqrt(DotProduct(v1,v1))
        dist_dn = sqrt(DotProduct(v2,v2))
        
        connections%dist(-1:3,iconn) = 0.d0
        connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connections%dist(0,iconn) = dist_up + dist_dn
        connections%dist(1:3,iconn) = v1 + v2
        connections%area(iconn) = area1 + area2
       
      endif
    enddo
  enddo

  deallocate(face_to_vertex)
  deallocate(cell_to_face)
  deallocate(face_to_cell)
  deallocate(vertex_to_cell)
  deallocate(dual_to_face)

  UGridComputeInternConnect => connections
#endif        

end function UGridComputeInternConnect

! ************************************************************************** !
!
! UGridGetCellFaceVertices: returns vertex ids of a given hex face
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine UGridGetCellFaceVertices(option,cell_type,iface,vertex_ids)

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscInt :: cell_type
  PetscInt :: iface
  PetscInt :: vertex_ids(*)
  
  select case(cell_type)
    case(HEX_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 6
          vertex_ids(4) = 5
        case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 7
          vertex_ids(4) = 6
        case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 4
          vertex_ids(3) = 8
          vertex_ids(4) = 7
        case(4)
          vertex_ids(1) = 4
          vertex_ids(2) = 1
          vertex_ids(3) = 5
          vertex_ids(4) = 8
        case(5)
          vertex_ids(1) = 1
          vertex_ids(2) = 4
          vertex_ids(3) = 3
          vertex_ids(4) = 2
        case(6)
          vertex_ids(1) = 5
          vertex_ids(2) = 6
          vertex_ids(3) = 7
          vertex_ids(4) = 8
      end select
    case default
      option%io_buffer = 'Cell type not recognized'
      call printErrMsg(option)
  end select

end subroutine UGridGetCellFaceVertices

! ************************************************************************** !
!
! UGridComputeCoord: Computes coordinates in x,y,z of unstructured grid cells
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine UGridComputeCoord(unstructured_grid,option, &
                             grid_x,grid_y,grid_z, &
                             x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: icell
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal :: centroid(3)

  do icell = 1, unstructured_grid%num_cells_ghosted
    do ivertex = 1, 8
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,icell) + 1
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    centroid = ComputeCentroid(HEX_TYPE,vertex_8)
    grid_x(icell) = centroid(1)
    grid_y(icell) = centroid(2)
    grid_z(icell) = centroid(3)
  enddo

  do ivertex = 1, unstructured_grid%num_vertices_local
    if (x_max < unstructured_grid%vertices(ivertex)%x) &
      x_max = unstructured_grid%vertices(ivertex)%x
    if (x_min > unstructured_grid%vertices(ivertex)%x) &
      x_min = unstructured_grid%vertices(ivertex)%x
    if (y_max < unstructured_grid%vertices(ivertex)%y) &
      y_max = unstructured_grid%vertices(ivertex)%y
    if (y_min > unstructured_grid%vertices(ivertex)%y) &
      y_min = unstructured_grid%vertices(ivertex)%y
    if (z_max < unstructured_grid%vertices(ivertex)%z) &
      z_max = unstructured_grid%vertices(ivertex)%z
    if (z_min > unstructured_grid%vertices(ivertex)%z) &
      z_min = unstructured_grid%vertices(ivertex)%z
  enddo
      
end subroutine UGridComputeCoord

! ************************************************************************** !
!
! UGridComputeVolumes: Computes volume of unstructured grid cells
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
subroutine UGridComputeVolumes(unstructured_grid,option,nL2G,volume)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  PetscInt :: nL2G(:)
  Vec :: volume
  

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal, pointer :: volume_p(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(volume,volume_p,ierr)

  do local_id = 1, unstructured_grid%num_cells_local
    ghosted_id = nL2G(local_id)
    do ivertex = 1, 8
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,ghosted_id) + 1
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    volume_p(local_id) = ComputeVolume(HEX_TYPE,vertex_8)
  enddo
      
  call VecRestoreArrayF90(volume,volume_p,ierr)

end subroutine UGridComputeVolumes

! ************************************************************************** !
!
! ComputeCentroid: Computes the centroid a grid cell
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
function ComputeCentroid(cell_type,vertices)

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  
  PetscReal :: ComputeCentroid(3)
  PetscInt :: ivertex
  
  ComputeCentroid = 0.d0
  select case(cell_type)
    case(HEX_TYPE)
      ! need something more sophisticated, but for now, just use average
      do ivertex = 1, 8
        ComputeCentroid(1) = ComputeCentroid(1) + vertices(ivertex)%x
        ComputeCentroid(2) = ComputeCentroid(2) + vertices(ivertex)%y
        ComputeCentroid(3) = ComputeCentroid(3) + vertices(ivertex)%z
      enddo
      ComputeCentroid = ComputeCentroid / 8.d0
  end select

end function ComputeCentroid

! ************************************************************************** !
!
! ComputeVolume: Computes the volume a grid cell
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
function ComputeVolume(cell_type,vertices)

  use Utility_module, only : DotProduct

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  
  PetscReal :: ComputeVolume
  PetscReal :: v(3)
  PetscReal :: l1, l2, l3
  
  ComputeVolume = 0.d0
  select case(cell_type)
    case(HEX_TYPE)
     ! need something more sophisticated, but for now, just use volume of 
     ! box
      v(1) = vertices(2)%x-vertices(1)%x
      v(2) = vertices(2)%y-vertices(1)%y
      v(3) = vertices(2)%z-vertices(1)%z
      l1 = DotProduct(v,v)
      v(1) = vertices(4)%x-vertices(1)%x
      v(2) = vertices(4)%y-vertices(1)%y
      v(3) = vertices(4)%z-vertices(1)%z
      l2 = DotProduct(v,v)
      v(1) = vertices(5)%x-vertices(1)%x
      v(2) = vertices(5)%y-vertices(1)%y
      v(3) = vertices(5)%z-vertices(1)%z
      l3 = DotProduct(v,v)
      ComputeVolume = l1*l2*l3
  end select

end function ComputeVolume

! ************************************************************************** !
!
! ComputePlane: Computes the plane intersected by 3 points
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine ComputePlane(plane,point1,point2,point3)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point1, point2, point3
  
  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: x3,y3,z3
  x1 = point1%x
  y1 = point1%y
  z1 = point1%z
  x2 = point2%x
  y2 = point2%y
  z2 = point2%z
  x3 = point3%x
  y3 = point3%y
  z3 = point3%z
  
  ! this grabbed from python script
  plane%A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
  plane%B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2)
  plane%C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
  plane%D = -1.*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1))

end subroutine ComputePlane

! ************************************************************************** !
!
! GetPlaneIntercept: Computes the intercept of a line with a plane
! author: Glenn Hammond
! date: 10/30/09
!
! ************************************************************************** !
subroutine GetPlaneIntercept(plane,point1,point2,intercept)

  implicit none
  
  type(plane_type) :: plane
  type(point_type) :: point1, point2
  type(point_type) :: intercept

  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: u
    
  x1 = point1%x
  y1 = point1%y
  z1 = point1%z
  x2 = point2%x
  y2 = point2%y
  z2 = point2%z
 
 
  u = (plane%A*x1 + plane%B*y1 + plane%C*z1 + plane%D) / &
      (plane%A*(x1-x2) + plane%B*(y1-y2) + plane%C*(z1-z2))

  intercept%x = point1%x + u*(point2%x-point1%x)
  intercept%y = point1%y + u*(point2%y-point1%y)
  intercept%z = point1%z + u*(point2%z-point1%z)

end subroutine GetPlaneIntercept

! ************************************************************************** !
!
! UGDMCreateJacobian: Creates a Jacobian matrix based on the unstructured
!                     grid dual
! author: Glenn Hammond
! date: 11/05/09
!
! ************************************************************************** !
subroutine UGDMCreateJacobian(unstructured_grid,ugdm,mat_type,J,option)

  use Option_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  MatType :: mat_type
  Mat :: J
  type(option_type) :: option

  PetscInt, allocatable :: d_nnz(:), o_nnz(:)
  PetscInt :: icell, ineighbor, neighbor_id
  PetscInt :: ndof_local
  PetscErrorCode :: ierr
  
  allocate(d_nnz(unstructured_grid%num_cells_local))
  allocate(o_nnz(unstructured_grid%num_cells_local))
  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0
  do icell = 1, unstructured_grid%num_cells_local
    do ineighbor = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      neighbor_id = unstructured_grid%cell_neighbors_local_ghosted(ineighbor,icell)
      if (neighbor_id > 0) then
        d_nnz(icell) = d_nnz(icell) + 1
      else
        o_nnz(icell) = o_nnz(icell) + 1
      endif
    enddo
  enddo

  ndof_local = unstructured_grid%num_cells_local*ugdm%ndof
  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*ugdm%ndof
        o_nnz = o_nnz*ugdm%ndof
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltogb,ierr)
      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltogb,ierr)
      case default
        option%io_buffer = 'MatType not recognized in UGDMCreateJacobian'
        call printErrMsg(option)
    end select
  else
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*ugdm%ndof
        call MatCreateSeqAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
      case(MATBAIJ)
        call MatCreateSeqBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
                             PETSC_NULL_INTEGER,d_nnz,J,ierr)
      case default
        option%io_buffer = 'MatType not recognized in UGDMCreateJacobian'
        call printErrMsg(option)
    end select
  endif

  deallocate(d_nnz)
  deallocate(o_nnz)
  
end subroutine UGDMCreateJacobian

! ************************************************************************** !
!
! UGDMCreateVector: Creates a global vector with PETSc ordering
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
subroutine UGDMCreateVector(unstructured_grid,ugdm,vec,vec_type,option)

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
      call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local* &
                        ugdm%ndof, &
                        PETSC_DETERMINE,vec,ierr)
      call VecSetLocalToGlobalMapping(vec,ugdm%mapping_ltog,ierr)
      call VecSetLocalToGlobalMappingBlock(vec,ugdm%mapping_ltogb,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
    case(LOCAL)
      call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%num_cells_ghosted* &
                        ugdm%ndof, &
                        vec,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
    case(NATURAL)
      call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local* &
                        ugdm%ndof, &
                        PETSC_DETERMINE,vec,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
  end select
    
end subroutine UGDMCreateVector

! ************************************************************************** !
!
! UGridMapIndices: maps global, local and natural indices of cells to each other
! author: Glenn Hammond
! date: 11/06/09
!
! ************************************************************************** !
subroutine UGridMapIndices(unstructured_grid,ugdm,nG2L,nL2G,nL2A,nG2A)

  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(ugdm_type) :: ugdm
  PetscInt, pointer :: nG2L(:)
  PetscInt, pointer :: nL2G(:)
  PetscInt, pointer :: nL2A(:)
  PetscInt, pointer :: nG2A(:)
  PetscErrorCode :: ierr
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id

  allocate(nG2L(unstructured_grid%num_cells_ghosted))
  allocate(nL2G(unstructured_grid%num_cells_local))
  allocate(nL2A(unstructured_grid%num_cells_local))
  allocate(nG2A(unstructured_grid%num_cells_ghosted))
  
  ! initialize ghosted to 0
  nG2L = 0

  call ISGetIndicesF90(ugdm%is_local_petsc,int_ptr,ierr)
  do local_id = 1, unstructured_grid%num_cells_local
    nL2G(local_id) = local_id
    nG2L(local_id) = local_id
    ! actually, nL2A is zero-based
    !nL2A(local_id) = int_ptr(local_id)+1
    nL2A(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(ugdm%is_local_petsc,int_ptr,ierr)
!zero-based  nL2A = nL2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%num_cells_local, &
                            nL2A,ierr)
!zero-based  nL2A = nL2A + 1

  call ISGetIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    nG2A(ghosted_id) = int_ptr(ghosted_id)+1
  enddo
  call ISRestoreIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  nG2A = nG2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%num_cells_ghosted, &
                            nG2A,ierr)
  nG2A = nG2A + 1

end subroutine UGridMapIndices

! ************************************************************************** !
!
! UnstructuredGridDestroy: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UnstructuredGridDestroy(unstructured_grid)

  implicit none
  
  type(unstructured_grid_type), pointer :: unstructured_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(unstructured_grid)) return
  
  if (associated(unstructured_grid%cell_vertices_0)) &
    deallocate(unstructured_grid%cell_vertices_0)
  nullify(unstructured_grid%cell_vertices_0)
  if (associated(unstructured_grid%vertices)) &
    deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)
  if (associated(unstructured_grid%hash)) &
    deallocate(unstructured_grid%hash)
  nullify(unstructured_grid%hash)
  if (associated(unstructured_grid%cell_ids_petsc)) &
    deallocate(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%cell_ids_petsc)
  if (associated(unstructured_grid%ghost_cell_ids_petsc)) &
    deallocate(unstructured_grid%ghost_cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  if (associated(unstructured_grid%cell_ids_natural)) &
    deallocate(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_natural)
  if (associated(unstructured_grid%ghost_cell_ids_natural)) &
    deallocate(unstructured_grid%ghost_cell_ids_natural)
  nullify(unstructured_grid%ghost_cell_ids_natural)
  if (associated(unstructured_grid%cell_neighbors_local_ghosted)) &
    deallocate(unstructured_grid%cell_neighbors_local_ghosted)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)

  if (unstructured_grid%ao_natural_to_petsc /= 0) &
    call AODestroy(unstructured_grid%ao_natural_to_petsc,ierr)

  deallocate(unstructured_grid)
  nullify(unstructured_grid)

end subroutine UnstructuredGridDestroy

! ************************************************************************** !
!
! UGDMDestroy: Deallocates a unstructured grid distributed mesh
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGDMDestroy(ugdm)

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
  call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltog,ierr)
  if (ugdm%mapping_ltogb /= 0) &
    call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltogb,ierr)
  call VecDestroy(ugdm%global_vec,ierr)
  call VecDestroy(ugdm%local_vec,ierr)
  deallocate(ugdm)
  nullify(ugdm)

end subroutine UGDMDestroy

end module Unstructured_Grid_module
