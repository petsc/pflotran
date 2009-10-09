module Unstructured_Grid_module

  use Connection_module
  
  implicit none

  private 
  
#include "definitions.h"

  type, public :: unstructured_grid_type
    PetscInt :: num_cells_global, num_cells_local
    PetscInt :: num_vertices_global, num_vertices_local
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash
    PetscInt, pointer :: cell_vertices(:,:)
    PetscReal, pointer :: vertex_coordinates(:,:)
  end type unstructured_grid_type

  type, public :: f90iptrwrap
    PetscInt, pointer, dimension(:) :: f90iptr
  end type f90iptrwrap
  
  public :: UnstructuredGridCreate, &
            UnstructuredGridCreateDM, &
            UnstGridComputeInternConnect, &
            UnstGridComputeBoundConnect, &
            UnstructGridGetGhostIdFromHash, &
            UnstructuredGridRead, &
            UnstructuredGridDecompose, &
            UnstructuredGridDestroy

contains

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
  unstructured_grid%nmax = 0
  unstructured_grid%nlmax = 0
  unstructured_grid%ngmax = 0
  nullify(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%vertex_coordinates)
  nullify(unstructured_grid%hash)
  unstructured_grid%num_hash = 100

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

  PetscInt :: icell, ivertex, idir, irank
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscInt :: status(MPI_STATUS_SIZE)
  
  input => InputCreate(86,filename)

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

  allocate(unstructured_grid%cell_vertices(8,unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices = 0

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
        unstructured_grid%cell_vertices(:,1:unstructured_grid%num_cells_local) = &
          temp_int_array(:,1:unstructured_grid%num_cells_local)
      else
print *, '0: ', num_to_read, ' cells sent'
        call MPI_Send(temp_int_array,num_to_read*8,MPI_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
print *, option%myrank,': ',unstructured_grid%num_cells_local, ' cells recv'
    call MPI_Recv(unstructured_grid%cell_vertices, &
                  unstructured_grid%num_cells_local*8,MPI_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status,ierr)
  endif


  unstructured_grid%num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = unstructured_grid%num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              unstructured_grid%num_vertices_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_vertices_local = &
                                 unstructured_grid%num_vertices_local + 1

  allocate(unstructured_grid%vertex_coordinates(3,unstructured_grid%num_vertices_local))
  unstructured_grid%vertex_coordinates = 0.d0

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
        unstructured_grid%vertex_coordinates(:,1:unstructured_grid%num_cells_local) = &
          temp_real_array(:,1:unstructured_grid%num_cells_local)
      else
        call MPI_Send(temp_real_array,num_to_read,MPI_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    call MPI_Recv(unstructured_grid%vertex_coordinates, &
                  unstructured_grid%num_vertices_local*3, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status,ierr)
  endif

  call InputDestroy(input)

end subroutine UnstructuredGridRead

! ************************************************************************** !
!
! UnstructuredGridDecompose: Decomposes an unstructured grid
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UnstructuredGridDecompose(unstructured_grid,dm,ndof,option)
  
  use Option_module
  
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
  DM :: dm
  PetscInt :: ndof
  type(option_type) :: option
  
#ifdef ENABLE_UNSTRUCTURED  
!  PetscInt, allocatable :: cell_distribution(:)
  PetscInt :: icell, ivertex, count, vertex_count
  PetscInt :: max_vertex_count
  PetscInt :: max_dual
  PetscInt :: stride
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscInt, allocatable :: cell_counts(:)
  PetscInt, pointer :: index_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: strided_indices(:)
  PetscInt, pointer :: ia_ptr(:), ja_ptr(:)
  PetscInt :: num_rows, num_cols, istart, iend, icol
  PetscTruth :: success
  PetscErrorCode :: ierr
  
  PetscViewer :: viewer
  Mat :: Adj_mat
  Mat :: Dual_mat
  MatPartitioning :: Part
  Vec :: elements_new
  Vec :: elements_old
  IS :: is_new
  IS :: is_num
  IS :: is_scatter
  VecScatter :: vec_scatter

  PetscInt :: global_offset
  
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

  max_vertex_count = 8
  max_dual = 6
  stride = max_vertex_count + max_dual + 1 + 1 ! another for -999

  allocate(local_vertices(max_vertex_count*unstructured_grid%num_cells_local))
  allocate(local_vertex_offset(unstructured_grid%num_cells_local+1))
  count = 0
  local_vertex_offset(1) = 0
  do icell = 1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices(ivertex,icell)
    enddo
    local_vertex_offset(icell+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
  num_common_vertices = 3 ! cells must share at least vertices

  global_offset = 0
  call MPI_Exscan(unstructured_grid%num_cells_local,global_offset,ONE_INTEGER,&
                  MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
  
  call MatCreateMPIAdj(option%mycomm,unstructured_grid%num_cells_local, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat,ierr)

#if 1
  call PetscViewerASCIIOpen(option%mycomm,'Adj.out',viewer,ierr)
  call MatView(Adj_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)

#if 1
  call PetscViewerASCIIOpen(option%mycomm,'Dual.out',viewer,ierr)
  call MatView(Dual_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call MatPartitioningCreate(option%mycomm,Part,ierr)
  call MatPartitioningSetAdjacency(Part,Dual_mat,ierr)
  call MatPartitioningSetFromOptions(Part,ierr)
  call MatPartitioningApply(Part,is_new,ierr)
  call MatPartitioningDestroy(Part,ierr)
  
  call PetscViewerASCIIOpen(option%mycomm,'is.out',viewer,ierr)
  call ISView(is_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  

  allocate(cell_counts(option%mycommsize))
  call ISPartitioningCount(is_new,option%mycommsize,cell_counts,ierr)
  
  call VecCreate(option%mycomm,elements_new,ierr)
  call VecSetSizes(elements_new, &
                   stride*cell_counts(option%myrank+1), &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_new,ierr)
  
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
    vec_ptr(count) = -(global_offset+icell)
    do ivertex = 1, max_vertex_count
      count = count + 1
      vertex_count = vertex_count + 1
        ! increment for 1-based ordering
      vec_ptr(count) = local_vertices(vertex_count) + 1
    enddo

    count = count + 1 
    vec_ptr(count) = -999  ! help differentiate

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
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call VecScatterCreate(elements_old,PETSC_NULL,elements_new,is_scatter,vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_new,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_new,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

  call PetscViewerASCIIOpen(option%mycomm,'elements_old.out',viewer,ierr)
  call VecView(elements_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  
  call VecDestroy(elements_old,ierr)
  
  call PetscViewerASCIIOpen(option%mycomm,'elements_new.out',viewer,ierr)
  call VecView(elements_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  

!  deallocate(cell_distribution)
  deallocate(local_vertices)
  deallocate(local_vertex_offset)
  
  stop
#endif
  
end subroutine UnstructuredGridDecompose

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
  if (associated(unstructured_grid%cell_vertices)) &
    deallocate(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%cell_vertices)
  if (associated(unstructured_grid%vertex_coordinates)) &
    deallocate(unstructured_grid%vertex_coordinates)
  nullify(unstructured_grid%vertex_coordinates)

  deallocate(unstructured_grid)
  nullify(unstructured_grid)

end subroutine UnstructuredGridDestroy

end module Unstructured_Grid_module
