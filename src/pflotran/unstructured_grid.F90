module Unstructured_Grid_module

  use Connection_module
  
  implicit none

  private 
  
#include "definitions.h"

  type, public :: unstructured_grid_type
    PetscInt :: num_cells
    PetscInt :: num_vertices
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

  unstructured_grid%num_cells = 0
  unstructured_grid%num_vertices = 0
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

  integer :: icell, ivertex, idir
  
  input => InputCreate(86,filename)

  card = 'Unstructured Grid'

  call InputReadFlotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,card)  

  call InputReadInt(input,option,unstructured_grid%num_cells)
  call InputErrorMsg(input,option,'number of cells',card)
  call InputReadInt(input,option,unstructured_grid%num_vertices)
  call InputErrorMsg(input,option,'number of vertices',card)

  allocate(unstructured_grid%cell_vertices(8,unstructured_grid%num_cells))
  unstructured_grid%cell_vertices = 0

  do icell = 1, unstructured_grid%num_cells
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)  
    do ivertex = 1, 8
      call InputReadInt(input,option, &
                        unstructured_grid%cell_vertices(ivertex,icell))
      call InputErrorMsg(input,option,'vertex id',card)
    enddo
  enddo

  allocate(unstructured_grid%vertex_coordinates(3,unstructured_grid%num_vertices))
  unstructured_grid%vertex_coordinates = 0.d0

  do ivertex = 1, unstructured_grid%num_vertices
    call InputReadFlotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)  
    do idir = 1, 3
      call InputReadDouble(input,option, &
                           unstructured_grid%vertex_coordinates(idir,ivertex))
      call InputErrorMsg(input,option,'vertex coordinate',card)
    enddo
  enddo

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
#include "finclude/petscviewer.h"

  interface
  end interface
  
  type(unstructured_grid_type) :: unstructured_grid
  DM :: dm
  PetscInt :: ndof
  type(option_type) :: option
  
#ifdef ENABLE_UNSTRUCTURED  
  PetscInt, allocatable :: cell_distribution(:)
  PetscInt :: icell, ivertex, count
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscInt, allocatable :: cell_counts(:)
  PetscInt, pointer :: index_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: strided_indices(:)
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
  
  ! cell distribution across processors (size = num_cores + 1)
  ! core i owns cells cell_distribution(i):cell_distribution(i+1), note
  ! the zero-based ordering
  allocate(cell_distribution(option%mycommsize+1))
  cell_distribution(1) = 0
  cell_distribution(2:) = unstructured_grid%num_cells
  num_local_cells = cell_distribution(option%myrank+1)- &
                    cell_distribution(option%myrank+2))
  
  ! provide ids of vertices for local cells
  allocate(local_vertices(8*unstructured_grid%num_cells))
  allocate(local_vertex_offset(unstructured_grid%num_cells+1))
  count = 0
  local_vertex_offset(1) = 0
  do icell = 1, unstructured_grid%num_cells
    do ivertex = 1, 8
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices(ivertex,icell)
    enddo
    local_vertex_offset(icell+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
  num_common_vertices = 3 ! cells must share at least vertices
  
  call MatCreateMPIAdj(option%mycomm,unstructured_grid%num_cells, &
                       unstructured_grid%num_vertices,local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat,ierr)

#if 0
  call PetscViewerASCIIOpen(option%mycomm,'Adj.out',viewer,ierr)
  call MatView(Adj_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call MatMeshToVertexGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)

#if 0
  call PetscViewerASCIIOpen(option%mycomm,'Dual.out',viewer,ierr)
  call MatView(Dual_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call MatPartitioningCreate(option%mycomm,Part,ierr)
  call MatPartitioningSetAdjacency(Part,Adj_mat,ierr)
  call MatPartitioningSetFromOptions(Part,ierr)
  call MatPartitioningAppy(Part,is_new,ierr)
  call MatPartitioningDestroy(Part,ierr)
  
  call PetscViewerASCIIOpen(option%mycomm,'is.out',viewer,ierr)
  call ISView(is_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  

  allocate(cell_counts(option%mycommsize))
  call ISPartitioningCount(is_new,option%mysize,cell_counts,ierr)
  
  call VecCreate(option%mycomm,elements_new,ierr)
  ! assuming hexahedron (8)
  call VecSetSizes(elements_new,8*cell_counts(option%myrank+1),PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_new)
  
  call ISPartitioningToNumbering(is_new,is_num,ierr)
  call ISDestroy(is_new,ierr)
  call ISGetIndicesF90(is_num,index_ptr,ierr)
  allocate(strided_indices(num_local_cells)
  do i=1, num_local_cells
    strided_indices(i) = 8*index_ptr(i)
  enddo
  call ISCreateBlock(option%mycomm,8,num_local_cells,strided_indices,is_scatter,ierr)
  call ISRestoreIndicesF90(is_num,index_ptr,ierr)
  deallocate(strided_indices)
  call ISDestroy(is_num,ierr)
  
  call VecCreateSeq(option%mycomm,8*num_local_cells,elements_old,ierr)
  call VecGetArrayF90(elements_old,vec_ptr,ierr)
  do i=1, 8*num_local_cells
    vec_ptr(i) = local_vertices(i)
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call VecScatterCreate(elements_old,PETSC_NULL,elements_new,is_scatter,vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_new,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_new,INSERT_VALUES,SCATTER_FORWARD,ierr)
  

#if 0
!  allocate(xadj_ptr)
!  nullify(xadj_ptr%f90iptr)
!  allocate(xadj_ptr%f90iptr(3))
!  xadj_ptr%f90iptr(1) = -1
!  xadj_ptr%f90iptr(2) = -2
!  xadj_ptr%f90iptr(3) = -3
!  call assign_c_int_ptr(xadj_cptr,xadj_ptr)
  call assign_c_int_ptr(xadj_cptr,xadj)

!  allocate(adjncy_ptr)
!  nullify(adjncy_ptr%f90iptr)
!  allocate(adjncy_ptr%f90iptr(3))
!  adjncy_ptr%f90iptr(1) = -10
!  adjncy_ptr%f90iptr(2) = -20
!  adjncy_ptr%f90iptr(3) = -30
!  call assign_c_int_ptr(adjncy_cptr,adjncy_ptr)
  call assign_c_int_ptr(adjncy_cptr,adjncy)

  call ParMETIS_V3_Mesh2Dual(cell_distribution,local_vertex_offset, &
                             local_vertices,index_format_flag, &
                             num_common_vertices,xadj_cptr,adjncy_cptr, &
                             option%mycomm)
!  call ParMETIS_V3_Mesh2Dual(idxtype *elmdist, idxtype *eptr, idxtype *eind, int *numflag,
!int *ncommonnodes, idxtype **xadj, idxtype **adjncy, MPI Comm *comm)

!  allocate(xadj_ptr)
!  nullify(xadj_ptr%f90ptr)
!  call assign_c_int_ptr(xadj_cptr,xadj_ptr)
  xadj => xadj_ptr%f90iptr
  deallocate(xadj_ptr)
  nullify(xadj_ptr)
!  
!  allocate(adjncy_ptr)
!  nullify(adjncy_ptr%f90ptr
!  call assign_c_int_ptr(adjncy_cptr,adjncy_ptr)
  adjncy => adjncy_ptr%f90iptr
  deallocate(adjncy_ptr)
  nullify(adjncy_ptr)
#endif
  
  deallocate(cell_distribution)
  deallocate(local_vertices)
  deallocate(local_vertex_offset)
  
  
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
