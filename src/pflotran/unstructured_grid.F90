module Unstructured_Grid_module

  use Connection_module
  
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
    PetscInt, pointer :: cell_vertices_nindex(:,:) ! vertices for each grid cell (natural index 0-based)
    PetscInt, pointer :: face_to_cell_locindex(:,:) !
    PetscInt, pointer :: face_to_vertex_nindex(:,:) !
    PetscInt, pointer :: cell_to_face_locindex(:,:)
    PetscInt, pointer :: vertex_ids_nindex(:)
    PetscInt, pointer :: cell_ids_natural(:) ! natural 1d right-hand i,j,k ordering
    PetscInt, pointer :: cell_ids_petsc(:) ! petsc ordering of cell ids
    PetscInt, pointer :: ghost_cell_ids_natural(:) ! natural ordering of ghost cell ids
    PetscInt, pointer :: ghost_cell_ids_petsc(:) ! petsc ordering of ghost cells ids
    PetscInt, pointer :: cell_neighbors_local_ghosted(:,:) ! local neighbors
    type(point_type), pointer :: vertices(:)
    type(point_type), pointer :: face_centroid(:)
    PetscReal, pointer :: face_area(:)
    AO :: ao_natural_to_petsc ! mapping of natural to petsc ordering
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

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: WEDGE_TYPE        = 2
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_CELL = 8
  !  PetscInt, parameter :: MAX_DUALS         = 6
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4
  !  PetscInt, parameter :: MAX_CELLS_SHARING_A_VERTEX = 16

  public :: UnstructuredGridCreate, &
            UnstructuredGridRead, &
#ifndef SAMR_HAVE_HDF5
            UnstructuredGridReadHDF5, &
#endif
#if defined(PARALLELIO_LIB)
            UnstructuredGridReadHDF5PIOLib, &
#endif
            UnstructuredGridDecompose, &
            UGridComputeInternConnect, &
            UGridPopulateConnection, &
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
  nullify(unstructured_grid%cell_vertices_nindex)
  nullify(unstructured_grid%face_to_cell_locindex)
  nullify(unstructured_grid%face_to_vertex_nindex)
  nullify(unstructured_grid%cell_to_face_locindex)
  nullify(unstructured_grid%vertex_ids_nindex)
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

  PetscInt :: icell, ivertex, idir, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename)

! Format of unstructured grid file
! Currently assumes hexahedron
! -----------------------------------------------------------------
! num_cells num_vertices  (integers)
! vert1 vert2 vert3 ... vert8  ! for cell 1 (integers)
! vert1 vert2 vert3 ... vert8  ! for cell 2
! ...
! ...
! vert1 vert2 vert3 ... vert8  ! for cell num_cells
! xcoord ycoord zcoord ! coordinates of vertex 1 (real)
! xcoord ycoord zcoord ! coordinates of vertex 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
! -----------------------------------------------------------------

!
! Note: When compiled with -DMIXED_FLAG the input file needs to be
!       modified. Each line corresponding to a cell now needs to contain
!       number of vertices as the first entry. Presnetly (4/25/2011)
!       only hexahedron and wedge cell-types supported
!
! -----------------------------------------------------------------
! num_cells num_vertices  (integers)
! nverts vert1 vert2 vert3 ... vert8  ! for cell 1 (integers)
! nverts vert1 vert2 vert3 ... vert8  ! for cell 2
! ...
! ...
! nverts vert1 vert2 vert3 ... vert8  ! for cell num_cells
! xcoord ycoord zcoord ! coordinates of vertex 1 (real)
! xcoord ycoord zcoord ! coordinates of vertex 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
! -----------------------------------------------------------------
!

  card = 'Unstructured Grid'

  call InputReadFlotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,card)  

  ! read num_cells
  call InputReadInt(input,option,unstructured_grid%num_cells_global)
  call InputErrorMsg(input,option,'number of cells',card)
  ! read num_vertices
  call InputReadInt(input,option,unstructured_grid%num_vertices_global)
  call InputErrorMsg(input,option,'number of vertices',card)

  ! divide cells across ranks
  unstructured_grid%num_cells_local = unstructured_grid%num_cells_global/ &
                                      option%mycommsize 
  num_cells_local_save = unstructured_grid%num_cells_local
  remainder = unstructured_grid%num_cells_global - &
              unstructured_grid%num_cells_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_cells_local = &
                                 unstructured_grid%num_cells_local + 1

  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL,unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices_0 = 0

  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(MAX_VERT_PER_CELL,num_cells_local_save+1))
    temp_int_array = -1
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do icell = 1, num_to_read
        ! read in the vertices defining the grid cell
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        num_vertices = MAX_VERT_PER_CELL
#ifdef MIXED_UMESH
        call InputReadInt(input,option,num_vertices)
        call InputErrorMsg(input,option,'num_vertices',card)
        if (num_vertices.gt.MAX_VERT_PER_CELL) then
          option%io_buffer = 'Cells verticies exceed maximum number of vertices'
          call printErrMsg(option)
        endif
        if((num_vertices.ne.8).and.(num_vertices.ne.6)) then
          write(*,*),num_vertices
          option%io_buffer = 'Only cells with 6 or 8 vertices supported'
          call printErrMsg(option)
        endif
#endif
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,icell))
          call InputErrorMsg(input,option,'vertex id',card)
        enddo
      enddo
      
      ! if the cells reside on io_rank
      if (irank == option%io_rank) then
print *, '0: ', unstructured_grid%num_cells_local, ' cells'
        unstructured_grid%cell_vertices_0(:,1:unstructured_grid%num_cells_local) = &
          temp_int_array(:,1:unstructured_grid%num_cells_local)
      else
        ! otherwise communicate to other ranks
print *, '0: ', num_to_read, ' cells sent'
        int_mpi = num_to_read*MAX_VERT_PER_CELL
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
print *, option%myrank,': ',unstructured_grid%num_cells_local, ' cells recv'
    int_mpi = unstructured_grid%num_cells_local*MAX_VERT_PER_CELL
    call MPI_Recv(unstructured_grid%cell_vertices_0,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif


  ! divide vertices across ranks
  unstructured_grid%num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = unstructured_grid%num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              unstructured_grid%num_vertices_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_vertices_local = &
                                 unstructured_grid%num_vertices_local + 1

  allocate(vertex_coordinates(3,unstructured_grid%num_vertices_local))
  vertex_coordinates = 0.d0

  ! just like above, but this time for vertex coordinates
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
  
  ! fill the vertices data structure
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

#ifndef SAMR_HAVE_HDF5

! ************************************************************************** !
!
! UnstructuredGridReadHDF5: Reads an unstructured grid from HDF5
! author: Gautam Bisht
! date: 04/25/11
!
! ************************************************************************** !
subroutine UnstructuredGridReadHDF5(unstructured_grid,filename,option)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use Input_module
  use Option_module

  implicit none

  type(unstructured_grid_type)   :: unstructured_grid
  type(option_type)              :: option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscMPIInt       :: hdf5_err
  PetscMPIInt       :: rank_mpi
  PetscInt          :: ndims
  PetscInt          :: istart, iend, ii, jj
  PetscInt          :: num_cells_local_save
  PetscInt          :: num_vertices_local_save
  PetscInt          :: remainder
  PetscInt,pointer  :: int_buffer(:,:)
  PetscReal,pointer :: double_buffer(:,:)
  PetscErrorCode    :: ierr

#if defined(PETSC_HAVE_HDF5)
  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id, grp_id2
  integer(HID_T) :: prop_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: data_space_id
  integer(HID_T) :: memory_space_id
  integer(HSIZE_T) :: num_data_in_file
  integer(HSIZE_T), allocatable :: dims_h5(:), max_dims_h5(:)
  integer(HSIZE_T) :: offset(2), length(2), stride(2), block(2), dims(2)
#endif

  ! Initialize FOTRAN predefined datatypes
  call h5open_f(hdf5_err)

  ! Setup file access property with parallel I/O access
  call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, hdf5_err)

#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm, MPI_INFO_NULL, hdf5_err)
#endif

  ! Open the file collectively
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdf5_err, prop_id)
  call h5pclose_f(prop_id, hdf5_err)
  
  !
  ! Domain/Cells
  !
  
  ! Open group
  group_name = "Domain"
  option%io_buffer = 'Opening group: ' // trim(group_name)
  call printMsg(option)

  ! Open dataset
  call h5dopen_f(file_id, "Domain/Cells", data_set_id, hdf5_err)

  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims, hdf5_err)
  if (ndims /= 2) then
    option%io_buffer='Dimension of Domain/Cells dataset in ' // filename // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%num_cells_global = dims_h5(2)
  unstructured_grid%num_cells_local = unstructured_grid%num_cells_global/ &
                                      option%mycommsize 
  num_cells_local_save = unstructured_grid%num_cells_local
  remainder = unstructured_grid%num_cells_global - &
              unstructured_grid%num_cells_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_cells_local = &
                                  unstructured_grid%num_cells_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(unstructured_grid%num_cells_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(unstructured_grid%num_cells_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
  ! Determine the length and offset of data to be read by each processor
  length(1) = dims_h5(1)
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart
  
  !
  rank_mpi = 2
  memory_space_id = -1
  
  ! Create data space for dataset
  call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)
  
  ! Select hyperslab
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(int_buffer(length(1), length(2)))
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, int_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL, &
                                            unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices_0 = 0
  
  do ii = 1, unstructured_grid%num_cells_local
    do jj = 2, int_buffer(1,ii) + 1
      unstructured_grid%cell_vertices_0(jj-1, ii) = int_buffer(jj, ii)
    enddo
  enddo
  
  call h5dclose_f(data_set_id, hdf5_err)
  
  deallocate(dims_h5)
  deallocate(max_dims_h5)

  !
  ! Domain/Vertices
  !
  
  ! Open dataset
  call h5dopen_f(file_id, "Domain/Vertices", data_set_id, hdf5_err)
  
  ! Get dataset's dataspace
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  
  ! Get number of dimensions and check
  call h5sget_simple_extent_ndims_f(data_space_id, ndims, hdf5_err)
  if (ndims /= 2) then
    option%io_buffer='Dimension of Domain/Vertices dataset in ' // filename // &
          ' is not equal to 2.'
    call printErrMsg(option)
  endif
  
  ! Allocate memory
  allocate(dims_h5(ndims))
  allocate(max_dims_h5(ndims))
  
  ! Get dimensions of dataset
  call h5sget_simple_extent_dims_f(data_space_id, dims_h5, max_dims_h5, &
                                   hdf5_err)
  
  ! Determine the number of cells each that will be saved on each processor
  unstructured_grid%num_vertices_global = dims_h5(2)
  unstructured_grid%num_vertices_local  = &
                                       unstructured_grid%num_vertices_global/ &
                                       option%mycommsize 
  num_cells_local_save = unstructured_grid%num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              unstructured_grid%num_vertices_local*option%mycommsize
  if (option%myrank < remainder) unstructured_grid%num_vertices_local = &
                                  unstructured_grid%num_vertices_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(unstructured_grid%num_vertices_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(unstructured_grid%num_vertices_local, iend, ONE_INTEGER_MPI, &
                MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  
  ! Determine the length and offset of data to be read by each processor
  length(1) = dims_h5(1)
  length(2) = iend-istart
  offset(1) = 0
  offset(2) = istart
  
  ! 
  rank_mpi = 2
  memory_space_id = -1
  
  ! Create data space for dataset
  call h5screate_simple_f(rank_mpi, length, memory_space_id, hdf5_err)
  
  ! Select hyperslab
  call h5dget_space_f(data_set_id, data_space_id, hdf5_err)
  call h5sselect_hyperslab_f(data_space_id, H5S_SELECT_SET_F, offset, length, &
                             hdf5_err)
  
  ! Initialize data buffer
  allocate(double_buffer(length(1), length(2)))
  
  ! Create property list
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
#endif
  
  ! Read the dataset collectively
  call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, double_buffer, &
                 dims_h5, hdf5_err, memory_space_id, data_space_id)
  
  call h5dclose_f(data_set_id, hdf5_err)
  !call h5gclose_f(grp_id, hdf5_err)
  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)

  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(unstructured_grid%num_vertices_local))
  do ii = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo
  
  
  deallocate(double_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)
  
  
end subroutine UnstructuredGridReadHDF5

#endif

! ************************************************************************** !
!
! UnstructuredGridReadHDF5PIOLib: Reads an unstructured grid from HDF5
! author: Gautam Bisht
! date: 05/13/11
!
! ************************************************************************** !
#if defined(PARALLELIO_LIB)

subroutine UnstructuredGridReadHDF5PIOLib(unstructured_grid, filename, &
                                          option)

#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

!#include "definitions.h"

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use Input_module
  use Option_module
  use HDF5_aux_module

  implicit none

  type(unstructured_grid_type)   :: unstructured_grid
  type(option_type)              :: option
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name

  PetscInt,pointer  :: int_buffer(:,:)
  PetscReal,pointer :: double_buffer(:,:)
  PetscInt          :: ii, jj
  PetscInt          :: dims(2), dataset_dims(2)

  character(len=MAXSTRINGLENGTH) :: cell_dataset_name = &
                                                       '/Domain/Cells'//CHAR(0)
  character(len=MAXSTRINGLENGTH) :: vert_dataset_name = &
                                                    '/Domain/Vertices'//CHAR(0)

  ! Read Domain/Cells
  call HDF5ReadDatasetInteger2D(filename, 
                                cell_dataset_name, &
                                NONUNIFORM_CONTIGUOUS_READ, &
                                option, &
                                int_buffer, &
                                dims, &
                                dataset_dims)

  ! Allocate array to store vertices for each cell
  unstructured_grid%num_cells_local  = dims(2)
  unstructured_grid%num_cells_global = dataset_dims(2)
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL, &
                                            unstructured_grid%num_cells_local))
  unstructured_grid%cell_vertices_0 = 0

  ! Fill the cell data structure
  do ii = 1, unstructured_grid%num_cells_local
    do jj = 2, int_buffer(1, ii) + 1
      unstructured_grid%cell_vertices_0(jj-1, ii) = int_buffer(jj, ii)
    enddo
  enddo

  ! Read Vertices
  call HDF5ReadDatasetReal2D(filename, &
                             vert_dataset_name, &
                             NONUNIFORM_CONTIGUOUS_READ, &
                             option, &
                             double_buffer, &
                             dims, &
                             dataset_dims)

  unstructured_grid%num_vertices_local = dims(2)
  unstructured_grid%num_vertices_global= dataset_dims(2)
  allocate(unstructured_grid%vertices(unstructured_grid%num_vertices_local))

  ! fill the vertices data structure
  do ii = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo

end subroutine UnstructuredGridReadHDF5PIOLib

#endif

! ************************************************************************** !
!
! UnstructuredGridDecompose: Decomposes an unstructured grid across ranks
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
#include "finclude/petscdm.h" 
#include "finclude/petscdm.h90"
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
  PetscBool :: success
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
  PetscInt :: num_cells_local_new  !sp 
  PetscInt :: num_cells_local_old  !sp  
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, pointer :: int_array_pointer(:)
  
  PetscInt :: idual, dual_id
  PetscBool :: found


  
!  cell distribution across processors (size = num_cores + 1)
!  core i owns cells cell_distribution(i):cell_distribution(i+1), note
!  the zero-based indexing
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%num_cells_local,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)
  
  ! provide ids of vertices for local cells

  ! in order to redistributed vertex/cell data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  vertex_ids_offset = 1 + 1 ! +1 for -777
  max_vertex_count = MAX_VERT_PER_CELL
  dual_offset = vertex_ids_offset + max_vertex_count + 1 ! +1 for -888
  max_dual = MAX_DUALS
  stride = dual_offset+ max_dual + 1 ! +1 for -999

  ! Information for each cell is packed in a strided petsc vec
  ! The information is ordered within each stride as follows:
  ! -cell_N   ! global cell id (negative indicates 1-based)
  ! -777      ! separator between cell id and vertex ids for cell_N
  ! vertex1   ! in cell_N
  ! vertex2
  ! ...
  ! vertexN   
  ! -888      ! separator between vertex and dual ids
  ! dual1     ! dual ids between cell_N and others
  ! dual2
  ! ...
  ! dualN     
  ! -999    ! separator indicating end of information for cell_N
  
  ! the purpose of -777, -888, and -999 is to allow one to use cells of 
  ! various geometry.  Currently, the max # vertices = 8 and max # duals = 6.
  ! But this will be generalized in the future.
  
  allocate(local_vertices(max_vertex_count*unstructured_grid%num_cells_local))
  allocate(local_vertex_offset(unstructured_grid%num_cells_local+1))
  count = 0
  local_vertex_offset(1) = 0
  do icell = 1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
#ifdef MIXED_UMESH
         if(unstructured_grid%cell_vertices_0(ivertex, icell) < 0) exit
#endif
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices_0(ivertex,icell)
    enddo
    local_vertex_offset(icell+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
  num_common_vertices = 3 ! cells must share at least this number of vertices

  ! determine the global offset from 0 for cells on this rank
  unstructured_grid%global_offset = 0
  call MPI_Exscan(unstructured_grid%num_cells_local, &
                  unstructured_grid%global_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! create an adjacency matrix for calculating the duals (connnections)
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

  ! petsc will call parmetis to calculate the graph/dual
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

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is.out',viewer,ierr)
  call ISView(is_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)  
#endif

  ! calculate the number of local grid cells on each processor
  allocate(cell_counts(option%mycommsize))
  ! ISPartitioningCount takes a ISPartitioning and determines the number of  
  ! resulting elements on each (partition) process - petsc
  call ISPartitioningCount(is_new,option%mycommsize,cell_counts,ierr)
  !sp  unstructured_grid%num_cells_local = cell_counts(option%myrank+1)
  !sp need to distinguish between old and new cell counts  
  num_cells_local_new=cell_counts(option%myrank+1) !sp 
  num_cells_local_old=unstructured_grid%num_cells_local  !sp 
  deallocate(cell_counts)
  
  ! create a petsc vec to store all the information for each element
  ! based on the stride calculated above.  
  call VecCreate(option%mycomm,elements_natural,ierr)
  call VecSetSizes(elements_natural, &
                   stride*num_cells_local_new, &    !sp 
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_natural,ierr)
   
  ! calculate the global offsets in the new vector for each grid cell
  
  ! ISPartitioningToNumbering takes an ISPartitioning and on each processor 
  ! generates an IS that contains a new global node number for each index 
  ! based on the partitioning. - petsc
  call ISPartitioningToNumbering(is_new,is_num,ierr)
  call ISDestroy(is_new,ierr)
  call ISGetIndicesF90(is_num,index_ptr,ierr)

  
  !sp loading of strided_indices eliminated here 10/22/2010 

  ! Create a mapping of local indices to global strided
  call ISCreateBlock(option%mycomm,stride, &
                     num_cells_local_old, &   !sp 
                     index_ptr,PETSC_COPY_VALUES,is_scatter,ierr)
  call ISRestoreIndicesF90(is_num,index_ptr,ierr)
  call ISDestroy(is_num,ierr)

  
  ! create another strided vector with the old cell/element distribution
  call VecCreate(option%mycomm,elements_old,ierr)
  call VecSetSizes(elements_old, &
                   stride*num_cells_local_old,PETSC_DECIDE,ierr)  !sp 
  call VecSetFromOptions(elements_old,ierr)


  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                      ja_ptr,success,ierr)

  if (.not.success .or. &
       num_rows /= num_cells_local_old) then  !sp 
    print *, option%myrank, num_rows, success, num_cells_local_old   !sp 
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  call VecGetArrayF90(elements_old,vec_ptr,ierr)
  count = 0
  vertex_count = 0
  do icell=1, num_cells_local_old !sp 
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(unstructured_grid%global_offset+icell)
    count = count + 1
    ! add the separator
    vec_ptr(count) = -777  ! help differentiate
    ! add the vertex ids
    do ivertex = 1, max_vertex_count
      count = count + 1
      vertex_count = vertex_count + 1
      ! increment for 1-based ordering
!      vec_ptr(count) = local_vertices(vertex_count) + 1
      vec_ptr(count) = unstructured_grid%cell_vertices_0(ivertex,icell) + 1
    enddo


    count = count + 1 
    ! another vertex/dual separator
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
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
    ! final separator
    vec_ptr(count) = -999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call MatRestoreRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                          ja_ptr,success,ierr)
  call MatDestroy(Dual_mat,ierr)
 
#if GEH_DEBUG
  call printMsg(option,'Before element scatter')
#endif

  ! scatter all the cell data from the old decomposition (as read in in parallel)
  ! to the more parmetis-calculated decomposition
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
  
  
  !sp 
  ! num_cells_local may be different with new partitioning 
  !  also need to update global_offset 
  unstructured_grid%num_cells_local = num_cells_local_new  
  unstructured_grid%global_offset = 0
  call MPI_Exscan(unstructured_grid%num_cells_local, &
                  unstructured_grid%global_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  !sp end 

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

  ! now we unpack the decomposed cell data
 
  
  ! store the natural grid cell id for each local cell as read from the grid file
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    unstructured_grid%cell_ids_natural(icell) = abs(vec_ptr((icell-1)*stride+1))
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

  ! make a list of petsc ids for each local cell (you simply take the global 
  ! offset and add it to the local contiguous cell ids on each processor
  allocate(int_array(unstructured_grid%num_cells_local))
  do icell=1, unstructured_grid%num_cells_local
    int_array(icell) = icell+unstructured_grid%global_offset
  enddo
  
  ! make the arrays zero-based
  int_array = int_array - 1
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural - 1
  ! create an application ordering (mapping of natural to petsc ordering)
  call AOCreateBasic(option%mycomm,unstructured_grid%num_cells_local, &
                     unstructured_grid%cell_ids_natural,int_array, &
                     unstructured_grid%ao_natural_to_petsc,ierr)
  deallocate(int_array)
  ! make cell_ids_natural 1-based again
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural + 1

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'ao.out',viewer,ierr)
  call AOView(unstructured_grid%ao_natural_to_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! The below creates a list of cells ids for the duals and converts them
  ! to petsc ordering   
  
  ! count the number of cells and their duals  
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit ! here we hit the -999 at the end of the stride 
      count = count + 1
    enddo
  enddo     
               
  ! allocate and fill an array with the natural cell and dual ids
  allocate(int_array(count))
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    int_array(count) = unstructured_grid%cell_ids_natural(icell)
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit ! again we hit the -999
      count = count + 1
      int_array(count) = dual_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

#if GEH_DEBUG
  call printMsg(option,'Application ordering')
#endif

  ! convert the dual ids in int_array from natural to petsc numbering
  int_array = int_array - 1             
  call AOApplicationToPetsc(unstructured_grid%ao_natural_to_petsc,count, &
                            int_array,ierr)
  int_array = int_array + 1                

#if GEH_DEBUG
  call printMsg(option,'PETSc-ordered duals')
#endif

  ! load mapped petsc-ordered dual ids back into duplicated vector
  ! exactly the opposite operation of when we loaded the temporary int_array
  ! vector
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  allocate(unstructured_grid%cell_ids_petsc(unstructured_grid%num_cells_local))
  count = 0
  do icell=1, unstructured_grid%num_cells_local
    count = count + 1
    ! extract the petsc id for the cell
    unstructured_grid%cell_ids_petsc(icell) = int_array(count)
    ! store it in the elements_petsc vector too
    vec_ptr((icell-1)*stride+1) = int_array(count)
    do idual = 1, max_dual
      dual_id = vec_ptr2(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      ! store the petsc numbered duals in the vector also
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
  ! allocate a temporarily-sized array
  max_ghost_cell_count = max(unstructured_grid%num_cells_local,100)
  allocate(unstructured_grid%ghost_cell_ids_petsc(max_ghost_cell_count))
  ! loop over all duals and find the off-processor cells on the other
  ! end of a dual
  do icell=1, unstructured_grid%num_cells_local
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      found = PETSC_FALSE
      if (dual_id < 1) exit
      do icell2 = 1, unstructured_grid%num_cells_local
        if (dual_id == unstructured_grid%cell_ids_petsc(icell2)) then
          vec_ptr(idual + dual_offset + (icell-1)*stride) = icell2
          found = PETSC_TRUE
          exit
        endif
      enddo
      ! if not found, add it to the list of ghost cells
      if (.not.found) then
        !sp  but only if it is not already in the list 
        found = PETSC_FALSE 
        do icell2 = 1, ghost_cell_count  
          if (dual_id == unstructured_grid%ghost_cell_ids_petsc(icell2)) then 
            vec_ptr(idual + dual_offset + (icell-1)*stride) = -icell2
            found = PETSC_TRUE
            exit
          end if 
        end do 
        if( .not. found) then ! if here then not on processor and not in the ghost list  
          ghost_cell_count = ghost_cell_count + 1
          ! reallocate the ghost cell array if necessary
          if (ghost_cell_count > max_ghost_cell_count) then
            call reallocateIntArray(unstructured_grid%ghost_cell_ids_petsc, &
                max_ghost_cell_count)
            endif
            unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count) = dual_id
            ! flag the id as negative for ghost cell and set to current count
            vec_ptr(idual + dual_offset + (icell-1)*stride) = -ghost_cell_count
          end if 
        endif
     enddo
  enddo
  
  ! This is just a note to add an algorithm that adds all ghost cells and
  ! then removed duplicates.  This will perform much faster.  Do not remove
  ! this print statemetn
  call printMsg(option,'Glenn: Add code that adds all ghost cells and removed duplicates')
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if GEH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_local_unsorted.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
 

  unstructured_grid%num_ghost_cells = ghost_cell_count
  unstructured_grid%num_cells_ghosted = &
    unstructured_grid%num_cells_local + ghost_cell_count

  ! sort ghost cell ids
  allocate(int_array(ghost_cell_count))
  allocate(int_array2(ghost_cell_count))
  allocate(int_array3(ghost_cell_count))
  do icell = 1, ghost_cell_count
   int_array(icell) = unstructured_grid%ghost_cell_ids_petsc(icell)
   int_array3(icell)=abs(int_array(icell) ) 
   int_array2(icell) = icell
  enddo
  !sp end 

  ! convert to 0-based
  int_array2 = int_array2-1
  call PetscSortIntWithPermutation(ghost_cell_count,int_array3,int_array2,ierr)
  ! convert to 1-based
  int_array2 = int_array2+1


  ! resize ghost cell array down to ghost_cell_count
  deallocate(unstructured_grid%ghost_cell_ids_petsc)
  allocate(unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count))
  ! fill with the sorted ids
  do icell = 1, ghost_cell_count
    unstructured_grid%ghost_cell_ids_petsc(int_array2(icell)) = int_array(icell)
  enddo
    
  ! now, rearrange the ghost cell ids of the dual accordingly
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ! add num_cells_local to get in local numbering
  ! basically, the first num_cells_local are the non-ghosted.  After that
  ! the cells ids are ghosted (at the end of the local vector)
!  int_array2 = int_array2 + unstructured_grid%num_cells_local
  do icell=1, unstructured_grid%num_cells_local
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (icell-1)*stride)
      if (dual_id < 0) then
        vec_ptr(idual + dual_offset + (icell-1)*stride) = int_array2(-dual_id)  & 
          + unstructured_grid%num_cells_local
      endif       
    enddo
  enddo
  deallocate(int_array)
  deallocate(int_array2)
  deallocate(int_array3)
!  deallocate(int_array4)
        
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_local.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! load cell neighbors into array
  ! start first index at zero to store # duals for a cell
  ! bad, bad, bad - hardwired to 6!!!  Should be max_num_duals_per_cell
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
      ! flag ghosted cells in dual as negative
!sp      if (dual_id > unstructured_grid%num_cells_local) dual_id = -dual_id
      unstructured_grid%cell_neighbors_local_ghosted(idual,icell) = dual_id
    enddo
    ! set the # of duals in for the cell
    unstructured_grid%cell_neighbors_local_ghosted(0,icell) = count
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! make a list of local vertices
  max_int_count = unstructured_grid%num_cells_local
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
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
  count = 1
  !
  ! The following line commented out by GB, because it was assuming
  ! that int_array2(1) = 1; which may not be the case
  !
  !int_array4(1) = 1
  int_array4(int_array2(1)) = count
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

  allocate(unstructured_grid%vertex_ids_nindex(vertex_count))
  unstructured_grid%vertex_ids_nindex = int_array3(1:vertex_count)
  !do ivertex=1,vertex_count
  !  if(option%myrank.eq.0) write(*,*), 'int_array3(',ivertex,')=',int_array3(ivertex)
  !enddo

  ! now load all the vertices needed to define all the local cells
  ! on the processor
  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  ! allocate the array that will store the vertex ids for each cell.
  ! remember that max_vertex_count is the max # of vertices in a cell
  ! currently hardwired to 8.
  deallocate(unstructured_grid%cell_vertices_0)
  allocate(unstructured_grid%cell_vertices_0(0:max_vertex_count, &
                                             unstructured_grid%num_cells_local))
  allocate(unstructured_grid%cell_vertices_nindex(1:max_vertex_count, &
                                                  unstructured_grid%num_cells_local))
  ! initialize to -999 (for error checking later)
  unstructured_grid%cell_vertices_0 = -999
  ! initialized the 0 entry (which stores the # of vertices in each cell) to zero
  unstructured_grid%cell_vertices_0(0,:) = 0
  
  ! permute the local ids calculated earlier in the int_array4
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do icell=1, unstructured_grid%num_cells_local
    do ivertex = 1, max_vertex_count
      ! extract the original vertex id
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride)
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices_0(0,icell)+1
      unstructured_grid%cell_vertices_0(count,icell) = int_array4(vertex_id)-1
      unstructured_grid%cell_vertices_0(0,icell) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (icell-1)*stride) = int_array4(vertex_id)
      unstructured_grid%cell_vertices_nindex(count,icell) = &
          int_array3(unstructured_grid%cell_vertices_0(count,icell) + 1) -1
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

  ! now we need to work on aligning the original vertex coordinates with 
  ! the current ordering or permuted/rearranged ordering.

  ! IS for gather operation - need local numbering
  allocate(strided_indices(vertex_count))
  ! vertex_count = # of local vertices (I believe ghosted+non-ghosted)
  !sp 22/10/2010 new ISCreateBlock  not strided 
  do ivertex = 1, vertex_count
    ! *3 for 3 coordinates x,y,z 
!    strided_indices(ivertex) = 3*(ivertex-1)
    strided_indices(ivertex) = ivertex-1
  enddo
  deallocate(int_array3)
  !sp end 
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     strided_indices,PETSC_COPY_VALUES,is_gather,ierr)
  deallocate(strided_indices)

  ! create a parallel petsc vector with a stride of 3.
  call VecCreateMPI(option%mycomm,unstructured_grid%num_vertices_local*3, &
                    PETSC_DETERMINE,vertices_old,ierr)
  call VecSetBlockSize(vertices_old,3,ierr)

  ! create serial petsc vector with a stride of 3
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
! load up the coordinates
  call VecGetArrayF90(vertices_old,vec_ptr,ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    vec_ptr((ivertex-1)*3+1) = unstructured_grid%vertices(ivertex)%x
    vec_ptr((ivertex-1)*3+2) = unstructured_grid%vertices(ivertex)%y
    vec_ptr((ivertex-1)*3+3) = unstructured_grid%vertices(ivertex)%z
  enddo
  call VecRestoreArrayF90(vertices_old,vec_ptr,ierr)
  deallocate(unstructured_grid%vertices)
  nullify(unstructured_grid%vertices)

  ! geh - get back to commenting here
  ! IS for scatter - provide petsc global numbering
  allocate(strided_indices(vertex_count))
  do ivertex = 1, vertex_count
!sp 22/10/2010 
!    strided_indices(ivertex) = 3*(needed_vertices_petsc(ivertex)-1)
    strided_indices(ivertex) = (needed_vertices_petsc(ivertex)-1)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     strided_indices,PETSC_COPY_VALUES,is_scatter,ierr)
  deallocate(strided_indices)

  ! resize vertex array to new size
  unstructured_grid%num_vertices_local = vertex_count
  allocate(unstructured_grid%vertices(vertex_count))
  do ivertex = 1, vertex_count
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
    unstructured_grid%vertices(ivertex)%id = needed_vertices_petsc(ivertex) !sp  
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

  unstructured_grid%nlmax = unstructured_grid%num_cells_local
  unstructured_grid%ngmax = unstructured_grid%num_cells_local + &
       unstructured_grid%num_ghost_cells

#endif
  
end subroutine UnstructuredGridDecompose

! ************************************************************************** !
!
! UnstructuredGridCreateUGDM: Constructs mappings / scatter contexts for PETSc DM 
!                                                  object
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
#include "finclude/petscdm.h"  
#include "finclude/petscdm.h90"
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
    !sp 22/10/2010 not strided 
!sp    int_array(icell) = (icell-1)*ndof+istart
    int_array(icell) = (icell-1)+istart
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_petsc,ierr)
  deallocate(int_array)
  
#if GEH_DEBUG  
  call PetscViewerASCIIOpen(option%mycomm,'is_local_petsc.out',viewer,ierr)
  call ISView(ugdm%is_local_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! IS for local numbering of ghosts cells
  allocate(int_array(unstructured_grid%num_ghost_cells))
  do icell = 1, unstructured_grid%num_ghost_cells
    ! sp 
!    int_array(icell) = (icell+unstructured_grid%num_cells_local-1)*ndof
    int_array(icell) = (icell+unstructured_grid%num_cells_local-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_local,ierr)
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
    !sp 
!    int_array(icell) = (unstructured_grid%ghost_cell_ids_petsc(icell)-1)*ndof
    int_array(icell) = (unstructured_grid%ghost_cell_ids_petsc(icell)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_ghost_cells, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosts_petsc,ierr)
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
!sp 
!    int_array(icell) = (icell-1)*ndof
    int_array(icell) = (icell-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_local,ierr)
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
    !sp 
!    int_array(icell) = (icell-1)*ndof
    int_array(icell) = (icell-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_local,ierr)
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
    !sp 
    int_array(unstructured_grid%num_cells_local+icell) = &
      (unstructured_grid%ghost_cell_ids_petsc(icell)-1)
!      (unstructured_grid%ghost_cell_ids_petsc(icell)-1)*ndof
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_petsc,ierr)
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
                       int_array,PETSC_COPY_VALUES,is_tmp,ierr) !sp 
  deallocate(int_array)
  call AOPetscToApplicationIS(unstructured_grid%ao_natural_to_petsc, &
                              is_tmp,ierr)
  ! remap for ndof > 1
  allocate(int_array(unstructured_grid%num_cells_local))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr)
  do icell = 1, unstructured_grid%num_cells_local
    !sp 
!    int_array(icell) = int_ptr(icell)*ndof
    int_array(icell) = int_ptr(icell)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr)
  call ISDestroy(is_tmp,ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_natural,ierr)
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
                                   scatter_ltol, option)

  use Connection_module
  use Option_module
  use Utility_module, only : DotProduct, CrossProduct

  implicit none

  type(connection_set_type), pointer :: UGridComputeInternConnect
  type(option_type) :: option
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(unstructured_grid_type) :: unstructured_grid
  VecScatter :: scatter_ltol 

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
  PetscBool :: found
  PetscBool :: match_found
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
  PetscInt :: nfaces, nfaces2, nvertices, nvertices2, cell_type, cell_type2
  PetscInt :: face_type
  PetscBool:: face_found, vertex_found
  
  PetscReal :: v1(3), v2(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  
  type(plane_type) :: plane1, plane2
  type(point_type) :: point1, point2, point3, point4
  type(point_type) :: point_up, point_dn
  type(point_type) :: intercept1, intercept2


  character(len=MAXSTRINGLENGTH) :: string  

  type(connection_set_type), pointer :: connections

  !sp 
  PetscReal, pointer :: vec_p(:) !sp 
  Vec :: local_vec1 !sp 
  Vec :: local_vec2 !sp 
  PetscInt :: ivert 
  PetscInt :: vert_id  
  PetscErrorCode :: ierr 
  PetscInt, pointer :: cell_vertices_0(:,:) 
  PetscInt :: max_vertex_count
  PetscInt, allocatable :: ltog(:)  
  PetscInt, allocatable :: gtol(:)  
  !sp end 
  

  !sp extend cell_vertices_0 to include ghosted cells 
 
  allocate( ltog( unstructured_grid%num_vertices_local) )
  do ivert=1, unstructured_grid%num_vertices_local 
   ltog(ivert) = unstructured_grid%vertices(ivert)%id 
  end do 

  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_local ,   &
         local_vec1, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_ghosted, &
         local_vec2, ierr)

  max_vertex_count=8 
  allocate( cell_vertices_0(0:max_vertex_count, unstructured_grid%num_cells_ghosted) ) 
  cell_vertices_0=-999 

    ! first the number of vertices per cell (ivertex=0)
    ivertex=0 
    call VecGetArrayF90(local_vec1,vec_p,ierr)
    do icell = 1, unstructured_grid%num_cells_local
     vec_p(icell) =  unstructured_grid%cell_vertices_0(ivertex,icell)  
    enddo
    call VecRestoreArrayF90(local_vec1,vec_p,ierr)

    call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

    call VecGetArrayF90(local_vec2,vec_p,ierr)

    do icell = 1, unstructured_grid%num_cells_ghosted
      cell_vertices_0(ivertex,icell)= vec_p(icell)  
    enddo
    call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! now the cell ids 
  do ivertex=1, max_vertex_count  
    call VecGetArrayF90(local_vec1,vec_p,ierr)
    do icell = 1, unstructured_grid%num_cells_local
     vert_id = unstructured_grid%cell_vertices_0(ivertex,icell)  + 1 
     vec_p(icell) = ltog(vert_id)
    enddo
    call VecRestoreArrayF90(local_vec1,vec_p,ierr)

    call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

    call VecGetArrayF90(local_vec2,vec_p,ierr)

    do icell = 1, unstructured_grid%num_cells_ghosted 
      do ivert=1, unstructured_grid%num_vertices_local 
       vert_id=ltog(ivert) 
       if( vert_id == vec_p(icell) )  exit 
      end do  
      if( vert_id .eq. vec_p(icell) ) cell_vertices_0(ivertex,icell)= ivert  - 1  
    enddo
    call VecRestoreArrayF90(local_vec2,vec_p,ierr) 
  end do 

  deallocate( unstructured_grid%cell_vertices_0) 
  allocate( unstructured_grid%cell_vertices_0(0:max_vertex_count,unstructured_grid%num_cells_ghosted) ) 
  do icell = 1, unstructured_grid%num_cells_ghosted 
   do ivertex=0, max_vertex_count  
      unstructured_grid%cell_vertices_0(ivertex,icell)= cell_vertices_0(ivertex,icell) 
   end do 
  end do 

  call VecDestroy(local_vec1,ierr) 
  call VecDestroy(local_vec2,ierr) 

  !sp end 
   

  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(MAX_VERT_PER_FACE,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  face_to_vertex = -999
  allocate(cell_to_face(MAX_DUALS,unstructured_grid%num_cells_ghosted))
  cell_to_face = -999
  allocate(face_to_cell(2,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  face_to_cell = -999
  allocate(vertex_to_cell(0:MAX_CELLS_SHARING_A_VERTEX,unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  
  allocate(unstructured_grid%face_to_vertex_nindex(MAX_VERT_PER_FACE,&
           MAX_DUALS*unstructured_grid%num_cells_ghosted))
  unstructured_grid%face_to_vertex_nindex = -999
  allocate(unstructured_grid%face_to_cell_locindex(1,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  unstructured_grid%face_to_cell_locindex = -999
  allocate(unstructured_grid%cell_to_face_locindex(MAX_DUALS,&
       unstructured_grid%num_cells_ghosted))
  unstructured_grid%cell_to_face_locindex = -999

  face_count = 0
  do icell = 1, unstructured_grid%num_cells_ghosted
    ! Determine number of faces and cell-type of the current cell
    select case(unstructured_grid%cell_vertices_0(0,icell))
      case(8)
        nfaces = 6
        cell_type = HEX_TYPE
      case(6)
        nfaces = 5
        cell_type = WEDGE_TYPE
      case default
        option%io_buffer = 'Cell type not recognized'
        call printErrMsg(option)
     end select
     do iface = 1, nfaces
      face_count = face_count + 1
      cell_to_face(iface,icell) = face_count
      face_to_cell(1,face_count) = icell
      unstructured_grid%face_to_cell_locindex(1,face_count) = icell
      call UGridGetCellFaceVertices(option,cell_type,iface,vertex_ids4)
      ! Determine number of vertices forming the face
      nvertices = 4
      if((cell_type.eq.WEDGE_TYPE).and.(iface.gt.3)) nvertices = 3
      do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices_0(vertex_ids4(ivertex),icell)+1
        if(face_to_vertex(ivertex,face_count).gt.0) then
          unstructured_grid%face_to_vertex_nindex(ivertex,face_count) = &
            unstructured_grid%vertex_ids_nindex(face_to_vertex(ivertex,face_count))
        endif
      enddo
    enddo
  enddo

#ifdef MIXED_UMESH
  !
  ! Remove duplicate faces:
  !
  ! A cell (cell_id) and Neighboring-Cell (cell_id2) will only share ONE face.
  ! Find the face that cell_id ane cell_id2 share and remove it.
  !
  ! Method:
  !        - Pick i-th face (iface) of cell_id and check if ALL the vertices of
  !          the iface are present in cell_id2. If all the vertices of iface are
  !          not present in cell_id2, move to the next face.
  !        - After finding the iface, now find iface2 in cell_id2 that
  !          corresponds to iface.
  !        - Check to ensure that atleast on face of cell_id is shared
  !          with cell_id2.
  !
  !
  !
  ! NOTE: For a cell_type = WEDGE_TYPE, faces 1-3 have 4 vertices; while
  !       faces 4-5 have 3 vertices
  !
  do icell = 1, unstructured_grid%num_cells_local
    ! Selet a cell and find number of vertices
    cell_id = icell
    select case(unstructured_grid%cell_vertices_0(0,icell))
      case(8)
        nfaces = 6
        cell_type = HEX_TYPE
      case(6)
        nfaces = 5
        cell_type = WEDGE_TYPE
    end select

    do icell2 = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      ! Selet a neighboring cell
      cell_id2 =  unstructured_grid%cell_neighbors_local_ghosted(icell2,icell)
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      ! Find the number of vertices for neighboring cell
      select case(unstructured_grid%cell_vertices_0(0,icell2))
        case(8)
          nfaces2 = 6
          cell_type2 = HEX_TYPE
        case(6)
          nfaces2 = 5
          cell_type2 = WEDGE_TYPE
      end select
      ! Initialize
      face_found = PETSC_FALSE
      do iface = 1, nfaces
        ! Select a face and find number of vertices forming the face
        face_id = cell_to_face(iface,cell_id)
        nvertices = 4
        if((cell_type.eq.WEDGE_TYPE).and.(iface.gt.3)) nvertices = 3
        do ivertex = 1,nvertices
          ! Select a vertex and initialize vertex_found
          vertex_id = face_to_vertex(ivertex,face_id) ! face_to_vertex is 1-based indexing
          vertex_found = PETSC_FALSE
          do ivertex2 = 1,unstructured_grid%cell_vertices_0(0,cell_id2)
            vertex_id2 = unstructured_grid%cell_vertices_0(ivertex2,cell_id2)+1 ! cell_vertices_0 is 0-based indexing
            if(vertex_id.eq.vertex_id2) then
              vertex_found = PETSC_TRUE
              exit
            endif
          enddo
          !
          ! If ivertex of iface of the Cell is not present as vertices of the
          ! Neighboring-Cell, then iface is not the shared face. Skip iterating
          ! over the remaing vertices of iface
          if(.not.vertex_found) exit
        enddo
        
        if(vertex_found) then
          ! All the vertices of iface are present in the Neighboring cells.
          ! Thus, iface is the shared face.
          face_found = PETSC_TRUE
          
          ! Now, we have to find iface2 that corresponds to iface
          do iface2 = 1, nfaces2
            face_id2 = cell_to_face(iface2,cell_id2)
            nvertices2 = 4
            if((cell_type.eq.WEDGE_TYPE).and.(iface.gt.3)) nvertices2 = 3
            
            ! Both iface and iface2 need to have same number of vertices
            if(nvertices.eq.nvertices2) then
              ! Count the number of vertices of iface which match vertices
              ! of iface2
              num_match = 0
              do ivertex = 1,nvertices
                vertex_id = face_to_vertex(ivertex,face_id)
                vertex_found = PETSC_TRUE
                
                do ivertex2 = 1,nvertices
                  vertex_id2 = face_to_vertex(ivertex2,face_id2)
                  if (vertex_id.eq.vertex_id2) then
                    vertex_found = PETSC_TRUE
                    num_match = num_match + 1
                    vertex_ids4(num_match) = vertex_id
                    exit
                  endif
                enddo
                !
                ! If vertex_id of face_id not found as one of vertices of face_id2,
                ! face_id2 is not shared between cells
                if(.not.vertex_found) exit
              enddo
              if(num_match.eq.nvertices) then
                ! remove duplicate face
                if (face_id2 > face_id) then
                  write(string,*) option%myrank, face_id2, ' -> ', face_id
                  option%io_buffer = 'Duplicated face removed:' // string //'\n'
                  !call printMsg(option)
                  cell_to_face(iface2,cell_id2) = face_id
                  !! flag as removed
                  face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
                  face_to_cell(2,face_id) = cell_id2
                else
                  write(string,*) option%myrank, face_id, ' -> ', face_id2
                  option%io_buffer = 'Duplicated face removed:' // string //'\n'
                  !call printMsg(option)
                  cell_to_face(iface,cell_id) = face_id2
                  !! flag as removed
                  face_to_cell(1,face_id) = -face_to_cell(1,face_id)
                  face_to_cell(2,face_id2) = cell_id
                endif
                exit
              endif
            endif
          enddo
          exit
        endif
      enddo ! iface-loop
      
      ! Check that one shared face was found between the Cell and Neighboring-Cell
      if(.not.face_found) then
        write(string,*),'rank=',option%myrank, 'icell_id',cell_id,'icell_id2',cell_id2
        option%io_buffer='No shared face found: ' // string // '\n'
        write(*,*),'No shared face found: rank=',option%myrank, 'icell_id',cell_id,'icell_id2',cell_id2
        call printErrMsg(option)
      endif
    enddo ! icell2-loop
  enddo  ! icell-loop
#endif ! #ifdef MIXED_UMESH

   

#ifndef MIXED_UMESH
  ! remove duplicate faces
  ! fill face ids
  do icell = 1, unstructured_grid%num_cells_local !sp   was num_cells_ghosted 
    cell_id = icell
    do icell2 = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      cell_id2 =  unstructured_grid%cell_neighbors_local_ghosted(icell2,icell)
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
              write(string,*) option%myrank, face_id2, ' -> ', face_id
              option%io_buffer = 'Duplicated face removed:' // string
              !call printMsg(option)
              cell_to_face(iface2,cell_id2) = face_id
              ! flag as removed
              face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
              face_to_cell(2,face_id) = cell_id2
            else
              write(string,*) option%myrank, face_id, ' -> ', face_id2
              option%io_buffer = 'Duplicated face removed:' // string
              !call printMsg(option)
              cell_to_face(iface,cell_id) = face_id2
              ! flag as removed
              face_to_cell(1,face_id) = -face_to_cell(1,face_id)
              face_to_cell(2,face_id2) = cell_id
            endif
          endif
        else
          match_found = PETSC_FALSE
        endif  
        if (.not.match_found) then
          option%io_buffer = 'Matching faces not found'
          call printErrMsg(option)
        endif
      endif    
    enddo
  enddo

#endif ! #ifndef MIXED_UMESH

  
  ! count up the # of faces
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) &
      face_count = face_count + 1
  enddo
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
        if(face_id<0) cycle
        if (face_id == temp_int(face_id2)) then
          found = PETSC_TRUE
          cell_to_face(iface2,cell_id) = face_id
          exit
        endif
      enddo
      if (.not.found) then
        option%io_buffer = 'Remapping of cell face id unsuccessful'
        call printErrMsg(option)
      endif
    enddo
  enddo
  deallocate(temp_int)
  
  
  do icell = 1, unstructured_grid%num_cells_ghosted
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,icell)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,icell)+1
      if( vertex_id <= 0) cycle 
      count = vertex_to_cell(0,vertex_id) + 1
      if( count .gt. MAX_CELLS_SHARING_A_VERTEX) then
        write(string,*) 'Vertex can be shared by at most by ',MAX_CELLS_SHARING_A_VERTEX, &
              ' cells. Rank = ', option%myrank, ' vertex_id = ', vertex_id, ' exceeds it.'
        option%io_buffer = string
        call printErrMsg(option)
      endif
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
!      if (dual_id < 0 .or. icell < dual_id) then
      if (dual_id > 0 .and. icell < dual_id) then !sp 
          nconn = nconn + 1
      endif
    enddo
  enddo


  connections => ConnectionCreate(nconn,option%nphase,INTERNAL_CONNECTION_TYPE)
  
  allocate(unstructured_grid%face_area(face_count))
  allocate(dual_to_face(nconn))
  dual_to_face = -999

  ! loop over connection again
  iconn = 0
  do icell = 1, unstructured_grid%num_cells_local
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,icell)
      dual_icell = unstructured_grid%cell_neighbors_local_ghosted(idual,icell)
      if (icell < abs(dual_icell)) then 
        iconn = iconn + 1
        ! find face
        found = PETSC_FALSE
        do iface = 1, unstructured_grid%cell_vertices_0(0,icell)
          face_id = cell_to_face(iface,icell)
          do icell2 = 1,2
            cell_id2 = face_to_cell(icell2,face_id)
            if (cell_id2 == abs(dual_icell)) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (found) exit
        enddo
        if (found) then
          dual_to_face(iconn) = face_id
        else
          write(string,*) option%myrank,icell,dual_icell 
          option%io_buffer = 'face not found in connection loop' // string 
          call printErrMsg(option)
        endif
        if(face_to_vertex(4,face_id).lt.0) then
          face_type = TRI_FACE_TYPE
        else
          face_type = QUAD_FACE_TYPE
        endif
        connections%id_up(iconn) = icell
        connections%id_dn(iconn) = abs(dual_icell)
        ! need to add the surface areas, distance, etc.
        point1 = unstructured_grid%vertices(face_to_vertex(1,face_id))
        point2 = unstructured_grid%vertices(face_to_vertex(2,face_id))
        point3 = unstructured_grid%vertices(face_to_vertex(3,face_id))
        if (face_type.eq.QUAD_FACE_TYPE) then
          point4 = unstructured_grid%vertices(face_to_vertex(4,face_id))
        else
          point4 = unstructured_grid%vertices(face_to_vertex(3,face_id))
        endif
        
        call ComputePlane(plane1,point1,point2,point3)
        if (face_type.eq.QUAD_FACE_TYPE) then
          call ComputePlane(plane2,point3,point4,point1)
        else
          call ComputePlane(plane2,point1,point2,point3)
        endif
        
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
        !area1 = 0.5d0*DotProduct(n1,n1)
        !area1 = area1*DotProduct(n1,n_up_dn)
        area1 = 0.5d0*sqrt(DotProduct(n1,n1))
        
        v1(1) = point1%x-point4%x
        v1(2) = point1%y-point4%y
        v1(3) = point1%z-point4%z
        v2(1) = point3%x-point4%x
        v2(2) = point3%y-point4%y
        v2(3) = point3%z-point4%z
        n2 = CrossProduct(v1,v2)
        area2 = 0.5d0*DotProduct(n2,n2)
        !area2 = area2*DotProduct(n2,n_up_dn)
        !area1 = 0.5*|(p2-p1)X(p3-p1)|.
        area2 = 0.5d0*sqrt(DotProduct(n2,n2))
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
        !connections%dist(1:3,iconn) = v1 + v2
        connections%dist(1:3,iconn) = (v1 + v2)/sqrt(DotProduct(v1+v2,v1+v2))
        connections%area(iconn) = area1 + area2
       
      endif
    enddo
  enddo
  
  ! Save area and centroid of faces
  unstructured_grid%cell_to_face_locindex = cell_to_face
  allocate(unstructured_grid%face_centroid(face_count))
  do iface = 1,face_count
    unstructured_grid%face_centroid(iface)%id = -999
  enddo
  
  do icell = 1, unstructured_grid%num_cells_local
    do iface = 1,MAX_DUALS
      face_id = cell_to_face(iface, icell)
      if( unstructured_grid%face_centroid(face_id)%id.eq.-999) then
        count = 0
        unstructured_grid%face_centroid(face_id)%x = 0.d0
        unstructured_grid%face_centroid(face_id)%y = 0.d0
        unstructured_grid%face_centroid(face_id)%z = 0.d0
        if(face_to_vertex(4,face_id).lt.0) then
          face_type = TRI_FACE_TYPE
        else
          face_type = QUAD_FACE_TYPE
        endif
        point1 = unstructured_grid%vertices(face_to_vertex(1,face_id))
        point2 = unstructured_grid%vertices(face_to_vertex(2,face_id))
        point3 = unstructured_grid%vertices(face_to_vertex(3,face_id))
        if (face_type.eq.QUAD_FACE_TYPE) then
          point4 = unstructured_grid%vertices(face_to_vertex(4,face_id))
        else
          point4 = unstructured_grid%vertices(face_to_vertex(3,face_id))
        endif
        v1(1) = point3%x-point2%x
        v1(2) = point3%y-point2%y
        v1(3) = point3%z-point2%z
        v2(1) = point1%x-point2%x
        v2(2) = point1%y-point2%y
        v2(3) = point1%z-point2%z
        n1 = CrossProduct(v1,v2)
        area1 = 0.5d0*sqrt(DotProduct(n1,n1))
        
        v1(1) = point1%x-point4%x
        v1(2) = point1%y-point4%y
        v1(3) = point1%z-point4%z
        v2(1) = point3%x-point4%x
        v2(2) = point3%y-point4%y
        v2(3) = point3%z-point4%z
        n2 = CrossProduct(v1,v2)
        area2 = 0.5d0*DotProduct(n2,n2)
        unstructured_grid%face_area(face_id) = area1 + area1
        
        do ivert = 1,MAX_VERT_PER_FACE
          vertex_id = face_to_vertex(ivert,face_id)
          if(vertex_id.ne.-999) then
            unstructured_grid%face_centroid(face_id)%x = &
              unstructured_grid%face_centroid(face_id)%x + unstructured_grid%vertices(vertex_id)%x
            unstructured_grid%face_centroid(face_id)%y = &
              unstructured_grid%face_centroid(face_id)%y + unstructured_grid%vertices(vertex_id)%y
            unstructured_grid%face_centroid(face_id)%z = &
              unstructured_grid%face_centroid(face_id)%z + unstructured_grid%vertices(vertex_id)%z
            count = count +1
          endif
        enddo
        unstructured_grid%face_centroid(face_id)%id = face_id
        unstructured_grid%face_centroid(face_id)%x  = &
          unstructured_grid%face_centroid(face_id)%x/count
        unstructured_grid%face_centroid(face_id)%y  = &
          unstructured_grid%face_centroid(face_id)%y/count
        unstructured_grid%face_centroid(face_id)%z  = &
          unstructured_grid%face_centroid(face_id)%z/count
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
! UGridPopulateConnection: Computes details of connection (area, dist, etc)
! author: Gautam Bisht
! date: 10/30/09
!
! ************************************************************************** !
subroutine UGridPopulateConnection(unstructured_grid, connection, iface, &
                                   iconn, ghosted_id)

  use Connection_module
  use Utility_module, only : DotProduct
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface,iface_cell
  PetscInt :: iconn
  PetscInt :: ghosted_id
  
  PetscErrorCode :: ierr
  
  PetscInt  :: face_id
  PetscInt  :: ivert,vert_id
  PetscReal :: v1(3),v2(3),n_dist(3), dist
  type(point_type) :: vertex_8(8)
  
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
    
      select case (unstructured_grid%cell_vertices_0(0, ghosted_id))
        case(8)
#if 0
          select case (iface)
            case(WEST_FACE)
              iface_cell = 4
            case(EAST_FACE)
              iface_cell = 2
            case(SOUTH_FACE)
              iface_cell = 1
            case(NORTH_FACE)
              iface_cell = 3
            case(BOTTOM_FACE)
              iface_cell = 5
            case(TOP_FACE)
              iface_cell = 6
          end select
#else
          iface_cell = iface
#endif
        case(6)
          write(*,*), 'UGridPopulateConnection: Add code for WEDGE_TYPE'
      end select
    
      ! Get face-centroid vector
      face_id = unstructured_grid%cell_to_face_locindex(iface_cell, ghosted_id)
      v1(1) = unstructured_grid%face_centroid(face_id)%x
      v1(2) = unstructured_grid%face_centroid(face_id)%y
      v1(3) = unstructured_grid%face_centroid(face_id)%z
      
      ! Compute cell centeroid
      v2 = 0d0
      do ivert = 1, unstructured_grid%cell_vertices_0(0, ghosted_id)
        vert_id = unstructured_grid%cell_vertices_0(ivert, ghosted_id) + 1
        vertex_8(ivert)%x = unstructured_grid%vertices(vert_id)%x
        vertex_8(ivert)%y = unstructured_grid%vertices(vert_id)%y
        vertex_8(ivert)%z = unstructured_grid%vertices(vert_id)%z
      enddo
      select case (unstructured_grid%cell_vertices_0(0, ghosted_id))
        case(8)
          v2 = ComputeCentroid(HEX_TYPE, vertex_8)
        case(6)
          v2 = ComputeCentroid(WEDGE_TYPE, vertex_8)
      end select

      ! Compute distance vector: cell_center - face_centroid
      v1(1) = v2(1) - v1(1)
      v1(2) = v2(2) - v1(2)
      v1(3) = v2(3) - v1(3)
      
      !
      dist = sqrt(DotProduct(v1, v1))
      n_dist = v1/dist
      connection%dist(0, iconn) = dist
      connection%dist(1, iconn) = n_dist(1)
      connection%dist(2, iconn) = n_dist(2)
      connection%dist(3, iconn) = n_dist(3)
      connection%area(iconn)    = unstructured_grid%face_area(face_id)
      
  end select
  
end subroutine UGridPopulateConnection

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
    case(WEDGE_TYPE)
      select case(iface)
        case(1)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 5
          vertex_ids(4) = 4
       case(2)
          vertex_ids(1) = 2
          vertex_ids(2) = 3
          vertex_ids(3) = 6
          vertex_ids(4) = 5
       case(3)
          vertex_ids(1) = 3
          vertex_ids(2) = 1
          vertex_ids(3) = 4
          vertex_ids(4) = 6
       case(4)
          vertex_ids(1) = 1
          vertex_ids(2) = 2
          vertex_ids(3) = 3
       case(5)
          vertex_ids(1) = 4
          vertex_ids(2) = 5
          vertex_ids(3) = 6
       case default
          option%io_buffer='Cell WEDGE_TYPE has only 5 faces'
          call printErrMsg(option)
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
! 11/2/10 Major rewrite to extend coordinates to ghost cells SP and GEH 
!
! ************************************************************************** !
subroutine UGridComputeCoord(unstructured_grid,option, &
                             scatter_ltol, & 
                             grid_x,grid_y,grid_z, &
                             x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  VecScatter :: scatter_ltol 
  
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: icell
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal :: centroid(3)
  PetscReal, pointer :: vec_p(:) !sp 
  Vec :: local_vec1 !sp 
  Vec :: local_vec2 !sp 
  PetscErrorCode :: ierr 

  do icell = 1, unstructured_grid%num_cells_local 
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,icell)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,icell) + 1
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    select case (unstructured_grid%cell_vertices_0(0,icell))
      case(8)
        centroid = ComputeCentroid(HEX_TYPE,vertex_8)
      case(6)
        centroid = ComputeCentroid(WEDGE_TYPE,vertex_8)
    end select
    grid_x(icell) = centroid(1)
    grid_y(icell) = centroid(2)
    grid_z(icell) = centroid(3)
  enddo

  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_local,   & 
         local_vec1, ierr) 
  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_ghosted, & 
         local_vec2, ierr) 

  ! x coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do icell = 1, unstructured_grid%num_cells_local
   vec_p(icell) = grid_x(icell)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do icell = 1, unstructured_grid%num_cells_ghosted
   grid_x(icell) = vec_p(icell)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! y coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do icell = 1, unstructured_grid%num_cells_local
   vec_p(icell) = grid_y(icell)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do icell = 1, unstructured_grid%num_cells_ghosted
   grid_y(icell) = vec_p(icell)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! z coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do icell = 1, unstructured_grid%num_cells_local
   vec_p(icell) = grid_z(icell)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do icell = 1, unstructured_grid%num_cells_ghosted
   grid_z(icell) = vec_p(icell)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)


  call VecDestroy(local_vec1,ierr) 
  call VecDestroy(local_vec2,ierr) 

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
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,ghosted_id) + 1
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
      select case (unstructured_grid%cell_vertices_0(0,ghosted_id))
        case(8)
          volume_p(local_id) = ComputeVolume(HEX_TYPE,vertex_8)
        case(6)
          volume_p(local_id) = ComputeVolume(WEDGE_TYPE,vertex_8)
      end select
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
    case(WEDGE_TYPE)
      ! need something more sophisticated, but for now, just use average
      do ivertex = 1, 6
        ComputeCentroid(1) = ComputeCentroid(1) + vertices(ivertex)%x
        ComputeCentroid(2) = ComputeCentroid(2) + vertices(ivertex)%y
        ComputeCentroid(3) = ComputeCentroid(3) + vertices(ivertex)%z
      enddo
      ComputeCentroid = ComputeCentroid / 6.d0
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

  use Utility_module, only : DotProduct, CrossProduct

  implicit none
  
  PetscInt :: cell_type
  type(point_type) :: vertices(*)
  
  PetscReal :: ComputeVolume
  PetscReal :: v(3)
  PetscReal :: l1, l2, l3
  PetscReal :: n1(3), area1, dz,v1(3),v2(3)
  
  ComputeVolume = 0.d0
  select case(cell_type)
    case(HEX_TYPE)
     ! need something more sophisticated, but for now, just use volume of 
     ! box
      v(1) = vertices(2)%x-vertices(1)%x
      v(2) = vertices(2)%y-vertices(1)%y
      v(3) = vertices(2)%z-vertices(1)%z
      l1 = sqrt(DotProduct(v,v))
      v(1) = vertices(4)%x-vertices(1)%x
      v(2) = vertices(4)%y-vertices(1)%y
      v(3) = vertices(4)%z-vertices(1)%z
      l2 = sqrt(DotProduct(v,v))
      v(1) = vertices(5)%x-vertices(1)%x
      v(2) = vertices(5)%y-vertices(1)%y
      v(3) = vertices(5)%z-vertices(1)%z
      l3 = sqrt(DotProduct(v,v))
      ComputeVolume = l1*l2*l3
    case(WEDGE_TYPE)
      v1(1) = vertices(3)%x-vertices(2)%x
      v1(2) = vertices(3)%y-vertices(2)%y
      v1(3) = vertices(3)%z-vertices(2)%z
      v2(1) = vertices(1)%x-vertices(2)%x
      v2(2) = vertices(1)%y-vertices(2)%y
      v2(3) = vertices(1)%z-vertices(2)%z
      n1 = CrossProduct(v1,v2)
      area1 = 0.5d0*sqrt(DotProduct(n1,n1))
      dz = (vertices(1)%z+vertices(2)%z+vertices(3)%z)/3.d0 - &
           (vertices(4)%z+vertices(5)%z+vertices(6)%z)/3.d0
      ComputeVolume = abs(dz)*area1
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
!  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*ugdm%ndof
        o_nnz = o_nnz*ugdm%ndof
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog,ugdm%mapping_ltog,ierr)
!        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltogb,ierr)
      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
!        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltogb,ugdm%mapping_ltogb,ierr)
      case default
        option%io_buffer = 'MatType not recognized in UGDMCreateJacobian'
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
!        option%io_buffer = 'MatType not recognized in UGDMCreateJacobian'
!        call printErrMsg(option)
!    end select
!  endif

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
  if (associated(unstructured_grid%cell_vertices_nindex))&
    deallocate(unstructured_grid%cell_vertices_nindex)
  nullify(unstructured_grid%cell_vertices_nindex)
  if (associated(unstructured_grid%face_to_cell_locindex))&
    deallocate(unstructured_grid%face_to_cell_locindex)
  nullify(unstructured_grid%face_to_cell_locindex)
  if (associated(unstructured_grid%face_to_vertex_nindex))&
    deallocate(unstructured_grid%face_to_vertex_nindex)
  nullify(unstructured_grid%face_to_vertex_nindex)

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
