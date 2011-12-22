module Unstructured_Grid_module

  use Connection_module
  use Unstructured_Cell_module
  
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
    PetscInt :: num_ghost_cells   ! number of ghost cells (only) on processor
    PetscInt :: num_vertices_global ! number of vertices in entire problem domain
    PetscInt :: num_vertices_local  ! number of vertices in local grid cells
    PetscInt :: global_offset ! offset in petsc ordering for the first cell on a processor???
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash
    PetscInt :: max_ndual_per_cell
    PetscInt :: max_nvert_per_cell
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
    PetscInt, pointer :: cell_ids_natural(:) ! natural 1d right-hand i,j,k ordering
    PetscInt, pointer :: cell_ids_petsc(:) ! petsc ordering of cell ids
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

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4
  !  PetscInt, parameter :: MAX_CELLS_SHARING_A_VERTEX = 16

  public :: UGridCreate, &
            UGridRead, &
#ifndef SAMR_HAVE_HDF5
            UGridReadHDF5, &
#endif
#if defined(PARALLELIO_LIB)
            UGridReadHDF5PIOLib, &
#endif
            UGridDecompose, &
            UGridComputeInternConnect, &
            UGridPopulateConnection, &
            UGridComputeCoord, &
            UGridComputeVolumes, &
            UGridMapIndices, &
            UGridDMCreateJacobian, &
            UGridDMCreateVector, &
            UGridGetCellFromPoint, &
            UGridGetCellsInRectangle, &
            UGridEnsureRightHandRule, &
            UGridMapSideSet, &
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
  ugdm%mapping_ltog = 0
  ugdm%mapping_ltogb = 0
  ugdm%global_vec = 0
  ugdm%local_vec = 0

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

  unstructured_grid%num_vertices_global = 0
  unstructured_grid%num_vertices_local = 0
  unstructured_grid%num_ghost_cells = 0
  unstructured_grid%global_offset = 0
  unstructured_grid%nmax = 0
  unstructured_grid%nlmax = 0
  unstructured_grid%ngmax = 0
  unstructured_grid%max_ndual_per_cell = 0
  unstructured_grid%max_nvert_per_cell = 0
  nullify(unstructured_grid%cell_type)
  nullify(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%cell_to_face_ghosted)
  nullify(unstructured_grid%vertex_ids_natural)
  nullify(unstructured_grid%vertices)
  nullify(unstructured_grid%hash)
  nullify(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
  nullify(unstructured_grid%connection_to_face)
  unstructured_grid%num_hash = 100
  unstructured_grid%ao_natural_to_petsc = 0

  UGridCreate => unstructured_grid
  
end function UGridCreate

  

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
! UGridRead: Reads an unstructured grid
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UGridRead(unstructured_grid,filename,option)

  use Input_module
  use Option_module
  use String_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card, word
  PetscInt :: num_cells_local_save
  PetscInt :: num_cells_local
  PetscInt :: num_vertices_local_save
  PetscInt :: num_vertices_local
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

  ! initial guess is 8 vertices per cell
  unstructured_grid%max_nvert_per_cell = 8

! Format of unstructured grid file
! type: H=hexahedron, T=tetrahedron, W=wedge, P=pyramid
! vertn(H) = 8
! vertn(T) = 4
! vertn(W) = 6
! vertn(P) = 5
! -----------------------------------------------------------------
! num_cells num_vertices  (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 1 (integers)
! type vert1 vert2 vert3 ... vertn  ! for cell 2
! ...
! ...
! type vert1 vert2 vert3 ... vertn  ! for cell num_cells
! xcoord ycoord zcoord ! coordinates of vertex 1 (real)
! xcoord ycoord zcoord ! coordinates of vertex 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of vertex num_vertices (real)
! -----------------------------------------------------------------

  card = 'Unstructured Grid'

  call InputReadFlotranString(input,option)
  string = 'unstructured grid'
  call InputReadStringErrorMsg(input,option,card)  

  ! read num_cells
  call InputReadInt(input,option,unstructured_grid%nmax)
  call InputErrorMsg(input,option,'number of cells',card)
  ! read num_vertices
  call InputReadInt(input,option,unstructured_grid%num_vertices_global)
  call InputErrorMsg(input,option,'number of vertices',card)

  ! divide cells across ranks
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices(unstructured_grid%max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -999

  ! for now, read all cells from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(unstructured_grid%max_nvert_per_cell, &
                            num_cells_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = -999
      num_to_read = num_cells_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do icell = 1, num_to_read
        ! read in the vertices defining the grid cell
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'element type',card)
        call StringToUpper(word)
        select case(word)
          case('H')
            num_vertices = 8
          case('W')
            num_vertices = 6
          case('P')
            num_vertices = 5
          case('T')
            num_vertices = 4
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,icell))
          call InputErrorMsg(input,option,'vertex id',card)
        enddo
      enddo
      
      ! if the cells reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_cells_local
        string = trim(adjustl(string)) // ' cells stored on p0'
        print *, trim(string)
#endif
        unstructured_grid%cell_vertices(:,1:num_cells_local) = &
          temp_int_array(:,1:num_cells_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' cells sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*unstructured_grid%max_nvert_per_cell
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_cells_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' cells received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_cells_local*unstructured_grid%max_nvert_per_cell
    call MPI_Recv(unstructured_grid%cell_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif


  ! divide vertices across ranks
  num_vertices_local = unstructured_grid%num_vertices_global/ &
                                         option%mycommsize 
  num_vertices_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                 num_vertices_local + 1

  allocate(vertex_coordinates(3,num_vertices_local))
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
        vertex_coordinates(:,1:num_vertices_local) = &
          temp_real_array(:,1:num_vertices_local)
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    int_mpi = num_vertices_local*3
    call MPI_Recv(vertex_coordinates, &
                  int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif
  
  ! fill the vertices data structure
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ivertex = 1, num_vertices_local
    unstructured_grid%vertices(ivertex)%id = 0
    unstructured_grid%vertices(ivertex)%x = vertex_coordinates(1,ivertex)
    unstructured_grid%vertices(ivertex)%y = vertex_coordinates(2,ivertex)
    unstructured_grid%vertices(ivertex)%z = vertex_coordinates(3,ivertex)
  enddo
  deallocate(vertex_coordinates)

  unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local

  call InputDestroy(input)

end subroutine UGridRead

#ifndef SAMR_HAVE_HDF5

! ************************************************************************** !
!
! UGridReadHDF5: Reads an unstructured grid from HDF5
! author: Gautam Bisht
! date: 04/25/11
!
! ************************************************************************** !
subroutine UGridReadHDF5(unstructured_grid,filename,option)

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
  PetscInt          :: num_cells_local
  PetscInt          :: num_cells_local_save
  PetscInt          :: num_vertices_local
  PetscInt          :: num_vertices_local_save
  PetscInt          :: remainder
  PetscInt,pointer  :: int_buffer(:,:)
  PetscReal,pointer :: double_buffer(:,:)
  PetscInt, parameter :: max_nvert_per_cell = 8  
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
  unstructured_grid%nmax = dims_h5(2)
  num_cells_local = unstructured_grid%nmax/option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%nmax - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                  num_cells_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_cells_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_cells_local, iend, ONE_INTEGER_MPI, &
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
  allocate(unstructured_grid%cell_vertices(max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -1
  
  do ii = 1, num_cells_local
    do jj = 2, int_buffer(1,ii) + 1
      unstructured_grid%cell_vertices(jj-1, ii) = int_buffer(jj, ii)
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
  num_vertices_local  = &
                                       unstructured_grid%num_vertices_global/ &
                                       option%mycommsize 
  num_cells_local_save = num_vertices_local
  remainder = unstructured_grid%num_vertices_global - &
              num_vertices_local*option%mycommsize
  if (option%myrank < remainder) num_vertices_local = &
                                  num_vertices_local + 1
  
  ! Find istart and iend
  istart = 0
  iend   = 0
  call MPI_Exscan(num_vertices_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  call MPI_Scan(num_vertices_local, iend, ONE_INTEGER_MPI, &
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
  allocate(unstructured_grid%vertices(num_vertices_local))
  do ii = 1, num_vertices_local
    unstructured_grid%vertices(ii)%id = 0
    unstructured_grid%vertices(ii)%x = double_buffer(1, ii)
    unstructured_grid%vertices(ii)%y = double_buffer(2, ii)
    unstructured_grid%vertices(ii)%z = double_buffer(3, ii)
  enddo
  
  
  deallocate(double_buffer)
  deallocate(dims_h5)
  deallocate(max_dims_h5)
  
  unstructured_grid%max_nvert_per_cell = max_nvert_per_cell
  unstructured_grid%nlmax = num_cells_local
  unstructured_grid%num_vertices_local = num_vertices_local
  
end subroutine UGridReadHDF5

#endif

! ************************************************************************** !
!
! UGridReadHDF5PIOLib: Reads an unstructured grid from HDF5
! author: Gautam Bisht
! date: 05/13/11
!
! ************************************************************************** !
#if defined(PARALLELIO_LIB)

subroutine UGridReadHDF5PIOLib(unstructured_grid, filename, &
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
  PetscInt, parameter :: max_nvert_per_cell = 8  

  character(len=MAXSTRINGLENGTH) :: cell_dataset_name = &
                                                       '/Domain/Cells'//CHAR(0)
  character(len=MAXSTRINGLENGTH) :: vert_dataset_name = &
                                                    '/Domain/Vertices'//CHAR(0)

  ! Read Domain/Cells
  call HDF5ReadDatasetInteger2D(filename, &
                                cell_dataset_name, &
                                NONUNIFORM_CONTIGUOUS_READ, &
                                option, &
                                int_buffer, &
                                dims, &
                                dataset_dims)

  ! Allocate array to store vertices for each cell
  num_cells_local  = dims(2)
  unstructured_grid%nmax = dataset_dims(2)
  allocate(unstructured_grid%cell_vertices(max_nvert_per_cell, &
                                             num_cells_local))
  unstructured_grid%cell_vertices = -1

  ! Fill the cell data structure
  do ii = 1, num_cells_local
    do jj = 2, int_buffer(1, ii) + 1
      unstructured_grid%cell_vertices(jj-1, ii) = int_buffer(jj, ii)
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

  unstructured_grid%max_nvert_per_cell = max_nvert_per_cell

end subroutine UGridReadHDF5PIOLib

#endif

! ************************************************************************** !
!
! UGridDecompose: Decomposes an unstructured grid across ranks
! author: Glenn Hammond
! date: 09/30/09
!
! ************************************************************************** !
subroutine UGridDecompose(unstructured_grid,option)
  
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
  
  PetscInt :: local_id, local_id2
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  PetscInt :: count, vertex_count
  PetscInt :: vertex_offset, global_vertex_offset
  PetscInt :: stride
  PetscInt, allocatable :: local_vertices(:)
  PetscInt, allocatable :: local_vertex_offset(:)
  PetscInt :: index_format_flag, num_common_vertices
  PetscInt, allocatable :: cell_counts(:)
  PetscInt, pointer :: index_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
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
  Vec :: elements_local
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
  PetscInt :: min_value
  PetscInt :: num_cells_local_new
  PetscInt :: num_cells_local_old  
  PetscInt :: global_offset_old
  PetscInt :: global_offset_new
  PetscInt, allocatable :: int_array(:)
  PetscInt, allocatable :: int_array2(:)
  PetscInt, allocatable :: int_array3(:)
  PetscInt, allocatable :: int_array4(:)
  PetscInt, allocatable :: int_array5(:)
  PetscInt, allocatable :: needed_vertices_petsc(:)
  PetscInt, pointer :: int_array_pointer(:)
  
  PetscInt :: idual, dual_id
  PetscInt :: iflag
  PetscBool :: found
  
!  cell distribution across processors (size = num_cores + 1)
!  core i owns cells cell_distribution(i):cell_distribution(i+1), note
!  the zero-based indexing
!  allocate(cell_distribution(option%mycommsize+1))
!  call MPI_Scan(unstructured_grid%nlmax,
!  cell_distribution(1) = 0
!  cell_distribution(2:) = unstructured_grid%num_cells
!  num_local_cells = cell_distribution(option%myrank+1)- &
!                    cell_distribution(option%myrank+2)

  num_cells_local_old = unstructured_grid%nlmax

  ! recalculate maximum number of vertices for any given cell
  temp_int = 0
  min_value = 2 ! min value should be either 0 or 1 after global reduction
  do local_id = 1, num_cells_local_old
    vertex_count = 0
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point, cell vertex can be 0
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) exit
      if (unstructured_grid%cell_vertices(ivertex,local_id) < min_value) then
        min_value = unstructured_grid%cell_vertices(ivertex,local_id)
      endif
      vertex_count = vertex_count+1
    enddo
    if (vertex_count > temp_int) temp_int = vertex_count
  enddo
  call MPI_Allreduce(temp_int,unstructured_grid%max_nvert_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(min_value,index_format_flag, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN,option%mycomm,ierr)

  ! let's make it Fortran indexing
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! at this point we may be zero-based
      if (unstructured_grid%cell_vertices(ivertex,local_id) < 0) then
        ! change no_value (-999) to '0'
        unstructured_grid%cell_vertices(ivertex,local_id) = 0
      else
        if (index_format_flag == 0) then
          ! let's make it Fortran indexing
          unstructured_grid%cell_vertices(ivertex,local_id) = &
            unstructured_grid%cell_vertices(ivertex,local_id) + 1
        endif
      endif
    enddo
  enddo

#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_nvert_per_cell
  option%io_buffer = 'Maximum number of vertices per cell: ' // adjustl(string)
  call printMsg(option)
  write(string,*) index_format_flag
  option%io_buffer = 'Vertex indexing starts at: ' // adjustl(string)
  call printMsg(option)
  if (index_format_flag == 0) then
    option%io_buffer = 'Changing vertex indexing to 1-based.'
    call printMsg(option)
  endif
#endif

  num_cells_local_old = unstructured_grid%nlmax 
  allocate(local_vertices(unstructured_grid%max_nvert_per_cell* &
                          num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  local_vertices = 0
  local_vertex_offset = 0
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      if (unstructured_grid%cell_vertices(ivertex,local_id) == 0) exit
      count = count + 1
      ! local vertices must be zero-based for MatCreateMPIAdj; thus subtract 1
      local_vertices(count) = &
        unstructured_grid%cell_vertices(ivertex,local_id) - 1
    enddo
    local_vertex_offset(local_id+1) = count 
  enddo
  num_common_vertices = 3 ! cells must share at least this number of vertices

  ! determine the global offset from 0 for cells on this rank
  global_offset_old = 0
  call MPI_Exscan(num_cells_local_old,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! create an adjacency matrix for calculating the duals (connnections)
#if UGRID_DEBUG
  call printMsg(option,'Adjacency matrix')
#endif

  call MatCreateMPIAdj(option%mycomm,num_cells_local_old, &
                       unstructured_grid%num_vertices_global, &
                       local_vertex_offset, &
                       local_vertices,PETSC_NULL_INTEGER,Adj_mat,ierr)

  ! do not free local_vertices; MatAdjDestroy will do it
  ! do not free local_vertex_offset; MatAdjDestroy will do it

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Adj.out',viewer,ierr)
  call MatView(Adj_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  call printMsg(option,'Dual matrix')
#endif

  ! petsc will call parmetis to calculate the graph/dual
#if defined(PETSC_HAVE_PARMETIS)
  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)
#else
  option%io_buffer = 'Must compile with Parmetis in order to use unstructured grids.'
  call printErrMsg(option)
#endif
  call MatDestroy(Adj_mat,ierr)
  if (allocated(local_vertices)) deallocate(local_vertices)
  if (allocated(local_vertex_offset)) deallocate(local_vertex_offset)
  
#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'Dual.out',viewer,ierr)
  call MatView(Dual_mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

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
  call PetscViewerASCIIOpen(option%mycomm,'is.out',viewer,ierr)
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

  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                      ja_ptr,success,ierr)

  if (.not.success .or. num_rows /= num_cells_local_old) then
    print *, option%myrank, num_rows, success, num_cells_local_old
    option%io_buffer = 'Error getting IJ row indices from dual matrix'
    call printErrMsg(option)
  endif

  ! calculate maximum number of connections for any given cell
  unstructured_grid%max_ndual_per_cell = 0
  do local_id = 1, num_cells_local_old
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) &
      unstructured_grid%max_ndual_per_cell = num_cols
  enddo
  temp_int = unstructured_grid%max_ndual_per_cell
  call MPI_Allreduce(temp_int,unstructured_grid%max_ndual_per_cell, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  
#if UGRID_DEBUG
  write(string,*) unstructured_grid%max_ndual_per_cell
  option%io_buffer = 'Maximum number of duals per cell: ' // adjustl(string)
  call printMsg(option)
#endif
  
  call MatRestoreRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                          ja_ptr,success,ierr)
  
  ! in order to redistributed vertex/cell data among ranks, I package it
  ! in a crude way within a strided petsc vec and pass it.  The stride 
  ! determines the size of each cells "packaged" data 
  vertex_ids_offset = 1 + 1 ! +1 for -777
  dual_offset = vertex_ids_offset + unstructured_grid%max_nvert_per_cell + 1 ! +1 for -888
  stride = dual_offset+ unstructured_grid%max_ndual_per_cell + 1 ! +1 for -999999

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
  ! -999999   ! separator indicating end of information for cell_N
  
  ! the purpose of -777, -888, and -999999 is to allow one to use cells of 
  ! various geometry.  Currently, the max # vertices = 8 and max # duals = 6.
  ! But this will be generalized in the future.
    
  ! create a petsc vec to store all the information for each element
  ! based on the stride calculated above.  
  call VecCreate(option%mycomm,elements_natural,ierr)
  call VecSetSizes(elements_natural, &
                   stride*num_cells_local_new, &
                   PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_natural,ierr)
   
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
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_old_to_new.out', &
                            viewer,ierr)
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! create another strided vector with the old cell/element distribution
  call VecCreate(option%mycomm,elements_old,ierr)
  call VecSetSizes(elements_old,stride*num_cells_local_old,PETSC_DECIDE,ierr)
  call VecSetFromOptions(elements_old,ierr)


  ! 0 = 0-based indexing
  ! MagGetRowIJF90 returns row and column pointers for compressed matrix data
  call MatGetRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                      ja_ptr,success,ierr)

  call VecGetArrayF90(elements_old,vec_ptr,ierr)
  count = 0
  vertex_count = 0
  do local_id = 1, num_cells_local_old
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(global_offset_old+local_id)
    count = count + 1
    ! add the separator
    vec_ptr(count) = -777  ! help differentiate
    ! add the vertex ids
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      count = count + 1
      vertex_count = vertex_count + 1
      ! increment for 1-based ordering
      vec_ptr(count) = unstructured_grid%cell_vertices(ivertex,local_id)
    enddo


    count = count + 1 
    ! another vertex/dual separator
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
    num_cols = iend-istart+1
    if (num_cols > unstructured_grid%max_ndual_per_cell) then
      option%io_buffer = &
        'Number of columns in Dual matrix is larger then max_ndual_per_cell.'
      call printErrMsgByRank(option)
    endif
    do icol = 1, unstructured_grid%max_ndual_per_cell
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
    vec_ptr(count) = -999999  ! help differentiate
  enddo
  call VecRestoreArrayF90(elements_old,vec_ptr,ierr)
  
  call MatRestoreRowIJF90(Dual_mat,0,PETSC_FALSE,PETSC_FALSE,num_rows,ia_ptr, &
                          ja_ptr,success,ierr)
  call MatDestroy(Dual_mat,ierr)
 
#if UGRID_DEBUG
  call printMsg(option,'Before element scatter')
#endif

  ! scatter all the cell data from the old decomposition (as read in in 
  ! parallel) to the more parmetis-calculated decomposition
  call VecScatterCreate(elements_old,PETSC_NULL,elements_natural,is_scatter, &
                        vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_natural, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_natural, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scatter,ierr)

#if UGRID_DEBUG
  call printMsg(option,'After element scatter')
  call PetscViewerASCIIOpen(option%mycomm,'elements_old.out',viewer,ierr)
  call VecView(elements_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  call VecDestroy(elements_old,ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_natural.out',viewer,ierr)
  call VecView(elements_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! update global offset based on new partitioning
  global_offset_new = 0
  call MPI_Exscan(num_cells_local_new,global_offset_new, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  allocate(unstructured_grid%cell_ids_natural(num_cells_local_new))
  unstructured_grid%cell_ids_natural = 0
  
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
    unstructured_grid%cell_ids_natural(local_id) = &
      abs(vec_ptr((local_id-1)*stride+1))
  enddo
  call VecRestoreArrayF90(elements_natural,vec_ptr,ierr)

  ! make a list of petsc ids for each local cell (you simply take the global 
  ! offset and add it to the local contiguous cell ids on each processor
  allocate(int_array(num_cells_local_new))
  do local_id = 1, num_cells_local_new
    int_array(local_id) = local_id+global_offset_new
  enddo
  
  ! make the arrays zero-based
  int_array = int_array - 1
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural - 1
  ! create an application ordering (mapping of natural to petsc ordering)
  call AOCreateBasic(option%mycomm,num_cells_local_new, &
                     unstructured_grid%cell_ids_natural,int_array, &
                     unstructured_grid%ao_natural_to_petsc,ierr)
  deallocate(int_array)
  ! make cell_ids_natural 1-based again
  unstructured_grid%cell_ids_natural = unstructured_grid%cell_ids_natural + 1

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'ao.out',viewer,ierr)
  call AOView(unstructured_grid%ao_natural_to_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! The below creates a list of cells ids for the duals and converts them
  ! to petsc ordering   
  
  ! count the number of cells and their duals  
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    do idual = 1, unstructured_grid%max_ndual_per_cell
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit ! here we hit the 0 at the end of last dual
      count = count + 1
    enddo
  enddo     
               
  ! allocate and fill an array with the natural cell and dual ids
  allocate(int_array(count))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    int_array(count) = unstructured_grid%cell_ids_natural(local_id)
    do idual = 1, unstructured_grid%max_ndual_per_cell
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
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
  call AOApplicationToPetsc(unstructured_grid%ao_natural_to_petsc,count, &
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
  allocate(unstructured_grid%cell_ids_petsc(num_cells_local_new))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    ! extract the petsc id for the cell
    unstructured_grid%cell_ids_petsc(local_id) = int_array(count)
    ! store it in the elements_petsc vector too
    vec_ptr((local_id-1)*stride+1) = int_array(count)
    do idual = 1, unstructured_grid%max_ndual_per_cell
!geh      dual_id = vec_ptr2(idual + dual_offset + (local_id-1)*stride)
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
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
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc.out',viewer,ierr)
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
    do idual = 1, unstructured_grid%max_ndual_per_cell
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
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
  call PetscViewerASCIIOpen(option%mycomm, &
                            'elements_local_dual_unsorted.out',viewer,ierr)
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
    allocate(unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count))
    unstructured_grid%ghost_cell_ids_petsc = 0

    unstructured_grid%ghost_cell_ids_petsc(1:ghost_cell_count) = &
      int_array3(1:ghost_cell_count)

#if UGRID_DEBUG
  call printMsg(option,'  Remappping ghost ids')
#endif

    ! remap of duals of ghost cells
    call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
    do local_id=1, num_cells_local_new
      do idual = 1, unstructured_grid%max_ndual_per_cell
        ! dual_id is now the negative of the local unsorted ghost cell id
        dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
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

    deallocate(int_array_pointer)
    nullify(int_array_pointer)
    deallocate(int_array2)
    deallocate(int_array3)
    deallocate(int_array4)
    deallocate(int_array5)
  endif

  unstructured_grid%nlmax = num_cells_local_new
  unstructured_grid%num_ghost_cells = ghost_cell_count
  unstructured_grid%ngmax = &
    num_cells_local_new + ghost_cell_count


#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_local_dual.out', &
                            viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

#if UGRID_DEBUG
  call printMsg(option,'Resizing natural cell id array')
#endif

  ! Resize cell_ids_natural to include ghosted cells
  allocate(int_array(unstructured_grid%nlmax))
  int_array(:) = unstructured_grid%cell_ids_natural(:)
  deallocate(unstructured_grid%cell_ids_natural)
  allocate(unstructured_grid%cell_ids_natural(unstructured_grid%ngmax))
  unstructured_grid%cell_ids_natural(:) = -999
  unstructured_grid%cell_ids_natural(1:unstructured_grid%nlmax) = int_array(:)
  deallocate(int_array)
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  do local_id=1, unstructured_grid%nlmax
    do idual = 1, unstructured_grid%max_ndual_per_cell
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit
      if (dual_id > unstructured_grid%nlmax) then
        unstructured_grid%cell_ids_natural(dual_id) = &
          vec_ptr2(idual + dual_offset + (local_id-1)*stride)
      endif       
    enddo
  enddo
  if (minval(unstructured_grid%cell_ids_natural) < 1) then
    write(string,*) minval( unstructured_grid%cell_ids_natural)
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
  allocate(unstructured_grid%cell_neighbors_local_ghosted( &
           0:unstructured_grid%max_ndual_per_cell,unstructured_grid%nlmax))
  unstructured_grid%cell_neighbors_local_ghosted = 0
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do local_id=1, unstructured_grid%nlmax
    count = 0
    do idual = 1, unstructured_grid%max_ndual_per_cell
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      ! flag ghosted cells in dual as negative
      !geh: these negative dual ids are used later in UGridDMCreateJacobian() 
      !     to specify off processor connectity in the Jacobian
      if (dual_id > unstructured_grid%nlmax) dual_id = -dual_id
      unstructured_grid%cell_neighbors_local_ghosted(idual,local_id) = dual_id
    enddo
    ! set the # of duals in for the cell
    unstructured_grid%cell_neighbors_local_ghosted(0,local_id) = count
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! need to create a local ghosted vector in which we can collect element info
  ! including ghost cells
  call VecCreateSeq(PETSC_COMM_SELF,stride*unstructured_grid%ngmax, &
                    elements_local,ierr)
  call VecSetBlockSize(elements_local,stride,ierr)
  allocate(int_array(unstructured_grid%ngmax))
  int_array(1:unstructured_grid%nlmax) = &
    unstructured_grid%cell_ids_petsc(:)
  if (unstructured_grid%num_ghost_cells > 0) then
    int_array(unstructured_grid%nlmax+1:unstructured_grid%ngmax) = &
      unstructured_grid%ghost_cell_ids_petsc(:)
  endif
  int_array = int_array-1
  call ISCreateBlock(option%mycomm,stride,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,is_scatter,ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    int_array(ghosted_id) = ghosted_id-1
  enddo
  call ISCreateBlock(option%mycomm,stride,unstructured_grid%ngmax, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_local_to_ghost.out',viewer,ierr)
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'is_gather_elem_local_to_ghost.out',viewer,ierr)
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
  call printMsg(option,'Scatter/gathering local ghosted vertices')
#endif

  ! make a list of local vertices
  max_int_count = 2*unstructured_grid%ngmax
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
  ! note that the vertices are still in natural numbering
  call VecGetArrayF90(elements_local,vec_ptr,ierr)
  do local_id=1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride)
      if (vertex_id < 1) exit
      vertex_count = vertex_count + 1
      if (vertex_count > max_int_count) then
        call reallocateIntArray(int_array_pointer,max_int_count)
      endif
      vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride) = vertex_count
      int_array_pointer(vertex_count) = vertex_id
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr)

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

  allocate(unstructured_grid%vertex_ids_natural(vertex_count))
  unstructured_grid%vertex_ids_natural = int_array3(1:vertex_count)

  ! now load all the vertices needed to define all the local cells
  ! on the processor
  allocate(needed_vertices_petsc(vertex_count))
  needed_vertices_petsc(1:vertex_count) = int_array3(1:vertex_count)

  ! allocate the array that will store the vertex ids for each cell.
  ! remember that max_nvert_per_cell is the max # of vertices in a cell
  ! currently hardwired to 8.
  deallocate(unstructured_grid%cell_vertices)
  allocate(unstructured_grid%cell_vertices( &
             0:unstructured_grid%max_nvert_per_cell,unstructured_grid%ngmax))
  unstructured_grid%cell_vertices = 0
  
  ! permute the local ids calculated earlier in the int_array4
  call VecGetArrayF90(elements_local,vec_ptr,ierr)
  do ghosted_id=1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%max_nvert_per_cell
      ! extract the original vertex id
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride)
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices(0,ghosted_id)+1
      unstructured_grid%cell_vertices(count,ghosted_id) = &
        int_array4(vertex_id)
      unstructured_grid%cell_vertices(0,ghosted_id) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (ghosted_id-1)*stride) = &
        int_array4(vertex_id)
    enddo
  enddo
  call VecRestoreArrayF90(elements_local,vec_ptr,ierr)
  deallocate(int_array2)
  deallocate(int_array3)
  deallocate(int_array4)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'elements_vert_local' // trim(adjustl(string)) // '.out'
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call VecView(elements_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  call VecDestroy(elements_local,ierr)

  ! now we need to work on aligning the original vertex coordinates with 
  ! the current ordering or permuted/rearranged ordering.

  ! IS for gather operation - need local numbering
  allocate(int_array(vertex_count))
  ! vertex_count = # of local vertices (I believe ghosted+non-ghosted)
  do ivertex = 1, vertex_count
    int_array(ivertex) = ivertex-1
  enddo

  ! include cell ids (use block ids, not indices)
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_gather,ierr)
  deallocate(int_array)

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

  ! IS for scatter - provide petsc global numbering
  allocate(int_array(vertex_count))
  do ivertex = 1, vertex_count
    int_array(ivertex) = (needed_vertices_petsc(ivertex)-1)
  enddo
  ! include cell ids
  call ISCreateBlock(option%mycomm,3,vertex_count, &
                     int_array,PETSC_COPY_VALUES,is_scatter,ierr)
  deallocate(int_array)

  ! resize vertex array to new size
  unstructured_grid%num_vertices_local = vertex_count
  allocate(unstructured_grid%vertices(vertex_count))
  do ivertex = 1, vertex_count
    unstructured_grid%vertices(ivertex)%x = 0.d0
    unstructured_grid%vertices(ivertex)%y = 0.d0
    unstructured_grid%vertices(ivertex)%z = 0.d0
  enddo

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter_vert_old_to_new.out',viewer,ierr)
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
  call PetscViewerASCIIOpen(option%mycomm,'is_gather_vert_old_to_new.out',viewer,ierr)
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

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'vertex_coord_old.out',viewer,ierr)
  call VecView(vertices_old,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecDestroy(vertices_old,ierr)


  call VecGetArrayF90(vertices_new,vec_ptr,ierr)
  do ivertex = 1, unstructured_grid%num_vertices_local
    unstructured_grid%vertices(ivertex)%id = needed_vertices_petsc(ivertex)
    unstructured_grid%vertices(ivertex)%x = vec_ptr((ivertex-1)*3+1)
    unstructured_grid%vertices(ivertex)%y = vec_ptr((ivertex-1)*3+2)
    unstructured_grid%vertices(ivertex)%z = vec_ptr((ivertex-1)*3+3)
  enddo
  call VecRestoreArrayF90(vertices_new,vec_ptr,ierr)
  
#if UGRID_DEBUG
  write(string,*) option%myrank
  string = 'vertex_coord_new' // trim(adjustl(string)) // '.out'
  call PetscViewerASCIIOpen(PETSC_COMM_SELF,trim(string),viewer,ierr)
  call VecView(vertices_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecDestroy(vertices_new,ierr)

#if UGRID_DEBUG
  call printMsg(option,'Setting cell types')
#endif

  allocate(unstructured_grid%cell_type(unstructured_grid%ngmax))
  do ghosted_id = 1, unstructured_grid%ngmax
    ! Determine number of faces and cell-type of the current cell
    select case(unstructured_grid%cell_vertices(0,ghosted_id))
      case(8)
        unstructured_grid%cell_type(ghosted_id) = HEX_TYPE
      case(6)
        unstructured_grid%cell_type(ghosted_id) = WEDGE_TYPE
      case(5)
        unstructured_grid%cell_type(ghosted_id) = PYR_TYPE
      case(4)
        unstructured_grid%cell_type(ghosted_id) = TET_TYPE
      case default
        option%io_buffer = 'Cell type not recognized: '
        call printErrMsg(option)
    end select
  enddo
  
  unstructured_grid%global_offset = global_offset_new  

end subroutine UGridDecompose

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

  interface
  end interface
  
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
  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
                    PETSC_DETERMINE,ugdm%global_vec,ierr)
  call VecSetBlockSize(ugdm%global_vec,ndof,ierr)
  ! create local vec
  call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax*ndof, &
                    ugdm%local_vec,ierr)
  call VecSetBlockSize(ugdm%local_vec,ndof,ierr)
  
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

  call VecCreateMPI(option%mycomm,unstructured_grid%nlmax*ndof, &
                    PETSC_DETERMINE,vec_tmp,ierr)
  call VecSetBlockSize(vec_tmp,ndof,ierr)
  call VecScatterCreate(ugdm%global_vec,ugdm%is_local_petsc,vec_tmp, &
                        ugdm%is_local_natural,ugdm%scatter_gton,ierr)
  call VecDestroy(vec_tmp,ierr)

#if UGRID_DEBUG
  string = 'scatter_gton' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecScatterView(ugdm%scatter_gton,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

end subroutine UGridCreateUGDM

! ************************************************************************** !
!
! UGridComputeInternConnect: computes internal connectivity of an
!                            unstructured grid
! author: Glenn Hammond
! date: 10/21/09
!
! ************************************************************************** !
function UGridComputeInternConnect(unstructured_grid,grid_x,grid_y,grid_z, &
                                   scatter_ltol,option)

  use Connection_module
  use Option_module
  use Utility_module, only : DotProduct, CrossProduct

  implicit none

  type(connection_set_type), pointer :: UGridComputeInternConnect
  type(option_type) :: option
  PetscReal :: grid_x(*), grid_y(*), grid_z(*)
  type(unstructured_grid_type) :: unstructured_grid
  VecScatter :: scatter_ltol 

  type(connection_set_type), pointer :: connections
  PetscInt :: nconn, iconn
  PetscInt :: idual, dual_id

  PetscInt, allocatable :: face_to_vertex(:,:)
  PetscInt, allocatable :: cell_to_face(:,:)
  PetscInt, allocatable :: face_to_cell(:,:)
  PetscInt, allocatable :: vertex_to_cell(:,:)
  PetscInt, allocatable :: temp_int(:)
  PetscInt, allocatable :: temp_int_2d(:,:)
  PetscBool, allocatable :: local_boundary_face(:)
  PetscInt :: num_match
  PetscInt :: found_count
  PetscBool :: found
  PetscBool :: match_found
  PetscInt :: face_count
  PetscInt :: count, i
  PetscInt :: iface, iface2, iside
  PetscInt :: face_id, face_id2
  PetscInt :: ghosted_id, ghosted_id2
  PetscInt :: local_id, local_id2
  PetscInt :: cell_id, cell_id2
  PetscInt :: dual_local_id
  PetscInt :: ivertex, ivertex2
  PetscInt :: vertex_id, vertex_id2
  PetscInt :: vertex_ids4(4)
  PetscInt :: nfaces, nfaces2, nvertices, nvertices2, cell_type, cell_type2
  PetscInt :: face_type, face_type2
  PetscBool:: face_found, vertex_found
  
  PetscReal :: v1(3), v2(3), v3(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: vcross(3), magnitude
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  PetscInt :: ivert
  
  type(plane_type) :: plane1, plane2
  type(point_type) :: point1, point2, point3, point4
  type(point_type) :: point_up, point_dn
  type(point_type) :: intercept1, intercept2, intercept

  character(len=MAXSTRINGLENGTH) :: string  
  
  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL* &
           unstructured_grid%ngmax))
  face_to_vertex = 0
  allocate(cell_to_face(MAX_FACE_PER_CELL, &
                        unstructured_grid%ngmax))
  cell_to_face = 0
  allocate(face_to_cell(2,MAX_FACE_PER_CELL* &
                        unstructured_grid%ngmax))
  face_to_cell = 0
  allocate(vertex_to_cell(0:MAX_CELLS_SHARING_A_VERTEX, &
                          unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  allocate(unstructured_grid%face_to_vertex_natural(MAX_VERT_PER_FACE, &
           MAX_FACE_PER_CELL*unstructured_grid%ngmax))
  unstructured_grid%face_to_vertex_natural = 0

  face_count = 0
  do ghosted_id = 1, unstructured_grid%ngmax
    cell_type = unstructured_grid%cell_type(ghosted_id)
    nfaces = UCellGetNFaces(cell_type)
    do iface = 1, nfaces
      face_count = face_count + 1
      cell_to_face(iface,ghosted_id) = face_count
      face_to_cell(1,face_count) = ghosted_id
      call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                      vertex_ids4)
      do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices(vertex_ids4(ivertex),ghosted_id)
          if (face_to_vertex(ivertex,face_count) > 0) then
            unstructured_grid%face_to_vertex_natural(ivertex,face_count) = &
              unstructured_grid%vertex_ids_natural(face_to_vertex(ivertex,face_count))
          endif
      enddo
    enddo
  enddo

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
  
  do local_id = 1, unstructured_grid%nlmax
    ! Selet a cell and find number of vertices
    cell_id = local_id
    ! cell_type is ghosted, but local cells are in the first nlmax entries
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type)
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      ! Select a neighboring cell
      ! ghosted neighbors have a negative id
      cell_id2 = &
        abs(unstructured_grid%cell_neighbors_local_ghosted(idual,local_id))
      cell_type2 = unstructured_grid%cell_type(cell_id2)
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      ! Find the number of vertices for neighboring cell
      nfaces2 = UCellGetNFaces(cell_type2)
      ! Initialize
      face_found = PETSC_FALSE
      do iface = 1, nfaces
        ! Select a face and find number of vertices forming the face
        face_id = cell_to_face(iface,cell_id)
        nvertices = UCellGetNFaceVertices(cell_type,iface)
        do ivertex = 1, nvertices
          ! Select a vertex and initialize vertex_found
          vertex_id = face_to_vertex(ivertex,face_id) ! face_to_vertex is 1-based indexing
          vertex_found = PETSC_FALSE
          do ivertex2 = 1, unstructured_grid%cell_vertices(0,cell_id2)
            vertex_id2 = unstructured_grid%cell_vertices(ivertex2,cell_id2) 
            if (vertex_id == vertex_id2) then
              vertex_found = PETSC_TRUE
              exit
            endif
          enddo
          !
          ! If ivertex of iface of the Cell is not present as vertices of the
          ! Neighboring-Cell, then iface is not the shared face. Skip iterating
          ! over the remaing vertices of iface
          if (.not.vertex_found) exit
        enddo
        
        if (vertex_found) then
          ! All the vertices of iface are present in the Neighboring cells.
          ! Thus, iface is the shared face.
          face_found = PETSC_TRUE
          
          ! Now, we have to find iface2 that corresponds to iface
          do iface2 = 1, nfaces2
            face_id2 = cell_to_face(iface2,cell_id2)
            !geh nvertices2 = 4
            !gehcomment: I believe that cell_type and iface on next line shoudl be the "2" versions
            !geh if ((cell_type == WEDGE_TYPE).and.(iface.gt.3)) nvertices2 = 3
            nvertices2 = UCellGetNFaceVertices(cell_type2,iface2)
            ! Both iface and iface2 need to have same number of vertices
            if (nvertices == nvertices2) then
              ! Count the number of vertices of iface which match vertices
              ! of iface2
              num_match = 0
              do ivertex = 1,nvertices
                vertex_id = face_to_vertex(ivertex,face_id)
                vertex_found = PETSC_FALSE ! gehbug - used to be PETSC_TRUE
                
                do ivertex2 = 1, nvertices2 ! gehbug - used to be nvertices
                  vertex_id2 = face_to_vertex(ivertex2,face_id2)
                  if (vertex_id == vertex_id2) then
                    vertex_found = PETSC_TRUE 
                    num_match = num_match + 1
                    vertex_ids4(num_match) = vertex_id
                    exit
                  endif
                enddo
                !
                ! If vertex_id of face_id not found as one of vertices of face_id2,
                ! face_id2 is not shared between cells
                if (.not.vertex_found) exit
              enddo
              if (num_match == nvertices) then
                ! remove duplicate face
                !geh: I believe that face_id2 will always be removed
                if (face_id2 > face_id) then
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id2, ' -> ', face_id
                  option%io_buffer = 'Duplicated face removed:' // string
                  call printMsg(option)
#endif
                  cell_to_face(iface2,cell_id2) = face_id
                  ! flag face_id2 as removed
                  face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
                  ! add cell_id2 to face_ids list
                  face_to_cell(2,face_id) = cell_id2
                else
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id, ' -> ', face_id2
                  option%io_buffer = 'Duplicated face removed:' // string
                  call printMsg(option)
#endif
                  cell_to_face(iface,cell_id) = face_id2
                  ! flag face_id as removed  
                  face_to_cell(1,face_id) = -face_to_cell(1,face_id)
                  ! add cell_id to face_ids2 list
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
      if (.not.face_found) then
        write(string,*) option%myrank
        string = '(' // trim(adjustl(string)) // ')'
        write(*,'(a,'' local_id = '',i3,'' natural_id = '',i3,/ &
                  &''  vertices: '',8i3)') &
                   trim(string), &
                   cell_id,unstructured_grid%cell_ids_natural(cell_id), &
                   (unstructured_grid%vertex_ids_natural( &
                     unstructured_grid%cell_vertices(ivertex,cell_id)), &
                     ivertex=1,unstructured_grid%cell_vertices(0,cell_id))
        write(*,'(a,'' local_id2 = '',i3,'' natural_id2 = '',i3,/ &
                       &''  vertices2: '',8i3 &
                       &)') &
                   trim(string), &
                   cell_id2,unstructured_grid%cell_ids_natural(cell_id2), &
                   (unstructured_grid%vertex_ids_natural( &
                     unstructured_grid%cell_vertices(ivertex2,cell_id2)), &
                     ivertex2=1,unstructured_grid%cell_vertices(0,cell_id2))
        option%io_buffer='No shared face found.'
        call printErrMsgByRank(option)
      endif
    enddo ! idual-loop
  enddo  ! local_id-loop

  ! count up the # of faces
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) &
      face_count = face_count + 1
  enddo
  allocate(unstructured_grid%face_to_vertex(MAX_VERT_PER_FACE,face_count))
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      unstructured_grid%face_to_vertex(:,face_count) = face_to_vertex(:,iface)
    endif
  enddo
  deallocate(face_to_vertex)
  ! reallocate face_to_cell to proper size
  allocate(temp_int_2d(2,face_count))
  allocate(temp_int(size(face_to_cell,2)))
  temp_int = 0
  face_count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      face_count = face_count + 1
      temp_int_2d(:,face_count) = face_to_cell(:,iface)
      temp_int(iface) = face_count
    endif
  enddo
  deallocate(face_to_cell)
  allocate(face_to_cell(2,face_count))
  face_to_cell = temp_int_2d
  deallocate(temp_int_2d)

  
  ! remap faces in cells using temp_int from above
  do iface = 1, size(face_to_cell,2)
    face_id = iface
    do i = 1,2
      cell_id = face_to_cell(i,face_id)
      ! check for exterior face
      if (cell_id < 1) cycle
      found = PETSC_FALSE
      cell_type = unstructured_grid%cell_type(cell_id)
      nfaces = UCellGetNFaces(cell_type)
      do iface2 = 1, nfaces
        face_id2 = cell_to_face(iface2,cell_id)
        if (face_id < 0) cycle
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
  
  do ghosted_id = 1, unstructured_grid%ngmax
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      if ( vertex_id <= 0) cycle 
      count = vertex_to_cell(0,vertex_id) + 1
      if (count > MAX_CELLS_SHARING_A_VERTEX) then
        write(string,*) 'Vertex can be shared by at most by ',MAX_CELLS_SHARING_A_VERTEX, &
              ' cells. Rank = ', option%myrank, ' vertex_id = ', vertex_id, ' exceeds it.'
        option%io_buffer = string
        call printErrMsg(option)
      endif
      vertex_to_cell(count,vertex_id) = ghosted_id
      vertex_to_cell(0,vertex_id) = count
    enddo
  enddo
  
  nconn = 0
  do local_id = 1, unstructured_grid%nlmax
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_id = unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      ! count all ghosted connections (dual_id < 0)
      ! only count connection with cells of larger ids to avoid double counts
!geh: we need to cound all local connection, but just once (local_id < dual_id) and all
!      ghosted connections (dual_id < 0)
      if (dual_id < 0 .or. local_id < dual_id) then
!geh: Nope      if (dual_id > 0 .and. local_id < dual_id) then !sp 
        nconn = nconn + 1
      endif
    enddo
  enddo


  connections => ConnectionCreate(nconn,option%nphase,INTERNAL_CONNECTION_TYPE)
  
  allocate(unstructured_grid%face_area(face_count))
  allocate(unstructured_grid%connection_to_face(nconn))
  unstructured_grid%connection_to_face = 0

  ! loop over connection again
  iconn = 0
  do local_id = 1, unstructured_grid%nlmax
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_local_id = &
        unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      ! abs(dual_local_id) to accommodate connections to ghost cells where 
      ! the dual id is < 0.
      if (local_id < abs(dual_local_id)) then 
        iconn = iconn + 1
        ! find face
        found = PETSC_FALSE
        do iface = 1, unstructured_grid%cell_vertices(0,local_id)
          face_id = cell_to_face(iface,local_id)
          do iside = 1,2
            cell_id2 = face_to_cell(iside,face_id)
            if (cell_id2 == abs(dual_local_id)) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (found) exit
        enddo
        if (found) then
          unstructured_grid%connection_to_face(iconn) = face_id
        else
          write(string,*) option%myrank,local_id,dual_local_id 
          option%io_buffer = 'face not found in connection loop' // string 
          call printErrMsg(option)
        endif
        face_type = &
          UCellGetFaceType(unstructured_grid%cell_type(local_id),iface)
        found = PETSC_FALSE
        do iface2 = 1, unstructured_grid%cell_vertices(0,cell_id2)
          if (cell_to_face(iface,local_id) == &
              cell_to_face(iface2,cell_id2)) then
            found = PETSC_TRUE
            exit
          endif
        enddo
        if (found) then
          face_type2 = &
            UCellGetFaceType(unstructured_grid%cell_type(cell_id2), &
                                                                 iface2)
          if (face_type /= face_type2) then
            write(string,*) option%myrank, local_id, cell_id2 
            option%io_buffer = 'face types do not match' // string 
            call printErrMsg(option)
          endif
        else
          write(string,*) option%myrank, iface, cell_id2
          option%io_buffer = 'global face not found' // string 
          call printErrMsg(option)
        endif
        connections%id_up(iconn) = local_id
        connections%id_dn(iconn) = abs(dual_local_id)
        ! need to add the surface areas, distance, etc.
        point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
        point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
        point3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
        if (face_type == QUAD_FACE_TYPE) then
          point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(4,face_id))
        endif
        
        call UCellComputePlane(plane1,point1,point2,point3)
       
        point_up%x = grid_x(local_id)
        point_up%y = grid_y(local_id)
        point_up%z = grid_z(local_id)
        point_dn%x = grid_x(abs(dual_local_id))
        point_dn%y = grid_y(abs(dual_local_id))
        point_dn%z = grid_z(abs(dual_local_id))
        v1(1) = point_dn%x-point_up%x
        v1(2) = point_dn%y-point_up%y
        v1(3) = point_dn%z-point_up%z
        n_up_dn = v1 / sqrt(DotProduct(v1,v1))
        call UCellGetPlaneIntercept(plane1,point_up,point_dn,intercept1)
        
        v1(1) = point3%x-point2%x
        v1(2) = point3%y-point2%y
        v1(3) = point3%z-point2%z
        v2(1) = point1%x-point2%x
        v2(2) = point1%y-point2%y
        v2(3) = point1%z-point2%z
        !geh: area = 0.5 * |v1 x v2|
        vcross = CrossProduct(v1,v2)
        !geh: but then we have to project the area onto the vector between
        !     the cell centers (n_up_dn)
        magnitude = sqrt(DotProduct(vcross,vcross))
        n1 = vcross/magnitude
        area1 = 0.5d0*magnitude
        area1 = dabs(area1*DotProduct(n1,n_up_dn))
        !geh: The below does not project onto the vector between cell centers.
        !gehbug area1 = 0.5d0*sqrt(DotProduct(n1,n1))
        
        if (face_type == QUAD_FACE_TYPE) then
          call UCellComputePlane(plane2,point3,point4,point1)
          call UCellGetPlaneIntercept(plane2,point_up,point_dn,intercept2)
          v1(1) = point1%x-point4%x
          v1(2) = point1%y-point4%y
          v1(3) = point1%z-point4%z
          v2(1) = point3%x-point4%x
          v2(2) = point3%y-point4%y
          v2(3) = point3%z-point4%z
          magnitude = sqrt(DotProduct(vcross,vcross))
          n2 = vcross/magnitude
          area2 = 0.5d0*magnitude
          area2 = dabs(area2*DotProduct(n2,n_up_dn))
        else 
          area2 = 0.d0
        endif

        if (face_type == QUAD_FACE_TYPE) then
          intercept%x = 0.5d0*(intercept1%x + intercept2%x)
          intercept%y = 0.5d0*(intercept1%y + intercept2%y)
          intercept%z = 0.5d0*(intercept1%z + intercept2%z)
        else
          intercept%x = intercept1%x
          intercept%y = intercept1%y
          intercept%z = intercept1%z
        endif
        
        !geh: this is very crude, but for now use average location of intercept
        v1(1) = intercept%x-point_up%x
        v1(2) = intercept%y-point_up%y
        v1(3) = intercept%z-point_up%z
        v2(1) = point_dn%x-intercept%x
        v2(2) = point_dn%y-intercept%y
        v2(3) = point_dn%z-intercept%z
        dist_up = sqrt(DotProduct(v1,v1))
        dist_dn = sqrt(DotProduct(v2,v2))
        
        connections%dist(-1:3,iconn) = 0.d0
        connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connections%dist(0,iconn) = dist_up + dist_dn
        v3 = v1 + v2
        connections%dist(1:3,iconn) = v3/sqrt(DotProduct(v3,v3))
        connections%area(iconn) = area1 + area2
       
      endif
    enddo
  enddo
  
  ! Save area and centroid of faces
  allocate(unstructured_grid%face_centroid(face_count))
  do iface = 1,face_count
    unstructured_grid%face_centroid(iface)%id = 0
  enddo
  
  do local_id = 1, unstructured_grid%nlmax
    do iface = 1,MAX_FACE_PER_CELL
      face_id = cell_to_face(iface, local_id)
      if (face_id == 0) cycle
      if ( unstructured_grid%face_centroid(face_id)%id == 0) then
        count = 0
        unstructured_grid%face_centroid(face_id)%x = 0.d0
        unstructured_grid%face_centroid(face_id)%y = 0.d0
        unstructured_grid%face_centroid(face_id)%z = 0.d0
        if (unstructured_grid%face_to_vertex(4,face_id) == 0) then
          face_type = TRI_FACE_TYPE
        else
          face_type = QUAD_FACE_TYPE
        endif
        point1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
        point2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
        point3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
        if (face_type == QUAD_FACE_TYPE) then
          point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(4,face_id))
        else
          point4 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
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
        area2 = 0.5d0*sqrt(DotProduct(n2,n2))
        unstructured_grid%face_area(face_id) = area1 + area2
        
        do ivert = 1,MAX_VERT_PER_FACE
          vertex_id = unstructured_grid%face_to_vertex(ivert,face_id)
          if (vertex_id.ne.0) then
            unstructured_grid%face_centroid(face_id)%x = &
              unstructured_grid%face_centroid(face_id)%x + &
              unstructured_grid%vertices(vertex_id)%x
            unstructured_grid%face_centroid(face_id)%y = &
              unstructured_grid%face_centroid(face_id)%y + &
              unstructured_grid%vertices(vertex_id)%y
            unstructured_grid%face_centroid(face_id)%z = &
              unstructured_grid%face_centroid(face_id)%z + &
              unstructured_grid%vertices(vertex_id)%z
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

  allocate(unstructured_grid%face_to_cell_ghosted(size(face_to_cell,1), &
                                                  size(face_to_cell,2)))
  unstructured_grid%face_to_cell_ghosted = face_to_cell
  allocate(unstructured_grid%cell_to_face_ghosted(size(cell_to_face,1), &
                                                  size(cell_to_face,2)))
  unstructured_grid%cell_to_face_ghosted(:,:) = cell_to_face(:,:)

  deallocate(cell_to_face)
  deallocate(face_to_cell)
  deallocate(vertex_to_cell)

  UGridComputeInternConnect => connections

end function UGridComputeInternConnect

! ************************************************************************** !
!
! UGridPopulateConnection: Computes details of connection (area, dist, etc)
! author: Gautam Bisht
! date: 10/30/09
!
! ************************************************************************** !
subroutine UGridPopulateConnection(unstructured_grid, connection, iface_cell, &
                                   iconn, ghosted_id, option)

  use Connection_module
  use Utility_module, only : DotProduct
  use Option_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface_cell
  PetscInt :: iconn
  PetscInt :: ghosted_id
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  PetscInt  :: face_id
  PetscInt  :: ivert,vert_id
  PetscReal :: v1(3),v2(3),n_dist(3), dist
  type(point_type) :: vertex_8(8)
  type(plane_type) :: plane
  type(point_type) :: point, vertex1, vertex2, vertex3, intercept
  character(len=MAXWORDLENGTH) :: word
  
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      if (iface_cell == 0) then
        write(word,*) ghosted_id
        option%io_buffer = 'Face id undefined for cell ' // &
          trim(adjustl(word)) // &
          ' in boundary condition.  Should this be a source/sink?'
        call printErrMsgByRank(option)
      endif
      ! Compute cell centeroid
      v2 = 0.d0
      do ivert = 1, unstructured_grid%cell_vertices(0, ghosted_id)
        vert_id = unstructured_grid%cell_vertices(ivert, ghosted_id)
        vertex_8(ivert)%x = unstructured_grid%vertices(vert_id)%x
        vertex_8(ivert)%y = unstructured_grid%vertices(vert_id)%y
        vertex_8(ivert)%z = unstructured_grid%vertices(vert_id)%z
      enddo
      v2 = UCellComputeCentroid(unstructured_grid%cell_type(ghosted_id), &
                                vertex_8)
! Instead of connecting centroid with face center, calculate the shortest
! distance between the centroid and face and use that distance - geh
#if 0
      ! Get face-centroid vector
      face_id = unstructured_grid%cell_to_face_ghosted(iface_cell, ghosted_id)
      v1(1) = unstructured_grid%face_centroid(face_id)%x
      v1(2) = unstructured_grid%face_centroid(face_id)%y
      v1(3) = unstructured_grid%face_centroid(face_id)%z
      
#endif
      !TODO(geh): add support for a quad face
      !TODO(geh): replace %face_to_vertex array with function that returns vertices
      !           based on cell type and iface
      face_id = unstructured_grid%cell_to_face_ghosted(iface_cell, ghosted_id)
      vertex1 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(1,face_id))
      vertex2 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(2,face_id))
      vertex3 = unstructured_grid%vertices(unstructured_grid%face_to_vertex(3,face_id))
      call UCellComputePlane(plane,vertex1,vertex2,vertex3)
      point%x = v2(1)
      point%y = v2(2)
      point%z = v2(3)
      call UCellProjectPointOntoPlane(plane,point,intercept)
      
      ! Compute distance vector: cell_center - face_centroid
      v1(1) = v2(1) - intercept%x
      v1(2) = v2(2) - intercept%y
      v1(3) = v2(3) - intercept%z
      
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

  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal :: centroid(3)
  PetscErrorCode :: ierr 

  do ghosted_id = 1, unstructured_grid%ngmax 
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    centroid = UCellComputeCentroid(unstructured_grid%cell_type(ghosted_id), &
                                    vertex_8)
    grid_x(ghosted_id) = centroid(1)
    grid_y(ghosted_id) = centroid(2)
    grid_z(ghosted_id) = centroid(3)
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
subroutine UGridComputeVolumes(unstructured_grid,option,volume)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  Vec :: volume
  

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal, pointer :: volume_p(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(volume,volume_p,ierr)

  do local_id = 1, unstructured_grid%nlmax
    ! ghosted_id = local_id on unstructured grids
    ghosted_id = local_id
    do ivertex = 1, unstructured_grid%cell_vertices(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices(ivertex,ghosted_id)
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
    volume_p(local_id) = UCellComputeVolume(unstructured_grid%cell_type( &
                           ghosted_id),vertex_8)
  enddo
      
  call VecRestoreArrayF90(volume,volume_p,ierr)

end subroutine UGridComputeVolumes

! ************************************************************************** !
!
! UGridEnsureRightHandRule: Rearranges order of vertices within each cell
!                           so that when the right hand rule is appied to a
!                           face, the thumb points away from teh centroid
! author: Glenn Hammond
! date: 10/24/11
!
! ************************************************************************** !
subroutine UGridEnsureRightHandRule(unstructured_grid,x,y,z,nL2A,option)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  PetscReal :: x(:), y(:), z(:)
  PetscInt :: nL2A(:)
  type(option_type) :: option

  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(point_type) :: point, point1, point2, point3
  type(plane_type) :: plane1
  PetscReal :: distance
  PetscInt :: cell_vertex_ids_before(8), cell_vertex_ids_after(8)
  PetscInt :: face_vertex_ids(4)
  PetscInt :: num_vertices, iface, cell_type, num_faces, face_type, i

  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id
    cell_type = unstructured_grid%cell_type(local_id)
    num_vertices = UCellGetNVertices(cell_type)
    cell_vertex_ids_before(1:num_vertices) = &
      unstructured_grid%cell_vertices(1:num_vertices,ghosted_id)
    cell_vertex_ids_after = cell_vertex_ids_before
    ! point is the centroid of cell
    point%x = x(ghosted_id)
    point%y = y(ghosted_id)
    point%z = z(ghosted_id)
    num_faces = UCellGetNFaces(cell_type)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface)
      call UCellGetFaceVertices(option,cell_type,iface,face_vertex_ids)
      ! Only need to use the first three vertices to produce a plane.  Will
      ! assume that the plane represents the entire face
      point1 = &
        unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(1)))
      point2 = &
        unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(2)))
      point3 = &
        unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(3)))
      call UCellComputePlane(plane1,point1,point2,point3)
      distance = UCellComputeDistanceFromPlane(plane1,point) 
      if (distance > 0.d0) then
        ! need to swap so that distance is negative (point lies below plane)
        write(option%io_buffer,'(''Cell '',i3,'' of type "'',a, &
              & ''" with vertices: '',8i6, &
              & '' violates right hand rule at face "'',a, &
              & ''" for face vertices '', &
              & '' based on face vertices '',3i3)') &
              nL2A(local_id), &
              trim(UCellTypeToWord(cell_type)), &
          (unstructured_grid%vertex_ids_natural(cell_vertex_ids_before(i)), &
             i = 1,8), &
              trim(UCellFaceTypeToWord(face_type)), &
          (face_vertex_ids(i),i = 1,3) 
        call printErrMsgByRank(option)
      endif
    enddo
  enddo

end subroutine UGridEnsureRightHandRule

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
  PetscInt :: ndof_local
  PetscErrorCode :: ierr
  
  allocate(d_nnz(unstructured_grid%nlmax))
  allocate(o_nnz(unstructured_grid%nlmax))
  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0
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

  ndof_local = unstructured_grid%nlmax*ugdm%ndof
!  if (option%mycommsize > 1) then
    select case(mat_type)
      case(MATAIJ)
        d_nnz = d_nnz*ugdm%ndof
        o_nnz = o_nnz*ugdm%ndof
        call MatCreateMPIAIJ(option%mycomm,ndof_local,ndof_local, &
                             PETSC_DETERMINE,PETSC_DETERMINE, &
                             PETSC_NULL_INTEGER,d_nnz, &
                             PETSC_NULL_INTEGER,o_nnz,J,ierr)
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog, &
                                        ugdm%mapping_ltog,ierr)
        call MatSetLocalToGlobalMappingBlock(J,ugdm%mapping_ltogb, &
                                             ugdm%mapping_ltogb,ierr)
      case(MATBAIJ)
        call MatCreateMPIBAIJ(option%mycomm,ugdm%ndof,ndof_local,ndof_local, &
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
      call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
                        ugdm%ndof, &
                        PETSC_DETERMINE,vec,ierr)
      call VecSetLocalToGlobalMapping(vec,ugdm%mapping_ltog,ierr)
      call VecSetLocalToGlobalMappingBlock(vec,ugdm%mapping_ltogb,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
    case(LOCAL)
      call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%ngmax* &
                        ugdm%ndof, &
                        vec,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
    case(NATURAL)
      call VecCreateMPI(option%mycomm,unstructured_grid%nlmax* &
                        ugdm%ndof, &
                        PETSC_DETERMINE,vec,ierr)
      call VecSetBlockSize(vec,ugdm%ndof,ierr)
  end select
    
end subroutine UGridDMCreateVector

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

  allocate(nG2L(unstructured_grid%ngmax))
  allocate(nL2G(unstructured_grid%nlmax))
  allocate(nL2A(unstructured_grid%nlmax))
  allocate(nG2A(unstructured_grid%ngmax))
  
  ! initialize ghosted to 0
  nG2L = 0

  call ISGetIndicesF90(ugdm%is_local_petsc,int_ptr,ierr)
  do local_id = 1, unstructured_grid%nlmax
    nL2G(local_id) = local_id
    nG2L(local_id) = local_id
    ! actually, nL2A is zero-based
    !nL2A(local_id) = int_ptr(local_id)+1
    nL2A(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(ugdm%is_local_petsc,int_ptr,ierr)
!zero-based  nL2A = nL2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%nlmax, &
                            nL2A,ierr)
!zero-based  nL2A = nL2A + 1

  call ISGetIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  do ghosted_id = 1, unstructured_grid%ngmax
    nG2A(ghosted_id) = int_ptr(ghosted_id)+1
  enddo
  call ISRestoreIndicesF90(ugdm%is_ghosted_petsc,int_ptr,ierr)
  nG2A = nG2A - 1
  call AOPetscToApplication(unstructured_grid%ao_natural_to_petsc, &
                            unstructured_grid%ngmax, &
                            nG2A,ierr)
  nG2A = nG2A + 1

end subroutine UGridMapIndices

! ************************************************************************** !
!
! UGridGetCellFromPoint: Returns the cell that encompasses a point in space
! author: Glenn Hammond
! date: 10/24/09
!
! ************************************************************************** !
subroutine UGridGetCellFromPoint(x,y,z,unstructured_grid,option,icell)

  use Option_module

  implicit none
  
  PetscReal :: x, y, z
  PetscInt :: icell
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
  PetscInt :: cell_type, num_faces, iface, face_type
  PetscInt :: vertex_ids(4)
  type(plane_type) :: plane1, plane2
  type(point_type) :: point, point1, point2, point3, point4
  PetscInt :: local_id, ghosted_id
  PetscReal :: distance
  PetscBool :: inside
  
  icell = 0
  
  point%x = x
  point%y = y
  point%z = z
  
  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type(ghosted_id)
    num_faces = UCellGetNFaces(cell_type)
 
    ! vertices should be ordered counter-clockwise so that a cross product
    ! of the two vectors v1-v2 and v1-v3 points outward.
    ! if the distance from the point to the planes making up the faces is always
    ! negative using counter-clockwise ordering, the point is within the volume
    ! encompassed by the faces.
    inside = PETSC_TRUE
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface)
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)
      point1 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(1),ghosted_id))
      point2 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(2),ghosted_id))
      point3 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(3),ghosted_id))
      call UCellComputePlane(plane1,point1,point2,point3)
      distance = UCellComputeDistanceFromPlane(plane1,point)
      if (distance > 0.d0) then
        inside = PETSC_FALSE
        exit
      endif
      if (face_type == QUAD_FACE_TYPE) then
        point4 = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(4),ghosted_id))
        call UCellComputePlane(plane2,point3,point4,point1)
        distance = UCellComputeDistanceFromPlane(plane2,point)
        if (distance > 0.d0) then
          inside = PETSC_FALSE
          exit
        endif
      endif
    enddo
    
    if (inside) then
      icell = local_id
      exit
    endif

  enddo
  
end subroutine UGridGetCellFromPoint

! ************************************************************************** !
!
! UGridGetCellsInRectangle: Returns the cell that encompasses a point in space
! author: Glenn Hammond
! date: 10/24/09
!
! ************************************************************************** !
subroutine UGridGetCellsInRectangle(x_min,x_max,y_min,y_max,z_min,z_max, &
                                    unstructured_grid,option,num_cells, &
                                    cell_ids,cell_face_ids)
  use Option_module
  use Utility_module, only : reallocateIntArray
  
  implicit none
                  
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  PetscInt :: num_cells
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: cell_face_ids(:)
  
  PetscInt :: cell_type, num_faces, iface, face_type
  PetscInt :: vertex_ids(4)
  PetscInt :: num_vertices, ivertex
  PetscInt :: local_id, ghosted_id
  type(point_type) :: point
  
  PetscReal :: x_min_adj, x_max_adj, y_min_adj, y_max_adj, z_min_adj, z_max_adj
  PetscReal :: pert
  PetscBool :: in_rectangle
  
  PetscInt, pointer :: temp_cell_array(:), temp_face_array(:)
  PetscInt :: temp_array_size
  
  temp_array_size = 100
  allocate(temp_cell_array(temp_array_size))
  allocate(temp_face_array(temp_array_size))
  temp_cell_array = 0
  temp_face_array = 0
  
  ! enlarge box slightly
  pert = max(1.d-8*(x_max-x_min),1.d-8)
  x_min_adj = x_min - pert 
  x_max_adj = x_max + pert 
  pert = max(1.d-8*(y_max-y_min),1.d-8)
  y_min_adj = y_min - pert 
  y_max_adj = y_max + pert 
  pert = max(1.d-8*(z_max-z_min),1.d-8)
  z_min_adj = z_min - pert 
  z_max_adj = z_max + pert 
  
  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type(ghosted_id)
    num_faces = UCellGetNFaces(cell_type)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface)
      num_vertices = UCellGetNFaceVertices(cell_type,iface)
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)
      in_rectangle = PETSC_TRUE
      do ivertex = 1, num_vertices
        point = unstructured_grid%vertices(unstructured_grid%cell_vertices(vertex_ids(ivertex),ghosted_id))
        if (point%x < x_min_adj .or. &
            point%x > x_max_adj .or. &
            point%y < y_min_adj .or. &
            point%y > y_max_adj .or. &
            point%z < z_min_adj .or. &
            point%z > z_max_adj) then
          in_rectangle = PETSC_FALSE
          exit
        endif
      enddo
     
      if (in_rectangle) then
        num_cells = num_cells + 1
        if (num_cells > temp_array_size) then
          call reallocateIntArray(temp_cell_array,temp_array_size)
          temp_array_size = temp_array_size / 2 ! convert back for next call
          call reallocateIntArray(temp_face_array,temp_array_size)
        endif
        temp_cell_array(num_cells) = local_id
        temp_face_array(num_cells) = iface
      endif

    enddo
  enddo
  
  allocate(cell_ids(num_cells))
  allocate(cell_face_ids(num_cells))
  cell_ids = temp_cell_array(1:num_cells)
  cell_face_ids = temp_face_array(1:num_cells)
  deallocate(temp_cell_array)
  nullify(temp_cell_array)
  deallocate(temp_face_array)
  nullify(temp_face_array)
  
end subroutine UGridGetCellsInRectangle

#if 0
! ************************************************************************** !
!
! UGridReadSideSet: Reads an unstructured grid sideset
! author: Glenn Hammond
! date: 12/16/11
!
! ************************************************************************** !
subroutine UGridReadSideSet(unstructured_grid1,sideset,filename,option)

  use Input_module
  use Option_module
  use String_module
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid1
  type(unstructured_sideset_type) :: sideset
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card, word
  PetscInt :: num_faces_local_save
  PetscInt :: num_faces_local
  PetscInt :: num_to_read
  PetscInt, parameter :: max_nvert_per_face = 4
  PetscInt, allocatable :: temp_int_array(:,:)

  PetscInt :: iface, ivertex, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename)

! Format of sideset file
! type: T=triangle, Q=quadrilateral
! vertn(Q) = 4
! vertn(T) = 3
! -----------------------------------------------------------------
! num_faces  (integer)
! type vert1 vert2 ... vertn  ! for face 1 (integers)
! type vert1 vert2 ... vertn  ! for face 2
! ...
! ...
! type vert1 vert2 ... vertn  ! for face num_faces
! -----------------------------------------------------------------

  card = 'Unstructured Sideset'

  call InputReadFlotranString(input,option)
  string = 'unstructured sideset'
  call InputReadStringErrorMsg(input,option,card)  

  ! read num_faces
  call InputReadInt(input,option,sideset%nmax)
  call InputErrorMsg(input,option,'number of faces',card)

  ! divide faces across ranks
  num_faces_local = sideset%nmax/option%mycommsize 
  num_faces_local_save = num_faces_local
  remainder = sideset%nmax - num_faces_local*option%mycommsize
  if (option%myrank < remainder) num_faces_local = &
                                 num_faces_local + 1

  ! allocate array to store vertices for each faces
  allocate(sideset%face_vertices(max_nvert_per_face, &
                                 num_faces_local))
  sideset%face_vertices = -999

  ! for now, read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(max_nvert_per_face, &
                            num_faces_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = -999
      num_to_read = num_faces_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do iface = 1, num_to_read
        ! read in the vertices defining the cell face
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,card)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face type',card)
        call StringToUpper(word)
        select case(word)
          case('Q')
            num_vertices = 4
          case('T')
            num_vertices = 3
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,iface))
          call InputErrorMsg(input,option,'vertex id',card)
        enddo
      enddo
      
      ! if the faces reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_faces_local
        string = trim(adjustl(string)) // ' faces stored on p0'
        print *, trim(string)
#endif
        sideset%face_vertices(:,1:num_faces_local) = &
          temp_int_array(:,1:num_faces_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' faces sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*max_nvert_per_face
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_faces_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' faces received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_faces_local*max_nvert_per_face
    call MPI_Recv(sideset%face_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif

!  unstructured_grid%nlmax = num_faces_local
!  unstructured_grid%num_vertices_local = num_vertices_local

  call InputDestroy(input)

end subroutine UGridReadSideSet
#endif
! ************************************************************************** !
!
! UGridMapSideSet: Maps a global boundary side set to the faces of local 
!                  ghosted cells
! author: Glenn Hammond
! date: 12/16/11
!
! ************************************************************************** !
subroutine UGridMapSideSet(unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)

  use Option_module

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"

  type(unstructured_grid_type) :: unstructured_grid
  PetscInt :: face_vertices(:,:)
  PetscInt :: n_ss_faces
  character(len=MAXWORDLENGTH) :: region_name
  type(option_type) :: option
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)
  
  Mat :: Mat_vert_to_face 
  Vec :: Vertex_vec, Face_vec
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4,1)
  PetscReal :: real_array4(4)
  PetscInt, allocatable :: boundary_faces(:)
  PetscInt, allocatable :: temp_int(:,:)
  PetscInt :: boundary_face_count
  PetscInt :: mapped_face_count
  PetscInt :: nfaces, nvertices
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: local_id
  PetscInt :: cell_type
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ivertex, cell_id
  PetscErrorCode :: ierr
    
  ! fill matrix with boundary faces of local cells
  ! count up the number of boundary faces
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    nfaces = UCellGetNFaces(unstructured_grid%cell_type(local_id))
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
      endif
    enddo
  enddo
  allocate(boundary_faces(boundary_face_count))
  boundary_faces = 0
  ! assume 4 vertices per face for simplicity
  call MatCreateSeqAIJ(PETSC_COMM_SELF,boundary_face_count, &
                       unstructured_grid%num_vertices_local,4, &
                       PETSC_NULL_INTEGER,Mat_vert_to_face,ierr)
  call MatZeroEntries(Mat_vert_to_face,ierr)
  real_array4 = 1.d0
  boundary_face_count = 0
  do local_id = 1, unstructured_grid%nlmax
    cell_type = unstructured_grid%cell_type(local_id)
    nfaces = UCellGetNFaces(cell_type)
    do iface = 1, nfaces
      face_id = unstructured_grid%cell_to_face_ghosted(iface,local_id)
      if (unstructured_grid%face_to_cell_ghosted(2,face_id) < 1) then
        ! boundary face, since not connected to 2 cells
        boundary_face_count = boundary_face_count + 1
        boundary_faces(boundary_face_count) = face_id
        call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                        int_array4)
        do ivertex = 1, nvertices
          int_array4_0(ivertex,1) = &
            unstructured_grid%cell_vertices(int_array4(ivertex),local_id)-1
        enddo
        call MatSetValues(Mat_vert_to_face,1,boundary_face_count-1, &
                          nvertices,int_array4_0,real_array4, &
                          INSERT_VALUES,ierr)
      endif
    enddo
  enddo
  call MatAssemblyBegin(Mat_vert_to_face,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Mat_vert_to_face,MAT_FINAL_ASSEMBLY,ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  string = 'Mat_vert_to_face_' // trim(region_name) // '_' // &
           trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call MatView(Mat_vert_to_face,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%num_vertices_local, &
                    Vertex_vec,ierr)
  call VecZeroEntries(Vertex_vec,ierr)
  
  call VecGetArrayF90(Vertex_vec,vec_ptr,ierr)
  do iface = 1, n_ss_faces
    do ivertex = 1, size(face_vertices,1)
      if (face_vertices(ivertex,iface) > 0) then
        vec_ptr(face_vertices(ivertex,iface)) = 1.d0
      endif
    enddo
  enddo
  call VecRestoreArrayF90(Vertex_vec,vec_ptr,ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  string = 'Vertex_vec_' // trim(region_name) // '_' // &
           trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(Vertex_vec,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  
  call VecCreateSeq(PETSC_COMM_SELF,boundary_face_count,Face_vec,ierr)
  call MatMult(Mat_vert_to_face,Vertex_vec,Face_vec,ierr)

#if UGRID_DEBUG
  write(string,*) option%myrank
  string = adjustl(string)
  string = 'Face_vec_' // trim(region_name) // '_' // &
           trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr)
  call VecView(Face_vec,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
    
  allocate(temp_int(MAX_FACE_PER_CELL,boundary_face_count))
  temp_int = 0
  
  mapped_face_count = 0
  call VecGetArrayF90(Face_vec,vec_ptr,ierr)
  do iface = 1, boundary_face_count
    face_id = boundary_faces(iface)
    if (vec_ptr(iface) > 2.d0) then ! 3 or more vertices in sideset
      ! need to ensure that the right number of vertices are included
      cell_id = unstructured_grid%face_to_cell_ghosted(1,face_id)
      cell_type = unstructured_grid%cell_type(cell_id)
      nfaces = UCellGetNFaces(cell_type)
      nvertices = 0
      do iface2 = 1, nfaces
        face_id2 = unstructured_grid%cell_to_face_ghosted(iface2,cell_id)
        if (face_id == face_id2) then
          nvertices = UCellGetNFaceVertices(cell_type,iface2)
          exit
        endif
      enddo
      if (nvertices == 0) then ! the case if not found 
        option%io_buffer = 'Face not found in UGridMapSideSet'
        call printErrMsgByRank(option)
      endif
      if (abs(nvertices - vec_ptr(iface)) < 0.5d0) then
        mapped_face_count = mapped_face_count + 1
        temp_int(1,mapped_face_count) = cell_id
        temp_int(2,mapped_face_count) = iface2
      endif
    endif
  enddo
  call VecRestoreArrayF90(Face_vec,vec_ptr,ierr)
  deallocate(boundary_faces)
  
  allocate(cell_ids(mapped_face_count))
  allocate(face_ids(mapped_face_count))
  
  cell_ids(:) = temp_int(1,1:mapped_face_count)
  face_ids(:) = temp_int(2,1:mapped_face_count)
  deallocate(temp_int)
  
  call MatDestroy(Mat_vert_to_face,ierr)
  call VecDestroy(Vertex_vec,ierr)
  call VecDestroy(Face_vec,ierr)
  
end subroutine UGridMapSideSet

! ************************************************************************** !
!
! UGridDestroy: Deallocates a unstructured grid
! author: Glenn Hammond
! date: 11/01/09
!
! ************************************************************************** !
subroutine UGridDestroy(unstructured_grid)

  implicit none
  
  type(unstructured_grid_type), pointer :: unstructured_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(unstructured_grid)) return

  if (associated(unstructured_grid%cell_type)) &
    deallocate(unstructured_grid%cell_type)
  nullify(unstructured_grid%cell_type)
  if (associated(unstructured_grid%cell_vertices)) &
    deallocate(unstructured_grid%cell_vertices)
  nullify(unstructured_grid%cell_vertices)
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
  if (associated(unstructured_grid%cell_neighbors_local_ghosted)) &
    deallocate(unstructured_grid%cell_neighbors_local_ghosted)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
  if (associated(unstructured_grid%face_to_cell_ghosted))&
    deallocate(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_cell_ghosted)
  if (associated(unstructured_grid%face_to_vertex_natural))&
    deallocate(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex_natural)
  if (associated(unstructured_grid%face_to_vertex))&
    deallocate(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%face_to_vertex)
  if (associated(unstructured_grid%connection_to_face)) &
    deallocate(unstructured_grid%connection_to_face)
  nullify(unstructured_grid%connection_to_face)

  if (unstructured_grid%ao_natural_to_petsc /= 0) &
    call AODestroy(unstructured_grid%ao_natural_to_petsc,ierr)

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
  call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltog,ierr)
  if (ugdm%mapping_ltogb /= 0) &
    call ISLocalToGlobalMappingDestroy(ugdm%mapping_ltogb,ierr)
  call VecDestroy(ugdm%global_vec,ierr)
  call VecDestroy(ugdm%local_vec,ierr)
  deallocate(ugdm)
  nullify(ugdm)

end subroutine UGridDMDestroy

end module Unstructured_Grid_module
