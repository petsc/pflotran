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
    !geh begin added
    !geh end added
    PetscInt, pointer :: cell_type(:)
    PetscInt, pointer :: cell_type_ghosted(:)
    PetscInt, pointer :: cell_vertices_0(:,:) ! vertices for each grid cell (zero-based)
    PetscInt, pointer :: cell_vertices_natural(:,:) ! vertices for each grid cell (natural index 0-based)
    PetscInt, pointer :: face_to_cell_ghosted(:,:) !
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

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: WEDGE_TYPE        = 2
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_CELL = 8
  !  PetscInt, parameter :: MAX_DUALS         = 6
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
  nullify(unstructured_grid%cell_type)
  nullify(unstructured_grid%cell_type_ghosted)
  nullify(unstructured_grid%cell_vertices_0)
  nullify(unstructured_grid%cell_vertices_natural)
  nullify(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%cell_to_face_ghosted)
  nullify(unstructured_grid%vertex_ids_natural)
  nullify(unstructured_grid%vertices)
  nullify(unstructured_grid%hash)
  nullify(unstructured_grid%cell_ids_natural)
  nullify(unstructured_grid%cell_ids_petsc)
  nullify(unstructured_grid%ghost_cell_ids_natural)
  nullify(unstructured_grid%ghost_cell_ids_petsc)
  nullify(unstructured_grid%cell_neighbors_local_ghosted)
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
  num_cells_local = unstructured_grid%num_cells_global/ &
                                      option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%num_cells_global - &
              num_cells_local*option%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1

  ! allocate array to store vertices for each cell
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL,num_cells_local))
  unstructured_grid%cell_vertices_0 = -1

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
#ifdef GLENN
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'element type',card)
        call StringToUpper(word)
        select case(word)
          case('H')
            num_vertices = 8
          case('T')
            num_vertices = 4
          case('W')
            num_vertices = 6
        end select
#else
        call InputReadInt(input,option,num_vertices)
        call InputErrorMsg(input,option,'num_vertices',card)
        if (num_vertices > MAX_VERT_PER_CELL) then
          option%io_buffer = 'Cells verticies exceed maximum number of vertices'
          call printErrMsg(option)
        endif
        if ((num_vertices /= 8).and.(num_vertices /= 6)) then
          write(*,*),num_vertices
          option%io_buffer = 'Only cells with 6 or 8 vertices supported'
          call printErrMsg(option)
        endif
#endif
#endif
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
        unstructured_grid%cell_vertices_0(:,1:num_cells_local) = &
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
        int_mpi = num_to_read*MAX_VERT_PER_CELL
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
    int_mpi = num_cells_local*MAX_VERT_PER_CELL
    call MPI_Recv(unstructured_grid%cell_vertices_0,int_mpi, &
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

  unstructured_grid%num_cells_local = num_cells_local
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
  num_cells_local = unstructured_grid%num_cells_global/ &
                                      option%mycommsize 
  num_cells_local_save = num_cells_local
  remainder = unstructured_grid%num_cells_global - &
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
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL, &
                                            num_cells_local))
  unstructured_grid%cell_vertices_0 = -1
  
  do ii = 1, num_cells_local
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
  
  unstructured_grid%num_cells_local = num_cells_local
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
  unstructured_grid%num_cells_global = dataset_dims(2)
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL, &
                                            num_cells_local))
  unstructured_grid%cell_vertices_0 = -1

  ! Fill the cell data structure
  do ii = 1, num_cells_local
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
  
#ifdef ENABLE_UNSTRUCTURED  
  PetscInt :: local_id, local_id2
  PetscInt :: ghosted_id
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
  PetscInt :: global_offset_old ! geh
  PetscInt :: global_offset_new ! geh
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
  
  num_cells_local_old = unstructured_grid%num_cells_local  !sp 
  allocate(local_vertices(max_vertex_count*num_cells_local_old))
  allocate(local_vertex_offset(num_cells_local_old+1))
  count = 0
  local_vertex_offset(1) = 0
  do local_id = 1, num_cells_local_old
    do ivertex = 1, max_vertex_count
      !geh: removed the #ifdef
      if (unstructured_grid%cell_vertices_0(ivertex,local_id) < 0) exit
      count = count + 1
      local_vertices(count) = unstructured_grid%cell_vertices_0(ivertex,local_id)
    enddo
    local_vertex_offset(local_id+1) = count 
  enddo
  index_format_flag = 0 ! C-style indexing
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
  call MatMeshToCellGraph(Adj_mat,num_common_vertices,Dual_mat,ierr)
  call MatDestroy(Adj_mat,ierr)
  !TODO(geh): why do I deallocate here when it says that it is not necessary above? And
  !           why does the code not crash?
  deallocate(local_vertices)
  deallocate(local_vertex_offset)
  
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
  num_cells_local_new=cell_counts(option%myrank+1) !sp 
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

  ! Create a mapping of local indices to global strided (use block ids, not indices)
  call ISCreateBlock(option%mycomm,stride, &
                     num_cells_local_old, &   !sp 
                     index_ptr,PETSC_COPY_VALUES,is_scatter,ierr)
  call ISRestoreIndicesF90(is_num,index_ptr,ierr)
  call ISDestroy(is_num,ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_scatter_elem_old_to_new.out',viewer,ierr)
  call ISView(is_scatter,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
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
  do local_id = 1, num_cells_local_old !sp 
    count = count + 1
    ! set global cell id
    ! negate to indicate cell id with 1-based numbering (-0 = 0)
    vec_ptr(count) = -(global_offset_old+local_id)
    count = count + 1
    ! add the separator
    vec_ptr(count) = -777  ! help differentiate
    ! add the vertex ids
    do ivertex = 1, max_vertex_count
      count = count + 1
      vertex_count = vertex_count + 1
      ! increment for 1-based ordering
!      vec_ptr(count) = local_vertices(vertex_count) + 1
      vec_ptr(count) = unstructured_grid%cell_vertices_0(ivertex,local_id) + 1
    enddo


    count = count + 1 
    ! another vertex/dual separator
    vec_ptr(count) = -888  ! help differentiate

    ! add the dual ids
    istart = ia_ptr(local_id)
    iend = ia_ptr(local_id+1)-1
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
 
#if UGRID_DEBUG
  call printMsg(option,'Before element scatter')
#endif

  ! scatter all the cell data from the old decomposition (as read in in parallel)
  ! to the more parmetis-calculated decomposition
  call VecScatterCreate(elements_old,PETSC_NULL,elements_natural,is_scatter,vec_scatter,ierr)
  call ISDestroy(is_scatter,ierr)
  call VecScatterBegin(vec_scatter,elements_old,elements_natural,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scatter,elements_old,elements_natural,INSERT_VALUES,SCATTER_FORWARD,ierr)
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
  
  
  !sp 
  ! num_cells_local may be different with new partitioning 
  !  also need to update global_offset 
  global_offset_new = 0
  call MPI_Exscan(num_cells_local_new,global_offset_new, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  !sp end 

  !geh: this was not being deallocated
  deallocate(unstructured_grid%cell_vertices_0)
  allocate(unstructured_grid%cell_vertices_0(MAX_VERT_PER_CELL,num_cells_local_new))
  unstructured_grid%cell_vertices_0 = 0
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
 
  
  ! store the natural grid cell id for each local cell as read from the grid file
  call VecGetArrayF90(elements_natural,vec_ptr,ierr)
  do local_id = 1, num_cells_local_new
    unstructured_grid%cell_ids_natural(local_id) = abs(vec_ptr((local_id-1)*stride+1))
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
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit ! here we hit the -999 at the end of the stride 
      count = count + 1
    enddo
  enddo     
               
  ! allocate and fill an array with the natural cell and dual ids
  allocate(int_array(count))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    int_array(count) = unstructured_grid%cell_ids_natural(local_id)
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit ! again we hit the -999
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
  call VecGetArrayF90(elements_natural,vec_ptr2,ierr)
  allocate(unstructured_grid%cell_ids_petsc(num_cells_local_new))
  count = 0
  do local_id=1, num_cells_local_new
    count = count + 1
    ! extract the petsc id for the cell
    unstructured_grid%cell_ids_petsc(local_id) = int_array(count)
    ! store it in the elements_petsc vector too
    vec_ptr((local_id-1)*stride+1) = int_array(count)
    do idual = 1, max_dual
      dual_id = vec_ptr2(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      ! store the petsc numbered duals in the vector also
      vec_ptr(idual + dual_offset + (local_id-1)*stride) = int_array(count)
    enddo
  enddo                
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  call VecRestoreArrayF90(elements_natural,vec_ptr2,ierr)
  deallocate(int_array)
  call VecDestroy(elements_natural,ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! make a list of ghosted ids in petsc numbering
  
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ghost_cell_count = 0
  ! allocate a temporarily-sized array
  ! geh: this assumes that the number of ghost cells will not exceed the number of local
  !      and 100 is used to ensure that if this is not true, the array is still large
  !      enough
  max_ghost_cell_count = max(num_cells_local_new,100)
  allocate(unstructured_grid%ghost_cell_ids_petsc(max_ghost_cell_count))
  ! loop over all duals and find the off-processor cells on the other
  ! end of a dual
  do local_id=1, num_cells_local_new
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      found = PETSC_FALSE
      if (dual_id < 1) exit
      !TODO(geh): add back in check based on global offset to determine whether a cell
      !            is ghosted or not, must faster (see changeset b97f7c29bc38)
      do local_id2 = 1, num_cells_local_new
        if (dual_id == unstructured_grid%cell_ids_petsc(local_id2)) then
          vec_ptr(idual + dual_offset + (local_id-1)*stride) = local_id2
          found = PETSC_TRUE
          exit
        endif
      enddo
      ! if not found, add it to the list of ghost cells
      if (.not.found) then
        !sp  but only if it is not already in the list 
        found = PETSC_FALSE 
        do local_id2 = 1, ghost_cell_count  
          if (dual_id == unstructured_grid%ghost_cell_ids_petsc(local_id2)) then 
            vec_ptr(idual + dual_offset + (local_id-1)*stride) = -local_id2
            found = PETSC_TRUE
            exit
          end if 
        end do 
        if ( .not. found) then ! if here then not on processor and not in the ghost list  
          ghost_cell_count = ghost_cell_count + 1
          ! reallocate the ghost cell array if necessary
          if (ghost_cell_count > max_ghost_cell_count) then
            call reallocateIntArray(unstructured_grid%ghost_cell_ids_petsc, &
                max_ghost_cell_count)
            endif
            unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count) = dual_id
            ! flag the id as negative for ghost cell and set to current count
            vec_ptr(idual + dual_offset + (local_id-1)*stride) = -ghost_cell_count
          end if 
        endif
     enddo
  enddo
  !TODO(geh): remove the check for duplicates above and add an algorithm to remove 
  !           duplicates afterwards.  As I mention below, it will be faster.
 
  ! This is just a note to add an algorithm that adds all ghost cells and
  ! then removed duplicates.  This will perform much faster.  Do not remove
  ! this print statemetn
  call printMsg(option,'Glenn: Add code that adds all ghost cells and removed duplicates')
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_local_unsorted.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
 

  unstructured_grid%num_ghost_cells = ghost_cell_count
  unstructured_grid%num_cells_ghosted = &
    num_cells_local_new + ghost_cell_count

  ! sort ghost cell ids
  allocate(int_array(ghost_cell_count))
  allocate(int_array2(ghost_cell_count))
!geh: no don't need this  allocate(int_array3(ghost_cell_count))
  do ghosted_id = 1, ghost_cell_count
   int_array(ghosted_id) = unstructured_grid%ghost_cell_ids_petsc(ghosted_id)
!   int_array3(ghosted_id)=abs(int_array(ghosted_id) ) 
   int_array2(ghosted_id) = ghosted_id
  enddo
  !sp end 

  if (minval(int_array) < 0) then
    write(string,*) minval(int_array)
    option%io_buffer = 'Negative ghost cell value: ' // trim(adjustl(string))
    call printErrMsgByRank(option)
  endif

  ! convert to 0-based
  int_array2 = int_array2-1
!  call PetscSortIntWithPermutation(ghost_cell_count,int_array3,int_array2,ierr)
  call PetscSortIntWithPermutation(ghost_cell_count,int_array,int_array2,ierr)
  ! convert to 1-based
  int_array2 = int_array2+1


  ! resize ghost cell array down to ghost_cell_count
  deallocate(unstructured_grid%ghost_cell_ids_petsc)
  allocate(unstructured_grid%ghost_cell_ids_petsc(ghost_cell_count))
  ! fill with the sorted ids
  do ghosted_id = 1, ghost_cell_count
    unstructured_grid%ghost_cell_ids_petsc(int_array2(ghosted_id)) = &
      int_array(ghosted_id)
  enddo
    
  ! now, rearrange the ghost cell ids of the dual accordingly
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  ! add num_cells_local to get in local numbering
  ! basically, the first num_cells_local are the non-ghosted.  After that
  ! the cells ids are ghosted (at the end of the local vector)
  do local_id=1, num_cells_local_new
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 0) then
        vec_ptr(idual + dual_offset + (local_id-1)*stride) = int_array2(-dual_id)  & 
          + num_cells_local_new
      endif       
    enddo
  enddo
  deallocate(int_array)
  deallocate(int_array2)
        
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_local.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  ! load cell neighbors into array
  ! start first index at zero to store # duals for a cell
  allocate(unstructured_grid%cell_neighbors_local_ghosted(0:max_dual, &
                                                          num_cells_local_new))
  unstructured_grid%cell_neighbors_local_ghosted = 0

  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do local_id=1, num_cells_local_new
    count = 0
    do idual = 1, max_dual
      dual_id = vec_ptr(idual + dual_offset + (local_id-1)*stride)
      if (dual_id < 1) exit
      count = count + 1
      ! flag ghosted cells in dual as negative
!TODO(geh): delete the line below if -dual_id does not matter
!geh: but it does matter as ghosted neighbors in duals should have an id < 0.
!sp      if (dual_id > num_cells_local_new) dual_id = -dual_id
      unstructured_grid%cell_neighbors_local_ghosted(idual,local_id) = dual_id
    enddo
    ! set the # of duals in for the cell
    unstructured_grid%cell_neighbors_local_ghosted(0,local_id) = count
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)

  ! make a list of local vertices
  max_int_count = num_cells_local_new
  allocate(int_array_pointer(max_int_count))
  int_array_pointer = 0
  vertex_count = 0
  ! yep - load them all into a petsc vector
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do local_id=1, num_cells_local_new
    do ivertex = 1, max_vertex_count
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

  allocate(unstructured_grid%vertex_ids_natural(vertex_count))
  unstructured_grid%vertex_ids_natural = int_array3(1:vertex_count)
  !do ivertex=1,vertex_count
  !  if (option%myrank == 0) write(*,*), 'int_array3(',ivertex,')=',int_array3(ivertex)
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
                                             num_cells_local_new))
  allocate(unstructured_grid%cell_vertices_natural(1:max_vertex_count, &
                                                  num_cells_local_new))
  ! initialize to -999 (for error checking later)
  unstructured_grid%cell_vertices_0 = -999
  ! initialized the 0 entry (which stores the # of vertices in each cell) to zero
  unstructured_grid%cell_vertices_0(0,:) = 0
  
  ! permute the local ids calculated earlier in the int_array4
  call VecGetArrayF90(elements_petsc,vec_ptr,ierr)
  do local_id=1, num_cells_local_new
    do ivertex = 1, max_vertex_count
      ! extract the original vertex id
      vertex_id = vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride)
      if (vertex_id < 1) exit
      count = unstructured_grid%cell_vertices_0(0,local_id)+1
      unstructured_grid%cell_vertices_0(count,local_id) = int_array4(vertex_id)-1
      unstructured_grid%cell_vertices_0(0,local_id) = count
      ! load the permuted value back into the petsc vector
      vec_ptr(ivertex + vertex_ids_offset + (local_id-1)*stride) = int_array4(vertex_id)
!TODO(geh): no need to store natural cell vertices for each cell -- remove
      unstructured_grid%cell_vertices_natural(count,local_id) = &
          int_array3(unstructured_grid%cell_vertices_0(count,local_id) + 1) -1
    enddo
  enddo
  call VecRestoreArrayF90(elements_petsc,vec_ptr,ierr)
  deallocate(int_array2)
  deallocate(int_array3)
  deallocate(int_array4)

#if UGRID_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'elements_petsc_vert_local.out',viewer,ierr)
  call VecView(elements_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  call VecDestroy(elements_petsc,ierr)

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
    unstructured_grid%vertices(ivertex)%id = needed_vertices_petsc(ivertex) !sp  
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

  unstructured_grid%nlmax = num_cells_local_new
  unstructured_grid%ngmax = num_cells_local_new + &
       unstructured_grid%num_ghost_cells
  unstructured_grid%num_cells_ghosted = &
    num_cells_local_new + unstructured_grid%num_ghost_cells

#ifdef GLENN
  allocate(unstructured_grid%cell_type(unstructured_grid%nlmax))
  do local_id = 1, unstructured_grid%nlmax
    ! Determine number of faces and cell-type of the current cell
    select case(unstructured_grid%cell_vertices_0(0,local_id))
      case(8)
        unstructured_grid%cell_type(local_id) = HEX_TYPE
      case(6)
        unstructured_grid%cell_type(local_id) = WEDGE_TYPE
      case(4)
        unstructured_grid%cell_type(local_id) = TET_TYPE
      case default
        option%io_buffer = 'Cell type not recognized: '
        call printErrMsg(option)
    end select
  enddo
#endif
  
  unstructured_grid%num_cells_local = num_cells_local_new  
  unstructured_grid%global_offset = global_offset_new  

#endif
  
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
  
#ifdef ENABLE_UNSTRUCTURED  
  PetscInt, pointer :: int_ptr(:)
  PetscInt :: istart, iend
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
  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local*ndof, &
                    PETSC_DETERMINE,ugdm%global_vec,ierr)
  call VecSetBlockSize(ugdm%global_vec,ndof,ierr)
  ! create local vec
  call VecCreateSeq(PETSC_COMM_SELF,unstructured_grid%num_cells_ghosted*ndof, &
                    ugdm%local_vec,ierr)
  call VecSetBlockSize(ugdm%local_vec,ndof,ierr)
  
  ! IS for global numbering of local, non-ghosted cells
  call VecGetOwnershipRange(ugdm%global_vec,istart,iend,ierr)
  allocate(int_array(unstructured_grid%num_cells_local))
  do local_id = 1, unstructured_grid%num_cells_local
    int_array(local_id) = (local_id-1)+istart
  enddo

  ! arguments for ISCreateBlock():
  ! option%mycomm  - the MPI communicator
  ! ndof  - number of elements in each block
  ! unstructured_grid%num_cells_local  - the length of the index set
  !                                      (the number of blocks
  ! int_array  - the list of integers, one for each block and count
  !              of block not indices
  ! PETSC_COPY_VALUES  - see PetscCopyMode, only PETSC_COPY_VALUES and
  !                      PETSC_OWN_POINTER are supported in this routine
  ! ugdm%is_local_petsc - the new index set
  ! ierr - PETScErrorCode
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
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
    int_array(ghosted_id) = (ghosted_id+unstructured_grid%num_cells_local-1)
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
    int_array(ghosted_id) = (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
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
  allocate(int_array(unstructured_grid%num_cells_local))
  do local_id = 1, unstructured_grid%num_cells_local
    int_array(local_id) = (local_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_local,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_local_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
  
  ! IS for ghosted numbering of local ghosted cells
  allocate(int_array(unstructured_grid%num_cells_ghosted))
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    int_array(ghosted_id) = (ghosted_id-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_ghosted_local,ierr)
  deallocate(int_array)
  
#if UGRID_DEBUG
  string = 'is_ghosted_local' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_ghosted_local,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif
             
  ! IS for petsc numbering of local ghosted cells
  allocate(int_array(unstructured_grid%num_cells_ghosted))
  do local_id = 1, unstructured_grid%num_cells_local
    int_array(local_id) = istart+(local_id-1)
  enddo
  do ghosted_id = 1,unstructured_grid%num_ghost_cells
    int_array(unstructured_grid%num_cells_local+ghosted_id) = &
      (unstructured_grid%ghost_cell_ids_petsc(ghosted_id)-1)
  enddo
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_ghosted, &
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
  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local, &
                    PETSC_DETERMINE,vec_tmp,ierr)
  call VecGetOwnershipRange(vec_tmp,istart,iend,ierr)
  call VecDestroy(vec_tmp,ierr)
  allocate(int_array(unstructured_grid%num_cells_local))
  do local_id = 1, unstructured_grid%num_cells_local 
    int_array(local_id) = (local_id-1)+istart
  enddo
  call ISCreateGeneral(option%mycomm,unstructured_grid%num_cells_local, &
                       int_array,PETSC_COPY_VALUES,is_tmp,ierr) 
  deallocate(int_array)
  call AOPetscToApplicationIS(unstructured_grid%ao_natural_to_petsc, &
                              is_tmp,ierr)
  ! remap for ndof > 1  !geh: no longer need to accommodate ndof > 1, but leave
  ! alone for now.
  allocate(int_array(unstructured_grid%num_cells_local))
  call ISGetIndicesF90(is_tmp,int_ptr,ierr)
  do local_id = 1, unstructured_grid%num_cells_local
    int_array(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_tmp,int_ptr,ierr)
  call ISDestroy(is_tmp,ierr)
  call ISCreateBlock(option%mycomm,ndof,unstructured_grid%num_cells_local, &
                     int_array,PETSC_COPY_VALUES,ugdm%is_local_natural,ierr)
  deallocate(int_array)

#if UGRID_DEBUG
  string = 'is_local_natural' // trim(ndof_word) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call ISView(ugdm%is_local_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif

  call VecCreateMPI(option%mycomm,unstructured_grid%num_cells_local*ndof, &
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

  type(connection_set_type), pointer :: connections
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
  PetscInt :: count, i
  PetscInt :: iface, iface2
  PetscInt :: face_id, face_id2
  PetscInt :: ghosted_id, ghosted_id2
  PetscInt :: local_id, local_id2
  PetscInt :: cell_id, cell_id2
  PetscInt :: dual_local_id
  PetscInt :: ivertex, ivertex2
  PetscInt :: vertex_id, vertex_id2
  PetscInt :: vertex_ids4(4)
  PetscInt :: nfaces, nfaces2, nvertices, nvertices2, cell_type, cell_type2
  PetscInt :: face_type
  PetscBool:: face_found, vertex_found
  
  PetscReal :: v1(3), v2(3), v3(3), n1(3), n2(3), n_up_dn(3)
  PetscReal :: vcross(3), magnitude
  PetscReal :: area1, area2
  PetscReal :: dist_up, dist_dn
  
  type(plane_type) :: plane1, plane2
  type(point_type) :: point1, point2, point3, point4
  type(point_type) :: point_up, point_dn
  type(point_type) :: intercept1, intercept2, intercept

  character(len=MAXSTRINGLENGTH) :: string  

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
  !sp end 
  
  !sp extend cell_vertices_0 to include ghosted cells 
 
  allocate(ltog( unstructured_grid%num_vertices_local))
  do ivert=1, unstructured_grid%num_vertices_local 
    ltog(ivert) = unstructured_grid%vertices(ivert)%id 
  end do 

  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_local ,   &
         local_vec1, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_ghosted, &
         local_vec2, ierr)

  max_vertex_count = 8 
  allocate(cell_vertices_0(0:max_vertex_count, &
                           unstructured_grid%num_cells_ghosted)) 
  cell_vertices_0 = -999 

  ! first the number of vertices per cell (ivertex=0)
  ivertex = 0 
  call VecGetArrayF90(local_vec1,vec_p,ierr)
  do local_id = 1, unstructured_grid%num_cells_local
    vec_p(local_id) = unstructured_grid%cell_vertices_0(ivertex,local_id)  
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &
                        INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)

  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    cell_vertices_0(ivertex,ghosted_id)= vec_p(ghosted_id)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! now the cell ids 
  do ivertex=1, max_vertex_count  
    call VecZeroEntries(local_vec1,ierr)
    call VecGetArrayF90(local_vec1,vec_p,ierr)
    do local_id = 1, unstructured_grid%num_cells_local
     vert_id = unstructured_grid%cell_vertices_0(ivertex,local_id)  + 1 
     if (vert_id > 0) vec_p(local_id) = ltog(vert_id)
    enddo
    call VecRestoreArrayF90(local_vec1,vec_p,ierr)

    call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

    call VecGetArrayF90(local_vec2,vec_p,ierr)

    do ghosted_id = 1, unstructured_grid%num_cells_ghosted 
      do ivert=1, unstructured_grid%num_vertices_local 
        vert_id = ltog(ivert) 
        if (vert_id == vec_p(ghosted_id))  exit 
      end do  
      if (vert_id == vec_p(ghosted_id)) cell_vertices_0(ivertex,ghosted_id) = ivert  - 1  
    enddo
    call VecRestoreArrayF90(local_vec2,vec_p,ierr) 
  end do 
  deallocate(ltog)

  deallocate( unstructured_grid%cell_vertices_0) 
  allocate(unstructured_grid%cell_vertices_0(0:max_vertex_count,unstructured_grid%num_cells_ghosted) ) 
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted 
    do ivertex=0, max_vertex_count  
      unstructured_grid%cell_vertices_0(ivertex,ghosted_id)= cell_vertices_0(ivertex,ghosted_id) 
    end do 
  end do

  call VecDestroy(local_vec1,ierr) 
  call VecDestroy(local_vec2,ierr) 

  !sp end 

  ! store the cell type
  ! gb: Moved from the end of subroutine UGridDecompose(), because when
  !     running with multiple processors, information about ghosted cells
  !     is not set by the end of UGridDecompse()
  allocate(unstructured_grid%cell_type_ghosted(unstructured_grid%ngmax))
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    ! Determine number of faces and cell-type of the current cell
    select case(unstructured_grid%cell_vertices_0(0,ghosted_id))
      case(8)
        unstructured_grid%cell_type_ghosted(ghosted_id) = HEX_TYPE
      case(6)
        unstructured_grid%cell_type_ghosted(ghosted_id) = WEDGE_TYPE
      case(4)
        unstructured_grid%cell_type_ghosted(ghosted_id) = TET_TYPE
      case default
        option%io_buffer = 'Cell type not recognized: '
        call printErrMsg(option)
    end select
  enddo


  ! create mappings of [cells,faces,vertices] to [cells,faces,vertices]
  allocate(face_to_vertex(MAX_VERT_PER_FACE,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  face_to_vertex = -999
  allocate(cell_to_face(MAX_DUALS,unstructured_grid%num_cells_ghosted))
  cell_to_face = -999
  allocate(face_to_cell(2,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  face_to_cell = -999
  allocate(vertex_to_cell(0:MAX_CELLS_SHARING_A_VERTEX,unstructured_grid%num_vertices_local))
  vertex_to_cell = 0

  allocate(unstructured_grid%face_to_vertex_natural(MAX_VERT_PER_FACE, &
           MAX_DUALS*unstructured_grid%num_cells_ghosted))
  unstructured_grid%face_to_vertex_natural = -999
  allocate(unstructured_grid%face_to_cell_ghosted(1,MAX_DUALS*unstructured_grid%num_cells_ghosted))
  unstructured_grid%face_to_cell_ghosted = -999
  allocate(unstructured_grid%cell_to_face_ghosted(MAX_DUALS,&
       unstructured_grid%num_cells_ghosted))
  unstructured_grid%cell_to_face_ghosted = -999

  face_count = 0
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    cell_type = unstructured_grid%cell_type_ghosted(ghosted_id)
    ! Determine number of faces and cell-type of the current cell
    select case(cell_type)
      case(HEX_TYPE)
        nfaces = 6
      case(WEDGE_TYPE)
        nfaces = 5
      case(TET_TYPE)
        nfaces = 4
      case default
        option%io_buffer = 'Cell type not recognized'
        call printErrMsg(option)
    end select
    do iface = 1, nfaces
      face_count = face_count + 1
      cell_to_face(iface,ghosted_id) = face_count
      face_to_cell(1,face_count) = ghosted_id
      unstructured_grid%face_to_cell_ghosted(1,face_count) = ghosted_id
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids4)
      ! Determine number of vertices forming the face
      select case(cell_type)
        case(HEX_TYPE) 
          nvertices = 4
        case(WEDGE_TYPE)
          if (iface > 3) then
            nvertices = 3
          else
            nvertices = 4
          endif
        case(TET_TYPE)
          nvertices = 3
      end select
      do ivertex = 1, nvertices
        face_to_vertex(ivertex,face_count) = &
          unstructured_grid%cell_vertices_0(vertex_ids4(ivertex),ghosted_id)+1
          if (face_to_vertex(ivertex,face_count) > 0) then
            unstructured_grid%face_to_vertex_natural(ivertex,face_count) = &
              unstructured_grid%vertex_ids_natural(face_to_vertex(ivertex,face_count))
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
  do local_id = 1, unstructured_grid%num_cells_local
    ! Selet a cell and find number of vertices
    cell_id = local_id
    cell_type = unstructured_grid%cell_type_ghosted(local_id)
    select case(cell_type)
      case(HEX_TYPE)
        nfaces = 6
      case(WEDGE_TYPE)
        nfaces = 5
      case(TET_TYPE)
        nfaces = 4
    end select
    do local_id2 = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      ! Selet a neighboring cell
      cell_id2 =  unstructured_grid%cell_neighbors_local_ghosted(local_id2,local_id)
      cell_type2 = unstructured_grid%cell_type_ghosted(local_id2)
      ! If cell-id is neighbor is lower, skip it
      if (cell_id2 <= cell_id) cycle
      ! Find the number of vertices for neighboring cell
      select case(cell_type2)
        case(HEX_TYPE)
          nfaces2 = 6
        case(WEDGE_TYPE)
          nfaces2 = 5
        case(TET_TYPE)
          nfaces2 = 4
      end select      
      ! Initialize
      face_found = PETSC_FALSE
      do iface = 1, nfaces
        ! Select a face and find number of vertices forming the face
        face_id = cell_to_face(iface,cell_id)
        select case(cell_type)
          case(HEX_TYPE) 
            nvertices = 4
          case(WEDGE_TYPE)
            if (iface > 3) then
              nvertices = 3
            else
              nvertices = 4
            endif
          case(TET_TYPE)
            nvertices = 3
        end select
        do ivertex = 1, nvertices
          ! Select a vertex and initialize vertex_found
          vertex_id = face_to_vertex(ivertex,face_id) ! face_to_vertex is 1-based indexing
          vertex_found = PETSC_FALSE
          do ivertex2 = 1, unstructured_grid%cell_vertices_0(0,cell_id2)
            vertex_id2 = unstructured_grid%cell_vertices_0(ivertex2,cell_id2)+1 ! cell_vertices_0 is 0-based indexing
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
            select case(cell_type2)
              case(HEX_TYPE) 
                nvertices2 = 4
              case(WEDGE_TYPE)
                if (iface2 > 3) then
                  nvertices2 = 3
                else
                  nvertices2 = 4
                endif
              case(TET_TYPE)
                nvertices2 = 3
            end select            
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
                if (face_id2 > face_id) then
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id2, ' -> ', face_id
                  option%io_buffer = 'Duplicated face removed:' // string
                  call printMsg(option)
#endif
                  cell_to_face(iface2,cell_id2) = face_id
                  !! flag as removed
                  face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
                  face_to_cell(2,face_id) = cell_id2
                else
#ifdef UGRID_DEBUG                
                  write(string,*) option%myrank, face_id, ' -> ', face_id2
                  option%io_buffer = 'Duplicated face removed:' // string
                  call printMsg(option)
#endif
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
      if (.not.face_found) then
        write(string,*),'rank=',option%myrank, 'local_id_id',cell_id,'local_id_id2',cell_id2
        option%io_buffer='No shared face found: ' // string // '\n'
        write(*,*),'No shared face found: rank=',option%myrank, 'local_id_id',cell_id,'local_id_id2',cell_id2
        call printErrMsg(option)
      endif
    enddo ! local_id2-loop
  enddo  ! local_id-loop
#endif ! #ifdef MIXED_UMESH

   

#ifndef MIXED_UMESH
  ! remove duplicate faces
  ! fill face ids
  do local_id = 1, unstructured_grid%num_cells_local !sp   was num_cells_ghosted 
    cell_id = local_id
    do local_id2 = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      cell_id2 =  unstructured_grid%cell_neighbors_local_ghosted(local_id2,local_id)
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
#ifdef UGRID_DEBUG
              write(string,*) option%myrank, face_id2, ' -> ', face_id
              option%io_buffer = 'Duplicated face removed:' // string
              call printMsg(option)
#endif
              cell_to_face(iface2,cell_id2) = face_id
              ! flag as removed
              face_to_cell(1,face_id2) = -face_to_cell(1,face_id2)
              face_to_cell(2,face_id) = cell_id2
            else
#ifdef UGRID_DEBUG
              write(string,*) option%myrank, face_id, ' -> ', face_id2
              option%io_buffer = 'Duplicated face removed:' // string
              call printMsg(option)
#endif
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
  allocate(unstructured_grid%face_to_vertex(MAX_VERT_PER_FACE,face_count))
  count = 0
  do iface = 1, size(face_to_cell,2)
    if (face_to_cell(1,iface) > 0) then
      count = count + 1
      unstructured_grid%face_to_vertex(:,count) = face_to_vertex(:,iface)
    endif
  enddo
  deallocate(face_to_vertex)
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
    do i = 1,2
      cell_id = face_to_cell(i,face_id)
      ! check for exterior face
      if (cell_id < 0) cycle
      found = PETSC_FALSE
      do iface2 = 1, 6
        face_id2 = cell_to_face(iface2,cell_id)
        if (face_id<0) cycle
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
  
  
  do ghosted_id = 1, unstructured_grid%num_cells_ghosted
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,ghosted_id)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,ghosted_id)+1
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
  do local_id = 1, unstructured_grid%num_cells_local
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_id = unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      ! count all ghosted connections (dual_id < 0)
      ! only count connection with cells of larger ids to avoid double counts
!      if (dual_id < 0 .or. local_id < dual_id) then
      if (dual_id > 0 .and. local_id < dual_id) then !sp 
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
  do local_id = 1, unstructured_grid%num_cells_local
    do idual = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      dual_local_id = unstructured_grid%cell_neighbors_local_ghosted(idual,local_id)
      if (local_id < abs(dual_local_id)) then 
        iconn = iconn + 1
        ! find face
        found = PETSC_FALSE
        do iface = 1, unstructured_grid%cell_vertices_0(0,local_id)
          face_id = cell_to_face(iface,local_id)
          do local_id2 = 1,2
            cell_id2 = face_to_cell(local_id2,face_id)
            if (cell_id2 == abs(dual_local_id)) then
              found = PETSC_TRUE
              exit
            endif
          enddo
          if (found) exit
        enddo
        if (found) then
          dual_to_face(iconn) = face_id
        else
          write(string,*) option%myrank,local_id,dual_local_id 
          option%io_buffer = 'face not found in connection loop' // string 
          call printErrMsg(option)
        endif
        if (unstructured_grid%face_to_vertex(4,face_id) < 0) then
          face_type = TRI_FACE_TYPE
        else
          face_type = QUAD_FACE_TYPE
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
  unstructured_grid%cell_to_face_ghosted = cell_to_face
  allocate(unstructured_grid%face_centroid(face_count))
  do iface = 1,face_count
    unstructured_grid%face_centroid(iface)%id = -999
  enddo
  
  do local_id = 1, unstructured_grid%num_cells_local
    do iface = 1,MAX_DUALS
      face_id = cell_to_face(iface, local_id)
      if (face_id == -999) cycle
      if ( unstructured_grid%face_centroid(face_id)%id == -999) then
        count = 0
        unstructured_grid%face_centroid(face_id)%x = 0.d0
        unstructured_grid%face_centroid(face_id)%y = 0.d0
        unstructured_grid%face_centroid(face_id)%z = 0.d0
        if (unstructured_grid%face_to_vertex(4,face_id) < 0) then
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
          if (vertex_id.ne.-999) then
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
  
  deallocate(cell_to_face)
  deallocate(face_to_cell)
  deallocate(vertex_to_cell)
  deallocate(dual_to_face)

#endif        
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
                                   iconn, ghosted_id)

  use Connection_module
  use Utility_module, only : DotProduct
  
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface_cell
  PetscInt :: iconn
  PetscInt :: ghosted_id
  
  PetscErrorCode :: ierr
  
  PetscInt  :: face_id
  PetscInt  :: ivert,vert_id
  PetscReal :: v1(3),v2(3),n_dist(3), dist
  type(point_type) :: vertex_8(8)
  type(plane_type) :: plane
  type(point_type) :: point, vertex1, vertex2, vertex3, intercept
  
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
      ! Compute cell centeroid
      v2 = 0.d0
      do ivert = 1, unstructured_grid%cell_vertices_0(0, ghosted_id)
        vert_id = unstructured_grid%cell_vertices_0(ivert, ghosted_id) + 1
        vertex_8(ivert)%x = unstructured_grid%vertices(vert_id)%x
        vertex_8(ivert)%y = unstructured_grid%vertices(vert_id)%y
        vertex_8(ivert)%z = unstructured_grid%vertices(vert_id)%z
      enddo
      v2 = UCellComputeCentroid(unstructured_grid%cell_type_ghosted(ghosted_id), &
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
      ! TODO(geh): add support for a quad face
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
                             nL2G, &
                             grid_x,grid_y,grid_z, &
                             x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  VecScatter :: scatter_ltol 
  PetscInt :: nL2G(:)
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: local_id
  PetscInt :: ivertex
  PetscInt :: vertex_id
  type(point_type) :: vertex_8(8)
  PetscReal :: centroid(3)
  PetscReal, pointer :: vec_p(:) !sp 
  Vec :: local_vec1 !sp 
  Vec :: local_vec2 !sp 
  PetscErrorCode :: ierr 

  do local_id = 1, unstructured_grid%num_cells_local 
    do ivertex = 1, unstructured_grid%cell_vertices_0(0,local_id)
      vertex_id = unstructured_grid%cell_vertices_0(ivertex,local_id) + 1
      vertex_8(ivertex)%x = &
        unstructured_grid%vertices(vertex_id)%x
      vertex_8(ivertex)%y = &
        unstructured_grid%vertices(vertex_id)%y
      vertex_8(ivertex)%z = &
        unstructured_grid%vertices(vertex_id)%z
    enddo
#ifdef GLENN
    ! TODO(geh): check if nL2G is working correctly
    centroid = UCellComputeCentroid(unstructured_grid%cell_type(local_id),vertex_8)
#else
    select case (unstructured_grid%cell_vertices_0(0,local_id))
      case(8)
        centroid = UCellComputeCentroid(HEX_TYPE,vertex_8)
      case(6)
        centroid = UCellComputeCentroid(WEDGE_TYPE,vertex_8)
    end select
#endif    
    grid_x(local_id) = centroid(1)
    grid_y(local_id) = centroid(2)
    grid_z(local_id) = centroid(3)
  enddo

  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_local,   & 
         local_vec1, ierr) 
  call VecCreateSeq(PETSC_COMM_SELF, unstructured_grid%num_cells_ghosted, & 
         local_vec2, ierr) 

  ! x coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do local_id = 1, unstructured_grid%num_cells_local
   vec_p(local_id) = grid_x(local_id)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do local_id = 1, unstructured_grid%num_cells_ghosted
   grid_x(local_id) = vec_p(local_id)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! y coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do local_id = 1, unstructured_grid%num_cells_local
   vec_p(local_id) = grid_y(local_id)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do local_id = 1, unstructured_grid%num_cells_ghosted
   grid_y(local_id) = vec_p(local_id)
  enddo
  call VecRestoreArrayF90(local_vec2,vec_p,ierr)

  ! z coordinate
  call VecGetArrayF90(local_vec1,vec_p,ierr)  
  do local_id = 1, unstructured_grid%num_cells_local
   vec_p(local_id) = grid_z(local_id)
  enddo
  call VecRestoreArrayF90(local_vec1,vec_p,ierr)

  call VecScatterBegin(scatter_ltol,local_vec1,local_vec2,  & 
                         INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scatter_ltol,local_vec1,local_vec2,    &  
                         INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(local_vec2,vec_p,ierr)
  do local_id = 1, unstructured_grid%num_cells_ghosted
   grid_z(local_id) = vec_p(local_id)
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
#ifdef GLENN
    volume_p(local_id) = UCellComputeVolume(unstructured_grid%cell_type_ghosted(nL2G(local_id)),vertex_8)
#else
      select case (unstructured_grid%cell_vertices_0(0,ghosted_id))
        case(8)
          volume_p(local_id) = UCellComputeVolume(HEX_TYPE,vertex_8)
        case(6)
          volume_p(local_id) = UCellComputeVolume(WEDGE_TYPE,vertex_8)
      end select
#endif      
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
subroutine UGridEnsureRightHandRule(unstructured_grid,x,y,z,option)

  use Option_module
  
  implicit none

  type(unstructured_grid_type) :: unstructured_grid
  PetscReal :: x(:), y(:), z(:) 
  type(option_type) :: option

  PetscInt :: local_id
  PetscInt :: ghosted_id
  type(point_type) :: point, point1, point2, point3
  type(plane_type) :: plane1
  PetscReal :: distance
  PetscInt :: cell_vertex_ids_before(8), cell_vertex_ids_after(8)
  PetscInt :: face_vertex_ids(4)
  PetscInt :: num_vertices, iface, cell_type, num_faces, face_type

  do local_id = 1, unstructured_grid%nlmax
    ghosted_id = local_id
    cell_type = unstructured_grid%cell_type(local_id)
    num_vertices = UCellGetNVertices(cell_type)
    cell_vertex_ids_before(1:num_vertices) = &
      unstructured_grid%cell_vertices_0(1:num_vertices,ghosted_id)+1
    cell_vertex_ids_after = cell_vertex_ids_before
    point%x = x(ghosted_id)
    point%y = y(ghosted_id)
    point%z = z(ghosted_id)
    num_faces = UCellGetNFaces(cell_type)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface)
      call UCellGetFaceVertices(option,cell_type,iface,face_vertex_ids)
      point1 = unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(1)))
      point2 = unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(2)))
      point3 = unstructured_grid%vertices(cell_vertex_ids_before(face_vertex_ids(3)))
      call UCellComputePlane(plane1,point1,point2,point3)
      distance = UCellComputeDistanceFromPlane(plane1,point) 
      if (distance > 0.d0) then
        ! need to swap so that distance is negative (point lies below plane)
        write(option%io_buffer,'(''Cell '',i3,'' with vertices: '',3i3, &
                               & '' violates right hand rule for face vertices.'')') &
                               ghosted_id, face_vertex_ids(1), &
                               face_vertex_ids(2), face_vertex_ids(3)
        call printErrMsg(option)
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
  
  allocate(d_nnz(unstructured_grid%num_cells_local))
  allocate(o_nnz(unstructured_grid%num_cells_local))
  d_nnz = 1 ! start 1 since diagonal connection to self
  o_nnz = 0
  do local_id = 1, unstructured_grid%num_cells_local
    do ineighbor = 1, unstructured_grid%cell_neighbors_local_ghosted(0,local_id)
      neighbor_id = unstructured_grid%cell_neighbors_local_ghosted(ineighbor,local_id)
      if (neighbor_id > 0) then
        d_nnz(local_id) = d_nnz(local_id) + 1
      else
        o_nnz(local_id) = o_nnz(local_id) + 1
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
        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltog,ugdm%mapping_ltog,ierr)
!        call MatSetLocalToGlobalMapping(J,ugdm%mapping_ltogb,ugdm%mapping_ltogb,ierr)
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
  
  icell = -999
  
  point%x = x
  point%y = y
  point%z = z
  
  do local_id = 1, unstructured_grid%num_cells_local
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type_ghosted(ghosted_id)
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
      point1 = unstructured_grid%vertices(unstructured_grid%cell_vertices_0(vertex_ids(1),ghosted_id)+1)
      point2 = unstructured_grid%vertices(unstructured_grid%cell_vertices_0(vertex_ids(2),ghosted_id)+1)
      point3 = unstructured_grid%vertices(unstructured_grid%cell_vertices_0(vertex_ids(3),ghosted_id)+1)
      call UCellComputePlane(plane1,point1,point2,point3)
      distance = UCellComputeDistanceFromPlane(plane1,point)
      if (distance > 0.d0) then
        inside = PETSC_FALSE
        exit
      endif
      if (face_type == QUAD_FACE_TYPE) then
        point4 = unstructured_grid%vertices(unstructured_grid%cell_vertices_0(vertex_ids(4),ghosted_id)+1)
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
  
  do local_id = 1, unstructured_grid%num_cells_local
    ghosted_id = local_id ! ghosted ids are same for first nlocal cells
    cell_type = unstructured_grid%cell_type_ghosted(ghosted_id)
    num_faces = UCellGetNFaces(cell_type)
    do iface = 1, num_faces
      face_type = UCellGetFaceType(cell_type,iface)
      num_vertices = UCellGetNFaceVertices(cell_type,iface)
      call UCellGetFaceVertices(option,cell_type,iface,vertex_ids)
      in_rectangle = PETSC_TRUE
      do ivertex = 1, num_vertices
        point = unstructured_grid%vertices(unstructured_grid%cell_vertices_0(vertex_ids(ivertex),ghosted_id)+1)
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
  if (associated(unstructured_grid%cell_type_ghosted)) &
    deallocate(unstructured_grid%cell_type_ghosted)
  nullify(unstructured_grid%cell_type_ghosted)
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
  if (associated(unstructured_grid%cell_vertices_natural))&
    deallocate(unstructured_grid%cell_vertices_natural)
  nullify(unstructured_grid%cell_vertices_natural)
  if (associated(unstructured_grid%face_to_cell_ghosted))&
    deallocate(unstructured_grid%face_to_cell_ghosted)
  nullify(unstructured_grid%face_to_cell_ghosted)
  if (associated(unstructured_grid%face_to_vertex_natural))&
    deallocate(unstructured_grid%face_to_vertex_natural)
  nullify(unstructured_grid%face_to_vertex_natural)
  if (associated(unstructured_grid%face_to_vertex))&
    deallocate(unstructured_grid%face_to_vertex)
  nullify(unstructured_grid%face_to_vertex)

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
