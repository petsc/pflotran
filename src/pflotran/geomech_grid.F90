#ifdef GEOMECH
module Geomechanics_Grid_module

  use Geomechanics_Grid_Aux_module
  use Unstructured_Cell_module
  use PFLOTRAN_Constants_module

  implicit none

  private 
  
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#if defined(SCORPIO)
  include "scorpiof.h"
#endif

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: GeomechGridRead, &
            CopySubsurfaceGridtoGeomechGrid, &
            GeomechGridLocalizeRegions, &
            GeomechGridVecGetArrayF90, &
            GeomechGridVecRestoreArrayF90, &
            GeomechGridCopyIntegerArrayToVec, &
            GeomechGridCopyVecToIntegerArray 

contains

! ************************************************************************** !
!
! GeomechGridRead: Reads a geomechanics grid
! author: Satish Karra, LANL
! date: 05/22/13
!
! ************************************************************************** !
subroutine GeomechGridRead(geomech_grid,filename,option)

  use Input_module
  use Option_module
  use String_module
  
  implicit none
  
  type(geomech_Grid_type) :: geomech_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: num_elems_local_save
  PetscInt :: num_elems_local
  PetscInt :: num_nodes_local_save
  PetscInt :: num_nodes_local
  PetscInt :: num_to_read
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscReal, allocatable :: node_coordinates(:,:)

  PetscInt :: ielem, inode, idir, irank, num_nodes
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename,option)

  ! initial guess is 8 nodes per elem
  geomech_grid%max_nnode_per_elem = 8

! Format of geomech grid file
! type: H=hexahedron, T=tetrahedron, W=wedge, P=pyramid
! nodeN(H) = 8
! nodeN(T) = 4
! nodeN(W) = 6
! nodeN(P) = 5
! -----------------------------------------------------------------
! num_elems num_nodes  (integers)
! type node1 node2 node3 ... nodeN  ! for elem 1 (integers)
! type node1 node2 node3 ... nodeN  ! for elem 2
! ...
! ...
! type node1 node2 node3 ... nodeN  ! for elem num_elems
! xcoord ycoord zcoord ! coordinates of node 1 (real)
! xcoord ycoord zcoord ! coordinates of node 2 (real)
! ...
! xcoord ycoord zcoord ! coordinates of node num_nodes (real)
! -----------------------------------------------------------------

  hint = 'Geomechanics Grid'

  call InputReadFlotranString(input,option)
  string = 'geomechanics grid'
  call InputReadStringErrorMsg(input,option,hint)  

  ! read num_elems
  call InputReadInt(input,option,geomech_grid%nmax_elem)
  call InputErrorMsg(input,option,'number of elems',hint)
  ! read num_nodes
  call InputReadInt(input,option,geomech_grid%nmax_node)
  call InputErrorMsg(input,option,'number of nodes',hint)

  ! divide elems across ranks
  num_elems_local = geomech_grid%nmax_elem/option%mycommsize 
  num_elems_local_save = num_elems_local
  remainder = geomech_grid%nmax_elem - &
              num_elems_local*option%mycommsize
  if (option%myrank < remainder) num_elems_local = &
                                 num_elems_local + 1

  ! allocate array to store nodes for each elem
  allocate(geomech_grid%elem_nodes(geomech_grid%max_nnode_per_elem, &
                                             num_elems_local))
  geomech_grid%elem_nodes = -999

  ! for now, read all elems from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(geomech_grid%max_nnode_per_elem, &
                            num_elems_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = -999
      num_to_read = num_elems_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do ielem = 1, num_to_read
        ! read in the nodes defining the grid elem
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'element type',hint)
        call StringToUpper(word)
        select case(word)
          case('H')
            num_nodes = 8
          case('W')
            num_nodes = 6
          case('P')
            num_nodes = 5
          case('T')
            num_nodes = 4
          case('Q')
            num_nodes = 4
        end select
        do inode = 1, num_nodes
          call InputReadInt(input,option,temp_int_array(inode,ielem))
          call InputErrorMsg(input,option,'node id',hint)
        enddo
      enddo
      
      ! if the elems reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_elems_local
        string = trim(adjustl(string)) // ' elems stored on p0'
        print *, trim(string)
#endif
        geomech_grid%elem_nodes(:,1:num_elems_local) = &
          temp_int_array(:,1:num_elems_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' elems sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*geomech_grid%max_nnode_per_elem
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_elems_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' elems received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    int_mpi = num_elems_local*geomech_grid%max_nnode_per_elem
    call MPI_Recv(geomech_grid%elem_nodes,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif


  ! divide nodes across ranks
  num_nodes_local = geomech_grid%nmax_node/ &
                                         option%mycommsize 
  num_nodes_local_save = num_nodes_local
  remainder = geomech_grid%nmax_node - &
              num_nodes_local*option%mycommsize
  if (option%myrank < remainder) num_nodes_local = &
                                 num_nodes_local + 1

  allocate(node_coordinates(3,num_nodes_local))
  node_coordinates = 0.d0

  ! just like above, but this time for node coordinates
  if (option%myrank == option%io_rank) then
    allocate(temp_real_array(3,num_nodes_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      num_to_read = num_nodes_local_save
      if (irank < remainder) num_to_read = num_to_read + 1
      do inode = 1, num_to_read
        call InputReadFlotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        do idir = 1, 3
          call InputReadDouble(input,option,temp_real_array(idir,inode))
          call InputErrorMsg(input,option,'node coordinate',hint)
        enddo
      enddo
      
      if (irank == option%io_rank) then
        node_coordinates(:,1:num_nodes_local) = &
          temp_real_array(:,1:num_nodes_local)
      else
        int_mpi = num_to_read*3
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_real_array)
  else
    int_mpi = num_nodes_local*3
    call MPI_Recv(node_coordinates, &
                  int_mpi, &
                  MPI_DOUBLE_PRECISION,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif
  
  ! fill the nodes data structure
  allocate(geomech_grid%nodes(num_nodes_local))
  do inode = 1, num_nodes_local
    geomech_grid%nodes(inode)%id = 0
    geomech_grid%nodes(inode)%x = node_coordinates(1,inode)
    geomech_grid%nodes(inode)%y = node_coordinates(2,inode)
    geomech_grid%nodes(inode)%z = node_coordinates(3,inode)
  enddo
  deallocate(node_coordinates)

  geomech_grid%nlmax_elem = num_elems_local
  geomech_grid%nlmax_node = num_nodes_local

  call InputDestroy(input)

end subroutine GeomechGridRead

! ************************************************************************** !
!
! CopySubsurfaceGridtoGeomechGrid: Subroutine to copy subsurface grid info.
! to geomechanics grid
! author: Satish Karra, LANL
! date: 05/30/13
!
! ************************************************************************** !
subroutine CopySubsurfaceGridtoGeomechGrid(ugrid,geomech_grid,option)
                                        
  use Unstructured_Grid_Aux_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Gauss_module
  
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
  
  type(unstructured_grid_type), pointer      :: ugrid
  type(geomech_grid_type), pointer           :: geomech_grid
  type(option_type), pointer                 :: option
  PetscInt                                   :: local_id
  PetscInt                                   :: ghosted_id
  PetscInt                                   :: vertex_count
  PetscInt                                   :: ivertex
  PetscInt                                   :: vertex_id
  PetscInt                                   :: count
  PetscInt, allocatable                      :: int_array(:)
  PetscInt, allocatable                      :: int_array2(:)
  PetscInt, allocatable                      :: int_array3(:)
  PetscInt, allocatable                      :: int_array4(:)
  PetscErrorCode                             :: ierr
  character(len=MAXSTRINGLENGTH)             :: string, string1
  PetscInt                                   :: global_offset_old
  PetscInt                                   :: global_offset
  Mat                                        :: Rank_Mat
  PetscReal                                  :: rank
  PetscViewer                                :: viewer
  PetscReal, pointer                         :: vec_ptr(:) 
  PetscInt                                   :: istart,iend
  PetscBool                                  :: vertex_found
  PetscInt                                   :: int_rank
  PetscInt                                   :: vertex_count2
  IS                                         :: is_rank 
  IS                                         :: is_rank_new
  IS                                         :: is_natural
  IS                                         :: is_ghost_petsc
  PetscReal                                  :: max_val
  PetscInt                                   :: row
  PetscScalar, allocatable                   :: val(:)
  PetscInt                                   :: ncols
  PetscInt, allocatable                      :: cols(:) 
  AO                                         :: ao_natural_to_petsc_nodes
  PetscInt                                   :: nlmax_node
  PetscInt, pointer                          :: int_ptr(:)
  type(point_type), pointer                  :: vertices(:)

  call printMsg(option,'GEOMECHANICS: Subsurface unstructured grid will ' // &
                  'be used for geomechanics.')
  
#ifdef GEOMECH_DEBUG
  call printMsg(option,'Copying unstructured grid to geomechanics grid')
#endif
  
  geomech_grid%global_offset_elem = ugrid%global_offset
  geomech_grid%nmax_elem = ugrid%nmax
  geomech_grid%nlmax_elem = ugrid%nlmax
  geomech_grid%nmax_node = ugrid%num_vertices_global
  
  allocate(geomech_grid%elem_ids_natural(geomech_grid%nlmax_elem))
  do local_id = 1, geomech_grid%nlmax_elem
    geomech_grid%elem_ids_natural(local_id) = ugrid%cell_ids_natural(local_id)
  enddo
  
  allocate(geomech_grid%elem_ids_petsc(size(ugrid%cell_ids_petsc)))
  geomech_grid%elem_ids_petsc = ugrid%cell_ids_petsc
  geomech_grid%ao_natural_to_petsc = ugrid%ao_natural_to_petsc
  geomech_grid%max_ndual_per_elem = ugrid%max_ndual_per_cell
  geomech_grid%max_nnode_per_elem = ugrid%max_nvert_per_cell
  geomech_grid%max_elem_sharing_a_node = ugrid%max_cells_sharing_a_vertex
  geomech_grid%nlmax_node = ugrid%num_vertices_natural

#ifdef GEOMECH_DEBUG
  call printMsg(option,'Removing ghosted elements (cells)')
#endif
  
  allocate(geomech_grid%elem_type(geomech_grid%nlmax_elem))
  do local_id = 1, geomech_grid%nlmax_elem
    geomech_grid%elem_type(local_id) = ugrid%cell_type(local_id)
  enddo

#ifdef GEOMECH_DEBUG
  call printMsg(option,'Reordering nodes (vertices)')
#endif

  ! Read all the cells including ghosted ones and remove the ghosted ones
  ! Find the nodes which are on the local cells
  vertex_count = 0
  ! First calculate number of nodes on local domain (without ghosted elements)
  do local_id = 1, geomech_grid%nlmax_elem
    vertex_count = vertex_count + ugrid%cell_vertices(0,local_id)
  enddo
  
  count = 0
  allocate(int_array(vertex_count))
  do local_id = 1, geomech_grid%nlmax_elem
    do ivertex = 1, ugrid%cell_vertices(0,local_id)
      count = count + 1
      int_array(count) = ugrid%cell_vertices(ivertex,local_id)
    enddo
  enddo

  ! Sort the vertex ids
  allocate(int_array2(vertex_count))
  do ivertex = 1, vertex_count
    int_array2(ivertex) = ivertex 
  enddo
  int_array2 = int_array2 - 1
  call PetscSortIntWithPermutation(vertex_count,int_array,int_array2,ierr)
  int_array2 = int_array2+1
    
  ! Remove duplicates
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
  
  allocate(geomech_grid%node_ids_ghosted_natural(vertex_count))
  do ivertex = 1, vertex_count
    geomech_grid%node_ids_ghosted_natural(ivertex) = ugrid% &
      vertex_ids_natural(int_array3(ivertex))  
  enddo

#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ids_ghosted_natural' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id = 1, vertex_count
    write(86,'(i5)') geomech_grid%node_ids_ghosted_natural(local_id)
  enddo  
  close(86)
#endif 

  allocate(geomech_grid%elem_nodes( &
            0:geomech_grid%max_nnode_per_elem,geomech_grid%nlmax_elem))
  geomech_grid%elem_nodes = 0
  
  count = 0
  do local_id = 1, geomech_grid%nlmax_elem
    geomech_grid%elem_nodes(0,local_id) = ugrid%cell_vertices(0,local_id)
    do ivertex = 1, geomech_grid%elem_nodes(0,local_id)
      count = count + 1     
      geomech_grid%elem_nodes(ivertex,local_id) = int_array4(count)
    enddo
  enddo
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_elem_nodes' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id = 1, geomech_grid%nlmax_elem
    write(86,'(i5)') geomech_grid%elem_nodes(0,local_id)
    do ivertex = 1, geomech_grid%max_nnode_per_elem
      write(86,'(i5)') geomech_grid%elem_nodes(ivertex,local_id)
    enddo
  enddo  
  close(86)
#endif   
  
  deallocate(int_array2)
  deallocate(int_array4)
    
  allocate(geomech_grid%nodes(vertex_count))
  do ivertex = 1, vertex_count
    geomech_grid%nodes(ivertex)%id = geomech_grid%node_ids_ghosted_natural(ivertex)
    geomech_grid%nodes(ivertex)%x = ugrid%vertices(int_array3(ivertex))%x
    geomech_grid%nodes(ivertex)%y = ugrid%vertices(int_array3(ivertex))%y
    geomech_grid%nodes(ivertex)%z = ugrid%vertices(int_array3(ivertex))%z    
  enddo
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_coordinates' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ivertex = 1, vertex_count
    write(86,'(i5)') geomech_grid%nodes(ivertex)%id
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%x
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%y
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%z
  enddo  
  close(86)
#endif  

  geomech_grid%ngmax_node = vertex_count

  deallocate(int_array3)
    
  call MatCreateAIJ(option%mycomm,PETSC_DECIDE,ONE_INTEGER, &
                    geomech_grid%nmax_node,option%mycommsize, &
                    option%mycommsize,PETSC_NULL_INTEGER, &
                    option%mycommsize,PETSC_NULL_INTEGER,Rank_Mat,ierr)
  
  call MatZeroEntries(Rank_Mat,ierr)
  
  rank = option%myrank + 1
  do ivertex = 1, geomech_grid%ngmax_node
    call MatSetValue(Rank_Mat,geomech_grid%node_ids_ghosted_natural(ivertex)-1, &
                     option%myrank,rank,INSERT_VALUES,ierr)
  enddo
  call MatAssemblyBegin(Rank_Mat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Rank_Mat,MAT_FINAL_ASSEMBLY,ierr)
  
#ifdef GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_Rank_Mat.out',viewer,ierr)
  call MatView(Rank_Mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  allocate(val(option%mycommsize))
  allocate(cols(option%mycommsize))
  call MatGetOwnershipRange(Rank_Mat,istart,iend,ierr)
  allocate(int_array(iend-istart))
  count = 0
  do row = istart, iend-1
    call MatGetRow(Rank_Mat,row,ncols,cols,val,ierr)
      max_val = 0.d0
      do local_id = 1, ncols
        max_val = max(max_val,val(local_id))
      enddo
    count = count + 1
    int_array(count) = int(max_val)
    call MatRestoreRow(Rank_Mat,row,ncols,cols,val,ierr)
  enddo
  deallocate(val)
  deallocate(cols)
  
  ! Change rank to start from 0
  int_array = int_array - 1

#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ranks' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do row = 1, count
    write(86,'(i5)') int_array(count)
  enddo  
  close(86)
#endif    
    
  ! Find the local nodes on a process
  do int_rank = 0, option%mycommsize
    vertex_count = 0
    do local_id = 1, count
      if (int_array(local_id) == int_rank) then
        vertex_count = vertex_count + 1
      endif
    enddo
    call MPI_Allreduce(vertex_count,vertex_count2, &
                       ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                       option%mycomm,ierr)
    if (option%myrank == int_rank) geomech_grid%nlmax_node = vertex_count2
    if (geomech_grid%nlmax_node > geomech_grid%ngmax_node) then
      option%io_buffer = 'Error: nlmax_node cannot be greater than' // &
                         ' ngmax_node.'
      call printErrMsg(option)
    endif
  enddo
  
  ! Add a check on nlmax_node to see if there are too many processes 
  call MPI_Allreduce(geomech_grid%nlmax_node,nlmax_node, &
                     ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                     option%mycomm,ierr)
                     
  if (nlmax_node < 1) then
    option%io_buffer = 'Error: Too many processes for the size of the domain.'
    call printErrMsg(option)
  endif
    
  call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                       is_rank,ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_rank_nodes.out', &
                            viewer,ierr)
  call ISView(is_rank,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  deallocate(int_array)
  
  allocate(int_array(count))
  global_offset_old = 0
  call MPI_Exscan(count,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  do local_id = 1, count
    int_array(local_id) = (local_id-1) + global_offset_old
  enddo
  call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                       is_natural,ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_natural.out', &
                            viewer,ierr)
  call ISView(is_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif       

  ! Find the global_offset for vertices on this rank
  global_offset_old = 0
  call MPI_Exscan(geomech_grid%nlmax_node,global_offset_old, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
                  
  geomech_grid%global_offset = global_offset_old
                  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  print *, 'Number of local vertices on process' // &
  trim(adjustl(string)) // ' is:', geomech_grid%nlmax_node
  print *, 'Global offset of vertices on process' // &
  trim(adjustl(string)) // ' is:', global_offset_old  
  print *, 'Number of ghosted vertices on process' // &
  trim(adjustl(string)) // ' is:', geomech_grid%ngmax_node  
#endif

  call ISPartitioningToNumbering(is_rank,is_rank_new,ierr)
  call ISDestroy(is_rank,ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_rank_nodes_new.out', &
                            viewer,ierr)
  call ISView(is_rank_new,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  !Create an application ordering
  call AOCreateBasicIS(is_natural,is_rank_new,ao_natural_to_petsc_nodes,ierr)
  call ISDestroy(is_rank_new,ierr)
  call ISDestroy(is_natural,ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_ao_natural_to_petsc_nodes.out', &
                            viewer,ierr)
  call AOView(ao_natural_to_petsc_nodes,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   

  ! Create a new IS with local PETSc numbering
  allocate(int_array(geomech_grid%nlmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = (local_id-1) + geomech_grid%global_offset
  enddo
  call ISCreateGeneral(option%mycomm,geomech_grid%nlmax_node, &
                       int_array,PETSC_COPY_VALUES,is_natural,ierr)
  deallocate(int_array)
  
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_petsc.out', &
                            viewer,ierr)
  call ISView(is_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   
  
  ! Rewrite the IS with natural numbering
  call AOPetscToApplicationIS(ao_natural_to_petsc_nodes,is_natural,ierr)
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_natural_after_ordering.out', &
                            viewer,ierr)
  call ISView(is_natural,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif    

  ! Get the local indices (natural)
  allocate(geomech_grid%node_ids_local_natural(geomech_grid%nlmax_node))
  call ISGetIndicesF90(is_natural,int_ptr,ierr)
  do local_id = 1, geomech_grid%nlmax_node
    geomech_grid%node_ids_local_natural(local_id) = int_ptr(local_id)
  enddo
  call ISRestoreIndicesF90(is_natural,int_ptr,ierr)
  call ISDestroy(is_natural,ierr)
  
  ! Changing to 1-based
  geomech_grid%node_ids_local_natural = geomech_grid%node_ids_local_natural + 1

  ! Find the natural ids of ghost nodes (vertices)
  if (geomech_grid%ngmax_node - geomech_grid%nlmax_node > 0) then  
    vertex_count = 0 
    allocate(int_array2(geomech_grid%ngmax_node-geomech_grid%nlmax_node))
    do ivertex = 1, geomech_grid%ngmax_node
      do local_id = 1, geomech_grid%nlmax_node
        vertex_found = PETSC_FALSE
        if (geomech_grid%node_ids_ghosted_natural(ivertex) == &
            geomech_grid%node_ids_local_natural(local_id)) then
          vertex_found = PETSC_TRUE
          exit
        endif
       enddo
       if (.not.vertex_found) then
         vertex_count = vertex_count + 1
         int_array2(vertex_count) = geomech_grid%node_ids_ghosted_natural(ivertex)
       endif
    enddo
  endif
  
  if (vertex_count /= geomech_grid%ngmax_node - geomech_grid%nlmax_node) then
    option%io_buffer = 'Error in number of ghost nodes!'
    call printErrMsg(option)
  endif
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ids_ghosts_natural' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  if (allocated(int_array2)) then
    do local_id = 1, vertex_count
      write(86,'(i5)') int_array2(local_id)
    enddo
  else
    write(86,*) 'There are no ghost nodes (vertices) on this process.'
  endif
  close(86)
#endif 
    
  geomech_grid%num_ghost_nodes = geomech_grid%ngmax_node - &
                                   geomech_grid%nlmax_node  
    
  ! Changing the index to 0 based
  if (allocated(int_array2))  int_array2 = int_array2 - 1    
    
  ! Create a new IS with local PETSc numbering
  call ISCreateGeneral(option%mycomm,geomech_grid%num_ghost_nodes, &
                       int_array2,PETSC_COPY_VALUES,is_ghost_petsc,ierr)
  
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_ghost_natural.out', &
                            viewer,ierr)
  call ISView(is_ghost_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   
  
  ! Rewrite the IS with natural numbering
  call AOApplicationtoPetscIS(ao_natural_to_petsc_nodes,is_ghost_petsc,ierr)
#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_ghost_petsc.out', &
                            viewer,ierr)
  call ISView(is_ghost_petsc,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif     
 
  if (allocated(int_array2)) then
    allocate(geomech_grid%ghosted_node_ids_natural(vertex_count))
    ! Change back to 1-based
    geomech_grid%ghosted_node_ids_natural = int_array2 + 1
    deallocate(int_array2)
  endif
  
  ! Get the petsc indices for ghost nodes
  if (geomech_grid%num_ghost_nodes  > 0) then
    allocate(geomech_grid%ghosted_node_ids_petsc(geomech_grid%num_ghost_nodes))
    call ISGetIndicesF90(is_ghost_petsc,int_ptr,ierr)
    do ghosted_id = 1, geomech_grid%num_ghost_nodes
      geomech_grid%ghosted_node_ids_petsc(ghosted_id) = int_ptr(ghosted_id) 
    enddo
    call ISRestoreIndicesF90(is_ghost_petsc,int_ptr,ierr)
  endif
  call ISDestroy(is_ghost_petsc,ierr)
  
  
  ! Changing back to 1-based
  if (geomech_grid%num_ghost_nodes > 0) &
  geomech_grid%ghosted_node_ids_petsc = geomech_grid%ghosted_node_ids_petsc + 1

  geomech_grid%ao_natural_to_petsc_nodes = ao_natural_to_petsc_nodes
  
  ! The following is for re-ordering of local ghosted numbering such that 
  ! the first nlmax_node values are local nodes and the rest are ghost nodes
  allocate(int_array(geomech_grid%ngmax_node))
  allocate(int_array2(geomech_grid%ngmax_node))
  do ivertex = 1, geomech_grid%nlmax_node
    int_array(ivertex) = geomech_grid%node_ids_local_natural(ivertex)
  enddo
  if (geomech_grid%num_ghost_nodes > 0) then
    do ivertex = geomech_grid%nlmax_node+1, geomech_grid%ngmax_node
      int_array(ivertex) = geomech_grid% &
                      ghosted_node_ids_natural(ivertex-geomech_grid%nlmax_node)
    enddo
  endif
  
  do ivertex = 1, geomech_grid%ngmax_node
    int_array2(ivertex) = ivertex
  enddo
  int_array2 = int_array2 - 1
  call PetscSortIntWithPermutation(geomech_grid%ngmax_node,int_array, &
                                   int_array2,ierr)
  int_array2 = int_array2+1
   
  do local_id = 1, geomech_grid%nlmax_elem
    do ivertex = 1, geomech_grid%elem_nodes(0,local_id)
      geomech_grid%elem_nodes(ivertex,local_id) = &
        int_array2(geomech_grid%elem_nodes(ivertex,local_id))
    enddo
  enddo
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_elem_nodes_reordered' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id = 1, geomech_grid%nlmax_elem
    write(86,'(i5)') geomech_grid%elem_nodes(0,local_id)
    do ivertex = 1, geomech_grid%elem_nodes(0,local_id)
      write(86,'(i5)') geomech_grid%elem_nodes(ivertex,local_id)
    enddo
  enddo  
  close(86)
#endif      
  
  allocate(vertices(geomech_grid%ngmax_node))
  do ghosted_id = 1, geomech_grid%ngmax_node
    vertices(ghosted_id)%id = geomech_grid%nodes(ghosted_id)%id
    vertices(ghosted_id)%x = geomech_grid%nodes(ghosted_id)%x
    vertices(ghosted_id)%y = geomech_grid%nodes(ghosted_id)%y
    vertices(ghosted_id)%z = geomech_grid%nodes(ghosted_id)%z
  enddo
  
  do ghosted_id = 1, geomech_grid%ngmax_node
    geomech_grid%nodes(int_array2(ghosted_id))%id = vertices(ghosted_id)%id
    geomech_grid%nodes(int_array2(ghosted_id))%x = vertices(ghosted_id)%x
    geomech_grid%nodes(int_array2(ghosted_id))%y = vertices(ghosted_id)%y
    geomech_grid%nodes(int_array2(ghosted_id))%z = vertices(ghosted_id)%z    
  enddo

#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_coordinates_reordered' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ivertex = 1, vertex_count
    write(86,'(i5)') geomech_grid%nodes(ivertex)%id
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%x
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%y
    write(86,'(1pe12.4)') geomech_grid%nodes(ivertex)%z
  enddo  
  close(86)
#endif  
   
  deallocate(int_array)
  deallocate(vertices)
  
  allocate(int_array(geomech_grid%ngmax_node))
  
  int_array = geomech_grid%node_ids_ghosted_natural
  
  do ghosted_id = 1, geomech_grid%ngmax_node
    geomech_grid%node_ids_ghosted_natural(int_array2(ghosted_id)) = &
                                          int_array(ghosted_id)
  enddo
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ids_ghosted_natural_reordered' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, geomech_grid%ngmax_node
    write(86,'(i5)') geomech_grid%node_ids_ghosted_natural(ghosted_id)
  enddo  
  close(86)
#endif   
  
  deallocate(int_array)
  deallocate(int_array2)
  
  allocate(geomech_grid%gauss_node(geomech_grid%nlmax_elem))
  do local_id = 1, geomech_grid%nlmax_elem
    call GaussInitialize(geomech_grid%gauss_node(local_id))  
    geomech_grid%gauss_node(local_id)%Eletype = geomech_grid%Elem_type(local_id)
    geomech_grid%gauss_node(local_id)%dim = THREE_DIM_GRID
    geomech_grid%gauss_node(local_id)%NGPTS = FOUR_INTEGER
    if (geomech_grid%gauss_node(local_id)%Eletype == PYR_TYPE) &
      geomech_grid%gauss_node(local_id)%NGPTS = FIVE_INTEGER
    call GaussCalculatePoints(geomech_grid%gauss_node(local_id))
  enddo
  
  allocate(int_array(geomech_grid%ngmax_node))
  do local_id = 1, geomech_grid%nlmax_node
    int_array(local_id) = local_id + geomech_grid%global_offset
  enddo  
  do ghosted_id = geomech_grid%nlmax_node+1, geomech_grid%ngmax_node
    int_array(ghosted_id) = geomech_grid%ghosted_node_ids_petsc(ghosted_id - geomech_grid%nlmax_node)
  enddo  
  allocate(geomech_grid%node_ids_ghosted_petsc(geomech_grid%ngmax_node))
  geomech_grid%node_ids_ghosted_petsc = int_array
  deallocate(int_array)

#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ids_ghosted_petsc' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ghosted_id = 1, geomech_grid%ngmax_node
    write(86,'(i5)') geomech_grid%node_ids_ghosted_petsc(ghosted_id)
  enddo  
  close(86)
#endif    
  
end subroutine CopySubsurfaceGridtoGeomechGrid


! ************************************************************************** !
!
! GeomechGridLocalizeRegions: Resticts regions to vertices local 
!                                    to processor for geomech grid when
!                                    the region is defined by a list of 
!                                    vertex ids
! author: Satish Karra
! date: 06/13/13
!
! ************************************************************************** !
subroutine GeomechGridLocalizeRegions(grid,region_list,option)

  use Option_module
  use Geomechanics_Region_module

  implicit none
  
  type(gm_region_list_type), pointer      :: region_list
  type(geomech_grid_type), pointer        :: grid
  type(option_type)                       :: option
  
  type(gm_region_type), pointer           :: region
  character(len=MAXSTRINGLENGTH)          :: string
  
  
  
  
  region => region_list%first
  do
    if (.not.(associated(region))) exit
    
    if (.not.(associated(region%vertex_ids))) then
      option%io_buffer = 'GeomechGridLocalizeRegions: define region only ' // &
                         'by list of vertices is currently implemented: ' //  &
                          trim(region%name)
      call printErrMsg(option)     
    else
      call GeomechGridLocalizeRegFromVertIDs(grid,region,option)
    endif
    
    if (region%num_verts == 0 .and. associated(region%vertex_ids)) then
      deallocate(region%vertex_ids)
      nullify(region%vertex_ids)
    endif
    
    region => region%next
  
  enddo
  

end subroutine GeomechGridLocalizeRegions

! ************************************************************************** !
!
! GeomechGridLocalizeRegFromVertIDs: Resticts regions to vertices local 
!                                    to processor for geomech grid when
!                                    the region is defined by a list of 
!                                    vertex ids
! author: Satish Karra
! date: 06/13/13
!
! ************************************************************************** !
subroutine GeomechGridLocalizeRegFromVertIDs(geomech_grid,geomech_region, &
                                             option)


  use Option_module
  use Geomechanics_Region_module
  
  implicit none
  
#include "finclude/petsclog.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscmat.h"


  type(geomech_grid_type)                 :: geomech_grid
  type(gm_region_type)                    :: geomech_region
  type(option_type)                       :: option
 
  Vec                                     :: vec_vertex_ids,vec_vertex_ids_loc
  IS                                      :: is_from, is_to
  VecScatter                              :: vec_scat
  PetscErrorCode                          :: ierr
  PetscViewer                             :: viewer
  PetscInt                                :: ii,jj,kk,count
  PetscInt                                :: istart,iend
  PetscInt                                :: ghosted_id,local_id
  PetscInt                                :: natural_id
  PetscInt, pointer                       :: tmp_int_array(:)
  PetscScalar, pointer                    :: v_loc_p(:)
  PetscScalar, pointer                    :: tmp_scl_array(:)
  character(len=MAXSTRINGLENGTH)          :: string,string1



  if (associated(geomech_region%vertex_ids)) then
    call VecCreateMPI(option%mycomm,geomech_grid%nlmax_node,PETSC_DECIDE, &
                      vec_vertex_ids,ierr)
    call VecCreateMPI(option%mycomm,geomech_grid%nlmax_node,PETSC_DECIDE,&
                      vec_vertex_ids_loc,ierr)
    call VecZeroEntries(vec_vertex_ids,ierr)
    
    allocate(tmp_int_array(geomech_region%num_verts))
    allocate(tmp_scl_array(geomech_region%num_verts))

    count = 0
    do ii = 1, geomech_region%num_verts
      count = count + 1
      tmp_int_array(count) = geomech_region%vertex_ids(ii) - 1 ! Change to zero-based numbering
      tmp_scl_array(count) = 1.d0
    enddo
    
#ifdef GEOMECH_DEBUG
    call PetscViewerASCIIOpen(option%mycomm,'vec_vertex_ids_bef.out', &
                              viewer,ierr)
    call VecView(vec_vertex_ids,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)  
#endif    

    call VecSetValues(vec_vertex_ids,geomech_region%num_verts,tmp_int_array, &
                      tmp_scl_array,ADD_VALUES,ierr)
    
    deallocate(tmp_int_array)
    deallocate(tmp_scl_array)

    call VecAssemblyBegin(vec_vertex_ids,ierr)
    call VecAssemblyEnd(vec_vertex_ids,ierr)
    
#ifdef GEOMECH_DEBUG
    call PetscViewerASCIIOpen(option%mycomm,'vec_vertex_ids_aft.out', &
                              viewer,ierr)
    call VecView(vec_vertex_ids,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)  
#endif       
    
  endif
  
  allocate(tmp_int_array(geomech_grid%nlmax_node))
  count = 0
  do ghosted_id = 1, geomech_grid%ngmax_node
    local_id = geomech_grid%nG2L(ghosted_id)
    if (local_id < 1) cycle
    count = count + 1
    natural_id = geomech_grid%nG2A(ghosted_id)
    tmp_int_array(count) = natural_id
  enddo

  tmp_int_array = tmp_int_array - 1
  call ISCreateBlock(option%mycomm,1,geomech_grid%nlmax_node, &
                     tmp_int_array,PETSC_COPY_VALUES,is_from,ierr)

  call VecGetOwnershipRange(vec_vertex_ids_loc,istart,iend,ierr)
  do ii = 1,geomech_grid%nlmax_node
    tmp_int_array(ii) = ii + istart
  enddo
 
  ! is_from is natural_numbering
  ! is_to is PETSc_numbering
 
  tmp_int_array = tmp_int_array - 1 
  call ISCreateBlock(option%mycomm,1,geomech_grid%nlmax_node,&
                     tmp_int_array,PETSC_COPY_VALUES,is_to,ierr)

  deallocate(tmp_int_array)
  
  call VecScatterCreate(vec_vertex_ids,is_from,vec_vertex_ids_loc,is_to, &
                        vec_scat,ierr)

  call ISDestroy(is_from,ierr)
  call ISDestroy(is_to,ierr)
  
  call VecScatterBegin(vec_scat,vec_vertex_ids,vec_vertex_ids_loc, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(vec_scat,vec_vertex_ids,vec_vertex_ids_loc, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterDestroy(vec_scat,ierr)
  
#if GEOMECH_DEBUG
    call PetscViewerASCIIOpen(option%mycomm,'vec_vertex_ids_loc.out', &
                              viewer,ierr)
    call VecView(vec_vertex_ids_loc,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)
#endif

  call VecGetArrayF90(vec_vertex_ids_loc,v_loc_p,ierr)
  count = 0
  do ii = 1,geomech_grid%nlmax_node
    if (v_loc_p(ii) == 1) count = count + 1
  enddo
  
  geomech_region%num_verts = count
  if (count > 0) then
    allocate(tmp_int_array(count))
    count = 0
    do ii = 1,geomech_grid%nlmax_node
      if (v_loc_p(ii) == 1) then
        count = count + 1
        tmp_int_array(count) = ii
      endif
    enddo
    
    deallocate(geomech_region%vertex_ids)
    allocate(geomech_region%vertex_ids(geomech_region%num_verts))
    geomech_region%vertex_ids = tmp_int_array
    deallocate(tmp_int_array)
  endif
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  write(string1,*) geomech_region%name
  string = 'vec_region_' // trim(adjustl(string1)) //  &
           '_mapped' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ii = 1,geomech_region%num_verts
    write(86,'(i5)') geomech_region%vertex_ids(ii)
  enddo  
  close(86)
#endif     
  
  call VecRestoreArrayF90(vec_vertex_ids_loc,v_loc_p,ierr)
  
  call VecDestroy(vec_vertex_ids,ierr)
  call VecDestroy(vec_vertex_ids_loc,ierr)
   
end subroutine GeomechGridLocalizeRegFromVertIDs

! ************************************************************************** !
!
! GeomechGridCopyIntegerArrayToVec: Copies values from an integer array into a 
!                                 PETSc Vec
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechGridCopyIntegerArrayToVec(grid,array,vector,num_values)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(geomech_grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call GeomechGridVecGetArrayF90(grid,vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call GeomechGridVecRestoreArrayF90(grid, vector,vec_ptr,ierr)
  
end subroutine GeomechGridCopyIntegerArrayToVec

! ************************************************************************** !
!
! GeomechGridVecGetArrayF90: Returns pointer to geomech veretex-based vector 
!                            values
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechGridVecGetArrayF90(grid,vec,f90ptr,ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_grid_type) :: grid
  Vec :: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(vec,f90ptr,ierr)

end subroutine GeomechGridVecGetArrayF90

! ************************************************************************** !
!
! GeomechGridVecRestoreArrayF90: Restores pointer to geomech vector values
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechGridVecRestoreArrayF90(grid,vec,f90ptr,ierr)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_grid_type) :: grid
  Vec :: vec
  PetscReal, pointer :: f90ptr(:)
  PetscErrorCode :: ierr

  call VecRestoreArrayF90(vec,f90ptr,ierr)
  
end subroutine GeomechGridVecRestoreArrayF90

! ************************************************************************** !
!
! GeomechGridCopyVecToIntegerArray: Copies values from a PETSc Vec to an  
!                                  integer array
! author: Satish Karra, LANL
! date: 06/17/13
!
! ************************************************************************** !
subroutine GeomechGridCopyVecToIntegerArray(grid,array,vector,num_values)

  implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type(geomech_grid_type) :: grid
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call GeomechGridVecGetArrayF90(grid,vector,vec_ptr,ierr)
  do i = 1,num_values
    if (vec_ptr(i) > 0.d0) then
      array(i) = int(vec_ptr(i)+1.d-4)
    else
      array(i) = int(vec_ptr(i)-1.d-4)
    endif
  enddo
  call GeomechGridVecRestoreArrayF90(grid,vector,vec_ptr,ierr)
  
end subroutine GeomechGridCopyVecToIntegerArray

end module Geomechanics_Grid_module
#endif 
! GEOMECH