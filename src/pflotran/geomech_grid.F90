#ifdef GEOMECH
module Geomech_Grid_module

  use Geomech_Grid_Aux_module
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

  !  PetscInt, parameter :: HEX_TYPE          = 1
  !  PetscInt, parameter :: TET_TYPE          = 2
  !  PetscInt, parameter :: WEDGE_TYPE        = 3
  !  PetscInt, parameter :: PYR_TYPE          = 4
  !  PetscInt, parameter :: TRI_FACE_TYPE     = 1
  !  PetscInt, parameter :: QUAD_FACE_TYPE    = 2
  !  PetscInt, parameter :: MAX_VERT_PER_FACE = 4

  public :: GeomechGridRead, &
            CopySubsurfaceGridtoGeomechGrid 

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
  use Geomech_Grid_Aux_module
  use Option_module
  
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
  type(geomech_Grid_type), pointer           :: geomech_grid
  type(option_type), pointer                 :: option
  PetscInt                                   :: local_id
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
  Mat                                        :: Rank_Mat
  PetscReal                                  :: rank
  PetscViewer                                :: viewer
  Vec                                        :: node_rank
  PetscReal, pointer                         :: vec_ptr(:) 
  PetscInt                                   :: istart,iend
  PetscBool                                  :: vertex_found
  PetscInt                                   :: int_rank
  PetscInt                                   :: vertex_count2
  IS                                         :: is_rank 

  
  call printMsg(option,'GEOMECHANICS: Subsurface unstructured grid will ' // &
                  'be used for geomechanics.')
  
#ifdef GEOMECH_DEBUG
  call printMsg(option,'Copying unstructured grid to geomechanics grid')
#endif
  
  geomech_grid%global_offset = ugrid%global_offset
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
  
  allocate(geomech_grid%node_ids_natural(vertex_count))
  do ivertex = 1, vertex_count
    geomech_grid%node_ids_natural(ivertex) = ugrid% &
      vertex_ids_natural(int_array3(ivertex))  
  enddo

#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  string = 'geomech_node_ids_natural' // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do local_id = 1, vertex_count
    write(86,'(i5)') geomech_grid%node_ids_natural(local_id)
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
    geomech_grid%nodes(ivertex)%id = geomech_grid%node_ids_natural(ivertex)
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
    call MatSetValue(Rank_Mat,geomech_grid%node_ids_natural(ivertex)-1, &
                     option%myrank,rank,INSERT_VALUES,ierr)
  enddo
  call MatAssemblyBegin(Rank_Mat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Rank_Mat,MAT_FINAL_ASSEMBLY,ierr)
  
#ifdef GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_Rank_Mat.out',viewer,ierr)
  call MatView(Rank_Mat,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  

  call VecCreate(option%mycomm,node_rank,ierr)
  call VecSetSizes(node_rank,PETSC_DECIDE,geomech_grid%nmax_node,ierr)
  call VecSetFromOptions(node_rank,ierr)
  call MatGetRowMax(Rank_Mat,node_rank,PETSC_NULL_INTEGER,ierr)
  
  ! Change rank to start from 0
  call VecGetOwnershipRange(node_rank,istart,iend,ierr)
  call VecGetArrayF90(node_rank,vec_ptr,ierr)
  do local_id = 1,iend-istart
    vec_ptr(local_id) = vec_ptr(local_id) - 1
  enddo
  call VecRestoreArrayF90(node_rank,vec_ptr,ierr)

#ifdef GEOMECH_DEBUG
  string = 'geomech_node_assigned_rank.out'
  call PetscViewerASCIIOpen(option%mycomm,trim(string),viewer,ierr)
  call VecView(node_rank,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif   
  
  ! Find the local nodes on a process
  call VecGetOwnershipRange(node_rank,istart,iend,ierr)
  call VecGetArrayF90(node_rank,vec_ptr,ierr)
  do int_rank = 0, option%mycommsize
    vertex_count = 0
    do local_id = 1, iend-istart
      if (vec_ptr(local_id) == int_rank) then
        vertex_count = vertex_count + 1
      endif
    enddo
    call MPI_Allreduce(vertex_count,vertex_count2, &
                       ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                       option%mycomm,ierr)
    if (vertex_count2> geomech_grid%ngmax_node) then
      option%io_buffer = 'Error: nlmax_node cannot be greater than' // &
                         ' ngmax_node.'
      call printErrMsg(option)
    endif
    if (option%myrank == int_rank) geomech_grid%nlmax_node = vertex_count2
  enddo
  
  allocate(int_array(geomech_grid%nlmax_node))
  count = 0
  do local_id = 1, iend-istart
    int_array(local_id) = int(vec_ptr(local_id)) 
    count = count + 1
  enddo
    
  call VecRestoreArrayF90(node_rank,vec_ptr,ierr)
  
  call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES,is_rank,ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'geomech_is_rank_nodes.out', &
                            viewer,ierr)
  call ISView(is_rank,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)
#endif  
  deallocate(int_array)
      
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
  print *, 'Number of ghosted of vertices on process' // &
  trim(adjustl(string)) // ' is:', geomech_grid%ngmax_node  
#endif


  vertex_count = 0 
  allocate(int_array3(geomech_grid%ngmax_node-geomech_grid%nlmax_node))
  do ivertex = 1, geomech_grid%ngmax_node
    do local_id = 1, geomech_grid%nlmax_node
      vertex_found = PETSC_FALSE
      if (geomech_grid%node_ids_natural(ivertex) == int_array2(local_id)) then
        vertex_found = PETSC_TRUE
        exit
      endif
     enddo
     if (.not.vertex_found) then
       vertex_count = vertex_count + 1
       int_array3(vertex_count) = geomech_grid%node_ids_natural(ivertex)
     endif
  enddo

end subroutine CopySubsurfaceGridtoGeomechGrid

end module Geomech_Grid_module
#endif 
! GEOMECH