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
! date: 05/23/13
!
! ************************************************************************** !
subroutine CopySubsurfaceGridtoGeomechGrid(ugrid,geomech_grid,option)
                                        
  use Unstructured_Grid_Aux_module
  use Geomech_Grid_Aux_module
  use Option_module
  
  implicit none
  
  type(unstructured_grid_type), pointer      :: ugrid
  type(geomech_Grid_type), pointer           :: geomech_grid
  type(option_type), pointer                 :: option
  
  if (option%myrank == option%io_rank) then
    write(*,*),'GEOMECHANICS: Subsurface unstructured grid will be used for '// &
               'geomechanics.'
  endif  
  
  geomech_grid%global_offset = ugrid%global_offset
  geomech_grid%nmax_elem = ugrid%nmax
  geomech_grid%nlmax_elem = ugrid%nlmax
  geomech_grid%nmax_node = ugrid%num_vertices_global
  geomech_grid%nlmax_node = ugrid%num_vertices_local
  allocate(geomech_grid%elem_ids_natural(size(ugrid%cell_ids_natural)))
  geomech_grid%elem_ids_natural = ugrid%cell_ids_natural
  allocate(geomech_grid%elem_ids_petsc(size(ugrid%cell_ids_petsc)))
  geomech_grid%elem_ids_petsc = ugrid%cell_ids_petsc
  geomech_grid%ao_natural_to_petsc = ugrid%ao_natural_to_petsc
  geomech_grid%max_ndual_per_elem = ugrid%max_ndual_per_cell
  geomech_grid%max_nnode_per_elem = ugrid%max_nvert_per_cell
  geomech_grid%max_elem_sharing_a_node = ugrid%max_cells_sharing_a_vertex
  allocate(geomech_grid%elem_type(size(ugrid%cell_type)))
  geomech_grid%elem_type = ugrid%cell_type
  allocate(geomech_grid%elem_nodes(size(ugrid%cell_vertices,ONE_INTEGER), &
           size(ugrid%cell_vertices,TWO_INTEGER)))
  allocate(geomech_grid%nodes(size(ugrid%vertices)))
  geomech_grid%nodes = ugrid%vertices
  geomech_grid%elem_nodes = ugrid%cell_vertices
  
end subroutine CopySubsurfaceGridtoGeomechGrid

end module Geomech_Grid_module
#endif 
! GEOMECH