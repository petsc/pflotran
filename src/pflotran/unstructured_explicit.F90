module Unstructured_Explicit_module
  
  use Geometry_module

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

  type, public :: unstructured_explicit_type
    PetscInt, pointer :: cell_ids(:)
    PetscReal, pointer :: cell_volumes(:)
    type(point3d_type), pointer :: cell_centroids(:)
    PetscInt, pointer :: connections(:,:)
    PetscReal, pointer :: face_areas(:)
    type(point3d_type), pointer :: face_centroids(:)
  end type unstructured_explicit_type

  public :: ExplicitUGridCreate, &
            ExplicitUGridRead, &
            ExplicitUGridDecompose, &
            ExplicitUGridSetInternConnect, &
            ExplicitUGridSetCellCentroids, &
            ExplicitUGridComputeVolumes, &
            ExplicitUGridSetBoundaryConnect, &
            ExplicitUGridSetConnections, &
            ExplicitUGridDestroy

contains

! ************************************************************************** !
!
! ExplicitUGridCreate: Creates an explicit unstructured grid object
! author: Glenn Hammond
! date: 05/14/12
!
! ************************************************************************** !
function ExplicitUGridCreate()

  implicit none
  
  type(unstructured_explicit_type), pointer :: ExplicitUGridCreate

  type(unstructured_explicit_type), pointer :: explicit_grid

  allocate(explicit_grid)

  nullify(explicit_grid%cell_ids)
  nullify(explicit_grid%cell_volumes)
  nullify(explicit_grid%cell_centroids)
  nullify(explicit_grid%connections)
  nullify(explicit_grid%face_areas)
  nullify(explicit_grid%face_centroids)

  ExplicitUGridCreate => explicit_grid
  
end function ExplicitUGridCreate

! ************************************************************************** !
!
! ExplicitUGridRead: Reads an explicit unstructured grid
! author: Glenn Hammond
! date: 05/14/12
!
! ************************************************************************** !
subroutine ExplicitUGridRead(explicit_grid,filename,option)

  use Input_module
  use Option_module
  use String_module
  
  implicit none
  
  type(unstructured_explicit_type) :: explicit_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: card, word
  PetscInt :: fileid, icell, iconn, id_up, id_dn
  
  PetscInt :: num_cells
  PetscInt :: num_connections
  
  fileid = 86
  input => InputCreate(fileid,filename,option)

! Format of explicit unstructured grid file
! id_, id_up_, id_dn_ = integer
! x_, y_, z_, area_, volume_ = real
! -----------------------------------------------------------------
! CELLS <integer>    integer = # cells (N)
! id_1 x_1 y_1 z_1 volume_1
! id_2 x_2 y_2 z_2 volume_2
! ...
! ...
! id_N x_N y_N z_N volume_N
! CONNECTIONS <integer>   integer = # connections (M)
! id_up_1 id_dn_1 x_1 y_1 z_1 area_1
! id_up_2 id_dn_2 x_2 y_2 z_2 area_2
! ...
! ...
! id_up_M id_dn_M x_M y_M z_M area_M
! -----------------------------------------------------------------


  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)
  
    select case(word)
      case('CELLS')
        card = 'Explicit Unstructured Grid CELLS'
        call InputReadInt(input,option,num_cells)
        call InputErrorMsg(input,option,'number of cells',card)
        allocate(explicit_grid%cell_ids(num_cells))
        explicit_grid%cell_ids = 0
        allocate(explicit_grid%cell_volumes(num_cells))
        explicit_grid%cell_volumes = 0
        allocate(explicit_grid%cell_centroids(num_cells))
        do icell = 1, num_cells
          explicit_grid%cell_centroids(icell)%x = 0.d0
          explicit_grid%cell_centroids(icell)%y = 0.d0
          explicit_grid%cell_centroids(icell)%z = 0.d0
        enddo
        do icell = 1, num_cells
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)  
          call InputReadInt(input,option,explicit_grid%cell_ids(icell))
          call InputErrorMsg(input,option,'cell id',card)
          call InputReadDouble(input,option, &
                               explicit_grid%cell_centroids(icell)%x)
          call InputErrorMsg(input,option,'cell x coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%cell_centroids(icell)%y)
          call InputErrorMsg(input,option,'cell y coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%cell_centroids(icell)%z)
          call InputErrorMsg(input,option,'cell z coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%cell_volumes(icell))
          call InputErrorMsg(input,option,'cell volume',card)
        enddo
      case('CONNECTIONS')
        card = 'Explicit Unstructured Grid CONNECTIONS'
        call InputReadInt(input,option,num_connections)
        call InputErrorMsg(input,option,'number of connections',card)
        allocate(explicit_grid%connections(2,num_connections))
        explicit_grid%connections = 0
        allocate(explicit_grid%face_areas(num_connections))
        explicit_grid%face_areas = 0    
        allocate(explicit_grid%face_centroids(num_connections))
        do iconn = 1, num_connections
          explicit_grid%face_centroids(iconn)%x = 0.d0
          explicit_grid%face_centroids(iconn)%y = 0.d0
          explicit_grid%face_centroids(iconn)%z = 0.d0
        enddo
        do iconn = 1, num_connections
          call InputReadFlotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)  
          call InputReadInt(input,option, &
                            explicit_grid%connections(1,iconn))
          call InputErrorMsg(input,option,'cell id upwind',card)
          call InputReadInt(input,option, &
                            explicit_grid%connections(2,iconn))
          call InputErrorMsg(input,option,'cell id downwind',card)
          call InputReadDouble(input,option, &
                               explicit_grid%face_centroids(iconn)%x)
          call InputErrorMsg(input,option,'face x coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%face_centroids(iconn)%y)
          call InputErrorMsg(input,option,'face y coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%face_centroids(iconn)%z)
          call InputErrorMsg(input,option,'face z coordinate',card)
          call InputReadDouble(input,option, &
                               explicit_grid%face_areas(iconn))
          call InputErrorMsg(input,option,'face area',card)
        enddo
      case default
        option%io_buffer = 'Keyword: ' // trim(word) // &
                           ' not recognized while reading explicit ' // &
                           'unstructured grid.'
        call printErrMsg(option)
    end select
  enddo

  call InputDestroy(input)

end subroutine ExplicitUGridRead

! ************************************************************************** !
!
! ExplicitUGridDecompose: Decomposes an explicit unstructured grid across 
!                         ranks
! author: Glenn Hammond
! date: 05/17/12
!
! ************************************************************************** !
subroutine ExplicitUGridDecompose(explicit_grid, num_ghost_cells, &
                                  global_offset, nmax, nlmax, ngmax, &
                                  cell_ids_natural, cell_ids_petsc, &
                                  ghost_cell_ids_petsc, ao_natural_to_petsc, &
                                  option)
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
  
  type(unstructured_explicit_type) :: explicit_grid
  PetscInt :: num_ghost_cells
  PetscInt :: global_offset
  PetscInt :: nmax, nlmax, ngmax
  PetscInt, pointer :: cell_ids_natural(:)
  PetscInt, pointer :: cell_ids_petsc(:)
  PetscInt, pointer :: ghost_cell_ids_petsc(:)
  AO :: ao_natural_to_petsc
  type(option_type) :: option
  
  PetscInt :: num_cells_local_new
  PetscInt :: global_offset_new
  PetscInt :: local_id
  PetscInt, allocatable :: int_array(:)
  PetscErrorCode :: ierr

  
  num_cells_local_new = size(explicit_grid%cell_ids)
  global_offset_new = 0
  
  allocate(cell_ids_natural(num_cells_local_new))
  cell_ids_natural = explicit_grid%cell_ids

  ! make a list of petsc ids for each local cell (you simply take the global 
  ! offset and add it to the local contiguous cell ids on each processor
  allocate(int_array(num_cells_local_new))
  do local_id = 1, num_cells_local_new
    int_array(local_id) = local_id+global_offset_new
  enddo
  
  ! make the arrays zero-based
  int_array = int_array - 1
  cell_ids_natural = cell_ids_natural - 1
  ! create an application ordering (mapping of natural to petsc ordering)
  call AOCreateBasic(option%mycomm,num_cells_local_new, &
                     cell_ids_natural,int_array, &
                     ao_natural_to_petsc,ierr)
  deallocate(int_array)
  ! make cell_ids_natural 1-based again
  cell_ids_natural = cell_ids_natural + 1
  
  allocate(cell_ids_petsc(num_cells_local_new))
  cell_ids_petsc = cell_ids_natural
  
  nmax = num_cells_local_new
  nlmax = nmax
  ngmax = nmax

  num_ghost_cells = 0
  global_offset = global_offset_new
  
end subroutine ExplicitUGridDecompose

! ************************************************************************** !
!
! ExplicitUGridSetCellCentroids: Sets the centroid of each grid cell
! author: Glenn Hammond
! date: 05/17/12
!
! ************************************************************************** !
subroutine ExplicitUGridSetCellCentroids(explicit_grid,x,y,z, &
                                         x_min,x_max,y_min,y_max,z_min,z_max)

  use Option_module

  implicit none
  
  type(unstructured_explicit_type) :: explicit_grid
  PetscReal :: x(:), y(:), z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

  PetscInt :: icell
  
  do icell = 1, size(x)
    x(icell) = explicit_grid%cell_centroids(icell)%x
    y(icell) = explicit_grid%cell_centroids(icell)%y
    z(icell) = explicit_grid%cell_centroids(icell)%z
  enddo
  
  x_min = minval(x)
  x_max = maxval(x)
  y_min = minval(y)
  y_max = maxval(y)
  z_min = minval(z)
  z_max = maxval(z)
      
end subroutine ExplicitUGridSetCellCentroids
      
! ************************************************************************** !
!
! ExplicitUGridSetInternConnect: Sets up the internal connectivity within  
!                                the connectivity object
! author: Glenn Hammond
! date: 05/17/12
!
! ************************************************************************** !
function ExplicitUGridSetInternConnect(explicit_grid,option)

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: ExplicitUGridSetInternConnect
  
  type(unstructured_explicit_type) :: explicit_grid
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id_up, id_dn
  PetscReal :: v(3), v_up(3), v_dn(3)
  PetscReal :: distance
  
  
  num_connections = size(explicit_grid%connections,2)
  connections => ConnectionCreate(num_connections,INTERNAL_CONNECTION_TYPE)
  
  do iconn = 1, num_connections
    id_up = explicit_grid%connections(1,iconn)
    id_dn = explicit_grid%connections(2,iconn)
    connections%id_up(iconn) = id_up
    connections%id_dn(iconn) = id_dn
    
    v_up(1) = explicit_grid%face_centroids(iconn)%x - &
              explicit_grid%cell_centroids(id_up)%x
    v_up(2) = explicit_grid%face_centroids(iconn)%y - &
              explicit_grid%cell_centroids(id_up)%y
    v_up(3) = explicit_grid%face_centroids(iconn)%z - &
              explicit_grid%cell_centroids(id_up)%z

    v_dn(1) = explicit_grid%cell_centroids(id_dn)%x - &
              explicit_grid%face_centroids(iconn)%x
    v_dn(2) = explicit_grid%cell_centroids(id_dn)%y - &
              explicit_grid%face_centroids(iconn)%y
    v_dn(3) = explicit_grid%cell_centroids(id_dn)%z - &
              explicit_grid%face_centroids(iconn)%z

    v = v_up + v_dn
    distance = sqrt(DotProduct(v,v))
    connections%dist(-1,iconn) = sqrt(DotProduct(v_up,v_up))/distance
    connections%dist(0,iconn) = distance
    connections%dist(1:3,iconn) = v/distance
    connections%area(iconn) = explicit_grid%face_areas(iconn)
  enddo
  
  ExplicitUGridSetInternConnect => connections

end function ExplicitUGridSetInternConnect

! ************************************************************************** !
!
! ExplicitUGridComputeVolumes: Sets the volume of each grid cell
! author: Glenn Hammond
! date: 05/17/12
!
! ************************************************************************** !
subroutine ExplicitUGridComputeVolumes(explicit_grid,option,volume)

  use Option_module

  implicit none
  
  type(unstructured_explicit_type) :: explicit_grid
  type(option_type) :: option
  Vec :: volume
  
  PetscInt :: icell
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayF90(volume,vec_ptr,ierr)
  do icell = 1, size(explicit_grid%cell_volumes)
    vec_ptr(icell) = explicit_grid%cell_volumes(icell)
  enddo
  call VecRestoreArrayF90(volume,vec_ptr,ierr)
  
end subroutine ExplicitUGridComputeVolumes

! ************************************************************************** !
!
! ExplicitUGridSetBoundaryConnect: Sets up the boundary connectivity within  
!                                  the connectivity object
! author: Glenn Hammond
! date: 05/18/12
!
! ************************************************************************** !
function ExplicitUGridSetBoundaryConnect(explicit_grid,cell_ids, &
                                         face_centroids,face_areas,option)

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: ExplicitUGridSetBoundaryConnect

  type(unstructured_explicit_type) :: explicit_grid
  PetscInt :: cell_ids(:)
  type(point3d_type) :: face_centroids(:)
  PetscReal :: face_areas(:)
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id
  PetscReal :: v(3)
  PetscReal :: distance
  
  
  num_connections = size(cell_ids)
  connections => ConnectionCreate(num_connections,BOUNDARY_CONNECTION_TYPE)
  
  do iconn = 1, num_connections
    id = cell_ids(iconn)
    connections%id_dn(iconn) = id
    
    v(1) = explicit_grid%cell_centroids(id)%x - &
           face_centroids(iconn)%x
    v(2) = explicit_grid%cell_centroids(id)%y - &
           face_centroids(iconn)%y
    v(3) = explicit_grid%cell_centroids(id)%z - &
           face_centroids(iconn)%z

    distance = sqrt(DotProduct(v,v))
    connections%dist(-1,iconn) = 0.d0
    connections%dist(0,iconn) = distance
    connections%dist(1:3,iconn) = v/distance
    connections%area(iconn) = face_areas(iconn)
  enddo
  
  ExplicitUGridSetBoundaryConnect => connections

end function ExplicitUGridSetBoundaryConnect

! ************************************************************************** !
!
! ExplicitUGridSetConnections: Sets up the connectivity for a region
! author: Glenn Hammond
! date: 05/18/12
!
! ************************************************************************** !
function ExplicitUGridSetConnections(explicit_grid,cell_ids,connection_type, &
                                     option)

  use Utility_module
  use Connection_module
  use Option_module

  implicit none
  
  type(connection_set_type), pointer :: ExplicitUGridSetConnections

  type(unstructured_explicit_type) :: explicit_grid
  PetscInt :: cell_ids(:)
  PetscInt :: connection_type
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connections
  PetscInt :: num_connections
  PetscInt :: iconn
  PetscInt :: id
  
  num_connections = size(cell_ids)
  connections => ConnectionCreate(num_connections,connection_type)
  
  do iconn = 1, num_connections
    id = cell_ids(iconn)
    connections%id_dn(iconn) = id
  enddo
  
  ExplicitUGridSetConnections => connections

end function ExplicitUGridSetConnections

! ************************************************************************** !
!
! ExplicitUGridDestroy: Deallocates an explicit unstructured grid object
! author: Glenn Hammond
! date: 05/14/12
!
! ************************************************************************** !
subroutine ExplicitUGridDestroy(explicit_grid)

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(unstructured_explicit_type), pointer :: explicit_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(explicit_grid)) return

  call DeallocateArray(explicit_grid%cell_ids)
  call DeallocateArray(explicit_grid%cell_volumes)
  if (associated(explicit_grid%cell_centroids)) &
    deallocate(explicit_grid%cell_centroids)
  nullify(explicit_grid%cell_centroids)
  call DeallocateArray(explicit_grid%connections)
  call DeallocateArray(explicit_grid%face_areas)
  if (associated(explicit_grid%face_centroids)) &
    deallocate(explicit_grid%face_centroids)
  nullify(explicit_grid%face_centroids)
  
  deallocate(explicit_grid)
  nullify(explicit_grid)

end subroutine ExplicitUGridDestroy

end module Unstructured_Explicit_module
