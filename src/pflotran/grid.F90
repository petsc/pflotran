module Grid_module

  use Structured_Grid_module
  use Unstructured_Grid_module
  use Connection_module
 
  implicit none

  private
 
#include "definitions.h"

  type, public :: grid_type 
  
    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured, unstructured, etc.)
    
    PetscInt :: nmax   ! Total number of nodes in global domain
    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
   
    !nL2G :  not collective, local processor: local  =>  ghosted local  
    !nG2L :  not collective, local processor:  ghosted local => local  
    !nG2N :  collective,  ghosted local => global index , used for   
    !                     matsetvaluesblocked ( not matsetvaluesblockedlocal)  
    !nL2A :   collective, local => natural index, used for initialization   
    !                              and source/sink setup  
    PetscInt, pointer :: nL2G(:), nG2L(:), nL2A(:)
    PetscInt, pointer :: nG2A(:)
    
    PetscReal, pointer :: x(:), y(:), z(:) ! coordinates of ghosted grid cells
    
    PetscReal :: x_min_global, x_max_global, y_min_global, y_max_global, z_min_global, z_max_global
    PetscReal :: x_min_local, x_max_local, y_min_local, y_max_local, z_min_local, z_max_local

    PetscInt, pointer :: hash(:,:,:)
    PetscInt :: num_hash_bins

    type(structured_grid_type), pointer :: structured_grid
    type(unstructured_grid_type), pointer :: unstructured_grid
    
    type(connection_set_list_type), pointer :: internal_connection_set_list

  end type grid_type


  public :: GridCreate, &
            GridDestroy, &
            GridComputeInternalConnect, &
            GridMapIndices, &
            GridComputeSpacing, &
            GridComputeCoordinates, &
            GridComputeVolumes, &
            GridLocalizeRegions, &
            GridPopulateConnection, &
            GridCopyIntegerArrayToPetscVec, &
            GridCopyRealArrayToPetscVec, &
            GridCopyPetscVecToIntegerArray, &
            GridCopyPetscVecToRealArray, &
            GridCreateNaturalToGhostedHash, &
            GridDestroyHashTable, &
            GridGetLocalGhostedIdFromHash, &
            GridVecGetArrayF90, &
            GridVecRestoreArrayF90
contains

! ************************************************************************** !
!
! GridCreate: Creates a structured or unstructured grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function GridCreate()

  implicit none
  
  type(grid_type), pointer :: GridCreate
  
  type(grid_type), pointer :: grid
  
  allocate(grid)
  grid%ctype = ''
  grid%itype = 0

  nullify(grid%structured_grid)
  nullify(grid%unstructured_grid)

  nullify(grid%internal_connection_set_list)

  nullify(grid%nL2G)
  nullify(grid%nG2L)
  nullify(grid%nL2A)
  nullify(grid%nG2A)

  nullify(grid%x)
  nullify(grid%y)
  nullify(grid%z)

  grid%x_min_global = 1.d20
  grid%x_max_global = -1.d20
  grid%y_min_global = 1.d20
  grid%y_max_global = -1.d20
  grid%z_min_global = 1.d20
  grid%z_max_global = -1.d20

  grid%x_min_local = 1.d20
  grid%x_max_local = -1.d20
  grid%y_min_local = 1.d20
  grid%y_max_local = -1.d20
  grid%z_min_local = 1.d20
  grid%z_max_local = -1.d20

  grid%nmax = 0
  grid%nlmax = 0 
  grid%ngmax = 0
  
  nullify(grid%hash)
  grid%num_hash_bins = 1000

  GridCreate => grid

end function GridCreate

! ************************************************************************** !
!
! GridComputeInternalConnect: computes internal connectivity of a grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
subroutine GridComputeInternalConnect(grid,option)

  use Connection_module
  use Option_module    
    
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  type(connection_set_type), pointer :: connection_set
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      connection_set => &
        StructGridComputeInternConnect(grid%structured_grid,option)
    case(UNSTRUCTURED_GRID) 
      connection_set => &
        UnstGridComputeInternConnect(grid%unstructured_grid,option)
  end select
  
  allocate(grid%internal_connection_set_list)
  call ConnectionInitList(grid%internal_connection_set_list)
  call ConnectionAddToList(connection_set,grid%internal_connection_set_list)
  
end subroutine GridComputeInternalConnect

! ************************************************************************** !
!
! GridPopulateConnection: computes connectivity coupler to a grid
! author: Glenn Hammond
! date: 11/09/07
!
! ************************************************************************** !
subroutine GridPopulateConnection(grid,connection,iface,iconn,cell_id_local)

  use Connection_module
  use Structured_Grid_module
  
  implicit none
 
  type(grid_type) :: grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: cell_id_local
  
  PetscInt :: cell_id_ghosted
  
  cell_id_ghosted = grid%nL2G(cell_id_local)
  ! Use ghosted index to access dx, dy, dz because we have
  ! already done a global-to-local scatter for computing the
  ! interior node connections.
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructGridPopulateConnection(grid%structured_grid,connection, &
                                        iface,iconn,cell_id_ghosted)
    case(UNSTRUCTURED_GRID)
  end select

end subroutine GridPopulateConnection

! ************************************************************************** !
!
! GridMapIndices: maps global, local and natural indices of cells 
!                 to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridMapIndices(grid)

  implicit none
  
  type(grid_type) :: grid
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructuredGridMapIndices(grid%structured_grid,grid%nG2L,grid%nL2G, &
                                    grid%nL2A,grid%nG2A)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine GridMapIndices

! ************************************************************************** !
!
! GridComputeSpacing: Computes grid spacing (only for structured grid
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine GridComputeSpacing(grid)

  implicit none
  
  type(grid_type) :: grid
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructuredGridComputeSpacing(grid%structured_grid,grid%nG2A, &
                                        grid%nG2L)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine GridComputeSpacing

! ************************************************************************** !
!
! GridComputeCoordinates: Computes x,y,z coordinates of grid cells
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine GridComputeCoordinates(grid,origin_global,option)

  use Option_module
  
  implicit none

  type(grid_type) :: grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  
  PetscErrorCode :: ierr  
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      allocate(grid%x(grid%ngmax))
      grid%x = 0.d0
      allocate(grid%y(grid%ngmax))
      grid%y = 0.d0
      allocate(grid%z(grid%ngmax))
      grid%z = 0.d0
      call StructuredGridComputeCoord(grid%structured_grid,option, &
                                      origin_global, &
                                      grid%x,grid%y,grid%z, &
                                      grid%x_min_local,grid%x_max_local, &
                                      grid%y_min_local,grid%y_max_local, &
                                      grid%z_min_local,grid%z_max_local)
    case(UNSTRUCTURED_GRID)
  end select

  if((grid%itype==STRUCTURED_GRID).and.(grid%structured_grid%p_samr_patch==0)) then
     ! compute global max/min from the local max/in
     call MPI_Allreduce(grid%x_min_local,grid%x_min_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MIN,PETSC_COMM_WORLD,ierr)
     call MPI_Allreduce(grid%y_min_local,grid%y_min_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MIN,PETSC_COMM_WORLD,ierr)
     call MPI_Allreduce(grid%z_min_local,grid%z_min_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MIN,PETSC_COMM_WORLD,ierr)
     call MPI_Allreduce(grid%x_max_local,grid%x_max_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MAX,PETSC_COMM_WORLD,ierr)
     call MPI_Allreduce(grid%y_max_local,grid%y_max_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MAX,PETSC_COMM_WORLD,ierr)
     call MPI_Allreduce(grid%z_max_local,grid%z_max_global,ONE_INTEGER, &
          MPI_DOUBLE_PRECISION,MPI_MAX,PETSC_COMM_WORLD,ierr)
  endif
end subroutine GridComputeCoordinates

! ************************************************************************** !
!
! GridComputeVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine GridComputeVolumes(grid,volume,option)

  use Option_module
  
  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  type(grid_type) :: grid
  type(option_type) :: option
  Vec :: volume
  
  select case(grid%itype)
    case(STRUCTURED_GRID)
      call StructuredGridComputeVolumes(grid%structured_grid,option, &
                                        grid%nL2G,volume)
    case(UNSTRUCTURED_GRID)
  end select

end subroutine GridComputeVolumes

! ************************************************************************** !
!
! GridLocalizeRegions: Resticts regions to cells local to processor
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine GridLocalizeRegions(grid,region_list,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type) :: option
  
  type(region_type), pointer :: region
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, allocatable :: temp_int_array(:)
  PetscInt :: i, j, k, count, local_count, local_ghosted_id, local_id
  PetscInt :: i_min, i_max, j_min, j_max, k_min, k_max
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max
  PetscReal :: shift
  PetscErrorCode :: ierr
  
  region => region_list%first
  do
  
    if (.not.associated(region)) exit
    
    if (.not.associated(region%cell_ids)) then
      if (region%i1 > 0 .and. region%i2 > 0 .and. &
          region%j1 > 0 .and. region%j2 > 0 .and. &
          region%k1 > 0 .and. region%k2 > 0) then

        ! convert indexing from global (entire domain) to local processor
        region%i1 = region%i1 - grid%structured_grid%nxs
        region%i2 = region%i2 - grid%structured_grid%nxs
        region%j1 = region%j1 - grid%structured_grid%nys
        region%j2 = region%j2 - grid%structured_grid%nys
        region%k1 = region%k1 - grid%structured_grid%nzs
        region%k2 = region%k2 - grid%structured_grid%nzs
          
        ! clip region to within local processor domain
        region%i1 = max(region%i1,1)
        region%i2 = min(region%i2,grid%structured_grid%nlx)
        region%j1 = max(region%j1,1)
        region%j2 = min(region%j2,grid%structured_grid%nly)
        region%k1 = max(region%k1,1)
        region%k2 = min(region%k2,grid%structured_grid%nlz)
         
        count = 0  
        if (region%i1 <= region%i2 .and. &
            region%j1 <= region%j2 .and. &
            region%k1 <= region%k2) then
          region%num_cells = (region%i2-region%i1+1)* &
                             (region%j2-region%j1+1)* &
                             (region%k2-region%k1+1)
          allocate(region%cell_ids(region%num_cells))
          if (region%iface /= 0) then
            allocate(region%faces(region%num_cells))
            region%faces = region%iface
          endif
          region%cell_ids = 0
          do k=region%k1,region%k2
            do j=region%j1,region%j2
              do i=region%i1,region%i2
                count = count + 1
                region%cell_ids(count) = &
                       i + (j-1)*grid%structured_grid%nlx + &
                       (k-1)*grid%structured_grid%nlxy
              enddo
            enddo
          enddo
!          if (region%num_cells > 0) then
!            region%coordinates(1)%x = grid%x(region%cell_ids(ONE_INTEGER))
!            region%coordinates(1)%y = grid%y(region%cell_ids(ONE_INTEGER))
!            region%coordinates(1)%z = grid%z(region%cell_ids(ONE_INTEGER))
!          endif
        else
          region%num_cells = 0
        endif
     
        if (count /= region%num_cells) &
          call printErrMsg(option,"Mismatch in number of cells in block region")

      else if (associated(region%coordinates)) then
        if (size(region%coordinates) == ONE_INTEGER) then
          if (region%coordinates(ONE_INTEGER)%x >= grid%x_min_global .and. &
              region%coordinates(ONE_INTEGER)%x <= grid%x_max_global .and. &
              region%coordinates(ONE_INTEGER)%y >= grid%y_min_global .and. &
              region%coordinates(ONE_INTEGER)%y <= grid%y_max_global .and. &
              region%coordinates(ONE_INTEGER)%z >= grid%z_min_global .and. &
              region%coordinates(ONE_INTEGER)%z <= grid%z_max_global) then
            select case(grid%itype)
              case(STRUCTURED_GRID)
                call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                                    region%coordinates(ONE_INTEGER)%x, &
                                                    region%coordinates(ONE_INTEGER)%y, &
                                                    region%coordinates(ONE_INTEGER)%z, &
                                                    i,j,k)
                if (i > 0 .and. j > 0 .and. k > 0) then
                  region%num_cells = 1
                  allocate(region%cell_ids(region%num_cells))
                  if (region%iface /= 0) then
                    allocate(region%faces(region%num_cells))
                    region%faces = region%iface
                  endif
                  region%cell_ids = 0
                  region%cell_ids(1) = i + (j-1)*grid%structured_grid%nlx + &
                                      (k-1)*grid%structured_grid%nlxy
                else
                  region%num_cells = 0
                endif
                call MPI_Allreduce(region%num_cells,count,ONE_INTEGER,MPI_INTEGER,MPI_SUM, &
                                   PETSC_COMM_WORLD,ierr)   

! the next test as designed will only work on a uniform grid
                if(grid%structured_grid%p_samr_patch==0) then
                   if (count /= 1) then
                      write(string,*) 'Region: (coord)', &
                           region%coordinates(ONE_INTEGER)%x, &
                           region%coordinates(ONE_INTEGER)%y, &
                           region%coordinates(ONE_INTEGER)%z, &
                           ' not found in global domain.', count
                      call printErrMsg(option,string)
                   endif
                endif
            end select
          endif
        else if (size(region%coordinates) == 2) then
          x_min = min(region%coordinates(ONE_INTEGER)%x, &
                      region%coordinates(TWO_INTEGER)%x)
          x_max = max(region%coordinates(ONE_INTEGER)%x, &
                      region%coordinates(TWO_INTEGER)%x)
          y_min = min(region%coordinates(ONE_INTEGER)%y, &
                      region%coordinates(TWO_INTEGER)%y)
          y_max = max(region%coordinates(ONE_INTEGER)%y, &
                      region%coordinates(TWO_INTEGER)%y)
          z_min = min(region%coordinates(ONE_INTEGER)%z, &
                      region%coordinates(TWO_INTEGER)%z)
          z_max = max(region%coordinates(ONE_INTEGER)%z, &
                      region%coordinates(TWO_INTEGER)%z)
                      
          ! shift box slightly inward
          shift = 1.d-8*(x_max-x_min)
          x_min = x_min+shift            
          x_max = x_max-shift
          shift = 1.d-8*(y_max-y_min)
          y_min = y_min+shift            
          y_max = y_max-shift
          shift = 1.d-8*(z_max-z_min)
          z_min = z_min+shift            
          z_max = z_max-shift
                      
          ! ensure overlap
          if (x_min <= grid%x_max_local .and. &
              x_max >= grid%x_min_local .and. &
              y_min <= grid%y_max_local .and. &
              y_max >= grid%y_min_local .and. &
              z_min <= grid%z_max_local .and. &
              z_max >= grid%z_min_local) then
              
            ! get I,J,K bounds
            select case(grid%itype)
              case(STRUCTURED_GRID)
                ! local, non-ghosted i,j,k's are returned
                call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                                    max(x_min,grid%x_min_local), &
                                                    max(y_min,grid%y_min_local), &
                                                    max(z_min,grid%z_min_local), &
                                                    i_min,j_min,k_min)
                call StructGridGetIJKFromCoordinate(grid%structured_grid, &
                                                    min(x_max,grid%x_max_local), &
                                                    min(y_max,grid%y_max_local), &
                                                    min(z_max,grid%z_max_local), &
                                                    i_max,j_max,k_max)
                if (i_min > 0 .and. j_min > 0 .and. k_min > 0 .and. &
                    i_max > 0 .and. j_max > 0 .and. k_max > 0) then
                  region%num_cells = (i_max-i_min+1)*(j_max-j_min+1)*(k_max-k_min+1)
                  allocate(region%cell_ids(region%num_cells))
                  if (region%iface /= 0) then
                    allocate(region%faces(region%num_cells))
                    region%faces = region%iface
                  endif
                  region%cell_ids = 0
                  count = 0
                  do k = k_min, k_max
                    do j = j_min, j_max
                      do i = i_min, i_max
                        count = count+1
                        region%cell_ids(count) = i + (j-1)*grid%structured_grid%nlx + &
                                            (k-1)*grid%structured_grid%nlxy
                      enddo
                    enddo
                  enddo
                else
                  write(string,*) 'GridLocalizeRegions, between two points'
                  call printErrMsg(option,string)
                endif
            end select
          endif
        endif    
      endif 
    else
#if 0
      allocate(temp_int_array(region%num_cells))
      temp_int_array = 0
      if (grid%itype == STRUCTURED_GRID) then
        do count=1,region%num_cells
          i = mod(region%cell_ids(count),grid%structured_grid%nx) - &
                grid%structured_grid%nxs
          j = mod((region%cell_ids(count)-1)/grid%structured_grid%nx, &
                  grid%structured_grid%ny)+1 - &
                grid%structured_grid%nys
          k = ((region%cell_ids(count)-1)/grid%structured_grid%nxy)+1 - &
                grid%structured_grid%nzs
          if (i > 0 .and. i <= grid%structured_grid%nlx .and. &
              j > 0 .and. j <= grid%structured_grid%nly .and. &
              k > 0 .and. k <= grid%structured_grid%nlz) then
            temp_int_array(local_count) = &
                i + (j-1)*grid%structured_grid%nlx + &
                (k-1)*grid%structured_grid%nlxy
            local_count = local_count + 1
          endif
        enddo
      else
        call GridCreateNaturalToGhostedHash(grid,option)
        do count=1,region%num_cells
          local_ghosted_id = UnstructGridGetGhostIdFromHash( &
                                                    grid%unstructured_grid, &
                                                    region%cell_ids(count))
          if (local_ghosted_id > -1) then
            local_id = grid%nG2L(local_ghosted_id)
            if (local_id > -1) then
              temp_int_array(local_count) = local_id
              local_count = local_count + 1
            endif
          endif
        enddo
        call GridDestroyHashTable(grid)
      endif
      if (local_count /= region%num_cells) then
        deallocate(region%cell_ids)
        allocate(region%cell_ids(local_count))
        region%cell_ids(1:local_count) = temp_int_array(1:local_count)
      endif
      deallocate(temp_int_array)
#endif        
    endif
    
    if (region%num_cells == 0 .and. associated(region%cell_ids)) &
      deallocate(region%cell_ids)
    if (region%num_cells == 0 .and. associated(region%faces)) &
      deallocate(region%faces)
    region => region%next
    
  enddo

end subroutine GridLocalizeRegions

! ************************************************************************** !
!
! GridCopyIntegerArrayToPetscVec: Copies values from an integer array into a 
!                                 PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyIntegerArrayToPetscVec(array,vector,num_values)

  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyIntegerArrayToPetscVec

! ************************************************************************** !
!
! GridCopyRealArrayToPetscVec: Copies values from an integer array into a 
!                              PETSc Vec
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyRealArrayToPetscVec(array,vector,num_values)

  implicit none
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
    
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  vec_ptr(1:num_values) = array(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyRealArrayToPetscVec

! ************************************************************************** !
!
! GridCopyPetscVecToIntegerArray: Copies values from a PETSc Vec to an  
!                                 integer array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyPetscVecToIntegerArray(array,vector,num_values)

  implicit none

#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  
  PetscInt :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  do i=1,num_values
    array(i) = int(vec_ptr(i)+1.d-4)
  enddo
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyPetscVecToIntegerArray

! ************************************************************************** !
!
! GridCopyPetscVecToRealArray: Copies values from a PETSc Vec to an integer 
!                              array
! author: Glenn Hammond
! date: 12/18/07
!
! ************************************************************************** !
subroutine GridCopyPetscVecToRealArray(array,vector,num_values)

  implicit none
  
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
    
  PetscReal :: array(:)
  Vec :: vector
  PetscInt :: num_values
  
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(vector,vec_ptr,ierr)
  array(1:num_values) = vec_ptr(1:num_values)
  call VecRestoreArrayF90(vector,vec_ptr,ierr)
  
end subroutine GridCopyPetscVecToRealArray

! ************************************************************************** !
!
! GridCreateNaturalToGhostedHash: Creates a hash table for looking up the  
!                                 local ghosted id of a natural id, if it 
!                                 exists
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridCreateNaturalToGhostedHash(grid,option)

  use Option_module
  use Logging_module  
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_ghosted_id, natural_id
  PetscInt :: num_in_hash, num_ids_per_hash, hash_id, id, ierr
  PetscInt :: max_num_ids_per_hash = 0
  PetscInt, pointer :: hash(:,:,:), temp_hash(:,:,:)

  if (associated(grid%hash)) return

  call PetscLogEventBegin(logging%event_hash_create, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
                          
  ! initial guess of 10% of ids per hash
  ! must be at least 5 so that reallocation (*1.2) works below
  num_ids_per_hash = max(grid%nlmax/(grid%num_hash_bins/10),5)

  allocate(hash(2,0:num_ids_per_hash,grid%num_hash_bins))
  hash(:,:,:) = 0

  ! recall that natural ids are zero-based
  do local_ghosted_id = 1, grid%ngmax
    natural_id = grid%nG2A(local_ghosted_id)+1
    hash_id = mod(natural_id,grid%num_hash_bins)+1 
    num_in_hash = hash(1,0,hash_id)
    num_in_hash = num_in_hash+1
    if (num_in_hash > max_num_ids_per_hash) max_num_ids_per_hash = num_in_hash
    ! if a hash runs out of space reallocate
    if (num_in_hash > num_ids_per_hash) then 
      allocate(temp_hash(2,0:num_ids_per_hash,0:grid%num_hash_bins))
      ! copy old hash
      temp_hash(1:2,0:num_ids_per_hash,grid%num_hash_bins) = &
                             hash(1:2,0:num_ids_per_hash,grid%num_hash_bins)
      deallocate(hash)
      ! recompute hash 20% larger
      num_ids_per_hash = int(dble(num_ids_per_hash)*1.2)
      allocate(hash(1:2,0:num_ids_per_hash,grid%num_hash_bins))
      ! copy old to new
      do hash_id = 1, grid%num_hash_bins
        do id = 1, temp_hash(1,0,hash_id)
          hash(1:2,id,hash_id) = temp_hash(1:2,id,hash_id)
        enddo
        hash(1,0,hash_id) = temp_hash(1,0,hash_id)
      enddo
      deallocate(temp_hash)
    endif
    hash(1,0,hash_id) = num_in_hash
    hash(1,num_in_hash,hash_id) = natural_id
    hash(2,num_in_hash,hash_id) = local_ghosted_id
  enddo

  grid%hash => hash
  
!  call GridPrintHashTable(grid)
  call mpi_allreduce(max_num_ids_per_hash,num_in_hash,ONE_INTEGER,MPI_INTEGER, &
                     MPI_MAX,PETSC_COMM_WORLD,ierr)
  if (option%myrank == 0) print *, 'max_num_ids_per_hash:', num_in_hash

  call PetscLogEventEnd(logging%event_hash_create, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

end subroutine GridCreateNaturalToGhostedHash

! ************************************************************************** !
!
! GetLocalIdFromNaturalId: Returns the local id corresponding to a natural
!                          id or 0, if the natural id is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalIdFromNaturalId(grid,natural_id)

  implicit none

  type(grid_type) :: grid

  PetscInt :: natural_id, local_id
  
  do local_id = 1, grid%nlmax
    if (natural_id == grid%nL2A(local_id)+1) then
      GridGetLocalIdFromNaturalId = local_id
      return
    endif
  enddo
  GridGetLocalIdFromNaturalId = 0

end function GridGetLocalIdFromNaturalId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromNatId: Returns the local ghosted id corresponding 
!                                 to a natural id or 0, if the natural id 
!                                 is off-processor
! WARNING: Extremely inefficient for large jobs
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromNatId(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: local_ghosted_id
  
  do local_ghosted_id = 1, grid%ngmax
    if (natural_id == grid%nG2A(local_ghosted_id)+1) then
      GridGetLocalGhostedIdFromNatId = local_ghosted_id
      return 
    endif
  enddo
  GridGetLocalGhostedIdFromNatId = 0

end function GridGetLocalGhostedIdFromNatId

! ************************************************************************** !
!
! GridGetLocalGhostedIdFromHash: Returns the local ghosted id of a natural 
!                                id, if it exists.  Otherwise 0 is returned
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
PetscInt function GridGetLocalGhostedIdFromHash(grid,natural_id)

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  
  PetscInt :: hash_id, id

  GridGetLocalGhostedIdFromHash = 0
  hash_id = mod(natural_id,grid%num_hash_bins)+1 
  do id = 1, grid%hash(1,0,hash_id)
    if (grid%hash(1,id,hash_id) == natural_id) then
      GridGetLocalGhostedIdFromHash = grid%hash(2,id,hash_id)
      return
    endif
  enddo

end function GridGetLocalGhostedIdFromHash

! ************************************************************************** !
!
! GridDestroyHashTable: Deallocates the hash table
! author: Glenn Hammond
! date: 03/07/07
!
! ************************************************************************** !
subroutine GridDestroyHashTable(grid)

  implicit none

  type(grid_type), pointer :: grid
  
  if (associated(grid%hash)) deallocate(grid%hash)
  nullify(grid%hash)
  grid%num_hash_bins = 100

end subroutine GridDestroyHashTable

! ************************************************************************** !
!
! UnstructGridPrintHashTable: Prints the hashtable for viewing
! author: Glenn Hammond
! date: 03/09/07
!
! ************************************************************************** !
subroutine GridPrintHashTable(grid)

  implicit none

  type(grid_type) :: grid
  
  PetscInt :: ihash, id, fid

  fid = 87 
  open(fid,file='hashtable.dat',action='write')
  do ihash=1,grid%num_hash_bins
    write(fid,'(a4,i3,a,i5,a2,x,200(i6,x))') 'Hash',ihash,'(', &
                         grid%hash(1,0,ihash), &
                         '):', &
                         (grid%hash(1,id,ihash),id=1,grid%hash(1,0,ihash))
  enddo
  close(fid)

end subroutine GridPrintHashTable

! ************************************************************************** !
!
! GridDestroy: Deallocates a grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine GridDestroy(grid)

  implicit none
  
  type(grid_type), pointer :: grid
    
  if (.not.associated(grid)) return
      
  if (associated(grid%nL2G)) deallocate(grid%nL2G)
  nullify(grid%nL2G)
  if (associated(grid%nG2L)) deallocate(grid%nG2L)
  nullify(grid%nG2L)
  if (associated(grid%nL2A)) deallocate(grid%nL2A)
  nullify(grid%nL2A)
  if (associated(grid%nG2A)) deallocate(grid%nG2A)
  nullify(grid%nG2A)

  if (associated(grid%x)) deallocate(grid%x)
  nullify(grid%x)
  if (associated(grid%y)) deallocate(grid%y)
  nullify(grid%y)
  if (associated(grid%z)) deallocate(grid%z)
  nullify(grid%z)
  
  if (associated(grid%hash)) call GridDestroyHashTable(grid)
  
  call UnstructuredGridDestroy(grid%unstructured_grid)    
  call StructuredGridDestroy(grid%structured_grid)
                                           
  call ConnectionDestroyList(grid%internal_connection_set_list)

end subroutine GridDestroy

subroutine GridVecGetArrayF90(grid, vec, f90ptr, ierr)

  implicit none

  type(grid_type) :: grid
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  integer :: ierr

  if (.not.associated(grid%structured_grid)) then
     call VecGetArrayF90(vec, f90ptr, ierr)
  else
     call StructuredGridVecGetArrayF90(grid%structured_grid, vec, f90ptr, ierr)
  endif

end subroutine GridVecGetArrayF90

subroutine GridVecRestoreArrayF90(grid, vec, f90ptr, ierr)

  implicit none

  type(grid_type) :: grid
  Vec:: vec
  PetscReal, pointer :: f90ptr(:)
  integer :: ierr

  if (.not.associated(grid%structured_grid)) then
     call VecRestoreArrayF90(vec, f90ptr, ierr)
  else
     call StructGridVecRestoreArrayF90(grid%structured_grid, vec, f90ptr, ierr)
  endif

end subroutine GridVecRestoreArrayF90

end module Grid_module
