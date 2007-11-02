module Grid_module

  use Structured_Grid_module
  use Unstructured_Grid_module
  use Connection_module
 
  implicit none

  private
 
#include "definitions.h"

! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "petscreldefs.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"

  type, public :: grid_type
  
    logical :: is_structured
    
    integer :: nmax   ! Total number of nodes in global domain
    integer :: nlmax  ! Total number of non-ghosted nodes in local domain.
    integer :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
#if 1    
    !nL2G :  not collective, local processor: local  =>  ghosted local  
    !nG2L :  not collective, local processor:  ghosted local => local  
    !nG2N :  collective,  ghosted local => global index , used for   
    !                     matsetvaluesblocked ( not matsetvaluesblockedlocal)  
    !nL2A :   collective, local => natural index, used for initialization   
    !                              and source/sink setup  
    integer, pointer :: nL2G(:), nG2L(:), nL2A(:),nG2N(:)
    integer, pointer :: nG2A(:)
#endif
    
    integer, pointer :: ibconn(:)
    
#if 1  
    real*8, pointer :: x(:), y(:), z(:), delz(:) 
#endif    
    
    integer :: igeom
    type(structured_grid_type), pointer :: structured_grid
    type(unstructured_grid_type), pointer :: unstructured_grid
    
    type(connection_list_type), pointer :: internal_connection_list, &
                                           boundary_connection_list

  end type grid_type


  public :: createGrid, &
            destroyGrid, &
            computeInternalConnectivity, &
            computeBoundaryConnectivity, &
            createPetscVector, &
            createJacobian, &
            createColoring, &
            DMGlobalToLocal, &
            DMGlobalToNatural, &
            mapGridIndices, &
            createDMs, &
            computeGridSpacing, &
            computeGridCoordinates, &
            computeGridCellVolumes, &
            computeBoundaryConnectivity2, &
            localizeRegions
  
contains

! ************************************************************************** !
!
! createGrid: Creates a structured or unstructured grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
function createGrid(igeom_)

  implicit none
  
  integer :: igeom_
  
  type(grid_type), pointer :: createGrid
  
  allocate(createGrid)
  call initGrid(createGrid)
  createGrid%igeom = igeom_
  
  if (igeom_ > 0) then
    createGrid%is_structured = .true.
    allocate(createGrid%structured_grid)
    call initStructuredGrid(createGrid%structured_grid)
    createGrid%structured_grid%igeom = igeom_
  else
    createGrid%is_structured = .false.
    allocate(createGrid%unstructured_grid)
    call initUnstructuredGrid(createGrid%unstructured_grid)
  endif

end function createGrid

! ************************************************************************** !
!
! initGrid: Initializes the abstract grid object
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine initGrid(grid)

  implicit none
  
  type(grid_type) :: grid
  
  grid%igeom = 0
  nullify(grid%structured_grid)
  nullify(grid%unstructured_grid)

end subroutine initGrid

! ************************************************************************** !
!
! createDMs: creates distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine createDMs(grid,option)
      
  use Option_module    
      
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  if (grid%is_structured) then
    call createStructuredDMs(grid%structured_grid,option)
    grid%nlmax = grid%structured_grid%nlmax
    grid%ngmax = grid%structured_grid%ngmax
  else 
    call createUnstructuredDMs(grid%unstructured_grid,option)
  endif

  ! allocate coordinate arrays  
  allocate(grid%x(grid%ngmax))
  allocate(grid%y(grid%ngmax))
  allocate(grid%z(grid%ngmax))

end subroutine createDMs

! ************************************************************************** !
!
! createPetscVector: Creates a global PETSc vector
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createPetscVector(grid,dm_index,vector,vector_type)

  implicit none
  
  type(grid_type) :: grid
  integer :: dm_index
  Vec :: vector
  integer :: vector_type
  
  if (grid%is_structured) then
    call createPetscVectorFromDA(grid%structured_grid,dm_index,vector, &
                                 vector_type)
  else
  endif
  
end subroutine createPetscVector
! ************************************************************************** !
!
! computeInternalConnectivity: computes internal connectivity of a grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
subroutine computeInternalConnectivity(grid,option)

  use Connection_module
  use Option_module    
    
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  type(connection_type), pointer :: connection
  
  if (grid%is_structured) then
    connection => &
      computeStructInternalConnect(grid%structured_grid,option)
  else 
    connection => &
      computeUnstructInternalConnect(grid%unstructured_grid,option)
  endif
  
  allocate(grid%internal_connection_list)
  call initConnectionList(grid%internal_connection_list)
  call addConnectionToList(connection,grid%internal_connection_list)
  
end subroutine computeInternalConnectivity

! ************************************************************************** !
!
! computeBoundaryConnectivity: computes boundary connectivity of a grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine computeBoundaryConnectivity(grid,option)

  use Connection_module
  use Option_module    

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  type(connection_type), pointer :: connection
  
  if (grid%is_structured) then
    connection => &
      computeStructBoundaryConnect(grid%structured_grid,option,grid%ibconn, &
                                   grid%nL2G)
  else 
    connection => &
      computeUnstructBoundaryConnect(grid%unstructured_grid,option)
  endif

  allocate(grid%boundary_connection_list)
  call initConnectionList(grid%boundary_connection_list)  
  call addConnectionToList(connection,grid%boundary_connection_list)

end subroutine computeBoundaryConnectivity

! ************************************************************************** !
!
! computeBoundaryConnectivity2: computes boundary connectivity of a grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine computeBoundaryConnectivity2(grid,option,boundary_condition_list)

  use Connection_module
  use Option_module  
  use Coupler_module  

  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_list_type) :: boundary_condition_list
  
  type(connection_type), pointer :: connection
  
  if (grid%is_structured) then
    connection => &
      computeStructBoundaryConnect2(grid%structured_grid,option,grid%ibconn, &
                                   grid%nL2G,boundary_condition_list)
  else 
    connection => &
      computeUnstructBoundaryConnect(grid%unstructured_grid,option)
  endif

  allocate(grid%boundary_connection_list)
  call initConnectionList(grid%boundary_connection_list)  
  call addConnectionToList(connection,grid%boundary_connection_list)

end subroutine computeBoundaryConnectivity2

! ************************************************************************** !
!
! mapGridIndices: maps global, local and natural indices of cells 
!                 to each other
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine mapGridIndices(grid)

  implicit none
  
  type(grid_type) :: grid
  
  if (grid%is_structured) then
    call mapStructuredGridIndices(grid%structured_grid,grid%nG2L,grid%nL2G, &
                                  grid%nL2A,grid%nG2A,grid%nG2N)
  else
  endif

end subroutine mapGridIndices

! ************************************************************************** !
!
! computeGridSpacing: Computes grid spacing (only for structured grid
! author: Glenn Hammond
! date: 10/26/07
!
! ************************************************************************** !
subroutine computeGridSpacing(grid)

  implicit none
  
  type(grid_type) :: grid
  
  if (grid%is_structured) then
    call computeStructuredGridSpacing(grid%structured_grid,grid%nL2A)
  endif
  
end subroutine computeGridSpacing

! ************************************************************************** !
!
! computeGridCoordinates: Computes x,y,z coordinates of grid cells
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine computeGridCoordinates(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  if (grid%is_structured) then
    call computeStructuredGridCoordinates(grid%structured_grid,option, &
                                          grid%x,grid%y,grid%z)
  else
  endif

end subroutine computeGridCoordinates

! ************************************************************************** !
!
! computeGridCellVolumes: Computes the volumes of cells in structured grid
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
subroutine computeGridCellVolumes(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  if (grid%is_structured) then
    call computeStructuredCellVolumes(grid%structured_grid,option, &
                                      grid%nL2G)
  else
  endif

end subroutine computeGridCellVolumes

! ************************************************************************** !
!
! createJacobian: Creates Jacobian matrix associated with grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createJacobian(grid,option)

  use Option_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  if (grid%is_structured) then
    call createStructuredGridJacobian(grid%structured_grid,option)
  else
  endif

end subroutine createJacobian

! ************************************************************************** !
!
! createColoring: Creates ISColoring for grid
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine createColoring(grid,option,coloring)

  use Option_module
  
  implicit none

#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
  
  type(grid_type) :: grid
  type(option_type) :: option
  ISColoring :: coloring
  
  if (grid%is_structured) then
    call createStructuredGridColoring(grid%structured_grid,option,coloring)
  else
  endif
  
end subroutine createColoring

! ************************************************************************** !
!
! DMGlobalToLocal: Performs global to local communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DMGlobalToLocal(grid,global_vec,local_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: local_vec
  integer :: dm_index
  
  if (grid%is_structured) then
    call DMStructGlobalToLocal(grid%structured_grid,global_vec,local_vec, &
                               dm_index)
  else
  endif
  
end subroutine DMGlobalToLocal
  
! ************************************************************************** !
!
! DMGlobalToNatural: Performs global to natural communication with DM
! author: Glenn Hammond
! date: 10/24/07
!
! ************************************************************************** !
subroutine DMGlobalToNatural(grid,global_vec,natural_vec,dm_index)

  implicit none
  
  type(grid_type) :: grid
  Vec :: global_vec
  Vec :: natural_vec
  integer :: dm_index
  
  if (grid%is_structured) then
    call DMStructGlobalToNatural(grid%structured_grid,global_vec,natural_vec, &
                               dm_index)
  else
  endif
  
end subroutine DMGlobalToNatural

! ************************************************************************** !
!
! localizeRegions: Resticts regions to cells local to processor
! author: Glenn Hammond
! date: 10/29/07
!
! ************************************************************************** !
subroutine localizeRegions(region_list,grid,option)

  use Option_module
  use Region_module

  implicit none
  
  type(region_list_type), pointer :: region_list
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  type(region_type), pointer :: region
  integer, allocatable :: temp_int_array(:)
  integer :: i, j, k, count, local_count, local_ghosted_id, local_id
  
  region => region_list%first
  do
  
    if (.not.associated(region)) exit
    
    if (.not.associated(region%cell_ids) .and. &
        region%i1 > 0 .and. region%i2 > 0 .and. &
        region%j1 > 0 .and. region%j2 > 0 .and. &
        region%k1 > 0 .and. region%k2 > 0) then
        
      ! clip region to within local processor domain
      region%i1 = max(region%i1,grid%structured_grid%nxs+1)
      region%i2 = min(region%i2,grid%structured_grid%nxe)
      region%j1 = max(region%j1,grid%structured_grid%nys+1)
      region%j2 = min(region%j2,grid%structured_grid%nye)
      region%k1 = max(region%k1,grid%structured_grid%nzs+1)
      region%k2 = min(region%k2,grid%structured_grid%nze)
        
      region%num_cells = (region%i2-region%i1+1)* &
                         (region%j2-region%j1+1)* &
                         (region%k2-region%k1+1)
                         
      allocate(region%cell_ids(region%num_cells))
      region%cell_ids = 0
        
      count = 0  
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
      if (count /= region%num_cells) &
        call printErrMsg(option,"Mismatch in number of cells in block region")
    else
      allocate(temp_int_array(region%num_cells))
      temp_int_array = 0
      if (grid%is_structured) then
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
        do count=1,region%num_cells
          local_ghosted_id = GetLocalGhostedIdFromHash(grid%unstructured_grid, &
                                                       region%cell_ids(count))
          if (local_ghosted_id > -1) then
            local_id = grid%nG2L(local_ghosted_id)
            if (local_id > -1) then
              temp_int_array(local_count) = local_id
              local_count = local_count + 1
            endif
          endif
        enddo
      endif
      if (local_count /= region%num_cells) then
        deallocate(region%cell_ids)
        allocate(region%cell_ids(local_count))
        region%cell_ids(1:local_count) = temp_int_array(1:local_count)
      endif
      deallocate(temp_int_array)
    endif
    
    if (region%num_cells == 0) deallocate(region%cell_ids)
    region => region%next
    
  enddo

end subroutine localizeRegions

! ************************************************************************** !
!
! destroyGrid: Deallocates a grid
! author: Glenn Hammond
! date: 11/01/07
!
! ************************************************************************** !
subroutine destroyGrid(grid)

  implicit none
  
  type(grid_type), pointer :: grid
    
  if (.not.associated(grid)) return
      
  if (associated(grid%nL2G)) deallocate(grid%nL2G)
  nullify(grid%nL2G)
  if (associated(grid%nG2L)) deallocate(grid%nG2L)
  nullify(grid%nG2L)
  if (associated(grid%nL2A)) deallocate(grid%nL2A)
  nullify(grid%nL2A)
  if (associated(grid%nG2N)) deallocate(grid%nG2N)
  nullify(grid%nG2N)
  if (associated(grid%nG2A)) deallocate(grid%nG2A)
  nullify(grid%nG2A)

  if (associated(grid%ibconn)) deallocate(grid%ibconn)
  nullify(grid%ibconn)

  if (associated(grid%x)) deallocate(grid%x)
  nullify(grid%x)
  if (associated(grid%y)) deallocate(grid%y)
  nullify(grid%y)
  if (associated(grid%z)) deallocate(grid%z)
  nullify(grid%z)
  if (associated(grid%delz)) deallocate(grid%delz)
  nullify(grid%delz)
  
  call destroyUnstructuredGrid(grid%unstructured_grid)    
  call destroyStructuredGrid(grid%structured_grid)
                                           
  call destroyConnectionList(grid%internal_connection_list)
  call destroyConnectionList(grid%boundary_connection_list)

end subroutine destroyGrid
  
end module Grid_module
