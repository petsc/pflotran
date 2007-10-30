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
            computeGridCellVolumes
  
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
                                  grid%nL2A,grid%nL2A,grid%nG2N)
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
    call computeStructuredCellVolumes(grid%structured_grid,option,grid%nL2G)
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
  
end module Grid_module
