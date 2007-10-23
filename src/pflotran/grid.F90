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
  
    !nL2G :  not collective, local processor: local  =>  ghosted local  
    !nG2L :  not collective, local processor:  ghosted local => local  
    !nG2N :  collective,  ghosted local => global index , used for   
    !                     matsetvaluesblocked ( not matsetvaluesblockedlocal)  
    !nL2A :   collective, local => natural index, used for initialization   
    !                              and source/sink setup  
    integer, pointer :: nL2G(:), nG2L(:), nL2A(:),nG2N(:)
    integer, pointer :: nG2A(:)
    
    integer, pointer :: ibconn(:)
    
    real*8, pointer :: x(:), y(:), z(:), delz(:) 
    real*8 :: x_max, x_min, y_max, y_min, z_max, z_min
    
    integer :: igeom
    
    Vec :: volume  ! Volume of a cell in the grid
    
    type(structured_grid_type), pointer :: structured_grid
    type(unstructured_grid_type), pointer :: unstructured_grid
    
    type(connection_list_type), pointer :: internal_connection_list, &
                                           boundary_connection_list

  end type grid_type


  public :: createGrid, &
            createStructuredDMs, &
            computeInternalConnectivity, &
            computeBoundaryConnectivity
  
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
! setNPXYZ: Initializes the decomposition of a structrued grid
! author: Glenn Hammond
! date: 10/23/07
!
! ************************************************************************** !
subroutine setNPXYZ(grid,npx,npy,npz)

  implicit none
  
  type(grid_type) :: grid
  integer :: npx, npy, npz
  
  if (grid%is_structured) then
    call setStructuredNPXYZ(grid%structured_grid,npx,npy,npz)
  else
    print *, 'ERROR: NPX,NPY,NPZ not supported for unstructured grid'
    stop
  endif

end subroutine setNPXYZ

! ************************************************************************** !
!
! createDMs: creates distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine createDMs(solution,grid)
      
  use Solution_module    
      
  implicit none
  
  type(solution_type) :: solution
  type(grid_type) :: grid
  
  if (grid%is_structured) then
    call createStructuredDMs(solution,grid%structured_grid)
  else 
    call createUnstructuredDMs(solution,grid%unstructured_grid)
  endif

end subroutine createDMs

! ************************************************************************** !
!
! computeInternalConnectivity: computes internal connectivity of a grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
subroutine computeInternalConnectivity(solution,grid)

  use Connection_module
  use Solution_module    
    
  implicit none
  
  type(solution_type) :: solution
  type(grid_type) :: grid
  
  type(connection_type), pointer :: connection
  
  if (grid%is_structured) then
    connection => &
      computeStructInternalConnect(solution,grid%structured_grid)
  else 
    connection => &
      computeUnstructInternalConnect(solution,grid%unstructured_grid)
  endif
  
  call addConnectionToList(connection,grid%internal_connection_list)
  
end subroutine computeInternalConnectivity

! ************************************************************************** !
!
! computeBoundaryConnectivity: computes boundary connectivity of a grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
subroutine computeBoundaryConnectivity(solution,grid)

  use Connection_module
  use Solution_module    

  implicit none
  
  type(solution_type) :: solution
  type(grid_type) :: grid
  
  type(connection_type), pointer :: connection
  
  if (grid%is_structured) then
    connection => &
      computeStructBoundaryConnect(solution,grid%structured_grid,grid%ibconn, &
                                   grid%nL2G)
  else 
    connection => &
      computeUnstructBoundaryConnect(solution,grid%unstructured_grid)
  endif
  
  call addConnectionToList(connection,grid%boundary_connection_list)

end subroutine computeBoundaryConnectivity

end module Grid_module
