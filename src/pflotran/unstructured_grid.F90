module Unstructured_Grid_module

  use Connection_module
  
  implicit none

  private 

#include "definitions.h"

  type, public :: unstructured_grid_type
    integer :: num_cells
  end type
  
  public :: initUnstructuredGrid, &
            createUnstructuredDMs, &
            computeUnstructInternalConnect, &
            computeUnstructBoundaryConnect


contains

! ************************************************************************** !
!
! initUnstructuredGrid: Initializes an unstructured grid object
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine initUnstructuredGrid(grid)

  implicit none
  
  type(unstructured_grid_type) :: grid

end subroutine initUnstructuredGrid

! ************************************************************************** !
!
! createStructuredDMs: creates unstructured distributed, parallel meshes/grids
! author: Glenn Hammond
! date: 10/22/07
!
! ************************************************************************** !
subroutine createUnstructuredDMs(solution,unstructured_grid)
      
  use Solution_module
      
  implicit none
  
  type(solution_type) :: solution
  type(unstructured_grid_type) :: unstructured_grid
  
end subroutine createUnstructuredDMs
  
! ************************************************************************** !
!
! computeUnstructInternalConnect: computes internal connectivity of an  
!                                 unstructured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function computeUnstructInternalConnect(solution,unstructured_grid)

  use Connection_module
  use Solution_module
  
  implicit none
  
  type(connection_type), pointer :: computeUnstructInternalConnect
  type(solution_type) :: solution
  type(unstructured_grid_type) :: unstructured_grid

  nullify(computeUnstructInternalConnect)
  
end function computeUnstructInternalConnect

! ************************************************************************** !
!
! computeUnstructBoundaryConnect: computes boundary connectivity of an 
!                                 unstructured grid
! author: Glenn Hammond
! date: 10/15/07
!
! ************************************************************************** !
function computeUnstructBoundaryConnect(solution,unstructured_grid)

  use Connection_module
  use Solution_module
  
  implicit none

  type(connection_type), pointer :: computeUnstructBoundaryConnect  
  type(solution_type) :: solution
  type(unstructured_grid_type) :: unstructured_grid
  
  nullify(computeUnstructBoundaryConnect)
  
end function computeUnstructBoundaryConnect

end module Unstructured_Grid_module
