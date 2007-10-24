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
subroutine createUnstructuredDMs(unstructured_grid,option)
      
  use Option_module
      
  implicit none
  
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
end subroutine createUnstructuredDMs
  
! ************************************************************************** !
!
! computeUnstructInternalConnect: computes internal connectivity of an  
!                                 unstructured grid
! author: Glenn Hammond
! date: 10/17/07
!
! ************************************************************************** !
function computeUnstructInternalConnect(unstructured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none
  
  type(connection_type), pointer :: computeUnstructInternalConnect
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option

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
function computeUnstructBoundaryConnect(unstructured_grid,option)

  use Connection_module
  use Option_module
  
  implicit none

  type(connection_type), pointer :: computeUnstructBoundaryConnect  
  type(unstructured_grid_type) :: unstructured_grid
  type(option_type) :: option
  
  nullify(computeUnstructBoundaryConnect)
  
end function computeUnstructBoundaryConnect

end module Unstructured_Grid_module
