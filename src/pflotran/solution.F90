module Solution_module

  use Grid_module
  use Option_module

  implicit none
  
private

  type, public :: solution_type

    type(grid_type), pointer :: grid
    type(option_type), pointer :: option

  end type solution_type

  public :: createSolution
  
contains
  
! ************************************************************************** !
!
! createSolution: Allocates and initializes a new Solution object
! author: Glenn Hammond
! date: 10/25/07
!
! ************************************************************************** !
function createSolution()

  implicit none
  
  type(solution_type), pointer :: createSolution
  
  allocate(createSolution)
  createSolution%option => createOption()
  nullify(createSolution%grid)
  
end function createSolution  
  
end module Solution_module
