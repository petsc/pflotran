module Solution_module

  use Grid_module
  use Option_module
  use Region_module

  implicit none
  
private

  type, public :: solution_type

    type(grid_type), pointer :: grid
    type(option_type), pointer :: option
    type(region_list_type), pointer :: regions

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
  
  type(solution_type), pointer :: solution
  
  allocate(solution)
  solution%option => createOption()
  nullify(solution%grid)
  allocate(solution%regions)
  call initRegionList(solution%regions)
  
  createSolution => solution
  
end function createSolution  
  
end module Solution_module
