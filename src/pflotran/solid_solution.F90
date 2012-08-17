module Solid_Solution_module

  use Mineral_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Solid_Solution_Aux_module
  
  implicit none
  
  private 

#include "definitions.h"

  public :: SolidSolutionRead
            
contains

! ************************************************************************** !
!
! SolidSolutionRead: Reads solid solution reactions
! author: Glenn Hammond
! date: 08/16/12
!
! ************************************************************************** !
subroutine SolidSolutionRead(reaction,input,option)

  use Option_module
  use String_module
  use Input_module
  use Utility_module
  
  implicit none
  
  type(reaction_type) :: reaction
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name  
  character(len=MAXWORDLENGTH) :: card
  type(solid_solution_type), pointer :: solid_solution, prev_solid_solution
           
  nullify(prev_solid_solution)
  do
    call InputReadFlotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
          
    solid_solution => SolidSolutionCreate()
    call InputReadWord(input,option,solid_solution%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SOLID_SOLUTIONS')   
    if (.not.associated(reaction%solid_solution_list)) then
      reaction%solid_solution_list => solid_solution
    endif
    if (associated(prev_solid_solution)) then
      prev_solid_solution%next => solid_solution
    endif
    prev_solid_solution => solid_solution
    nullify(solid_solution)
  enddo

end subroutine SolidSolutionRead

end module Solid_Solution_module
