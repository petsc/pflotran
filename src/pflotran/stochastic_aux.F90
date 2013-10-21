module Stochastic_Aux_module

  use Simulation_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "finclude/petscsys.h"

  type, public :: stochastic_type
    PetscInt :: num_groups
    PetscInt :: num_realizations
    PetscInt :: num_local_realizations
    PetscInt, pointer :: realization_ids(:)
    type(simulation_type), pointer :: simulation
  end type stochastic_type
  
  public :: StochasticCreate, &
            StochasticRead, &
            StochasticDestroy
  
contains

! ************************************************************************** !
!
! StochasticCreate: Create stochastic simulation object
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
function StochasticCreate()

  implicit none
  
  type(stochastic_type), pointer :: StochasticCreate
  
  type(stochastic_type), pointer :: stochastic
  
  allocate(stochastic)
  stochastic%num_realizations = 0
  stochastic%num_groups = 0
  nullify(stochastic%simulation)
  StochasticCreate => stochastic

end function StochasticCreate

! ************************************************************************** !
!
! StochasticRead: Reads debugging data from the input file
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticRead(stochastic,input,option)

  use Option_module
  use Input_Aux_module
  
  implicit none
    
  type(stochastic_type) :: stochastic
  type(input_type) :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','STOCHASTIC')   
      
    select case(trim(keyword))
    
      case('NUM_REALIZATIONS')
        call InputReadInt(input,option,stochastic%num_realizations)
        call InputErrorMsg(input,option,'num_realizations','STOCHASTIC')
      case('NUM_GROUPS')
        call InputReadInt(input,option,stochastic%num_groups)
        call InputErrorMsg(input,option,'num_groups','STOCHASTIC')

    end select 
  
  enddo  

end subroutine StochasticRead

! ************************************************************************** !
!
! StochasticDestroy: Deallocates a stochastic object
! author: Glenn Hammond
! date: 02/04/09
!
! ************************************************************************** !
subroutine StochasticDestroy(stochastic)

  implicit none
  
  type(stochastic_type), pointer :: stochastic
  
  if (.not.associated(stochastic)) return
  
  if (associated(stochastic%simulation)) then
    call SimulationDestroy(stochastic%simulation)
    nullify(stochastic%simulation)
  endif
  
  if (associated(stochastic%realization_ids)) &
    deallocate(stochastic%realization_ids)
  nullify(stochastic%realization_ids)
  
  deallocate(stochastic)
  nullify(stochastic)

end subroutine StochasticDestroy

end module Stochastic_Aux_module
