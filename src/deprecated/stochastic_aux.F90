module Stochastic_Aux_module

#include "finclude/petscsys.h"
  use petscsys
  use Simulation_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
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

function StochasticCreate()
  ! 
  ! Create stochastic simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 

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

subroutine StochasticRead(stochastic,input,option)
  ! 
  ! Reads debugging data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 

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

subroutine StochasticDestroy(stochastic)
  ! 
  ! Deallocates a stochastic object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/09
  ! 

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
