module Option_Flow_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"


  type, public :: flow_option_type 
  
    PetscReal :: inf_rel_update_tol
    PetscReal :: inf_scaled_res_tol   
    PetscBool :: check_post_convergence
  
  end type flow_option_type
  
  public :: OptionFlowCreate, &
            OptionFlowInitAll, &
            OptionFlowInitRealization, &
            OptionFlowDestroy

contains

! ************************************************************************** !

function OptionFlowCreate()
  ! 
  ! Allocates and initializes a new Option object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type), pointer :: OptionFlowCreate
  
  type(flow_option_type), pointer :: option
  
  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide 
  ! whether the member needs initialization once for all stochastic 
  ! simulations or initialization for every realization (e.g. within multiple 
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionFlowInitAll(option)
  OptionFlowCreate => option
  
end function OptionFlowCreate

! ************************************************************************** !

subroutine OptionFlowInitAll(option)
  ! 
  ! Initializes all option variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type) :: option
  
  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)
  
  call OptionFlowInitRealization(option)

end subroutine OptionFlowInitAll

! ************************************************************************** !

subroutine OptionFlowInitRealization(option)
  ! 
  ! Initializes option variables specific to a single
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type) :: option
  
  ! These variables should be initialized once at the beginning of every 
  ! PFLOTRAN realization or simulation of a single realization
    
  option%check_post_convergence = PETSC_FALSE
  option%inf_rel_update_tol = UNINITIALIZED_DOUBLE
  option%inf_scaled_res_tol = UNINITIALIZED_DOUBLE 
  
end subroutine OptionFlowInitRealization

! ************************************************************************** !

subroutine OptionFlowDestroy(option)
  ! 
  ! Deallocates an option
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  type(flow_option_type), pointer :: option
  
  if (.not.associated(option)) return
  ! all kinds of stuff needs to be added here.

  ! all the below should be placed somewhere other than option.F90
  deallocate(option)
  nullify(option)
  
end subroutine OptionFlowDestroy

end module Option_Flow_module
