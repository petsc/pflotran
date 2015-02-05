module Argonne_Mixed_Potential_module
  
  implicit none
  
  private

#include "finclude/petscsys.h"

  public :: AMPInit, &
            AMPRun, &
            AMPDestroy
  
contains

! ************************************************************************** !

subroutine AMPInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  ! 

  implicit none
  
end subroutine AMPInit

! ************************************************************************** !

subroutine AMPRun(time,conc)
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  !

  implicit none
  
  PetscReal :: time
  PetscReal :: conc(:)
  
end subroutine AMPRun

! ************************************************************************** !

subroutine AMPDestroy()
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/15
  ! 

  implicit none
  
end subroutine AMPDestroy

end module Argonne_Mixed_Potential_module
