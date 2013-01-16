module Realization_Base_class

  implicit none

  private

#include "definitions.h"
  type, public :: realization_base_type

    PetscInt :: id
    
  end type realization_base_type
  
  public :: RealizationBaseInit

contains

! ************************************************************************** !
!
! RealizationBaseInit: Initializes variables/objects in base realization class
! author: Glenn Hammond
! date: 01/16/13
!
! ************************************************************************** !
subroutine RealizationBaseInit(realization)

  implicit none
  
  class(realization_base_type) :: realization
  
  realization%id = 0

end subroutine RealizationBaseInit
  
end module Realization_Base_class
