#ifdef GEOMECH

module Geomechanics_Realization_module

  use Geomechanics_Discretization_module
  
  implicit none
  
private


#include "definitions.h"

  type, public :: geomech_realization_type

    PetscInt :: id
    type(geomech_discretization_type), pointer :: discretization

  end type geomech_realization_type
  


contains

! ************************************************************************** !
!
! GeomechRealizCreateDiscretization: This subroutine creates a grid for 
! geomechanics
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
!subroutine GeomechRealizCreateDiscretization(geomechanics_realization)




!end subroutine GeomechRealizCreateDiscretization

end module Geomechanics_Realization_module
#endif
! GEOMECH
