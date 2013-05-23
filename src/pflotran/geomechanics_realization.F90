#ifdef GEOMECH

module Geomechanics_Realization_module

  use Geomechanics_Discretization_module
  use Input_module
  use Patch_module
  use Option_module
  use Level_module
  
  implicit none
  
private


#include "definitions.h"

  type, public :: geomech_realization_type

    PetscInt :: id
    type(geomech_discretization_type), pointer :: discretization
    type(input_type), pointer :: input
    type(patch_type), pointer :: patch
    type(option_type), pointer :: option
    type(level_list_type), pointer :: level_list

  end type geomech_realization_type

public :: GeomechRealizCreate, &
          GeomechRealizDestroy

contains

! ************************************************************************** !
!
! GeomechRealizCreate: This subroutine creates realization for geomechanics
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
function GeomechRealizCreate(option)

  implicit none

  type(geomech_realization_type), pointer :: GeomechRealizCreate
  type(geomech_realization_type), pointer :: geomech_realization
  type(option_type), pointer              :: option
  
  geomech_realization%id = 0
  if (associated(option)) then
    geomech_realization%option => option
  else
    geomech_realization%option => OptionCreate()
  endif
  
  nullify(geomech_realization%input)
  geomech_realization%discretization => GeomechDiscretizationCreate()
  geomech_realization%level_list     => LevelCreateList()
  
  nullify(geomech_realization%patch)

  GeomechRealizCreate => geomech_realization
  
end function GeomechRealizCreate

! ************************************************************************** !
!
! GeomechRealizDestroy: This subroutine deallocates geomechanics realization
! author: Satish Karra, LANL
! date: 05/23/13
!
! ************************************************************************** !
subroutine GeomechRealizDestroy(geomech_realization)

  implicit none
  
  type(geomech_realization_type), pointer :: geomech_realization
  
  if(.not.associated(geomech_realization)) return
  
  call LevelDestroyList(geomech_realization%level_list)
  call GeomechDiscretizationDestroy(geomech_realization%discretization)
  

end subroutine GeomechRealizDestroy



end module Geomechanics_Realization_module
#endif
! GEOMECH
