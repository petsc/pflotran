#ifdef GEOMECH

module Geomechanics_Realization_module

  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Input_module
  use Option_module
  use Level_module
  use Geomechanics_Material_module
  
  implicit none
  
private


#include "definitions.h"

  type, public :: geomech_realization_type

    PetscInt :: id
    type(geomech_discretization_type), pointer :: discretization
    type(input_type), pointer :: input
    type(geomech_patch_type), pointer :: geomech_patch
    type(option_type), pointer :: option
    type(geomech_material_property_type), &
                           pointer :: geomech_material_properties
    type(geomech_material_property_ptr_type), &
                           pointer :: geomech_material_property_array(:)


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
  
  nullify(geomech_realization%geomech_material_properties)
  nullify(geomech_realization%geomech_material_property_array)
  
  nullify(geomech_realization%geomech_patch)

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
  
  call GeomechDiscretizationDestroy(geomech_realization%discretization)
  

end subroutine GeomechRealizDestroy



end module Geomechanics_Realization_module
#endif
! GEOMECH
