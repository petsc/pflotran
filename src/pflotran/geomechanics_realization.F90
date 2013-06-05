#ifdef GEOMECH

module Geomechanics_Realization_module

  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Input_module
  use Option_module
  use Geomechanics_Material_module
  use Waypoint_module
  
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
    type(waypoint_list_type), pointer  :: waypoints

  end type geomech_realization_type

public :: GeomechRealizCreate, &
          GeomechRealizDestroy, &
          GeomechRealizCreateDiscretization

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
  
  allocate(geomech_realization)
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
! GeomechRealizCreateDiscretization: Creates grid
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine GeomechRealizCreateDiscretization(realization)

  use Geomech_Grid_Aux_module
  
  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  type(geomech_realization_type) :: realization
  
  type(geomech_discretization_type), pointer :: discretization
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  discretization => realization%discretization
  grid => discretization%grid
  option => realization%option
  
  call GeomechDiscretizationCreateDMs(discretization,option)

  
end subroutine GeomechRealizCreateDiscretization

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
  
  if (associated(geomech_realization%geomech_material_property_array)) &
    deallocate(geomech_realization%geomech_material_property_array)
  nullify(geomech_realization%geomech_material_property_array)
  call GeomechanicsMaterialPropertyDestroy(geomech_realization% &
                                           geomech_material_properties)
  call GeomechDiscretizationDestroy(geomech_realization%discretization)
  

end subroutine GeomechRealizDestroy



end module Geomechanics_Realization_module
#endif
! GEOMECH
