module Material_module
 
  use Region_module
 
  implicit none

  private

#include "definitions.h"
! Apparently the PETSc authors believe that Fortran 90 modules should ensure
! that PETSC_AVOID_DECLARATIONS and PETSC_AVOID_MPIF_H are defined when the
! PETSc header files are included.  I can get around this, though, by making
! the definitions in these headers private.
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscda.h"
#include "include/finclude/petscda.h90"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petsclog.h"
 
  type, public :: material_type
    integer :: id
    character(len=MAXNAMELENGTH) :: name
    real*8 :: permeability(3,3)
    real*8 :: permeability_pwr
    real*8 :: porosity
    real*8 :: tortuosity
    integer :: ithrm
    integer :: icap
    type(material_type), pointer :: next
  end type material_type
  
  type, public :: thermal_property_type
    integer :: id
    real*8 :: rock_density
    real*8 :: spec_heat
    real*8 :: therm_cond_dry
    real*8 :: therm_cond_wet
    real*8 :: pore_compress
    real*8 :: pore_expansivity
    real*8 :: tort_bin_diff
    real*8 :: vap_air_diff_coef
    real*8 :: exp_binary_diff
    real*8 :: enh_binary_diff_coef
    type(thermal_property_type), pointer :: next
  end type thermal_property_type
  
  type, public :: saturation_function_type
    integer :: id
    character(len=MAXNAMELENGTH) :: ctype
    integer :: itype
    real*8, pointer :: Sr(:)
    real*8 :: m
    real*8 :: lambda
    real*8 :: alpha
    real*8 :: pcwmax
    real*8 :: betac
    real*8 :: power
    integer :: ihist 
    type(saturation_function_type), pointer :: next
  end type saturation_function_type
  
  public :: MaterialCreate, ThermalPropertyCreate, SaturationFunctionCreate, &
            MaterialDestroy, ThermalPropertyDestroy, &
            SaturationFunctionDestroy, &
            MaterialAddToList, ThermalAddPropertyToList, &
            SaturationFunctionAddToList, &
            MaterialGetPtrFromList
  
contains

! ************************************************************************** !
!
! MaterialCreate: Creates a material
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialCreate()
  
  implicit none

  type(material_type), pointer :: MaterialCreate
  
  type(material_type), pointer :: material
  
  allocate(material)
  material%id = 0
  material%name = ""
  material%permeability = 0.d0
  material%permeability_pwr = 0.d0
  material%porosity = 0.d0
  material%tortuosity = 0.d0
  material%ithrm = 0
  material%icap = 0
  nullify(material%next)
  MaterialCreate => material

end function MaterialCreate

! ************************************************************************** !
!
! MaterialCreate: Creates a thermal property
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function ThermalPropertyCreate()
  
  implicit none

  type(thermal_property_type), pointer :: ThermalPropertyCreate
  
  type(thermal_property_type), pointer :: thermal_property
  
  allocate(thermal_property)
  thermal_property%id = 0
  thermal_property%rock_density = 0.d0
  thermal_property%spec_heat = 0.d0
  thermal_property%therm_cond_dry = 0.d0
  thermal_property%therm_cond_wet = 0.d0
  thermal_property%pore_compress = 0.d0
  thermal_property%pore_expansivity = 0.d0
  thermal_property%tort_bin_diff = 0.d0
  thermal_property%vap_air_diff_coef = 0.d0
  thermal_property%exp_binary_diff = 0.d0
  thermal_property%enh_binary_diff_coef = 0.d0
  nullify(thermal_property%next)
  ThermalPropertyCreate => thermal_property

end function ThermalPropertyCreate

! ************************************************************************** !
!
! SaturationFunctionCreate: Creates a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function SaturationFunctionCreate(option)
  
  use Option_module
  
  implicit none

  type(saturation_function_type), pointer :: SaturationFunctionCreate
  type(option_type) :: option
  
  type(saturation_function_type), pointer :: saturation_function
  
  allocate(saturation_function)
  saturation_function%id = 0
  saturation_function%ctype = ""
  saturation_function%itype = 0
  allocate(saturation_function%Sr(option%ndof))
  saturation_function%Sr = 0.d0
  saturation_function%m = 0.d0
  saturation_function%lambda = 0.d0
  saturation_function%alpha = 0.d0
  saturation_function%pcwmax = 0.d0
  saturation_function%betac = 0.d0
  saturation_function%power = 0.d0
  saturation_function%ihist = 0
  nullify(saturation_function%next)
  SaturationFunctionCreate => saturation_function

end function SaturationFunctionCreate

! ************************************************************************** !
!
! MaterialAddToList: Adds a thermal property to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine MaterialAddToList(material,list)

  implicit none
  
  type(material_type), pointer :: material
  type(material_type), pointer :: list

  type(material_type), pointer :: cur_material
  
  if (associated(list)) then
    cur_material => list
    ! loop to end of list
    do
      if (.not.associated(cur_material%next)) exit
      cur_material => cur_material%next
    enddo
    cur_material%next => material
  else
    list => material
  endif
  
end subroutine MaterialAddToList

! ************************************************************************** !
!
! ThermalAddPropertyToList: Adds a thermal property to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine ThermalAddPropertyToList(thermal_property,list)

  implicit none
  
  type(thermal_property_type), pointer :: thermal_property
  type(thermal_property_type), pointer :: list

  type(thermal_property_type), pointer :: cur_thermal_property
  
  if (associated(list)) then
    cur_thermal_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_thermal_property%next)) exit
      cur_thermal_property => cur_thermal_property%next
    enddo
    cur_thermal_property%next => thermal_property
  else
    list => thermal_property
  endif
  
end subroutine ThermalAddPropertyToList

! ************************************************************************** !
!
! SaturationFunctionAddToList: Adds a saturation function to linked list
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine SaturationFunctionAddToList(saturation_function,list)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  type(saturation_function_type), pointer :: list

  type(saturation_function_type), pointer :: cur_saturation_function
  
  if (associated(list)) then
    cur_saturation_function => list
    ! loop to end of list
    do
      if (.not.associated(cur_saturation_function%next)) exit
      cur_saturation_function => cur_saturation_function%next
    enddo
    cur_saturation_function%next => saturation_function
  else
    list => saturation_function
  endif
  
end subroutine SaturationFunctionAddToList

! ************************************************************************** !
!
! MaterialGetPtrFromList: Returns a pointer to the material matching 
!                         material_name
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
function MaterialGetPtrFromList(material_name,material_list)

  use Fileio_module

  implicit none
  
  type(material_type), pointer :: MaterialGetPtrFromList
  character(len=MAXNAMELENGTH) :: material_name
  type(material_type), pointer :: material_list

  type(material_type), pointer :: material
    
  nullify(MaterialGetPtrFromList)
  material => material_list
  
  do 
    if (.not.associated(material)) exit
    if (fiStringCompare(material%name,material_name, &
                        len_trim(material_name))) then
      MaterialGetPtrFromList => material
      return
    endif
    material => material%next
  enddo
  
end function MaterialGetPtrFromList

! ************************************************************************** !
!
! MaterialDestroy: Destroys a material
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine MaterialDestroy(material)

  implicit none
  
  type(material_type), pointer :: material
  
  if (.not.associated(material)) return
  
  call MaterialDestroy(material%next)
  
  deallocate(material)
  nullify(material)
  
end subroutine MaterialDestroy

! ************************************************************************** !
!
! SaturationFunctionDestroy: Destroys a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine ThermalPropertyDestroy(thermal_property)

  implicit none
  
  type(thermal_property_type), pointer :: thermal_property
  
  if (.not.associated(thermal_property)) return
  
  call ThermalPropertyDestroy(thermal_property%next)

  deallocate(thermal_property)
  nullify(thermal_property)
  
end subroutine ThermalPropertyDestroy

! ************************************************************************** !
!
! SaturationFunctionDestroy: Destroys a saturation function
! author: Glenn Hammond
! date: 11/02/07
!
! ************************************************************************** !
recursive subroutine SaturationFunctionDestroy(saturation_function)

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  
  if (.not.associated(saturation_function)) return
  
  call SaturationFunctionDestroy(saturation_function%next)
    
  if (associated(saturation_function%Sr)) deallocate(saturation_function%Sr)
  nullify(saturation_function%Sr)
    
  deallocate(saturation_function)
  nullify(saturation_function)
  
end subroutine SaturationFunctionDestroy

end module Material_module
