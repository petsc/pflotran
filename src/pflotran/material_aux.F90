module Material_Aux_class
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  PetscInt, parameter, public :: perm_xx_index = 1
  PetscInt, parameter, public :: perm_yy_index = 2
  PetscInt, parameter, public :: perm_zz_index = 3
  PetscInt, parameter, public :: perm_xy_index = 4
  PetscInt, parameter, public :: perm_yz_index = 5
  PetscInt, parameter, public :: perm_xz_index = 6

!  PetscInt, public :: soil_thermal_conductivity_index
!  PetscInt, public :: soil_heat_capacity_index
  PetscInt, public :: soil_compressibility_index
  PetscInt, public :: soil_reference_pressure_index
  PetscInt, public :: max_material_index
 
  type, public :: material_auxvar_type
    PetscReal :: volume
    PetscReal :: porosity
    PetscReal :: tortuosity
    PetscReal :: soil_particle_density
    PetscReal, pointer :: permeability(:)
    PetscReal, pointer :: sat_func_prop(:)
    PetscReal, pointer :: soil_properties(:) ! den, therm. cond., heat cap.
!    procedure(SaturationFunction), nopass, pointer :: SaturationFunction
  contains
    procedure, public :: PermeabilityTensorToScalar => &
                           MaterialPermTensorToScalar
  end type material_auxvar_type
  
  type, public :: material_parameter_type
    PetscReal, pointer :: soil_residual_saturation(:,:)
    PetscReal, pointer :: soil_heat_capacity(:) ! MJ/kg rock-K
    PetscReal, pointer :: soil_thermal_conductivity(:,:) ! W/M
  end type material_parameter_type  
  
  type, public :: material_type
    PetscInt :: num_aux
    type(material_parameter_type), pointer :: material_parameter
    class(material_auxvar_type), pointer :: auxvars(:)
  end type material_type
  
  public :: MaterialAuxCreate, &
            MaterialAuxVarInit, &
            MaterialAuxVarCopy, &
            MaterialAuxVarStrip, &
            MaterialAuxDestroy
  
contains

! ************************************************************************** !

function MaterialAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  type(material_type), pointer :: MaterialAuxCreate
  
  type(material_type), pointer :: aux

  allocate(aux)
  nullify(aux%auxvars)
  allocate(aux%material_parameter)
  nullify(aux%material_parameter%soil_residual_saturation)
  nullify(aux%material_parameter%soil_heat_capacity)
  nullify(aux%material_parameter%soil_thermal_conductivity)
  aux%num_aux = 0

  MaterialAuxCreate => aux
  
end function MaterialAuxCreate

! ************************************************************************** !

subroutine MaterialAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%volume = -999.d0
  auxvar%porosity = -999.d0
  auxvar%tortuosity = -999.d0
  auxvar%soil_particle_density = -999.d0
  if (option%iflowmode /= NULL_MODE) then
    allocate(auxvar%permeability(3))
    auxvar%permeability = -999.d0
  else
    nullify(auxvar%permeability)
  endif
  nullify(auxvar%sat_func_prop)
  if (max_material_index > 0) then
    allocate(auxvar%soil_properties(max_material_index))
    ! initialize these to zero for now
    auxvar%soil_properties = 0.d0
  else
    nullify(auxvar%soil_properties)
  endif
  
end subroutine MaterialAuxVarInit

! ************************************************************************** !

subroutine MaterialAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option
  
  auxvar2%volume = auxvar%volume
  auxvar2%porosity = auxvar%porosity
  auxvar2%tortuosity = auxvar%tortuosity
  if (associated(auxvar%permeability)) then
    auxvar2%permeability = auxvar%permeability
  endif
  if (associated(auxvar%sat_func_prop)) then
    auxvar2%sat_func_prop = auxvar%sat_func_prop
  endif
  if (associated(auxvar%soil_properties)) then
    auxvar2%soil_properties = auxvar%soil_properties
  endif

end subroutine MaterialAuxVarCopy

! ************************************************************************** !

subroutine MaterialPermTensorToScalar(material_auxvar,dist, &
                                      scalar_permeability)
  ! 
  ! Transforms a diagonal permeability tensor to a scalar through a dot 
  ! product.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(out) :: scalar_permeability
  
  scalar_permeability = &
            material_auxvar%permeability(perm_xx_index)*dabs(dist(1))+ &
            material_auxvar%permeability(perm_yy_index)*dabs(dist(2))+ &
            material_auxvar%permeability(perm_zz_index)*dabs(dist(3))

end subroutine MaterialPermTensorToScalar

! ************************************************************************** !

subroutine MaterialAuxVarStrip(auxvar)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Utility_module, only : DeallocateArray
  
  implicit none

  class(material_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%permeability)
  call DeallocateArray(auxvar%sat_func_prop)
  call DeallocateArray(auxvar%soil_properties)
  
end subroutine MaterialAuxVarStrip

! ************************************************************************** !

subroutine MaterialAuxDestroy(aux)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/11
  ! 

  implicit none

  type(material_type), pointer :: aux
  
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MaterialAuxVarStrip(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
    
  if (associated(aux%material_parameter)) then
    if (associated(aux%material_parameter%soil_residual_saturation)) &
      deallocate(aux%material_parameter%soil_residual_saturation)
    nullify(aux%material_parameter%soil_residual_saturation)
    if (associated(aux%material_parameter%soil_heat_capacity)) &
      deallocate(aux%material_parameter%soil_heat_capacity)
    nullify(aux%material_parameter%soil_heat_capacity)
    if (associated(aux%material_parameter%soil_thermal_conductivity)) &
      deallocate(aux%material_parameter%soil_thermal_conductivity)
    nullify(aux%material_parameter%soil_thermal_conductivity)
    deallocate(aux%material_parameter)
  endif
  nullify(aux%material_parameter)
  
  deallocate(aux)
  nullify(aux)

end subroutine MaterialAuxDestroy

end module Material_Aux_class
