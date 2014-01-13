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

  PetscInt, public :: rock_density_index
  PetscInt, public :: rock_thermal_conductivity_index
  PetscInt, public :: rock_heat_capacity_index
  PetscInt, public :: rock_compressibility_index
 
  type, public :: material_auxvar_type
    PetscReal :: volume
    PetscReal :: porosity
    PetscReal :: tortuosity
    PetscReal, pointer :: permeability(:)
    PetscReal, pointer :: sat_func_prop(:)
    PetscReal, pointer :: rock_properties(:) ! den, therm. cond., heat cap.
!    procedure(SaturationFunction), nopass, pointer :: SaturationFunction
  contains
    procedure, public :: PermeabilityTensorToScalar => &
                           MaterialPermTensorToScalar
  end type material_auxvar_type
  
  type, public :: material_parameter_type
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: dencpr(:) ! MJ/kg rock-K
  end type material_parameter_type  
  
  type, public :: material_type
    PetscInt :: num_aux
    type(material_parameter_type), pointer :: material_parameter
    class(material_auxvar_type), pointer :: aux_vars(:)
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
  nullify(aux%aux_vars)
  allocate(aux%material_parameter)
  nullify(aux%material_parameter%sir)
  nullify(aux%material_parameter%dencpr)
  aux%num_aux = 0

  MaterialAuxCreate => aux
  
end function MaterialAuxCreate

! ************************************************************************** !

subroutine MaterialAuxVarInit(aux_var,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%volume = -999.d0
  aux_var%porosity = -999.d0
  aux_var%tortuosity = -999.d0
  if (option%iflowmode /= NULL_MODE) then
    allocate(aux_var%permeability(3))
    aux_var%permeability = -999.d0
  else
    nullify(aux_var%permeability)
  endif
  nullify(aux_var%sat_func_prop)
  nullify(aux_var%rock_properties)
  
end subroutine MaterialAuxVarInit

! ************************************************************************** !

subroutine MaterialAuxVarCopy(aux_var,aux_var2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option
  
  aux_var2%volume = aux_var%volume
  aux_var2%porosity = aux_var%porosity
  aux_var2%tortuosity = aux_var%tortuosity
  if (associated(aux_var%permeability)) then
    aux_var2%permeability = aux_var%permeability
  endif
  if (associated(aux_var%sat_func_prop)) then
    aux_var2%sat_func_prop = aux_var%sat_func_prop
  endif
  if (associated(aux_var%rock_properties)) then
    aux_var2%rock_properties = aux_var%rock_properties
  endif

end subroutine MaterialAuxVarCopy

! ************************************************************************** !

subroutine MaterialPermTensorToScalar(material_aux_var,dist, &
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
  
  class(material_auxvar_type) :: material_aux_var
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(out) :: scalar_permeability
  
  scalar_permeability = &
            material_aux_var%permeability(perm_xx_index)*dabs(dist(1))+ &
            material_aux_var%permeability(perm_yy_index)*dabs(dist(2))+ &
            material_aux_var%permeability(perm_zz_index)*dabs(dist(3))

end subroutine MaterialPermTensorToScalar

! ************************************************************************** !

subroutine MaterialAuxVarStrip(aux_var)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Utility_module, only : DeallocateArray
  
  implicit none

  class(material_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%permeability)
  call DeallocateArray(aux_var%sat_func_prop)
  call DeallocateArray(aux_var%rock_properties)
  
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
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call MaterialAuxVarStrip(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
    
  if (associated(aux%material_parameter)) then
    if (associated(aux%material_parameter%sir)) &
      deallocate(aux%material_parameter%sir)
    nullify(aux%material_parameter%sir)
    if (associated(aux%material_parameter%dencpr)) &
      deallocate(aux%material_parameter%dencpr)
    nullify(aux%material_parameter%dencpr)
    deallocate(aux%material_parameter)
  endif
  nullify(aux%material_parameter)
  
  deallocate(aux)
  nullify(aux)

end subroutine MaterialAuxDestroy

end module Material_Aux_class
