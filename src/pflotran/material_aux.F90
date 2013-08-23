module Material_Aux_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
 
  type, public :: material_auxvar_type
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: dencpr(:) ! MJ/kg rock-K
  end type material_auxvar_type
  
  type, public :: material_parameter_type
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: dencpr(:) ! MJ/kg rock-K
  end type material_parameter_type

  type, public :: material_type
    PetscInt :: num_aux
    type(material_parameter_type), pointer :: material_parameter
    type(material_auxvar_type), pointer :: aux_vars(:)    
  end type material_type  
  
  public :: MaterialAuxCreate, &
            MaterialAuxDestroy
  
contains

! ************************************************************************** !
!
! MaterialAuxCreate: Allocate and initialize auxiliary object
! author: Glenn Hammond
! date: 03/02/11
!
! ************************************************************************** !
function MaterialAuxCreate()

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
!
! MaterialAuxVarInit: Initialize auxiliary object
! author: Glenn Hammond
! date: 03/02/11
!
! ************************************************************************** !
subroutine MaterialAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(material_auxvar_type) :: aux_var
  type(option_type) :: option
  
end subroutine MaterialAuxVarInit

! ************************************************************************** !
!
! MaterialAuxVarCopy: Copies an auxiliary variable
! author: Glenn Hammond
! date: 03/02/11
!
! ************************************************************************** !  
subroutine MaterialAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(material_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

end subroutine MaterialAuxVarCopy

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a material auxiliary object
! author: Glenn Hammond
! date: 03/02/11
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(material_auxvar_type) :: aux_var
  
end subroutine AuxVarDestroy

! ************************************************************************** !
!
! MaterialAuxDestroy: Deallocates a material auxiliary object
! author: Glenn Hammond
! date: 03/02/11
!
! ************************************************************************** !
subroutine MaterialAuxDestroy(aux)

  implicit none

  type(material_type), pointer :: aux
  
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%aux_vars(iaux))
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

end module Material_Aux_module
