#ifdef GEOMECH

module Geomechanics_Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: geomech_global_auxvar_type
    PetscReal, pointer :: disp_vector(:)   ! [m]
    PetscReal, pointer :: strain(:,:)      ! dimensionless
    PetscReal, pointer :: stress(:,:)      ! [N]
  end type geomech_global_auxvar_type
  
  type, public :: geomech_global_type
    PetscInt :: num_aux
    type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  end type geomech_global_type
  
  interface GeomechGlobalAuxVarDestroy
    module procedure GeomechGlobalAuxVarSingleDestroy
    module procedure GeomechGlobalAuxVarArrayDestroy
  end interface GeomechGlobalAuxVarDestroy

  public :: GeomechGlobalAuxCreate, &
            GeomechGlobalAuxDestroy, &
            GeomechGlobalAuxVarInit, &
            GeomechGlobalAuxVarCopy, &
            GeomechGlobalAuxVarDestroy, &
            GeomechGlobalAuxVarStrip

contains

! ************************************************************************** !
!
! GeomechGlobalAuxCreate: Creates a geomech global aux
! author: Satish Karra, LANL
! date: 06/14/13
!
! ************************************************************************** !
function GeomechGlobalAuxCreate()

  implicit none
  
  type(geomech_global_type), pointer       :: GeomechGlobalAuxCreate
  
  type(geomech_global_type), pointer       :: aux

  allocate(aux) 
  aux%num_aux = 0
  nullify(aux%aux_vars)

  GeomechGlobalAuxCreate => aux
  
end function GeomechGlobalAuxCreate

! ************************************************************************** !
!
! GeomechGlobalAuxVarInit: Initializes a geomech global aux
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechGlobalAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(geomech_global_auxvar_type)       :: aux_var
  type(option_type)                      :: option
  
  allocate(aux_var%disp_vector(option%ngeomechdof))
  allocate(aux_var%strain(option%ngeomechdof,option%ngeomechdof))
  allocate(aux_var%stress(option%ngeomechdof,option%ngeomechdof))
  aux_var%disp_vector = 0.d0
  aux_var%strain = 0.d0
  aux_var%stress = 0.d0
  
end subroutine GeomechGlobalAuxVarInit

! ************************************************************************** !
!
! GeomechGlobalAuxVarCopy: Copies a geomech global aux to another
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !  
subroutine GeomechGlobalAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(geomech_global_auxvar_type)      :: aux_var, aux_var2
  type(option_type)                     :: option

  aux_var%disp_vector = aux_var2%disp_vector
  aux_var%strain = aux_var2%strain
  aux_var%stress = aux_var2%stress
  
end subroutine GeomechGlobalAuxVarCopy

! ************************************************************************** !
!
! GeomechGlobalAuxVarSingleDestroy: Destroys a geomech global aux
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechGlobalAuxVarSingleDestroy(aux_var)

  implicit none

  type(geomech_global_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call GeomechGlobalAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)

end subroutine GeomechGlobalAuxVarSingleDestroy
  
! ************************************************************************** !
!
! GeomechGlobalAuxVarArrayDestroy: Destroys an array of geomech global auxvar
!                                  type
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechGlobalAuxVarArrayDestroy(aux_vars)

  implicit none

  type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call GeomechGlobalAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)

end subroutine GeomechGlobalAuxVarArrayDestroy
  
! ************************************************************************** !
!
! GeomechGlobalAuxVarStrip: Strips a geomech global auxvar
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechGlobalAuxVarStrip(aux_var)

  use Utility_module, only: DeallocateArray

  implicit none

  type(geomech_global_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%disp_vector)
  call DeallocateArray(aux_var%strain)
  call DeallocateArray(aux_var%stress)

end subroutine GeomechGlobalAuxVarStrip

! ************************************************************************** !
!
! GeomechGlobalAuxDestroy: Destroys a geomech global type
! author: Satish Karra, LANL
! date: 06/14/13
! 
! ************************************************************************** !
subroutine GeomechGlobalAuxDestroy(aux)

  implicit none

  type(geomech_global_type), pointer       :: aux
  
  if (.not.associated(aux)) return
  
  call GeomechGlobalAuxVarDestroy(aux%aux_vars)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeomechGlobalAuxDestroy

end module Geomechanics_Global_Aux_module
#endif
