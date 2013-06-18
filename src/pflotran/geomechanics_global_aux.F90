#ifdef GEOMECH

module Geomechanics_Global_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: geomech_global_auxvar_type
    PetscReal, pointer :: disp_x(:)   ! [m]
    PetscReal, pointer :: disp_y(:)   ! [m]
    PetscReal, pointer :: disp_z(:)   ! [m]
  end type geomech_global_auxvar_type
  
  type, public :: geomech_global_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(geomech_global_auxvar_type), pointer :: aux_vars(:)
    type(geomech_global_auxvar_type), pointer :: aux_vars_bc(:)
    type(geomech_global_auxvar_type), pointer :: aux_vars_ss(:)
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
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)

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
  
  nullify(aux_var%disp_x)
  nullify(aux_var%disp_y)
  nullify(aux_var%disp_z)  

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

  if (associated(aux_var%disp_x) .and. &
      associated(aux_var2%disp_x)) then
    aux_var2%disp_x = aux_var%disp_x
  endif
  if (associated(aux_var%disp_y) .and. &
      associated(aux_var2%disp_y)) then
    aux_var2%disp_y = aux_var%disp_y
  endif
  if (associated(aux_var%disp_z) .and. &
      associated(aux_var2%disp_z)) then
    aux_var2%disp_z = aux_var%disp_z
  endif
  
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
  
  call DeallocateArray(aux_var%disp_x)
  call DeallocateArray(aux_var%disp_y)
  call DeallocateArray(aux_var%disp_z)  

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
  call GeomechGlobalAuxVarDestroy(aux%aux_vars_bc)
  call GeomechGlobalAuxVarDestroy(aux%aux_vars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeomechGlobalAuxDestroy

end module Geomechanics_Global_Aux_module
#endif
