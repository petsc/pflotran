#ifdef SURFACE_FLOW

module Surface_Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"

  type, public :: surface_global_auxvar_type
    PetscInt :: istate
    PetscReal, pointer :: head(:)   ! [m]
    PetscReal, pointer :: temp(:)   ! [C]
    PetscReal, pointer :: den_kg(:) ! [kg/m^3]
  end type surface_global_auxvar_type
  
  type, public :: surface_global_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(surface_global_auxvar_type), pointer :: aux_vars(:)
    type(surface_global_auxvar_type), pointer :: aux_vars_bc(:)
    type(surface_global_auxvar_type), pointer :: aux_vars_ss(:)
  end type surface_global_type
  
  interface SurfaceGlobalAuxVarDestroy
    module procedure SurfaceGlobalAuxVarSingleDestroy
    module procedure SurfaceGlobalAuxVarArrayDestroy
  end interface SurfaceGlobalAuxVarDestroy

  public :: SurfaceGlobalAuxCreate, &
            SurfaceGlobalAuxDestroy, &
            SurfaceGlobalAuxVarInit, &
            SurfaceGlobalAuxVarCopy, &
            SurfaceGlobalAuxVarDestroy, &
            SurfaceGlobalAuxVarStrip

contains

! ************************************************************************** !

function SurfaceGlobalAuxCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Option_module

  implicit none
  
  type(surface_global_type), pointer :: SurfaceGlobalAuxCreate
  
  type(surface_global_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)

  SurfaceGlobalAuxCreate => aux
  
end function SurfaceGlobalAuxCreate

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarInit(aux_var,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Option_module

  implicit none
  
  type(surface_global_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%istate = 0

  allocate(aux_var%head(option%nphase))
  aux_var%head = 0.d0
  allocate(aux_var%temp(ONE_INTEGER))
  aux_var%temp = option%reference_temperature
  allocate(aux_var%den_kg(option%nphase))
  aux_var%den_kg = 0.d0

end subroutine SurfaceGlobalAuxVarInit

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarCopy(aux_var,aux_var2,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Option_module

  implicit none
  
  type(surface_global_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%istate = aux_var%istate
  aux_var2%head = aux_var%head
  aux_var2%temp = aux_var%temp
  aux_var2%den_kg = aux_var%den_kg

end subroutine SurfaceGlobalAuxVarCopy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarSingleDestroy(aux_var)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  implicit none

  type(surface_global_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call SurfaceGlobalAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)

end subroutine SurfaceGlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarArrayDestroy(aux_vars)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  implicit none

  type(surface_global_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call SurfaceGlobalAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)

end subroutine SurfaceGlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarStrip(aux_var)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(surface_global_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%head)
  call DeallocateArray(aux_var%temp)
  call DeallocateArray(aux_var%den_kg)


end subroutine SurfaceGlobalAuxVarStrip

! ************************************************************************** !

subroutine SurfaceGlobalAuxDestroy(aux)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 02/28/13
  ! 

  implicit none

  type(surface_global_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  call SurfaceGlobalAuxVarDestroy(aux%aux_vars)
  call SurfaceGlobalAuxVarDestroy(aux%aux_vars_bc)
  call SurfaceGlobalAuxVarDestroy(aux%aux_vars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine SurfaceGlobalAuxDestroy

end module Surface_Global_Aux_module
#endif
