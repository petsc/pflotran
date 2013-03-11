#ifdef SURFACE_FLOW

module Surface_TH_Aux_module

  implicit none

  private

#include "definitions.h"

  type, public :: Surface_TH_auxvar_type
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
  end type Surface_TH_auxvar_type

  type, public :: Surface_TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(Surface_TH_auxvar_type), pointer :: aux_vars(:)
    type(Surface_TH_auxvar_type), pointer :: aux_vars_bc(:)
    type(Surface_TH_auxvar_type), pointer :: aux_vars_ss(:)
  end type Surface_TH_type

  public :: SurfaceTHAuxCreate, &
            SurfaceTHAuxDestroy, &
            SurfaceTHAuxVarCompute, &
            SurfaceTHAuxVarInit, &
            SurfaceTHAuxVarCopy

contains

! ************************************************************************** !
!> This routine creates an empty object
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
function SurfaceTHAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(Surface_TH_type), pointer :: SurfaceTHAuxCreate
  
  type(Surface_TH_type), pointer :: aux

  allocate(aux)
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  nullify(aux%aux_vars_ss)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  SurfaceTHAuxCreate => aux
  
end function SurfaceTHAuxCreate

! ************************************************************************** !
!> This routine initilizes an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(Surface_TH_auxvar_type) :: aux_var
  type(option_type) :: option

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%pc = 0.d0
  aux_var%vis = 0.d0

end subroutine SurfaceTHAuxVarInit

! ************************************************************************** !
!> This routine makes a copy of an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(Surface_TH_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
  aux_var2%vis = aux_var%vis

end subroutine SurfaceTHAuxVarCopy

! ************************************************************************** !
!> This routine computes values for auxiliary variables.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxVarCompute(xx,aux_var,global_aux_var, &
                                  !saturation_function,por,perm,&
                                  option)

  use Option_module
  use Surface_Global_Aux_module
  use Water_EOS_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: xx(option%nflowdof)
  type(Surface_TH_auxvar_type) :: aux_var
  type(surface_global_auxvar_type) :: global_aux_var
  PetscReal :: por, perm

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  
  global_aux_var%den_kg = 0.d0

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  kr = 0.d0
 
  global_aux_var%head(1) = xx(1)
  global_aux_var%temp(1) = xx(2)
!  global_aux_var%temp(1) = option%reference_temperature
 

!***************  Liquid phase properties **************************

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

  call wateos_noderiv(global_aux_var%temp(1),pw,dw_kg,dw_mol,hw,option%scale,ierr)
  global_aux_var%den_kg(1) = dw_kg
  
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale
  
end subroutine SurfaceTHAuxVarCompute

! ************************************************************************** !
!> This routine deallocates an object.
!!
!> @author
!! Gautam Bisht, LBNL
!!
!! date: 03/07/13
! ************************************************************************** !
subroutine SurfaceTHAuxDestroy(aux)

  implicit none

  type(Surface_TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) deallocate(aux%aux_vars_bc)
  nullify(aux%aux_vars_bc)
  if (associated(aux%aux_vars_ss)) deallocate(aux%aux_vars_ss)
  nullify(aux%aux_vars_ss)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

  deallocate(aux)
  nullify(aux)  

  end subroutine SurfaceTHAuxDestroy

end module Surface_TH_Aux_module

#endif
