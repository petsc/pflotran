module General_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: general_auxvar_type
    PetscReal, pointer :: xmol(:,:)
    PetscReal, pointer :: dsat_dp(:,:)
    PetscReal, pointer :: dden_dp(:,:)
    PetscReal, pointer :: dsat_dt(:)
    PetscReal, pointer :: dden_dt(:)
    PetscReal, pointer :: kvr(:)
    PetscReal, pointer :: dkvr_dp(:)
  end type general_auxvar_type
  
  type, public :: general_parameter_type
    PetscReal, pointer :: sir(:,:)
  end type general_parameter_type
  
  type, public :: general_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc
    type(general_parameter_type), pointer :: general_parameter
    type(general_auxvar_type), pointer :: aux_vars(:)
    type(general_auxvar_type), pointer :: aux_vars_bc(:)
  end type general_type

  public :: GeneralAuxCreate, GeneralAuxDestroy, &
            GeneralAuxVarCompute, GeneralAuxVarInit, &
            GeneralAuxVarCopy

contains


! ************************************************************************** !
!
! GeneralAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
function GeneralAuxCreate()

  use Option_module

  implicit none
  
  type(general_type), pointer :: GeneralAuxCreate
  
  type(general_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)
#ifdef GLENN
  nullify(aux%matrix_buffer)
#endif

  GeneralAuxCreate => aux
  
end function GeneralAuxCreate

! ************************************************************************** !
!
! GeneralAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 01/10/10
!
! ************************************************************************** !
subroutine GeneralAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: aux_var
  type(option_type) :: option
  
end subroutine GeneralAuxVarInit

! ************************************************************************** !
!
! GeneralAuxVarCopy: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine GeneralAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(general_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

end subroutine GeneralAuxVarCopy
  
! ************************************************************************** !
!
! GeneralAuxVarCompute: Computes auxilliary variables for each grid cell
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine GeneralAuxVarCompute(x,aux_var,global_aux_var,&
                                saturation_function,por,perm,option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Saturation_Function_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(general_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var

  PetscReal :: por, perm

#if 0
  PetscInt :: iphase
  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  
  global_aux_var%sat = 0.d0
  global_aux_var%den = 0.d0
  global_aux_var%den_kg = 0.d0
  aux_var%kvr = 0.d0
  kr = 0.d0
 
  global_aux_var%pres = x(1)
  global_aux_var%temp = option%reference_temperature
 
  aux_var%pc = option%reference_pressure - global_aux_var%pres(1)

!***************  Liquid phase properties **************************
  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (aux_var%pc > 0.d0) then
  if (aux_var%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(global_aux_var%pres(1),global_aux_var%sat(1),kr, &
                                   ds_dp,dkr_dp, &
                                   saturation_function, &
                                   por,perm,option)
  else
    iphase = 1
    aux_var%pc = 0.d0
    global_aux_var%sat = 1.d0  
    kr = 1.d0    
    pw = global_aux_var%pres(1)
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call wateos(global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw, &
              hw_dp,hw_dt,option%scale,ierr)

! may need to compute dpsat_dt to pass to VISW
  call psat(global_aux_var%temp(1),sat_pressure,ierr)
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call VISW(global_aux_var%temp(1),pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr) 
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  aux_var%kvr = kr/visl
  
!  aux_var%vis = visl
!  aux_var%dvis_dp = dvis_dp
!  aux_var%kr = kr
!  aux_var%dkr_dp = dkr_dp
  aux_var%dsat_dp = ds_dp

  aux_var%dden_dp = dw_dp
  
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

#endif

end subroutine GeneralAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a general auxilliary object
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(general_auxvar_type) :: aux_var
  
end subroutine AuxVarDestroy

! ************************************************************************** !
!
! GeneralAuxDestroy: Deallocates a general auxilliary object
! author: Glenn Hammond
! date: 01/05/10
!
! ************************************************************************** !
subroutine GeneralAuxDestroy(aux)

  implicit none

  type(general_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%aux_vars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%aux_vars(iaux))
    enddo  
    deallocate(aux%aux_vars)
  endif
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%aux_vars_bc(iaux))
    enddo  
    deallocate(aux%aux_vars_bc)
  endif
  nullify(aux%aux_vars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)

end subroutine GeneralAuxDestroy

end module General_Aux_module
