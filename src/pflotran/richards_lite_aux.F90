module Richards_Lite_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: richards_lite_auxvar_type
    PetscReal :: pres
    PetscReal :: temp
    PetscReal :: sat
    PetscReal :: den
    PetscReal :: den_kg
    PetscReal :: avgmw
    PetscReal :: pc
!    PetscReal :: vis
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
    PetscReal :: dkvr_dp
  end type richards_lite_auxvar_type
  
  type, public :: richards_lite_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscTruth :: aux_vars_up_to_date
    PetscTruth :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc
    type(richards_lite_auxvar_type), pointer :: aux_vars(:)
    type(richards_lite_auxvar_type), pointer :: aux_vars_bc(:)
  end type richards_lite_type

  public :: RichardsLiteAuxCreate, RichardsLiteAuxDestroy, &
            RichardsLiteAuxVarCompute, RichardsLiteAuxVarInit, &
            RichardsLiteAuxVarCopy

contains


! ************************************************************************** !
!
! RichardsLiteAuxVarCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function RichardsLiteAuxCreate()

  use Option_module

  implicit none
  
  type(richards_lite_type), pointer :: RichardsLiteAuxCreate
  
  type(richards_lite_type), pointer :: aux

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

  RichardsLiteAuxCreate => aux
  
end function RichardsLiteAuxCreate

! ************************************************************************** !
!
! RichardsLiteAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsLiteAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(richards_lite_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%pres = 0.d0
  aux_var%temp = 0.d0
  aux_var%sat = 0.d0
  aux_var%den = 0.d0
  aux_var%den_kg = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%pc = 0.d0
!  aux_var%kr = 0.d0
!  aux_var%dkr_dp = 0.d0
!  aux_var%vis = 0.d0
!  aux_var%dvis_dp = 0.d0
  aux_var%kvr = 0.d0
  aux_var%dsat_dp = 0.d0
  aux_var%dden_dp = 0.d0
  aux_var%dkvr_dp = 0.d0

end subroutine RichardsLiteAuxVarInit

! ************************************************************************** !
!
! RichardsLiteAuxVarCopy: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsLiteAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(richards_lite_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%den_kg = aux_var%den_kg
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%pc = aux_var%pc
!  aux_var2%kr = aux_var%kr
!  aux_var2%dkr_dp = aux_var%dkr_dp
!  aux_var2%vis = aux_var%vis
!  aux_var2%dvis_dp = aux_var%dvis_dp
  aux_var2%kvr = aux_var%kvr
  aux_var2%dsat_dp = aux_var%dsat_dp
  aux_var2%dden_dp = aux_var%dden_dp
  aux_var2%dkvr_dp = aux_var%dkvr_dp

end subroutine RichardsLiteAuxVarCopy
  
! ************************************************************************** !
!
! RichardsLiteAuxVarCompute: Computes auxilliary variables for each grid cell
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsLiteAuxVarCompute(x,aux_var,iphase,saturation_function, &
                                     por,perm,option)

  use Option_module
  use water_eos_module
  use Material_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(richards_lite_auxvar_type) :: aux_var
  PetscInt :: iphase
  PetscReal :: por, perm

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  
  aux_var%sat = 0.d0
  aux_var%den = 0.d0
  aux_var%den_kg = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%kvr = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)
  aux_var%temp = 25.d0
 
  aux_var%pc = option%pref - aux_var%pres

!***************  Liquid phase properties **************************
  !geh aux_var%avgmw = option%fmwh2o  ! hardwire for comparison with old code
  aux_var%avgmw = 18.0153d0

  pw = option%pref
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (aux_var%pc > 0.d0) then
  if (aux_var%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(aux_var%pres,aux_var%sat,kr, &
                                   ds_dp,dkr_dp, &
                                   saturation_function, &
                                   por,perm,option)
  else
    iphase = 1
    aux_var%pc = 0.d0
    aux_var%sat = 1.d0  
    kr = 1.d0    
    pw = aux_var%pres
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call wateos(aux_var%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
              option%scale,ierr)

! may need to compute dpsat_dt to pass to VISW
  call psat(aux_var%temp,sat_pressure,ierr)
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call VISW(aux_var%temp,pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr) 
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  aux_var%den = dw_mol
  aux_var%den_kg = dw_kg
  aux_var%kvr = kr/visl
  
!  aux_var%vis = visl
!  aux_var%dvis_dp = dvis_dp
!  aux_var%kr = kr
!  aux_var%dkr_dp = dkr_dp
  aux_var%dsat_dp = ds_dp

  aux_var%dden_dp = dw_dp
  
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

end subroutine RichardsLiteAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(richards_lite_auxvar_type) :: aux_var
  

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! RichardsLiteAuxDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsLiteAuxDestroy(aux)

  implicit none

  type(richards_lite_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call AuxVarDestroy(aux%aux_vars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call AuxVarDestroy(aux%aux_vars_bc(iaux))
  enddo  
  
  if (associated(aux%aux_vars)) deallocate(aux%aux_vars)
  nullify(aux%aux_vars)
  if (associated(aux%aux_vars_bc)) deallocate(aux%aux_vars_bc)
  nullify(aux%aux_vars_bc)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
    
end subroutine RichardsLiteAuxDestroy

end module Richards_Lite_Aux_module
