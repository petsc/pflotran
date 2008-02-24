module Richards_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: richards_auxvar_type
    PetscReal :: pres
    PetscReal :: temp
    PetscReal :: sat
    PetscReal :: den
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
!    PetscReal :: vis
!    PetscReal :: dvis_dp
!    PetscReal :: kr
!    PetscReal :: dkr_dp
    PetscReal :: kvr
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dt
    PetscReal :: dh_dp
    PetscReal :: dh_dt
    PetscReal :: du_dp
    PetscReal :: du_dt
    PetscReal, pointer :: xmol(:)
    PetscReal, pointer :: diff(:)
  end type richards_auxvar_type
  
  type, public :: richards_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    logical :: aux_vars_up_to_date
    logical :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc
    type(richards_auxvar_type), pointer :: aux_vars(:)
    type(richards_auxvar_type), pointer :: aux_vars_bc(:)
  end type richards_type

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy

contains


! ************************************************************************** !
!
! RichardsAuxCreate: Allocate and initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
function RichardsAuxCreate()

  use Option_module

  implicit none
  
  type(richards_type), pointer :: RichardsAuxCreate
  
  type(richards_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = .false.
  aux%inactive_cells_exist = .false.
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  RichardsAuxCreate => aux
  
end function RichardsAuxCreate

! ************************************************************************** !
!
! RichardsAuxVarInit: Initialize auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%pres = 0.d0
  aux_var%temp = 0.d0
  aux_var%sat = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%pc = 0.d0
!  aux_var%kr = 0.d0
!  aux_var%dkr_dp = 0.d0
!  aux_var%vis = 0.d0
!  aux_var%dvis_dp = 0.d0
  aux_var%kvr = 0.d0
  aux_var%dsat_dp = 0.d0
  aux_var%dden_dp = 0.d0
  aux_var%dden_dt = 0.d0
  aux_var%dkvr_dp = 0.d0
  aux_var%dkvr_dt = 0.d0
  aux_var%dh_dp = 0.d0
  aux_var%dh_dt = 0.d0
  aux_var%du_dp = 0.d0
  aux_var%du_dt = 0.d0    
  allocate(aux_var%xmol(option%nspec))
  aux_var%xmol = 0.d0
  allocate(aux_var%diff(option%nspec))
  aux_var%diff = 0.d0

end subroutine RichardsAuxVarInit

! ************************************************************************** !
!
! RichardsAuxVarCopy: Copies an auxilliary variable
! author: Glenn Hammond
! date: 12/13/07
!
! ************************************************************************** !  
subroutine RichardsAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(richards_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var2%pres = aux_var%pres
  aux_var2%temp = aux_var%temp
  aux_var2%sat = aux_var%sat
  aux_var2%den = aux_var%den
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
!  aux_var2%kr = aux_var%kr
!  aux_var2%dkr_dp = aux_var%dkr_dp
!  aux_var2%vis = aux_var%vis
!  aux_var2%dvis_dp = aux_var%dvis_dp
  aux_var2%kvr = aux_var%kvr
  aux_var2%dsat_dp = aux_var%dsat_dp
  aux_var2%dden_dp = aux_var%dden_dp
  aux_var2%dden_dt = aux_var%dden_dt
  aux_var2%dkvr_dp = aux_var%dkvr_dp
  aux_var2%dkvr_dt = aux_var%dkvr_dt
  aux_var2%dh_dp = aux_var%dh_dp
  aux_var2%dh_dt = aux_var%dh_dt
  aux_var2%du_dp = aux_var%du_dp
  aux_var2%du_dt = aux_var%du_dt  
  aux_var2%xmol = aux_var%xmol
  aux_var2%diff = aux_var%diff

end subroutine RichardsAuxVarCopy

! ************************************************************************** !
!
! RichardsAuxVarCompute: Computes auxilliary variables for each grid cell
! author: Glenn Hammond
! date: 02/22/08
!
! ************************************************************************** !
subroutine RichardsAuxVarCompute(x,aux_var,iphase,saturation_function,option)

  use Option_module
  use water_eos_module
  use Material_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: aux_var
  PetscInt ::iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  
  aux_var%sat = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%den = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%xmol = 0.d0
  aux_var%kvr = 0.d0
  aux_var%diff = 0.d0
  kr = 0.d0
 
  aux_var%pres = x(1)  
  aux_var%temp = x(2)
 
  aux_var%pc = option%pref - aux_var%pres
  aux_var%xmol(1) = 1.d0
  if (option%nspec > 1) aux_var%xmol(2:option%nspec) = x(3:option%nspec+1)   

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
                                   option)
    dpw_dp = 0
  else
    iphase = 1
    aux_var%pc = 0.d0
    aux_var%sat = 1.d0  
    kr = 1.d0    
    pw = aux_var%pres
    dpw_dp = 1.d0
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call wateos(aux_var%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
              option%scale,ierr)

! may need to compute dpsat_dt to pass to VISW
  call psat(aux_var%temp,sat_pressure,dpsat_dt,ierr)
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call VISW(aux_var%temp,pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr) 
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif
 
  aux_var%den = dw_mol
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale
  aux_var%diff(1:option%nspec) = option%difaq
  aux_var%kvr = kr/visl
  
!  aux_var%vis = visl
!  aux_var%dvis_dp = dvis_dp
!  aux_var%kr = kr
!  aux_var%dkr_dp = dkr_dp
  aux_var%dsat_dp = ds_dp
  aux_var%dden_dt = dw_dt

  aux_var%dden_dp = dw_dp
  
  aux_var%dkvr_dt = -kr/(visl*visl)*(dvis_dt+dvis_dpsat*dpsat_dt)
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    aux_var%dh_dp = hw_dp
    aux_var%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    aux_var%dh_dp = 0.d0
    aux_var%du_dp = 0.d0
  endif
  aux_var%dh_dt = hw_dt
  aux_var%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

end subroutine RichardsAuxVarCompute

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(richards_auxvar_type) :: aux_var
  
  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
  nullify(aux_var%xmol)
  if (associated(aux_var%diff))deallocate(aux_var%diff)
  nullify(aux_var%diff)

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! RichardsAuxDestroy: Deallocates a richards auxilliary object
! author: Glenn Hammond
! date: 02/14/08
!
! ************************************************************************** !
subroutine RichardsAuxDestroy(aux)

  implicit none

  type(richards_type), pointer :: aux
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
    
end subroutine RichardsAuxDestroy

end module Richards_Aux_module
