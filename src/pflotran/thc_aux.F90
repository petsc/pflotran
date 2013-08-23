module THC_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  type, public :: thc_auxvar_type
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
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
#ifdef ICE
    PetscReal :: sat_ice
    PetscReal :: sat_gas
    PetscReal :: dsat_dt
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_gas_dp
    PetscReal :: dsat_ice_dt
    PetscReal :: dsat_gas_dt
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dt
    PetscReal :: u_ice
    PetscReal :: du_ice_dt
#endif
  end type thc_auxvar_type

  type, public :: thc_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
#ifdef ICE
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:)
#endif
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type thc_parameter_type
  
  type, public :: thc_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(thc_parameter_type), pointer :: thc_parameter
    type(thc_auxvar_type), pointer :: aux_vars(:)
    type(thc_auxvar_type), pointer :: aux_vars_bc(:)
    type(thc_auxvar_type), pointer :: aux_vars_ss(:)
  end type thc_type


  public :: THCAuxCreate, THCAuxDestroy, &
            THCAuxVarCompute, THCAuxVarInit, &
            THCAuxVarCopy

#ifdef ICE
  public :: THCAuxVarComputeIce
#endif

contains


! ************************************************************************** !
!
! THCAuxCreate: Allocate and initialize auxiliary object
! author: ???
! date: 02/14/08
!
! ************************************************************************** !
function THCAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(thc_type), pointer :: THCAuxCreate
  
  type(thc_type), pointer :: aux

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
  allocate(aux%thc_parameter)
  nullify(aux%thc_parameter%sir)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%thc_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%thc_parameter%diffusion_activation_energy(option%nphase))
  aux%thc_parameter%diffusion_coefficient = 1.d-9
  aux%thc_parameter%diffusion_activation_energy = 0.d0
 
  THCAuxCreate => aux
  
end function THCAuxCreate

! ************************************************************************** !
!
! THCAuxVarInit: Initialize auxiliary object
! author: ???
! date: 02/14/08
!
! ************************************************************************** !
subroutine THCAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(thc_auxvar_type) :: aux_var
  type(option_type) :: option
  

  aux_var%avgmw = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%pc = 0.d0
!  aux_var%kr = 0.d0
!  aux_var%dkr_dp = 0.d0
  aux_var%vis = 0.d0
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
  allocate(aux_var%xmol(option%nflowspec))
  aux_var%xmol = 0.d0
  allocate(aux_var%diff(option%nflowspec))
  aux_var%diff = 1.d-9
#ifdef ICE
  aux_var%sat_ice = 0.d0
  aux_var%sat_gas = 0.d0
  aux_var%dsat_dt = 0.d0
  aux_var%dsat_ice_dp = 0.d0
  aux_var%dsat_gas_dp = 0.d0
  aux_var%dsat_ice_dt = 0.d0
  aux_var%dsat_gas_dt = 0.d0
  aux_var%den_ice = 0.d0
  aux_var%dden_ice_dp = 0.d0
  aux_var%dden_ice_dt = 0.d0
  aux_var%u_ice = 0.d0
  aux_var%du_ice_dt = 0.d0
#endif

end subroutine THCAuxVarInit

! ************************************************************************** !
!
! THCAuxVarCopy: Copies an auxiliary variable
! author: ???
! date: 12/13/07
!
! ************************************************************************** !  
subroutine THCAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(thc_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

! aux_var2%pres = aux_var%pres
! aux_var2%temp = aux_var%temp
! aux_var2%den = aux_var%den
! aux_var2%den_kg = aux_var%den_kg
    
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
!  aux_var2%kr = aux_var%kr
!  aux_var2%dkr_dp = aux_var%dkr_dp
  aux_var2%vis = aux_var%vis
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
#ifdef ICE
  aux_var2%sat_ice = aux_var%sat_ice 
  aux_var2%sat_gas = aux_var%sat_gas
  aux_var2%dsat_dt = aux_var%dsat_dt
  aux_var2%dsat_ice_dp = aux_var%dsat_ice_dp
  aux_var2%dsat_gas_dp = aux_var%dsat_gas_dp
  aux_var2%dsat_ice_dt = aux_var%dsat_ice_dt
  aux_var2%dsat_gas_dt = aux_var%dsat_gas_dt
  aux_var2%den_ice = aux_var%den_ice
  aux_var2%dden_ice_dp = aux_var%dden_ice_dp
  aux_var2%dden_ice_dt = aux_var%dden_ice_dt
  aux_var2%u_ice = aux_var%u_ice
  aux_var2%du_ice_dt = aux_var%du_ice_dt
#endif

end subroutine THCAuxVarCopy

! ************************************************************************** !
!
! THCAuxVarCompute: Computes auxiliary variables for each grid cell
! author: ???
! date: 02/22/08
!
! ************************************************************************** !
subroutine THCAuxVarCompute(x,aux_var,global_aux_var, &
                            iphase,saturation_function,por,perm,option)

  use Option_module
  use Global_Aux_module
  use Water_EOS_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(thc_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: por, perm
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  
! aux_var%den = 0.d0
! aux_var%den_kg = 0.d0
  global_aux_var%sat = 0.d0
  global_aux_var%den = 0.d0
  global_aux_var%den_kg = 0.d0

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%xmol = 0.d0
  aux_var%kvr = 0.d0
  aux_var%diff = 0.d0
  kr = 0.d0
 
! aux_var%pres = x(1)  
! aux_var%temp = x(2)
  global_aux_var%pres = x(1)  
  global_aux_var%temp = x(2)
 
! aux_var%pc = option%reference_pressure - aux_var%pres
  aux_var%pc = option%reference_pressure - global_aux_var%pres(1)
  aux_var%xmol(1) = 1.d0
  if (option%nflowspec > 1) aux_var%xmol(2:option%nflowspec) = x(3:option%nflowspec+1)   

!***************  Liquid phase properties **************************
  aux_var%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (aux_var%pc > 0.d0) then
  if (aux_var%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(global_aux_var%pres(1),global_aux_var%sat(1), &
                                   kr,ds_dp,dkr_dp, &
                                   saturation_function, &
                                   por,perm, &
                                   option)
    dpw_dp = 0.d0
  else
    iphase = 1
    aux_var%pc = 0.d0
    global_aux_var%sat(1) = 1.d0  
    kr = 1.d0    
!   pw = aux_var%pres
    pw = global_aux_var%pres(1)
    dpw_dp = 1.d0
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call wateos(global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
              option%scale,ierr)

! may need to compute dpsat_dt to pass to VISW
  call psat(global_aux_var%temp(1),sat_pressure,dpsat_dt,ierr)
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call VISW(global_aux_var%temp(1),pw,sat_pressure,visl,dvis_dt,dvis_dp,ierr)
  
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

! aux_var%den = dw_mol
! aux_var%den_kg = dw_kg
  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale
  aux_var%kvr = kr/visl
  
  aux_var%vis = visl
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
  
end subroutine THCAuxVarCompute

! ************************************************************************** !
! 
! THCAuxVarComputeIce: Computes auxillary variables for each grid cell when
!                      ice and vapor phases are present
! author: Satish Karra
! Date: 11/16/11
!
! ************************************************************************** !

#ifdef ICE
subroutine THCAuxVarComputeIce(x, aux_var, global_aux_var, iphase, &
                               saturation_function, por, perm, option)

!sk: Not sure if we need por, perm

  use Option_module
  use Global_Aux_module
  use Water_EOS_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(thc_auxvar_type) :: aux_var
  type(global_auxvar_type) :: global_aux_var
  PetscReal :: por, perm
  PetscInt :: iphase

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: kr, ds_dp, dkr_dp, dkr_dt
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: ice_saturation, gas_saturation
  PetscReal :: dsl_temp
  PetscReal :: dsg_pl, dsg_temp
  PetscReal :: dsi_pl, dsi_temp
  PetscReal :: den_ice, dden_ice_dT, dden_ice_dP
  PetscReal :: u_ice, du_ice_dT
  PetscBool :: out_of_table_flag
  PetscReal :: p_th
  
  out_of_table_flag = PETSC_FALSE
 
  global_aux_var%sat = 0.d0
  global_aux_var%den = 0.d0
  global_aux_var%den_kg = 0.d0

  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%avgmw = 0.d0
  aux_var%xmol = 0.d0
  aux_var%kvr = 0.d0
  aux_var%diff = 0.d0
   
  global_aux_var%pres = x(1)  
  global_aux_var%temp = x(2)
  
  ! Check if the capillary pressure is less than -100MPa
  
  if (global_aux_var%pres(1) - option%reference_pressure < -1.d8 + 1.d0) then
    global_aux_var%pres(1) = -1.d8 + option%reference_pressure + 1.d0
  endif

 
  aux_var%pc = option%reference_pressure - global_aux_var%pres(1)
  aux_var%xmol(1) = 1.d0
  if (option%nflowspec > 1) aux_var%xmol(2:option%nflowspec) = x(3:option%nflowspec+1)   

!***************  Liquid phase properties **************************
  aux_var%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (aux_var%pc > 1.d0) then
    iphase = 3
    dpw_dp = 0.d0
  else
    iphase = 1
    aux_var%pc = 0.d0
    pw = global_aux_var%pres(1)
    dpw_dp = 1.d0
  endif  
  
  call CapillaryPressureThreshold(saturation_function,p_th,option)


  call SaturationFunctionComputeIce(global_aux_var%pres(1), & 
                                    global_aux_var%temp(1), ice_saturation, &
                                    global_aux_var%sat(1), gas_saturation, &
                                    kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                    dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                    saturation_function, p_th, option)


  call wateos(global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
              option%scale,ierr)

!  call wateos_flag (global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw, &
!                     hw_dp,hw_dt,option%scale,out_of_table_flag,ierr)
  
!  if (out_of_table_flag) then  
!    option%out_of_table = PETSC_TRUE                 
!  endif

!  call wateos_simple(global_aux_var%temp(1), pw, dw_kg, dw_mol, dw_dp, &
!                         dw_dt, hw, hw_dp, hw_dt, ierr)
                         
  call psat(global_aux_var%temp(1), sat_pressure, dpsat_dt, ierr)
  
  call VISW(global_aux_var%temp(1), pw, sat_pressure, visl, dvis_dt, &
            dvis_dp, ierr)

!  call VISW_temp(global_aux_var%temp(1),visl,dvis_dt,ierr)
!  dvis_dp = 0.d0
  
  dvis_dpsat = -dvis_dp 
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale
  aux_var%kvr = kr/visl
  aux_var%vis = visl
  aux_var%dsat_dp = ds_dp
  aux_var%dden_dt = dw_dt
  aux_var%dden_dp = dw_dp
  aux_var%dkvr_dt = -kr/(visl*visl)*(dvis_dt + dvis_dpsat*dpsat_dt) + dkr_dt/visl
  aux_var%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  aux_var%dh_dp = hw_dp
  aux_var%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* &
                  option%scale
  aux_var%dh_dt = hw_dt
  aux_var%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

  aux_var%sat_ice = ice_saturation
  aux_var%sat_gas = gas_saturation
  aux_var%dsat_dt = dsl_temp
  aux_var%dsat_ice_dp = dsi_pl
  aux_var%dsat_gas_dp = dsg_pl
  aux_var%dsat_ice_dt = dsi_temp
  aux_var%dsat_gas_dt = dsg_temp
  
! Calculate the density, internal energy and derivatives for ice
  call DensityIce(global_aux_var%temp(1), global_aux_var%pres(1), &
                  den_ice, dden_ice_dT, dden_ice_dP)

  call InternalEnergyIce(global_aux_var%temp(1), u_ice, du_ice_dT)

  aux_var%den_ice = den_ice
  aux_var%dden_ice_dt = dden_ice_dT
  aux_var%dden_ice_dp = dden_ice_dP
  aux_var%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  aux_var%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 

end subroutine THCAuxVarComputeIce
#endif

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a thc auxiliary object
! author: ???
! date: 02/14/08
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(thc_auxvar_type) :: aux_var
  
  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
  nullify(aux_var%xmol)
  if (associated(aux_var%diff))deallocate(aux_var%diff)
  nullify(aux_var%diff)

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! THCAuxDestroy: Deallocates a thc auxiliary object
! author: ???
! date: 02/14/08
!
! ************************************************************************** !
subroutine THCAuxDestroy(aux)

  implicit none

  type(thc_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call AuxVarDestroy(aux%aux_vars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call AuxVarDestroy(aux%aux_vars_bc(iaux))
  enddo  
  do iaux = 1, aux%num_aux_ss
    call AuxVarDestroy(aux%aux_vars_ss(iaux))
  enddo  
  
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
  if (associated(aux%thc_parameter)) then
    if (associated(aux%thc_parameter%diffusion_coefficient)) &
      deallocate(aux%thc_parameter%diffusion_coefficient)
    nullify(aux%thc_parameter%diffusion_coefficient)
    if (associated(aux%thc_parameter%diffusion_activation_energy)) &
      deallocate(aux%thc_parameter%diffusion_activation_energy)
    nullify(aux%thc_parameter%diffusion_activation_energy)
    if (associated(aux%thc_parameter%dencpr)) deallocate(aux%thc_parameter%dencpr)
    nullify(aux%thc_parameter%dencpr)
    if (associated(aux%thc_parameter%ckwet)) deallocate(aux%thc_parameter%ckwet)
    nullify(aux%thc_parameter%ckwet)
    if (associated(aux%thc_parameter%ckdry)) deallocate(aux%thc_parameter%ckdry)
    nullify(aux%thc_parameter%ckdry)
    if (associated(aux%thc_parameter%alpha)) deallocate(aux%thc_parameter%alpha)
    nullify(aux%thc_parameter%alpha)
#ifdef ICE
    if (associated(aux%thc_parameter%ckfrozen)) deallocate(aux%thc_parameter%ckfrozen)
    nullify(aux%thc_parameter%ckfrozen)
    if (associated(aux%thc_parameter%alpha_fr)) deallocate(aux%thc_parameter%alpha_fr)
    nullify(aux%thc_parameter%alpha_fr)
#endif
    if (associated(aux%thc_parameter%sir)) deallocate(aux%thc_parameter%sir)
    nullify(aux%thc_parameter%sir)
  endif
  nullify(aux%thc_parameter)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine THCAuxDestroy

end module THC_Aux_module
