module THMC_Aux_module

  implicit none
  
  private 

#include "definitions.h"

  type, public :: thmc_auxvar_type
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
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
    PetscReal :: stress(3,3)
    PetscReal :: gradient(3,3)
    PetscReal :: Minv(3,3)           ! this matrix depends only on the grid and is precomputed and stored
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
  end type thmc_auxvar_type

  type, public :: thmc_parameter_type
    PetscReal, pointer :: rock_den(:)
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: youngs_modulus(:) !Elasticity parameters
    PetscReal, pointer :: poissons_ratio(:)
#ifdef ICE
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:)
#endif
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type thmc_parameter_type
  
  type, public :: thmc_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)

    PetscBool :: aux_vars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc
    type(thmc_parameter_type), pointer :: thmc_parameter
    type(thmc_auxvar_type), pointer :: aux_vars(:)
    type(thmc_auxvar_type), pointer :: aux_vars_bc(:)
  end type thmc_type

  public :: THMCAuxCreate, THMCAuxDestroy, &
            THMCAuxVarCompute, THMCAuxVarInit, &
            THMCAuxVarCopy

#ifdef ICE
  public :: THMCAuxVarComputeIce
#endif

contains


! ************************************************************************** !
!
! THMCAuxCreate: Allocate and initialize auxiliary object
! author:
! date: 3/2/12
!
! ************************************************************************** !
function THMCAuxCreate(option)

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(thmc_type), pointer :: THMCAuxCreate
  
  type(thmc_type), pointer :: aux

  allocate(aux) 
  aux%aux_vars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  nullify(aux%aux_vars)
  nullify(aux%aux_vars_bc)
  aux%n_zero_rows = 0
  allocate(aux%thmc_parameter)
  nullify(aux%thmc_parameter%sir)
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%thmc_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%thmc_parameter%diffusion_activation_energy(option%nphase))
  aux%thmc_parameter%diffusion_coefficient = 1.d-9
  aux%thmc_parameter%diffusion_activation_energy = 0.d0
 ! aux%thmc_parameter%youngs_modulus = 2.d11 ! in Pa (check this)
 ! aux%thmc_parameter%poissons_ratio = 0.3

  THMCAuxCreate => aux
  
end function THMCAuxCreate

! ************************************************************************** !
!
! THMCAuxVarInit: Initialize auxiliary object
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCAuxVarInit(aux_var,option)

  use Option_module

  implicit none
  
  type(thmc_auxvar_type) :: aux_var
  type(option_type) :: option
  
  aux_var%avgmw = 0.d0
  aux_var%h = 0.d0
  aux_var%u = 0.d0
  aux_var%pc = 0.d0
  aux_var%vis = 0.d0
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
  aux_var%gradient = 0.d0
  aux_var%stress = 0.d0
  aux_var%Minv = 0.d0
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

end subroutine THMCAuxVarInit

! ************************************************************************** !
!
! THMCAuxVarCopy: Copies an auxiliary variable
! author:
! date: 3/2/12
!
! ************************************************************************** !  
subroutine THMCAuxVarCopy(aux_var,aux_var2,option)

  use Option_module

  implicit none
  
  type(thmc_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option
  
  aux_var2%avgmw = aux_var%avgmw
  aux_var2%h = aux_var%h
  aux_var2%u = aux_var%u
  aux_var2%pc = aux_var%pc
  aux_var2%vis = aux_var%vis
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
  aux_var2%gradient = aux_var%gradient
  aux_var2%stress = aux_var%stress
  aux_var2%Minv = aux_var%Minv
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

end subroutine THMCAuxVarCopy

! ************************************************************************** !
!
! THMCAuxVarCompute: Computes auxiliary variables for each grid cell
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCAuxVarCompute(x,aux_var,global_aux_var, &
                            iphase,saturation_function,por,perm,option)

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(thmc_auxvar_type) :: aux_var
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
 
  global_aux_var%pres = x(1)  
  global_aux_var%temp = x(2)
 
  aux_var%pc = option%reference_pressure - global_aux_var%pres(1)
  aux_var%xmol(1) = 1.d0
  if (option%nflowspec > 1) aux_var%xmol(2:option%nflowspec) = x(3:option%nflowspec+1)   

  global_aux_var%displacement(:) = x(option%nflowdof-option%nmechdof+1:option%nflowdof)

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
    pw = global_aux_var%pres(1)
    dpw_dp = 1.d0
  endif  

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

  global_aux_var%den = dw_mol
  global_aux_var%den_kg = dw_kg
  
  aux_var%h = hw
  aux_var%u = aux_var%h - pw / dw_mol * option%scale
  aux_var%kvr = kr/visl
  
  aux_var%vis = visl
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
  
end subroutine THMCAuxVarCompute

! ************************************************************************** !
! 
! THMCAuxVarComputeIce: Computes auxillary variables for each grid cell when
!                      ice and vapor phases are present
! author:
! date: 3/2/12
!
! ************************************************************************** !

#ifdef ICE
subroutine THMCAuxVarComputeIce(x, aux_var, global_aux_var, iphase, &
                               saturation_function, por, perm, option)

!sk: Not sure if we need por, perm

  use Option_module
  use Global_Aux_module
  use water_eos_module
  use Saturation_Function_module  
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(thmc_auxvar_type) :: aux_var
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
  PetscReal :: u_ice, du_ice_dT, pth
 
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
  
  pth = 1.d8 - 1

  call SaturationFunctionComputeIce(global_aux_var%pres(1), & 
                                    global_aux_var%temp(1), ice_saturation, &
                                    global_aux_var%sat(1), gas_saturation, &
                                    kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                    dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                    saturation_function, pth, option)


!  call wateos(global_aux_var%temp(1),pw,dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt, &
!              option%scale,ierr)
              
!  print *, 'wateos:', dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt

  call wateos_simple(global_aux_var%temp(1), pw, dw_kg, dw_mol, dw_dp, &
                         dw_dt, hw, hw_dp, hw_dt, ierr)
                         
!  print *, 'wateos_simple', dw_kg,dw_mol,dw_dp,dw_dt,hw,hw_dp,hw_dt

   

  call psat(global_aux_var%temp(1), sat_pressure, dpsat_dt, ierr)
  
  call VISW(global_aux_var%temp(1), pw, sat_pressure, visl, dvis_dt, &
            dvis_dp, ierr)
  
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

end subroutine THMCAuxVarComputeIce
#endif

! ************************************************************************** !
!
! AuxVarDestroy: Deallocates a thmc auxiliary object
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine AuxVarDestroy(aux_var)

  implicit none

  type(thmc_auxvar_type) :: aux_var
  
  if (associated(aux_var%xmol)) deallocate(aux_var%xmol)
  nullify(aux_var%xmol)
  if (associated(aux_var%diff))deallocate(aux_var%diff)
  nullify(aux_var%diff)

end subroutine AuxVarDestroy

! ************************************************************************** !
!
! THMCAuxDestroy: Deallocates a thmc auxiliary object
! author:
! date: 3/2/12
!
! ************************************************************************** !
subroutine THMCAuxDestroy(aux)

  implicit none

  type(thmc_type), pointer :: aux
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
  if (associated(aux%thmc_parameter)) then
    if (associated(aux%thmc_parameter%diffusion_coefficient)) &
      deallocate(aux%thmc_parameter%diffusion_coefficient)
    nullify(aux%thmc_parameter%diffusion_coefficient)
    if (associated(aux%thmc_parameter%diffusion_activation_energy)) &
      deallocate(aux%thmc_parameter%diffusion_activation_energy)
    nullify(aux%thmc_parameter%diffusion_activation_energy)
    if (associated(aux%thmc_parameter%dencpr)) deallocate(aux%thmc_parameter%dencpr)
    nullify(aux%thmc_parameter%dencpr)
    if (associated(aux%thmc_parameter%rock_den)) deallocate(aux%thmc_parameter%rock_den)
    nullify(aux%thmc_parameter%rock_den)
    if (associated(aux%thmc_parameter%ckwet)) deallocate(aux%thmc_parameter%ckwet)
    nullify(aux%thmc_parameter%ckwet)
    if (associated(aux%thmc_parameter%ckdry)) deallocate(aux%thmc_parameter%ckdry)
    nullify(aux%thmc_parameter%ckdry)
    if (associated(aux%thmc_parameter%alpha)) deallocate(aux%thmc_parameter%alpha)
    nullify(aux%thmc_parameter%alpha)
    if (associated(aux%thmc_parameter%youngs_modulus)) deallocate(aux%thmc_parameter%youngs_modulus)
    nullify(aux%thmc_parameter%youngs_modulus)
    if (associated(aux%thmc_parameter%poissons_ratio)) deallocate(aux%thmc_parameter%poissons_ratio)
    nullify(aux%thmc_parameter%poissons_ratio)
#ifdef ICE
    if (associated(aux%thmc_parameter%ckfrozen)) deallocate(aux%thmc_parameter%ckfrozen)
    nullify(aux%thmc_parameter%ckfrozen)
    if (associated(aux%thmc_parameter%alpha_fr)) deallocate(aux%thmc_parameter%alpha_fr)
    nullify(aux%thmc_parameter%alpha_fr)
#endif
    if (associated(aux%thmc_parameter%sir)) deallocate(aux%thmc_parameter%sir)
    nullify(aux%thmc_parameter%sir)
  endif
  nullify(aux%thmc_parameter)

  deallocate(aux)
  nullify(aux)  

  end subroutine THMCAuxDestroy

end module THMC_Aux_module
