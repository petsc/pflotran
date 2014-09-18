module TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "finclude/petscsys.h"

  type, public :: TH_auxvar_type
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
    PetscReal :: transient_por
    PetscReal :: Dk_eff
    PetscReal :: Ke
    PetscReal :: Ke_fr
    PetscReal :: dKe_dp
    PetscReal :: dKe_dt
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dt
    ! ice
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
    PetscReal :: den_gas
    PetscReal :: dden_gas_dt
    PetscReal :: u_gas
    PetscReal :: du_gas_dt
    PetscReal :: mol_gas
    PetscReal :: dmol_gas_dt
    ! For DallAmico model
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dt
    ! For surface-flow
    PetscBool :: surf_wat
    PetscReal :: P_min
    PetscReal :: P_max
    PetscReal :: coeff_for_cubic_approx(4)
    PetscReal :: coeff_for_deriv_cubic_approx(4)
    PetscReal :: range_for_linear_approx(4)
    PetscReal :: dlinear_slope_dT
    PetscBool :: bcflux_default_scheme

  end type TH_auxvar_type

  type, public :: TH_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:) ! exponent frozen
    PetscReal, pointer :: sir(:,:)
    PetscReal, pointer :: diffusion_coefficient(:)
    PetscReal, pointer :: diffusion_activation_energy(:)
  end type TH_parameter_type
  
  type, public :: TH_type
    PetscInt :: n_zero_rows
    PetscInt, pointer :: zero_rows_local(:), zero_rows_local_ghosted(:)
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(TH_parameter_type), pointer :: TH_parameter
    type(TH_auxvar_type), pointer :: auxvars(:)
    type(TH_auxvar_type), pointer :: auxvars_bc(:)
    type(TH_auxvar_type), pointer :: auxvars_ss(:)
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6

  public :: THAuxCreate, THAuxDestroy, &
            THAuxVarComputeNoFreezing, THAuxVarInit, &
            THAuxVarCopy

  public :: THAuxVarComputeFreezing

contains

! ************************************************************************** !

function THAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(TH_type), pointer :: THAuxCreate
  
  type(TH_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  aux%n_zero_rows = 0

  allocate(aux%TH_parameter)
  nullify(aux%TH_parameter%dencpr)
  nullify(aux%TH_parameter%ckdry)
  nullify(aux%TH_parameter%ckwet)
  nullify(aux%TH_parameter%alpha)
  nullify(aux%TH_parameter%ckfrozen)
  nullify(aux%TH_parameter%alpha_fr)
  nullify(aux%TH_parameter%sir)
  nullify(aux%TH_parameter%diffusion_coefficient)
  nullify(aux%TH_parameter%diffusion_activation_energy)
  
  nullify(aux%zero_rows_local)
  nullify(aux%zero_rows_local_ghosted)

  allocate(aux%TH_parameter%diffusion_coefficient(option%nphase))
  allocate(aux%TH_parameter%diffusion_activation_energy(option%nphase))
  aux%TH_parameter%diffusion_coefficient = 1.d-9
  aux%TH_parameter%diffusion_activation_energy = 0.d0
 
  THAuxCreate => aux
  
end function THAuxCreate

! ************************************************************************** !

subroutine THAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Option_module

  implicit none
  
  type(TH_auxvar_type) :: auxvar
  type(option_type) :: option
  

  auxvar%avgmw = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%pc = 0.d0
!  auxvar%kr = 0.d0
!  auxvar%dkr_dp = 0.d0
  auxvar%vis = 0.d0
!  auxvar%dvis_dp = 0.d0
  auxvar%kvr = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dden_dp = 0.d0
  auxvar%dden_dt = 0.d0
  auxvar%dkvr_dp = 0.d0
  auxvar%dkvr_dt = 0.d0
  auxvar%dh_dp = 0.d0
  auxvar%dh_dt = 0.d0
  auxvar%du_dp = 0.d0
  auxvar%du_dt = 0.d0    
  auxvar%transient_por = 0.d0
  auxvar%Dk_eff = 0.d0
  auxvar%Ke = 0.d0
  auxvar%Ke_fr = 0.d0
  auxvar%dKe_dp = 0.d0
  auxvar%dKe_dt = 0.d0
  auxvar%dKe_fr_dp = 0.d0
  auxvar%dKe_fr_dt = 0.d0
  allocate(auxvar%xmol(option%nflowspec))
  auxvar%xmol = 0.d0
  allocate(auxvar%diff(option%nflowspec))
  auxvar%diff = 1.d-9
  ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if not used!
  auxvar%sat_ice = 0.d0
  auxvar%sat_gas = 0.d0
  auxvar%dsat_dt = 0.d0
  auxvar%dsat_ice_dp = 0.d0
  auxvar%dsat_gas_dp = 0.d0
  auxvar%dsat_ice_dt = 0.d0
  auxvar%dsat_gas_dt = 0.d0
  auxvar%den_ice = 0.d0
  auxvar%dden_ice_dp = 0.d0
  auxvar%dden_ice_dt = 0.d0
  auxvar%u_ice = 0.d0
  auxvar%du_ice_dt = 0.d0
  auxvar%den_gas = 0.d0
  auxvar%dden_gas_dt = 0.d0
  auxvar%u_gas = 0.d0
  auxvar%du_gas_dt = 0.d0
  auxvar%mol_gas = 0.d0
  auxvar%dmol_gas_dt = 0.d0
  auxvar%pres_fh2o = 0.d0
  auxvar%dpres_fh2o_dp = 0.d0
  auxvar%dpres_fh2o_dt = 0.d0
  auxvar%surf_wat = PETSC_FALSE
  auxvar%P_min = 0.d0
  auxvar%P_max = 0.d0
  auxvar%coeff_for_cubic_approx(:) = 0.d0
  auxvar%coeff_for_deriv_cubic_approx(:) = 0.d0
  auxvar%range_for_linear_approx(:) = 0.d0
  auxvar%dlinear_slope_dT = 0.d0
  auxvar%bcflux_default_scheme = PETSC_FALSE

end subroutine THAuxVarInit

! ************************************************************************** !

subroutine THAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module

  implicit none
  
  type(TH_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

! auxvar2%pres = auxvar%pres
! auxvar2%temp = auxvar%temp
! auxvar2%den = auxvar%den
! auxvar2%den_kg = auxvar%den_kg
    
  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
!  auxvar2%kr = auxvar%kr
!  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%vis = auxvar%vis
!  auxvar2%dvis_dp = auxvar%dvis_dp
  auxvar2%kvr = auxvar%kvr
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dkvr_dp = auxvar%dkvr_dp
  auxvar2%dkvr_dt = auxvar%dkvr_dt
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dt = auxvar%dh_dt
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dt = auxvar%du_dt  
  auxvar2%transient_por = auxvar%transient_por
  auxvar2%Dk_eff = auxvar%Dk_eff
  auxvar2%Ke = auxvar%Ke
  auxvar2%Ke_fr = auxvar%Ke_fr
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_fr_dp = auxvar%dKe_fr_dp
  auxvar2%dKe_dt = auxvar%dKe_dt
  auxvar2%dKe_fr_dt = auxvar%dKe_fr_dt
  auxvar2%xmol = auxvar%xmol
  auxvar2%diff = auxvar%diff
  if (option%use_th_freezing) then
     auxvar2%sat_ice = auxvar%sat_ice 
     auxvar2%sat_gas = auxvar%sat_gas
     auxvar2%dsat_dt = auxvar%dsat_dt
     auxvar2%dsat_ice_dp = auxvar%dsat_ice_dp
     auxvar2%dsat_gas_dp = auxvar%dsat_gas_dp
     auxvar2%dsat_ice_dt = auxvar%dsat_ice_dt
     auxvar2%dsat_gas_dt = auxvar%dsat_gas_dt
     auxvar2%den_ice = auxvar%den_ice
     auxvar2%dden_ice_dp = auxvar%dden_ice_dp
     auxvar2%dden_ice_dt = auxvar%dden_ice_dt
     auxvar2%u_ice = auxvar%u_ice
     auxvar2%du_ice_dt = auxvar%du_ice_dt
     auxvar2%pres_fh2o = auxvar%pres_fh2o
     auxvar2%dpres_fh2o_dp = auxvar%dpres_fh2o_dp
     auxvar2%dpres_fh2o_dt = auxvar%dpres_fh2o_dt
     auxvar2%den_gas = auxvar%den_gas
     auxvar2%dden_gas_dt = auxvar%dden_gas_dt
     auxvar2%u_gas = auxvar%u_gas
     auxvar2%du_gas_dt = auxvar%du_gas_dt
     auxvar2%mol_gas = auxvar%mol_gas
     auxvar2%dmol_gas_dt = auxvar%dmol_gas_dt
  endif
  auxvar2%surf_wat = auxvar%surf_wat
  auxvar2%P_min = auxvar%P_min
  auxvar2%P_max = auxvar%P_max
  auxvar2%coeff_for_cubic_approx(:) = auxvar%coeff_for_cubic_approx(:)
  auxvar2%coeff_for_deriv_cubic_approx(:) = auxvar%coeff_for_deriv_cubic_approx(:)
  auxvar2%range_for_linear_approx(:) = auxvar%range_for_linear_approx(:)
  auxvar2%dlinear_slope_dT = auxvar%dlinear_slope_dT
  auxvar2%bcflux_default_scheme = auxvar%bcflux_default_scheme

end subroutine THAuxVarCopy

! ************************************************************************** !

subroutine THAuxVarComputeNoFreezing(x,auxvar,global_auxvar, &
                                     material_auxvar, &
                                     iphase,saturation_function, &
                                     th_parameter, ithrm, &
                                     option)
  ! 
  ! Computes auxiliary variables for each grid cell
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module  
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  class(material_auxvar_type) :: material_auxvar

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp, dvis_dpsat
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry

! auxvar%den = 0.d0
! auxvar%den_kg = 0.d0
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%xmol = 0.d0
  auxvar%kvr = 0.d0
  auxvar%diff = 0.d0
  kr = 0.d0
 
! auxvar%pres = x(1)  
! auxvar%temp = x(2)
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
 
! auxvar%pc = option%reference_pressure - auxvar%pres
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)
  auxvar%xmol(1) = 1.d0
  if (option%nflowspec > 1) auxvar%xmol(2:option%nflowspec) = x(3:option%nflowspec+1)   

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
!  if (auxvar%pc > 0.d0) then
  if (auxvar%pc > 1.d0) then
    iphase = 3
    call SaturationFunctionCompute(auxvar%pc,global_auxvar%sat(1), &
                                   kr,ds_dp,dkr_dp, &
                                   saturation_function, &
                                   material_auxvar%porosity, &
                                   material_auxvar%permeability(perm_xx_index), &
                                   option)
    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    global_auxvar%sat(1) = 1.d0  
    kr = 1.d0    
!   pw = auxvar%pres
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  

!  call wateos_noderiv(option%temp,pw,dw_kg,dw_mol,hw,option%scale,ierr)
  call EOSWaterDensityEnthalpy(global_auxvar%temp,pw,dw_kg,dw_mol,hw, &
                               dw_dp,dw_dt,hw_dp,hw_dt,ierr)
  ! J/kmol -> whatever units
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  
! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,dpsat_dt,ierr)
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,dpsat_dt,visl, &
                         dvis_dt,dvis_dp,dvis_dpsat,ierr)
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

! auxvar%den = dw_mol
! auxvar%den_kg = dw_kg
  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  
  auxvar%vis = visl
!  auxvar%dvis_dp = dvis_dp
!  auxvar%kr = kr
!  auxvar%dkr_dp = dkr_dp
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dt = dw_dt

  auxvar%dden_dp = dw_dp
  
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt+dvis_dpsat*dpsat_dt)
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    auxvar%dh_dp = hw_dp
    auxvar%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    auxvar%dh_dp = 0.d0
    auxvar%du_dp = 0.d0
  endif

  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt
  
  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)

  !unfrozen soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  auxvar%Ke = Ke

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk_dry + (Dk - Dk_dry)*Ke

  ! Derivative of soil Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dp
  auxvar%dKe_dt = 0.d0

end subroutine THAuxVarComputeNoFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing(x, auxvar, global_auxvar, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   th_parameter, ithrm, &
                                   option)
  ! 
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  ! 

!sk: Not sure if we need por, perm

  use Option_module
  use Global_Aux_module
  
  use EOS_Water_module
  use Saturation_Function_module  
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  PetscReal :: x(option%nflowdof)
  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
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

  PetscReal :: p_g
  PetscReal :: p_sat
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dt
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K

  PetscReal :: Ke
  PetscReal :: Ke_fr
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice

  out_of_table_flag = PETSC_FALSE
 
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%xmol = 0.d0
  auxvar%kvr = 0.d0
  auxvar%diff = 0.d0
   
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
  
  ! Check if the capillary pressure is less than -100MPa
  
  if (global_auxvar%pres(1) - option%reference_pressure < -1.d8 + 1.d0) then
    global_auxvar%pres(1) = -1.d8 + option%reference_pressure + 1.d0
  endif

 
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)
  auxvar%xmol(1) = 1.d0
  if (option%nflowspec > 1) auxvar%xmol(2:option%nflowspec) = x(3:option%nflowspec+1)   

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 1.d0) then
    iphase = 3
    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    pw = global_auxvar%pres(1)
    dpw_dp = 1.d0
  endif  
  
  call CapillaryPressureThreshold(saturation_function,p_th,option)

  select case (option%ice_model)
    case (PAINTER_EXPLICIT)
      ! Model from Painter, Comp. Geosci. (2011)
      call SatFuncComputeIcePExplicit(global_auxvar%pres(1), & 
                                      global_auxvar%temp, ice_saturation, &
                                      global_auxvar%sat(1), gas_saturation, &
                                      kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                      saturation_function, p_th, option)    
    case (PAINTER_KARRA_IMPLICIT)
      ! Implicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKImplicit(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option)    
    case (PAINTER_KARRA_EXPLICIT)
      ! Explicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKExplicit(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option) 
    case (DALL_AMICO)
      ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
      call SatFuncComputeIceDallAmico(global_auxvar%pres(1), &
                                      global_auxvar%temp, &
                                      auxvar%pres_fh2o, &
                                      auxvar%dpres_fh2o_dp, &
                                      auxvar%dpres_fh2o_dt, &
                                      ice_saturation, &
                                      global_auxvar%sat(1), gas_saturation, &
                                      kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                      saturation_function, option)
    case (PAINTER_KARRA_EXPLICIT_NOCRYO)
      ! Explicit model from Painter & Karra, VJZ (2013) and removed cryosuction
      call SatFuncComputeIcePKExplicitNoCryo(global_auxvar%pres(1), & 
                                       global_auxvar%temp, ice_saturation, &
                                       global_auxvar%sat(1), gas_saturation, &
                                       kr, ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dt, &
                                       saturation_function, p_th, option) 
    case default
      option%io_buffer = 'THCAuxVarComputeIce: Ice model not recognized.'
      call printErrMsg(option)
  end select

  call EOSWaterDensityEnthalpy(global_auxvar%temp,pw,dw_kg,dw_mol,hw, &
                               dw_dp,dw_dt,hw_dp,hw_dt,ierr)
  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
                         
  call EOSWaterSaturationPressure(global_auxvar%temp, sat_pressure, &
                                  dpsat_dt, ierr)
  call EOSWaterViscosity(global_auxvar%temp, pw, sat_pressure, dpsat_dt, &
                         visl, dvis_dt,dvis_dp, dvis_dpsat, ierr)

  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  global_auxvar%den = dw_mol
  global_auxvar%den_kg = dw_kg
  
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%kvr = kr/visl
  auxvar%vis = visl
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dt = dw_dt
  auxvar%dden_dp = dw_dp
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity  
!  auxvar%dkvr_dt = -kr/(visl*visl)*(dvis_dt + dvis_dpsat*dpsat_dt) + dkr_dt/visl
  auxvar%dkvr_dt = -kr/(visl*visl)*dvis_dt + dkr_dt/visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp
  auxvar%dh_dp = hw_dp
  auxvar%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* &
                  option%scale
  auxvar%dh_dt = hw_dt
  auxvar%du_dt = hw_dt + pw/(dw_mol*dw_mol)*option%scale*dw_dt

  auxvar%sat_ice = ice_saturation
  auxvar%sat_gas = gas_saturation
  auxvar%dsat_dt = dsl_temp
  auxvar%dsat_ice_dp = dsi_pl
  auxvar%dsat_gas_dp = dsg_pl
  auxvar%dsat_ice_dt = dsi_temp
  auxvar%dsat_gas_dt = dsg_temp
  
  ! Calculate the density, internal energy and derivatives for ice
  call EOSWaterDensityIce(global_auxvar%temp, global_auxvar%pres(1), &
                          den_ice, dden_ice_dT, dden_ice_dP, ierr)

  call EOSWaterInternalEnergyIce(global_auxvar%temp, u_ice, du_ice_dT)

  auxvar%den_ice = den_ice
  auxvar%dden_ice_dt = dden_ice_dT
  auxvar%dden_ice_dp = dden_ice_dP
  auxvar%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  auxvar%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 

  ! Calculate the values and derivatives for density and internal energy
  call EOSWaterSaturationPressure(global_auxvar%temp, p_sat, ierr)

  p_g            = option%reference_pressure
  auxvar%den_gas = p_g/(IDEAL_GAS_CONST*(global_auxvar%temp + 273.15d0))*1.d-3 !in kmol/m3
  mol_g          = p_sat/p_g
  C_g            = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  auxvar%u_gas   = C_g*(global_auxvar%temp + 273.15d0)           ! in MJ/kmol
  auxvar%mol_gas = mol_g

  auxvar%dden_gas_dt = - p_g/(IDEAL_GAS_CONST*(global_auxvar%temp + 273.15d0)**2)*1.d-3
  dmolg_dt           = dpsat_dt/p_g
  auxvar%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)* &
                       (global_auxvar%temp + 273.15d0)
  auxvar%dmol_gas_dt = dmolg_dt

  if (option%ice_model == DALL_AMICO) then
    auxvar%den_ice = dw_mol
    auxvar%dden_ice_dt = auxvar%dden_dt
    auxvar%dden_ice_dp = auxvar%dden_dp
!    auxvar%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%du_ice_dt = auxvar%du_dt
  endif

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dp
  auxvar%dKe_dt = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)*auxvar%dsat_dt
  auxvar%dKe_fr_dt = alpha_fr*(global_auxvar%sat(1) + epsilon)**(alpha_fr - 1.d0)*auxvar%dsat_dt
  auxvar%dKe_fr_dp = alpha_fr*(global_auxvar%sat(1) + epsilon)**(alpha_fr - 1.d0)*auxvar%dsat_dp

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_auxvar_type) :: auxvar
  
  if (associated(auxvar%xmol)) deallocate(auxvar%xmol)
  nullify(auxvar%xmol)
  if (associated(auxvar%diff))deallocate(auxvar%diff)
  nullify(auxvar%diff)

end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine THAuxDestroy(aux)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_type), pointer :: aux
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  do iaux = 1, aux%num_aux
    call AuxVarDestroy(aux%auxvars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call AuxVarDestroy(aux%auxvars_bc(iaux))
  enddo  
  do iaux = 1, aux%num_aux_ss
    call AuxVarDestroy(aux%auxvars_ss(iaux))
  enddo  
  
  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) deallocate(aux%auxvars_ss)
  nullify(aux%auxvars_ss)
  if (associated(aux%zero_rows_local)) deallocate(aux%zero_rows_local)
  nullify(aux%zero_rows_local)
  if (associated(aux%zero_rows_local_ghosted)) deallocate(aux%zero_rows_local_ghosted)
  nullify(aux%zero_rows_local_ghosted)
  if (associated(aux%TH_parameter)) then
    if (associated(aux%TH_parameter%diffusion_coefficient)) &
      deallocate(aux%TH_parameter%diffusion_coefficient)
    nullify(aux%TH_parameter%diffusion_coefficient)
    if (associated(aux%TH_parameter%diffusion_activation_energy)) &
      deallocate(aux%TH_parameter%diffusion_activation_energy)
    nullify(aux%TH_parameter%diffusion_activation_energy)
    if (associated(aux%TH_parameter%dencpr)) deallocate(aux%TH_parameter%dencpr)
    nullify(aux%TH_parameter%dencpr)
    if (associated(aux%TH_parameter%ckwet)) deallocate(aux%TH_parameter%ckwet)
    nullify(aux%TH_parameter%ckwet)
    if (associated(aux%TH_parameter%ckdry)) deallocate(aux%TH_parameter%ckdry)
    nullify(aux%TH_parameter%ckdry)
    if (associated(aux%TH_parameter%alpha)) deallocate(aux%TH_parameter%alpha)
    nullify(aux%TH_parameter%alpha)
    ! ice
    if (associated(aux%TH_parameter%ckfrozen)) deallocate(aux%TH_parameter%ckfrozen)
    nullify(aux%TH_parameter%ckfrozen)
    if (associated(aux%TH_parameter%alpha_fr)) deallocate(aux%TH_parameter%alpha_fr)
    nullify(aux%TH_parameter%alpha_fr)

    if (associated(aux%TH_parameter%sir)) deallocate(aux%TH_parameter%sir)
    nullify(aux%TH_parameter%sir)
  endif
  nullify(aux%TH_parameter)
  
  deallocate(aux)
  nullify(aux)  

  end subroutine THAuxDestroy

end module TH_Aux_module
