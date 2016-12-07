module TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: th_itol_scaled_res = 1.d-5
  PetscReal, public :: th_itol_rel_update = UNINITIALIZED_DOUBLE

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
    PetscReal :: dsat_dt
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dkvr_dp
    PetscReal :: dkvr_dt
    PetscReal :: dh_dp
    PetscReal :: dh_dt
    PetscReal :: du_dp
    PetscReal :: du_dt
    PetscReal :: transient_por
    PetscReal :: Dk_eff
    PetscReal :: Ke
    PetscReal :: dKe_dp
    PetscReal :: dKe_dt
    ! for ice
    type(th_ice_type), pointer :: ice
    ! For surface-flow
    type(th_surface_flow_type), pointer :: surface
  end type TH_auxvar_type

  type, public :: th_ice_type
    PetscReal :: Ke_fr
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dt
    ! ice
    PetscReal :: sat_ice
    PetscReal :: sat_gas
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
  end type th_ice_type
  
  type, public :: th_surface_flow_type
    PetscBool :: surf_wat
    PetscReal :: P_min
    PetscReal :: P_max
    PetscReal :: coeff_for_cubic_approx(4)
    PetscReal :: coeff_for_deriv_cubic_approx(4)
    PetscReal :: range_for_linear_approx(4)
    PetscReal :: dlinear_slope_dT
    PetscBool :: bcflux_default_scheme
  end type th_surface_flow_type

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
            THAuxVarCopy, THAuxVarDestroy

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
#include "petsc/finclude/petscsys.h"
  use petscsys

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
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none
  
  type(TH_auxvar_type) :: auxvar
  type(option_type) :: option
  
  PetscReal :: uninit_value
  uninit_value     = UNINITIALIZED_DOUBLE

  auxvar%avgmw     = uninit_value
  auxvar%h         = uninit_value
  auxvar%u         = uninit_value
  auxvar%pc        = uninit_value
  !auxvar%kr       = uninit_value
  !auxvar%dkr_dp   = uninit_value
  auxvar%vis       = uninit_value
  !auxvar%dvis_dp  = uninit_value
  auxvar%kvr       = uninit_value
  auxvar%dsat_dp   = uninit_value
  auxvar%dsat_dt   = uninit_value
  auxvar%dden_dp   = uninit_value
  auxvar%dden_dt   = uninit_value
  auxvar%dkvr_dp   = uninit_value
  auxvar%dkvr_dt   = uninit_value
  auxvar%dh_dp     = uninit_value
  auxvar%dh_dt     = uninit_value
  auxvar%du_dp     = uninit_value
  auxvar%du_dt     = uninit_value    
  auxvar%transient_por = uninit_value
  auxvar%Dk_eff    = uninit_value
  auxvar%Ke        = uninit_value
  auxvar%dKe_dp    = uninit_value
  auxvar%dKe_dt    = uninit_value
  if (option%use_th_freezing) then
    allocate(auxvar%ice)
    auxvar%ice%Ke_fr     = uninit_value
    auxvar%ice%dKe_fr_dp = uninit_value
    auxvar%ice%dKe_fr_dt = uninit_value
    ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if not used!
    auxvar%ice%sat_ice       = uninit_value
    auxvar%ice%sat_gas       = uninit_value
    auxvar%ice%dsat_ice_dp   = uninit_value
    auxvar%ice%dsat_gas_dp   = uninit_value
    auxvar%ice%dsat_ice_dt   = uninit_value
    auxvar%ice%dsat_gas_dt   = uninit_value
    auxvar%ice%den_ice       = uninit_value
    auxvar%ice%dden_ice_dp   = uninit_value
    auxvar%ice%dden_ice_dt   = uninit_value
    auxvar%ice%u_ice         = uninit_value
    auxvar%ice%du_ice_dt     = uninit_value
    auxvar%ice%den_gas       = uninit_value
    auxvar%ice%dden_gas_dt   = uninit_value
    auxvar%ice%u_gas         = uninit_value
    auxvar%ice%du_gas_dt     = uninit_value
    auxvar%ice%mol_gas       = uninit_value
    auxvar%ice%dmol_gas_dt   = uninit_value
    auxvar%ice%pres_fh2o     = uninit_value
    auxvar%ice%dpres_fh2o_dp = uninit_value
    auxvar%ice%dpres_fh2o_dt = uninit_value
  else
    nullify(auxvar%ice)
  endif
  if (option%surf_flow_on) then
    allocate(auxvar%surface)
    auxvar%surface%surf_wat      = PETSC_FALSE
    auxvar%surface%P_min         = uninit_value
    auxvar%surface%P_max         = uninit_value
    auxvar%surface%coeff_for_cubic_approx(:)       = uninit_value
    auxvar%surface%coeff_for_deriv_cubic_approx(:) = uninit_value
    auxvar%surface%range_for_linear_approx(:)      = uninit_value
    auxvar%surface%dlinear_slope_dT                = uninit_value
    auxvar%surface%bcflux_default_scheme           = PETSC_FALSE
  else
    nullify(auxvar%surface)
  endif
  
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
  auxvar2%dsat_dt = auxvar%dsat_dt
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
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_dt = auxvar%dKe_dt
  if (associated(auxvar%ice)) then
    auxvar2%ice%Ke_fr = auxvar%ice%Ke_fr
    auxvar2%ice%dKe_fr_dp = auxvar%ice%dKe_fr_dp
    auxvar2%ice%dKe_fr_dt = auxvar%ice%dKe_fr_dt
    auxvar2%ice%sat_ice = auxvar%ice%sat_ice 
    auxvar2%ice%sat_gas = auxvar%ice%sat_gas
    auxvar2%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp
    auxvar2%ice%dsat_gas_dp = auxvar%ice%dsat_gas_dp
    auxvar2%ice%dsat_ice_dt = auxvar%ice%dsat_ice_dt
    auxvar2%ice%dsat_gas_dt = auxvar%ice%dsat_gas_dt
    auxvar2%ice%den_ice = auxvar%ice%den_ice
    auxvar2%ice%dden_ice_dp = auxvar%ice%dden_ice_dp
    auxvar2%ice%dden_ice_dt = auxvar%ice%dden_ice_dt
    auxvar2%ice%u_ice = auxvar%ice%u_ice
    auxvar2%ice%du_ice_dt = auxvar%ice%du_ice_dt
    auxvar2%ice%pres_fh2o = auxvar%ice%pres_fh2o
    auxvar2%ice%dpres_fh2o_dp = auxvar%ice%dpres_fh2o_dp
    auxvar2%ice%dpres_fh2o_dt = auxvar%ice%dpres_fh2o_dt
    auxvar2%ice%den_gas = auxvar%ice%den_gas
    auxvar2%ice%dden_gas_dt = auxvar%ice%dden_gas_dt
    auxvar2%ice%u_gas = auxvar%ice%u_gas
    auxvar2%ice%du_gas_dt = auxvar%ice%du_gas_dt
    auxvar2%ice%mol_gas = auxvar%ice%mol_gas
    auxvar2%ice%dmol_gas_dt = auxvar%ice%dmol_gas_dt
  endif
  if (associated(auxvar%surface)) then
    auxvar2%surface%surf_wat = auxvar%surface%surf_wat
    auxvar2%surface%P_min = auxvar%surface%P_min
    auxvar2%surface%P_max = auxvar%surface%P_max
    auxvar2%surface%coeff_for_cubic_approx(:) = &
      auxvar%surface%coeff_for_cubic_approx(:)
    auxvar2%surface%coeff_for_deriv_cubic_approx(:) = &
      auxvar%surface%coeff_for_deriv_cubic_approx(:)
    auxvar2%surface%range_for_linear_approx(:) = &
      auxvar%surface%range_for_linear_approx(:)
    auxvar2%surface%dlinear_slope_dT = auxvar%surface%dlinear_slope_dT
    auxvar2%surface%bcflux_default_scheme = &
      auxvar%surface%bcflux_default_scheme
  endif

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
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp, hw_dt
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dt
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: aux(1)

! auxvar%den = 0.d0
! auxvar%den_kg = 0.d0
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%kvr = 0.d0
  kr = 0.d0
 
! auxvar%pres = x(1)  
! auxvar%temp = x(2)
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
 
! auxvar%pc = option%reference_pressure - auxvar%pres
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)

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

  ! may need to compute dpsat_dt to pass to VISW
  call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,dpsat_dt,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dt,ierr)
  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,dpsat_dt,visl, &
                           dvis_dt,dvis_dp,ierr)
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,dpsat_dt,aux, &
                              visl,dvis_dt,dvis_dp,ierr)
  endif
  ! J/kmol -> whatever units
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
  
!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
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
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
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
#include "petsc/finclude/petscsys.h"
  use petscsys

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
  PetscReal :: dvis_dt, dvis_dp
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
  auxvar%kvr = 0.d0
   
  global_auxvar%pres = x(1)  
  global_auxvar%temp = x(2)
  
  ! Check if the capillary pressure is less than -100MPa
  
  if (global_auxvar%pres(1) - option%reference_pressure < -1.d8 + 1.d0) then
    global_auxvar%pres(1) = -1.d8 + option%reference_pressure + 1.d0
  endif

 
  auxvar%pc = option%reference_pressure - global_auxvar%pres(1)

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
                                      auxvar%ice%pres_fh2o, &
                                      auxvar%ice%dpres_fh2o_dp, &
                                      auxvar%ice%dpres_fh2o_dt, &
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

  call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dt,ierr)
  call EOSWaterEnthalpy(global_auxvar%temp,pw,hw,hw_dp,hw_dt,ierr)
  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dt = hw_dt * option%scale
                         
  call EOSWaterSaturationPressure(global_auxvar%temp, sat_pressure, &
                                  dpsat_dt, ierr)
  call EOSWaterViscosity(global_auxvar%temp, pw, sat_pressure, dpsat_dt, &
                         visl, dvis_dt,dvis_dp, ierr)

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

  auxvar%ice%sat_ice = ice_saturation
  auxvar%ice%sat_gas = gas_saturation
  auxvar%dsat_dt = dsl_temp
  auxvar%ice%dsat_ice_dp = dsi_pl
  auxvar%ice%dsat_gas_dp = dsg_pl
  auxvar%ice%dsat_ice_dt = dsi_temp
  auxvar%ice%dsat_gas_dt = dsg_temp
  
  ! Calculate the density, internal energy and derivatives for ice
  call EOSWaterDensityIce(global_auxvar%temp, global_auxvar%pres(1), &
                          den_ice, dden_ice_dT, dden_ice_dP, ierr)

  call EOSWaterInternalEnergyIce(global_auxvar%temp, u_ice, du_ice_dT)

  auxvar%ice%den_ice = den_ice
  auxvar%ice%dden_ice_dt = dden_ice_dT
  auxvar%ice%dden_ice_dp = dden_ice_dP
  auxvar%ice%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  auxvar%ice%du_ice_dt = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K 

  ! Calculate the values and derivatives for density and internal energy
  call EOSWaterSaturationPressure(global_auxvar%temp, p_sat, ierr)

  p_g            = option%reference_pressure
  auxvar%ice%den_gas = p_g/(IDEAL_GAS_CONSTANT*(global_auxvar%temp + 273.15d0))*1.d-3 !in kmol/m3
  mol_g          = p_sat/p_g
  C_g            = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  auxvar%ice%u_gas   = C_g*(global_auxvar%temp + 273.15d0)           ! in MJ/kmol
  auxvar%ice%mol_gas = mol_g

  auxvar%ice%dden_gas_dt = - p_g/(IDEAL_GAS_CONSTANT*(global_auxvar%temp + 273.15d0)**2)*1.d-3
  dmolg_dt           = dpsat_dt/p_g
  auxvar%ice%du_gas_dt   = C_g + (C_wv*dmolg_dt*FMWH2O - C_a*dmolg_dt*FMWAIR)* &
                       (global_auxvar%temp + 273.15d0)
  auxvar%ice%dmol_gas_dt = dmolg_dt

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(ithrm)
  alpha_fr = th_parameter%alpha_fr(ithrm)
  Dk = th_parameter%ckwet(ithrm)
  Dk_dry = th_parameter%ckdry(ithrm)
  Dk_ice = th_parameter%ckfrozen(ithrm)

  !Soil Kersten number
  Ke = (global_auxvar%sat(1) + epsilon)**(alpha)
  Ke_fr = (auxvar%ice%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%ice%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  auxvar%Dk_eff = Dk*Ke + Dk_ice*Ke_fr + (1.d0 - Ke - Ke_fr)*Dk_dry

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dt = alpha*(global_auxvar%sat(1) + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dt
  auxvar%ice%dKe_fr_dt = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dt
  auxvar%ice%dKe_fr_dp = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dp

  if (option%ice_model == DALL_AMICO) then
    auxvar%ice%den_ice = dw_mol
    auxvar%ice%dden_ice_dt = auxvar%dden_dt
    auxvar%ice%dden_ice_dp = auxvar%dden_dp
!    auxvar%ice%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%ice%du_ice_dt = auxvar%du_dt

    auxvar%ice%sat_gas       = 0.d0
    auxvar%ice%dsat_gas_dp   = 0.d0
    auxvar%ice%dsat_gas_dt   = 0.d0
    auxvar%ice%den_gas       = 0.d0
    auxvar%ice%dden_gas_dt   = 0.d0
    auxvar%ice%u_gas         = 0.d0
    auxvar%ice%du_gas_dt     = 0.d0
    auxvar%ice%mol_gas       = 0.d0
    auxvar%ice%dmol_gas_dt   = 0.d0
  endif

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !

subroutine THAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a TH auxiliary object
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  implicit none

  type(TH_auxvar_type) :: auxvar
  
  if (associated(auxvar%ice)) deallocate(auxvar%ice)
  nullify(auxvar%ice)
  if (associated(auxvar%surface)) deallocate(auxvar%surface)
  nullify(auxvar%surface)
  
end subroutine THAuxVarDestroy

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
    call THAuxVarDestroy(aux%auxvars(iaux))
  enddo  
  do iaux = 1, aux%num_aux_bc
    call THAuxVarDestroy(aux%auxvars_bc(iaux))
  enddo  
  do iaux = 1, aux%num_aux_ss
    call THAuxVarDestroy(aux%auxvars_ss(iaux))
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
