module PM_TOWG_Aux_module

  use PFLOTRAN_Constants_module

  use PM_Base_Aux_module
  use AuxVars_TOWG_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  !global variable to TOWG
  PetscReal, public :: towg_window_epsilon = 1.d-4
  !PetscReal, public :: towg_fmw_comp(3) = & ! initialised after EOSread
  !           [UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE,UNINITIALIZED_DOUBLE]
  PetscReal, pointer, public :: towg_fmw_comp(:)
  PetscInt, pointer, public :: towg_dof_to_primary_variable(:,:)

  PetscInt, public :: towg_debug_cell_id = UNINITIALIZED_INTEGER
  PetscReal, parameter, public :: towg_pressure_scale = 1.d0

  PetscBool, public :: towg_isothermal = PETSC_FALSE
  !PO:needs to add input for towg_no_gas and towg_no_oil in pm_towg%Read
  !   towg_no_oil currently supported only for TOWG_IMMISCIBLE 
  !   and TOWG_TODD_LONGSTAFF. To have it working also for BLACK_OIL and SOLV.
  !   must swap the orger of primary vars for TOWG_LIQ_GAS_STATE,
  !   Po,Sg,Xg^G -> Po,Xg^G,Sg
  PetscBool, public :: towg_no_gas = PETSC_FALSE
  PetscBool, public :: towg_no_oil = PETSC_FALSE 
  PetscInt, public :: towg_miscibility_model = UNINITIALIZED_INTEGER

  !list of TOWG paramters
  PetscInt, parameter, public :: TOWG_PREV_TS = 1
  PetscInt, parameter, public :: TOWG_PREV_IT = 2

  !available miscibility models
  PetscInt, parameter, public :: TOWG_IMMISCIBLE = 1
  PetscInt, parameter, public :: TOWG_TODD_LONGSTAFF = 2
  PetscInt, parameter, public :: TOWG_BLACK_OIL = 3
  PetscInt, parameter, public :: TOWG_SOLVENT_TL = 4

  ! thermodynamic state of fluid ids - for BLACK OIL 
  PetscInt, parameter, public :: TOWG_NULL_STATE = 0
  PetscInt, parameter, public :: TOWG_LIQ_OIL_STATE = 1
  PetscInt, parameter, public :: TOWG_LIQ_GAS_STATE = 2
  PetscInt, parameter, public :: TOWG_THREE_PHASE_STATE = 3
  PetscInt, parameter, public :: TOWG_ANY_STATE = 4

  ! Primary DOF indices - 
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_DOF = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_2PH_DOF = 2
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_3PH_DOF = 3
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_DOF = 3
  PetscInt, parameter, public :: TOWG_X_OIL_IN_GAS_DOF = 3
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_DOF = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_DOF = 5 
  !towg_energy_dof assigned TOWG_3CMPS_ENERGY_DOF or TOWG_SOLV_TL_ENERGY_DOF
  !in PMTOWGCreate 
  PetscInt, public :: towg_energy_dof = UNINITIALIZED_INTEGER 
  

  ! Equation indices -   
  PetscInt, parameter, public :: TOWG_LIQ_EQ_IDX = 1
  PetscInt, parameter, public :: TOWG_OIL_EQ_IDX = 2
  PetscInt, parameter, public :: TOWG_GAS_EQ_IDX = 3
  PetscInt, parameter, public :: TOWG_SOLV_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_3CMPS_ENERGY_EQ_IDX = 4
  PetscInt, parameter, public :: TOWG_SOLV_TL_ENERGY_EQ_IDX = 5
  PetscInt, public :: towg_energy_eq_idx = UNINITIALIZED_INTEGER

  !Indices used to map aux_real for condition values 
  PetscInt, parameter, public :: TOWG_OIL_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: TOWG_OIL_SATURATION_INDEX = 2 
  PetscInt, parameter, public :: TOWG_GAS_SATURATION_INDEX = 3  
  PetscInt, parameter, public :: TOWG_SOLV_SATURATION_INDEX = 4
  PetscInt, parameter, public :: TOWG_X_GAS_IN_OIL_INDEX = 5
  PetscInt, parameter, public :: TOWG_X_GAS_IN_GAS_INDEX = 6
  PetscInt, parameter, public :: TOWG_TEMPERATURE_INDEX = 7
  PetscInt, parameter, public :: TOWG_LIQUID_FLUX_INDEX = 8
  PetscInt, parameter, public :: TOWG_OIL_FLUX_INDEX = 9
  PetscInt, parameter, public :: TOWG_GAS_FLUX_INDEX = 10
  PetscInt, parameter, public :: TOWG_SOLV_FLUX_INDEX = 11
  PetscInt, parameter, public :: TOWG_ENERGY_FLUX_INDEX = 12
  PetscInt, parameter, public :: TOWG_LIQ_CONDUCTANCE_INDEX = 13
  PetscInt, parameter, public :: TOWG_OIL_CONDUCTANCE_INDEX = 14
  PetscInt, parameter, public :: TOWG_GAS_CONDUCTANCE_INDEX = 15
  PetscInt, parameter, public :: TOWG_MAX_INDEX = 15

  !Indices used to map aux_int_var for condition values
  PetscInt, parameter, public :: TOWG_STATE_INDEX = 1

  !flags to identify type of auxvar update
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TOWG_UPDATE_FOR_BOUNDARY = 2

  ! it might be required for thermal diffusion terms and tough conv criteria
  type, public :: towg_parameter_type
     !  PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
     !  PetscReal :: newton_inf_scaled_res_tol
     !  PetscBool :: check_post_converged
  end type towg_parameter_type

  !if required, could add other intermediate classes:
  type, public, extends(pm_base_aux_type) :: pm_towg_aux_type 
    type(towg_parameter_type), pointer :: parameter
    class(auxvar_towg_type), pointer :: auxvars(:,:)
    class(auxvar_towg_type), pointer :: auxvars_bc(:)
    class(auxvar_towg_type), pointer :: auxvars_ss(:)
  contains
    !add bound-procedure
    procedure, public :: Init => InitTOWGAuxVars
    !procedure, public :: Perturb => PerturbTOilIms
  end type pm_towg_aux_type

  interface TOWGAuxVarStrip
    module procedure TOWGAuxVarArray1Strip
    module procedure TOWGAuxVarArray2Strip
  end interface TOWGAuxVarStrip


  !pointing to null() function
  procedure(TOWGAuxVarComputeDummy), pointer :: TOWGAuxVarCompute => null()
  procedure(TOWGAuxVarPerturbDummy), pointer :: TOWGAuxVarPerturb => null()

  abstract interface
    subroutine TOWGAuxVarComputeDummy(x,auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
      use Option_module
      use AuxVars_TOWG_module
      use Global_Aux_module
      use EOS_Water_module
      use Characteristic_Curves_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      class(characteristic_curves_type) :: characteristic_curves
      PetscReal :: x(option%nflowdof)
      class(auxvar_towg_type) :: auxvar
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      PetscInt :: natural_id
    end subroutine TOWGAuxVarComputeDummy

    subroutine TOWGAuxVarPerturbDummy(auxvar,global_auxvar, &
                                      material_auxvar, &
                                      characteristic_curves,natural_id, &
                                      option)
      use AuxVars_TOWG_module
      use Option_module
      use Characteristic_Curves_module
      use Global_Aux_module
      use Material_Aux_class
      implicit none
      type(option_type) :: option
      PetscInt :: natural_id
      type(auxvar_towg_type) :: auxvar(0:)
      type(global_auxvar_type) :: global_auxvar
      class(material_auxvar_type) :: material_auxvar
      class(characteristic_curves_type) :: characteristic_curves
    end subroutine 

  end interface
 
  public :: TOWGAuxCreate, &
            TOWGAuxDestroy, &
            TOWGAuxVarStrip, &
            TOWGAuxVarCompute, &
            TOWGImsAuxVarComputeSetup, &
            TOWGAuxVarPerturb
  !          TOilImsAuxVarPerturb, TOilImsAuxDestroy, &
  !          TOilImsAuxVarStrip

contains

! ************************************************************************** !

function TOWGAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object for TOWG
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/05/16
  ! 

  use Option_module
  use EOS_Oil_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
    
  class(pm_towg_aux_type), pointer :: TOWGAuxCreate

  class(pm_towg_aux_type), pointer :: aux

  allocate(towg_fmw_comp(option%nflowspec)) 

  !in TOWG the gas FMW must be defined in the input deck
  if ( Uninitialized(EOSGasGetFMW()) ) then
    option%io_buffer = 'TOWG: gas FMW not initialised. ' // &
                       'Define its value in the the input deck' // &
                       ' or add EOS GAS card to default to FMWAIR' 
    call printErrMsg(option)    
  endif
  
  towg_fmw_comp(1) = FMWH2O
  towg_fmw_comp(2) = EOSOilGetFMW()
  towg_fmw_comp(3) = EOSGasGetFMW()

  !need to add an EOS for solvent, before adding solvent 
  !if ( towg_miscibility_model == TOWG_SOLVENT_TL ) then
  !  towg_fmw_comp(4) = 
  !end if   

  allocate( towg_dof_to_primary_variable(1:option%nflowdof,1:3) )

  select case(towg_miscibility_model)
    case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF,TOWG_BLACK_OIL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_X_GAS_IN_OIL_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, &
                TOWG_X_GAS_IN_GAS_INDEX, TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_GAS_SATURATION_INDEX,TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
    case(TOWG_SOLVENT_TL)
      towg_dof_to_primary_variable(1:option%nflowdof,1:3) = &
        reshape([TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_X_GAS_IN_OIL_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_GAS_SATURATION_INDEX, &
                TOWG_X_GAS_IN_GAS_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX, &
                TOWG_OIL_PRESSURE_INDEX, TOWG_OIL_SATURATION_INDEX, &
                TOWG_GAS_SATURATION_INDEX, TOWG_SOLV_SATURATION_INDEX, &
                TOWG_TEMPERATURE_INDEX], &
                shape(towg_dof_to_primary_variable))
  end select

  !allocate here to define this is a pm_toil_ims_aux_type
  allocate(aux) 

  call PMBaseAuxInit(aux) 

  !nullify here and not in the parent class, because auxvars are mode dependent
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  allocate(aux%parameter)
  
  !PO - to allocate when supporting diffusion the oil and gas phase
  !allocate(aux%parameter%diffusion_coefficient(option%nphase)) 
  !aux%parameter%diffusion_coefficient(option%liquid_phase) = & 
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%oil_phase) = &
  !                                                 UNINITIALIZED_DOUBLE
  !aux%parameter%diffusion_coefficient(option%gas_phase) = 2.13d-5
  !aux%parameter%newton_inf_scaled_res_tol = 1.d-50
  !aux%parameter%check_post_converged = PETSC_FALSE


  TOWGAuxCreate => aux
  
end function TOWGAuxCreate

! ************************************************************************** !

subroutine InitTOWGAuxVars(this,grid,num_bc_connection, &
                              num_ss_connection,option)
  ! 
  ! Initialize pm_towg_auxvars  
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/07/16
  ! 

  use Option_module
  use Grid_module 

  implicit none

  class(pm_towg_aux_type) :: this
  PetscInt :: num_bc_connection
  PetscInt :: num_ss_connection  
  type(grid_type) :: grid
  type(option_type) :: option

  PetscInt :: ghosted_id, iconn, local_id
  PetscInt :: idof 
 
  allocate(this%auxvars(0:option%nflowdof,grid%ngmax)) 
  do ghosted_id = 1, grid%ngmax
    do idof = 0, option%nflowdof
      call this%auxvars(idof,ghosted_id)%Init(option)
    enddo
  enddo

  this%num_aux = grid%ngmax

  if (num_bc_connection > 0) then
    allocate(this%auxvars_bc(num_bc_connection))
    do iconn = 1, num_bc_connection
      call this%auxvars_bc(iconn)%Init(option)
    enddo
  endif
  this%num_aux_bc = num_bc_connection

  if (num_ss_connection > 0) then
    allocate(this%auxvars_ss(num_ss_connection))
    do iconn = 1, num_ss_connection
      call this%auxvars_ss(iconn)%Init(option)
    enddo
  endif
  this%num_aux_ss = num_ss_connection

  call PMBaseAuxSetup(this,grid,option)

end subroutine InitTOWGAuxVars

! ************************************************************************** !

subroutine TOWGImsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  ! 
  ! Computes auxiliary variables for each grid cell for TOWGIms
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/02/16
  ! 

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Oil_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  
  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  class(auxvar_towg_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar ! passing this for salt conc.
                                            ! not currenty used  
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id !only for debugging/print out - currently not used 

  PetscInt :: lid, oid, gid
  PetscReal :: cell_pressure, wat_sat_pres
  PetscReal :: krl, visl, dkrl_Se
  PetscReal :: kro, viso, dkro_Se
  PetscReal :: krg, visg, dkrg_Se
  PetscReal :: dummy
  !PetscReal :: Uoil_J_kg, Hoil_J_kg
  PetscErrorCode :: ierr

  ! from SubsurfaceSetFlowMode
  ! option%liquid_phase = 1           ! liquid_pressure
  ! option%oil_phase = 2              ! oil_pressure
  ! option%gas_phase = 3              ! gas_pressure
  lid = option%liquid_phase
  oid = option%oil_phase
  gid = option%gas_phase

  auxvar%effective_porosity = 0.d0
  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mobility = 0.d0
  auxvar%H = 0.d0
  auxvar%U = 0.d0
  auxvar%temp = 0.0d0

  !assing auxvars given by the solution variables
  auxvar%pres(oid) = x(TOWG_OIL_PRESSURE_DOF)
  auxvar%sat(oid) = x(TOWG_OIL_SATURATION_DOF)
  auxvar%sat(gid) = x(TOWG_GAS_SATURATION_3PH_DOF)
  auxvar%temp = x(towg_energy_dof)

  auxvar%sat(lid) = 1.d0 - auxvar%sat(oid) - auxvar%sat(gid)

  !below pc_ow /= 0
  !compute capillary presssure water/oil (pc_ow)
  !call characteristic_curves%saturation_function% &
  !           CapillaryPressure(auxvar%sat(lid),auxvar%pc(lid),option)
  !auxvar%pres(lid) = auxvar%pres(oid) - auxvar%pc(lid)

  !Assumptions below on capillary pressure for comparison with TOUGH2-EOS8:
  ! pc_ow = 0, only pc_gw /= 0 and computed considering water saturation only
  ! this results in pc_gw = pc_go
  
  !assuming no capillary pressure between oil and water: pc_ow = 0
  ! pc_ow = 0.0d0
  auxvar%pres(lid) = auxvar%pres(oid) 

  !compute capillary pressure gas/water (pc_gw),
  ! pc_go = pc_gw when pc_ow = 0
  call characteristic_curves%saturation_function% &
             CapillaryPressure(auxvar%sat(lid),auxvar%pc(lid),option)
  
  auxvar%pres(gid) = auxvar%pres(lid) + auxvar%pc(lid)

  auxvar%pc(oid) = 0.0d0


  cell_pressure = max(auxvar%pres(lid),auxvar%pres(oid),auxvar%pres(gid))


  ! calculate effective porosity as a function of pressure
  if (option%iflag /= TOWG_UPDATE_FOR_BOUNDARY) then
    auxvar%effective_porosity = material_auxvar%porosity_base

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                auxvar%effective_porosity,dummy)
    endif
    if (option%iflag /= TOWG_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

  ! UPDATE THERMODYNAMIC PROPERTIES FOR BOTH PHASES!!!

  ! Liquid phase thermodynamic properties
  ! using cell_pressure (which is the max press)? or %pres(lid)?
  call EOSWaterDensity(auxvar%temp,cell_pressure, &
                       auxvar%den_kg(lid),auxvar%den(lid),ierr)
  call EOSWaterEnthalpy(auxvar%temp,cell_pressure,auxvar%H(lid),ierr)
  auxvar%H(lid) = auxvar%H(lid) * 1.d-6 ! J/kmol -> MJ/kmol

  ! MJ/kmol comp                  ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
  auxvar%U(lid) = auxvar%H(lid) - (cell_pressure / auxvar%den(lid) * 1.d-6)

  ! ADD HERE BRINE dependency. Two options (see mphase)
  ! - salinity constant in space and time (passed in option%option%m_nacl)
  ! - salt can be trasnported by RT (sequential coupling) and passed 
  !   and passed with global_auxvar%m_nacl 
  !  ! Assign salinity 
  !  m_na=option%m_nacl; m_cl=m_na; m_nacl=m_na 
  !  if (option%ntrandof > 0) then
  !    m_na = global_auxvar%m_nacl(1)
  !    m_cl = global_auxvar%m_nacl(2)
  !    m_nacl = m_na
  !    if (m_cl > m_na) m_nacl = m_cl
  !  endif    
  !
  !  ! calculate density for pure water
  !  call EOSWaterDensityEnthalpy(t,pw,dw_kg,dw_mol,hw,ierr)
  !  !..................
  !  xm_nacl = m_nacl*FMWNACL
  !  xm_nacl = xm_nacl/(1.D3 + xm_nacl)
  !  ! corrects water densit previously calculated as pure water
  !  call EOSWaterDensityNaCl(t,p,xm_nacl,dw_kg)  
  !  ! water viscosity dependence on salt concetration, but no derivatives
  !  !  call EOSWaterViscosityNaCl(t,p,xm_nacl,visl)
  !  call EOSWaterViscosity(t,pw,sat_pressure,0.d0,visl,dvdt,dvdp,dvdps,ierr)

  call EOSOilDensityEnergy(auxvar%temp,auxvar%pres(oid),&
                           auxvar%den(oid),auxvar%H(oid), &
                           auxvar%U(oid),ierr)

  auxvar%den_kg(oid) = auxvar%den(oid) * EOSOilGetFMW()

  auxvar%H(oid) = auxvar%H(oid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(oid) = auxvar%U(oid) * 1.d-6 ! J/kmol -> MJ/kmol

  !compute gas properties (default is air - but methane can be set up)
  call EOSGasDensityEnergy(auxvar%temp,auxvar%pres(gid),auxvar%den(gid), &
                           auxvar%H(gid),auxvar%U(gid),ierr)

  auxvar%den_kg(gid) = auxvar%den(gid) * EOSGasGetFMW()
  auxvar%H(gid) = auxvar%H(gid) * 1.d-6 ! J/kmol -> MJ/kmol
  auxvar%U(gid) = auxvar%U(gid) * 1.d-6 ! J/kmol -> MJ/kmol


  ! compute water mobility (rel. perm / viscostiy)
  call characteristic_curves%liq_rel_perm_function% &
         RelativePermeability(auxvar%sat(lid),krl,dkrl_Se,option)
                            
  call EOSWaterSaturationPressure(auxvar%temp, wat_sat_pres,ierr)                   

  ! use cell_pressure; cell_pressure - psat calculated internally
  call EOSWaterViscosity(auxvar%temp,cell_pressure,wat_sat_pres,visl,ierr)

  auxvar%mobility(lid) = krl/visl


  ! compute oil mobility (rel. perm / viscostiy)
  call characteristic_curves%oil_rel_perm_function% &
         RelativePermeability(auxvar%sat(lid),kro,dkro_Se,option)

  call EOSOilViscosity(auxvar%temp,auxvar%pres(oid), &
                       auxvar%den(oid), viso, ierr)

  auxvar%mobility(oid) = kro/viso

  ! compute gas mobility (rel. perm / viscosity)
  call characteristic_curves%gas_rel_perm_function% &
         RelativePermeability(auxvar%sat(lid),krg,dkrg_Se,option)

  !currently only a viscosity model for air or constant value   
  call EOSGasViscosity(auxvar%temp,auxvar%pres(gid), &
                       auxvar%pres(gid),auxvar%den(gid),visg,ierr)

  auxvar%mobility(gid) = krg/visg


end subroutine TOWGImsAuxVarCompute

! ************************************************************************** !

subroutine TOWGImsTLAuxVarPerturb(auxvar,global_auxvar, &
                                  material_auxvar, &
                                  characteristic_curves,natural_id, &
                                  option)
  ! 
  ! Calculates auxiliary variables for perturbed system
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 12/27/16
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(auxvar_towg_type) :: auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves
     
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), &
               pert(option%nflowdof), x_pert_save(option%nflowdof)

  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: tempreal
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: min_perturbation = 1.d-10
  PetscInt :: idof

  x(TOWG_OIL_PRESSURE_DOF) = auxvar(ZERO_INTEGER)%pres(option%oil_phase)
  x(TOWG_OIL_SATURATION_DOF) = auxvar(ZERO_INTEGER)%sat(option%oil_phase)
  x(TOWG_GAS_SATURATION_3PH_DOF) = auxvar(ZERO_INTEGER)%sat(option%gas_phase)
  x(TOWG_3CMPS_ENERGY_DOF) = auxvar(ZERO_INTEGER)%temp
  
  pert(TOWG_OIL_PRESSURE_DOF) = &
    perturbation_tolerance*x(TOWG_OIL_PRESSURE_DOF)+min_perturbation
  pert(TOWG_3CMPS_ENERGY_DOF) = &
    perturbation_tolerance*x(TOWG_3CMPS_ENERGY_DOF)+min_perturbation
  if (x(TOWG_OIL_SATURATION_DOF) > 0.5d0) then 
    pert(TOWG_OIL_SATURATION_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_OIL_SATURATION_DOF) = perturbation_tolerance
  endif
  if (x(TOWG_GAS_SATURATION_3PH_DOF) > 0.5d0) then 
    pert(TOWG_GAS_SATURATION_3PH_DOF) = -1.d0 * perturbation_tolerance
  else
    pert(TOWG_GAS_SATURATION_3PH_DOF) = perturbation_tolerance
  endif

  ! TOWG_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = TOWG_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    auxvar(idof)%pert = pert(idof)
    x_pert = x
    x_pert(idof) = x(idof) + pert(idof)
    x_pert_save = x_pert
    call TOWGAuxVarCompute(x_pert,auxvar(idof),global_auxvar, &
                           material_auxvar, &
                           characteristic_curves,natural_id,option)
  enddo

  auxvar(TOWG_OIL_PRESSURE_DOF)%pert = &
     auxvar(TOWG_OIL_PRESSURE_DOF)%pert / towg_pressure_scale
 
end subroutine TOWGImsTLAuxVarPerturb

! ************************************************************************** !
subroutine TOWGImsAuxVarComputeSetup()

  implicit none

  TOWGAuxVarCompute => TOWGImsAuxVarCompute
  TOWGAuxVarPerturb => TOWGImsTLAuxVarPerturb

end subroutine TOWGImsAuxVarComputeSetup

! ************************************************************************** !

subroutine TOWGAuxDestroy(aux)
  ! 
  ! Deallocates a towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_towg_aux_type), pointer :: aux
  PetscInt :: iaux, idof
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars) ) then
    call TOWGAuxVarStrip(aux%auxvars)
    deallocate(aux%auxvars)
  end if 
  nullify(aux%auxvars) 

  if (associated(aux%auxvars_bc) ) then
    call TOWGAuxVarStrip(aux%auxvars_bc)
    deallocate(aux%auxvars_bc)
  end if
  nullify(aux%auxvars_bc)

  if ( associated(aux%auxvars_ss) ) then
    call TOWGAuxVarStrip(aux%auxvars_ss)
    deallocate(aux%auxvars_ss)
  end if
  nullify(aux%auxvars_ss)

  call PMBaseAuxStrip(aux)

  if (associated(aux%parameter)) then
    deallocate(aux%parameter) 
    !to add paramter strip when introducing diff in oil and gas phases
  end if
  nullify(aux%parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine TOWGAuxDestroy

! ************************************************************************** !

! ************************************************************************** !

subroutine  TOWGAuxVarArray1Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !here we can pass by pointer, we could destroy the array within the routine
  !but we don't to be consistent with TOilImsAuxVarArray2Strip 
  !class(auxvar_towg_type), pointer :: auxvars(:)
  type(auxvar_towg_type) :: auxvars(:)

  PetscInt :: iaux

  !print *, "den oil bc/ss pass = ", auxvars(1)%den(2)
  
  do iaux = 1, size(auxvars)
    call auxvars(iaux)%Strip
  enddo  

end subroutine TOWGAuxVarArray1Strip

! ************************************************************************** !

subroutine TOWGAuxVarArray2Strip(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! this could be generalised for different modes 
  ! using class(*) (unlimited polymorphic)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 

  use AuxVars_TOWG_module

  implicit none

  !can't use class due to gfortran (4.8.4) bug (values passed are not correct)
  !bug fixed in gfortran 6.2.0
  !cannot use type(...) with pointer attribute, therefore we deallocate and 
  !nullify pointer outide this routine 
  !because the compiler does not allow to specify lower 0-bound in auxvar
  !type(auxvar_towg_type), pointer :: auxvars(:,:)
  !class(auxvar_towg_type) :: auxvars(0:,:)
  type(auxvar_towg_type) :: auxvars(0:,:)

  PetscInt :: iaux, idof

  do iaux = 1, size(auxvars,2)
    do idof = 1, size(auxvars,1)
      call auxvars(idof-1,iaux)%Strip()
    enddo
  enddo  

end subroutine TOWGAuxVarArray2Strip

! ************************************************************************** !

end module PM_TOWG_Aux_module
